//! Code related to calculating acceleration between bodies using the Barnes-Hut
//! approach. It groups source masses. O(N Log(N)) complexity. Faster than naive
//! n-body; slower than Fast Multipole (FMM). Using in place of FMM due to its relative
//! ease of implementation.
//!
//! todo: If you can make this general enough, publish it as a standalone library.

// todo: You can use rayon more throughout this, e.g. during tree construction.

use lin_alg::f64::Vec3;
use rayon::prelude::*;

use crate::{
    accel::{acc_newton_inner, MondFn},
    units::A0_MOND,
    Body,
};

const MAX_BODIES_PER_NODE: usize = 1;

#[derive(Clone, Debug)]
/// A cubical bounding box. length=width=depth.
pub struct Cube {
    pub center: Vec3,
    pub width: f64,
    // These mins and maxes are derivative of center/width.
    pub x_min: f64,
    pub x_max: f64,
    pub y_min: f64,
    pub y_max: f64,
    pub z_min: f64,
    pub z_max: f64,
}

impl Cube {
    /// Construct minimum limits that encompass all bodies. Run these each time the bodies change,
    /// or perhaps use a pad and do it at a coarser interval.
    ///
    /// The pad allows us to keep the same cube for multiple timesteps, but taking into acacount
    /// potential movement of bodies outside the cube between these updates.
    ///
    /// The z offset is intended for the case where the Z coordinate for all particles is 0.
    /// This prevents the divisions straddling the points, doubling the number of nodes.
    pub fn from_bodies(bodies: &[Body], pad: f64,  z_offset: bool) -> Option<Self> {
        if bodies.is_empty() {
            return None;
        }

        let mut x_min = f64::MAX;
        let mut x_max = f64::MIN;
        let mut y_min = f64::MAX;
        let mut y_max = f64::MIN;
        let mut z_min = f64::MAX;
        let mut z_max = f64::MIN;

        for body in bodies {
            let p = &body.posit;
            x_min = x_min.min(p.x);
            x_max = x_max.max(p.x);
            y_min = y_min.min(p.y);
            y_max = y_max.max(p.y);
            z_min = z_min.min(p.z);
            z_max = z_max.max(p.z);
        }

        x_min -= pad;
        x_max += pad;
        y_min -= pad;
        y_max += pad;
        z_min -= pad;
        z_max += pad;

        if z_offset {
            z_max += 1e-5;
        }

        let x_size = x_max - x_min;
        let y_size = y_max - y_min;
        let z_size = z_max - z_min;

        // Coerce to a cube.
        let mut width = x_size.max(y_size).max(z_size);
        let width_div2 = width / 2.;

        let center = Vec3::new(
            (x_max + x_min) / 2.,
            (y_max + y_min) / 2.,
            (z_max + z_min) / 2.,
        );

        Some(Self::new(center, width, width_div2))
    }

    pub fn new(center: Vec3, width: f64, width_div2: f64) -> Self {
        let x_min = center.x - width_div2;
        let x_max = center.x + width_div2;
        let y_min = center.y - width_div2;
        let y_max = center.y + width_div2;
        let z_min = center.z - width_div2;
        let z_max = center.z + width_div2;

        Self {
            center,
            width,
            x_min,
            x_max,
            y_min,
            y_max,
            z_min,
            z_max,
        }
    }

    /// Divide this into equal-area octants.
    pub(crate) fn divide_into_octants(&self) -> [Self; 8] {
        let width = self.width / 2.;
        let wd2 = self.width / 4.; // short for brevity below.

        // Every combination of + and - for the center offset.
        // The order matters, due to the binary index logic used when partitioning bodies into octants.
        [
            Self::new(self.center + Vec3::new(-wd2, -wd2, -wd2), width, wd2),
            Self::new(self.center + Vec3::new(wd2, -wd2, -wd2), width, wd2),
            Self::new(self.center + Vec3::new(-wd2, wd2, -wd2), width, wd2),
            Self::new(self.center + Vec3::new(wd2, wd2, -wd2), width, wd2),
            Self::new(self.center + Vec3::new(-wd2, -wd2, wd2), width, wd2),
            Self::new(self.center + Vec3::new(wd2, -wd2, wd2), width, wd2),
            Self::new(self.center + Vec3::new(-wd2, wd2, wd2), width, wd2),
            Self::new(self.center + Vec3::new(wd2, wd2, wd2), width, wd2),
        ]
    }
}

#[derive(Debug)]
// todo: Implement with flat structure.
pub struct Node {
    pub id: usize, // todo: Remove A/R. Ideally, we use indices.
    pub bounding_box: Cube,
    /// Node indices in the tree. We use this to guide the transversal process while finding
    /// relevant nodes for a given target body.
    pub children: Vec<usize>,
    mass: f64,
    center_of_mass: Vec3,
}

#[derive(Debug)]
/// A recursive tree. Each node can be subdivided  Terminates with `NodeType::NodeTerminal`.
pub struct Tree {
    pub nodes: Vec<Node>,
}

impl Tree {
    /// Constructs a tree. Call this externaly using all bodies, once per time step.
    /// It creates the entire tree, branching until each cell has `MAX_BODIES_PER_NODE` or fewer
    /// bodies.
    ///
    /// We partially transverse it as-required while calculating the force on a given target.
    pub fn new(
        bodies: &[Body],
        bb: &Cube,
    ) -> Self {
        // Convert &[Body] to &[&Body], and remove the target body.
        let body_refs: Vec<&Body> = bodies
            .iter()
            .collect();

        let mut nodes = Vec::new();

        let mut current_node_i: usize = 0;
        populate_nodes(&mut nodes, &body_refs, bb, &mut current_node_i);

        Self { nodes }
    }


    /// Recursive function called when getting leaves.
    fn leaves_inner<'a>(&'a self, leaves: &mut Vec<&'a Node>, node_i: usize, posit_target: Vec3, id_target: usize, θ: f64) {
        let mut node = &self.nodes[node_i];

        let dist = (posit_target - node.center_of_mass).magnitude();

        // todo: Is this where you should prevent self-interaction? Likely.

        if dist < 1e-8 {
            // todo: Better way, using id_target?
            println!("Attempting to avoid self-interaction");
            return;
        }

        // The source is far, or we have few items (e.g. 1) items in this region; this node is
        // terminal. (A leaf).
        if node.children.len() <= MAX_BODIES_PER_NODE || node.bounding_box.width / dist < θ {
            leaves.push(node);
        } else {
            // The source is near; go deeper.
            for child_i in &node.children {
                self.leaves_inner(leaves, child_i, posit_target, id_target, θ);
            }
        }
    }

    /// Get all leaves relevant to a given target. We use this to create a coarser
    /// version of the tree, containing only the nodes we need to calculate acceleration
    /// on a specific target.
    pub(crate) fn leaves(&self, posit_target: Vec3, id_target: usize, θ: f64) -> Vec<&Node> {
        let mut result = Vec::new();

        if self.nodes.is_empty() {
            return result;
        }

        self.leaves_inner(&mut result, 0, posit_target, id_target, θ);

        result
    }
}

/// Compute center of mass as a position, and mass value.
fn center_of_mass(bodies: &[&Body]) -> (Vec3, f64) {
    let mut mass = 0.;
    let mut center_of_mass = Vec3::new_zero();

    for body in bodies {
        mass += body.mass;
        // Weight the center by mass.
        center_of_mass += body.posit * body.mass;
    }
    if mass.abs() > f64::EPSILON {
        // Remove the weighting's effect on magnitude of the center.
        center_of_mass /= mass;
    }

    (center_of_mass, mass)
}


/// Recursively populate nodes of the tree. This contains the bulk of the Barnes Hut algorithm.
/// This is effectively the `Tree` constructor, but is separate from it so we can recur.
fn populate_nodes(nodes: &mut Vec<Node>, bodies: &[&Body], bb: &Cube, current_node_i: &mut usize) {
    // todo: Is this at the core? At least the initial one can be cached.
    //todo: I suspect the fix is a more fundamental change than that.
    let (center_of_mass, mass) = center_of_mass(bodies);

    if bodies.len() <= MAX_BODIES_PER_NODE {
        nodes.push(Node {
            id: *current_node_i,
            bounding_box: bb.clone(), // todo: Way around this?
            mass,
            center_of_mass,
            children: Vec::new(),
        });
            *current_node_i += 1;
    } else {
        // Populate up to 8 children.
        let octants = bb.divide_into_octants();
        let bodies_by_octant = partition(bodies, bb);

        let node_i_this =  *current_node_i;

        let mut children = Vec::new();
        for bodies in &bodies_by_octant {
            if !bodies.is_empty() {
                // todo QC off-by-one etc.

                // todo: Wrong.
                *current_node_i += 1;
                children.push(*current_node_i);
            }
        }

        nodes.push(Node {
            id: node_i_this,
            bounding_box: bb.clone(), // todo: Way around this?
            mass,
            center_of_mass,
            children,
        });

        for (i, octant) in octants.into_iter().enumerate() {
            if !bodies_by_octant[i].is_empty() {
                populate_nodes(
                    nodes,
                    bodies_by_octant[i].as_slice(),
                    &octant,
                    current_node_i,
                )
            }
        }
        // octants
        //     .into_par_iter()
        //     .zip(bodies_by_octant.into_par_iter())
        //     .filter(|(_, bodies)| !bodies.is_empty())
        //     .for_each(|(octant, bodies)| {
        //         populate_nodes(nodes, bodies.as_slice(), &octant, posit_target, θ);
        //     });
        // }
    }
}

/// Partition bodies into each of the 8 octants.
fn partition<'a>(bodies: &[&'a Body], bb: &Cube) -> [Vec<&'a Body>; 8] {
    let mut result: [Vec<&Body>; 8] = Default::default();

    for body in bodies {
        let mut index = 0;
        if body.posit.x > bb.center.x {
            index |= 0b001;
        }
        if body.posit.y > bb.center.y {
            index |= 0b010;
        }
        if body.posit.z > bb.center.z {
            index |= 0b100;
        }

        result[index].push(body);
    }

    result
}

/// Calculate Newtonian acceleration using the Barnes Hut algorithm.
pub fn acc_newton_bh(
    posit_target: Vec3,
    id_target: usize,
    tree: &Tree,
    θ: f64,
    mond: Option<MondFn>,
    softening_factor_sq: f64,
) -> Vec3 {
    // todo: Put back the part checking for self interaction.
    tree.leaves(posit_target, id_target, θ)
        // tree.nodes
        // .par_chunks(32) // Process in chunks
        // .map(|chunk| {
        //     chunk.iter().map(|leaf| {
        .par_iter()
        // .iter()
        .map(|leaf| {
            let acc_diff = leaf.center_of_mass - posit_target;
            let dist = acc_diff.magnitude();
            let acc_dir = acc_diff / dist; // Unit vec

            let mut acc = acc_newton_inner(acc_dir, leaf.mass, dist, softening_factor_sq);

            if let Some(mond_fn) = mond {
                let x = acc.magnitude() / A0_MOND;
                acc = acc / mond_fn.μ(x);
            }

            acc
        })
        .reduce(Vec3::new_zero, |acc, elem| acc + elem)
    // .fold(Vec3::new_zero(), |acc, elem| acc + elem)
}