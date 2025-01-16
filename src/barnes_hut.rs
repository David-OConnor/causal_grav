//! Code related to calculating acceleration between bodies using the Barnes-Hut
//! approach. It groups source masses. O(N Log(N)) complexity. Faster than naive
//! n-body; slower than Fast Multipole (FMM). Using in place of FMM due to its relative
//! ease of implementation.
//!
//! todo: If you can make this general enough, publish it as a standalone library.

// todo: You can use rayon more throughout this, e.g. during tree construction.

use std::{fmt, fmt::Formatter};

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
    pub fn from_bodies(bodies: &[Body], pad: f64, z_offset: bool) -> Option<Self> {
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
    // pub(crate) fn divide_into_octants(&self) -> [Self; 8] { // todo: TS stack overflow
    pub(crate) fn divide_into_octants(&self) -> Vec<Self> {
        let width = self.width / 2.;
        let wd2 = self.width / 4.; // short for brevity below.

        // Every combination of + and - for the center offset.
        // The order matters, due to the binary index logic used when partitioning bodies into octants.
        // [
        vec![
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
    /// We use `id` while building the tree, then sort by it, replacing with index.
    /// Once complete, `id` == index in `Tree::nodes`.
    pub id: usize,
    pub bounding_box: Cube,
    /// Node indices in the tree. We use this to guide the transversal process while finding
    /// relevant nodes for a given target body.
    pub children: Vec<usize>,
    mass: f64,
    center_of_mass: Vec3,
}

impl fmt::Display for Node {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "Id: {}, Width: {:.3}, Ch: {:?}",
            self.id, self.bounding_box.width, self.children
        )
    }
}

#[derive(Debug)]
/// A recursive tree. Each node can be subdivided  Terminates with `NodeType::NodeTerminal`.
pub struct Tree {
    /// Order matters; we index this by `Node::children`.
    pub nodes: Vec<Node>,
}

impl Tree {
    /// Constructs a tree. Call this externaly using all bodies, once per time step.
    /// It creates the entire tree, branching until each cell has `MAX_BODIES_PER_NODE` or fewer
    /// bodies.
    ///
    /// We partially transverse it as-required while calculating the force on a given target.
    pub fn new(bodies: &[Body], bb: &Cube) -> Self {
        // Convert &[Body] to &[&Body].
        let body_refs: Vec<&Body> = bodies.iter().collect();

        let mut nodes = Vec::new();

        let mut current_node_i: usize = 0;
        // populate_nodes(&mut nodes, &body_refs, bb, &mut current_node_i);

        // Stack to simulate recursion: Each entry contains (bodies, bounding box, parent_id, child_index).
        let mut stack = Vec::new();
        stack.push((body_refs.to_vec(), bb.clone(), None));

        while let Some((bodies_, bb, parent_id)) = stack.pop() {
            let (center_of_mass, mass) = center_of_mass(&bodies_);

            if bodies_.len() <= MAX_BODIES_PER_NODE {
                // Create a leaf node.
                let node_id = current_node_i;
                nodes.push(Node {
                    id: node_id,
                    bounding_box: bb.clone(),
                    mass,
                    center_of_mass,
                    children: Vec::new(),
                });
                current_node_i += 1;

                if let Some(pid) = parent_id {
                    // nodes[pid].children.push(node_id);
                    // todo: Odd: Rust is requesting an explicit type...
                    let a: &mut Node = &mut nodes[pid];
                    a.children.push(node_id);
                }
            } else {
                // Create an internal node and push its ID.
                let node_id = current_node_i;
                current_node_i += 1;

                nodes.push(Node {
                    id: node_id,
                    bounding_box: bb.clone(),
                    mass,
                    center_of_mass,
                    children: Vec::new(),
                });

                if let Some(pid) = parent_id {
                    // nodes[pid].children.push(node_id);
                    // todo: Odd: Rust is requesting an explicit type...
                    let a: &mut Node = &mut nodes[pid];
                    a.children.push(node_id);
                }

                // Divide into octants and partition bodies.
                let octants = bb.divide_into_octants();
                let bodies_by_octant = partition(&bodies_, &bb);

                // Add each octant with bodies to the stack.
                for (i, octant) in octants.into_iter().enumerate() {
                    if !bodies_by_octant[i].is_empty() {
                        stack.push((bodies_by_octant[i].clone(), octant, Some(node_id)));
                    }
                }
            }
        }

        // Now that nodes are populated, rearrange so index == `id`. We will then index by `children`.
        nodes.sort_by(|l, r| l.id.partial_cmp(&r.id).unwrap());

        Self { nodes }
    }

    /// Get all leaves relevant to a given target. We use this to create a coarser
    /// version of the tree, containing only the nodes we need to calculate acceleration
    /// on a specific target.
    pub(crate) fn leaves(&self, posit_target: Vec3, id_target: usize, θ: f64) -> Vec<&Node> {
        let mut result = Vec::new();

        if self.nodes.is_empty() {
            return result;
        }

        let node_i = 0;

        let mut stack = Vec::new();
        stack.push(node_i);

        while let Some(current_node_i) = stack.pop() {
            let node = &self.nodes[current_node_i];

            if node.children.len() <= MAX_BODIES_PER_NODE {
                result.push(node);
                continue;
            }

            let dist = (posit_target - node.center_of_mass).magnitude();

            // Avoid self-interaction based on distance or id_target.
            // todo: Use id_target, if able.
            if dist < 1e-8 {
                continue;
            }

            if node.bounding_box.width / dist < θ {
                result.push(node);
            } else {
                // The source is near; add children to the stack to go deeper.
                for &child_i in &node.children {
                    stack.push(child_i);
                }
            }
        }

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

/// Partition bodies into each of the 8 octants.
// fn partition<'a>(bodies: &[&'a Body], bb: &Cube) -> [Vec<&'a Body>; 8] { // todo TS stack overflow
fn partition<'a>(bodies: &[&'a Body], bb: &Cube) -> Vec<Vec<&'a Body>> {
    // let mut result: [Vec<&Body>; 8] = Default::default();
    let mut result: Vec<Vec<&Body>> = vec![Vec::new(); 8];

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
        .par_iter()
        .filter_map(|leaf| {
            let acc_diff = leaf.center_of_mass - posit_target;
            let dist = acc_diff.magnitude();

            // todo: Not sure why we get 0 here when we check it in `leaves()`.
            // todo: QC this.
            if dist < 1e-8 {
                return None;
            }

            let acc_dir = acc_diff / dist; // Unit vec

            let mut acc = acc_newton_inner(acc_dir, leaf.mass, dist, softening_factor_sq);

            if let Some(mond_fn) = mond {
                let x = acc.magnitude() / A0_MOND;
                acc = acc / mond_fn.μ(x);
            }

            Some(acc)
        })
        .reduce(Vec3::new_zero, |acc, elem| acc + elem)
}
