//! Code related to calculating acceleration between bodies using the Barnes-Hut
//! approach. It groups source masses. O(N Log(N)) complexity. Faster than naive
//! n-body; slower than Fast Multipole (FMM). Using in place of FMM due to its relative
//! ease of implementation.
//!
//! todo: If you can make this general enough, publish it as a standalone library.

// todo: Should we store Trees recursively, or with a top level struct that has a
// Vec of terminal nodes, and a vec of non-terminal ones?


// todo: You can use rayon more throughout this, e.g. during tree construction.


use lin_alg::f64::Vec3;

use rayon::prelude::*;

use crate::{accel::acc_newton_inner, Body};
use crate::accel::MondFn;
use crate::units::A0_MOND;

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
    /// Construct minimum limits that encompass all bodies.
    /// The z offset is intended for the case where the Z coordinate for all particles is 0.
    /// This prevents the divisions straddling the points, doubling the number of nodes.
    pub fn from_bodies(bodies: &[Body], z_offset: bool) -> Option<Self> {
        // todo: You could also accept a pad, and run this
        // periodically instead of eacy time.
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

        if z_offset {
            z_max += 1e-5;
        }

        let x_size = x_max - x_min;
        let y_size = y_max - y_min;
        let z_size = z_max - z_min;

        // Coerce to a cube.
        let mut width = x_size.max(y_size).max(z_size);
        let width_div2 = width / 2.;

        let mut center = Vec3::new(
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
    pub fn divide_into_octants(&self) -> [Self; 8] {
        let width = self.width / 2.;
        let wd2 = self.width / 4.; // short for brevity below.

        // Every combination of + and - for the center offset.
        [
            Self::new(self.center + Vec3::new(wd2, wd2, wd2), width, wd2),
            Self::new(self.center + Vec3::new(-wd2, wd2, wd2), width, wd2),
            Self::new(self.center + Vec3::new(wd2, -wd2, wd2), width, wd2),
            Self::new(self.center + Vec3::new(wd2, wd2, -wd2), width, wd2),
            Self::new(self.center + Vec3::new(-wd2, -wd2, wd2), width, wd2),
            Self::new(self.center + Vec3::new(wd2, -wd2, -wd2), width, wd2),
            Self::new(self.center + Vec3::new(-wd2, wd2, -wd2), width, wd2),
            Self::new(self.center + Vec3::new(-wd2, -wd2, -wd2), width, wd2),
        ]
    }

    // todo: Use or rem A/R.
    pub fn contains(&self, posit: Vec3) -> bool {
        self.x_min <= posit.x
            && posit.x <= self.x_max
            && self.y_min <= posit.y
            && posit.y <= self.y_max
            && self.z_min <= posit.z
            && posit.z <= self.z_max
    }
}

#[derive(Debug)]
enum NodeType {
    NonTerminal(Vec<Box<Tree>>),
    Terminal
}

#[derive(Debug)]
/// A recursive tree. Each node can be subdivided  Terminates with `NodeType::NodeTerminal`.
pub struct Tree {
    // todo: Consider replacing this with `children: Vec<Box<Tree>>`. Terminal means children is empty.
    // data: Box<NodeType>,
    children: Vec<Box<Tree>>,
    pub bounding_box: Cube, // todo temp pub?
    /// We use mass and center-of-mass to calculate Newtonian acceleration
    /// with an acted-on body.
    mass: f64,
    center_of_mass: Vec3,
}

// todo: Consider converting this recursive struct to an iterative approach
// todo: with an explicit stack to reduce function call overhead.
impl Tree {
    /// Constructs a tree. Call this externaly using all bodies.
    pub fn new(
        bodies: &[Body],
        posit_acted_on: Vec3,
        id_acted_on: usize,
        θ: f64,
        z_offset: bool,
    ) -> Self {
        let bb = Cube::from_bodies(bodies, z_offset).unwrap();

        // Convert &[Body] to &[&Body], and remove the acted-on body.
        let body_refs: Vec<&Body> = bodies.iter().enumerate().filter(|(i, b)| i != *id_acted_on).collect();

        let body_ids: Vec<usize> = (0..body_refs.len()).collect();
        Self::new_internal(&body_refs, &body_ids, &bb, posit_acted_on, θ)
    }

    /// Constructs a tree. This can be called externally, but has a slightly lower-level API, requiring
    /// body IDs and a bounding box to be manually specified, for use during recursion.
    /// We assume that all bodies passed are inside the bounding box.
    /// `body_indices` must correspond to `bodies`.
    pub fn new_internal(
        bodies: &[&Body],
        body_ids: &[usize],
        bb: &Cube,
        posit_acted_on: Vec3,
        // id_acted_on: usize,
        θ: f64,
    ) -> Self {
        let (center_of_mass, mass) = center_of_mass(bodies);

        let data = match bodies.len() {
            // 0 | 1 => NodeType::Terminal,
            0 | 1 => Vec::new(),
            _ => {
                let dist = (posit_acted_on - center_of_mass).magnitude();

                // The actor is far; group all bodies in this bounding box.
                if bb.width / dist < θ {
                    // NodeType::Terminal
                    Vec::new()
                } else {
                    // The actor is close; continue subdividing.
                    let octants = bb.divide_into_octants();

                    // Populate up to 8 children.
                    let mut children = Vec::new();
                    for octant in &octants {
                        let mut bodies_this_octant = Vec::new();
                        let mut body_indices_this_octant = Vec::new();

                        for (i, body) in bodies.iter().enumerate() {
                            // Todo: Use a more efficient method, perhaps, where you position
                            // todo each body using > logic.
                            if octant.contains(body.posit) {
                                bodies_this_octant.push(*body);
                                body_indices_this_octant.push(body_ids[i]);
                            }
                        }
                        if !bodies_this_octant.is_empty() {
                            children.push(Box::new(Tree::new_internal(
                                bodies_this_octant.as_slice(),
                                &body_indices_this_octant,
                                octant,
                                posit_acted_on,
                                // id_acted_on,
                                θ,
                            )))
                        }
                    }

                    // NodeType::NonTerminal(children)
                    children
                }
            }
        };

        Self {
            children: data,
            bounding_box: bb.clone(), // todo: Re
            mass,
            center_of_mass,
        }
    }

    /// For debugging only?
    /// Get all the terminal nodes (i.e. that contain 0 or 1 bodies, and no sub-nodes containing bodies)
    pub fn get_leaves(&self) -> Vec<&Self> {
        let mut result = Vec::new();


        // match self.data.as_ref() {
        if !self.children.is_empty() {
            // NodeType::NonTerminal(nodes) => {
            // Recur
            // for node in nodes {
            for node in &self.children {
                result.extend(node.get_leaves());
            }
            // }
            // NodeType::Terminal(_) => {
            // NodeType::Terminal => {
        } else {
            // Terminate recursion.
            result.push(self);
        }

        result
    }
    //
    // fn is_terminal(&self) -> bool {
    //     self.children.as_ref().is_empty()
    //     // match self.data.as_ref() {
    //     //     // NodeType::Terminal(_) => true,
    //     //     NodeType::Terminal => true,
    //     //     _ => false,
    //     // }
    // }
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

/// Calculate Newtonian acceleration using the Barnes Hut algorithm.
pub fn acc_newton_bh(
    posit_acted_on: Vec3,
    id_acted_on: usize,
    bodies_other: &[Body],
    mond: Option<MondFn>,
    θ: f64,
    softening_factor_sq: f64,
) -> Vec3 {
    // todo: The tree building is taking too long.
    let tree = Tree::new(bodies_other, posit_acted_on, id_acted_on, θ, true);

    // Self-interaction is prevented. in the tree construction.
    tree.get_leaves()
        // .par_iter()
        .iter()
        .map(|leaf| {
            let acc_diff = leaf.center_of_mass - posit_acted_on;
            let dist = acc_diff.magnitude();
            let acc_dir = acc_diff / dist; // Unit vec

            let mut acc = acc_newton_inner(acc_dir, leaf.mass, dist, softening_factor_sq);

            if let Some(mond_fn) = mond {
                let x = acc.magnitude() / A0_MOND;
                acc = acc / mond_fn.μ(x);
            }

            acc
        })
        // .reduce(Vec3::new_zero, |acc, elem| acc + elem)
        .fold(Vec3::new_zero(), |acc, elem| acc + elem)
}
