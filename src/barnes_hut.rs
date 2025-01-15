//! Code related to calculating acceleration between bodies using the Barnes-Hut
//! approach. It groups source masses. O(N Log(N)) complexity. Faster than naive
//! n-body; slower than Fast Multipole (FMM). Using in place of FMM due to its relative
//! ease of implementation.
//!
//! todo: If you can make this general enough, publish it as a standalone library.

// todo: Should we store Trees recursively, or with a top level struct that has a
// Vec of terminal nodes, and a vec of non-terminal ones?

use egui::IdMap;
use lin_alg::f64::Vec3;

use crate::{accel::acc_newton_inner, Body};

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
    Terminal(NodeTerminal),
}

#[derive(Debug)]
/// A terminal node, containing 0 or 1 body.
struct NodeTerminal {
    /// Inner: Body index.
    body: Option<usize>,
}

#[derive(Debug)]
/// A recursive tree. Terminates with `NodeType::NodeTerminal`.
pub struct Tree {
    data: Box<NodeType>,
    pub bounding_box: Cube, // todo temp pub?
    /// We use mass and center-of-mass to calculate Newtonian acceleration
    /// with an acted-on body.
    mass: f64,
    center_of_mass: Vec3,
}

impl Tree {
    /// Constructs a tree. Call this externaly using all bodies.
    pub fn new(bodies:  &[Body], id_acted_on: usize, z_offset: bool) -> Self {
        let bb = Cube::from_bodies(bodies, z_offset).unwrap();
        let body_ids: Vec<usize> = (0..bodies.len()).collect();

        // Convert &[Body] to &[&Body]
        let body_refs: Vec<&Body> = bodies.iter().collect();
        Self::new_internal(&body_refs, &body_ids, &bb, id_acted_on)
    }

    /// Constructs a tree. This can be called externally, but has a slightly lower-level API, requiring
    /// body IDs and a bounding box to be manually specified, for use during recursion.
    /// We assume that all bodies passed are inside the bounding box.
    /// `body_indices` must correspond to `bodies`.
    // pub fn new_internal(bodies: &[Body], body_ids: &[usize], bb: &Cube, id_acted_on: usize) -> Self {
    pub fn new_internal(bodies: &[&Body], body_ids: &[usize], bb: &Cube, id_acted_on: usize) -> Self {
        let data = match bodies.len() {
            0 => NodeType::Terminal(NodeTerminal { body: None }),
            1 => NodeType::Terminal(NodeTerminal {
                body: Some(body_ids[0]),
            }),
            _ => {
                let octants = bb.divide_into_octants();

                // Populate up to 8 children.
                let mut children = Vec::new();
                for octant in &octants {
                    let mut bodies_this_octant = Vec::new();
                    let mut body_indices_this_octant = Vec::new();

                    for (i, body) in bodies.iter().enumerate() {
                        if i == id_acted_on {
                            continue // Exclude the body acted on from the tree.
                        }
                        // Todo: Use a more efficient method, perhaps, where you position
                        // todo each body using > logic.
                        if octant.contains(body.posit) {
                            // bodies_this_octant.push(body.clone());
                            bodies_this_octant.push(*body);
                            body_indices_this_octant.push(body_ids[i]);
                        }
                    }
                    if !bodies_this_octant.is_empty() {
                        children.push(Box::new(Tree::new_internal(
                            bodies_this_octant.as_slice(),
                            &body_indices_this_octant,
                            octant,
                            id_acted_on
                        )))
                    }
                }

                NodeType::NonTerminal(children)
            }
        };

        // todo: Temp removed here.?
        // let (center_of_mass, mass) = center_of_mass(bodies);
        let (center_of_mass, mass) = (Vec3::new_zero(), 0.);

        Self {
            data: Box::new(data),
            bounding_box: bb.clone(),
            mass,
            center_of_mass,
        }
    }

    /// Get all bosides
    // fn get_bodies(&self) -> Vec<&Body> {
    fn get_bodies(&self) -> Vec<usize> {
        let mut result = Vec::new();
        match self.data.as_ref() {
            NodeType::NonTerminal(nodes) => {
                // Recur
                for node in nodes {
                    result.extend(node.get_bodies());
                }
            }
            NodeType::Terminal(node) => {
                // Terminate recursion.
                if let Some(body) = &node.body {
                    result.push(*body);
                }
            }
        }

        result
    }

    /// For debugging only?
    /// Get all the terminal nodes (i.e. that contain 0 or 1 bodies, and no sub-nodes containing bodies)
    pub fn get_leaves(&self) -> Vec<&Self> {
        let mut result = Vec::new();

        match self.data.as_ref() {
            NodeType::NonTerminal(nodes) => {
                // Recur
                for node in nodes {
                    result.extend(node.get_leaves());
                }
            }
            NodeType::Terminal(node) => {
                // Terminate recursion.
                result.push(self);
            }
        }

        result
    }

    // /// Calculates the position, and value of the center-of-mass. Run this once all sub-nodes of this are populated.
    // fn update_mass(&mut self) {
    //     (self.center_of_mass, self.mass) = center_of_mass(&self.get_bodies());
    // }

    fn is_terminal(&self) -> bool {
        match self.data.as_ref() {
            NodeType::Terminal(_) => true,
            _ => false,
        }
    }

    /// Gets the next node that should count towards the force calculation for the current particle.
    ///
    /// Whether a node is or isn't sufficiently far away from a body,
    /// depends on the quotient s/d,
    /// where s is the width of the region represented by the internal node,
    /// and d is the distance between the body and the node's center of mass.
    /// The node is sufficiently far away when this ratio is smaller than a threshold value θ.
    /// The parameter θ determines the accuracy of the simulation;
    /// larger values of θ increase the speed of the simulation but decreases its accuracy.
    /// If θ = 0, no internal node is treated as a single body and the algorithm degenerates to a direct-sum algorithm.
    fn collect_relevant_nodes<'a>(&'a self, posit_acted_on: Vec3, θ: f64, results: &mut Vec<&'a Tree>) {
        let dist = (posit_acted_on - self.center_of_mass).magnitude();
        let s = self.bounding_box.width;

        if s / dist < θ {
            // If the node is far enough, add it to the results.
            results.push(self);
            return;
        }

        // Otherwise, traverse into child nodes if this is a non-terminal node.
        if let NodeType::NonTerminal(ref children) = *self.data {
            for child in children {
                child.collect_relevant_nodes(posit_acted_on, θ, results);
            }
        }
    }
}

#[derive(Clone, Copy, PartialEq)]
enum NodeTypeAlt {
    NonTerminal,
    Terminal,
}

pub struct LeafAlt {
    bb: Cube,
    children_nonterminal: Vec<usize>,
    children_terminal: Vec<usize>,
}

pub struct NodeTerminalAlt {
    bb: Cube,
    body: Option<usize>, // Index
}

/// An alternate approach, where we store nodes in a flat manner, and index them.
pub struct TreeAlt {
    // limits: GridLimits,
    // bodies: Vec<Body>,
    bodies: Vec<usize>, // Body index
    /// All non-terminal nodes.
    leaves: Vec<LeafAlt>,
    /// All terminal nodes. None means no body.
    nodes_terminal: Vec<NodeTerminalAlt>, // Body indices
}

impl TreeAlt {
    pub fn new(bodies: &Vec<Body>) -> Self {
        Self {
            bodies: Vec::new(),         // todo temp
            leaves: Vec::new(),         // todo temp
            nodes_terminal: Vec::new(), // todo temp
        }
    }
}

/// Compute center of mass as a position, and mass value.
// fn center_of_mass(bodies: &[&Body]) -> (Vec3, f64) {
fn center_of_mass(bodies: &[Body]) -> (Vec3, f64) {
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
pub fn calc_acc_bh(
    posit_acted_on: Vec3,
    id_acted_on: usize,
    bodies_other: &mut [Body],
    θ: f64,
    softening_factor_sq: f64,
) -> Vec3 {
    let tree = Tree::new(bodies_other, id_acted_on, true);
    let mut relevant_nodes = Vec::new();

    tree.collect_relevant_nodes(posit_acted_on, θ, &mut relevant_nodes);

    let mut result = Vec3::new_zero();

    // Compute acceleration contributions from each relevant node.
    for node in relevant_nodes {
        let acc_diff = node.center_of_mass - posit_acted_on;
        let dist = acc_diff.magnitude();
        let acc_dir = acc_diff / dist; // Unit vec

        // Skip self-interaction.
        if dist < f64::EPSILON {
            continue;
        }

        let acc_contribution = acc_newton_inner(acc_dir, node.mass, dist,softening_factor_sq);
        result += acc_contribution;
    }


    result
}
