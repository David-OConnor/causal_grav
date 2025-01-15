//! Code related to calculating acceleration between bodies using the Barnes-Hut
//! approach. It groups source masses. O(N Log(N)) complexity. Faster than naive
//! n-body; slower than Fast Multipole (FMM). Using in place of FMM due to its relative
//! ease of implementation.
//!
//! todo: If you can make this general enough, publish it as a standalone library.

// todo: Should we store Trees recursively, or with a top level struct that has a
// Vec of terminal nodes, and a vec of non-terminal ones?

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
    pub fn from_bodies(bodies: &[Body]) -> Option<Self> {
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

        let x_size = x_max - x_min;
        let y_size = y_max - y_min;
        let z_size = z_max - z_min;

        // Coerce to a cube.
        let mut width = x_size.max(y_size).max(z_size);
        let width_div2 = width / 2.;

        // todo: QC this...
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

    // /// Uses binary logic.
    // /// Note: The centers and sizes can be computed dyanmically; we use parameters as a cahe.
    // fn from_octant(octant: usize, center_x: f64, center_y: f64, center_z: f64,
    //                half_size_x: f64, half_size_y: f64, half_size_z: f64) -> Self {
    //     assert!(octant < 8);
    //
    //     // todo: This isn't right. On the right track, but wrong.
    //     Self {
    //         x_min: center_x + ((octant >> 0) & 1) as f64,
    //         x_max: center_x + ((octant >> 0) & 1) as f64,
    //         y_min: center_y + ((octant >> 1) & 1) as f64,
    //         y_max: center_y + ((octant >> 1) & 1) as f64,
    //         z_min: center_z + ((octant >> 2) & 1) as f64,
    //         z_max: center_z + ((octant >> 2) & 1) as f64,
    //     }
    // }

    /// Divide this into equal-area octants.
    pub fn divide_into_octants(&self) -> [Self; 8] {
        let center_x = (self.x_max + self.x_min) / 2.;
        let center_y = (self.y_max + self.y_min) / 2.;
        let center_z = (self.z_max + self.z_min) / 2.;

        // let half_size_x = (self.x_max - self.x_min) / 2.;
        // let half_size_y = (self.y_max - self.y_min) / 2.;
        // let half_size_z = (self.z_max - self.z_min) / 2.;
        //
        // [0, 1, 2, 3, 5, 6, 7, 8].map(|i| self.from_octant(i, center_x, center_y, center_z, half_size_x, half_size_y, half_size_z))

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
            // Each combination of (min, center) and (center, max).
            // Self {
            //     width,
            //     x_min: self.x_min,
            //     x_max: self.center.x,
            //     y_min: self.y_min,
            //     y_max: self.center.y,
            //     z_min: self.z_min,
            //     z_max: self.center.z,
            // },
            // Self {
            //     width,
            //     x_min: self.center.x,
            //     x_max: self.x_max,
            //     y_min: self.y_min,
            //     y_max: self.center.y,
            //     z_min: self.z_min,
            //     z_max: self.center.z,
            // },
            // Self {
            //     width,
            //     x_min: self.x_min,
            //     x_max: self.center.x,
            //     y_min: self.center.y,
            //     y_max: self.y_max,
            //     z_min: self.z_min,
            //     z_max: self.center.z,
            // },
            // Self {
            //     width,
            //     x_min: self.center.x,
            //     x_max: self.x_max,
            //     y_min: self.center.y,
            //     y_max: self.y_max,
            //     z_min: self.z_min,
            //     z_max: self.center.z,
            // },
            // Self {
            //     width,
            //     x_min: self.x_min,
            //     x_max: self.center.x,
            //     y_min: self.y_min,
            //     y_max: self.center.y,
            //     z_min: self.center.z,
            //     z_max: self.z_min,
            // },
            // Self {
            //     width,
            //     x_min: self.center.x,
            //     x_max: self.x_max,
            //     y_min: self.y_min,
            //     y_max: self.center.y,
            //     z_min: self.center.z,
            //     z_max: self.z_min,
            // },
            // Self {
            //     width,
            //     x_min: self.x_min,
            //     x_max: self.center.x,
            //     y_min: self.center.y,
            //     y_max: self.y_max,
            //     z_min: self.center.z,
            //     z_max: self.z_min,
            // },
            // Self {
            //     width,
            //     x_min: self.center.x,
            //     x_max: self.x_max,
            //     y_min: self.center.y,
            //     y_max: self.y_max,
            //     z_min: self.center.z,
            //     z_max: self.z_min,
            // },
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

    // // todo: Use or rem A/R.
    // pub fn intersects(&self, other: &Self) -> bool {
    //     false // todo: Implement A/R.
    //
    // }
}

#[derive(Debug)]
enum NodeType {
    // NonTerminal([Box<Tree>; 8]),
    NonTerminal(Vec<Box<Tree>>),
    Terminal(NodeTerminal),
}

#[derive(Debug)]
/// A terminal node, containing 0 or 1 body.
struct NodeTerminal {
    /// Inner: Body index.
    // pub body: Option<Body>,
    body: Option<usize>,
}

#[derive(Debug)]
/// A recursive tree. Terminates with `NodeType::NodeTerminal`.
pub struct Tree {
    data: Box<NodeType>,
    pub bounding_box: Cube, // todo temp pub?
    /// We cache mass ands its center. We use this to calculate Newtonian acceleration
    /// with a destination body.
    mass: f64, // Calculated as needed.
    center_of_mass: Vec3,   // Calculated as needed.
}

impl Tree {
    /// Constructs a tree. Call this externaly using all bodies, and a bounding box that
    /// encompasses them.
    /// We assume that all bodies passed are inside the bounding box.
    /// `body_indices` must correspond to `bodies`.
    pub fn new(bodies: &[Body], body_indices: &[usize], bb: &Cube) -> Self {
        let data = match bodies.len() {
            0 => NodeType::Terminal(NodeTerminal { body: None }),
            1 => NodeType::Terminal(NodeTerminal {
                body: Some(body_indices[0]),
            }),
            _ => {
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
                            bodies_this_octant.push(body.clone()); // todo: Clone...
                            body_indices_this_octant.push(body_indices[i]);
                        }
                    }
                    if !bodies_this_octant.is_empty() {
                        children.push(Box::new(Tree::new(
                            &bodies_this_octant,
                            &body_indices_this_octant,
                            octant,
                        )))
                    }
                }

                NodeType::NonTerminal(children)
            }
        };

        let (center_of_mass, mass) = center_of_mass(bodies);

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
    pub fn get_all_children(&self) -> Vec<&Self> {
        let mut result = Vec::new();

        match self.data.as_ref() {
            NodeType::NonTerminal(nodes) => {
                // Recur
                for node in nodes {
                    result.extend(node.get_all_children());
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
    fn find_relevant_node(&self, posit_acted_on: Vec3, θ: f64) -> &Self {
        let dist = (posit_acted_on - self.center_of_mass).magnitude();
        let s = self.bounding_box.width;

        if s / dist < θ {
            return self;
        }

        // // If region size / distance < THETA, treat as single mass
        // if (internal.region_half_size * 2.0) / dist < θ {
        //     // Compute force from center-of-mass
        //     pairwise_force_cm(body, internal.center_of_mass, internal.total_mass)
        // } else {
        //     // Otherwise, descend into children
        //     let mut force = Vec3 { x: 0.0, y: 0.0, z: 0.0 };
        //     for child_opt in internal.children.iter() {
        //         if let Some(child) = child_opt {
        //             force = force + compute_force_single(body, child);
        //         }
        //     }
        //     force
        // }

        // // If not far enough away, add children to the stack.
        // for (quadrant, child) in node.children.iter().enumerate() {
        //     match child {
        //         Some(child) => self.stack.push((child, bb.child(quadrant))),
        //         None => (),
        //     }
        // }

        self
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
    bodies_other: &mut [Body],
    id_acted_on: usize,
    θ: f64,
) -> Vec3 {
    let bb = Cube::from_bodies(bodies_other).unwrap();
    let body_ids: Vec<usize> = (0..bodies_other.len()).collect();

    // todo: You must address self-interaction; likely during tree construction.
    // todo: Check examples.
    let tree = Tree::new(bodies_other, &body_ids, &bb);

    let mut result = Vec3::new_zero();

    // Transverse our tree, choosing a grouping coarseness based on distance.
    // let dist = (posit_acted_on - self.center_of_mass).magnitude();
    // let s = self.bounding_box.width;
    //
    // if s / dist < θ {
    //     return self
    // }

    // Iterate over bodies acting on our target.
    for (i, body) in bodies_other.iter().enumerate() {
        if i == id_acted_on {
            continue; // self-interaction
        }

        let node_this_body = tree.find_relevant_node(posit_acted_on, θ);

        // todo: How do we know if a body has already been taken into account?
        // todo: mark the tree node as used??

        // match tree.data.as_ref() {
        //     NodeType::NonTerminal(nodes) => {
        //         // Recur
        //         for node in nodes {
        //             result.extend(node.get_bodies());
        //         }
        //     }
        //     NodeType::Terminal(node) => {
        //         // Terminate recursion.
        //         if let Some(body) = &node.body {
        //             result.push(body);
        //         }
        //     }
        // }
        //
        // let v = if tree.bounding_box.width / dist > θ {
        //     let mass = tree.mass;
        //     let posit = tree.center_of_mass;
        // };

        // for group in groups {
        //     let group.posit - body.posit; // todo: QC order.
        //     let acc_dir = acc_diff.magnitude();
        //     body.acc += acc_newton_inner(acc_dir, group.mass);
        // }
    }
    result
}
