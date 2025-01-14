//! Code related to calculating acceleration between bodies using the Barnes-Hut
//! approach. It groups source masses. O(N Log(N)) complexity. Faster than naive
//! n-body; slower than Fast Multipole (FMM). Using in place of FMM due to its relative
//! ease of implementation.

// todo: Should we store Trees recursively, or with a top level struct that has a
// Vec of terminal nodes, and a vec of non-terminal ones?

use lin_alg::f64::Vec3;

use crate::{accel::acc_newton_inner, Body};

#[derive(Clone)]
pub struct BoundingBox {
    // pub is Temp
    pub x_min: f64,
    pub x_max: f64,
    pub y_min: f64,
    pub y_max: f64,
    pub z_min: f64,
    pub z_max: f64,
}

impl BoundingBox {
    /// Construct minimum limits that encompass all bodies.
    pub fn new(bodies: &[Body]) -> Option<Self> {
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

        Some(Self {
            x_min,
            x_max,
            y_min,
            y_max,
            z_min,
            z_max,
        })
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

        // Each combination of (min, center) and (center, max).
        [
            Self {
                x_min: self.x_min,
                x_max: center_x,
                y_min: self.y_min,
                y_max: center_y,
                z_min: self.z_min,
                z_max: center_z,
            },
            Self {
                x_min: center_x,
                x_max: self.x_max,
                y_min: self.y_min,
                y_max: center_y,
                z_min: self.z_min,
                z_max: center_z,
            },
            Self {
                x_min: self.x_min,
                x_max: center_x,
                y_min: center_y,
                y_max: self.y_max,
                z_min: self.z_min,
                z_max: center_z,
            },
            Self {
                x_min: center_x,
                x_max: self.x_max,
                y_min: center_y,
                y_max: self.y_max,
                z_min: self.z_min,
                z_max: center_z,
            },
            Self {
                x_min: self.x_min,
                x_max: center_x,
                y_min: self.y_min,
                y_max: center_y,
                z_min: center_z,
                z_max: self.z_min,
            },
            Self {
                x_min: center_x,
                x_max: self.x_max,
                y_min: self.y_min,
                y_max: center_y,
                z_min: center_z,
                z_max: self.z_min,
            },
            Self {
                x_min: self.x_min,
                x_max: center_x,
                y_min: center_y,
                y_max: self.y_max,
                z_min: center_z,
                z_max: self.z_min,
            },
            Self {
                x_min: center_x,
                x_max: self.x_max,
                y_min: center_y,
                y_max: self.y_max,
                z_min: center_z,
                z_max: self.z_min,
            },
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

enum NodeType {
    Leaf([Box<Tree>; 8]),
    Terminal(NodeTerminal),
}

#[derive(Default)]
/// A terminal node, containing 0 or 1 body.
struct NodeTerminal {
    pub body: Option<Body>,
    // pub bb: BoundingBox,
}

// todo: Leaf for the name?
struct Tree {
    pub node_type: Box<NodeType>,
    pub bb: BoundingBox,
    /// We cache mass ands its center. We use this to calculate Newtonian acceleration
    /// with a destination body.
    mass: f64, // Calculated
    center_of_mass: Vec3, // Calculated
}

impl Tree {
    pub fn new(bodies: &[Body], bb: &BoundingBox) -> Self {
        let node_type = match bodies.len() {
            0 => NodeType::Terminal(NodeTerminal { body: None }),
            1 => NodeType::Terminal(NodeTerminal {
                body: Some(bodies[0].clone()),
            }), // todo: Clone...
            _ => {
                let octants = bb.divide_into_octants();

                let mut bodies_by_octant = Vec::new();
                for octant in &octants {
                    let mut bodies_this_octant = Vec::new();
                    for body in bodies {
                        if octant.contains(body.posit) {
                            bodies_this_octant.push(body.clone()); // todo: Clone...
                        }
                    }
                    bodies_by_octant.push(bodies_this_octant);
                }

                NodeType::Leaf([
                    Box::new(Tree::new(&bodies_by_octant[0], &octants[0])),
                    Box::new(Tree::new(&bodies_by_octant[1], &octants[1])),
                    Box::new(Tree::new(&bodies_by_octant[2], &octants[2])),
                    Box::new(Tree::new(&bodies_by_octant[3], &octants[3])),
                    Box::new(Tree::new(&bodies_by_octant[4], &octants[4])),
                    Box::new(Tree::new(&bodies_by_octant[5], &octants[5])),
                    Box::new(Tree::new(&bodies_by_octant[6], &octants[6])),
                    Box::new(Tree::new(&bodies_by_octant[7], &octants[7])),
                ])
            }
        };

        Self {
            node_type: Box::new(node_type),
            bb: bb.clone(),                   // todo: Avoid clone?
            mass: 0.,                         // Calculated later,
            center_of_mass: Vec3::new_zero(), // calculated later.
        }
    }

    fn get_bodies(&self) -> Vec<&Body> {
        let mut result = Vec::new();
        match self.node_type.as_ref() {
            NodeType::Leaf(nodes) => {
                // Recur
                for node in nodes {
                    result.extend(node.get_bodies());
                }
            }
            NodeType::Terminal(node) => {
                // Terminate recursion.
                if let Some(body) = &node.body {
                    result.push(body);
                }
            }
        }

        result
    }

    /// Calculates the position, and value of the center-of-mass. Run this once all sub-nodes of this are populated.
    fn calc_mass(&mut self) {
        (self.center_of_mass, self.mass) = center_of_mass(&self.get_bodies());
    }
}

#[derive(Clone, Copy, PartialEq)]
enum NodeTypeAlt {
    Leaf,
    Terminal,
}

pub struct LeafAlt {
    bb: BoundingBox,
    children_leaves: Vec<usize>,
    children_terminal: Vec<usize>,
}

pub struct NodeTerminalAlt {
    bb: BoundingBox,
    body: Option<usize>, // Index
}

/// An alternate approach, where we store nodes in a flat manner, and index them.
pub struct TreeAlt {
    // limits: GridLimits,
    bodies: Vec<Body>,
    /// All non-terminal nodes.
    leaves: Vec<LeafAlt>,
    /// All terminal nodes. None means no body.
    nodes_terminal: Vec<NodeTerminalAlt>, // Body indices
}

impl TreeAlt {
    pub fn new(bodies: &Vec<Body>) -> Self {
        Self {
            bodies: bodies.clone(),     // todo: Way to avoid this, and instead use refs?
            leaves: Vec::new(),         // todo temp
            nodes_terminal: Vec::new(), // todo temp
        }
    }
}

/// Compute center of mass as a position, and mass value.
pub fn center_of_mass(bodies: &[&Body]) -> (Vec3, f64) {
    // todo pub temp
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
