//! Acceleration calculations for n-body using Fast Multipole Methods (FMM) in conjunction with tree-(?)
//! Initial code from ChatGPT

// todo: Update this GPT basis

use lin_alg::f64::Vec3;

use crate::Body;

/// A node in the octree, storing either:
/// - A list of bodies (leaf node), OR
/// - 8 children (internal node).
///
/// We also store multipole expansion information here, for example:
/// - The total mass in the node
/// - The center of mass
/// - Higher-order moments (not detailed here)
#[derive(Debug)]
pub struct OctreeNode {
    pub center: Vec3,      // Center of the node in space
    pub half_width: f64,   // Half the width of this octree node
    pub mass: f64,         // Total mass (for M0 moment, or center-of-mass calc)
    pub com: Vec3,         // Center of mass (M1 / M0)
    pub bodies: Vec<Body>, // Bodies in this node (if leaf)
    pub children: Option<[Box<OctreeNode>; 8]>,
    // Higher-order multipole coefficients would go here in a full FMM implementation
}

impl OctreeNode {
    /// Create an empty octree node
    pub fn new(center: Vec3, half_width: f64) -> Self {
        Self {
            center,
            half_width,
            mass: 0.0,
            com: Vec3::new_zero(),
            bodies: Vec::new(),
            children: None,
        }
    }
}

/// Build an octree from a list of bodies.
/// This is a recursive function that subdivides space until
/// each node contains <= some number of bodies (leaf).
pub fn build_octree(
    bodies: &[Body],
    center: Vec3,
    half_width: f64,
    max_bodies: usize,
) -> OctreeNode {
    let mut node = OctreeNode::new(center, half_width);

    // Collect bodies that belong to this node
    let mut in_node = Vec::new();
    for b in bodies {
        // Check if b is in the bounding cube
        if within_bounds(&b.posit, &center, half_width) {
            in_node.push(b.clone());
        }
    }

    // If there are too many bodies, subdivide
    if in_node.len() > max_bodies {
        // Create the 8 children
        let mut children = Vec::with_capacity(8);
        for child_center in octant_centers(center, half_width) {
            let child_node = build_octree(&in_node, child_center, half_width / 2.0, max_bodies);
            children.push(Box::new(child_node));
        }
        node.children = Some(children.try_into().unwrap());
    } else {
        // Make this a leaf node
        node.bodies = in_node;
    }

    // Compute multipole expansions (we only do mass and COM here for brevity).
    let (mass, com) = compute_mass_com(&node);
    node.mass = mass;
    node.com = com;

    node
}

/// Compute total mass and center of mass for this node (and its children).
fn compute_mass_com(node: &OctreeNode) -> (f64, Vec3) {
    let mut total_mass = 0.0;
    let mut com = Vec3::new_zero();

    match &node.children {
        Some(children) => {
            for child in children {
                let (child_mass, child_com) = (child.mass, child.com);
                total_mass += child_mass;
                com += child_com * child_mass;
            }
        }
        None => {
            // Leaf node => sum directly
            for b in &node.bodies {
                total_mass += b.mass;
                com += b.posit * b.mass;
            }
        }
    }

    if total_mass > 0.0 {
        com *= 1.0 / total_mass;
    }
    (total_mass, com)
}

/// Helper: check if a position is within a node's bounding cube
fn within_bounds(pos: &Vec3, center: &Vec3, half_width: f64) -> bool {
    (pos.x >= center.x - half_width)
        && (pos.x < center.x + half_width)
        && (pos.y >= center.y - half_width)
        && (pos.y < center.y + half_width)
        && (pos.z >= center.z - half_width)
        && (pos.z < center.z + half_width)
}

/// Return the centers of each of the 8 octants.
fn octant_centers(center: Vec3, half_width: f64) -> Vec<Vec3> {
    let hw = half_width / 2.0;
    let mut centers = Vec::with_capacity(8);
    for dx in [-hw, hw].iter() {
        for dy in [-hw, hw].iter() {
            for dz in [-hw, hw].iter() {
                centers.push(Vec3::new(center.x + dx, center.y + dy, center.z + dz));
            }
        }
    }
    centers
}

/// Gravitational constant (for demonstration)
const G: f64 = 1.0;

/// A typical "theta" parameter used in FMM to decide whether
/// to use a node's multipole expansion or recurse.
const THETA: f64 = 0.5;

/// Top-level function that, for demonstration, returns the total gravitational
/// acceleration on the **first body** in `bodies` using an FMM approach.
///
/// In a real application, you would compute acceleration for **each** body,
/// or in parallel for all bodies. Here we just show how you might compute
/// it for one particular body.
pub fn acc_fmm(bodies: &[Body]) -> Vec3 {
    // Safety check
    if bodies.is_empty() {
        return Vec3::new_zero();
    }
    let target_body = &bodies[0];

    // 1. Build octree around the bounding region of bodies.
    //    For simplicity, pick a "center" and "half_width" that can encompass all bodies.
    //    In a real scenario, you might compute bounding boxes from the data.
    let center = Vec3::new_zero();
    let half_width = 100.0;
    let root = build_octree(bodies, center, half_width, 8);

    // 2. Walk the tree to accumulate acceleration
    compute_acceleration(&root, target_body, Vec3::new_zero())
}

/// Recursively computes the acceleration contribution from this node (and sub-nodes).
///
/// - If the node is "far enough" from the target body (based on theta),
///   use the node's multipole expansion (just mass & COM in this demo).
/// - Otherwise, descend (if the node is internal) or do direct summation (if leaf).
fn compute_acceleration(node: &OctreeNode, target: &Body, mut acc: Vec3) -> Vec3 {
    // If this node has no mass or is empty, return early
    if node.mass <= f64::EPSILON {
        return acc;
    }

    // Distance from target to the node center of mass
    let r_vec = node.com - target.posit;
    let dist_sq = r_vec.magnitude_squared();

    // If leaf or acceptable opening criterion => use multipole approximation
    if node.children.is_none() || use_multipole_approx(dist_sq, node.half_width) {
        // Use the node's mass & COM
        if dist_sq > 1e-12 {
            let dist = dist_sq.sqrt();
            let inv_dist3 = 1.0 / (dist_sq * dist);
            // Acceleration contribution: G * m_node / r^3 * r_vec
            let f = G * node.mass * inv_dist3;
            acc += r_vec * f;
        }
    } else {
        // Otherwise, descend if children exist
        if let Some(children) = &node.children {
            for child in children {
                acc = compute_acceleration(child, target, acc);
            }
        } else {
            // Node is a leaf with a small number of bodies => direct summation
            for b in &node.bodies {
                if std::ptr::eq(b, target) {
                    continue; // skip self
                }
                let r_vec = b.posit - target.posit;
                let dist_sq = r_vec.magnitude_squared();
                if dist_sq > 1e-12 {
                    let dist = dist_sq.sqrt();
                    let inv_dist3 = 1.0 / (dist_sq * dist);
                    let f = G * b.mass * inv_dist3;
                    acc = acc + r_vec * f;
                }
            }
        }
    }

    acc
}

/// Simple opening criterion for demonstration:
/// If (width of node / distance from node center to target) < THETA => use multipole
fn use_multipole_approx(dist_sq: f64, half_width: f64) -> bool {
    // approximate diameter of node
    let size = 2.0 * half_width;
    // compare ratio
    if dist_sq <= 1e-12 {
        // Avoid division by near zero => must do direct summation
        false
    } else {
        let dist = dist_sq.sqrt();
        (size / dist) < THETA
    }
}
