use std::sync::Arc;
use lin_alg::f64::Vec3;
use crate::barnes_hut::center_of_mass;
use crate::Body;


/// Barnes–Hut node that can either be an internal node (with children)
/// or a leaf that directly stores a small number of bodies.
pub enum Node {
    Leaf(LeafNode),
    Internal(InternalNode),
}

/// A leaf node storing a list of bodies.
pub struct LeafNode {
    pub bodies: Vec<Arc<Body>>,
    pub center_of_mass: Vec3,
    pub total_mass: f64,
    /// Center of this node's region in space
    pub region_center: Vec3,
    /// Half the size of this node's cubic region
    pub region_half_size: f64,
}

/// An internal node storing aggregated information and octant children.
pub struct InternalNode {
    pub center_of_mass: Vec3,
    pub total_mass: f64,
    pub region_center: Vec3,
    pub region_half_size: f64,
    pub children: Vec<Option<Box<Node>>>, // 8 children for an octree
}

// Threshold for splitting a leaf into an internal node
const LEAF_CAPACITY: usize = 4;

/// Builds a Barnes–Hut tree from a list of bodies and a known bounding cube.
///
/// * `bodies` - the slice of bodies to include in this node
/// * `region_center` - the center of the cube region
/// * `region_half_size` - half the length of a side of the cube region
pub fn build_tree(
    bodies: &[Arc<Body>],
    region_center: Vec3,
    region_half_size: f64,
) -> Node {
    if bodies.len() <= LEAF_CAPACITY {
        // Create a leaf node directly
        let (com, mass) = center_of_mass(bodies);
        Node::Leaf(LeafNode {
            bodies: bodies.to_vec(),
            center_of_mass: com,
            total_mass: mass,
            region_center,
            region_half_size,
        })
    } else {
        // Subdivide into 8 child octants
        let mut children_bodies: [Vec<Arc<Body>>; 8] = Default::default();
        for body in bodies {
            let octant_idx = find_octant_index(body.posit, region_center);
            children_bodies[octant_idx].push(body.clone());
        }

        let half = region_half_size * 0.5;
        let mut children: Vec<Option<Box<Node>>> = vec![None; 8];
        for (i, child_bodies) in children_bodies.into_iter().enumerate() {
            if !child_bodies.is_empty() {
                let child_center = compute_child_center(region_center, region_half_size, i);
                let child_node = build_tree(&child_bodies, child_center, half);
                children[i] = Some(Box::new(child_node));
            }
        }

        // Compute center of mass from children
        let mut total_mass = 0.0;
        let mut weighted_pos = Vec3 { x: 0.0, y: 0.0, z: 0.0 };
        for child in children.iter() {
            if let Some(node) = child {
                match node.as_ref() {
                    Node::Leaf(leaf) => {
                        total_mass += leaf.total_mass;
                        weighted_pos = weighted_pos + leaf.center_of_mass * leaf.total_mass;
                    }
                    Node::Internal(internal) => {
                        total_mass += internal.total_mass;
                        weighted_pos =
                            weighted_pos + internal.center_of_mass * internal.total_mass;
                    }
                }
            }
        }
        let center_of_mass = if total_mass > 0.0 {
            weighted_pos * (1.0 / total_mass)
        } else {
            region_center
        };

        Node::Internal(InternalNode {
            center_of_mass,
            total_mass,
            region_center,
            region_half_size,
            children,
        })
    }
}


/// Recursively compute the gravitational force on `body` by traversing the tree.
fn compute_force_single(body: &Body, node: &Node, θ: f64) -> Vec3 {
    match node {
        Node::Leaf(leaf) => {
            // Sum over each body in this leaf
            let mut force = Vec3 { x: 0.0, y: 0.0, z: 0.0 };
            for other in leaf.bodies.iter() {
                // Skip self-interaction if the pointers match
                if std::ptr::eq(&**other as *const Body, body as *const Body) {
                    continue;
                }
                force = force + pairwise_force(body, other);
            }
            force
        }
        Node::Internal(internal) => {
            // Check if node is "far enough"
            let dx = internal.center_of_mass - body.posit;
            let dist_sq = dx.magnitude_squared();
            let dist = dist_sq.sqrt();

            // If region size / distance < THETA, treat as single mass
            if (internal.region_half_size * 2.0) / dist < θ {
                // Compute force from center-of-mass
                pairwise_force_cm(body, internal.center_of_mass, internal.total_mass)
            } else {
                // Otherwise, descend into children
                let mut force = Vec3 { x: 0.0, y: 0.0, z: 0.0 };
                for child_opt in internal.children.iter() {
                    if let Some(child) = child_opt {
                        force = force + compute_force_single(body, child);
                    }
                }
                force
            }
        }
    }
}


/// Determine in which octant of the parent region a position should go.
/// We'll number octants as bits in x/y/z. E.g.:
/// 0 => -x, -y, -z
/// 1 => +x, -y, -z
/// 2 => -x, +y, -z
/// 3 => +x, +y, -z
/// 4 => -x, -y, +z
/// 5 => +x, -y, +z
/// 6 => -x, +y, +z
/// 7 => +x, +y, +z
fn find_octant_index(pos: Vec3, center: Vec3) -> usize {
    let mut idx = 0;
    if pos.x >= center.x {
        idx |= 1;
    }
    if pos.y >= center.y {
        idx |= 2;
    }
    if pos.z >= center.z {
        idx |= 4;
    }
    idx
}

/// Compute the center of one of the eight children given parent region info.
fn compute_child_center(parent_center: Vec3, parent_half_size: f64, octant_idx: usize) -> Vec3 {
    let half = parent_half_size * 0.5;
    let mut x = parent_center.x;
    let mut y = parent_center.y;
    let mut z = parent_center.z;

    if (octant_idx & 1) != 0 {
        x += half;
    } else {
        x -= half;
    }
    if (octant_idx & 2) != 0 {
        y += half;
    } else {
        y -= half;
    }
    if (octant_idx & 4) != 0 {
        z += half;
    } else {
        z -= half;
    }

    Vec3 { x, y, z }
}