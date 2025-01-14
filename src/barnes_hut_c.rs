//! https://github.com/Katsutoshii/barnes-hut-rs/tree/master/src


use crate::Body;
use lin_alg::f64::Vec3;
use crate::barnes_hut::BoundingBox;

fn in_bounds(v: Vec3, min: Vec3, max: Vec3) -> bool {
    v.x >= min.x && v.x <= max.x &&
        v.y >= min.y && v.y <= max.y &&
        v.z >= min.z && v.z <= max.z
}


const EPSILON: f64 = 1e-4;

/// Computes the l2 norm of a 2d vector represented by x1, y1, x2, y2
fn l2(x1: f64, y1: f64, x2: f64, y2: f64) -> f64 {
    let dx: f64 = x2 - x1;
    let dy: f64 = y2 - y1;
    (dx * dx + dy * dy).sqrt()
}

/// Definition of the mass quadtree
#[derive(Debug)]
pub struct MassOcttree {
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub m: f64,
    pub children: Vec<Option<Self>>,
}

/// Implementation for the mass quadtree
impl MassOcttree {
    /// Constructs a child with no children
    pub fn empty() -> Self {
        Self {
            x: 0.,
            y: 0.,
            z: 0.,
            m: 0.,
            children: vec![None, None, None, None, None, None, None, None]
        }
    }

    // Constructs a new child under a node
    pub fn new_child(&mut self, quadrant: usize, x: f64, y: f64, z: f64, m: f64) {
        // println!("New child ({}, {}, {}) under ({}, {}, {}) in quad {}", x, y, m, self.x, self.y, self.m, quadrant);
        self.children[quadrant] = Some(Self {
            x,
            y,
            z,
            m,
            children: vec![None, None, None, None, None, None, None, None]
        })
    }

    /// Constructs a quadtree for the given bounds and list of points
    pub fn new(r: &Vec<Vec3>, m: &Vec<f64>, bb: BoundingBox) -> Self {
        let mut root = Self::empty();
        for i in 0..r.len() {
            root.insert(r[i].x, r[i].y, m[i], bb);
        }
        root
    }

    // Updates the center of mass
    pub fn update_com(&mut self, x: f64, y: f64, m: f64) {
        let total_m: f64 = self.m + m;
        self.x = (self.m * self.x + m * x) / total_m;
        self.y = (self.m * self.y + m * y) / total_m;
        self.m = total_m;
    }

    /// Inserts a point into the quadtree.
    pub fn insert(&mut self, x: f64, y: f64, m: f64, bb: BoundingBox) {
        // Edge cases: if inserting empty objects or inserting the first element of the tree
        if m == 0. { return }
        if self.m == 0. { self.x = x; self.y = y; self.m = m; return }

        // Find the parent to insert this node under
        let mut parent: &mut Self = self;
        let mut parent_bb: BoundingBox = bb;
        let mut quadrant: usize = parent_bb.quadrant(x, y);
        while let Some(_) = &mut parent.children[quadrant] {
            // Update the parent's center of mass
            parent.update_com(x, y, m);

            // Update the bounding box while searching for new parents deeper in the tree
            parent_bb = parent_bb.child(quadrant);
            parent = parent.children[quadrant].as_mut().unwrap();

            // Compute quadrant for next ieration
            quadrant = parent_bb.quadrant(x, y);
        }

        // Leaves must be re-inserted
        if parent.is_leaf() {
            let (px, py, pz, pm) = (parent.x, parent.y, parent.z, parent.m);

            // Edge case: if the parent is too close to the child, don't insert as child
            if (px - x).abs() < EPSILON && (py - y).abs() < EPSILON && (pz - z).abs() < EPSILON { return }

            // Find the center of mass between the two
            parent.update_com(x, y, m);
            let (cx, cy, cz, cm) = (parent.x, parent.y, parent.z, parent.m);

            // Then split until the parent and child are in separate cells
            let mut parent_quadrant = parent_bb.quadrant(px, py);
            while quadrant == parent_quadrant {
                // Create the cell containing both
                parent.new_child(quadrant, cx, cy, cz, cm);
                parent = parent.children[quadrant].as_mut().unwrap();

                // Split the center and continue down
                parent_bb = parent_bb.child(quadrant);
                quadrant = parent_bb.quadrant(x, y);
                parent_quadrant = parent_bb.quadrant(px, py);
            }
            // Once the quadrants are different, insert the parent into its quadrant
            parent.new_child(parent_quadrant, px, py, pz, pm);
        }

        // Insert the new child in the correct quadrant
        parent.new_child(quadrant, x, y, z, m);
    }

    /// Checks if this node is a leaf
    pub fn is_leaf(&self) -> bool {
        for child in &self.children {
            if child.is_some() {
                return false
            }
        }
        true
    }
}

/// Iterator for iterating over all nearby nodes of the tree
pub struct MassQuadtreeIterator<'a> {
    x: f64,
    y: f64,
    theta: f64,
    stack: Vec<(&'a MassOcttree, BoundingBox3D)>
}

/// Implementation of the constructor for the mass quadtree iterator.
impl<'a> MassQuadtreeIterator<'a> {
    /// Constructs a new iterator with the stack initialized to the root.
    pub fn new(x: f64, y: f64, z: f64, z: f64, theta: f64, tree: &'a MassOcttree, bb: BoundingBox3D) -> Self {
        Self {
            x,
            y,
            z,
            theta,
            stack: vec![(tree, bb)]
        }
    }
}

/// Implements the iterator
impl<'a> Iterator for MassQuadtreeIterator<'a> {
    type Item = &'a MassOcttree;

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
    fn next(&mut self) -> Option<&'a MassOcttree> {
        while !self.stack.is_empty() {
            let (node, bb) = self.stack.pop()?;

            let d: f64 = l2(node.x, node.y, self.x, self.y);
            let s: f64 = bb.width();
            if s / d < self.theta || node.is_leaf() { return Some(node) }

            // If not far enough away, add children to the stack.
            for (quadrant, child) in node.children.iter().enumerate() {
                match child {
                    Some(child) => self.stack.push((child, bb.child(quadrant))),
                    None => (),
                }
            }
        }
        None
    }
}



/// Runs a single timestep of the simulation using the Barnes-Hut algorithm.
pub fn nbody_barnes_hut(sim: &mut NBodySimulation3D, dt: f64, theta: f64) {
    let bb: BoundingBox= BoundingBox { min_x, max_x, min_y, max_y, };
    let quadtree: MassQuadtree = MassQuadtree::new(&sim.r, &sim.m, bb);
    // println!("\n\nQuadtree: {:?}", quadtree);

    // For each point
    for i in 0..sim.n {
        sim.a[i] = Vec3::new_zero();
        // println!("r[i] = ({}, {})", sim.rx[i], sim.ry[i]);

        let quadtree_iter =
            MassQuadtreeIterator::new(sim.r[i].x, sim.r[i].y, theta, &quadtree, bb);

        // Get all points that are close enough to treat as individuals
        for node in quadtree_iter {
            let d = Vec3 {
                x: node.x - sim.r[i].x,
                y: node.y - sim.r[i].y,
                z: 0.,
            };
            let d_sqrd: f64 = d.l2_sqrd();
            if d_sqrd < sim.config.min_dist_sqrd {
                continue;
            }

            // if i == 0 { println!("Node: ({}, {}, {})", node.x, node.y, node.m); }

            let inv_d_cubed: f64 = 1. / d_sqrd.powf(3.);
            sim.a[i] += d * node.m * inv_d_cubed;
        }
        // if i == 0 { println!(); }
    }

    sim.integrate(dt);
}