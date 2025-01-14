//! Adapted from https://gitlab.scd31.com/stephen/nbody_barnes_hut/-/blob/master/src/barnes_hut_3d.rs?ref_type=heads

use crate::Body;
use lin_alg::f64::Vec3;
use crate::barnes_hut::BoundingBox;

struct OctNode {
    children_base: usize,
    children_mask: u8,
    quad_size_squared: f64,
    mass: f64,
}

impl OctNode {
    fn new(min: Vec3, max: Vec3) -> Self {
        let s = (min - max).magnitude_squared();

        OctNode {
            children_base: 0,
            children_mask: 0,
            mass: 0.0,
            quad_size_squared: s,
        }
    }
}

pub struct OctTree {
    nodes: Vec<OctNode>,
    min_coords: Vec<Vec3>,
    max_coords: Vec<Vec3>,
    avg_weighted_positions: Vec<Vec3>,
    center_of_masses: Vec<Vec3>,
    theta_squared: f64,
}

impl OctTree {
    /// Constructs a tree from an array of particles.
    /// Also takes in theta, which - roughly speaking - describes how accurate the simulation should be.
    /// For more information, check [here](https://en.wikipedia.org/wiki/Barnes%E2%80%93Hut_simulation).
    ///
    /// > "Whether a node is or isn't sufficiently far away from a body, depends on the quotient s / d, where s is the width of the region represented by the internal node, and d is the distance between the body and the node's center of mass. The node is sufficiently far away when this ratio is smaller than a threshold value θ. The parameter θ determines the accuracy of the simulation; larger values of θ increase the speed of the simulation but decreases its accuracy. If θ = 0, no internal node is treated as a single body and the algorithm degenerates to a direct-sum algorithm."
    ///
    /// A value of 0.5 is quite common.
    pub fn new(bodies: &[&Body], theta: f64) -> Self {
        let limits = BoundingBox::new(bodies);

        let mut tree = OctTree::construct_empty_tree(min, max, theta, bodies.len());

        tree.add_particles_to_specific_node(bodies, 0);
        //calculate center of masses
        for (i, n) in &mut tree.nodes.iter().enumerate() {
            tree.center_of_masses
                .push(tree.avg_weighted_positions[i] / n.mass);
        }

        tree
    }

    pub fn update_theta(&mut self, theta: f64) {
        assert!(theta >= 0.0);
        assert!(theta <= 1.0);

        self.theta_squared = theta * theta;
    }

    /// Calculates the force on a particle.
    ///
    /// Takes in the position of the particle, any attributes you wish to pass in, and a closure
    ///
    /// The closure will need to take in the distance squared, the mass, the distance vector, and the other attributes which you passed in.
    /// It will need to return a force (or acceleration) vector
    /// Please note that the distance vector is *not* normalized! It is your responsibility to do this.
    pub fn calc_forces_on_particle<T: Copy, F: Fn(f64, f64, Vec3, T) -> Vec3 + Copy>(
        &self,
        particle_position: Vec3,
        other_attributes: T,
        calc: F,
    ) -> Vec3 {
        self.barnes_hut_specific_node(particle_position, 0, other_attributes, calc)
    }

    fn barnes_hut_specific_node<T: Copy, F: Fn(f64, f64, Vec3, T) -> Vec3 + Copy>(
        &self,
        particle_position: Vec3,
        cur_node: usize,
        other_attributes: T,
        calc: F,
    ) -> Vec3 {
        let c_node = &self.nodes[cur_node];

        let s_squared = c_node.quad_size_squared;
        let dis_vec = self.center_of_masses[cur_node] - particle_position;
        let d_squared = dis_vec.magnitude_squared();
        let theta_squared = self.theta_squared;

        if c_node.children_mask == 0 || s_squared < theta_squared * d_squared {
            if d_squared.abs() < 0.00001 {
                return Vec3::new_zero();
            }

            calc(d_squared, c_node.mass, dis_vec, other_attributes)
        } else {
            let mut c = c_node.children_base;
            //Divide and conquer
            let mut sum: Vec3 = Vec3::new_zero();
            for i in 0..8 {
                if c_node.children_mask & (1 << i) != 0 {
                    sum +=
                        self.barnes_hut_specific_node(particle_position, c, other_attributes, calc);
                    c += 1;
                }
            }

            sum
        }
    }

    fn construct_empty_tree(min: Vec3, max: Vec3, theta: f64, n: usize) -> Self {
        let mut o = OctTree {
            nodes: Vec::with_capacity(n),
            min_coords: Vec::with_capacity(n),
            max_coords: Vec::with_capacity(n),
            avg_weighted_positions: Vec::with_capacity(n),
            center_of_masses: Vec::with_capacity(n),
            theta_squared: theta * theta,
        };
        o.nodes.push(OctNode::new(min, max)); // Create our root node
        o.min_coords.push(min);
        o.max_coords.push(max);
        o.avg_weighted_positions.push(Vec3::new_zero());
        o
    }

    fn add_particles_to_specific_node(&mut self, bodies: &[&Body], cur_node: usize) {
        if bodies.is_empty() {
            return;
        }

        if bodies.len() == 1 {
            self.nodes[cur_node].mass += bodies[0].mass;
            self.avg_weighted_positions[cur_node] += bodies[0].posit * bodies[0].mass;
            return;
        }

        let center_point = (self.min_coords[cur_node] + self.max_coords[cur_node]) / 2.0;

        let mut particle_octs: [Vec<&Body>; 8] = [
            Vec::with_capacity(bodies.len() / 4),
            Vec::with_capacity(bodies.len() / 4),
            Vec::with_capacity(bodies.len() / 4),
            Vec::with_capacity(bodies.len() / 4),
            Vec::with_capacity(bodies.len() / 4),
            Vec::with_capacity(bodies.len() / 4),
            Vec::with_capacity(bodies.len() / 4),
            Vec::with_capacity(bodies.len() / 4),
        ];

        for &particle in bodies {
            // Update ourself
            self.nodes[cur_node].mass += particle.mass;
            self.avg_weighted_positions[cur_node] += particle.posit * particle.mass;

            let octant_index = OctTree::get_id_from_center(center_point, particle.posit);
            particle_octs[octant_index].push(particle);
        }

        //Recurse
        self.nodes[cur_node].children_base = self.nodes.len();
        for (i, particle_oct) in particle_octs.iter().enumerate() {
            if !particle_oct.is_empty() {
                let (min, max) = OctTree::get_octant_bounding_box_from_id(
                    self.min_coords[cur_node],
                    self.max_coords[cur_node],
                    center_point,
                    i,
                );

                self.nodes.push(OctNode::new(min, max));
                self.min_coords.push(min);
                self.max_coords.push(max);
                self.avg_weighted_positions.push(Vec3::new_zero());
                self.nodes[cur_node].children_mask |= 1 << i;
            }
        }

        let mut c = self.nodes[cur_node].children_base;
        for particle_oct in particle_octs.iter() {
            if !particle_oct.is_empty() {
                self.add_particles_to_specific_node(particle_oct, c);
                c += 1;
            }
        }
    }

    fn get_id_from_center(center: Vec3, point: Vec3) -> usize {
        let offset = point - center;
        //We can look at the sign of the components of the offset to figure out which octant it is in.
        let x_offset = if offset.x > 0.0 { 0 } else { 1 };
        let y_offset = if offset.y > 0.0 { 0 } else { 1 };
        let z_offset = if offset.z > 0.0 { 0 } else { 1 };
        x_offset + y_offset * 2 + z_offset * 4 // basic binary stuff here
    }

    fn get_octant_bounding_box_from_id(
        min: Vec3,
        max: Vec3,
        center: Vec3,
        idx: usize,
    ) -> (Vec3, Vec3) {
        let mut min_coord: Vec3 = Vec3::new_zero();
        let mut max_coord: Vec3 = Vec3::new_zero();
        if (idx & 1) != 0 {
            min_coord.x = min.x;
            max_coord.x = center.x;
        } else {
            min_coord.x = center.x;
            max_coord.x = max.x;
        }

        if (idx & 2) != 0 {
            min_coord.y = min.y;
            max_coord.y = center.y;
        } else {
            min_coord.y = center.y;
            max_coord.y = max.y;
        }

        if (idx & 4) != 0 {
            min_coord.z = min.z;
            max_coord.z = center.z;
        } else {
            min_coord.z = center.z;
            max_coord.z = max.z;
        }

        (min_coord, max_coord)
    }
}

