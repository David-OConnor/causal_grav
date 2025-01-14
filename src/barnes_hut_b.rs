//! https://github.com/BaGreal2/gravitation-particles/tree/main/src

use lin_alg::f64::Vec3;
use rand::Rng;
use crate::barnes_hut::BoundingBox;
use crate::Body;
use crate::units::G;

#[derive(Clone)]
pub struct OctTree {
    bounds: Cube,
    children: Option<[Box<OctTree>; 8]>,
    body: Option<Body>,
    mass: f64,
    m_center_pos: Vec3,
}

impl OctTree {
    pub fn new(bounds: Cube) -> Self {
        let copy_bounds = bounds.clone();
        Self {
            bounds,
            children: None,
            body: None,
            mass: 0.0,
            m_center_pos: Vec3::new(
                                        copy_bounds.top_left_fwd_pos.x + copy_bounds.w / 2.0,
                                        copy_bounds.top_left_fwd_pos.y + copy_bounds.h / 2.0,
                                        copy_bounds.top_left_fwd_pos.z + copy_bounds.d / 2.0,
            )
        }
    }

    fn is_divided(&self) -> bool {
        !self.children.is_none()
    }

    fn subdivide(&mut self) {
        let (x, y, z) = (self.bounds.top_left_fwd_pos.x, self.bounds.top_left_fwd_pos.y, self.bounds.top_left_fwd_pos.z);
        let (w, h, d) = (self.bounds.w, self.bounds.h, self.bounds.d);

        let top_left_f = Cube::new(Vec3::new(x, y, z), w / 2.0, h / 2.0, d/2.);
        let top_right_f = Cube::new(Vec3::new(x + w / 2.0, y, z), w / 2.0, h / 2.0, d/2.);
        let bottom_left_f = Cube::new(Vec3::new(x, y + h / 2.0, z), w / 2.0, h / 2.0, d/2.);
        let bottom_right_f = Cube::new(Vec3::new(x + w / 2.0, y + h / 2.0, z), w / 2.0, h / 2.0, d/2.);

        let top_left_a = Cube::new(Vec3::new(x, y, z), w / 2.0, h / 2.0, d/2.);
        let top_right_a = Cube::new(Vec3::new(x + w / 2.0, y, z + d/2.), w / 2.0, h / 2.0, d/2.);
        let bottom_left_a = Cube::new(Vec3::new(x, y + h / 2.0, z + d/2.), w / 2.0, h / 2.0, d/2.);
        let bottom_right_a = Cube::new(Vec3::new(x + w / 2.0, y + h / 2.0, z + d/2.), w / 2.0, h / 2.0, d/2.);
        
        self.children = Some([
            Box::new(OctTree::new(top_left_f)),
            Box::new(OctTree::new(top_right_f)),
            Box::new(OctTree::new(bottom_left_f)),
            Box::new(OctTree::new(bottom_right_f)),
            Box::new(OctTree::new(top_left_a)),
            Box::new(OctTree::new(top_right_a)),
            Box::new(OctTree::new(bottom_left_a)),
            Box::new(OctTree::new(bottom_right_a)),
        ]);
    }

    pub fn insert(&mut self, body: &Body) {
        if !self.bounds.contains(body) {
            return;
        }

        if self.body.is_none() {
            self.body = Some(body.clone());
        } else {
            if !self.is_divided() {
                self.subdivide();
            }
            // self.children
            //     .as_mut()
            //     .unwrap()
            //     .par_iter_mut()
            //     .for_each(|leaf| {
            //         leaf.as_mut().insert(particle);
            //     });
            for leaf in self.children.as_mut().unwrap().as_mut() {
                leaf.as_mut().insert(body);
            }
            self.update_mass();
        }
    }

    pub fn calculate_force(&mut self, body: &mut Body) {
        if !self.is_divided() {
            if let Some(existent_particle) = &self.body {
                if existent_particle.index != body.index {
                    let attraction_force =
                        body.get_attraction_force(&self.body.as_ref().unwrap());
                    // let attraction_force = self.particle.unwrap().get_attraction_force(particle);
                    body.net_force += attraction_force;
                }
            }
            return;
        }

        let ratio = self.bounds.w / body.get_distance_to(&self.m_center_pos);
        if ratio < 0.5 {
            let attraction_force = body.get_attraction_force(&Body::new(
                self.m_center_pos,
                Vec3::new(0.0, 0.0),
                self.mass,
                1.0,
                1000000,
            ));
            body.net_force += attraction_force;
            return;
        }

        for leaf in self.children.as_mut().unwrap() {
            leaf.calculate_force(body);
        }
    }

    fn update_mass(&mut self) {
        if !self.is_divided() {
            if self.body.is_none() {
                return;
            }
            self.mass = self.body.as_ref().unwrap().mass;
            self.m_center_pos = self.body.as_ref().unwrap().pos;
            return;
        }
        let mut mass_sum = 0.0;
        let mut center_x= 0.0;
        let mut center_y = 0.0;
        let mut center_z = 0.0;

        for leaf in self.children.as_mut().unwrap() {
            leaf.as_mut().update_mass();
            mass_sum += leaf.as_ref().mass;

            center_x += leaf.as_ref().m_center_pos.x * leaf.as_ref().mass;
            center_y += leaf.as_ref().m_center_pos.y * leaf.as_ref().mass;
            center_z += leaf.as_ref().m_center_pos.z * leaf.as_ref().mass;
        }
        self.mass = mass_sum;

        center_x /= mass_sum;
        center_y /= mass_sum;
        center_z /= mass_sum;

        self.m_center_pos = Vec3::new(center_x, center_y, center_z);
    }

    pub fn query(&self, bb: &BoundingBox) -> Vec<Body> {
        let mut results = Vec::new();
        if !self.bounds.intersects(bb) {
            return results;
        }

        if let Some(body) = &self.body {
            if bb.contains(body) {
                results.push(body.clone());
            }
        }
        if self.is_divided() {
            for leaf in self.children.as_ref().unwrap().iter() {
                results.extend(leaf.query(bb));
            }
        }

        results
    }
}
