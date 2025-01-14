//! https://github.com/DeadlockCode/barnes-hut/tree/improved/src


use std::ops::Range;
use lin_alg::f64::Vec3;
use crate::barnes_hut::BoundingBox;
use crate::Body;

#[derive(Clone)]
pub struct Node {
    pub children: usize,
    pub next: usize,
    pub pos: Vec3,
    pub mass: f64,
    pub oct: Oct,
    pub bodies: Range<usize>,
}

impl Node {
    pub fn new(next: usize, quad: Oct, bodies: Range<usize>) -> Self {
        Self {
            children: 0,
            next,
            pos: Vec3::new_zero(),
            mass: 0.0,
            oct: quad,
            bodies,
        }
    }

    pub fn is_leaf(&self) -> bool {
        self.children == 0
    }

    pub fn is_branch(&self) -> bool {
        self.children != 0
    }

    pub fn is_empty(&self) -> bool {
        self.mass == 0.0
    }
}

pub struct OctTree {
    pub t_sq: f64,
    pub e_sq: f64,
    pub leaf_capacity: usize,
    pub nodes: Vec<Node>,
    pub parents: Vec<usize>,
}

impl OctTree {
    pub const ROOT: usize = 0;

    pub fn new(theta: f64, epsilon: f64, leaf_capacity: usize) -> Self {
        Self {
            t_sq: theta * theta,
            e_sq: epsilon * epsilon,
            leaf_capacity,
            nodes: Vec::new(),
            parents: Vec::new(),
        }
    }

    pub fn clear(&mut self) {
        self.nodes.clear();
        self.parents.clear();
    }

    pub fn subdivide(&mut self, node: usize, bodies: &mut [Body], range: Range<usize>) {
        let center = self.nodes[node].oct.center;

        let mut split = [range.start, 0, 0, 0, range.end];

        let predicate = |body: &Body| body.posit.y < center.y;
        split[2] = split[0] + bodies[split[0]..split[4]].partition(predicate);

        let predicate = |body: &Body| body.posit.x < center.x;
        split[1] = split[0] + bodies[split[0]..split[2]].partition(predicate);
        split[3] = split[2] + bodies[split[2]..split[4]].partition(predicate);

        self.parents.push(node);
        let children = self.nodes.len();
        self.nodes[node].children = children;

        let nexts = [
            children + 1,
            children + 2,
            children + 3,
            self.nodes[node].next,
        ];
        let quads = self.nodes[node].oct.subdivide();
        for i in 0..4 {
            let bodies = split[i]..split[i + 1];
            self.nodes.push(Node::new(nexts[i], quads[i], bodies));
        }
    }

    pub fn propagate(&mut self) {
        for &node in self.parents.iter().rev() {
            let i = self.nodes[node].children;

            self.nodes[node].pos = self.nodes[i].pos
                + self.nodes[i + 1].pos
                + self.nodes[i + 2].pos
                + self.nodes[i + 3].pos;

            self.nodes[node].mass = self.nodes[i].mass
                + self.nodes[i + 1].mass
                + self.nodes[i + 2].mass
                + self.nodes[i + 3].mass;
        }
        for node in &mut self.nodes {
            node.pos /= node.mass.max(f64::MIN_POSITIVE);
        }
    }

    pub fn build(&mut self, bodies: &mut [Body]) {
        self.clear();

        let quad = BoundingBox::new(bodies);
        self.nodes.push(Node::new(0, quad, 0..bodies.len()));

        let mut node = 0;
        while node < self.nodes.len() {
            let range = self.nodes[node].bodies.clone();
            if range.len() > self.leaf_capacity {
                self.subdivide(node, bodies, range);
            } else {
                for i in range {
                    self.nodes[node].pos += bodies[i].posit * bodies[i].mass;
                    self.nodes[node].mass += bodies[i].mass;
                }
            }
            node += 1;
        }

        self.propagate();
    }

    pub fn acc(&self, pos: Vec3, bodies: &[Body]) -> Vec3 {
        let mut acc = Vec3::new_zero();

        let mut node = Self::ROOT;
        loop {
            let n = self.nodes[node].clone();

            let d = n.pos - pos;
            let d_sq = d.mag_sq();

            if n.oct.size * n.oct.size < d_sq * self.t_sq {
                let denom = (d_sq + self.e_sq) * d_sq.sqrt();
                acc += d * (n.mass / denom);

                if n.next == 0 {
                    break;
                }
                node = n.next;
            } else if n.is_leaf() {
                for i in n.bodies {
                    let body = &bodies[i];
                    let d = body.posit - pos;
                    let d_sq = d.mag_sq();

                    let denom = (d_sq + self.e_sq) * d_sq.sqrt();
                    acc += d * (body.mass / denom).min(f64::MAX);
                }

                if n.next == 0 {
                    break;
                }
                node = n.next;
            } else {
                node = n.children;
            }
        }

        acc
    }
}