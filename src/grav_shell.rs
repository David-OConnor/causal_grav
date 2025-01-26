use lin_alg::f64::Vec3;

use crate::{gaussian::GaussianShell, units::C};

// Find a value of C, given spacing and amplitude, that provides a good balance between distribution
// uniformity, and sharp edges.
//
// A higher coefficient results in a more uniform distribution, at the cost of responsiveness at the edges.
// 0.55: Sharper falloff. 0.6: More uniform distribution.

// todo: Adjust this approach A/R.

pub const COEFF_C: f64 = 0.6;
// pub const COEFF_C: f64 = 0.55;

pub const AMP_SCALER: f64 = 0.6649; // Based on COEFF = 0.6. Found from trial + error using `gauss_spacing.py`.
                                    // pub const AMP_SCALER: f64 = 0.7253; // Based on COEFF = 0.55. Found from trial + error using `gauss_spacing.py`.

#[derive(Debug, Clone)]
/// Represents gravitational potential, as a shell. This allows for gravitational force to have finite speed,
/// and act locally. We combine gaussians to achieve a uniform-like distribution.
///
/// See the S. Carlip "Abberation and the speed of gravity" for ideas on how to model the apparent lack of
/// abberation in real models, in conjunction with a finite speed of gravity. Should these velocity-dependent
/// abberation-cancelling terms be included at shell init from data on the source body alone, or do
/// they develop as the shell propogates, e.g. thorugh interaction with other shells?
///
/// First step: Add first and second-order (velocity and acc) terms from the body's init condition.
pub struct GravShell {
    pub source_id: usize,
    pub center: Vec3,
    pub radius: f64,
    pub src_mass: f64,
    /// Body velocity at creation. Experimenting with velocity-dependent effects, which may be required.
    /// todo: Alternate model: Maybe the wave doesn't get this information at the start, but it gets
    /// todo it on the way from interaction with other waves.
    pub body_vel: Vec3,
    /// Is this also required?
    /// *photon-rocket* accel doesn't affect it, but gravitational accel?
    pub body_acc: Vec3,
}

impl GravShell {
    /// Expand the radius at C, in one timestep.
    pub fn iter_t(&mut self, dt: f64) {
        self.radius += C * dt;
    }

    pub fn value(&self, posit: Vec3, gauss_c: f64) -> f64 {
        let gauss = GaussianShell {
            center: self.center,
            radius: self.radius,
            a: self.src_mass,
            c: gauss_c,
        };

        gauss.value(posit)
    }
}

// pub const MAX_SHELL_R: f64 = 50.; // todo: Adjust this approach A/R.
pub const MAX_SHELL_R: f64 = 20.;
