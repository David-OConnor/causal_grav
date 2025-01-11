#![allow(non_ascii_idents)]

//! This module contains acceleration calculations.

use lin_alg::f64::Vec3;

use crate::{
    gaussian::{GaussianShell, AMP_SCALER},
    units::{A0_MOND, G},
    Body, GravShell, SOFTENING_FACTOR_SQ,
};
// /// Calculate the force acting on a body, given the local environment of gravity shells intersecting it.
// pub fn acc_shells(
//     body: &Body,
//     rays: &[GravRay],
//     shells: &[GravShell],
//     emitter_id: usize,
//     dt: f64,
//     shell_c: f64,
// ) -> Vec3 {
//     let rect = SampleRect::new(body.posit, RAY_SAMPLE_WIDTH);
//     let properties = rect.measure_properties(rays, shells, emitter_id, dt, shell_c);
//
//     // todo: Is this too indirect?
//
//     properties.acc_shell
// }

/// A helper function, where the inputs are precomputed.
fn acc_newton_inner(acc_dir: Vec3, src_mass: f64, r: f64) -> Vec3 {
    acc_dir * G * src_mass / (r.powi(2) + SOFTENING_FACTOR_SQ)
}

pub fn calc_acc_shell(shells: &[GravShell], posit: Vec3, id_acted_on: usize, shell_c: f64) -> Vec3 {
    let mut result = Vec3::new_zero();

    // todo: Once you have more than one body acting on a target, you need to change this, so you get
    // todo exactly 0 or 1 shells per other body.
    for shell in shells {
        if shell.emitter_id == id_acted_on {
            continue;
        }

        let gauss = GaussianShell {
            center: shell.center,
            radius: shell.radius,
            a: 1., // Keep at one; we handle magnitude below.
            c: shell_c,
        };

        let acc_diff = shell.center - posit;
        let r = acc_diff.magnitude();
        let acc_dir = acc_diff / r; // Unit vec

        // todo: Experimenting using our consistent unit system
        // let mass_kg = shell.src_mass * SOLAR_MASS;
        // let r_m = r * KPC;

        result += acc_newton_inner(acc_dir, shell.src_mass * gauss.value(posit), r);
    }

    result * AMP_SCALER
}

#[derive(Clone, Copy, PartialEq)]
pub enum MondFn {
    /// Famaey & Binney. More realistic fits than the standard one. `x` is a_Newton / a_0.
    Simple,
    /// Sanders & Noordermeer. `x` is a_Newton / a_0.
    Standard,
}

impl MondFn {
    pub fn μ(&self, x: f64) -> f64 {
        match self {
            Self::Simple => x / (1. + x),
            Self::Standard => x / (1. + x.powi(2)).sqrt(),
        }
    }
}

/// An instantaneous acceleration computation. Either Newtonian, or Newtonian modified with MOND.
/// `mond_params` are `(a, a_0)`.
pub fn acc_newton(
    posit: Vec3,
    bodies_other: &[Body],
    id_acted_on: usize,
    mond: Option<MondFn>,
) -> Vec3 {
    let mut result = Vec3::new_zero();

    for (i, body_src) in bodies_other.iter().enumerate() {
        if i == id_acted_on {
            continue; // self-interaction.
        }

        let acc_diff = body_src.posit - posit;
        let r = acc_diff.magnitude();
        let acc_dir = acc_diff / r; // Unit vec

        let mut acc = acc_newton_inner(acc_dir, body_src.mass, r);

        if let Some(mond_fn) = mond {
            // todo: This may not be correct. r may be the wrong arbgument?
            // todo: Also, div/0 error.
            let x = acc.magnitude() / A0_MOND;
            acc = acc / mond_fn.μ(x)
        }

        result += acc;
    }

    result
}

/// Finds the gravitomagnetic vector potential, analagous to magnetism in Maxwell's equations for EM.
pub fn gravitomagnetic_force(bodies: &[Body]) -> Vec3 {
    // todo: Is this from motion of masses, or rotation? A fn for each?
    Vec3::new_zero()
}
