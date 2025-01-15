#![allow(non_ascii_idents)]

//! This module contains acceleration calculations.

use lin_alg::f64::Vec3;

use crate::{
    gaussian::{GaussianShell, AMP_SCALER},
    units::{A0_MOND, G},
    Body, GravShell,
};

use rayon::prelude::*;

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
/// `acc_dir` is a unit vector.
pub(crate) fn acc_newton_inner(acc_dir: Vec3, src_mass: f64, dist: f64, softening_factor_sq: f64) -> Vec3 {
    acc_dir * G * src_mass / (dist.powi(2) + softening_factor_sq)
}

pub fn calc_acc_shell(shells: &[GravShell], posit: Vec3, id_acted_on: usize, shell_c: f64, softening_factor_sq: f64) -> Vec3 {
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
        let dist = acc_diff.magnitude();
        let acc_dir = acc_diff / dist; // Unit vec

        // todo: Experimenting using our consistent unit system
        // let mass_kg = shell.src_mass * SOLAR_MASS;
        // let r_m = r * KPC;

        result += acc_newton_inner(acc_dir, shell.src_mass * gauss.value(posit), dist, softening_factor_sq);
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

/// An instantaneous acceleration computation, from all actors, on a single body acted on.
/// Either Newtonian, or Newtonian modified with MOND.
/// `mond_params` are `(a, a_0)`.
pub fn acc_newton(
    posit_acted_on: Vec3,
    id_acted_on: usize,
    bodies_other: &[Body],
    mond: Option<MondFn>,
    softening_factor_sq: f64,
) -> Vec3 {
    let mut result = Vec3::new_zero();

    // Iterate over bodies acting on our target.
    // for (i, body_actor) in bodies_other.iter().enumerate() {
    for (i, body_actor) in bodies_other.iter().enumerate() {
        if i == id_acted_on {
            continue; // self-interaction.
        }

        let acc_diff = body_actor.posit - posit_acted_on;
        let r = acc_diff.magnitude();
        let acc_dir = acc_diff / r; // Unit vec

        let mut acc = acc_newton_inner(acc_dir, body_actor.mass, r, softening_factor_sq);

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

/// Uses Rayon. Restructured to use its required functional approach.
/// ~10x slower than the non-Rayon approach, from initial tests.
pub fn acc_newton_parallel(
    posit_acted_on: Vec3,
    id_acted_on: usize,
    bodies_other: &[Body],
    mond: Option<MondFn>,
    softening_factor_sq: f64,
) -> Vec3 {
    // Compute the result in parallel and then sum the contributions.
    bodies_other
        .par_iter()
        .enumerate()
        .filter_map(|(i, body_actor)| {
            if i == id_acted_on {
                return None; // Skip self-interaction.
            }

            let acc_diff = body_actor.posit - posit_acted_on;
            let r = acc_diff.magnitude();
            let acc_dir = acc_diff / r; // Unit vector.

            let mut acc = acc_newton_inner(acc_dir, body_actor.mass, r, softening_factor_sq);

            if let Some(mond_fn) = mond {
                let x = acc.magnitude() / A0_MOND;
                acc = acc / mond_fn.μ(x);
            }

            Some(acc)
        })
        .reduce(Vec3::new_zero, |acc, elem| acc + elem) // Sum the contributions.
}

/// Finds the gravitomagnetic vector potential, analagous to magnetism in Maxwell's equations for EM.
pub fn gravitomagnetic_force(bodies: &[Body]) -> Vec3 {
    // todo: Is this from motion of masses, or rotation? A fn for each?
    Vec3::new_zero()
}
