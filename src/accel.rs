#![allow(non_ascii_idents)]

//! This module contains acceleration calculations.

use lin_alg::f64::Vec3;

use crate::{
    gaussian::{GaussianShell, AMP_SCALER},
    units::{KPC, SOLAR_MASS},
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
fn acc_newton_simple(acc_dir: Vec3, src_mass: f64, r: f64) -> Vec3 {
    acc_dir * src_mass / (r.powi(2) + SOFTENING_FACTOR_SQ)
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
        let mass_kg = shell.src_mass * SOLAR_MASS;
        let r_kpc = r * KPC;

        result += acc_newton_simple(acc_dir, shell.src_mass * gauss.value(posit), r);
    }

    result * AMP_SCALER
}

/// An instantaneous acceleration computation. Either Newtonian, or Newtonian modified with MOND.
/// `mond_params` are `(a, a_0)`.
pub fn acc_newton(
    posit: Vec3,
    bodies_other: &[Body],
    id_acted_on: usize,
    mond_a0: Option<f64>,
) -> Vec3 {
    let mut result = Vec3::new_zero();

    for (i, body_src) in bodies_other.iter().enumerate() {
        if i == id_acted_on {
            continue; // self-interaction.
        }

        let acc_diff = body_src.posit - posit;
        let r = acc_diff.magnitude();
        let acc_dir = acc_diff / r; // Unit vec

        // result += acc_newton_simple(acc_dir,  body_src.mass, r);
        // todo: Experimenting using our consistent unit system
        let mass_kg = body_src.mass * SOLAR_MASS;
        let r_kpc = r * KPC;
        result += acc_newton_simple(acc_dir, mass_kg, r_kpc);

        // if let Some(a_0) = mond_a0 {
        //     let x = a / a_0;
        //     let μ = x / (1. + x.powi(2)).sqrt();
        //
        //     accel_newton *= μ;
        // }

        // result += accel_newton;
    }

    if let Some(a_0) = mond_a0 {
        let acc_newton_mag = result.magnitude();

        if acc_newton_mag < 1e-14 {
            // If it's basically zero, no reason to do anything special
            return result;
        }

        fn μ_simple(x: f64) -> f64 {
            x / (1. + x.powi(2)).sqrt()
        }

        // todo: QC all this ChatGPT slop.
        let mut a_mond = acc_newton_mag;
        // Fixed-point iteration:
        //   a_{k+1} = aN / mu(a_k/a0)
        for _ in 0..50 {
            let x = a_mond / a_0;
            let mu_val = μ_simple(x);
            if mu_val.abs() < 1e-30 {
                break; // avoid division by zero
            }
            let new_a = acc_newton_mag / mu_val;
            // Check for convergence
            if ((new_a - a_mond).abs() / a_mond) < 1e-12 {
                a_mond = new_a;
                break;
            }
            a_mond = new_a;
        }

        // 4) Rescale the direction of the Newtonian acceleration to have magnitude a_mond
        //    That is our final (toy) MOND acceleration vector
        let scale = a_mond / acc_newton_mag;

        return result * scale;
    }

    result
}

/// Finds the gravitomagnetic vector potential, analagous to magnetism in Maxwell's equations for EM.
pub fn gravitomagnetic_force(bodies: &[Body]) -> Vec3 {
    // todo: Is this from motion of masses, or rotation? A fn for each?
    Vec3::new_zero()
}
