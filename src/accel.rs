#![allow(non_ascii_idents)]

//! This module contains acceleration calculations.

use lin_alg::f64::Vec3;

use crate::{gaussian::{GaussianShell, AMP_SCALER}, Body, GravShell, SOFTENING_FACTOR_SQ};

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

        let acc_dir = shell.center - posit;

        let acc_mag = acc_dir.magnitude();

        result += acc_dir * shell.src_mass * gauss.value(posit) / acc_mag.powi(3);
    }

    result * AMP_SCALER
}

/// An instantaneous acceleration computation. Either Newtonian, or Newtonian modified with MOND.
/// `mond_params` are `(a, a_0)`.
pub fn acc_newton(posit: Vec3, bodies_other: &[Body], id_acted_on: usize, mond_a0: Option<f64>) -> Vec3 {
    let mut result = Vec3::new_zero();

    for (i, body_src) in bodies_other.iter().enumerate() {
        if i == id_acted_on {
            continue; // self-interaction.
        }

        let acc_dir = body_src.posit - posit;

        let mut accel_newton = acc_dir * body_src.mass / (acc_dir.magnitude().powi(3) + SOFTENING_FACTOR_SQ);

        // if let Some(a_0) = mond_a0 {
        //     let x = a / a_0;
        //     let μ = x / (1. + x.powi(2)).sqrt();
        //
        //     accel_newton *= μ;
        // }

        result += accel_newton;
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
        let scale = a_mond /acc_newton_mag;

        return result * scale;
    }

    result
}

/// Finds the gravitomagnetic vector potential, analagous to magnetism in Maxwell's equations for EM.
pub fn gravitomagnetic_force(bodies: &[Body]) -> Vec3 {
    // todo: Is this from motion of masses, or rotation? A fn for each?
    Vec3::new_zero()
}
