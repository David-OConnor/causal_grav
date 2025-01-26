#![allow(non_ascii_idents)]

//! This module contains acceleration calculations.

use lin_alg::f64::Vec3;
use rayon::prelude::*;

use crate::{
    grav_shell::{GravShell, AMP_SCALER},
    units::{A0_MOND, C, G},
    Body,
};

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

/// The most fundamental part of Newtonian acceleration calculation.
/// `acc_dir` is a unit vector.
pub fn acc_newton_inner(acc_dir: Vec3, src_mass: f64, dist: f64, softening_factor_sq: f64) -> Vec3 {
    acc_dir * G * src_mass / (dist.powi(2) + softening_factor_sq)
}

/// This optionally applies MOND to our basic Newton acceleration.
pub fn acc_newton_inner_with_mond(
    acc_dir: Vec3,
    mass: f64,
    dist: f64,
    mond: Option<MondFn>,
    softening_factor_sq: f64,
) -> Vec3 {
    let mut acc = acc_newton_inner(acc_dir, mass, dist, softening_factor_sq);

    if let Some(mond_fn) = mond {
        let x = acc.magnitude() / A0_MOND;
        acc /= mond_fn.μ(x);
    }
    return acc;
}

pub fn calc_acc_shell(
    shells: &[GravShell],
    posit: Vec3,
    id_target: usize,
    shell_c: f64,
    softening_factor_sq: f64,
) -> Vec3 {
    // todo: Once you have more than one body acting on a target, you need to change this, so you get
    // todo exactly 0 or 1 shells per other body.

    // todo: DRY with non-shell.
    shells.par_iter().filter_map(|shell| {
        if shell.source_id == id_target {
            return None; // Skip self-interaction.
        }

        let source_posit = shell.center;
        let t_since_creation = shell.radius / C;
        // let source_posit = shell.center +  shell.body_vel * t_since_creation;
        // let source_posit = shell.center + shell.body_vel * t_since_creation + shell.body_acc / 2. * t_since_creation.powi(2);

        let acc_diff = source_posit - posit;
        let dist = acc_diff.magnitude();
        let acc_dir = acc_diff / dist; // Unit vec

        Some(acc_newton_inner(
            acc_dir,
            shell.value(posit, shell_c),
            dist,
            softening_factor_sq,
        ))
    })
        .reduce(Vec3::new_zero, |acc, elem| acc + elem) // Sum the contributions.
    * AMP_SCALER
}

/// An instantaneous acceleration computation, from all sources, on a single target.
/// Either Newtonian, or Newtonian modified with MOND.
/// `mond_params` are `(a, a_0)`.
///
/// Uses Rayon for parallel execution. The functional approach is required for use with Rayon.
pub fn acc_newton(
    posit_target: Vec3,
    id_target: usize,
    bodies_other: &[Body],
    mond: Option<MondFn>,
    softening_factor_sq: f64,
) -> Vec3 {
    // Compute the result in parallel and then sum the contributions.
    bodies_other
        .par_iter()
        .enumerate()
        .filter_map(|(i, body_source)| {
            if i == id_target {
                return None; // Skip self-interaction.
            }

            let acc_diff = body_source.posit - posit_target;
            let dist = acc_diff.magnitude();
            let acc_dir = acc_diff / dist; // Unit vector.

            Some(acc_newton_inner_with_mond(
                acc_dir,
                body_source.mass,
                dist,
                mond,
                softening_factor_sq,
            ))
        })
        .reduce(Vec3::new_zero, |acc, elem| acc + elem) // Sum the contributions.
}

/// Finds the gravitomagnetic vector potential, analagous to magnetism in Maxwell's equations for EM.
pub fn gravitomagnetic_force(bodies: &[Body]) -> Vec3 {
    // todo: Is this from motion of masses, or rotation? A fn for each?
    Vec3::new_zero()
}
