use rayon::prelude::*;

use crate::{
    accel, barnes_hut,
    barnes_hut::{BhConfig, Tree},
    Body, ForceModel, GravShell,
};

pub fn integrate_rk4(
    bodies: &mut [Body],
    shells: &[GravShell],
    dt: f64,
    gauss_c: f64,
    tree: &Tree,
    bh_config: &BhConfig,
    force_model: ForceModel,
    softening_factor_sq: f64,
) {
    // todo: You must update this acc using your tree approach.
    let acc = |id_target, posit_target| match force_model {
        ForceModel::Newton => {
            // accel::acc_newton(posit, id, &bodies_other, None, softening_factor_sq)
            barnes_hut::acc_newton_bh(
                posit_target,
                id_target,
                tree,
                bh_config,
                None,
                softening_factor_sq,
            )
        }
        ForceModel::GaussShells => accel::calc_acc_shell(
            shells,
            posit_target,
            id_target,
            gauss_c,
            softening_factor_sq,
        ),
        ForceModel::Mond(mond_fn) => {
            // accel::acc_newton(posit_target, id_target, &bodies_other, Some(mond_fn), softening_factor_sq)
            barnes_hut::acc_newton_bh(
                posit_target,
                id_target,
                tree,
                bh_config,
                None,
                softening_factor_sq,
            )
        }
    };

    bodies.par_iter_mut().enumerate().for_each(|(id, body)| {
        // Step 1: Calculate the k-values for position and velocity
        body.accel = acc(id, body.posit);

        let k1_v = body.accel * dt;
        let k1_pos = body.vel * dt;

        let body_pos_k2 = body.posit + k1_pos * 0.5;
        let k2_v = acc(id, body_pos_k2) * dt;
        let k2_pos = (body.vel + k1_v * 0.5) * dt;

        let body_pos_k3 = body.posit + k2_pos * 0.5;
        let k3_v = acc(id, body_pos_k3) * dt;
        let k3_pos = (body.vel + k2_v * 0.5) * dt;

        let body_pos_k4 = body.posit + k3_pos;
        let k4_v = acc(id, body_pos_k4) * dt;
        let k4_pos = (body.vel + k3_v) * dt;

        // Step 2: Update position and velocity using weighted average of k-values
        body.vel += (k1_v + k2_v * 2. + k3_v * 2. + k4_v) / 6.;
        body.posit += (k1_pos + k2_pos * 2. + k3_pos * 2. + k4_pos) / 6.;
    });
}
