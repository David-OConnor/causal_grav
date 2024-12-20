//! This module contains acceleration calculations.

use lin_alg::f64::Vec3;

use crate::{
    gaussian::{GaussianShell, AMP_SCALER},
    Body, GravRay, GravShell, SampleRect, RAY_SAMPLE_WIDTH,
};

/// Calculate the force acting on a body, given the local environment of gravity rays around it.
pub fn acc_rays(
    body: &mut Body,
    rays: &[GravRay],
    shells: &[GravShell],
    emitter_id: usize,
    shell_c: f64,
) -> Vec3 {
    // let ray_density_grad = density_gradient(body.posit, rays);

    let rect = SampleRect::new(body.posit, RAY_SAMPLE_WIDTH);
    let properties = rect.measure_properties(rays, shells, emitter_id, shell_c);

    // println!("Prop: {:?}", properties);

    // todo: The rate constant must take into account the angle of the rays; dividing by dt
    // todo is likely only to work in the limit of infinitely small d_theta for ray emission, etc.
    // todo: Divide by a rate constant.
    let rate_const = 460.;
    // let rate_const = 1. / dt; // todo: we can pre-calc this, but not a big deal.

    properties.ray_net_direction * properties.ray_density * rate_const
}

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

pub fn calc_acc_shell(shells: &[GravShell], posit: Vec3, emitter_id: usize, shell_c: f64) -> Vec3 {
    let mut result = Vec3::new_zero();

    // todo: Once you have more than one body acting on a target, you need to change this, so you get
    // todo exactly 0 or 1 shells per other body.
    for shell in shells {
        if shell.emitter_id == emitter_id {
            continue;
        }

        let gauss = GaussianShell {
            center: shell.center,
            radius: shell.radius,
            a: 1., // Keep at one; we handle magnitude below.
            c: shell_c,
        };

        let acc_dir = shell.center - posit;

        result += acc_dir * shell.src_mass * gauss.value(posit) / acc_dir.magnitude().powi(3);
    }

    result * AMP_SCALER
}

/// An instantaneous-accel control.
pub fn calc_acc_inst(posit: Vec3, bodies_other: &[Body], emitter_id: usize) -> Vec3 {
    let mut result = Vec3::new_zero();

    for (i, body_src) in bodies_other.iter().enumerate() {
        if i == emitter_id {
            continue; // self-interaction.
        }

        let acc_dir = body_src.posit - posit;

        result += acc_dir * body_src.mass / acc_dir.magnitude().powi(3);

        // println!("RESULT: {:?}, Acc dir: {:?}", result, acc_dir);
    }

    result
}
