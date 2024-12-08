//! This module contains acceleration calculations.

use lin_alg::f64::Vec3;

use crate::{gaussian::GaussianShell, Body, GravRay, GravShell, SampleRect, RAY_SAMPLE_WIDTH};

/// Calculate the force acting on a body, given the local environment of gravity rays around it.
pub fn accel(
    body: &mut Body,
    rays: &[GravRay],
    shells: &[GravShell],
    emitter_id: usize,
    dt: f64,
) -> Vec3 {
    // let ray_density_grad = density_gradient(body.posit, rays);

    let rect = SampleRect::new(body.posit, RAY_SAMPLE_WIDTH);
    let properties = rect.measure_properties(rays, shells, emitter_id);

    // println!("Prop: {:?}", properties);

    // todo: The rate constant must take into account the angle of the rays; dividing by dt
    // todo is likely only to work in the limit of infinitely small d_theta for ray emission, etc.
    // todo: Divide by a rate constant.
    let rate_const = 460.;
    // let rate_const = 1. / dt; // todo: we can pre-calc this, but not a big deal.

    properties.ray_net_direction * properties.ray_density * rate_const

    // todo: We need to take into account the destination body's mass, not just
    // todo for inertia, but for attraction... right?

    // todo: We may be missing part of the puzzle.
    // let force = gradient * body.mass;

    // a = F / m
    // force / body.mass

    // body.V_acting_on = ray_density_grad;

    // ray_density_grad
}

/// Calculate the force acting on a body, given the local environment of gravity shells intersecting it.
pub fn accel_shells(
    body: &Body,
    rays: &[GravRay],
    shells: &[GravShell],
    emitter_id: usize,
    dt: f64,
) -> Vec3 {
    let rect = SampleRect::new(body.posit, RAY_SAMPLE_WIDTH);
    let properties = rect.measure_properties(rays, shells, emitter_id);

    // todo: Is this too indirect?

    // println!("Prop: {:?}", properties);

    properties.acc_shell
}

// todo: Put back once you have an accel computation.
// fn integrate_rk4(bodies: &mut [Body], dt: f64) {
//     for body in bodies.iter_mut() {
//         // Step 1: Calculate the k-values for position and velocity
//         let k1_v = body.accel * dt;
//         let k1_posit = body.vel * dt;
//
//         let k2_v = compute_acceleration(body.posit + k1_posit * 0.5) * dt;
//         let k2_posit = (body.vel + k1_v * 0.5) * dt;
//
//         let k3_v = compute_acceleration(body.posit + k2_posit * 0.5) * dt;
//         let k3_posit = (body.vel + k2_v * 0.5) * dt;
//
//         let k4_v = compute_acceleration(body.posit + k3_posit) * dt;
//         let k4_posit = (body.vel + k3_v) * dt;
//
//         // Step 2: Update position and velocity using weighted average of k-values
//         body.vel += (k1_v + k2_v * 2. + k3_v * 2. + k4_v) / 6.;
//         body.posit += (k1_posit + k2_posit * 2. + k3_posit * 2. + k4_posit) / 6.;
//     }
// }

pub fn calc_acc_shell(shells: &[GravShell], posit: Vec3, emitter_id: usize) -> Vec3 {
    let mut shell_value = 0.;
    let mut shell_acc_sum = Vec3::new_zero();

    // todo: Once you have more than one body acting on a target, you need to change this, so you get
    // todo exactly 0 or 1 shells per other body.
    for shell in shells {
        if shell.emitter_id == emitter_id {
            continue;
        }

        let gauss = GaussianShell {
            center: shell.center,
            radius: shell.radius,
            a: 1.,  // Keep at one; scales magnitude
            c: 0.1, // todo: Exper with c.
        };
        // shell_value += shell.src_mass * gauss.value(posit) / (shell.center - posit).magnitude_squared();
        // todo: Experimenting.
        shell_value +=
            shell.src_mass * gauss.value(posit) / (shell.center - posit).magnitude().powi(3);

        // todo: QC what the acc dir from the shell is.
        let shell_acc_dir = (posit - shell.center).to_normalized();
        shell_acc_sum += shell_acc_dir * shell_value; // todo: QC order

        // if shell.intersects_rect(self) {
        //     println!("Intersects: {:?}", shell);
        //     acc_shell + (shell.center - center) * shell.src_mass / shell.radius.powi(2);
        // }
    }

    let compensator = 0.0001 / 25.;

    shell_acc_sum * compensator
}
