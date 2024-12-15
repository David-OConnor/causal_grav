use lin_alg::f64::Vec3;

use crate::{accel, Body, GravRay, GravShell};

pub fn integrate_rk4(bodies: &mut Vec<Body>, shells: &[GravShell], dt: f64, shell_c: f64) {
    let bodies_other = bodies.clone(); // todo: I don't like this. Avoids mut error.
    for (id, body) in bodies.iter_mut().enumerate() {
        // Step 1: Calculate the k-values for position and velocity
        // body.accel = accel::calc_acc_shell(shells, body.posit, id, dt, shell_c) * dt;
        body.accel = accel::calc_acc_inst(body.posit, &bodies_other, id) * dt;

        let k1_v = body.accel * dt;
        let k1_posit = body.vel * dt;

        let body_pos_k2 = body.posit + k1_posit * 0.5;
        // let k2_v = accel::calc_acc_shell(shells, body_pos_k2, id, dt, shell_c) * dt;
        let k2_v = accel::calc_acc_inst(body_pos_k2, &bodies_other, id) * dt;
        let k2_posit = (body.vel + k1_v * 0.5) * dt;

        let body_pos_k3 = body.posit + k2_posit * 0.5;
        // let k3_v = accel::calc_acc_shell(shells, body_pos_k3, id, dt, shell_c) * dt;
        let k3_v = accel::calc_acc_inst(body_pos_k3, &bodies_other, id)  * dt;
        let k3_posit = (body.vel + k2_v * 0.5) * dt;

        let body_pos_k4 = body.posit + k3_posit;
        // let k4_v = accel::calc_acc_shell(shells, body_pos_k4, id, dt, shell_c) * dt;
        let k4_v = accel::calc_acc_inst(body_pos_k4, &bodies_other, id)  * dt;
        let k4_posit = (body.vel + k3_v) * dt;

        // Step 2: Update position and velocity using weighted average of k-values
        body.vel += (k1_v + k2_v * 2. + k3_v * 2. + k4_v) / 6.;
        body.posit += (k1_posit + k2_posit * 2. + k3_posit * 2. + k4_posit) / 6.;
    }
}

// todo: DRY
// todo: Unless we want to apply accel etc to these rays, RK4 here is not required; accel
// todo is always 0, and v is always c.
pub fn integrate_rk4_ray(rays: &mut [GravRay], dt: f64) {
    // todo: Pending further exporation, no grav accel on rays.
    let a = Vec3::new_zero();

    for body in rays.iter_mut() {
        // Step 1: Calculate the k-values for position and velocity
        let k1_v = a;
        let k1_posit = body.vel * dt;

        let k2_v = a;
        let k2_posit = (body.vel + k1_v * 0.5) * dt;

        let k3_v = a;
        let k3_posit = (body.vel + k2_v * 0.5) * dt;

        let k4_v = a;
        let k4_posit = (body.vel + k3_v) * dt;

        // Step 2: Update position and velocity using weighted average of k-values
        // body.vel += (k1_v + k2_v * 2. + k3_v * 2. + k4_v) / 6.;

        // Vel is constant; C.
        body.posit += (k1_posit + k2_posit * 2. + k3_posit * 2. + k4_posit) / 6.;
    }
}
