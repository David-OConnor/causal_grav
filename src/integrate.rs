use lin_alg::f64::Vec3;

use crate::{accel, Body, GravRay, GravShell};

pub fn integrate_rk4(
    bodies: &mut Vec<Body>,
    shells: &[GravShell],
    dt: f64,
    gauss_c: f64,
    acc_inst: bool,
) {
    let bodies_other = bodies.clone(); // todo: I don't like this. Avoids mut error.

    let acc = |id, posit| {
        if acc_inst {
            accel::calc_acc_inst(posit, &bodies_other, id)
        } else {
            accel::calc_acc_shell(shells, posit, id, gauss_c)
        }
    };

    for (id, body) in bodies.iter_mut().enumerate() {
        // Step 1: Calculate the k-values for position and velocity
        // body.accel = acc(id, body.posit); // todo: Temp placed before fn.

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
