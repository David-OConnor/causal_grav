use crate::{accel, Body, ForceModel, GravShell};

pub fn integrate_rk4(
    bodies: &mut [Body],
    shells: &[GravShell],
    dt: f64,
    gauss_c: f64,
    force_model: ForceModel,
) {
    let bodies_other = bodies.to_owned(); // todo: I don't like this. Avoids mut error.

    let acc = |id, posit| match force_model {
        ForceModel::Newton => accel::acc_newton(posit, &bodies_other, id, None),
        ForceModel::GaussShells => accel::calc_acc_shell(shells, posit, id, gauss_c),
        ForceModel::Mond(mond_fn) => accel::acc_newton(posit, &bodies_other, id, Some(mond_fn)),
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
