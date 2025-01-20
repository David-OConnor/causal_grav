use lin_alg::f64::Vec3;

use crate::Body;

/// Compute acceleration, position, and velocity, using RK4.
pub fn integrate_rk4<F>(body: &mut Body, id: usize, acc: &F, dt: f64)
where
    F: Fn(usize, Vec3) -> Vec3,
{
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
}
