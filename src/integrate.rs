use lin_alg::f64::Vec3;

use crate::Body;

/// Compute acceleration, position, and velocity, using RK4.
/// The acc fn: (id, target posit, target charge) -> Acceleration.
/// Target charge is N/A for gravity mode.
pub fn integrate_rk4<F>(body_tgt: &mut Body, id_tgt: usize, acc: &F, dt: f64)
where
    F: Fn(usize, Vec3, f64) -> Vec3,
{
    // Step 1: Calculate the k-values for position and velocity
    body_tgt.accel = acc(id_tgt, body_tgt.posit, body_tgt.mass);

    let k1_v = body_tgt.accel * dt;
    let k1_pos = body_tgt.vel * dt;

    let body_pos_k2 = body_tgt.posit + k1_pos * 0.5;
    let k2_v = acc(id_tgt, body_pos_k2, body_tgt.mass) * dt;
    let k2_pos = (body_tgt.vel + k1_v * 0.5) * dt;

    let body_pos_k3 = body_tgt.posit + k2_pos * 0.5;
    let k3_v = acc(id_tgt, body_pos_k3, body_tgt.mass) * dt;
    let k3_pos = (body_tgt.vel + k2_v * 0.5) * dt;

    let body_pos_k4 = body_tgt.posit + k3_pos;
    let k4_v = acc(id_tgt, body_pos_k4, body_tgt.mass) * dt;
    let k4_pos = (body_tgt.vel + k3_v) * dt;

    // Step 2: Update position and velocity using weighted average of k-values
    body_tgt.vel += (k1_v + k2_v * 2. + k3_v * 2. + k4_v) / 6.;
    body_tgt.posit += (k1_pos + k2_pos * 2. + k3_pos * 2. + k4_pos) / 6.;
}
