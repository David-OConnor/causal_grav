//! Code related to the playback of computed snapshots.
//!
use lin_alg::f32::Vec3 as Vec3f32;

#[derive(Debug)]
pub struct SnapShot {
    pub time: f64,
    // To save memory, we store the snapshots as f32; we only need f64 precision
    // during the integration.
    pub body_posits: Vec<Vec3f32>,
    // todo: Determine if you want to store and show these.
    pub rays: Vec<Vec3f32>,
}
