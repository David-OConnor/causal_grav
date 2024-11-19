//! Code related to the playback of computed snapshots.
//! 
use lin_alg::f32::Vec3 as Vec3f32;


#[derive(Debug)]
struct SnapShot {
    pub time: f64,
    // todo: To save memory, you could store the snapshots as f32; we only need f64 precision
    // todo during the integration.
    pub body_posits: Vec<Vec3f32>,
}