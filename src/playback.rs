//! Code related to the playback of computed snapshots.
//!

use graphics::Entity;
use lin_alg::f32::{Quaternion, Vec3 as Vec3f32};

use crate::render::{BODY_COLOR, BODY_SHINYNESS, BODY_SIZE, RAY_COLORS, RAY_SHINYNESS, RAY_SIZE};

#[derive(Debug)]
pub struct SnapShot {
    pub time: usize,
    // To save memory, we store the snapshots as f32; we only need f64 precision
    // during the integration.
    pub body_posits: Vec<Vec3f32>,
    pub V_at_bodies: Vec<Vec3f32>,
    // todo: Determine if you want to store and show these.
    // todo: Store a posit and a velocity for rays A/R.
    // The usize is body id.
    pub rays: Vec<(Vec3f32, usize)>,
}

pub fn change_snapshot(entities: &mut Vec<Entity>, snapshot: &SnapShot) {
    *entities = Vec::new();

    for body_posit in &snapshot.body_posits {
        entities.push(Entity::new(
            0,
            *body_posit,
            Quaternion::new_identity(),
            BODY_SIZE,
            BODY_COLOR,
            BODY_SHINYNESS,
        ));
    }

    for (ray_posit, body_id) in &snapshot.rays {
        entities.push(Entity::new(
            0,
            *ray_posit,
            Quaternion::new_identity(),
            RAY_SIZE,
            RAY_COLORS[body_id % RAY_COLORS.len()],
            RAY_SHINYNESS,
        ));
    }
}
