//! Code related to the playback of computed snapshots.
//!

use graphics::Entity;
use lin_alg::{
    f32::{Quaternion, Vec3 as Vec3f32},
    f64::Vec3,
};

use crate::{
    render::{
        BODY_COLOR, BODY_SHINYNESS, BODY_SIZE, RAY_COLORS, RAY_SHINYNESS, RAY_SIZE, SHELL_OPACITY,
    },
    GravShell,
};

#[derive(Debug)]
pub struct SnapShot {
    pub time: usize,
    // To save memory, we store the snapshots as f32; we only need f64 precision
    // during the integration.
    pub body_posits: Vec<Vec3f32>,
    pub V_at_bodies: Vec<Vec3f32>,
    pub body_accs: Vec<Vec3f32>,
    // todo: Determine if you want to store and show these.
    // todo: Store a posit and a velocity for rays A/R.
    // The usize is body id.
    pub rays: Vec<(Vec3f32, usize)>,
    // todo: Compact form for shells, as above?
    pub shells: Vec<GravShell>,
}

pub fn vec_to_f32(v: Vec3) -> Vec3f32 {
    Vec3f32::new(v.x as f32, v.y as f32, v.z as f32)
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

    // todo: Draw an actual shell instead of a sphere.
    for shell in &snapshot.shells {
        let mut entity = Entity::new(
            1,
            vec_to_f32(shell.center),
            Quaternion::new_identity(),
            shell.radius as f32,
            RAY_COLORS[shell.emitter_id % RAY_COLORS.len()],
            RAY_SHINYNESS,
        );
        // entity.opacity = SHELL_OPACITY;
        entity.opacity = 0.; // todo temp
        entities.push(entity);
    }
}
