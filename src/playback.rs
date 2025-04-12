//! Code related to the playback of computed snapshots.
//!

use std::{
    fs::File,
    io,
    io::{ErrorKind, Read, Write},
    path::Path,
};

use barnes_hut::{Cube, Node};
use bincode::{
    config,
    error::{DecodeError, EncodeError},
    Decode, Encode,
};
use graphics::{Entity, UP_VEC};
use lin_alg::{
    f32::{Quaternion, Vec3 as Vec3f32},
    f64::Vec3,
};

use crate::{
    grav_shell::GravShell,
    render::{
        ARROW_COLOR, ARROW_SHINYNESS, BODY_COLOR, BODY_SHINYNESS, BODY_SIZE_MAX, BODY_SIZE_MIN,
        BODY_SIZE_SCALER, MESH_ARROW, MESH_CUBE, MESH_SPHERE, SHELL_COLOR, TREE_COLOR,
        TREE_CUBE_SCALE_FACTOR, TREE_SHINYNESS,
    },
};

#[derive(Debug, Encode, Decode)]
/// A compact version
pub struct GravShellSnapshot {
    center: Vec3f32,
    radius: f32,
    src_mass: f32,
}

impl GravShellSnapshot {
    pub fn new(shell: &GravShell) -> Self {
        Self {
            center: shell.center.into(),
            radius: shell.radius as f32,
            src_mass: shell.src_mass as f32,
        }
    }
}

#[derive(Debug, Encode, Decode, Default)]
pub struct SnapShot {
    pub time: f32,
    // To save memory, we store the snapshots as f32; we only need f64 precision
    // during the integration.
    pub body_posits: Vec<Vec3f32>,
    // pub V_at_bodies: Vec<Vec3f32>,
    pub body_accs: Vec<Vec3f32>,
    // todo: Velocity?
    // pub body_vels: Vec<Vec3f32>,
    // todo: Determine if you want to store and show these.
    // todo: Store a posit and a velocity for rays A/R.
    // The usize is body id.
    // pub rays: Vec<(Vec3f32, usize)>,
    // todo: Compact form for shells, as above?
    pub shells: Vec<GravShellSnapshot>,
    pub dt: f32,
    pub tree_cubes: Vec<Cube>, // todo: Custom type type f32, as above.
}

/// Body masses are separate from the snapshot, since it's invariant.
pub fn change_snapshot(entities: &mut Vec<Entity>, snapshot: &SnapShot, body_masses: &[f32]) {
    // todo: Shells, acc vecs A/R
    *entities = Vec::with_capacity(snapshot.body_posits.len() + snapshot.tree_cubes.len());

    for (i, posit) in snapshot.body_posits.iter().enumerate() {
        let entity_size = f32::clamp(
            BODY_SIZE_SCALER * body_masses[i],
            BODY_SIZE_MIN,
            BODY_SIZE_MAX,
        );
        entities.push(Entity::new(
            MESH_SPHERE,
            *posit,
            Quaternion::new_identity(),
            entity_size,
            BODY_COLOR,
            BODY_SHINYNESS,
        ));

        // entities.push(Entity::new(
        //     MESH_ARROW,
        //     *posit,
        //     Quaternion::from_unit_vecs(UP_VEC, snapshot.body_accs[i].to_normalized()),
        //     snapshot.body_accs[i].magnitude() * 0.2,
        //     ARROW_COLOR,
        //     ARROW_SHINYNESS,
        // ));
    }

    for cube in &snapshot.tree_cubes {
        entities.push(Entity::new(
            MESH_CUBE,
            Vec3f32::new(
                cube.center.x as f32,
                cube.center.y as f32,
                cube.center.z as f32,
            ),
            Quaternion::new_identity(),
            cube.width as f32 * TREE_CUBE_SCALE_FACTOR,
            TREE_COLOR,
            TREE_SHINYNESS,
        ));
    }

    // for (ray_posit, body_id) in &snapshot.rays {
    //     entities.push(Entity::new(
    //         0,
    //         *ray_posit,
    //         Quaternion::new_identity(),
    //         RAY_SIZE,
    //         RAY_COLORS[body_id % RAY_COLORS.len()],
    //         RAY_SHINYNESS,
    //     ));
    // }

    // todo: Draw an actual shell instead of a sphere.
    // todo: Add back once you sort out transparency.
    for shell in &snapshot.shells {
        // let center = Vec3f32::new(shell.center.x, shell.center.y, shell.center.z);

        // let entity = Entity::new(
        //     2,
        //     center,
        //     Quaternion::new_identity(), // todo: A/R?
        //     shell.radius,
        //     // todo: Custom colors A/R
        //     // todo: Color code by radius.
        //     SHELL_COLOR, // todo: Color-code A/R
        //     BODY_SHINYNESS,
        // );
        //
        // // entity.opacity = SHELL_OPACITY;
        // entities.push(entity);
    }
}
