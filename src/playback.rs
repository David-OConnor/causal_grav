//! Code related to the playback of computed snapshots.
//!

use std::{
    fs::File,
    io,
    io::{ErrorKind, Read, Write},
    path::Path,
};

use bincode::{
    config,
    error::{DecodeError, EncodeError},
    Decode, Encode,
};
use graphics::Entity;
use lin_alg::{
    f32::{Quaternion, Vec3 as Vec3f32_b},
    f64::Vec3,
};

use crate::{
    render::{
        BODY_COLOR, BODY_SHINYNESS, BODY_SIZE_MAX, BODY_SIZE_MIN, BODY_SIZE_SCALER, SHELL_COLOR,
    },
    GravShell,
};

pub const DEFAULT_SNAPSHOT_FILE: &str = "snapshot.cg";

/// We use a custom type, vice from lin_alg, so we can impl Encode and Decode.
#[derive(Clone, Copy, Debug, Encode, Decode)]
pub struct Vec3f32 {
    pub x: f32,
    pub y: f32,
    pub z: f32,
}

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
            center: vec3_to_f32(shell.center),
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
    // todo: Determine if you want to store and show these.
    // todo: Store a posit and a velocity for rays A/R.
    // The usize is body id.
    // pub rays: Vec<(Vec3f32, usize)>,
    // todo: Compact form for shells, as above?
    pub shells: Vec<GravShellSnapshot>,
    pub dt: f32,
}

pub fn vec3_to_f32(v: Vec3) -> Vec3f32 {
    Vec3f32 {
        x: v.x as f32,
        y: v.y as f32,
        z: v.z as f32,
    }
}

/// Body masses are separate from the snapshot, since it's invariant.
pub fn change_snapshot(entities: &mut Vec<Entity>, snapshot: &SnapShot, body_masses: &[f32]) {
    // *entities = Vec::with_capacity(entities.len() + snapshot.rays.len());
    *entities = Vec::with_capacity(entities.len());

    for (i, body_posit_) in snapshot.body_posits.iter().enumerate() {
        let body_posit = Vec3f32_b::new(body_posit_.x, body_posit_.y, body_posit_.z);
        // println!("Body mass: {:?}. Scaled size: {:?}", body_masses[i], BODY_SIZE_SCALER * body_masses[i]);

        let entity_size = f32::clamp(
            BODY_SIZE_SCALER * body_masses[i],
            BODY_SIZE_MIN,
            BODY_SIZE_MAX,
        );
        entities.push(Entity::new(
            0,
            body_posit,
            Quaternion::new_identity(),
            entity_size,
            BODY_COLOR,
            BODY_SHINYNESS,
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
        // let center = Vec3f32_b::new(shell.center.x, shell.center.y, shell.center.z);

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

/// Save to file, using Bincode.
pub fn save<T: Encode>(path: &Path, data: &T) -> io::Result<()> {
    let config = config::standard();

    let encoded: Vec<u8> = bincode::encode_to_vec(data, config).unwrap();

    let mut file = File::create(path)?;
    file.write_all(&encoded)?;
    Ok(())
}

/// Load from file, using Bincode.
pub fn load<T: Decode>(path: &Path) -> io::Result<T> {
    let config = config::standard();

    let mut file = File::open(path)?;
    let mut buffer = Vec::new();
    file.read_to_end(&mut buffer)?;
    let (decoded, _len) = match bincode::decode_from_slice(&buffer, config) {
        Ok(v) => v,
        Err(_) => {
            eprintln!("Error loading from file. Did the format change?");
            return Err(io::Error::new(ErrorKind::Other, "error loading"));
        }
    };
    Ok(decoded)
}
