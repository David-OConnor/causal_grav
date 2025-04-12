use std::{
    f64::consts::TAU,
    fs::File,
    io,
    io::{ErrorKind, Read, Write},
    path::Path,
};

use bincode::{config, Decode, Encode};
use lin_alg::f64::Vec3;
use rand::{rngs::ThreadRng, Rng};

use crate::{Body, State};

/// This function generates an interpolated value for the given `val` based on the
/// provided `data`. The `data` is a set of (x, y) pairs, where `x` is the input
/// and `y` is the corresponding output value.
pub fn interpolate(data: &[(f64, f64)], val: f64) -> Option<f64> {
    // Ensure there are at least two points for interpolation
    if data.len() < 2 {
        panic!("At least two data points are required for interpolation.");
    }

    // If `val` is outside the range of data, extrapolate
    if val < data[0].0 {
        let (x0, y0) = data[0];
        let (x1, y1) = data[1];
        return Some(y0 + (val - x0) / (x1 - x0) * (y1 - y0));
    } else if val > data[data.len() - 1].0 {
        let (x0, y0) = data[data.len() - 2];
        let (x1, y1) = data[data.len() - 1];
        return Some(y0 + (val - x0) / (x1 - x0) * (y1 - y0));
    }

    // Find the interval in which `val` falls
    for i in 0..data.len() - 1 {
        let (x0, y0) = data[i];
        let (x1, y1) = data[i + 1];

        if val >= x0 && val <= x1 {
            // if (x1 - x0).abs() < 1e-10 {
            //     return
            // }
            // Perform linear interpolation
            let t = (val - x0) / (x1 - x0);
            return Some(y0 + t * (y1 - y0));
        }
    }

    None
}

/// Defaults to `Config::dt_integration`, but becomes more precise when
/// bodies are close. This is a global DT, vice local only for those bodies.
#[allow(unused)]
pub fn calc_dt_dynamic(state: &State, bodies_other: &[Body]) -> f64 {
    let mut result = state.config.dt_integration_max;

    // todo: Consider cacheing the distances, so this second iteration can be reused.
    for (id_tgt, body_tgt) in state.bodies.iter().enumerate() {
        for (i_src, body_src) in bodies_other.iter().enumerate() {
            if i_src == id_tgt {
                continue; // self-interaction.
            }

            let dist = (body_src.posit - body_tgt.posit).magnitude();
            let rel_velocity = (body_src.vel - body_tgt.vel).magnitude();
            let dt = state.config.dynamic_dt_scaler * dist / rel_velocity;
            if dt < result {
                result = dt;
            }
        }
    }

    result
}

/// E.g. converting arcseconds to kpc, for galaxy radius.
pub fn scale_x_axis(data: &[(f64, f64)], scaler: f64) -> Vec<(f64, f64)> {
    data.iter().map(|(x, y)| (scaler * x, *y)).collect()
}

/// Combine separate radius and mass density or velocity data.
pub fn zip_data(r: &[f64], vals: Vec<f64>) -> Vec<(f64, f64)> {
    r.iter().copied().zip(vals).collect()
}

/// Save to file, using Bincode. We currently use this for preference files.
pub fn save<T: Encode>(path: &Path, data: &T) -> io::Result<()> {
    let config = config::standard();

    let encoded: Vec<u8> = bincode::encode_to_vec(data, config).unwrap();

    let mut file = File::create(path)?;
    file.write_all(&encoded)?;
    Ok(())
}

/// Load from file, using Bincode. We currently use this for preference files.
pub fn load<T: Decode<()>>(path: &Path) -> io::Result<T> {
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

pub fn volume_sphere(r: f64) -> f64 {
    const COEFF: f64 = 2. / 3. * TAU;
    r.powi(3) * COEFF
}

pub fn random_unit_vec(rng: &mut ThreadRng) -> Vec3 {
    let θ = rng.random_range(0.0..TAU);

    // Random phi for polar angle with area weighting
    let u: f64 = rng.random_range(-1.0..1.0); // Uniform random variable
    let ϕ = u.acos(); // Inverse cosine for area-preserving sampling

    // Convert spherical coordinates to Cartesian coordinates
    let x = ϕ.sin() * θ.cos();
    let y = ϕ.sin() * θ.sin();
    let z = ϕ.cos();

    Vec3::new(x, y, z)
}
