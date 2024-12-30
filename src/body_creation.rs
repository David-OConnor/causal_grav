//! This model creates distributions of bodies, e.g. ones that coarsely represent galaxies.

use std::f64::consts::TAU;

use lin_alg::f64::Vec3;
use rand::Rng;

use crate::Body;

/// Create n bodies, in circular orbit, at equal distances from each other.
pub fn make_bodies_balanced(num: usize, r: f64, mass_body: f64, mass_central: f64) -> Vec<Body> {
    let mut result = Vec::with_capacity(num);

    for i in 0..num {
        let theta = TAU / num as f64 * i as f64;
        let posit = Vec3::new(r * theta.cos(), r * theta.sin(), 0.0);

        // Velocity magnitude for circular orbit
        let v_mag = (mass_central / r).sqrt();

        // Velocity direction: perpendicular to the radius vector
        let v_x = -v_mag * theta.sin(); // Tangential velocity in x-direction
        let v_y = v_mag * theta.cos(); // Tangential velocity in y-direction

        let vel = Vec3::new(v_x, v_y, 0.0);

        result.push(Body {
            posit,
            vel,
            accel: Vec3::new_zero(),
            mass: mass_body,
        });
    }

    result
}

pub fn make_galaxy_coarse(num_bands: usize, bodies_per_band: usize) -> Vec<Body> {
    let mass_central = 1_000.;

    let mut result = vec![Body {
        posit: Vec3::new_zero(),
        vel: Vec3::new_zero(),
        accel: Vec3::new_zero(),
        mass: mass_central,
    }];

    // todo: Pass as params etc.

    let dist_spacing = 2.5;
    let mass = 6.;

    for i in 1..num_bands {
        let r = dist_spacing * i as f64;

        // if r < 4. {
        //     continue; // todo temp
        // }

        // todo: the central mass method here isn't correct, as it neglects the other masses.
        // todo: (Used to circularize orbits)

        let mut mass_central_adj = mass_central + (i - 1) as f64 * mass * bodies_per_band as f64
            - (num_bands - (i - 1)) as f64 * mass * bodies_per_band as f64;

        // mass_central_adj *= 2.0;

        result.extend_from_slice(&make_bodies_balanced(
            bodies_per_band,
            r,
            mass,
            mass_central_adj,
        ));

        // result
        //     .extend_from_slice(&make_bodies_balanced(3, 5., 10., mass_central));
        // result
        //     .extend_from_slice(&make_bodies_balanced(3, 10., 10., mass_central));
    }
    result
}

// todo: Unused
fn _make_bodies(
    n: usize,
    posit_max_dist: f64,
    mass_range: (f64, f64),
    ω_range: (f64, f64),
) -> Vec<Body> {
    let mut rng = rand::thread_rng();
    let mut result = Vec::with_capacity(n);

    for _ in 0..n {
        // Generate a random position within the maximum distance
        // todo: QC this posit gen.
        let r = rng.gen_range(0.0..posit_max_dist);
        let theta = rng.gen_range(0.0..TAU);
        let posit = Vec3::new(r * theta.cos(), r * theta.sin(), 0.0);

        // Generate a random mass within the range
        let mass = rng.gen_range(mass_range.0..mass_range.1);

        // Generate a velocity with an angular velocity in the specified range
        let ω = rng.gen_range(ω_range.0..ω_range.1);
        let vel = Vec3::new(-ω * posit.y, ω * posit.x, 0.0);

        // Create the body and add it to the result vector
        result.push(Body {
            posit,
            vel,
            accel: Vec3::new_zero(),
            mass,
        });
    }

    result
}
