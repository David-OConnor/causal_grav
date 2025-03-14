//! Hijacking this programs' infrastructure to do some electron modelling.

use std::{f64::consts::TAU, fmt, fmt::Formatter};

use lin_alg::f64::{Quaternion, Vec3};
use rand::Rng;

use crate::Body;

pub fn force_coulomb(
    acc_dir: Vec3,
    q_src: f64,
    q_tgt: f64,
    dist: f64,
    softening_factor_sq: f64,
) -> Vec3 {
    acc_dir * q_src * q_tgt / (dist.powi(2) + softening_factor_sq)
}

pub fn make_particles() -> Vec<Body> {
    // todo: Maybe don't make even at R; distribute spacially uniformly.
    let n_particles = 10_000;
    let mut rng = rand::rng();

    let mut result = Vec::with_capacity(n_particles);
    // let mut result = Vec::with_capacity(n_particles + 1);

    // Nucleus
    // result.push(Body {
    //     posit: Vec3::new_zero(),
    //     vel: Vec3::new_zero(),
    //     accel: Vec3::new_zero(),
    //     mass: 1.,
    // });

    // Implicit nuc for now.

    // todo: You may need to continuously generate particles, and use something
    // todo like flux as your psi equivalent.

    // todo: Re-use algos you have for this random dist.
    for _ in 0..n_particles {
        let r = rng.random_range(0.0..=20.);
        let θ = rng.random_range(0.0..TAU);

        let ϕ = {
            // Random phi for polar angle with area weighting
            let u: f64 = rng.random_range(-1.0..1.0); // Uniform random variable
            u.acos() // Inverse cosine for area-preserving sampling
        };

        // Convert spherical coordinates to Cartesian coordinates
        let x = r * ϕ.sin() * θ.cos();
        let y = r * ϕ.sin() * θ.sin();
        let z = r * ϕ.cos();

        let posit = Vec3::new(x, y, z);

        // Add a slight tangental velocity to avoid a singularity in the center.
        // This is an arbitraryish direction (?)
        let w = Vec3::new(1., 0., 0.); // arbitrary
        let tangent_vec = w.cross(posit).to_normalized();
        let rotator =
            Quaternion::from_axis_angle(posit.to_normalized(), rng.random_range(0.0..TAU));

        // Hmm... KE = 1/2 mv^2 classically
        // set v = sqrt(KE * 2/m) ?
        // If KE = -1/2 and m = 1: v = sqrt(-1) = i... or say KE = 1. v = 1 for n=1 Hydrogen?

        // let speed = rng.random_range(0.0..1.);
        let speed = 1.;
        let vel = rotator.rotate_vec(tangent_vec) * speed;

        // todo: "intertial" mass vs charge...

        result.push(Body {
            posit,
            vel,
            accel: Vec3::new_zero(),
            mass: 1.,
        });
    }

    result
}

#[derive(Debug)]
pub struct FieldProperties {
    pub avg_vel: Vec3,
    pub density: f64,
    pub flux: Vec3,
    pub divergence: f64,
    pub curl: Vec3,
    pub avg_accel: Vec3,
    pub accel_divergence: f64,
    pub accel_curl: Vec3,
}

impl fmt::Display for FieldProperties {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        writeln!(f, "V avg: {}", self.avg_vel)?;
        writeln!(f, "ρ: {}", self.density)?;
        writeln!(f, "Flux: {}", self.flux)?;
        writeln!(f, "∇·F: {:.3}", self.divergence)?;
        writeln!(f, "∇×F: {}", self.curl)?;
        writeln!(f, "Acc avg: {}", self.avg_accel)?;
        writeln!(f, "Acc ∇·F: {:.3}", self.accel_divergence)?;
        writeln!(f, "Acc ∇×F: {}", self.accel_curl)?;

        Ok(())
    }
}

impl FieldProperties {
    pub fn new(bodies: &[Body], center: Vec3, r: f64) -> Self {
        Self {
            avg_vel: average_velocity(bodies, center, r),
            density: get_density(bodies, center, r),
            flux: get_flux(bodies, center, r),
            divergence: get_divergence(bodies, center, r),
            curl: get_curl(bodies, center, r),
            avg_accel: average_acceleration(bodies, center, r),
            accel_divergence: get_accel_divergence(bodies, center, r),
            accel_curl: get_accel_curl(bodies, center, r),
        }
    }
}

/// Return the average velocity among the bodies within a distance `r` of `center`,
/// weighted by their mass. (Naive example.)
fn average_velocity(bodies: &[Body], center: Vec3, r: f64) -> Vec3 {
    let mut total_mass = 0.0;
    let mut momentum_sum = Vec3::new_zero();

    for b in bodies {
        let dist = (b.posit - center).magnitude();
        if dist <= r {
            total_mass += b.mass;
            momentum_sum = momentum_sum + (b.vel * b.mass);
        }
    }

    if total_mass > 0.0 {
        momentum_sum / total_mass
    } else {
        Vec3::new_zero()
    }
}

/// Get a naive local density by summing the mass of all bodies within distance `dx`.
/// Then we divide by the volume of the sphere (4/3 π dx^3).
pub fn get_density(bodies: &[Body], posit: Vec3, dx: f64) -> f64 {
    let mut total_mass = 0.0;
    for b in bodies {
        let dist = (b.posit - posit).magnitude();
        if dist <= dx {
            total_mass += b.mass;
        }
    }
    // For a spherical region:
    let volume = (4.0 / 3.0) * std::f64::consts::PI * dx.powi(3);
    if volume > 0.0 {
        total_mass / volume
    } else {
        0.0
    }
}

/// Get the mass flux (a vector) at a point by multiplying density by the *locally
/// averaged velocity*. Another naive approach.
pub fn get_flux(bodies: &[Body], posit: Vec3, dx: f64) -> Vec3 {
    // local density:
    let density = get_density(bodies, posit, dx);

    // local average velocity:
    let avg_vel = average_velocity(bodies, posit, dx);

    // mass flux = ρ * v  (vector)
    avg_vel * density
}

///  ∇ · F
/// Estimate the divergence of the velocity field at `posit` by finite differences.
///  div(v) = dVx/dx + dVy/dy + dVz/dz
pub fn get_divergence(bodies: &[Body], posit: Vec3, dx: f64) -> f64 {
    // We'll reuse `average_velocity` to sample the velocity at +/- dx around `posit`.

    let v_xp = average_velocity(
        bodies,
        posit
            + Vec3 {
                x: dx,
                y: 0.0,
                z: 0.0,
            },
        dx,
    );
    let v_xm = average_velocity(
        bodies,
        posit
            - Vec3 {
                x: dx,
                y: 0.0,
                z: 0.0,
            },
        dx,
    );

    let v_yp = average_velocity(
        bodies,
        posit
            + Vec3 {
                x: 0.0,
                y: dx,
                z: 0.0,
            },
        dx,
    );
    let v_ym = average_velocity(
        bodies,
        posit
            - Vec3 {
                x: 0.0,
                y: dx,
                z: 0.0,
            },
        dx,
    );

    let v_zp = average_velocity(
        bodies,
        posit
            + Vec3 {
                x: 0.0,
                y: 0.0,
                z: dx,
            },
        dx,
    );
    let v_zm = average_velocity(
        bodies,
        posit
            - Vec3 {
                x: 0.0,
                y: 0.0,
                z: dx,
            },
        dx,
    );

    let dvx_dx = (v_xp.x - v_xm.x) / (2.0 * dx);
    let dvy_dy = (v_yp.y - v_ym.y) / (2.0 * dx);
    let dvz_dz = (v_zp.z - v_zm.z) / (2.0 * dx);

    dvx_dx + dvy_dy + dvz_dz
}

/// Estimate the curl of the velocity field by finite differences of the local
/// average velocity. We sample v at +/- dx along each axis, then compute:
///
///   curl(v) = (dVz/dy - dVy/dz, dVx/dz - dVz/dx, dVy/dx - dVx/dy)
///
/// ∇ × F
pub fn get_curl(bodies: &[Body], posit: Vec3, dx: f64) -> Vec3 {
    // We'll reuse the same `dx` for both the sampling radius and the derivative step.
    // In an actual code base, you might want separate parameters for the "sampling radius"
    // vs. the "delta used in finite differences".
    //
    // Sample velocities at 6 surrounding points:
    let v_xp = average_velocity(
        bodies,
        posit
            + Vec3 {
                x: dx,
                y: 0.0,
                z: 0.0,
            },
        dx,
    );
    let v_xm = average_velocity(
        bodies,
        posit
            - Vec3 {
                x: dx,
                y: 0.0,
                z: 0.0,
            },
        dx,
    );
    let v_yp = average_velocity(
        bodies,
        posit
            + Vec3 {
                x: 0.0,
                y: dx,
                z: 0.0,
            },
        dx,
    );
    let v_ym = average_velocity(
        bodies,
        posit
            - Vec3 {
                x: 0.0,
                y: dx,
                z: 0.0,
            },
        dx,
    );
    let v_zp = average_velocity(
        bodies,
        posit
            + Vec3 {
                x: 0.0,
                y: 0.0,
                z: dx,
            },
        dx,
    );
    let v_zm = average_velocity(
        bodies,
        posit
            - Vec3 {
                x: 0.0,
                y: 0.0,
                z: dx,
            },
        dx,
    );

    // Partial derivatives of v_x, v_y, v_z
    let dvx_dy = (v_yp.x - v_ym.x) / (2.0 * dx);
    let dvx_dz = (v_zp.x - v_zm.x) / (2.0 * dx);
    let dvy_dx = (v_xp.y - v_xm.y) / (2.0 * dx);
    let dvy_dz = (v_zp.y - v_zm.y) / (2.0 * dx);
    let dvz_dx = (v_xp.z - v_xm.z) / (2.0 * dx);
    let dvz_dy = (v_yp.z - v_ym.z) / (2.0 * dx);

    // curl(v) = (dVz/dy - dVy/dz, dVx/dz - dVz/dx, dVy/dx - dVx/dy)
    Vec3 {
        x: dvz_dy - dvy_dz,
        y: dvx_dz - dvz_dx,
        z: dvy_dx - dvx_dy,
    }
}

/// Returns the mass-weighted average acceleration of all bodies within distance `r`.
fn average_acceleration(bodies: &[Body], center: Vec3, r: f64) -> Vec3 {
    let mut total_mass = 0.0;
    let mut accel_sum = Vec3::new_zero();

    for b in bodies {
        let dist = (b.posit - center).magnitude();
        if dist <= r {
            total_mass += b.mass;
            accel_sum = accel_sum + (b.accel * b.mass);
        }
    }

    if total_mass > 0.0 {
        accel_sum / total_mass
    } else {
        Vec3::new_zero()
    }
}

/// Divergence of acceleration field: ∇·a
pub fn get_accel_divergence(bodies: &[Body], posit: Vec3, dx: f64) -> f64 {
    let a_xp = average_acceleration(
        bodies,
        posit
            + Vec3 {
                x: dx,
                y: 0.0,
                z: 0.0,
            },
        dx,
    );
    let a_xm = average_acceleration(
        bodies,
        posit
            - Vec3 {
                x: dx,
                y: 0.0,
                z: 0.0,
            },
        dx,
    );

    let a_yp = average_acceleration(
        bodies,
        posit
            + Vec3 {
                x: 0.0,
                y: dx,
                z: 0.0,
            },
        dx,
    );
    let a_ym = average_acceleration(
        bodies,
        posit
            - Vec3 {
                x: 0.0,
                y: dx,
                z: 0.0,
            },
        dx,
    );

    let a_zp = average_acceleration(
        bodies,
        posit
            + Vec3 {
                x: 0.0,
                y: 0.0,
                z: dx,
            },
        dx,
    );
    let a_zm = average_acceleration(
        bodies,
        posit
            - Vec3 {
                x: 0.0,
                y: 0.0,
                z: dx,
            },
        dx,
    );

    // Partial derivatives via central difference
    let dax_dx = (a_xp.x - a_xm.x) / (2.0 * dx);
    let day_dy = (a_yp.y - a_ym.y) / (2.0 * dx);
    let daz_dz = (a_zp.z - a_zm.z) / (2.0 * dx);

    dax_dx + day_dy + daz_dz
}

/// Curl of acceleration field: ∇×a
pub fn get_accel_curl(bodies: &[Body], posit: Vec3, dx: f64) -> Vec3 {
    let a_xp = average_acceleration(
        bodies,
        posit
            + Vec3 {
                x: dx,
                y: 0.0,
                z: 0.0,
            },
        dx,
    );
    let a_xm = average_acceleration(
        bodies,
        posit
            - Vec3 {
                x: dx,
                y: 0.0,
                z: 0.0,
            },
        dx,
    );

    let a_yp = average_acceleration(
        bodies,
        posit
            + Vec3 {
                x: 0.0,
                y: dx,
                z: 0.0,
            },
        dx,
    );
    let a_ym = average_acceleration(
        bodies,
        posit
            - Vec3 {
                x: 0.0,
                y: dx,
                z: 0.0,
            },
        dx,
    );

    let a_zp = average_acceleration(
        bodies,
        posit
            + Vec3 {
                x: 0.0,
                y: 0.0,
                z: dx,
            },
        dx,
    );
    let a_zm = average_acceleration(
        bodies,
        posit
            - Vec3 {
                x: 0.0,
                y: 0.0,
                z: dx,
            },
        dx,
    );

    let dax_dy = (a_yp.x - a_ym.x) / (2.0 * dx);
    let dax_dz = (a_zp.x - a_zm.x) / (2.0 * dx);

    let day_dx = (a_xp.y - a_xm.y) / (2.0 * dx);
    let day_dz = (a_zp.y - a_zm.y) / (2.0 * dx);

    let daz_dx = (a_xp.z - a_xm.z) / (2.0 * dx);
    let daz_dy = (a_yp.z - a_ym.z) / (2.0 * dx);

    // curl(a) = (dAz/dy - dAy/dz, dAx/dz - dAz/dx, dAy/dx - dAx/dy)
    Vec3 {
        x: daz_dy - day_dz,
        y: dax_dz - daz_dx,
        z: day_dx - dax_dy,
    }
}
