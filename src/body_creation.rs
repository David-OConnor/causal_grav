//! This model creates distributions of bodies, e.g. ones that coarsely represent galaxies.

use std::f64::consts::TAU;

use lin_alg::{f64::Vec3, linspace};
use rand::Rng;

use crate::{util::interpolate, Body, C};

/// Create n bodies, in circular orbit, at equal distances from each other.
pub fn make_bodies_balanced(num: usize, r: f64, mass_body: f64, mass_central: f64) -> Vec<Body> {
    let mut result = Vec::with_capacity(num);

    for i in 0..num {
        let θ = TAU / num as f64 * i as f64;
        let posit = Vec3::new(r * θ.cos(), r * θ.sin(), 0.0);

        // Velocity magnitude for circular orbit
        let v_mag = (mass_central / r).sqrt();

        // Velocity direction: perpendicular to the radius vector
        let v_x = -v_mag * θ.sin(); // Tangential velocity in x-direction
        let v_y = v_mag * θ.cos(); // Tangential velocity in y-direction

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

/// Make a halo of dark matter, or a galaxy's central bulge
pub fn make_halo_bulge(radius: f64, n_bodies: usize, mass: f64) -> Vec<Body> {
    let mut result = Vec::with_capacity(n_bodies);
    let mut rng = rand::thread_rng();

    for _ in 0..n_bodies {
        let r = radius * rng.gen::<f64>().cbrt(); // Random radius scaled within [0, distribution_radius]
        let θ = rng.gen_range(0.0..TAU); // Random angle θ in [0, 2*pi]
        let ϕ = rng.gen_range(0.0..TAU / 2.); // Random angle phi in [0, pi]

        // Convert spherical coordinates to Cartesian coordinates
        let x = r * ϕ.sin() * θ.cos();
        let y = r * ϕ.sin() * θ.sin();
        let z = r * ϕ.cos();

        result.push(Body {
            posit: Vec3::new(x, y, z),
            vel: Vec3::new_zero(), // todo A/R
            accel: Vec3::new_zero(),
            mass,
        });
    }

    result
}

#[derive(Clone, Copy, PartialEq)]
pub enum GalaxyShape {
    GrandDesignSpiral,
    FlocculentSpiral,
    MultiArmSpiral,
    BarredSpiral,
    Lenticular,
    Elliptical,
}

/// todo: We assume a spiral galaxy for now
pub struct GalaxyDescrip {
    pub shape: GalaxyShape,
    /// See `properties` for what these units are
    pub mass_density: Vec<(f64, f64)>,
    /// r (kpc), v/c
    pub rotation_curve: Vec<(f64, f64)>,
    /// alpha (arcsec), mu (mac arcsec^-2)
    pub luminosity: Vec<(f64, f64)>,
    // todo: More A/R
    /// 0 means a circle. 1 is fully elongated.
    pub eccentricity: f64,
    pub arm_count: usize,
    /// Not a fundamental property; used to normalize mass density etc?
    /// todo: I'm not sure what this is
    pub r_s: f64,
}

fn ring_area(r: f64, dr: f64) -> f64 {
    let r_outer = r + dr / 2.;
    let r_inner = r - dr / 2.;
    let area_outer = r_outer.powi(2) * TAU / 2.;
    let area_inner = r_inner.powi(2) * TAU / 2.;

    area_outer - area_inner
}

impl GalaxyDescrip {
    /// See the `properties` module for info on distributions
    /// todo: Luminosity A/R
    pub fn make_bodies(&self) -> Vec<Body> {
        let mut result = Vec::with_capacity(69); // todo
        let mut rng = rand::thread_rng();

        // todo: Event along rings instead of random?

        // Note: Our distributions tend to be heavily biased towards low r, so if we extend
        // all the way to the end, we will likely leave out lots of values there.

        let r_last = self.mass_density.last().unwrap().0;

        let num_rings = 12;
        let dr = r_last / num_rings as f64;

        let r_all = linspace(0., r_last, num_rings);

        // todo: Maybe dynamically vary this.
        let mass_per_body = 6.;

        // todo: Split into disc + bulge, among other things.
        for r in r_all {
            // todo: Most of this is a hack.

            // This is ρ/ρ_0.
            let rho = interpolate(&self.mass_density, r).unwrap();

            // Calculate
            let disc_thickness = 30.; // todo: Important question. Disc thickness?
            let ring_volume = ring_area(r, dr) * disc_thickness;
            // let circum = TAU * r;

            // todo: Normalize?
            // rho is (normalized) M/L^3. So, M = rho * L^3
            let bodies_this_r = (1. * rho * ring_volume / mass_per_body) as usize;

            println!("N bodies: {:?}", bodies_this_r);

            // Multiply by C, because the curve is normalized to C.
            // todo: The fudge factor...
            let v_mag = 1_500. * C * interpolate(&self.rotation_curve, r).unwrap();

            for i in 0..bodies_this_r {
                // todo: Random, or even? Even is more uniform, which may be nice, but
                // todo it may cause resonances. Maybe even, but with a random offset per ring?
                let θ = rng.gen_range(0.0..TAU);
                // let θ = TAU / bodies_this_r as f64 * i as f64;

                // Apply eccentricity: Scale radius in x-direction
                let scale_x = 1.0 - self.eccentricity; // Eccentricity factor for x-axis

                let posit = Vec3::new(r * θ.cos() * scale_x, r * θ.sin(), 0.0);

                // todo: Does v need to be a function of theta due to eccentricity?

                // Velocity direction: perpendicular to the radius vector
                let v_x = -v_mag * θ.sin(); // Tangential velocity in x-direction
                let v_y = v_mag * θ.cos(); // Tangential velocity in y-direction

                let vel = Vec3::new(v_x, v_y, 0.0);

                result.push(Body {
                    posit,
                    vel,
                    accel: Vec3::new_zero(),
                    mass: mass_per_body,
                })
            }
        }

        result
    }
}

/// todo: Move specific galaxy creation to its own module A/R
#[derive(Clone, Copy, PartialEq)]
pub enum GalaxyModel {
    Ngc1560,
    Ngc3198,
    Ngc3115,
    Ngc3031,
    Ngc7331,
}

impl GalaxyModel {
    pub fn descrip(&self) -> GalaxyDescrip {
        match self {
            /// Ludwig, Figures 3 and 5. todo: Partial/rough
            Self::Ngc1560 => GalaxyDescrip {
                shape: GalaxyShape::FlocculentSpiral, // todo ?
                mass_density: vec![
                    (0.01, 1.),
                    (0.02, 1.),
                    (0.05, 1.),
                    (0.10, 1.),
                    (0.8, 0.8),
                    (1.0, 0.61),
                    (3.0, 0.05),
                    (10.0, 0.0),
                ],
                rotation_curve: vec![
                    (0., 0.),
                    (1., 0.00009),
                    (2., 0.00013),
                    (3., 0.00017),
                    (4., 0.00020),
                    (5., 0.00022),
                    (6., 0.00024),
                    (7., 0.000245),
                    (8., 0.00025),
                    (9., 0.000251),
                    (10., 0.000251),
                    (11., 0.000252),
                    (12., 0.000252),
                ],
                luminosity: vec![],
                eccentricity: 0.0, // todo temp
                arm_count: 2,
                r_s: 1.46e-6,
            },
            Self::Ngc3198 => GalaxyDescrip {
                shape: GalaxyShape::BarredSpiral,
                mass_density: vec![
                    (0.1, 1.),
                    (0.2, 0.95),
                    (0.5, 0.6),
                    (1.0, 0.38),
                    (3.0, 0.1),
                    (5.0, 0.07),
                    (7.0, 0.03),
                    (10.0, 0.008),
                ],
                rotation_curve: vec![
                    (0., 0.),
                    (2., 0.00035),
                    (4., 0.00045),
                    (6., 0.00048),
                    (8., 0.00051),
                    (10., 0.00050),
                    (12., 0.00049),
                    (14., 0.00049),
                    (16., 0.00049),
                    (18., 0.00048),
                    (20., 0.00048),
                    (22., 0.00048),
                    (24., 0.00047),
                    (26., 0.00048),
                    (28., 0.00049),
                    (30., 0.00049),
                ],
                luminosity: vec![],
                eccentricity: 0.,
                arm_count: 2,
                r_s: 1.2e-5,
            },
            Self::Ngc3115 => GalaxyDescrip {
                shape: GalaxyShape::Lenticular,
                mass_density: vec![],
                rotation_curve: vec![],
                luminosity: vec![],
                eccentricity: 0.,
                arm_count: 2,
                r_s: 6.97e-16,
            },
            _ => unimplemented!(), // todo
        }
    }

    pub fn make_bodies(&self) -> Vec<Body> {
        self.descrip().make_bodies()
    }
}
