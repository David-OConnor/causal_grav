//! This model creates distributions of bodies, e.g. ones that coarsely represent galaxies.

use std::f64::consts::TAU;

use lin_alg::f64::Vec3;
use rand::Rng;

use crate::{
    util::{interpolate, linspace},
    Body, C,
};

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

pub fn make_galaxy_coarse(
    num_bands: usize,
    bodies_per_band: usize,
    bodies_bulge: usize,
) -> Vec<Body> {
    let mass_central = 1_000.;

    let n_bodies_spiral = (num_bands - 1) * bodies_per_band + 1;
    let mut result = Vec::with_capacity(n_bodies_spiral + bodies_bulge);

    // A central mass
    result.push(Body {
        posit: Vec3::new_zero(),
        vel: Vec3::new_zero(),
        accel: Vec3::new_zero(),
        mass: mass_central,
    });

    let dist_spacing = 2.5;
    let mass = 6.;

    // The bulge
    result.extend(make_halo_bulge(4., bodies_bulge, mass));

    // todo: Pass as params etc.

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
    pub eccentricity: f64,
    pub arm_count: usize,
}

impl GalaxyDescrip {
    /// See the `properties` module for info on distributions
    /// todo: Luminosity A/R
    pub fn make_bodies(&self) -> Vec<Body> {
        let mut result = Vec::with_capacity(69); // todo
        let mut rng = rand::thread_rng();

        // todo: eccentricity.
        // todo: Event along rings instead of random?

        // Note: Our distributions tend to be heavily biased towards low r, so if we extend
        // all the way to the end, we will likely leave out lots of values there.

        let num_rings = 10;
        let rs = linspace(0., self.mass_density.last().unwrap().0, num_rings);

        // todo: Split into disc + bulge, among other things.
        for r in rs {
            // todo: Most of this is a hack.

            let rho = interpolate(&self.mass_density, r).unwrap();
            let circum = TAU * r;
            let bodies_this_r = (10. * rho * circum) as usize;

            // Multiply by C, because the curve is normalized to C.
            // todo: The fudge factor...
            let v_mag = 1_500. * C * interpolate(&self.rotation_curve, r).unwrap();

            for i in 0..bodies_this_r {
                // todo: Random, or even? Even is more uniform, which may be nice, but
                // todo it may cause resonances. Maybe even, but with a random offset per ring?
                let θ = rng.gen_range(0.0..TAU);
                // let θ = TAU / bodies_this_r as f64 * i as f64;

                let posit = Vec3::new(r * θ.cos(), r * θ.sin(), 0.0);

                // Velocity direction: perpendicular to the radius vector
                let v_x = -v_mag * θ.sin(); // Tangential velocity in x-direction
                let v_y = v_mag * θ.cos(); // Tangential velocity in y-direction

                let vel = Vec3::new(v_x, v_y, 0.0);

                // todo?
                let mass = 6.;

                result.push(Body {
                    posit,
                    vel,
                    accel: Vec3::new_zero(),
                    mass,
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
                eccentricity: 0.,
                arm_count: 2,
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
            },
            Self::Ngc3115 => GalaxyDescrip {
                shape: GalaxyShape::Lenticular,
                mass_density: vec![],
                rotation_curve: vec![],
                luminosity: vec![],
                eccentricity: 0.,
                arm_count: 2,
            },
            _ => unimplemented!(), // todo
        }
    }

    pub fn make_bodies(&self) -> Vec<Body> {
        self.descrip().make_bodies()
    }
}
