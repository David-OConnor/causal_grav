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
    /// todo: Consider removing `luminosity`; unused.
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
            Self::Ngc1560 => {
                // These rotation curve values are from Broeils.
                // X axis is ''. (Arcseconds?
                // todo: Use theta and i? i is always 80. Theta ranges from 20.1 to 22.7

                // Multiply arcsecond measurements by this to get distance (i.e. radius?)
                // in kpc. This applies to both rotation curve, and luminosity. (Note: Broils indicates
                // this using two different, equivalent conventions)
                let α_conv_factor = 0.01455;

                let luminosity = vec![
                    (0.0, 22.27),
                    (2.0, 22.30),
                    (4.0, 22.31),
                    (6.0, 22.29),
                    (8.0, 22.29),
                    (10.0, 22.27),
                    (12.0, 22.27),
                    (14.0, 22.27),
                    (16.0, 22.27),
                    (18.0, 22.30),
                    (20.0, 22.31),
                    (22.0, 22.32),
                    (24.0, 22.36),
                    (26.0, 22.35),
                    (28.0, 22.36),
                    (30.0, 22.42),
                    (32.0, 22.36),
                    (34.0, 22.36),
                    (36.0, 22.38),
                    (38.0, 22.36),
                    (40.0, 22.41),
                    (42.0, 22.41),
                    (44.0, 22.47),
                    (46.0, 22.47),
                    (48.0, 22.50),
                    (50.0, 22.51),
                    (52.0, 22.50),
                    (54.0, 22.54),
                    (56.0, 22.54),
                    (58.0, 22.57),
                    (60.0, 22.60),
                    (62.0, 22.63),
                    (64.0, 22.65),
                    (66.0, 22.68),
                    (68.0, 22.73),
                    (70.0, 22.75),
                    (72.0, 22.81),
                    (74.0, 22.81),
                    (76.0, 22.83),
                    (78.0, 22.89),
                    (80.0, 22.95),
                    (82.0, 22.92),
                    (84.0, 22.98),
                    (86.0, 23.01),
                    (88.0, 23.02),
                    (90.0, 23.03),
                    (92.0, 23.09),
                    (94.0, 23.11),
                    (96.0, 23.14),
                    (98.0, 23.11),
                    (100.0, 23.13),
                    (102.0, 23.20),
                    (104.0, 23.22),
                    (106.0, 23.23),
                    (108.0, 23.24),
                    (110.0, 23.24),
                    (112.0, 23.32),
                    (114.0, 23.36),
                    (116.0, 23.28),
                    (118.0, 23.34),
                    (120.0, 23.39),
                    (122.0, 23.43),
                    (124.0, 23.46),
                    (126.0, 23.46),
                    (128.0, 23.49),
                    (130.0, 23.52),
                    (132.0, 23.53),
                    (134.0, 23.57),
                    (136.0, 23.62),
                    (138.0, 23.60),
                    (140.0, 23.62),
                    (142.0, 23.63),
                    (144.0, 23.61),
                    (146.0, 23.67),
                    (148.0, 23.70),
                    (150.0, 23.73),
                    (152.0, 23.75),
                    (154.0, 23.79),
                    (156.0, 23.84),
                    (158.0, 23.80),
                    (160.0, 23.84),
                    (162.0, 23.90),
                    (164.0, 23.90),
                    (166.0, 23.95),
                    (168.0, 24.01),
                    (170.0, 24.02),
                    (172.0, 24.04),
                    (174.0, 24.04),
                    (176.0, 24.12),
                    (178.0, 24.13),
                    (180.0, 24.11),
                    (182.0, 24.18),
                    (184.0, 24.18),
                    (186.0, 24.25),
                    (188.0, 24.24),
                    (190.0, 24.23),
                    (192.0, 24.29),
                    (194.0, 24.30),
                    (196.0, 24.34),
                    (198.0, 24.34),
                    (200.0, 24.32),
                    (202.0, 24.36),
                    (204.0, 24.41),
                    (206.0, 24.44),
                    (208.0, 24.43),
                    (210.0, 24.41),
                    (212.0, 24.54),
                    (214.0, 24.55),
                    (216.0, 24.47),
                    (218.0, 24.51),
                    (220.0, 24.57),
                    (222.0, 24.56),
                    (224.0, 24.57),
                    (226.0, 24.59),
                    (228.0, 24.67),
                    (230.0, 24.73),
                    (232.0, 24.75),
                    (234.0, 24.71),
                    (236.0, 24.82),
                    (238.0, 24.86),
                    (240.0, 24.86),
                    (242.0, 24.91),
                    (244.0, 24.89),
                    (246.0, 24.99),
                    (248.0, 24.93),
                    (250.0, 24.97),
                    (253.0, 24.99),
                    (257.0, 25.05),
                    (261.0, 25.11),
                    (265.0, 25.08),
                    (269.0, 25.16),
                    (273.0, 25.18),
                    (277.0, 25.37),
                    (281.0, 25.33),
                    (285.0, 25.43),
                    (289.0, 25.37),
                    (293.0, 25.39),
                    (297.0, 25.44),
                    (301.0, 25.41),
                    (305.0, 25.61),
                    (309.0, 25.63),
                    (313.0, 25.72),
                    (317.0, 25.73),
                    (321.0, 25.75),
                    (325.0, 25.60),
                    (329.0, 25.80),
                    (333.0, 25.88),
                    (337.0, 25.81),
                    (341.0, 25.94),
                    (345.0, 25.99),
                    (349.0, 25.91),
                    (353.0, 26.13),
                ];

                // X axis: arc sec.
                let rot_curve_raw = vec![
                    (15., 4.5),
                    (30., 8.4),
                    (45., 14.0),
                    (60., 26.1),
                    (75., 28.7),
                    (90., 27.9),
                    (105., 32.1),
                    (120., 43.0),
                    (135., 47.8),
                    (150., 47.6),
                    (165., 49.8),
                    (180., 52.6),
                    (195., 56.4),
                    (210., 58.3),
                    (225., 58.8),
                    (240., 59.4),
                    (255., 59.1),
                    (270., 59.9),
                    (285., 61.9),
                    (300., 60.9),
                    (315., 60.6),
                    (330., 62.1),
                    (345., 64.5),
                    (360., 65.6),
                    (375., 67.1),
                    (390., 68.7),
                    (405., 70.4),
                    (420., 71.3),
                    (435., 72.0),
                    (450., 72.3),
                    (465., 73.2),
                    (480., 74.1),
                    (510., 74.4),
                    (540., 75.2),
                    (570., 76.6),
                ];

                let rot_curve_raw_corr = vec![
                    (15., 5.0),
                    (30., 8.9),
                    (45., 14.5),
                    (60., 24.6),
                    (75., 28.9),
                    (90., 27.8),
                    (105., 31.8),
                    (120., 42.8),
                    (135., 48.2),
                    (150., 48.4),
                    (165., 50.6),
                    (180., 53.5),
                    (195., 57.2),
                    (210., 59.1),
                    (225., 59.8),
                    (240., 60.3),
                    (255., 60.1),
                    (270., 62.1),
                    (285., 63.6),
                    (300., 62.0),
                    (315., 60.5),
                    (330., 63.1),
                    (345., 63.8),
                    (360., 66.1),
                    (375., 67.7),
                    (390., 70.4),
                    (405., 73.0),
                    (420., 74.2),
                    (435., 75.1),
                    (450., 75.2),
                    (465., 76.3),
                    (480., 77.2),
                    (510., 76.9),
                    (540., 77.5),
                    (570., 78.7),
                ];

                let rot_curve: Vec<(f64, f64)> = rot_curve_raw
                    .iter()
                    .map(|(x, y)| (α_conv_factor * x, *y))
                    .collect();
                let rot_curve_corr: Vec<(f64, f64)> = rot_curve_raw_corr
                    .iter()
                    .map(|(x, y)| (α_conv_factor * x, *y))
                    .collect();

                GalaxyDescrip {
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
                    // Broeils. "
                    // X axis: Arcsec (''). 1'' is 14.5pc at 3.0 Mpc"
                    // Y axis is mu; mag arcsec^-2
                    luminosity: vec![],
                    eccentricity: 0.0, // todo temp
                    // eccentricity: 0.18, // Broeils
                    arm_count: 2,
                    r_s: 1.46e-6,
                    // total H mass: 8.2e8 solar masses // Broeils
                    // Total blue luminosity: 3.45e8 solar luminosities // Broeils
                }
            }
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
