//! Data on specific galaxies.

use crate::{
    body_creation::{GalaxyDescrip, GalaxyShape},
    units::ARCSEC_CONV_FACTOR,
    util::scale_x_axis,
};

pub fn ngc_1560() -> GalaxyDescrip {
    // These rotation curve values are from Broeils.
    // X axis is ''.
    // todo: Use theta and i? i is always 80. Theta ranges from 20.1 to 22.7

    // X: arcsec(''), Y: μ (mag arcsec^2) - surface brightness profile.
    let luminosity_arcsec = vec![
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

    // X: arcsec (''). Y: km/s
    let rot_curve_arcsec = vec![
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

    // X: arcsec (''). Y: km/s
    let rot_curve_corr_arcsec = vec![
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

    let dist_from_earth = 2_990.; // Wikipedia, J2000 epoch, converted from Mly.

    // Convert the x values from arcsec ('') to kpc.
    let α_conv_factor = ARCSEC_CONV_FACTOR * dist_from_earth;
    let luminosity = scale_x_axis(&luminosity_arcsec, α_conv_factor);
    let rotation_curve = scale_x_axis(&rot_curve_arcsec, α_conv_factor);
    let rotation_curve_corr = scale_x_axis(&rot_curve_corr_arcsec, α_conv_factor);

    // solar masses
    let mass_total = 1.2e10; // todo temp. H mass: 8.2e8 solar masses.

    // Generate mass density from luminosity; we multiply by a mapping between light
    // and mass, and convert the tabular data in terms of arcsec^-2, to m.

    let mass_to_light_ratio = 35.; // Broeils.

    // todo: Come back to this; re-examine converting surface brightness to solar luminosity etc.
    // todo: You are likely missing a step.
    let mut mass_density = Vec::with_capacity(luminosity.len());
    for (i, (r, lum)) in luminosity.iter().enumerate() {
        // surface brightness profile (μ)
        // μ = mag / arcsec^2
        // mag = μ * arcsec^2
        // todo. Hmm: Why are we then, dividing?
        let mut arcsec_sq = luminosity_arcsec[i].0.powi(2);
        if arcsec_sq < 1e-14 {
            arcsec_sq = 1.; // avoid div by 0. // todo: Kludge
        }
        // mass_density.push((*r, mass_to_light_ratio * r / arcsec_sq))
    }

    // Let's try a simple technique where, to find mass density, we scale luminosity by the total mass.
    // todo: You really need to QC this!
    let mut lum_total = 0.;
    for (r, lum) in &luminosity {
        if *r < 1e-14 {
            continue; // todo kludege
        }
        lum_total += lum / r.powi(2);
    }

    let lum_scaler = mass_total / lum_total;
    for (r, lum) in &luminosity {
        if *r < 1e-14 {
            continue; // todo kludege
        }
        mass_density.push((*r, lum_scaler * lum / r.powi(2)));
    }

    // todo temp
    // for (r, lum) in &luminosity {
    //     println!("R: {r} lum: {lum}");
    // }

    // for (r, mass) in &mass_density{
    //     println!("R: {r} mass: {mass}");
    // }

    GalaxyDescrip {
        shape: GalaxyShape::FlocculentSpiral, // todo ?
        mass_density,
        rotation_curve,
        luminosity,
        eccentricity: 0.0, // todo temp
        // eccentricity: 0.18, // Broeils
        arm_count: 2,
        r_s: 1.46e-6, // todo?
        // total H mass: 8.2e8 solar masses // Broeils
        // Total blue luminosity: 3.45e8 solar luminosities // Broeils
        mass_total,
        mass_to_light_ratio,
        dist_from_earth,
        // gas-to-blue luminosity ratio
        //M_HI / L_B = 2.4
    }
}

pub fn ngc_3198() -> GalaxyDescrip {
    GalaxyDescrip {
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
        mass_total: 0.,
        mass_to_light_ratio: 0.,  // todo
        dist_from_earth: 14_410., // Wikipedia, J2000 epoch
    }
}

pub fn ngc_3115() -> GalaxyDescrip {
    GalaxyDescrip {
        shape: GalaxyShape::Lenticular,
        mass_density: vec![],
        rotation_curve: vec![],
        luminosity: vec![],
        eccentricity: 0.,
        arm_count: 2,
        r_s: 6.97e-16,
        mass_total: 0.,
        mass_to_light_ratio: 0., // todo
        dist_from_earth: 9_700., // Wikipedia, J2000 epoch.
    }
}
