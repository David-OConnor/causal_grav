//! Data on specific galaxies.
//!
//! [SPARC](http://astroweb.cwru.edu/SPARC/) has tabular .dat data files of mass density and rotation curves.

use crate::{
    body_creation::{mass_density_from_lum, GalaxyDescrip, GalaxyShape},
    units::ARCSEC_CONV_FACTOR,
    util::{scale_x_axis, zip_data},
};
// todo: Method to auto-parse from SPARC etc Rotmod dat files?

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
    let dist_from_earth = 3_270.; // Jacobs et al. (2009)

    // Convert the x values from arcsec ('') to kpc.
    let α_conv_factor = ARCSEC_CONV_FACTOR * dist_from_earth;
    let luminosity = scale_x_axis(&luminosity_arcsec, α_conv_factor);
    let rotation_curve = scale_x_axis(&rot_curve_arcsec, α_conv_factor);
    let rotation_curve_corr = scale_x_axis(&rot_curve_corr_arcsec, α_conv_factor);

    // solar masses
    let mass_disk = 1.2e10; // H mass: 8.2e8 solar masses.

    // Generate mass density from luminosity; we multiply by a mapping between light
    // and mass, and convert the tabular data in terms of arcsec^-2, to m.

    let mass_to_light_ratio = 35.; // Broeils.

    let mass_density = mass_density_from_lum(&luminosity, mass_disk, &luminosity_arcsec);

    GalaxyDescrip {
        shape: GalaxyShape::FlocculentSpiral, // todo ?
        mass_density_disk: mass_density,
        rotation_curve_disk: rotation_curve,
        luminosity_disk: luminosity,
        mass_density_bulge: vec![], // Thin disk
        rotation_curve_bulge: vec![],
        luminosity_bulge: vec![],
        eccentricity: 0.0, // todo temp
        // eccentricity: 0.18, // Broeils
        arm_count: 0,
        // Gentile (2024), section 6. Note: s_0 is 0.8e-24 g/cm^3.
        burkert_params: (5.6, 1.182e7),
        r_s: 1.46e-6, // todo? For nfw Halo?
        // total H mass: 8.2e8 solar masses // Broeils
        // Total blue luminosity: 3.45e8 solar luminosities // Broeils
        mass_disk,
        mass_bulge: 0., // Our data is from a thin disk model.
        mass_to_light_ratio,
        dist_from_earth,
        // gas-to-blue luminosity ratio
        //M_HI / L_B = 2.4
    }
}

/// Python lib: https://github.com/hsalas/rotation_curves/blob/master/data/ngc3198.dat
pub fn ngc_3198() -> GalaxyDescrip {
    let dist_from_earth = 47_000.;

    // Convert the x values from arcsec ('') to kpc.
    let α_conv_factor = ARCSEC_CONV_FACTOR * dist_from_earth;

    let rotation_curve = vec![];

    let luminosity = vec![];

    let mass_density = vec![];

    GalaxyDescrip {
        shape: GalaxyShape::BarredSpiral,
        mass_density_disk: mass_density,
        rotation_curve_disk: rotation_curve,
        luminosity_disk: luminosity,
        mass_density_bulge: vec![],
        rotation_curve_bulge: vec![],
        luminosity_bulge: vec![],
        eccentricity: 0.,
        arm_count: 0,
        burkert_params: (0., 0.),
        r_s: 1.2e-5,
        mass_bulge: 0.,
        mass_disk: 0.,
        mass_to_light_ratio: 0., // todo
        dist_from_earth,
    }
}

pub fn ngc_3115() -> GalaxyDescrip {
    GalaxyDescrip {
        shape: GalaxyShape::Lenticular,
        mass_density_disk: vec![],
        rotation_curve_disk: vec![],
        luminosity_disk: vec![],
        mass_density_bulge: vec![],
        rotation_curve_bulge: vec![],
        luminosity_bulge: vec![],
        eccentricity: 0.,
        arm_count: 2,
        burkert_params: (0., 0.),
        r_s: 6.97e-16,
        mass_disk: 0.,
        mass_bulge: 0.,
        mass_to_light_ratio: 0., // todo
        dist_from_earth: 9_700., // Wikipedia, J2000 epoch.
    }
}

/// SPARC Rotmod_ETG
pub fn ngc_2685() -> GalaxyDescrip {
    // NGC2685_disk.dat. kpc
    let radius_disk = vec![
        0.00000, 0.05074, 0.05620, 0.06166, 0.06791, 0.07493, 0.08196, 0.09054, 0.09913, 0.10928,
        0.12020, 0.13191, 0.14518, 0.16001, 0.17562, 0.19358, 0.21309, 0.23417, 0.25758, 0.28334,
        0.31144, 0.34266, 0.37701, 0.41447, 0.45662, 0.50189, 0.55185, 0.60727, 0.66815, 0.73528,
        0.80865, 0.88905, 0.97803, 1.07560, 1.18331, 1.30196, 1.43231, 1.57515, 1.73282, 1.90610,
        2.09656, 2.30653, 2.53679, 2.79047, 3.06990, 3.37666, 3.71464, 4.08618, 4.49441, 4.94400,
        5.43887, 5.98214, 6.58082, 7.23882, 7.96239, 8.75855, 9.63433, 10.59831, 11.65830,
        12.82366, 14.10610, 15.51655, 16.92701, 18.33746, 19.74792, 21.15837, 22.56882, 23.97928,
        25.38973,
    ];

    // At the disk radius indexies. MSUN/pc^2
    let density_disk_ = vec![
        2025.89108, 1939.86514, 1930.82150, 1921.82001, 1911.58398, 1900.13359, 1888.75180,
        1874.93325, 1861.21581, 1845.13362, 1827.96971, 1809.75695, 1789.33508, 1766.78323,
        1743.35144, 1716.78880, 1688.37521, 1658.21649, 1625.33820, 1589.92453, 1552.17086,
        1511.27283, 1467.52856, 1421.25035, 1370.92946, 1318.86387, 1263.70392, 1205.20599,
        1144.05741, 1080.22699, 1014.52534, 947.11242, 926.21031, 749.65245, 691.03713, 623.52740,
        546.30428, 487.00515, 429.64095, 371.31642, 327.45032, 287.02406, 242.94126, 196.28233,
        148.23835, 108.79154, 77.59057, 58.31954, 44.03678, 32.82585, 28.80156, 21.64792, 16.15162,
        12.36563, 10.59293, 9.29429, 7.65978, 5.52352, 3.87807, 2.45141, 1.29365, 1.13715, 0.68090,
        0.43749, 0.28110, 0.18062, 0.11605, 0.07457, 0.04791,
    ];

    // At the disk radius indexies. km/s
    let velocity_disk_ = vec![
        0.000000, 8.206425, 9.063075, 9.915770, 10.886032, 11.965243, 13.035349, 14.326730,
        15.605824, 17.097387, 18.691297, 20.397338, 22.333736, 24.468565, 26.638496, 29.071377,
        31.620722, 34.323400, 37.306805, 40.565222, 43.968341, 47.506242, 51.283070, 55.683579,
        60.024268, 64.348085, 69.218281, 74.254233, 79.416406, 84.811164, 90.311281, 95.878239,
        102.117083, 108.131218, 112.110334, 116.642601, 120.773710, 124.344712, 127.809843,
        130.869902, 133.457396, 136.050867, 138.539741, 141.619056, 141.577247, 140.278912,
        138.028615, 133.877826, 129.512849, 124.828333, 120.324483, 116.895271, 112.645003,
        108.022512, 103.473958, 99.909351, 96.976756, 94.034644, 90.582204, 86.943555, 82.730548,
        78.376724, 75.321101, 72.144864, 69.337426, 66.795913, 64.500683, 62.430336, 60.624714,
    ];

    // NGC2685_disk.dat. kpc.
    let radius_bulge = vec![];

    // At the disk radius indexies. MSUN/PC^2
    let density_bulge_ = vec![];

    // km/s
    let velocity_bulge_ = vec![];

    let mass_density_disk = zip_data(&radius_disk, density_disk_);
    let velocity_disk = zip_data(&radius_disk, velocity_disk_);
    let mass_density_bulge = zip_data(&radius_bulge, density_bulge_);
    let velocity_bulge = zip_data(&radius_bulge, velocity_bulge_);

    let mass_disk = 20.6046e9;
    let mass_bulge = 8.3897e9;

    // todo: Mass and rotation for the bulge.

    GalaxyDescrip {
        shape: GalaxyShape::LenticularRingSeyfertType2,
        mass_density_disk,
        rotation_curve_disk: velocity_disk,
        luminosity_disk: vec![],       // todo
        mass_density_bulge,
        rotation_curve_bulge: velocity_bulge,
        luminosity_bulge: vec![],       // todo
        eccentricity: 0.,         // todo
        arm_count: 0,             // todo
        burkert_params: (0., 0.), // todo
        r_s: 0.,                  // todo
        mass_disk,
        mass_bulge,
        mass_to_light_ratio: 0., // todo
        dist_from_earth: 14.79e3,     // Wikipedia
    }
}
