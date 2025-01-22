//! Data on specific galaxies.
//!
//! [SPARC](http://astroweb.cwru.edu/SPARC/) has tabular .dat data files of mass density and rotation curves.

use crate::{
    body_creation::{GalaxyDescrip, GalaxyShape, mass_density_from_lum},
    units::{ARCSEC_CONV_FACTOR, KPC_MYR_PER_KM_S},
    util::{scale_x_axis, zip_data},
};
// todo: Method to auto-parse from SPARC etc Rotmod dat files?

/// todo: Move specific galaxy creation to its own module A/R
#[derive(Clone, Copy, PartialEq, Default)]
pub enum GalaxyModel {
    #[default]
    Ngc1560,
    Ngc3198,
    Ngc3115,
    Ngc3031,
    Ngc7331,
    Ngc2685,
    Ngc2824,
    Ngc3626,
    Ugc6176,
    M31,
}

impl GalaxyModel {
    pub fn to_str(&self) -> String {
        match self {
            Self::Ngc1560 => "NGC 1560",
            Self::Ngc3198 => "NGC 3198",
            Self::Ngc3115 => "NGC 3115",
            Self::Ngc3031 => "NGC 3031",
            Self::Ngc7331 => "NGC 7331",
            Self::Ngc2685 => "NGC 2685",
            Self::Ngc2824 => "NGC 2824", // Shelest paper
            Self::Ngc3626 => "NGC 3626", // Shelest paper
            Self::Ugc6176 => "UGC 6176", // Shlest paper
            Self::M31 => "M31-NGC 224", // Andromeda
        }
            .to_owned()
    }

    pub fn descrip(&self) -> GalaxyDescrip {
        match self {
            /// Ludwig, Figures 3 and 5. todo: Partial/rough
            Self::Ngc1560 => ngc_1560(),
            Self::Ngc3198 => ngc_3198(),
            Self::Ngc3115 => ngc_3115(),
            Self::Ngc2685 => ngc_2685(),
            Self::Ngc2824 => ngc_2824(),
            Self::Ngc3626 => ngc_3636(),
            Self::Ugc6176 => ugc_6176(),
            Self::M31 => m31(),
            _ => unimplemented!(), // todo
        }
    }

    // pub fn make_bodies(&self, num_bodies_disk: usize, num_rings_disk: usize, num_bodies_disk: usize, num_rings_disk: usize,) -> Vec<Body> {
    //     self.descrip().make_bodies(num_bodies, num_rings)
    // }
}

/// A subset of GalaxyDescrip, from SPARC .dat file data. In units provided by SPARC, which
/// are not the same units we use internally.
struct SparcData {
    /// kpc
    pub r: Vec<f64>,
    /// X: r (kpc). Y:  M☉ / pc^2.
    pub mass_density_disk: Vec<f64>,
    /// X: r (kpc). Y: km/s
    pub velocity_disk: Vec<f64>,
    // /// Luminosity brightness profile. r (kpc), mu (mac arcsec^-2) -
    // pub luminosity_disk: Vec<(f64, f64)>,
    /// X: r (kpc). Y:  M☉ / pc^2. Note: This is only valid on the galactic plane. You must extrapolate
    /// to get a volume profile; the value you need is in M☉ / kpc^3?
    pub mass_density_bulge: Vec<f64>,
    /// X: r (kpc). Y: km/s
    pub velocity_bulge: Vec<f64>,
    pub mass_disk: f64,
    pub mass_bulge: f64,
}

impl SparcData {
    /// Handles unit conversions, and zipping radius with each param, since in the general case,
    /// velocity, mass, and luminosity data may not have the same radius indexes.
    fn galaxy_descrip(
        &self,
    ) -> (
        Vec<(f64, f64)>,
        Vec<(f64, f64)>,
        Vec<(f64, f64)>,
        Vec<(f64, f64)>,
    ) {
        // M☉/pc^2
        let density_disk_ = self.mass_density_disk.iter().map(|v| v * 1e6).collect(); // convert to M☉/kpc^2

        // kpc/MYR
        let velocity_disk = self
            .velocity_disk
            .iter()
            .map(|v| v * KPC_MYR_PER_KM_S)
            .collect();

        // At the disk radius indexies. M☉/pc^2
        let density_bulge_ = self.mass_density_bulge.iter().map(|v| v * 1e6).collect(); // convert to M☉/kpc^3;

        // kpc/MYR
        let velocity_bulge = self
            .velocity_bulge
            .iter()
            .map(|v| v * KPC_MYR_PER_KM_S)
            .collect();

        (
            zip_data(&self.r, density_disk_),
            zip_data(&self.r, velocity_disk),
            zip_data(&self.r, density_bulge_),
            zip_data(&self.r, velocity_bulge),
        )
    }
}

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

    let rotation_curve = scale_x_axis(&rot_curve_arcsec, α_conv_factor)
        .into_iter()
        .map(|(r, v)| (r, v * KPC_MYR_PER_KM_S))
        .collect();

    // let rotation_curve_corr = scale_x_axis(&rot_curve_corr_arcsec, α_conv_factor)
    //     .into_iter()
    //     .map(|(r, v)| (r, v * KPC_MYR_PER_KM_S))
    //     .collect();

    // let (mass_density_disk, rotation_curve_disk, mass_density_bulge, rotation_curve_bulge) =
    //     density_vel_from_sparc(
    //         &radius,
    //         &density_disk_,
    //         &velocity_disk_,
    //         &density_bulge_,
    //         &velocity_bulge_,
    //     );

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
        eccentricity: 0.18, // Broeils
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
    // NGC2685_disk.dat. kpc.  For disk and bulge.
    let radius = vec![
        0.00000, 0.05074, 0.05620, 0.06166, 0.06791, 0.07493, 0.08196, 0.09054, 0.09913, 0.10928,
        0.12020, 0.13191, 0.14518, 0.16001, 0.17562, 0.19358, 0.21309, 0.23417, 0.25758, 0.28334,
        0.31144, 0.34266, 0.37701, 0.41447, 0.45662, 0.50189, 0.55185, 0.60727, 0.66815, 0.73528,
        0.80865, 0.88905, 0.97803, 1.07560, 1.18331, 1.30196, 1.43231, 1.57515, 1.73282, 1.90610,
        2.09656, 2.30653, 2.53679, 2.79047, 3.06990, 3.37666, 3.71464, 4.08618, 4.49441, 4.94400,
        5.43887, 5.98214, 6.58082, 7.23882, 7.96239, 8.75855, 9.63433, 10.59831, 11.65830,
        12.82366, 14.10610, 15.51655, 16.92701, 18.33746, 19.74792, 21.15837, 22.56882, 23.97928,
        25.38973,
    ];

    // M☉/pc^2
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

    // M☉/pc^2
    let density_bulge_ = vec![
        29637.09388,
        29637.09388,
        29119.00646,
        28546.91739,
        27835.95108,
        27114.68817,
        26034.58429,
        25166.29367,
        24000.76381,
        22735.53888,
        21385.99481,
        19925.57520,
        18449.90519,
        17216.03165,
        15941.72231,
        14650.03410,
        13123.50973,
        11364.50947,
        9843.52233,
        8492.69425,
        7178.64917,
        5992.28305,
        5128.05857,
        4283.19400,
        3667.81694,
        2925.09072,
        2275.31105,
        1769.97708,
        1403.04449,
        967.60125,
        663.41545,
        473.96446,
        326.67018,
        217.20846,
        138.42412,
        84.27460,
        48.85553,
        26.88098,
        13.90071,
        6.73397,
        3.03608,
        1.26156,
        0.48155,
        0.16666,
        0.05179,
        0.01436,
        0.00349,
        0.00074,
        0.00014,
        0.00002,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
    ];

    // km/s
    let velocity_bulge_ = vec![
        0.000000, 56.773061, 64.501042, 71.702888, 79.234420, 87.507350, 95.315333, 103.639398,
        111.993003, 120.749348, 129.727592, 138.031938, 145.906989, 153.624711, 161.360114,
        169.621952, 178.742618, 186.237683, 192.964654, 198.777669, 203.353461, 206.522821,
        209.188005, 210.251855, 211.078492, 212.073347, 212.005465, 210.643292, 208.701420,
        206.235205, 201.292280, 197.279176, 191.666135, 185.813825, 179.443042, 172.799653,
        165.959388, 159.087032, 152.211996, 145.464481, 138.884213, 132.497063, 126.394519,
        120.532659, 114.961551, 109.589733, 104.486860, 99.640248, 95.042747, 90.694627, 86.470197,
        82.450350, 78.610541, 74.952634, 71.465933, 68.140401, 64.969576, 61.944461, 59.061324,
        56.313795, 53.692956, 51.194488, 49.015187, 47.092441, 45.379542, 43.840923, 42.448892,
        41.181553, 40.021351,
    ];

    let sparc_data = SparcData {
        r: radius,
        mass_density_disk: density_disk_,
        velocity_disk: velocity_disk_,
        mass_density_bulge: density_bulge_,
        velocity_bulge: velocity_bulge_,
        mass_disk: 20.6046e9,
        mass_bulge: 9.4320e9,
    };

    let (mass_density_disk, rotation_curve_disk, mass_density_bulge, rotation_curve_bulge) =
        sparc_data.galaxy_descrip();

    // todo: Mass and rotation for the bulge.

    GalaxyDescrip {
        shape: GalaxyShape::LenticularRingSeyfertType2,
        mass_density_disk,
        rotation_curve_disk,
        luminosity_disk: vec![], // todo
        mass_density_bulge,
        rotation_curve_bulge,
        luminosity_bulge: vec![], // todo
        eccentricity: 0.,         // todo
        arm_count: 0,             // todo
        burkert_params: (0., 0.), // todo
        r_s: 0.,                  // todo
        mass_disk: sparc_data.mass_disk,
        mass_bulge: sparc_data.mass_bulge,
        mass_to_light_ratio: 0.,  // todo
        dist_from_earth: 14.79e3, // Wikipedia
    }
}

/// SPARC Rotmod_ETG
/// NGC2824_disk.dat. kpc. For disk and bulge.
pub fn ngc_2824() -> GalaxyDescrip {
    let radius = vec![
        0.00000, 0.12479, 0.13823, 0.15167, 0.16703, 0.18431, 0.20159, 0.22270, 0.24382, 0.26878,
        0.29566, 0.32446, 0.35709, 0.39357, 0.43197, 0.47613, 0.52412, 0.57596, 0.63355, 0.69691,
        0.76603, 0.84282, 0.92729, 1.01945, 1.12312, 1.23447, 1.35734, 1.49365, 1.64340, 1.80851,
        1.98898, 2.18672, 2.40559, 2.64557, 2.91051, 3.20233, 3.52295, 3.87428, 4.26209, 4.68830,
        5.15675, 5.67319, 6.23955, 6.86351, 7.55082, 8.30532, 9.13662, 10.05048, 11.05457,
        12.16041, 13.37760,
    ];

    // M☉/pc^2
    let density_disk_ = vec![
        6238.98586, 5309.65746, 5218.22652, 5128.37001, 5027.56994, 4916.53615, 4807.95454,
        4678.49490, 4552.52112, 4408.00897, 4257.50647, 4101.95139, 3932.51506, 3751.41380,
        3569.78366, 3371.75383, 3168.94363, 2963.58823, 2750.98516, 2534.69479, 2318.08916,
        2055.79046, 1887.95677, 1665.83787, 1394.14145, 1191.40312, 1042.90608, 910.11888,
        773.23372, 633.74339, 495.12393, 370.43669, 268.85245, 192.62560, 147.33774, 114.89840,
        96.36382, 82.62548, 73.70789, 66.05626, 59.03560, 39.76619, 26.05637, 18.10982, 15.51365,
        11.76830, 10.51753, 8.11174, 6.36669, 4.49484, 4.44134,
    ];

    // At the disk radius indexies. km/s
    let velocity_disk_ = vec![
        0., 33.517614, 36.787605, 39.996933, 43.603965, 47.581888, 51.445770, 56.027996, 60.495887,
        65.630964, 70.949192, 76.463117, 82.466260, 88.778166, 95.181257, 102.123863, 109.186352,
        116.319076, 123.659274, 131.077472, 138.452257, 145.640479, 152.549992, 159.521405,
        165.708430, 170.410751, 174.592974, 178.756782, 182.648707, 185.866435, 187.791580,
        187.948079, 186.062472, 182.671018, 177.775649, 172.484678, 167.055254, 162.011551,
        157.382832, 153.755473, 151.275223, 148.862110, 143.914483, 137.866121, 131.945657,
        126.861749, 121.861989, 117.786450, 113.629696, 109.724267, 107.707366,
    ];

    // M☉/pc^2
    let density_bulge_ = vec![
        19929.70080,
        19929.70080,
        19066.07962,
        18169.49574,
        17433.27344,
        16491.24791,
        15400.27144,
        14345.66756,
        13132.79604,
        11892.65755,
        10478.35539,
        8884.86077,
        7311.13877,
        5941.92075,
        4585.60688,
        3315.97579,
        2536.46001,
        1730.24166,
        1067.57138,
        461.69833,
        204.88269,
        83.07041,
        30.77383,
        10.41623,
        3.07915,
        0.83166,
        0.19618,
        0.03951,
        0.00680,
        0.00098,
        0.00012,
        0.00001,
        0.,
        0.,
        0.,
        0.,
        0.,
        0.,
        0.,
        0.,
        0.,
        0.,
        0.,
        0.,
        0.,
        0.,
        0.,
        0.,
        0.,
        0.,
    ];

    // km/s
    let velocity_bulge_ = vec![
        0.0, 88.399177, 101.509195, 112.280335, 123.138977, 135.232161, 146.115845, 157.543821,
        168.533135, 179.951289, 191.930881, 202.089477, 210.782109, 219.682516, 224.190533,
        227.383737, 228.084656, 227.845947, 226.343141, 222.696708, 214.967925, 205.905089,
        196.866884, 188.352819, 178.575062, 171.446095, 163.536416, 155.689078, 148.472538,
        141.467596, 134.687368, 128.616150, 122.625640, 116.931754, 111.482717, 106.281837,
        101.330189, 96.626576, 92.125696, 87.838396, 83.753710, 79.850645, 76.140453, 72.597047,
        69.214173, 65.995437, 62.921535, 59.992739, 57.203311, 54.540358, 51.999958,
    ];

    let sparc_data = SparcData {
        r: radius,
        mass_density_disk: density_disk_,
        velocity_disk: velocity_disk_,
        mass_density_bulge: density_bulge_,
        velocity_bulge: velocity_bulge_,
        mass_disk: 31.1825e9,
        mass_bulge: 8.3897e9,
    };

    let (mass_density_disk, rotation_curve_disk, mass_density_bulge, rotation_curve_bulge) =
        sparc_data.galaxy_descrip();

    GalaxyDescrip {
        shape: GalaxyShape::Lenticular,
        mass_density_disk,
        rotation_curve_disk,
        luminosity_disk: vec![], // todo
        mass_density_bulge,
        rotation_curve_bulge,
        luminosity_bulge: vec![], // todo
        eccentricity: 0.,         // todo
        arm_count: 0,             // todo
        burkert_params: (0., 0.), // todo
        r_s: 0.,                  // todo
        mass_disk: sparc_data.mass_disk,
        mass_bulge: sparc_data.mass_bulge,
        mass_to_light_ratio: 0., // todo
        dist_from_earth: 0.,     // Not sure.
    }
}

// todo: Fill out
pub fn m31() -> GalaxyDescrip {
    let radius = vec![];

    // M☉/pc^2
    let density_disk_ = vec![];

    // At the disk radius indexies. km/s
    let velocity_disk_ = vec![];

    // M☉/pc^2
    let density_bulge_ = vec![];

    // km/s
    let velocity_bulge_ = vec![];

    let sparc_data = SparcData {
        r: radius,
        mass_density_disk: density_disk_,
        velocity_disk: velocity_disk_,
        mass_density_bulge: density_bulge_,
        velocity_bulge: velocity_bulge_,
        mass_disk: 31.1825e9,
        mass_bulge: 8.3897e9,
    };

    let (mass_density_disk, rotation_curve_disk, mass_density_bulge, rotation_curve_bulge) =
        sparc_data.galaxy_descrip();

    GalaxyDescrip {
        shape: GalaxyShape::Lenticular,
        mass_density_disk,
        rotation_curve_disk,
        luminosity_disk: vec![], // todo
        mass_density_bulge,
        rotation_curve_bulge,
        luminosity_bulge: vec![], // todo
        eccentricity: 0.,         // todo
        arm_count: 0,             // todo
        burkert_params: (0., 0.), // todo
        r_s: 0.,                  // todo
        mass_disk: sparc_data.mass_disk,
        mass_bulge: sparc_data.mass_bulge,
        mass_to_light_ratio: 0., // todo
        dist_from_earth: 0.,     // Not sure.
    }
}

// todo: Fill out
/// From SPARC, and paper
pub fn ngc_3636() -> GalaxyDescrip {
    let radius = vec![
        0.00000, 0.07216, 0.07994, 0.08771, 0.09659, 0.10658, 0.11657, 0.12879, 0.14100, 0.15543,
        0.17097, 0.18763, 0.20650, 0.22760, 0.24980, 0.27534, 0.30309, 0.33307, 0.36637, 0.40301,
        0.44298, 0.48739, 0.53624, 0.58953, 0.64948, 0.71387, 0.78493, 0.86375, 0.95035, 1.04583,
        1.15019, 1.26454, 1.39111, 1.52989, 1.68310, 1.85185, 2.03726, 2.24043, 2.46470, 2.71117,
        2.98206, 3.28071, 3.60823, 3.96905, 4.36651, 4.80283, 5.28355, 5.81202, 6.39267, 7.03215,
        7.73604, 8.50875, 9.36029, 10.29621, 11.32539, 12.45782, 13.70349, 15.07461, 16.58230,
        18.23986,
    ];

    // M☉/pc^2
    let density_disk_ = vec![
        1976.40336, 1923.06796, 1917.41069, 1911.77007, 1905.34395, 1898.14040, 1890.96407,
        1882.22984, 1873.53595, 1863.31312, 1852.36628, 1840.70889, 1827.58584, 1813.02958,
        1797.83244, 1780.51314, 1761.87700, 1741.96900, 1720.11271, 1696.38740, 1670.87830,
        1642.98449, 1612.83884, 1580.58315, 1545.06594, 1507.80630, 1467.73411, 1424.52387,
        1378.51847, 1329.51480, 1277.94302, 1223.72759, 1166.39968, 1094.73185, 1038.33195,
        991.10077, 904.99363, 800.70140, 690.35883, 559.59716, 453.60264, 368.36264, 293.95143,
        231.56642, 189.44170, 154.41006, 116.59384, 83.15255, 58.43534, 49.05413, 41.25490,
        29.31403, 23.26353, 17.61471, 13.59802, 8.68312, 5.39853, 2.66119, 1.91900, 0.83922,
    ];

    // At the disk radius indexies. km/s
    let velocity_disk_ = vec![
        0.000000, 9.833079, 10.870669, 11.901031, 13.071282, 14.378063, 15.674886, 17.247795,
        18.804476, 20.626984, 22.570611, 24.631977, 26.940055, 29.485366, 32.122622, 35.104202,
        38.296060, 41.682904, 45.364003, 49.342064, 53.598720, 58.191988, 63.124232, 68.338247,
        74.075058, 80.007684, 86.349486, 93.116151, 100.279133, 107.863446, 115.779528,
        124.080218, 132.876313, 141.857920, 151.035739, 161.445752, 172.849141, 183.414661,
        192.704414, 199.754085, 203.297757, 204.979085, 205.271037, 203.685342, 201.025874,
        198.796267, 196.223946, 191.499709, 184.400383, 176.639563, 171.300876, 166.267545,
        159.697726, 154.141116, 148.706218, 143.506151, 137.121326, 130.369985, 123.351035,
        118.086839,
    ];

    // M☉/pc^2
    let density_bulge_ = vec![
        68628.50634, 68628.50634, 66423.43908, 64286.61643, 61916.83667, 59574.53437, 57094.99441,
        54081.14089, 50473.91952, 46821.85189, 43379.37337, 38844.84654, 34240.46060, 29845.07327,
        26323.59338, 22475.73430, 19081.80530, 16203.81633, 12938.82242, 10551.93825, 8645.35359,
        6901.69597, 5213.48420, 3932.44190, 2917.67431, 2149.67731, 1641.42443, 1179.78967,
        855.69300, 561.28383, 334.08879, 152.65729, 46.26051, 12.49315, 2.94446, 0.59933,
        0.10425, 0.01534, 0.00185, 0.00018, 0.00002, 0.00000, 0.00000, 0.00000, 0.00000,
        0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
        0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,];

    // km/s
    let velocity_bulge_ = vec![
        0.000000, 114.569865, 130.417300, 145.114151, 159.809702, 175.233893, 189.386936,
        206.006674, 221.301723, 237.263683, 253.234833, 268.917507, 283.770990, 295.877263,
        306.345735, 316.659999, 324.776489, 331.162745, 337.497005, 340.898969, 341.208815,
        340.317195, 339.328442, 335.958417, 330.391159, 324.001759, 316.146018, 308.113147,
        299.411775, 290.357793, 280.892678, 270.810973, 259.920976, 248.315467, 236.915804,
        225.975411, 215.454260, 205.415586, 195.808481, 186.670890, 177.962590, 169.669189,
        161.785576, 154.256553, 147.068505, 140.229149, 133.697718, 127.474511, 121.547425,
        115.889168, 110.491165, 105.354727, 100.448222, 95.774129, 91.318827, 87.069466,
        83.017818, 79.152359, 75.468292, 71.957515,
    ];

    let sparc_data = SparcData {
        r: radius,
        mass_density_disk: density_disk_,
        velocity_disk: velocity_disk_,
        mass_density_bulge: density_bulge_,
        velocity_bulge: velocity_bulge_,
        mass_disk: 52.7180e9,
        mass_bulge: 21.9047e9,
    };

    let (mass_density_disk, rotation_curve_disk, mass_density_bulge, rotation_curve_bulge) =
        sparc_data.galaxy_descrip();

    GalaxyDescrip {
        shape: GalaxyShape::Lenticular, // todo
        mass_density_disk,
        rotation_curve_disk,
        luminosity_disk: vec![], // todo
        mass_density_bulge,
        rotation_curve_bulge,
        luminosity_bulge: vec![], // todo
        eccentricity: 0.,         // todo
        arm_count: 0,             // todo
        burkert_params: (0., 0.), // todo
        r_s: 0.,                  // todo
        mass_disk: sparc_data.mass_disk,
        mass_bulge: sparc_data.mass_bulge,
        mass_to_light_ratio: 0., // todo
        dist_from_earth: 0.,     // Not sure.
    }
}

// todo: Fill out
pub fn ugc_6176() -> GalaxyDescrip {
    let radius = vec![
        0.00000, 0.12385, 0.13718, 0.15052, 0.16576, 0.18291, 0.20006, 0.22102, 0.24198, 0.26674,
        0.29342, 0.32200, 0.35439, 0.39059, 0.42870, 0.47252, 0.52015, 0.57160, 0.62875, 0.69163,
        0.76022, 0.83643, 0.92027, 1.01172, 1.11461, 1.22512, 1.34706, 1.48234, 1.63095, 1.79481,
        1.97391, 2.17016, 2.38736, 2.62553, 2.88846, 3.17807, 3.49626, 3.84493, 4.22981, 4.65279,
        5.11768, 5.63021, 6.19228, 6.81151, 7.49361, 8.24240, 9.06741, 9.97434, 10.97082, 12.06828,
        13.27625,
    ];

    // M☉/pc^2
    let density_disk_ = vec![
        4510.40087, 3817.97191, 3750.05553, 3683.34730, 3608.56112, 3526.24014, 3445.79713,
        3349.96580, 3256.79965, 3150.02905, 3038.95662, 2924.29288, 2799.56195, 2666.44203,
        2533.14840, 2388.07738, 2239.79879, 2089.98161, 1935.24162, 1778.22701, 1608.69851,
        1439.30903, 1304.35081, 1172.16506, 1014.75891, 868.84963, 732.82096, 606.19384, 493.19026,
        392.83940, 313.77268, 240.66622, 192.58175, 155.24411, 129.00756, 110.41171, 100.32638,
        91.58306, 84.60872, 65.61682, 50.60753, 36.62816, 23.19611, 15.28323, 10.70077, 7.27470,
        4.47317, 2.51082, 1.53821, 1.01162, 0.41325,
    ];

    // At the disk radius indexies. km/s
    let velocity_disk_ = vec![
        0.000000, 33.877963, 37.112867, 40.278384, 43.813323, 47.686778, 51.455444, 55.918653,
        60.226403, 65.128268, 70.204495, 75.417094, 81.049706, 86.976212, 92.882054, 99.244169,
        105.668078, 112.074487, 118.590959, 125.093021, 131.370656, 137.271110, 142.876922,
        148.526649, 153.980662, 158.413294, 161.963867, 164.534078, 165.962947, 166.142473,
        165.160235, 163.011767, 159.444175, 155.527053, 151.186997, 146.860160, 143.075042,
        140.539389, 139.521069, 139.346424, 137.310678, 134.740138, 130.784647, 125.108250,
        119.454212, 114.204857, 108.966315, 103.488531, 98.026380, 93.108608, 88.763872,
    ];

    // M☉/pc^2
    let density_bulge_ = vec![
        17558.88885, 17558.88885, 17112.23442, 16476.87729, 15665.89842, 14704.51322,
        13745.19616, 12785.28093, 11512.53583, 10385.02567, 9031.23957, 7694.87133,
        6210.34791, 5067.65689, 3723.31146, 2547.43712, 1607.43208, 889.97162, 440.57540,
        195.28228, 80.38629, 29.98265, 10.13274, 3.10279, 0.81946, 0.19610, 0.04047, 0.00703,
        0.00103, 0.00013, 0.00001, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
        0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
        0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
    ];

    // km/s
    let velocity_bulge_ = vec![
        0.000000, 81.710823, 93.364608, 104.373964, 115.855284, 127.719002, 138.287674,
        149.738838, 160.663448, 171.392185, 182.111633, 191.960778, 200.137109, 206.841095,
        211.675565, 215.795167, 215.976056, 213.631785, 208.581629, 201.678924, 194.484106,
        185.361991, 178.230888, 169.341827, 160.993517, 153.619877, 146.510294, 139.687622,
        133.142175, 126.842946, 121.001574, 115.400799, 110.026105, 104.917082, 100.027984,
        95.361476, 90.918620, 86.698279, 82.659784, 78.813020, 75.148122, 71.646074,
        68.317097, 65.137786, 62.102510, 59.214471, 56.456385, 53.828541, 51.325732,
        48.936399, 46.657019,
    ];

    let sparc_data = SparcData {
        r: radius,
        mass_density_disk: density_disk_,
        velocity_disk: velocity_disk_,
        mass_density_bulge: density_bulge_,
        velocity_bulge: velocity_bulge_,
        mass_disk: 22.2694e9,
        mass_bulge: 6.7031e9,
    };

    let (mass_density_disk, rotation_curve_disk, mass_density_bulge, rotation_curve_bulge) =
        sparc_data.galaxy_descrip();

    GalaxyDescrip {
        shape: GalaxyShape::Lenticular,
        mass_density_disk,
        rotation_curve_disk,
        luminosity_disk: vec![], // todo
        mass_density_bulge,
        rotation_curve_bulge,
        luminosity_bulge: vec![], // todo
        eccentricity: 0.,         // todo
        arm_count: 0,             // todo
        burkert_params: (0., 0.), // todo
        r_s: 0.,                  // todo
        mass_disk: sparc_data.mass_disk,
        mass_bulge: sparc_data.mass_bulge,
        mass_to_light_ratio: 0., // todo
        dist_from_earth: 0.,     // Not sure.
    }
}
