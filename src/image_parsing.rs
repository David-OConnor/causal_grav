//! For estimating galaxy parameters based on telescope images.

use crate::body_creation::{GalaxyDescrip, GalaxyShape};

/// Dist is in kpc. `image_buf` is a bitmap.
pub fn examine_image(image_buf: &[u8], dist: f64) -> GalaxyDescrip {
    GalaxyDescrip {
        shape: GalaxyShape::BarredSpiral,
        mass_density_disk: Vec::new(),
        rotation_curve_disk: Vec::new(),
        luminosity_disk: Vec::new(),
        mass_density_bulge: Vec::new(),
        rotation_curve_bulge: Vec::new(),
        luminosity_bulge: Vec::new(),
        eccentricity: 0.,
        arm_count: 0,
        burkert_params: (0., 0.),
        r_s: 0.,
        mass_bulge: 0.,
        mass_disk: 0.,
        mass_to_light_ratio: 0.,
        dist_from_earth: 0.,
    }
}
