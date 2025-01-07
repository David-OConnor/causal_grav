//! Contains unit conversions, and similar.
//!
//!
//! // Let's define base units for our sim:

// Basic:
// Dist: kpc
// Time: Myr: (million years)
// Mass: M☉ x 10^10
//
// Derived:
// Speed: kpc / Myr
// Accel kpc / Myr^2
// Force: M☉ x 10^10 x kpc / Myr^2
// Energy: M☉ x 10^10 x (kpc / Myr)^2

use std::f64::consts::TAU;

// Cache, vice using `.to_radians`.
pub const ARCSEC_CONV_FACTOR: f64 = TAU / (360. * 3_600.);

// Multiply our native units by these to get SI units.

pub const KPC: f64 = 3.086e16; // KPC in meters
pub const KPC_KM: f64 = 3.086e13; // KPC in km
pub const MYR: f64 = 3.15576e13; // 1 million years in seconds
pub const SOLAR_MASS: f64 = 1.989e30; // Yikes on the size of this in SI.
pub const SOLAR_MASS_E10: f64 = 1.989e20;

// todo: Use these units A/R, or don't.
/// Inner values: The value.
#[derive(Clone, Copy, PartialEq)]
pub enum Length {
    // M(f64),
    // Km(f64),
    // Pc(f64),
    // Kpc(f64),
    M,
    Km,
    Pc,
    Kpc,
}

#[derive(Clone, Copy, PartialEq)]
pub enum Unit {
    Length((Length, f64)),
}

/// Measure dist is in Kpc; i.e. the distance from earth.
pub fn arcsec_to_km(val_arcsec: f64, measure_dist: f64) -> f64 {
    val_arcsec * measure_dist * ARCSEC_CONV_FACTOR
}
