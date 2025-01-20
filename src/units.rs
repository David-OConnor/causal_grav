//! Contains unit conversions, and a definition of the base units used throughout this program.

#![allow(unused)]

use std::f64::consts::TAU;

// We use this to convert angle to length, when multiplied by distance.
// Cache, vice using `.to_radians`.
pub const ARCSEC_CONV_FACTOR: f64 = TAU / (360. * 3_600.);

// Let's base units for our sim:

// Basic:
// Dist: kpc
// Time: Myr: (million years)
// Mass: M☉ (or M☉ × 10^10)?
//
// Derived:
// Velocity: kpc / Myr
// Accel: kpc / Myr^2
// Force: M☉ × kpc / Myr^2
// Energy: M☉ × (kpc / Mry)^2

// An alternate system:
// pc, M☉, km/s  (I've seen G written this way, but it uses both pc and km for distance.)

const M_PER_KPC: f64 = 3.085677581e19; // 1 KPC in meters. m/kpc
const S_PER_MYR: f64 = 3.15576e13; // 1 Myr in seconds. s/Mry
const KG_PER_SOLAR: f64 = 1.988_416e30; // 1 Solar mass in kg. kg/M☉

// If we end up using these scaled units instead of the native ones above.
// pub const SOLAR_MASS_E10: f64 = SOLAR_MASS/ 1.0e10; // 1 M☉ × 10^10 in kg

const G_SI: f64 = 6.67_430e-11; // m^3 / (kg s^2) or N m^2 / kg^2
                                // const G_PC_M_KMS: f64 = 4.3009172706e-3; // Wikipedia

// gravitational constant, in KPC^3 / (M☉ × Myr^2)
const KPC_PER_M: f64 = 1. / M_PER_KPC; // 1 m in KPC. Kpc/m
const MYR_PER_S: f64 = 1. / S_PER_MYR; // 1 s in Myr. Myr/s
const SOLAR_PER_KG: f64 = 1. / KG_PER_SOLAR; // 1 kg in M☉. M☉/kg

// We use this to convert rotational velocity data (km/s) to our base units of kpc/Myr
pub const KPC_MYR_PER_KM_S: f64 = KPC_PER_M * 1_000. / MYR_PER_S; // ~1.022e-3

// Confirmed: https://www.wolframalpha.com/input?i=gravitational+constant+to+kiloparsec%5E3%2F%28solar+mass*megayear%5E2%29&assumption=%22UnitClash%22+-%3E+%7B%22gravitational+constant%22%2C+%7B%22GravitationalConstant%22%2C+%22dflt%22%7D%7D&assumption=%7B%22F%22%2C+%22UnitsConversion2%22%2C+%22fromValue%22%7D+-%3E%221%22
pub const G: f64 =
    G_SI * (KPC_PER_M * KPC_PER_M * KPC_PER_M / (SOLAR_PER_KG * MYR_PER_S * MYR_PER_S)); // 4.4984e-12

// Let's see if we can get 4.3e-3 pc⋅M⊙−1⋅(km/s)2 (The above Wikipedia const) from SI units to validated
// our approach.

// todo: Consider a natural unit set where G = 1.

pub const A0_MOND: f64 = 1.2e-10 * KPC_PER_M / (MYR_PER_S * MYR_PER_S); // 3.87e-3

// Note: Setting this too high is problematic.
// todo: Maybe a different time unit?
// const C: f64 = 9.72e-12; // Rough; kpc/s^2.
// const C: f64 = 40.;
pub const C: f64 = 306.4; // KPC/Myr

// G_SI × 1/3.0856e16 × c^2 / SOLAR_MASS_INV
// = 4.45e-3. This checks out. Our approach of using inverses, and preserving multiplication/division order
// is validated; the above G should work.

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
