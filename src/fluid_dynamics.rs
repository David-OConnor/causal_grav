//! Code related to fluid dyanmics models, vice point masses.

use lin_alg::f64::Vec3;

/// Smoothed-particle hydrodynamics
pub struct SphPoint {
    pub posit: Vec3,
    pub vel: Vec3,
    pub accel: Vec3,
    /// Local density
    pub density: f64,
    pub pressure: f64,
    pub mass: f64,
    /// Smoothing length (influence radius)
    pub smoothing_length: f64,
    /// Indices of neighboring particles
    pub neighbors: Vec<usize>,
    // Potentially optional fields below
    /// For thermal dynamics.
    pub temp: f64,
    /// For multiphase simulations
    pub color: i8, // todo?
    /// To differentiate fluid and boundary particles.
    pub boundary_flags: i8,
    /// If using particle-specific viscosity.
    pub viscosity: f64,
}

// todo: Mesh methods.
