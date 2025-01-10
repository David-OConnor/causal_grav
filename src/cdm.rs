//! Related to Cold Dark Matter (CDM)

/// Generate a Berkert Halo. Generally gives good fites to rotation curves.
/// rho_0 is the central density. r_core is the core radius.
pub fn density_burkert(r: f64, rho_0: f64, r_core: f64) -> f64 {
    rho_0 * r_core.powi(3) / ((r + r_core) * (r.powi(2) + r_core.powi(2)))
}

/// Generate a NFW Halo.
/// rho_s is the characteristic density. r_csis the characteristic scale.
pub fn density_nfw(r: f64, rho_s: f64, r_s: f64) -> f64 {
    rho_s / ((r / r_s) * (1. + r/r_s).powi(2))
}