//! Gravitoelectromagnetism

use lin_alg::f64::Vec3;

pub struct FourPotential {
    /// ϕ
    pub scaler: f64,
    /// A
    pub vector: Vec3,
}

impl FourPotential {
    /// E = -∇ϕ - ∂A/∂t
    pub fn elec_field(&self) -> Vec3 {
        Vec3::new_zero() // todo
    }

    /// B = ∇ × A
    pub fn mag_field(&self) -> Vec3 {
        Vec3::new_zero() // todo
    }
}

// https://en.wikipedia.org/wiki/Gravitoelectromagnetism
//todo A/R the 4 GEM equations.
