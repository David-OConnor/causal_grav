//! For constructing and working with Gaussian rings.

use lin_alg::f64::Vec3;

#[derive(Debug)]
pub struct GaussianShell {
    pub center: Vec3,
    pub radius: f64,
    /// a, b, and c are the gaussian constants. A scales in magnitude
    pub a: f64,
    /// C controls width. Higher C means higher width. (Symmetric)
    pub c: f64,
}

impl GaussianShell {
    pub fn value(&self, posit: Vec3) -> f64 {
        // Find the effective "x" value, as used in the 1D gaussian function. This X scale is centered
        // on the radius, and is on an imaginary line along its radial direction.

        // Calculate the slope of the line between the point in question, and the sphere center.
        // let slope = (posit.y - self.center.y) / (posit.x - self.center.x);
        // The point on the circle is the intersection between this line and the circle.
        let v = posit - self.center;
        let mag_v = v.magnitude();

        // todo: QC this and related.
        let closest_pt_on_shell = self.center + v / mag_v * self.radius;

        let x_1d = (closest_pt_on_shell - posit).magnitude();
        // Note: We don't need to worry about this being negative, since the gaussian
        // is symmetric around 0.

        // Calculate the 1D gaussian value, using its standard definition.
        let b = 0.;

        self.a * (-(x_1d - b).powi(2) / (2. * self.c.powi(2))).exp()
    }
}
