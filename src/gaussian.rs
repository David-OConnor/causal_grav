//! For constructing and working with Gaussian rings.
//!
//! Notes:
//! - As C increases, the distribution becomes uniform, at the cost of more sloping edges. (Sloppier falloff)
//! - As the gaussian spacing becomes smaller, the distribution becomes uniform, at the cost of more computation required.
//! - As C decreases and/or spacing becomes larger, the distribution becomes wavier; less uniform.
//! - It appears that amplitude (`a`, and/or our gravity amplitude scaler) is not relevant to the shape; ignore it.
//! - We can select a constant, `COEFF_C`, to multiply spacing by, to balance uniformity vs falloff.
//! - We must scale amplitude down, according to a rule I haven't determined that depends only on `COEFF_C`.
//!
//! todo: Update: It appears in your calcs below and assumptions above, you may be missing part of the effect
//! of spacing on uniformity.
//!
//! Recommendation:
//!
//!
//! todo: Is there a function more suitable for this task than gaussians? Maybe boxes. This could potentially
//! todo: have perfectly sharp edges, but I can imagine some negative effect from these combining from different
//! todo: directions.

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

        // Attempt at a performance saver. // todo: QC that this doesn't introduce problems.
        if x_1d > self.c * 5. {
            // todo: Adjust a/r.
            // let test_r = self.a * (-x_1d.powi(2) / (2. * self.c.powi(2))).exp();
            // if test_r > 0.00001 {
            //     println!("Test R: {:?}", test_r);
            // }
            return 0.;
        }

        // Calculate the 1D gaussian value, using its standard definition.
        self.a * (-x_1d.powi(2) / (2. * self.c.powi(2))).exp()
    }
}
