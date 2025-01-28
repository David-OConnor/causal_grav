//! Holding ground for our unused ray code.

use std::f64::consts::TAU;

const MAX_RAY_DIST: f64 = 30.; // todo: Adjust this approach A/R.

pub const RAY_SIZE: f32 = 0.005;
pub const RAY_SHINYNESS: f32 = 2.;
pub const RAY_COLORS: [Color; 5] = [
    (1., 1.0, 0.2),
    (0.5, 0.7, 1.),
    (0., 0., 0.9),
    (0.9, 0.3, 0.2),
    (1., 0.0, 1.0),
];

// This cubed is d-volume
const RAY_SAMPLE_WIDTH: f64 = 0.8;

// We use this to calculate divergence of ray density, numerically.
const DX_RAY_GRADIENT: f64 = RAY_SAMPLE_WIDTH;

impl State {
    /// Remove rays that are far from the area of interest, for performance reasons.
    /// todo: This is crude; improve it to be more efficient.
    fn remove_far_rays(&mut self) {
        self.rays.retain(|ray| {
            ray.posit.x <= MAX_RAY_DIST
                && ray.posit.y <= MAX_RAY_DIST
                && ray.posit.z <= MAX_RAY_DIST
        });
    }
}

impl Body {
    /// Generate a ray traveling in a random direction.
    pub fn create_ray(&self, source_id: usize) -> GravRay {
        let mut rng = rand::thread_rng();

        let unit_vec = util::random_unit_vec(&rng);

        let unit_vec = Vec3::new(x, y, z);
        GravRay {
            source_id: source_id,
            posit: self.posit,
            vel: unit_vec * C,
            src_mass: self.mass,
        }
    }
}

/// We are treating the ray density gradient as a proxy for gravitational potential.
/// Calculate it numberically.
fn density_gradient(
    posit: Vec3,
    rays: &[GravRay],
    shells: &[GravShell],
    source_id: usize,
    dt: f64,
    gauss_c: f64,
) -> Vec3 {
    // let rect_this = SampleRect::new(posit, RAY_SAMPLE_WIDTH);

    let rect_x_prev = SampleRect::new(posit - Vec3::new(DX_RAY_GRADIENT, 0., 0.), RAY_SAMPLE_WIDTH);
    let rect_x_next = SampleRect::new(posit + Vec3::new(DX_RAY_GRADIENT, 0., 0.), RAY_SAMPLE_WIDTH);

    let rect_y_prev = SampleRect::new(posit - Vec3::new(0., DX_RAY_GRADIENT, 0.), RAY_SAMPLE_WIDTH);
    let rect_y_next = SampleRect::new(posit + Vec3::new(0., DX_RAY_GRADIENT, 0.), RAY_SAMPLE_WIDTH);

    let rect_z_prev = SampleRect::new(posit - Vec3::new(0., 0., DX_RAY_GRADIENT), RAY_SAMPLE_WIDTH);
    let rect_z_next = SampleRect::new(posit + Vec3::new(0., 0., DX_RAY_GRADIENT), RAY_SAMPLE_WIDTH);

    // let properties_this = rect_this.measure_properties(rays);

    let properties_x_prev = rect_x_prev.measure_properties(rays, shells, source_id, gauss_c);
    let properties_x_next = rect_x_next.measure_properties(rays, shells, source_id, gauss_c);
    let properties_y_prev = rect_y_prev.measure_properties(rays, shells, source_id, gauss_c);
    let properties_y_next = rect_y_next.measure_properties(rays, shells, source_id, gauss_c);
    let properties_z_prev = rect_z_prev.measure_properties(rays, shells, source_id, gauss_c);
    let properties_z_next = rect_z_next.measure_properties(rays, shells, source_id, gauss_c);

    Vec3::new(
        (properties_x_next.ray_density - properties_x_prev.ray_density) / 2.,
        (properties_y_next.ray_density - properties_y_prev.ray_density) / 2.,
        (properties_z_next.ray_density - properties_z_prev.ray_density) / 2.,
    )
    //
    // // todo: QC this etc
    // // Potential
    // let V_this = properties_this.ray_density;
    //
    // // todo: Is this the divergence property, or should we take the gradient using finite difference
    // // todo from two sample rect locations?
    // let gradient = Vec3::new_zero();
}

/// A rectangular prism for sampling properties (generally an "infinitessimal" volume.
struct SampleRect {
    pub start: Vec3,
    pub end: Vec3,
}

impl SampleRect {
    /// Create a rectangle of a giFven size, centered on a position
    /// todo: Accept d-volume instead of width?
    pub fn new(posit: Vec3, width: f64) -> Self {
        let w_div2 = width / 2.;
        Self {
            start: Vec3::new(posit.x - w_div2, posit.y - w_div2, posit.z - w_div2),
            end: Vec3::new(posit.x + w_div2, posit.y + w_div2, posit.z + w_div2),
        }
    }

    fn measure_properties(
        &self,
        rays: &[GravRay],
        shells: &[GravShell],
        source_id: usize,
        shell_c: f64,
    ) -> SampleProperties {
        let volume = {
            let dist_x = self.end.x - self.start.x;
            let dist_y = self.end.y - self.start.y;
            let dist_z = self.end.z - self.start.z;

            (dist_x.powi(2) + dist_y.powi(2) + dist_z.powi(2)).sqrt()
        };

        let mut ray_value = 0.;
        let mut ray_vel_sum = Vec3::new_zero();

        for ray in rays {
            // A method to avoid self-interaction.
            if ray.source_id == source_id {
                continue;
            }

            if ray.posit.x >= self.start.x
                && ray.posit.x <= self.end.x
                && ray.posit.y >= self.start.y
                && ray.posit.y <= self.end.y
                && ray.posit.z >= self.start.z
                && ray.posit.z <= self.end.z
            {
                ray_value += ray.src_mass;
                ray_vel_sum += ray.vel;
            }
        }

        // println!("\nVel sum: {:?}", vel_sum);
        // println!("Ray val: {ray_value}\n");

        // note: This is undoing the conversion to the rectangle.
        let sample_center = Vec3::new(
            (self.end.x + self.start.x) / 2.,
            (self.end.y + self.start.y) / 2.,
            (self.end.z + self.start.z) / 2.,
        );

        // let mut shell_inner = None;
        // let mut shell_outer = None;

        let acc_shell = accel::calc_acc_shell(shells, sample_center, source_id, shell_c);

        // todo: To calculate div and curl, we need multiple sets of rays.

        SampleProperties {
            ray_density: ray_value / volume,
            ray_net_direction: ray_vel_sum.to_normalized() * -1.,
            acc_shell,
            div: 0.,
            curl: 0.,
        }
    }
}

#[derive(Debug)]
/// Our model "particle" that travels outward from a source
struct GravRay {
    /// A quick and dirty way to prevent a body from interacting with itself.
    source_id: usize,
    posit: Vec3,
    /// The magnitude corresponds to source mass. Direction is outward
    /// at any angle from the source.
    vel: Vec3,
    // todo: We can assume there is no acc on this, right?
    /// We are currently applying a weight based on source mass; we could alternatively have more
    /// massive objects output more rays. This method here is likely more computationally efficient.
    src_mass: f64,
}


impl GravShell {
    /// Determine if the shell surface intersects a rectangular box.
    pub fn intersects_rect(&self, rect: &SampleRect) -> bool {
        // Iterate over all the vertices of the rectangular prism.
        let box_vertices = [
            Vec3::new(rect.start.x, rect.start.y, rect.start.z),
            Vec3::new(rect.start.x, rect.start.y, rect.end.z),
            Vec3::new(rect.start.x, rect.end.y, rect.start.z),
            Vec3::new(rect.start.x, rect.end.y, rect.end.z),
            Vec3::new(rect.end.x, rect.start.y, rect.start.z),
            Vec3::new(rect.end.x, rect.start.y, rect.end.z),
            Vec3::new(rect.end.x, rect.end.y, rect.start.z),
            Vec3::new(rect.end.x, rect.end.y, rect.end.z),
        ];

        for vertex in &box_vertices {
            let distance = (self.center.x - vertex.x).powi(2)
                + (self.center.y - vertex.y).powi(2)
                + (self.center.z - vertex.z).powi(2);

            let radius_squared = self.radius.powi(2);

            // Check if the distance to the vertex is exactly the radius.
            if (distance - radius_squared).abs() < f64::EPSILON {
                return true;
            }
        }

        false
    }
}



#[derive(Debug)]
struct SampleProperties {
    /// Is this an analog for V (potential) ?
    ray_density: f64,
    /// A unit vec, from the weighted average of ray velocities.
    ray_net_direction: Vec3,
    acc_shell: Vec3,
    /// Divergence
    div: f64,
    curl: f64,
}

/// Calculate the force acting on a body, given the local environment of gravity rays around it.
pub fn acc_rays(
    body: &mut Body,
    rays: &[GravRay],
    shells: &[GravShell],
    source_id: usize,
    shell_c: f64,
) -> Vec3 {
    // let ray_density_grad = density_gradient(body.posit, rays);

    let rect = SampleRect::new(body.posit, RAY_SAMPLE_WIDTH);
    let properties = rect.measure_properties(rays, shells, source_id, shell_c);

    // println!("Prop: {:?}", properties);

    // todo: The rate constant must take into account the angle of the rays; dividing by dt
    // todo is likely only to work in the limit of infinitely small d_theta for ray emission, etc.
    // todo: Divide by a rate constant.
    let rate_const = 460.;
    // let rate_const = 1. / dt; // todo: we can pre-calc this, but not a big deal.

    properties.ray_net_direction * properties.ray_density * rate_const
}
