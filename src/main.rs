#![allow(non_snake_case)]

use std::f64::consts::TAU;

use lin_alg::{f32::Vec3 as Vec3f32, f64::Vec3};
use rand::Rng;

use crate::{
    gaussian::GaussianShell,
    playback::{vec_to_f32, SnapShot},
    render::render,
};
// Add this to use the random number generator

mod gaussian;
mod playback;
mod render;
mod ui;
// Shower thought, from looking at this from a first person view: View things from the body's perspective.
// Can you make of it something like that?

// todo: What if you have total electric charge different from a whole number? Misc elec "particles"
// todo zipping around in the either, to be capture.

// todo: Next, try having the electrons add up to more than 1 charge.
// todo: Try to quantify the elec density, so you can compare it to Schrodinger.

const M_ELEC: f64 = 1.;
const Q_ELEC: f64 = -1.;

const MAX_RAY_DIST: f64 = 30.; // todo: Adjust this approach A/R.
const MAX_SHELL_R: f64 = 300.; // todo: Adjust this approach A/R.

const SNAPSHOT_RATIO: usize = 10;

const C: f64 = 10.;

// todo: A/R
// This cubed is d-volume
const RAY_SAMPLE_WIDTH: f64 = 0.8;

// We use this to calculate divergence of ray density, numerically.
const DX_RAY_GRADIENT: f64 = RAY_SAMPLE_WIDTH;

// Don't calculate force if particles are farther than this from each other. Computation saver.
// const MAX_DIST: f64 = 0.1;

pub struct Config {
    num_timesteps: usize,
    dt_integration: f64,
    // dt_pulse: f64,
    num_rays_per_iter: usize,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            // num_timesteps: 2_000,
            num_timesteps: 20_000,
            // dt_integration: 0.001,
            dt_integration: 0.0001,
            num_rays_per_iter: 200,
        }
    }
}

#[derive(Default)]
pub struct StateUi {
    snapshot_selected: usize,
}

#[derive(Default)]
struct State {
    config: Config,
    ui: StateUi,
    bodies: Vec<Body>,
    rays: Vec<GravRay>,
    shells: Vec<GravShell>,
    snapshots: Vec<SnapShot>,
}

impl State {
    /// Remove rays that are far from the area of interest, for performance reasons.
    /// todo: This is crude; improve it to be more efficient.
    fn remove_far_rays(&mut self) {
        self.rays.retain(|ray| {
            ray.posit.x <= MAX_RAY_DIST
                && ray.posit.y <= MAX_RAY_DIST
                && ray.posit.z <= MAX_RAY_DIST
        });

        self.shells.retain(|shell| shell.radius <= MAX_SHELL_R);
    }

    fn take_snapshot(&mut self, time: usize) {
        let mut accs = Vec::new();
        let mut accs_shell = Vec::new();

        for (i, body) in self.bodies.iter_mut().enumerate() {
            let acc = accel(
                body,
                &self.rays,
                &self.shells,
                i,
                self.config.dt_integration,
            );
            let acc_shell = accel_shells(
                body,
                &self.rays,
                &self.shells,
                i,
                self.config.dt_integration,
            );

            accs.push(acc);
            accs_shell.push(acc_shell);
        }

        self.snapshots.push(SnapShot {
            time,
            body_posits: self.bodies.iter().map(|b| vec_to_f32(b.posit)).collect(),
            V_at_bodies: self
                .bodies
                .iter()
                .map(|b| vec_to_f32(b.V_acting_on))
                .collect(),
            acc_at_bodies: accs_shell.iter().map(|a| vec_to_f32(*a)).collect(),
            rays: self
                .rays
                .iter()
                .map(|r| {
                    (
                        Vec3f32::new(r.posit.x as f32, r.posit.y as f32, r.posit.z as f32),
                        r.emitter_id,
                    )
                })
                .collect(),
            // shells: self.shells.clone(),
            shells: Vec::new(),
        })
    }
}

struct Body {
    posit: Vec3,
    vel: Vec3,
    accel: Vec3,
    mass: f64,
    /// We use this for debugging, testing etc. It goes in the snapshots.
    V_acting_on: Vec3,
}

impl Body {
    /// Generate a ray traveling in a random direction.
    pub fn create_ray(&self, emitter_id: usize) -> GravRay {
        let mut rng = rand::thread_rng();

        // todo: You could also generate a random quaternion.

        let theta = rng.gen_range(0.0..TAU); // Random angle in [0, 𝜏)
        let cos_phi = rng.gen_range(-1.0..1.0); // Uniform in [-1, 1]
        let sin_phi = (1.0_f64 - cos_phi * cos_phi).sqrt();

        let x = sin_phi * theta.cos();
        let y = sin_phi * theta.sin();
        let z = cos_phi;
        let unit_vec = Vec3::new(x, y, z);
        GravRay {
            emitter_id,
            posit: self.posit,
            vel: unit_vec * C,
            src_mass: self.mass,
        }
    }

    /// Generate a shell traveling outward.
    pub fn create_shell(&self, emitter_id: usize) -> GravShell {
        GravShell {
            emitter_id,
            center: self.posit,
            radius: 0.,
            src_mass: self.mass,
        }
    }
}

// // (radius, charge_in_radius)
// fn density(
//     elec_posits: &[Vec3f32],
//     q_per_elec: f64,
//     center_ref: Vec3f32,
// ) -> Vec<((f32, f32), f32)> {
//     // Must be ascending radii for the below code to work.
//     let per_elec = q_per_elec as f32;
//
//     // Equidistant ish, for now.
//     let mut result = vec![
//         ((0.0, 0.2), 0.),
//         ((0.2, 0.4), 0.),
//         ((0.4, 0.6), 0.),
//         ((0.8, 1.), 0.),
//         ((1., 1.2), 0.),
//         ((1.2, 1.4), 0.),
//         ((1.4, 1.6), 0.),
//         ((1.6, 1.8), 0.),
//         ((1.8, 2.0), 0.),
//         ((2.0, 9999.), 0.),
//     ];
//
//     for (r_bound, q) in &mut result {
//         for posit in elec_posits {
//             let mag = (center_ref - *posit).magnitude();
//             if mag < r_bound.1 && mag > r_bound.0 {
//                 *q += per_elec;
//             }
//         }
//     }
//
//     result
// }

/// We are treating the ray density gradient as a proxy for gravitational potential.
/// Calculate it numberically.
fn density_gradient(
    posit: Vec3,
    rays: &[GravRay],
    shells: &[GravShell],
    emitter_id: usize,
) -> Vec3 {
    // let rect_this = SampleRect::new(posit, RAY_SAMPLE_WIDTH);

    let rect_x_prev = SampleRect::new(posit - Vec3::new(DX_RAY_GRADIENT, 0., 0.), RAY_SAMPLE_WIDTH);
    let rect_x_next = SampleRect::new(posit + Vec3::new(DX_RAY_GRADIENT, 0., 0.), RAY_SAMPLE_WIDTH);

    let rect_y_prev = SampleRect::new(posit - Vec3::new(0., DX_RAY_GRADIENT, 0.), RAY_SAMPLE_WIDTH);
    let rect_y_next = SampleRect::new(posit + Vec3::new(0., DX_RAY_GRADIENT, 0.), RAY_SAMPLE_WIDTH);

    let rect_z_prev = SampleRect::new(posit - Vec3::new(0., 0., DX_RAY_GRADIENT), RAY_SAMPLE_WIDTH);
    let rect_z_next = SampleRect::new(posit + Vec3::new(0., 0., DX_RAY_GRADIENT), RAY_SAMPLE_WIDTH);

    // let properties_this = rect_this.measure_properties(rays);

    let properties_x_prev = rect_x_prev.measure_properties(rays, shells, emitter_id);
    let properties_x_next = rect_x_next.measure_properties(rays, shells, emitter_id);
    let properties_y_prev = rect_y_prev.measure_properties(rays, shells, emitter_id);
    let properties_y_next = rect_y_next.measure_properties(rays, shells, emitter_id);
    let properties_z_prev = rect_z_prev.measure_properties(rays, shells, emitter_id);
    let properties_z_next = rect_z_next.measure_properties(rays, shells, emitter_id);

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

/// Calculate the force acting on a body, given the local environment of gravity rays around it.
fn accel(
    body: &mut Body,
    rays: &[GravRay],
    shells: &[GravShell],
    emitter_id: usize,
    dt: f64,
) -> Vec3 {
    // let ray_density_grad = density_gradient(body.posit, rays);

    let rect = SampleRect::new(body.posit, RAY_SAMPLE_WIDTH);
    let properties = rect.measure_properties(rays, shells, emitter_id);

    // println!("Prop: {:?}", properties);

    // todo: The rate constant must take into account the angle of the rays; dividing by dt
    // todo is likely only to work in the limit of infinitely small d_theta for ray emission, etc.
    // todo: Divide by a rate constant.
    let rate_const = 460.;
    // let rate_const = 1. / dt; // todo: we can pre-calc this, but not a big deal.

    properties.ray_net_direction * properties.ray_density * rate_const

    // todo: We need to take into account the destination body's mass, not just
    // todo for inertia, but for attraction... right?

    // todo: We may be missing part of the puzzle.
    // let force = gradient * body.mass;

    // a = F / m
    // force / body.mass

    // body.V_acting_on = ray_density_grad;

    // ray_density_grad
}

/// Calculate the force acting on a body, given the local environment of gravity shells intersecting it.
fn accel_shells(
    body: &mut Body,
    rays: &[GravRay],
    shells: &[GravShell],
    emitter_id: usize,
    dt: f64,
) -> Vec3 {
    let rect = SampleRect::new(body.posit, RAY_SAMPLE_WIDTH);
    let properties = rect.measure_properties(rays, shells, emitter_id);

    // println!("Prop: {:?}", properties);

    properties.acc_shell
}

// todo: Put back once you have an accel computation.
// fn integrate_rk4(bodies: &mut [Body], dt: f64) {
//     for body in bodies.iter_mut() {
//         // Step 1: Calculate the k-values for position and velocity
//         let k1_v = body.accel * dt;
//         let k1_posit = body.vel * dt;
//
//         let k2_v = compute_acceleration(body.posit + k1_posit * 0.5) * dt;
//         let k2_posit = (body.vel + k1_v * 0.5) * dt;
//
//         let k3_v = compute_acceleration(body.posit + k2_posit * 0.5) * dt;
//         let k3_posit = (body.vel + k2_v * 0.5) * dt;
//
//         let k4_v = compute_acceleration(body.posit + k3_posit) * dt;
//         let k4_posit = (body.vel + k3_v) * dt;
//
//         // Step 2: Update position and velocity using weighted average of k-values
//         body.vel += (k1_v + k2_v * 2. + k3_v * 2. + k4_v) / 6.;
//         body.posit += (k1_posit + k2_posit * 2. + k3_posit * 2. + k4_posit) / 6.;
//     }
// }

// todo: DRY
// todo: Unless we want to apply accel etc to these rays, RK4 here is not required; accel
// todo is always 0, and v is always c.
fn integrate_rk4_ray(rays: &mut [GravRay], dt: f64) {
    // todo: Pending further exporation, no grav accel on rays.
    let a = Vec3::new_zero();

    for body in rays.iter_mut() {
        // Step 1: Calculate the k-values for position and velocity
        let k1_v = a;
        let k1_posit = body.vel * dt;

        let k2_v = a;
        let k2_posit = (body.vel + k1_v * 0.5) * dt;

        let k3_v = a;
        let k3_posit = (body.vel + k2_v * 0.5) * dt;

        let k4_v = a;
        let k4_posit = (body.vel + k3_v) * dt;

        // Step 2: Update position and velocity using weighted average of k-values
        // body.vel += (k1_v + k2_v * 2. + k3_v * 2. + k4_v) / 6.;

        // Vel is constant; C.
        body.posit += (k1_posit + k2_posit * 2. + k3_posit * 2. + k4_posit) / 6.;
    }
}

/// A rectangular prism for sampling properties (generally an "infinitessimal" volume.
struct SampleRect {
    pub start: Vec3,
    pub end: Vec3,
}

impl SampleRect {
    /// Create a rectangle of a given size, centered on a position
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
        emitter_id: usize,
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
            if ray.emitter_id == emitter_id {
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

        let sample_center = Vec3::new(
            (self.end.x + self.start.x) / 2.,
            (self.end.y + self.start.y) / 2.,
            (self.end.z + self.start.z) / 2.,
        );

        // let mut shell_inner = None;
        // let mut shell_outer = None;

        let mut shell_value = 0.;
        let mut shell_vel_sum = Vec3::new_zero();

        // todo: Once you have more than one body acting on a target, you need to change this, so you get
        // todo exactly 0 or 1 shells per other body.
        for shell in shells {
            if shell.emitter_id == emitter_id {
                continue;
            }

            let gauss = GaussianShell {
                center: shell.center,
                radius: shell.radius,
                a: 1.,
                c: 1., // todo: Exper with a and c.
            };
            shell_value += shell.src_mass * gauss.value(sample_center);
            // todo: QC what the acc dir from the shell is.
            let shell_acc_dir = (sample_center - shell.center).to_normalized();
            shell_vel_sum += shell_acc_dir * shell_value; // todo: QC order

            // if shell.intersects_rect(self) {
            //     println!("Intersects: {:?}", shell);
            //     acc_shell + (shell.center - center) * shell.src_mass / shell.radius.powi(2);
            // }
        }
        // todo: To calculate div and curl, we need multiple sets of rays.

        SampleProperties {
            ray_density: ray_value / volume,
            ray_net_direction: ray_vel_sum.to_normalized() * -1.,
            acc_shell: shell_vel_sum,
            div: 0.,
            curl: 0.,
        }
    }
}

#[derive(Debug)]
/// Our model "particle" that travels outward from a source
struct GravRay {
    /// A quick and dirty way to prevent a body from interacting with itself.
    emitter_id: usize,
    posit: Vec3,
    /// The magnitude corresponds to source mass. Direction is outward
    /// at any angle from the source.
    vel: Vec3,
    // todo: We can assume there is no acc on this, right?
    /// We are currently applying a weight based on source mass; we could alternatively have more
    /// massive objects output more rays. This method here is likely more computationally efficient.
    src_mass: f64,
}

#[derive(Debug, Clone)]
struct GravShell {
    emitter_id: usize,
    center: Vec3,
    radius: f64,
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

/// Entry point for computation; rename A/R.
fn build(state: &mut State) {
    // todo temp
    // let mut acc = Vec3::new_zero();
    // let mut acc_shell = Vec3::new_zero();
    // let mut counter = 0;

    for t in 0..state.config.num_timesteps {
        // Create a new set of rays.
        for (id, body) in state.bodies.iter().enumerate() {
            // todo temp RM
            for _ in 0..state.config.num_rays_per_iter {
                // state.rays.push(body.create_ray(id));
            }
            state.shells.push(body.create_shell(id));
        }

        state.remove_far_rays();

        // Update ray propogation
        integrate_rk4_ray(&mut state.rays, state.config.dt_integration);

        for shell in &mut state.shells {
            shell.radius += C * state.config.dt_integration;
        }

        // Update body motion.
        // integrate_rk4(&mut state.bodies, state.config.dt_integration);

        // Save the current state to a snapshot, for later playback.
        // Note: This can use a substantial amount of memory.

        if t % SNAPSHOT_RATIO == 0 {
            state.take_snapshot(t / SNAPSHOT_RATIO);
        }

        //     // todo temp
        //     if t + 30 >= state.config.num_timesteps {
        //         counter += 1;
        //         acc += accel(
        //             &mut state.bodies[1],
        //             &state.rays,
        //             &state.shells,
        //             1,
        //             state.config.dt_integration,
        //         );
        //         acc_shell += accel_shells(
        //             &mut state.bodies[1],
        //             &state.rays,
        //             &state.shells,
        //             1,
        //             state.config.dt_integration,
        //         );
        //     }
    }
    //
    // println!("SHELLS: {:?}", state.shells.len());
    //
    // // todo temp
    // acc = acc / counter as f64;
    // acc_shell = acc_shell / counter as f64;
    // println!("Acc net: {:?}. Shell: {:?}", acc, acc_shell);
}

fn main() {
    println!("Building snapshots...");
    let mut state = State::default();

    state.bodies = vec![
        Body {
            posit: Vec3::new_zero(),
            vel: Vec3::new_zero(),
            accel: Vec3::new_zero(),
            mass: 1.,
            V_acting_on: Vec3::new_zero(),
        },
        Body {
            posit: Vec3::new(3., 0., 0.),
            // vel: Vec3::new_zero(),
            vel: Vec3::new(0.01, 0., 0.), // todo tmep
            accel: Vec3::new_zero(),
            mass: 1.,
            V_acting_on: Vec3::new_zero(),
        },
    ];

    build(&mut state);

    // todo: Temp. Testing our approach
    // let acc = accel(&mut state.bodies[1], &state.rays, 1);

    // todo: We must discount etc the body's own rays. ANd possibly consider their
    // todo significance.

    // println!("Acc: {:?}", acc);

    println!("Complete. Rendering.");

    render(state);
}
