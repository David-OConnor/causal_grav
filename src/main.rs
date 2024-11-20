use std::f64::{consts::TAU, MAX};

use lin_alg::f64::Vec3;
use rand::Rng;

use crate::{playback::SnapShot, render::render};
// Add this to use the random number generator

mod playback;
mod render;
mod ui;
// todo: What if you have total electric charge different from a whole number? Misc elec "particles"
// todo zipping around in the either, to be capture.

// todo: Next, try having the electrons add up to more than 1 charge.
// todo: Try to quantify the elec density, so you can compare it to Schrodinger.

const M_ELEC: f64 = 1.;
const Q_ELEC: f64 = -1.;

// If two particles are closer to each other than this, don't count accel, or cap it.
const MIN_DIST: f64 = 0.1;

const MAX_RAY_DIST: f64 = 100.; // todo: Adjust this approach A/R.

// Don't calculate force if particles are farther than this from each other. Computation saver.
// const MAX_DIST: f64 = 0.1;

pub struct Config {
    num_timesteps: usize,
    dt_integration: f64,
    dt_pulse: f64,
    num_rays_per_iter: usize,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            num_timesteps: 1_000,
            dt_integration: 0.01,
            dt_pulse: 0.01,
            num_rays_per_iter: 40,
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
    snapshots: Vec<SnapShot>,
}

impl State {
    /// Remove rays that are far from the area of interest, for performance reasons.
    /// todo: This is crude; improve it to be more efficient.
    fn remove_far_rays(&mut self) {
        let mut removed_rays = Vec::new();

        for (i, ray) in self.rays.iter().enumerate() {
            if ray.posit.x > MAX_RAY_DIST
                || ray.posit.y > MAX_RAY_DIST
                || ray.posit.z > MAX_RAY_DIST
            {
                removed_rays.push(i);
                continue;
            }
        }

        for i in removed_rays {
            self.rays.remove(i);
        }
    }
}

struct Body {
    posit: Vec3,
    vel: Vec3,
    accel: Vec3,
    mass: f64,
}

impl Body {
    /// Generate a ray in a random direction.
    pub fn create_ray(&self) -> GravRay {
        let mut rng = rand::thread_rng();

        let theta = rng.gen_range(0.0..TAU); // Random angle in [0, 𝜏)
        let phi = rng.gen_range(0.0..TAU / 2.); // Random angle in [0, 𝜏/2]

        let x = phi.sin() * theta.cos();
        let y = phi.sin() * theta.sin();
        let z = phi.cos();
        let unit_vec = Vec3::new(x, y, z);

        GravRay {
            posit: self.posit,
            vel: unit_vec * self.mass,
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

/// Calculate the Coulomb or gravitational acceleration on a particle, from a single other particle.
fn accel(
    posit_acted_on: Vec3,
    posit_actor: Vec3,
    q_acted_on: f64,
    q_actor: f64,
    mass_acted_on: f64,
    nuc_actor: bool,
) -> Vec3 {
    let posit_diff = posit_acted_on - posit_actor;
    let dist = posit_diff.magnitude();

    let posit_diff_unit = posit_diff / dist;

    // Note: We factor out KC * Q_PROT, since they are present in every calculation.
    // Calculate the Coulomb force between nuclei.

    let f_mag = q_acted_on * q_actor / dist.powi(2);

    // todo: Solving the nuc-toss-out problem.
    if nuc_actor && dist < 0.1 {
        // println!("Min dist: {:?} Force: {:?}", dist, f_mag);
        return Vec3::new_zero();
    }

    if dist < MIN_DIST {
        println!("Min dist: {:?} Force: {:?}", dist, f_mag);
        return Vec3::new_zero();
    }

    posit_diff_unit * f_mag / mass_acted_on
}

fn integrate_rk4(bodies: &mut [Body], dt: f64) {
    for body in bodies.iter_mut() {
        // Step 1: Calculate the k-values for position and velocity
        let k1_v = body.accel * dt;
        let k1_posit = body.vel * dt;

        let k2_v = body.accel * dt;
        let k2_posit = (body.vel + k1_v * 0.5) * dt;

        let k3_v = body.accel * dt;
        let k3_posit = (body.vel + k2_v * 0.5) * dt;

        let k4_v = body.accel * dt;
        let k4_posit = (body.vel + k3_v) * dt;

        // Step 2: Update position and velocity using weighted average of k-values
        body.vel += (k1_v + k2_v * 2. + k3_v * 2. + k4_v) / 6.;
        body.posit += (k1_posit + k2_posit * 2. + k3_posit * 2. + k4_posit) / 6.;
    }
}

// todo: DRY
fn integrate_rk4_ray(rays: &mut [GravRay], dt: f64) {
    // todo: Pending further exporation, no grav accel on rays.
    let a = Vec3::new_zero();

    for body in rays.iter_mut() {
        // Step 1: Calculate the k-values for position and velocity
        let k1_v = a * dt;
        let k1_posit = body.vel * dt;

        let k2_v = a * dt;
        let k2_posit = (body.vel + k1_v * 0.5) * dt;

        let k3_v = a * dt;
        let k3_posit = (body.vel + k2_v * 0.5) * dt;

        let k4_v = a * dt;
        let k4_posit = (body.vel + k3_v) * dt;

        // Step 2: Update position and velocity using weighted average of k-values
        body.vel += (k1_v + k2_v * 2. + k3_v * 2. + k4_v) / 6.;
        body.posit += (k1_posit + k2_posit * 2. + k3_posit * 2. + k4_posit) / 6.;
    }
}

// fn run(n_elecs: usize, n_timesteps: usize, dt: f64) -> Vec<SnapShot> {
//     let mut snapshots = Vec::new();
//
//     let charge_per_elec = Q_ELEC / n_elecs as f64;
//
//     let nuc = Nucleus {
//         mass: 2_000., // todo temp?
//         charge: -Q_ELEC,
//         posit: Vec3::new_zero(),
//     };
//
//     let mut elecs = make_initial_elecs(n_elecs);
//
//     for snap_i in 0..n_timesteps {
//         let len = elecs.len();
//         for elec_acted_on in 0..len {
//             let mut a = Vec3::new_zero();
//
//             // Force of other elecs on this elec.
//             for elec_actor in 0..len {
//                 // continue; // todo removing elec-elec charge for now.
//
//                 // Split the slice to get mutable and immutable elements
//                 // let (acted_on, actor) = elecs.split_at_mut(i);
//
//                 if elec_acted_on == elec_actor {
//                     continue;
//                 }
//
//                 a += accel_coulomb(
//                     elecs[elec_acted_on].posit,
//                     elecs[elec_actor].posit,
//                     charge_per_elec,
//                     charge_per_elec,
//                     M_ELEC,
//                     false,
//                 );
//             }
//
//             // Nuc force on this elec.
//             a += accel_coulomb(
//                 elecs[elec_acted_on].posit,
//                 nuc.posit,
//                 charge_per_elec,
//                 -Q_ELEC,
//                 M_ELEC,
//                 true,
//             );
//             elecs[elec_acted_on].a = a;
//         }
//
//         integrate_rk4(&mut elecs, dt);
//
//         snapshots.push(SnapShot {
//             time: snap_i as f64 * dt,
//             elec_posits: elecs.iter().map(|e| Vec3f32::new(e.posit.x as f32, e.posit.y as f32, e.posit.z as f32)).collect(),
//             nuc_posits: vec![Vec3f32::new(nuc.posit.x as f32, nuc.posit.y as f32, nuc.posit.z as f32)]
//         })
//     }
//
//     snapshots
// }

/// A rectangular prism for sampling properties (generally an "infinitessimal" volume.
struct SampleRect {
    pub start: Vec3,
    pub end: Vec3,
}

impl SampleRect {
    fn measure_properties(&self, rays: &[GravRay]) -> SampleProperties {
        let volume = {
            let dist_x = self.end.x - self.start.x;
            let dist_y = self.end.y - self.start.y;
            let dist_z = self.end.z - self.start.z;

            (dist_x.powi(2) + dist_y.powi(2) + dist_z.powi(2)).sqrt()
        };

        let mut num_rays = 0;
        for ray in rays {
            if ray.posit.x >= self.start.x
                && ray.posit.x <= self.end.x
                && ray.posit.y >= self.start.y
                && ray.posit.y <= self.end.y
                && ray.posit.z >= self.start.z
                && ray.posit.z <= self.end.z
            {
                num_rays += 1;
            }
        }

        // todo: To calculate div and curl, we need multiple sets of rays.

        SampleProperties {
            charge_density: num_rays as f64 / volume,
            div: 0.,
            curl: 0.,
        }
    }
}

/// Our model "particle" that travels outward from a source
struct GravRay {
    posit: Vec3,
    /// The magnitude corresponds to source mass. Direction is outward
    /// at any angle from the source.
    vel: Vec3,
    // todo: We can assume there is no acc on this, right?
}

#[derive(Debug)]
struct SampleProperties {
    charge_density: f64,
    /// Divergence
    div: f64,
    curl: f64,
}

/// Entry point for computation; rename A/R.
fn build(state: &mut State) {
    for t in 0..state.config.num_timesteps {
        // Create a new set of rays.
        for body in &state.bodies {
            for _ in 0..state.config.num_rays_per_iter {
                state.rays.push(body.create_ray())
            }
        }

        state.remove_far_rays();

        // Update ray propogation
        integrate_rk4_ray(&mut state.rays, state.config.dt_integration)

        // Update body motion.

        // Save the current state to a snapshot, for later playback.
        // Note: This can use a substantial amount of memory.
    }
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
        },
        Body {
            posit: Vec3::new(1., 0., 0.),
            vel: Vec3::new_zero(),
            accel: Vec3::new_zero(),
            mass: 1.,
        },
    ];

    build(&mut state);

    println!("Complete. Rendering...");
    render(state);
    // let snapshots = run(200, 50_000, 0.001);
    // let snapshots = run(1, 1000_000, 0.0001);
}
