use rand::Rng; // Add this to use the random number generator
use std::f64::consts::TAU;

use lin_alg::f64::Vec3;

mod ui;
mod playback;
mod render;
// todo: What if you have total electric charge different from a whole number? Misc elec "particles"
// todo zipping around in the either, to be capture.


// todo: Next, try having the electrons add up to more than 1 charge.
// todo: Try to quantify the elec density, so you can compare it to Schrodinger.

const M_ELEC: f64 = 1.;
const Q_ELEC: f64 = -1.;

// If two particles are closer to each other than this, don't count accel, or cap it.
const MIN_DIST: f64 = 0.1;


// Don't calculate force if particles are farther than this from each other. Computation saver.
// const MAX_DIST: f64 = 0.1;


pub struct Config {
    num_timesteps: usize,
    dt_integration: f64,
    dt_pulse: f64,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            num_timesteps: 1_000,
            dt_integration: 0.01,
            dt_pulse: 0.01
        }
    }
}

#[derive(Default)]
pub struct StateUi {
    snapshot_selected: usize
}


#[derive(Default)]
struct State {
    config: Config,
    ui: StateUi,
    bodies: Vec<Body>,
    snapshots: Vec<SnapShot>
}

struct Body {
    posit: Vec3,
    vel: Vec3,
}

// (radius, charge_in_radius)
fn density(elec_posits: &[Vec3f32], q_per_elec: f64, center_ref: Vec3f32) -> Vec<((f32, f32), f32)> {
    // Must be ascending radii for the below code to work.
    let per_elec = q_per_elec as f32;

    // Equidistant ish, for now.
    let mut result = vec![
        ((0.0, 0.2), 0.),
        ((0.2, 0.4),0.),
        ((0.4, 0.6),0.),
        ((0.8, 1.), 0.),
        ((1., 1.2),0.),
        ((1.2, 1.4),0.),
        ((1.4, 1.6), 0.),
        ((1.6, 1.8), 0.),
        ((1.8, 2.0), 0.),
        ((2.0, 9999.), 0.),
    ];

    for (r_bound, q) in &mut result {
        for posit in elec_posits {
            let mag = (center_ref - *posit).magnitude();
            if mag < r_bound.1 && mag > r_bound.0 {
                *q += per_elec;
            }
        }
    }

    result
}

/// Calculate the Coulomb or gravitational acceleration on a particle, from a single other particle.
fn accel(
    posit_acted_on: Vec3,
    posit_actor: Vec3,
    q_acted_on: f64,
    q_actor: f64,
    mass_acted_on: f64,
    nuc_actor: bool
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

// fn integrate_rk4(bodies: &mut [Body], dt: f64) {
//     for elec in elecs.iter_mut() {
//         // Step 1: Calculate the k-values for position and velocity
//         let k1_v = elec.a * dt;
//         let k1_posit = elec.v * dt;
//
//         let k2_v = (elec.a) * dt;
//         let k2_posit = (elec.v + k1_v * 0.5) * dt;
//
//         let k3_v = (elec.a) * dt;
//         let k3_posit = (elec.v + k2_v * 0.5) * dt;
//
//         let k4_v = (elec.a) * dt;
//         let k4_posit = (elec.v + k3_v) * dt;
//
//         // Step 2: Update position and velocity using weighted average of k-values
//         elec.v += (k1_v + k2_v * 2. + k3_v * 2. + k4_v) / 6.;
//         elec.posit += (k1_posit + k2_posit * 2. + k3_posit * 2. + k4_posit) / 6.;
//     }
// }


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

/// Entry point for computation; rename A/R.
fn build(state: &mut State) {
    for t in 0..state.config.num_timesteps {

    }
}


fn main() {
    println!("Building snapshots...");
    let mut state = State::default();

    state.bodies = vec![
        Body {posit: Vec3::new_zero(), vel: Vec3::new_zero() },
        Body {posit: Vec3::new(1., 0., 0.), vel: Vec3::new_zero() },
    ];

    build(&mut state);

    println!("Complete. Rendering...");
    ui::render(state);
    // let snapshots = run(200, 50_000, 0.001);
    // let snapshots = run(1, 1000_000, 0.0001);
}
