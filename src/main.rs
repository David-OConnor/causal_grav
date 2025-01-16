#![allow(non_snake_case)]
#![allow(non_ascii_idents)]

use std::{path::PathBuf, time::Instant};

use lin_alg::f64::Vec3;
use rand::Rng;
use rayon::prelude::*;

use crate::{
    accel::MondFn,
    barnes_hut::{Cube, Tree},
    body_creation::{GalaxyDescrip, GalaxyModel},
    gaussian::{COEFF_C, MAX_SHELL_R},
    playback::{save, vec3_to_f32, GravShellSnapshot, SnapShot, DEFAULT_SNAPSHOT_FILE},
    render::render,
    units::{A0_MOND, C},
};

mod accel;
mod body_creation;
mod cdm;
mod fluid_dynamics;
mod fmm_gpt;
mod galaxy_data;
mod gaussian;
mod gem;
mod integrate;
mod playback;
mod properties;
mod render;
mod ui;
mod units;
mod util;
// mod fmm_gadget4;
// mod fmm_py;
mod barnes_hut;

// Shower thought, from looking at this from a first person view: View things from the body's perspective.
// Can you make of it something like that?

// todo: What if you have total electric charge different from a whole number? Misc elec "particles"
// todo zipping around in the either, to be capture.

// todo: Next, try having the electrons add up to more than 1 charge.
// todo: Try to quantify the elec density, so you can compare it to Schrodinger.

// todo: Soften gravity at short distances to prevent unstability?
// todo: F = G*m1*m2/(R^2 + eps^2)

// todo: Try adaptive time-steps, i.e. tighter steps whenever bodies get close to each other.

// todo: Use soft masses:  Bodies get a finite radius. Use plummer or hernquist potentials. (Similar to force softening?)

// todo: Thermal velocities? Init the system consistent with virial equilibrium.

// todo: Symplectit integrators: Leapfrog/WIdsom-Homan.

// Re MOND and our wavefront propogation. Something that could cause a reduced falloff of gravitation
// at distances:
//
// - A reduced wave propogation speed at distances (Causing the waves to be closer together; amplitude
// scales inverse proportionally to propogation speed C.) for a given generation rate
//
// Some sort of cumulative drag?? Can we estimate and test this?

const BOUNDING_BOX_PAD: f64 = 0.3;
const BB_GEN_RATIO: usize = 5;

pub struct Config {
    num_timesteps: usize,
    dt_integration_max: f64,
    /// Unit: MYR
    dt: f64, // Fixed.
    /// Lower values here lead to higher precision, and slower time evolution.
    dynamic_dt_scaler: f64,
    shell_creation_ratio: usize,
    // num_rays_per_iter: usize,
    gauss_c: f64,
    num_bodies_disk: usize, // todo: You may, in the future, not make this a constant.
    num_bodies_bulge: usize, // todo: You may, in the future, not make this a constant.
    /// When placing bodies.
    num_rings_disk: usize, // todo: You may, in the future, not make this a constant.
    num_rings_bulge: usize, // todo: You may, in the future, not make this a constant.
    softening_factor_sq: f64,
    snapshot_ratio: usize,
    /// Smaller values use less grouping; these are slower, but more accurate.
    barnes_hut_θ: f64,
}

impl Default for Config {
    fn default() -> Self {
        let dt = 2.0e-3;
        let shell_creation_ratio = 12;

        // Important: Shell spacing is only accurate if using non-dynamic DT.

        // In distance: t * d/t = d.
        let shell_spacing = dt * shell_creation_ratio as f64 * C;

        let num_bodies_disk = 300;
        let num_bodies_bulge = 100;

        // todo: Delegate this.
        let num_rings_disk = num_bodies_disk / 10;
        let num_rings_bulge = num_bodies_bulge / 10;

        Self {
            num_timesteps: 5_000,
            shell_creation_ratio,
            dt,
            dt_integration_max: 0.01, // Not used
            // dynamic_dt_scaler: 0.01,
            dynamic_dt_scaler: 0.1, // not used.
            // num_rays_per_iter: 200,
            gauss_c: shell_spacing * COEFF_C,
            num_bodies_disk,
            num_rings_disk,
            num_bodies_bulge,
            num_rings_bulge,
            softening_factor_sq: 0.01,
            snapshot_ratio: 4,
            barnes_hut_θ: 0.4,
        }
    }
}

#[derive(Copy, Clone, PartialEq)]
pub enum ForceModel {
    Newton,
    Mond(MondFn),
    GaussShells,
}

impl Default for ForceModel {
    fn default() -> Self {
        ForceModel::Newton
    }
}

pub struct StateUi {
    snapshot_selected: usize,
    force_model: ForceModel,
    building: bool,
    /// We include text input fields for user-typeable floats. Not required for int.
    dt_input: String,
    θ_input: String,
    // num_timesteps_input: String,
    add_halo: bool, // todo: A/R
    galaxy_model: GalaxyModel,
    /// For display in the UI. cached.
    galaxy_descrip: GalaxyDescrip,
}

impl Default for StateUi {
    fn default() -> Self {
        let galaxy_model = Default::default();

        Self {
            snapshot_selected: Default::default(),
            force_model: Default::default(),
            building: Default::default(),
            dt_input: Default::default(),
            θ_input: Default::default(),
            add_halo: Default::default(),
            galaxy_model,
            galaxy_descrip: galaxy_model.descrip(),
        }
    }
}

#[derive(Default)]
struct State {
    config: Config,
    ui: StateUi,
    bodies: Vec<Body>,
    // rays: Vec<GravRay>,
    shells: Vec<GravShell>,
    snapshots: Vec<SnapShot>,
    /// For rendering; separate from snapshots since it's invariant.
    body_masses: Vec<f32>,
    time_elapsed: f64,
}

impl State {
    fn remove_far_shells(&mut self) {
        self.shells.retain(|shell| shell.radius <= MAX_SHELL_R);
    }

    fn take_snapshot(&mut self, time: f64, dt: f64) {
        self.snapshots.push(SnapShot {
            time: time as f32,
            body_posits: self.bodies.iter().map(|b| vec3_to_f32(b.posit)).collect(),
            // V_at_bodies: self
            //     .bodies
            //     .iter()
            //     .map(|b| vec_to_f32(b.V_acting_on))
            //     .collect(),
            body_accs: self.bodies.iter().map(|b| vec3_to_f32(b.accel)).collect(),
            // rays: self
            //     .rays
            //     .iter()
            //     .map(|r| {
            //         (
            //             Vec3f32::new(r.posit.x as f32, r.posit.y as f32, r.posit.z as f32),
            //             r.emitter_id,
            //         )
            //     })
            //     .collect(),
            shells: self
                .shells
                .iter()
                .map(|s| GravShellSnapshot::new(s))
                .collect(),
            dt: dt as f32,
        })
    }
}

#[derive(Clone, Debug)]
struct Body {
    posit: Vec3,
    vel: Vec3,
    accel: Vec3,
    mass: f64,
    // /// We use this for debugging, testing etc. It goes in the snapshots.
    // V_source: Vec3,
}

impl Body {
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

#[derive(Debug, Clone)]
struct GravShell {
    emitter_id: usize,
    center: Vec3,
    radius: f64,
    src_mass: f64,
}

/// Entry point for computation; rename A/R.
fn build(state: &mut State, force_model: ForceModel) {
    state.ui.building = true;
    println!("Building...");

    // Allow gravity shells to propogate to reach a steady state, ideally.
    // todo: make this dynamic
    let mut farthest_r = 0.;
    for body in &state.bodies {
        let r = body.posit.magnitude();
        if r > farthest_r {
            farthest_r = r;
        }
    }
    // 2x: For the case of opposite sides of circle. More: A pad. May not be required.
    farthest_r *= 2.5;
    let integrate_start_t = farthest_r / C;

    println!(
        "T start integration: {:?} T: {:?}. Farthest r: {:.1}",
        integrate_start_t, state.time_elapsed, farthest_r
    );

    let mut start_time = Instant::now();
    let mut start_time_tree = Instant::now();
    let mut start_time_integ = Instant::now();
    // let mut start_time_general = Instant::now();

    let mut bb = Cube::from_bodies(&state.bodies, BOUNDING_BOX_PAD, true).unwrap();

    const BENCH_RATIO: usize = 1_000;

    for t in 0..state.config.num_timesteps {
        // if t % BENCH_RATIO == 0 {
        //     start_time_general = Instant::now();
        // }

        if force_model == ForceModel::GaussShells {
            if t % state.config.shell_creation_ratio == 0 {
                for (id, body) in state.bodies.iter().enumerate() {
                    state.shells.push(body.create_shell(id));
                }
            }
            for shell in &mut state.shells {
                shell.radius += C * state.config.dt;
            }
        }

        if force_model == ForceModel::GaussShells && t % state.config.shell_creation_ratio == 0 {
            state.remove_far_shells(); // Note grouped above due to a borrow problem.
        }

        // todo: Put back.
        // if t % BB_GEN_RATIO == 0 {
        bb = Cube::from_bodies(&state.bodies, BOUNDING_BOX_PAD, true).unwrap();
        // }

        if bb.width.is_nan() {
            eprintln!("Error: NaN");
            return;
        }

        // Calculate dt for this step, based on the closest/fastest rel velocity.
        // This affects motion integration only; not shell creation.
        // let dt = util::calc_dt_dynamic(state, &bodies_other);
        // todo: Static DT for now, or shells won't work.
        let dt = state.config.dt;

        if t % BENCH_RATIO == 0 {
            start_time_tree = Instant::now();
        }

        let mut tree = None;
        if force_model != ForceModel::GaussShells {
            tree = Some(Tree::new(&state.bodies, &bb));
        }

        if t % BENCH_RATIO == 0 {
            println!("Tree time: {}μs", start_time_tree.elapsed().as_micros());
            // start_time = Instant::now();
        }

        // Benchmarking, 100k bodies, 2025-01-16, theta = 0.4, BH algo.
        // Without rayon: Tree time: 51ms. N body time: 1,371ms
        // With rayon: Tree time: 51ms N body time: 144ms (Solid speedup)

        // if t % BENCH_RATIO == 0 {
        //     println!("N body time: {}μs", start_time.elapsed().as_micros());
        // }

        state.time_elapsed += dt;

        if t % BENCH_RATIO == 0 {
            start_time_integ = Instant::now();
        }
        if force_model != ForceModel::GaussShells || state.time_elapsed > integrate_start_t {
            // for body in &mut state.bodies {
            //     body.vel += body.accel * dt;
            //     body.posit += body.vel * dt;
            // }

            // Update body motion.
            integrate::integrate_rk4(
                &mut state.bodies,
                &state.shells,
                dt,
                state.config.gauss_c,
                &tree.as_ref().unwrap(),
                state.config.barnes_hut_θ,
                force_model,
                state.config.softening_factor_sq,
            );
        }

        if t % BENCH_RATIO == 0 {
            println!("Integ time: {}μs", start_time_integ.elapsed().as_micros());
        }

        // Save the current state to a snapshot, for later playback.
        if t % state.config.snapshot_ratio == 0 {
            state.take_snapshot(state.time_elapsed, dt);
        }

        // if t % BENCH_RATIO == 0 {
        //     println!("General time: {}μs", start_time_general.elapsed().as_micros());
        // }
    }

    state.ui.building = false;
    println!("Build complete.");
}

fn main() {
    let mut state = State::default();

    state.bodies = state.ui.galaxy_descrip.make_bodies(
        state.config.num_bodies_disk,
        state.config.num_rings_disk,
        state.config.num_bodies_bulge,
        state.config.num_rings_bulge,
    );
    state.body_masses = state.bodies.iter().map(|b| b.mass as f32).collect();

    state.ui.dt_input = state.config.dt.to_string();
    state.ui.θ_input = state.config.barnes_hut_θ.to_string();

    // todo: Don't auto-build, but we a graphics engine have a prob when we don't.
    state.take_snapshot(0., 0.); // Initial snapshot; t=0.

    let rotation_curve = properties::rotation_curve(&state.bodies, Vec3::new_zero(), C);
    let mass_density = properties::mass_density(&state.bodies, Vec3::new_zero());

    properties::plot_rotation_curve(&rotation_curve, &state.ui.galaxy_model.to_str());
    properties::plot_mass_density(&mass_density, &state.ui.galaxy_model.to_str());

    if let Err(e) = save(&PathBuf::from(DEFAULT_SNAPSHOT_FILE), &state.snapshots) {
        eprintln!("Error saving snapshots: {e}");
    }

    println!("a_0 (MOND): {:?}", A0_MOND);

    render(state);
}
