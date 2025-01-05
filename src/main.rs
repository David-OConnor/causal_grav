#![allow(non_snake_case)]
#![allow(non_ascii_idents)]

use std::path::PathBuf;

use lin_alg::f64::Vec3;
use rand::Rng;

use crate::{
    body_creation::GalaxyModel,
    gaussian::{COEFF_C, MAX_SHELL_R},
    playback::{save, vec3_to_f32, GravShellSnapshot, SnapShot, DEFAULT_SNAPSHOT_FILE},
    render::render,
};

mod accel;
mod body_creation;
mod fluid_dynamics;
mod gaussian;
mod integrate;
mod playback;
mod properties;
mod render;
mod ui;
mod util;
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

const SNAPSHOT_RATIO: usize = 4;

// Note: Setting this too high is problematic.
// todo: Maybe a different time unit?
// const C: f64 = 9.72e-12; // Rough; kpc/s^2.
const C: f64 = 40.;

pub struct Config {
    num_timesteps: usize,
    dt_integration_max: f64,
    /// Lower values here lead to higher precision, and slower time evolution.
    dynamic_dt_scaler: f64,
    shell_creation_ratio: usize,
    // num_rays_per_iter: usize,
    gauss_c: f64,
}

impl Default for Config {
    fn default() -> Self {
        let dt_integration_max = 0.01;
        let shell_creation_ratio = 12;

        // In distance: t * d/t = d.
        let shell_spacing = dt_integration_max * shell_creation_ratio as f64 * C;

        Self {
            // num_timesteps: 1_000_000,
            num_timesteps: 10_000,
            // num_timesteps: 1000,
            shell_creation_ratio,
            dt_integration_max,
            // dynamic_dt_scaler: 0.01,
            dynamic_dt_scaler: 0.1,
            // num_rays_per_iter: 200,
            gauss_c: shell_spacing * COEFF_C,
        }
    }
}

#[derive(Copy, Clone, PartialEq)]
pub enum ForceModel {
    Newton,
    Mond(f64), // inner is a placeholder for a coefficient
    GaussRings,
}

impl Default for ForceModel {
    fn default() -> Self {
        ForceModel::Newton
    }
}

#[derive(Default)]
pub struct StateUi {
    snapshot_selected: usize,
    force_model: ForceModel,
    building: bool,
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
    // /// Defaults to `Config::dt_integration`, but becomes more precise when
    // /// bodies are close. This is a global DT, vice local only for those bodies.
    // dt_dynamic: f64,
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
    // V_acting_on: Vec3,
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
    state.take_snapshot(0., 0.); // Initial snapshot; t=0.

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

    for t in 0..state.config.num_timesteps {
        // Create a new set of rays.
        if force_model == ForceModel::GaussRings && t % state.config.shell_creation_ratio == 0 {
            for (id, body) in state.bodies.iter().enumerate() {
                // for _ in 0..state.config.num_rays_per_iter {
                // state.rays.push(body.create_ray(id));
                // }
                state.shells.push(body.create_shell(id));
            }
        }

        // Update ray propogation
        // integrate::integrate_rk4_ray(&mut state.rays, state.config.dt_integration);

        // state.remove_far_rays();
        state.remove_far_shells();

        for shell in &mut state.shells {
            shell.radius += C * state.config.dt_integration_max;
        }

        // todo: C+P from integrate, so we can test acc vals.
        let bodies_other = state.bodies.clone(); // todo: I don't like this. Avoids mut error.

        let acc = |id, posit| match force_model {
            ForceModel::Newton => accel::acc_newton(posit, &bodies_other, id),
            ForceModel::GaussRings => {
                accel::calc_acc_shell(&state.shells, posit, id, state.config.gauss_c)
            }
            ForceModel::Mond(_) => unimplemented!(),
        };

        // Calculate dt for this step, based on the closest/fastest rel velocity.
        // This affects motion integration only; not shell creation.
        // todo: Separate fn
        let mut dt_min = state.config.dt_integration_max;
        // todo: Consider cacheing the distances, so this second iteration can be reused.
        for (id_acted_on, body) in &mut state.bodies.iter_mut().enumerate() {
            for (i, body_src) in bodies_other.iter().enumerate() {
                if i == id_acted_on {
                    continue; // self-interaction.
                }

                let dist = (body_src.posit - body.posit).magnitude();
                let rel_velocity = (body_src.vel - body.vel).magnitude();
                let dt = state.config.dynamic_dt_scaler * dist / rel_velocity;
                if dt < dt_min {
                    dt_min = dt;
                }
            }
        }

        for (id, body) in &mut state.bodies.iter_mut().enumerate() {
            body.accel = acc(id, body.posit);
        }

        state.time_elapsed += dt_min;

        // todo: Is
        if force_model != ForceModel::GaussRings || state.time_elapsed > integrate_start_t {
            // Update body motion.
            integrate::integrate_rk4(
                &mut state.bodies,
                &state.shells,
                dt_min,
                state.config.gauss_c,
                force_model,
            );
        }

        // Save the current state to a snapshot, for later playback.
        // Note: This can use a substantial amount of memory.

        if t % SNAPSHOT_RATIO == 0 {
            state.take_snapshot(state.time_elapsed, dt_min);
        }
        // println!("Shell ct: {:?}", state.shells.len());
    }

    state.ui.building = false;
}

fn main() {
    println!("Building snapshots...");
    let mut state = State::default();
    // state.dt_dynamic = state.config.dt_integration; // todo: Integrate this into State::default();

    let model = GalaxyModel::Ngc1560;
    state.bodies = model.make_bodies();

    state.body_masses = state.bodies.iter().map(|b| b.mass as f32).collect();

    state.ui.force_model = ForceModel::Newton;

    let fm = state.ui.force_model;
    // todo: Don't auto-build, but we a graphics engine have a prob when we don't.
    build(&mut state, fm);

    let rotation_curve = properties::rotation_curve(&state.bodies, Vec3::new_zero(), 80., C);
    let mass_density = properties::mass_density(&state.bodies, Vec3::new_zero(), 80.);

    properties::plot_rotation_curve(&rotation_curve, "NGC 3198");
    properties::plot_mass_density(&mass_density, "NGC 3198");

    println!("Complete. Rendering.");

    if let Err(e) = save(&PathBuf::from(DEFAULT_SNAPSHOT_FILE), &state.snapshots) {
        eprintln!("Error saving snapshots: {e}");
    }

    render(state);
}
