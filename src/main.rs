#![allow(non_snake_case)]
#![allow(non_ascii_idents)]

use lin_alg::f64::Vec3;
use rand::Rng;

use crate::{
    gaussian::{COEFF_C, MAX_SHELL_R},
    playback::{vec_to_f32, SnapShot},
    render::render,
};
// Add this to use the random number generator

mod accel;
mod body_creation;
mod gaussian;
mod integrate;
mod playback;
mod render;
mod ui;
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
const C: f64 = 40.;

pub struct Config {
    num_timesteps: usize,
    dt_integration: f64,
    shell_creation_ratio: usize,
    // num_rays_per_iter: usize,
    gauss_c: f64,
}

impl Default for Config {
    fn default() -> Self {
        let dt_integration = 0.001;
        let shell_creation_ratio = 6;

        // In distance: t * d/t = d.
        let shell_spacing = dt_integration * shell_creation_ratio as f64 * C;

        Self {
            num_timesteps: 15_000,
            shell_creation_ratio,
            dt_integration,
            // num_rays_per_iter: 200,
            gauss_c: shell_spacing * COEFF_C,
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
    // rays: Vec<GravRay>,
    shells: Vec<GravShell>,
    snapshots: Vec<SnapShot>,
    /// For rendering; separate from snapshots since it's invariant.
    body_masses: Vec<f32>,
    /// Defaults to `Config::dt_integration`, but becomes more precise when
    /// bodies are close.
    dt_dynamic: f64,
}

impl State {
    fn remove_far_shells(&mut self) {
        self.shells.retain(|shell| shell.radius <= MAX_SHELL_R);
    }

    fn take_snapshot(&mut self, time: usize) {
        self.snapshots.push(SnapShot {
            time,
            body_posits: self.bodies.iter().map(|b| vec_to_f32(b.posit)).collect(),
            // V_at_bodies: self
            //     .bodies
            //     .iter()
            //     .map(|b| vec_to_f32(b.V_acting_on))
            //     .collect(),
            body_accs: self.bodies.iter().map(|b| vec_to_f32(b.accel)).collect(),
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
            // shells: self.shells.clone(),
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
fn build(state: &mut State) {
    for t in 0..state.config.num_timesteps {
        // Create a new set of rays.
        if t % state.config.shell_creation_ratio == 0 {
            for (id, body) in state.bodies.iter().enumerate() {
                // for _ in 0..state.config.num_rays_per_iter {
                // state.rays.push(body.create_ray(id));
                // }
                state.shells.push(body.create_shell(id));
            }
        }

        // for (i, body) in state.bodies.iter_mut().enumerate() {
        // body.accel = accel::accel(
        //     body,
        //     &state.rays,
        //     &state.shells,
        //     i,
        //     state.config.dt_integration,
        // );

        // body.accel = accel::acc_shells(
        //     body,
        //     &state.rays,
        //     &state.shells,
        //     i,
        //     state.config.dt_integration,
        //     state.config.gauss_c,
        // );

        // body.accel = accel::calc_acc_shell(&state.shells, body.posit, i, state.config.gauss_c);
        // }

        // Update ray propogation
        // integrate::integrate_rk4_ray(&mut state.rays, state.config.dt_integration);

        // state.remove_far_rays();
        state.remove_far_shells();

        for shell in &mut state.shells {
            shell.radius += C * state.config.dt_integration;
        }

        let acc_inst = true;
        // let acc_inst = false;

        // todo: C+P from integrate, so we can test acc vals.
        let bodies_other = state.bodies.clone(); // todo: I don't like this. Avoids mut error.

        let acc = |id, posit| {
            if acc_inst {
                accel::calc_acc_inst(posit, &bodies_other, id)
            } else {
                accel::calc_acc_shell(&state.shells, posit, id, state.config.gauss_c)
            }
        };

        for (id, body) in &mut state.bodies.iter_mut().enumerate() {
            body.accel = acc(id, body.posit);

            // todo: Ideally, use distance.
            if body.accel.magnitude() > 69. {
                let dt = state.config.dt_integration / 2.;
                // todo: Split time stamp.
            }
        }

        state.take_snapshot(t / SNAPSHOT_RATIO); // Initial snapshot; t=0.

        // Allow waves to propogate to reach a steady state, ideally.
        if acc_inst || t > 1_000 {
            // Update body motion.
            integrate::integrate_rk4(
                &mut state.bodies,
                &state.shells,
                state.config.dt_integration,
                state.config.gauss_c,
                acc_inst,
            );
        }

        // Save the current state to a snapshot, for later playback.
        // Note: This can use a substantial amount of memory.

        if t % SNAPSHOT_RATIO == 0 {
            state.take_snapshot(t / SNAPSHOT_RATIO);
        }
        // println!("Shell ct: {:?}", state.shells.len());
    }
}

fn main() {
    println!("Building snapshots...");
    let mut state = State::default();
    state.dt_dynamic = state.config.dt_integration; // todo: Integrate this into State::default();

    state.bodies = body_creation::make_galaxy_coarse(4, 6);
    // state.bodies = body_creation::make_galaxy_coarse(10, 8);
    state.body_masses = state.bodies.iter().map(|b| b.mass as f32).collect();

    build(&mut state);

    println!("Complete. Rendering.");

    render(state);
}
