#![allow(non_snake_case)]
#![allow(non_ascii_idents)]

use std::{
    io,
    path::{Path, PathBuf},
    str::FromStr,
    sync::Arc,
    time::Instant,
};

use barnes_hut::{BhConfig, BodyModel, Cube, Node, Tree};
use bincode::{Decode, Encode};
#[cfg(feature = "cuda")]
use cudarc::{driver::CudaDevice, nvrtc::Ptx};
use galaxy_data::GalaxyModel;
use grav_shell::{GravShell, MAX_SHELL_R};
use lin_alg::f64::Vec3;
use rand::Rng;
use rayon::prelude::*;

use crate::{
    accel::{acc_newton_inner_with_mond, MondFn},
    body_creation::GalaxyDescrip,
    gaussian::GaussianShell,
    grav_shell::COEFF_C,
    integrate::integrate_rk4,
    playback::{vec3_to_f32, GravShellSnapshot, SnapShot},
    render::render,
    units::{A0_MOND, C},
};

mod accel;
mod body_creation;
mod cdm;
mod fluid_dynamics;
// mod fmm_gpt;
mod galaxy_data;
mod gaussian;
mod gem;
#[cfg(feature = "cuda")]
mod gpu;
mod grav_shell;
mod image_parsing;
mod integrate;
mod playback;
mod properties;
mod render;
mod ui;
mod units;
mod util;
// todo: Try a Galaxy filament simulation; large scale CDM theory. Can we get filaments without CDM?
// todo: - try an earth-perspective visualization and analysis. From the perspective of earth, validate these
// todo galaxies vs the images we get.

// Shower thought, from looking at this from a first person view: View things from the body's perspective.
// Can you make of it something like that?

// todo: What if you have total electric charge different from a whole number? Misc elec "particles"
// todo zipping around in the either, to be capture.

// todo: Next, try having the electrons add up to more than 1 charge.
// todo: Try to quantify the elec density, so you can compare it to Schrodinger.

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

// const BOUNDING_BOX_PAD: f64 = 0.3;
const BOUNDING_BOX_PAD: f64 = 0.;
const BB_GEN_RATIO: usize = 1;

const SAVE_FILE: &str = "config.grav";
const DEFAULT_SNAPSHOT_FILE: &str = "snapshot.grav";

const DISK_RING_PORTION: usize = 10;
const BULGE_RING_PORTION: usize = 5;

#[derive(Debug, Clone)]
pub enum ComputationDevice {
    Cpu,
    #[cfg(feature = "cuda")]
    Gpu(Arc<CudaDevice>),
}

// todo: Custom Bincode config that only contains the fields you customize directly.
#[derive(Encode, Decode)]
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
    // num_rings_disk: usize, // todo: You may, in the future, not make this a constant.
    // num_rings_bulge: usize, // todo: You may, in the future, not make this a constant.
    softening_factor_sq: f64,
    snapshot_ratio: usize,
    bh_config: BhConfig,
    /// Experimental tool; scale all published V magnitudes by this.
    v_scaler: f64,
    /// Use instantaneous Newtonian forces instead of tree code.
    skip_tree: bool,
}

impl Default for Config {
    fn default() -> Self {
        let dt = 2.0e-3;

        // Important: Shell spacing is only accurate if using non-dynamic DT.

        let num_bodies_disk = 600;
        let num_bodies_bulge = 100;

        let mut result = Self {
            num_timesteps: 5_000,
            shell_creation_ratio: 1,
            // shell_creation_ratio: 12,
            dt,
            dt_integration_max: 0.01, // Not used
            // dynamic_dt_scaler: 0.01,
            dynamic_dt_scaler: 0.1, // not used.
            // num_rays_per_iter: 200,
            gauss_c,
            num_bodies_disk,
            num_bodies_bulge,
            softening_factor_sq: 0.0001,
            // softening_factor_sq: 0.05,
            snapshot_ratio: 8,
            bh_config: BhConfig {
                // θ: 0.4,
                ..Default::default()
            },
            v_scaler: 1.0,
            skip_tree: false,
        };

        result.refresh_shell_ratio();
        result
    }
}

impl Config {
    /// Load the config.w
    pub fn load(&mut self, path: &Path) -> io::Result<Self> {
        util::load(path)
    }

    /// Run this whne dt, or shell creation ratio changes, to update the Gaussian width.
    pub fn refresh_shell_ratio(&mut self) {
        // In distance: t * d/t = d.
        let shell_spacing = self.dt * self.shell_creation_ratio as f64 * C;
        self.gauss_c = shell_spacing * COEFF_C;
    }
}

#[derive(Copy, Clone, PartialEq, Default)]
pub enum ForceModel {
    #[default]
    Newton,
    Mond(MondFn),
    GaussShells,
}

pub struct StateUi {
    snapshot_selected: usize,
    force_model: ForceModel,
    building: bool,
    /// We include text input fields for user-typeable floats. Not required for int.
    dt_input: String,
    θ_input: String,
    v_scaler_input: String,
    // num_timesteps_input: String,
    add_halo: bool, // todo: A/R
    galaxy_model: GalaxyModel,
    /// For display in the UI. cached.
    galaxy_descrip: GalaxyDescrip,
    draw_tree: bool,
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
            v_scaler_input: Default::default(),
            add_halo: Default::default(),
            galaxy_model,
            galaxy_descrip: galaxy_model.descrip(),
            draw_tree: false,
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
    fn refresh_bodies(&mut self) {
        self.bodies = self.ui.galaxy_descrip.make_bodies(
            self.config.num_bodies_disk,
            self.config.num_bodies_bulge,
            self.config.v_scaler,
        );
        self.body_masses = self.bodies.iter().map(|b| b.mass as f32).collect();

        self.time_elapsed = 0.;
        self.snapshots = Vec::new();
        self.take_snapshot(0., Vec::new()); // Initial snapshot; t=0.
        self.ui.snapshot_selected = 0;

        self.shells = Vec::new();

        let rotation_curve = properties::rotation_curve(&self.bodies, Vec3::new_zero(), C);
        let mass_density = properties::mass_density(&self.bodies, Vec3::new_zero());

        properties::plot_rotation_curve(&rotation_curve, &self.ui.galaxy_model.to_str());
        properties::plot_mass_density(&mass_density, &self.ui.galaxy_model.to_str());
    }

    fn remove_far_shells(&mut self) {
        self.shells.retain(|shell| shell.radius <= MAX_SHELL_R);
    }

    fn take_snapshot(&mut self, dt: f64, tree_nodes: Vec<Cube>) {
        self.snapshots.push(SnapShot {
            time: self.time_elapsed as f32,
            body_posits: self.bodies.iter().map(|b| vec3_to_f32(b.posit)).collect(),
            body_accs: self.bodies.iter().map(|b| vec3_to_f32(b.accel)).collect(),
            shells: self.shells.iter().map(GravShellSnapshot::new).collect(),
            dt: dt as f32,
            tree_cubes: tree_nodes,
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
    pub fn create_shell(&self, source_id: usize) -> GravShell {
        GravShell {
            source_id,
            center: self.posit,
            radius: 0.,
            src_mass: self.mass,
            body_vel: self.vel,   // todo: A/R
            body_acc: self.accel, // todo: A/R
        }
    }
}

impl BodyModel for Body {
    fn posit(&self) -> Vec3 {
        self.posit
    }

    fn mass(&self) -> f64 {
        self.mass
    }
}

/// Entry point for computation; rename A/R.
fn build(state: &mut State, force_model: ForceModel) {
    println!("Building...");
    state.ui.building = true;

    // We must refresh bodies prior to building, to reset their positions after the previous update.
    state.refresh_bodies();

    // Allow gravity shells to propogate to reach a steady state, ideally.
    // todo: make this dynamic
    let mut farthest_r = 0.;
    for body in &state.bodies {
        let r = body.posit.magnitude();
        if r > farthest_r {
            farthest_r = r;
        }
    }
    // 2x: For the case of opposite sides of circle.
    farthest_r *= 2.;
    let integrate_start_t = farthest_r / C;

    println!(
        "T start integration: {:?} T: {:?}. Farthest r: {:.1}",
        integrate_start_t, state.time_elapsed, farthest_r
    );

    let mut start_time_tree = Instant::now();
    let mut tree_time = 0;
    let mut start_time_integ = Instant::now();

    // todo: A/R.
    // let mut bb = Cube::from_bodies(&state.bodies, BOUNDING_BOX_PAD, true).unwrap();
    let mut bb = Cube::from_bodies(&state.bodies, BOUNDING_BOX_PAD, false).unwrap();

    const BENCH_RATIO: usize = 1_000;

    // Store a partially pre-allocated Vec of the nodes to be used by the tree, as an optimization.
    // This is faster than creating it every time.
    // todo: Make a guess at capacity based on bodies.
    // let mut nodes = Vec::with_capacity(69);

    for t in 0..state.config.num_timesteps {
        if force_model == ForceModel::GaussShells && t % state.config.shell_creation_ratio == 0 {
            state.remove_far_shells(); // Note grouped above due to a borrow problem.
        }

        if force_model == ForceModel::GaussShells {
            if t % state.config.shell_creation_ratio == 0 {
                for (id, body) in state.bodies.iter().enumerate() {
                    state.shells.push(body.create_shell(id));
                }
            }
            for shell in &mut state.shells {
                shell.iter_t(state.config.dt);
            }
        }

        let cfg = &state.config; // Code cleaner.

        if t % BB_GEN_RATIO == 0 && !cfg.skip_tree {
            bb = Cube::from_bodies(&state.bodies, BOUNDING_BOX_PAD, true).unwrap();
        }

        if bb.width.is_nan() {
            eprintln!("Error: NaN");
            return;
        }

        // Calculate dt for this step, based on the closest/fastest rel velocity.
        // This affects motion integration only; not shell creation.
        // let dt = util::calc_dt_dynamic(state, &bodies_other);
        // todo: Static DT for now, or shells won't work.
        let dt = cfg.dt;

        if t % BENCH_RATIO == 0 {
            start_time_tree = Instant::now();
        }

        let mut tree = None;
        if force_model != ForceModel::GaussShells && !cfg.skip_tree {
            tree = Some(Tree::new(&state.bodies, &bb, &cfg.bh_config));
        }

        if t % BENCH_RATIO == 0 && force_model != ForceModel::GaussShells && !cfg.skip_tree {
            tree_time = start_time_tree.elapsed().as_micros();
        }

        // Benchmarking, 100k bodies, 2025-01-16, theta = 0.4, BH algo.
        // Without rayon: Tree time: 51ms. N body time: 1,371ms
        // With rayon: Tree time: 51ms N body time: 144ms (Solid speedup)

        state.time_elapsed += dt;

        if t % BENCH_RATIO == 0 {
            start_time_integ = Instant::now();
        }

        let bodies_other = if cfg.skip_tree {
            Some(state.bodies.clone())
        } else {
            None
        };

        // This acceleration function acts on a target id and position.
        let acc = |id_target, posit_target| match force_model {
            ForceModel::Newton => {
                if cfg.skip_tree {
                    accel::acc_newton(
                        posit_target,
                        id_target,
                        &bodies_other.as_ref().unwrap(),
                        None,
                        cfg.softening_factor_sq,
                    )
                } else {
                    let acc_fn = |acc_dir, mass, dist| {
                        acc_newton_inner_with_mond(
                            acc_dir,
                            mass,
                            dist,
                            None,
                            cfg.softening_factor_sq,
                        )
                    };

                    barnes_hut::acc_bh(
                        posit_target,
                        id_target,
                        tree.as_ref().unwrap(),
                        &cfg.bh_config,
                        &acc_fn,
                    )
                }
            }
            ForceModel::GaussShells => accel::calc_acc_shell(
                &state.shells,
                posit_target,
                id_target,
                cfg.gauss_c,
                cfg.softening_factor_sq,
            ),
            ForceModel::Mond(mond_fn) => {
                if cfg.skip_tree {
                    accel::acc_newton(
                        posit_target,
                        id_target,
                        &bodies_other.as_ref().unwrap(),
                        Some(mond_fn),
                        cfg.softening_factor_sq,
                    )
                } else {
                    let acc_fn = |acc_dir, mass, dist| {
                        acc_newton_inner_with_mond(
                            acc_dir,
                            mass,
                            dist,
                            Some(mond_fn),
                            cfg.softening_factor_sq,
                        )
                    };

                    barnes_hut::acc_bh(
                        posit_target,
                        id_target,
                        tree.as_ref().unwrap(),
                        &cfg.bh_config,
                        &acc_fn,
                    )
                }
            }
        };

        if force_model != ForceModel::GaussShells || state.time_elapsed > integrate_start_t {
            // todo: COme back to skiping the first body. Setting the central body as immovable for now.
            // todo: While we have a central body...
            // Iterate, in parallel, over target bodies. The loop over source bodies, per target, is handled
            // by the acceleration function.
            state
                .bodies
                .par_iter_mut()
                .enumerate()
                // .skip(1) // Skip the central body
                .for_each(|(id_target, body_target)| {
                    integrate_rk4(body_target, id_target, &acc, dt);
                });
        }

        if t % BENCH_RATIO == 0 && force_model != ForceModel::GaussShells && !cfg.skip_tree {
            println!(
                "t: {}k, Tree time: {}μs Tree size: {} Integ time: {}μs",
                t / 1_000,
                tree_time,
                tree.as_ref().unwrap().nodes.len(),
                start_time_integ.elapsed().as_micros()
            );
        }

        // Save the current state to a snapshot, for later playback.
        if t % cfg.snapshot_ratio == 0 {
            let nodes: Vec<Cube> = if let Some(t) = &tree {
                if state.ui.draw_tree {
                    // Whole tree
                    // t.nodes.iter().map(|n| n.bounding_box.clone()).collect()

                    // Tree WRT a specific (arbitrary) body.
                    let leaves = t.leaves(state.bodies[0].posit, &state.config.bh_config);

                    leaves.iter().map(|n| n.bounding_box.clone()).collect()
                } else {
                    Vec::new()
                }
            } else {
                Vec::new()
            };
            state.take_snapshot(dt, nodes);
        }
    }

    state.ui.building = false;
    println!("Final V/c: {:.6}", state.bodies[0].vel.magnitude() / C); // todo temp
    println!("Build complete.");
}

fn main() {
    #[cfg(feature = "cuda")]
    let dev = {
        // This is compiled in `build_`.
        let cuda_dev = CudaDevice::new(0).unwrap();
        cuda_dev
            .load_ptx(Ptx::from_file("./cuda.ptx"), "cuda", &["newton_kernel"])
            .unwrap();

        // println!("Using the GPU for computations.");
        ComputationDevice::Gpu(cuda_dev)
    };

    #[cfg(not(feature = "cuda"))]
    let dev = ComputationDevice::Cpu;

    let mut state = State::default();
    if let Ok(cfg) = util::load(&PathBuf::from_str(SAVE_FILE).unwrap()) {
        state.config = cfg;
        println!("V scaler loaded: {:?}", state.config.v_scaler);
    }

    state.ui.dt_input = state.config.dt.to_string();
    state.ui.θ_input = state.config.bh_config.θ.to_string();
    state.ui.v_scaler_input = state.config.v_scaler.to_string();

    state.refresh_bodies();

    // if let Err(e) = util::save(&PathBuf::from(DEFAULT_SNAPSHOT_FILE), &state.snapshots) {
    //     eprintln!("Error saving snapshots: {e}");
    // }

    render(state);
}
