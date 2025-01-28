//! This module integraties this application with the graphics engine.

use std::f32::consts::TAU;

use graphics::{
    Camera, ControlScheme, DeviceEvent, EngineUpdates, InputSettings, LightType, Lighting,
    Mesh, PointLight, Scene, UiLayout, UiSettings,
};
use lin_alg::f32::{Quaternion, Vec3};

use crate::{playback::change_snapshot, ui::ui_handler, State};

type Color = (f32, f32, f32);

const WINDOW_TITLE: &str = "Causal gravity model";
const WINDOW_SIZE_X: f32 = 1_600.;
const WINDOW_SIZE_Y: f32 = 1_000.;
const BACKGROUND_COLOR: Color = (0., 0., 0.);

const RENDER_DIST: f32 = 200.;

pub const BODY_SIZE_SCALER: f32 = 1.0e-9; // We multiply mass by this.
pub const BODY_SIZE_MIN: f32 = 0.01;
pub const BODY_SIZE_MAX: f32 = 0.6;

pub const BODY_COLOR: Color = (1.0, 0.4, 0.4);
pub const BODY_SHINYNESS: f32 = 2.;

pub const SHELL_COLOR: Color = (1.0, 0.6, 0.2);
pub const SHELL_SHINYNESS: f32 = 2.;

pub const TREE_COLOR: Color = (0.4, 0.4, 1.0);
pub const TREE_SHINYNESS: f32 = 1.;

pub const ARROW_COLOR: Color = (0.2, 1.0, 0.6);
pub const ARROW_SHINYNESS: f32 = 1.;

// Allows individual cubes to be distinguished by creating gaps between them.
pub const TREE_CUBE_SCALE_FACTOR: f32 = 0.85;

// Keep this in sync with mesh init.
pub const MESH_SPHERE: usize = 0;
pub const MESH_CUBE: usize = 1;
pub const MESH_ARROW: usize = 2;

pub const SHELL_OPACITY: f32 = 0.01;

fn event_handler(
    _state: &mut State,
    _event: DeviceEvent,
    _scene: &mut Scene,
    _dt: f32,
) -> EngineUpdates {
    EngineUpdates::default()
}

/// This runs each frame. Currently, no updates.
fn render_handler(_state: &mut State, _scene: &mut Scene, _dt: f32) -> EngineUpdates {
    EngineUpdates::default()
}

/// Entry point to our render and event loop.
pub fn render(state: State) {
    // Initialize entities.
    let mut entities = Vec::new();
    change_snapshot(
        &mut entities,
        &state.snapshots[state.ui.snapshot_selected],
        &state.body_masses,
    );

    let scene = Scene {
        meshes: vec![
            Mesh::new_sphere(1., 12, 12),
            Mesh::new_box(1., 1., 1.),
            Mesh::new_arrow(1., 0.05, 8),
        ],
        entities,
        camera: Camera {
            fov_y: TAU / 8.,
            position: Vec3::new(0., 0., -60.),
            far: RENDER_DIST,
            near: 0.2, // todo: Adjust A/R
            // orientation: Quaternion::from_axis_angle(Vec3::new(1., 0., 0.), TAU / 16.),
            orientation: Quaternion::from_axis_angle(Vec3::new(1., 0., 0.), 0.),
            ..Default::default()
        },
        lighting: Lighting {
            ambient_color: [-1., 1., 1., 0.5],
            ambient_intensity: 0.02,
            point_lights: vec![
                // Light from above
                PointLight {
                    type_: LightType::Omnidirectional,
                    position: Vec3::new(20., 20., 100.),
                    diffuse_color: [0.3, 0.4, 0.5, 1.],
                    specular_color: [0.3, 0.4, 0.5, 1.],
                    diffuse_intensity: 1_000.,
                    specular_intensity: 4_000.,
                },
                // Light from below
                PointLight {
                    type_: LightType::Omnidirectional,
                    position: Vec3::new(-20., 20., -100.),
                    diffuse_color: [0.3, 0.4, 0.5, 1.],
                    specular_color: [0.3, 0.4, 0.5, 1.],
                    diffuse_intensity: 1_000.,
                    specular_intensity: 4_000.,
                },
            ],
        },
        background_color: BACKGROUND_COLOR,
        window_size: (WINDOW_SIZE_X, WINDOW_SIZE_Y),
        window_title: WINDOW_TITLE.to_owned(),
    };

    let input_settings = InputSettings {
        initial_controls: ControlScheme::FreeCamera,
        move_sens: 3.5,
        ..Default::default()
    };
    let ui_settings = UiSettings {
        layout: UiLayout::Top,
        icon_path: Some("./resources/icon.png".to_owned()),
    };

    graphics::run(
        state,
        scene,
        input_settings,
        ui_settings,
        render_handler,
        event_handler,
        ui_handler,
    );
}
