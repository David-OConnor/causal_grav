use std::f32::consts::TAU;
use egui::{Context, Slider, ThemePreference, TopBottomPanel};
use graphics::{Camera, ControlScheme, DeviceEvent, EngineUpdates, InputSettings, Lighting, LightType, PointLight, Scene, UiLayout, UiSettings};
use lin_alg::f32::{Quaternion, Vec3};
use crate::State;

type Color = (f32, f32, f32);

const WINDOW_TITLE: &str = "Causal gravity model";
const WINDOW_SIZE_X: f32 = 1_600.;
const WINDOW_SIZE_Y: f32 = 1_200.;
const RENDER_DIST: f32 = 200.;
const BACKGROUND_COLOR: Color = (0.5, 0.5, 0.5);

const SLIDER_WIDTH: f32 = 460.;
const SLIDER_WIDTH_ORIENTATION: f32 = 100.;

pub const ROW_SPACING: f32 = 22.;
pub const COL_SPACING: f32 = 30.;


/// This function draws the (immediate-mode) GUI.
/// [UI items](https://docs.rs/egui/latest/egui/struct.Ui.html#method.heading)
pub fn ui_handler(state: &mut State, ctx: &Context, scene: &mut Scene) -> EngineUpdates {

    let mut engine_updates = EngineUpdates::default();

    TopBottomPanel::top("0").show(ctx, |ui| {
        // todo: Wider values on larger windows?
        // todo: Delegate the numbers etc here to consts etc.
        ui.spacing_mut().slider_width = SLIDER_WIDTH;

        ui.horizontal(|ui| {
            ui.add_space(COL_SPACING);
            ui.label("Time:");
            ui.add(Slider::new(&mut state.ui.snapshot_selected, 0..=state.snapshots.len()));
        });

        ui.add_space(ROW_SPACING / 2.);

    });

    engine_updates
}