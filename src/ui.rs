use egui::{Context, Slider, TopBottomPanel};
use graphics::{EngineUpdates, Scene};

use crate::{playback::change_snapshot, State};

const SLIDER_WIDTH: f32 = 460.;
const SLIDER_WIDTH_ORIENTATION: f32 = 100.;

pub const ROW_SPACING: f32 = 22.;
pub const COL_SPACING: f32 = 30.;

/// This function draws the (immediate-mode) GUI.
/// [UI items](https://docs.rs/egui/latest/egui/struct.Ui.html)
pub fn ui_handler(state: &mut State, ctx: &Context, scene: &mut Scene) -> EngineUpdates {
    let mut engine_updates = EngineUpdates::default();

    TopBottomPanel::top("0").show(ctx, |ui| {
        ui.spacing_mut().slider_width = SLIDER_WIDTH;

        ui.horizontal(|ui| {
            ui.add_space(COL_SPACING);
            ui.label("Time:");

            let snapshot_prev = state.ui.snapshot_selected;
            ui.add(Slider::new(
                &mut state.ui.snapshot_selected,
                0..=state.snapshots.len() - 1,
            ));

            if state.ui.snapshot_selected != snapshot_prev {
                change_snapshot(
                    &mut scene.entities,
                    &state.snapshots[state.ui.snapshot_selected],
                );

                engine_updates.entities = true;
            }
        });

        ui.add_space(ROW_SPACING / 2.);

        ui.horizontal(|ui| {
            for (i, body_V) in state.snapshots[state.ui.snapshot_selected]
                .V_at_bodies
                .iter()
                .enumerate()
            {
                ui.label(&format!("V at Body {i}:"));
                ui.label(&format!("{:?}", body_V));
            }
        });
    });

    engine_updates
}
