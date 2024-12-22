use egui::{Context, Slider, TopBottomPanel};
use graphics::{EngineUpdates, Scene};

use crate::{playback::change_snapshot, State};

pub const ROW_SPACING: f32 = 22.;
pub const COL_SPACING: f32 = 30.;

/// This function draws the (immediate-mode) GUI.
/// [UI items](https://docs.rs/egui/latest/egui/struct.Ui.html)
pub fn ui_handler(state: &mut State, ctx: &Context, scene: &mut Scene) -> EngineUpdates {
    let mut engine_updates = EngineUpdates::default();

    TopBottomPanel::top("0").show(ctx, |ui| {
        ui.spacing_mut().slider_width = ui.available_width() - 200.;

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
                    &state.body_masses,
                );

                engine_updates.entities = true;
            }
        });

        ui.add_space(ROW_SPACING / 2.);

        ui.horizontal(|ui| {
            // for (i, body_V) in state.snapshots[state.ui.snapshot_selected]
            //     .V_at_bodies
            //     .iter()
            //     .enumerate()
            // {
            //     ui.label(&format!("V at Body {i}:"));
            //     ui.label(&format!("{:?}", body_V));
            // }

            ui.add_space(COL_SPACING);

            for (i, body_V) in state.snapshots[state.ui.snapshot_selected]
                .body_accs
                .iter()
                .enumerate()
            {
                ui.label(format!("Acc at Body {i}:"));
                ui.label(format!("{:?}", body_V));
            }
        });

        ui.horizontal(|ui| {
            for (i, body_p) in state.snapshots[state.ui.snapshot_selected]
                .body_posits
                .iter()
                .enumerate()
            {
                ui.label(format!("Posit Body {i}:"));
                ui.label(format!("{:?}", body_p));
            }
        })
    });

    engine_updates
}
