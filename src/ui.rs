use egui::{Color32, Context, RichText, Slider, TopBottomPanel, Ui};
use graphics::{EngineUpdates, Scene};

use crate::{
    build,
    playback::{change_snapshot, SnapShot},
    State,
};

pub const ROW_SPACING: f32 = 22.;
pub const COL_SPACING: f32 = 30.;

fn force_debug(snapshot: &SnapShot, ui: &mut Ui) {
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

        for (i, body_V) in snapshot.body_accs.iter().enumerate() {
            ui.label(format!("Acc at Body {i}:"));
            ui.label(format!("{:?}", body_V));
        }
    });

    ui.horizontal(|ui| {
        for (i, body_p) in snapshot.body_posits.iter().enumerate() {
            ui.label(format!("Posit Body {i}:"));
            ui.label(format!("{:?}", body_p));
        }
    });
}

/// This function draws the (immediate-mode) GUI.
/// [UI items](https://docs.rs/egui/latest/egui/struct.Ui.html)
pub fn ui_handler(state: &mut State, ctx: &Context, scene: &mut Scene) -> EngineUpdates {
    let mut engine_updates = EngineUpdates::default();

    TopBottomPanel::top("0").show(ctx, |ui| {
        let snapshot = if state.snapshots.is_empty() {
            &SnapShot::default()
        } else {
            &state.snapshots[state.ui.snapshot_selected]
        };

        ui.spacing_mut().slider_width = ui.available_width() - 240.;

        ui.horizontal(|ui| {
            ui.label("Snap:");

            let snapshot_prev = state.ui.snapshot_selected;
            ui.add(Slider::new(
                &mut state.ui.snapshot_selected,
                0..=state.snapshots.len() - 1,
            ));

            if state.ui.snapshot_selected != snapshot_prev {
                change_snapshot(&mut scene.entities, snapshot, &state.body_masses);

                engine_updates.entities = true;
            }

            if !state.snapshots.is_empty() {
                ui.add_space(COL_SPACING);
                ui.label(format!(
                    "t: {:.4}",
                    &state.snapshots[state.ui.snapshot_selected].time,
                ));
                ui.label(format!(
                    "dt: {:.6}",
                    &state.snapshots[state.ui.snapshot_selected].dt
                ));
            }
        });

        ui.add_space(ROW_SPACING / 2.);

        force_debug(snapshot, ui);

        ui.horizontal(|ui| {
            // todo: Prog bar
            if ui.button("Build").clicked() {
                build(state, state.ui.force_model);
            }

            if state.ui.building {
                ui.heading(RichText::new("Building...").color(Color32::ORANGE));
            }
        });
    });

    engine_updates
}
