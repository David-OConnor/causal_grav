use egui::{Color32, Context, RichText, Slider, TopBottomPanel, Ui};
use graphics::{EngineUpdates, Scene};

use crate::{
    accel::MondFn,
    body_creation::GalaxyModel,
    build,
    playback::{change_snapshot, SnapShot},
    ForceModel, State,
};

pub const ROW_SPACING: f32 = 10.;
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

    // This variable prevents mutliple borrow errors.
    let mut reset_snapshot = false;

    TopBottomPanel::top("0").show(ctx, |ui| {
        let snapshot = if state.snapshots.is_empty() {
            &SnapShot::default()
        } else {
            &state.snapshots[state.ui.snapshot_selected]
        };

        ui.spacing_mut().slider_width = ui.available_width() - 280.;

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

        ui.add_space(ROW_SPACING);

        // force_debug(snapshot, ui);

        ui.horizontal(|ui| {
            // todo: Prog bar
            if ui
                .button(RichText::new("Build").color(Color32::GOLD))
                .clicked()
            {
                build(state, state.ui.force_model);
            }

            if state.ui.building {
                ui.heading(RichText::new("Building...").color(Color32::ORANGE));
            }

            ui.add_space(COL_SPACING);

            ui.radio_value(&mut state.ui.force_model, ForceModel::Newton, "Newton");
            ui.radio_value(
                &mut state.ui.force_model,
                ForceModel::Mond(MondFn::Simple),
                "MOND simple",
            );
            ui.radio_value(
                &mut state.ui.force_model,
                ForceModel::Mond(MondFn::Standard),
                "MOND",
            );
            ui.radio_value(
                &mut state.ui.force_model,
                ForceModel::GaussShells,
                "Causal shells",
            );

            ui.add_space(COL_SPACING);

            // todo: Combo box for Galaxy model;
            // todo: Populate this as you add data.
            for model in [
                GalaxyModel::Ngc1560,
                GalaxyModel::Ngc3198,
                // GalaxyModel::Ngc3115,
                // GalaxyModel::Ngc3031,
                // GalaxyModel::Ngc7331,
                GalaxyModel::Ngc2685,
            ] {
                if ui
                    .radio_value(&mut state.ui.galaxy_model, model, model.to_str())
                    .changed()
                {
                    state.bodies = state.ui.galaxy_model.make_bodies();
                    state.ui.galaxy_descrip = model.descrip();
                    state.body_masses = state.bodies.iter().map(|b| b.mass as f32).collect();

                    state.snapshots = Vec::new();
                    state.take_snapshot(0., 0.); // Initial snapshot; t=0.

                    reset_snapshot = true;
                    engine_updates.entities = true;
                };
            }

            ui.add_space(COL_SPACING);

            ui.checkbox(&mut state.ui.add_halo, "Add halo");

            ui.add_space(COL_SPACING);

            ui.label("dt:");
            // ui.text_edit_singleline(&mut state.ui.dt_input);
            ui.add_sized(
                [60., Ui::available_height(ui)],
                egui::TextEdit::singleline(&mut state.ui.dt_input),
            );
            if ui.button("Save dt").clicked() {
                if let Ok(v) = state.ui.dt_input.parse() {
                    state.config.dt = v;
                }
            }

            ui.add_space(COL_SPACING);

            ui.label("Steps (x1000):");
            // todo: Width?
            let mut val = (state.config.num_timesteps / 1_000).to_string();
            if ui
                .add_sized(
                    [40., Ui::available_height(ui)],
                    egui::TextEdit::singleline(&mut val),
                )
                .changed()
            {
                if let Ok(v) = val.parse::<usize>() {
                    state.config.num_timesteps = v * 1_000;
                }
            }

            // todo: Other params like ring ratio.
        });
        ui.add_space(ROW_SPACING);

        ui.horizontal(|ui| {
            let desc = &state.ui.galaxy_descrip;
            ui.label(format!("Mass: {} ×10⁸ M☉", desc.mass_disk / 1.0e8));
            ui.add_space(COL_SPACING);
            ui.label(format!("M/L: {}", desc.mass_to_light_ratio)); // todo: Remove A/R
            ui.add_space(COL_SPACING);
            ui.label(format!("Dist: {} kpc", desc.dist_from_earth));
            ui.add_space(COL_SPACING);
            ui.label(format!("Eccentricity: {}", desc.eccentricity));
            ui.add_space(COL_SPACING);
        });

        ui.add_space(ROW_SPACING);
    });

    if reset_snapshot {
        change_snapshot(&mut scene.entities, &state.snapshots[0], &state.body_masses);
    }

    engine_updates
}
