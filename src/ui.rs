use std::{collections::HashMap, path::PathBuf, str::FromStr};

use barnes_hut::{Cube, Tree};
use egui::{Color32, ComboBox, Context, RichText, Slider, TopBottomPanel, Ui};
use graphics::{EngineUpdates, Entity, Scene};
use lin_alg::{
    f32::{Quaternion, Vec3},
    f64::Vec3 as Vec3F64,
    linspace,
};

use crate::{
    accel::MondFn,
    build,
    charge::{plot_field_properties, FieldProperties},
    galaxy_data::GalaxyModel,
    playback::{change_snapshot, SnapShot},
    render::{TREE_COLOR, TREE_CUBE_SCALE_FACTOR, TREE_SHINYNESS},
    util, ForceModel, State, BOUNDING_BOX_PAD, SAVE_FILE,
};

pub const ROW_SPACING: f32 = 10.;
pub const COL_SPACING: f32 = 30.;

fn int_field(val: &mut usize, label: &str, redraw_bodies: &mut bool, ui: &mut Ui) {
    ui.label(label);
    let mut val_str = val.to_string();
    if ui
        .add_sized(
            [60., Ui::available_height(ui)],
            egui::TextEdit::singleline(&mut val_str),
        )
        .changed()
    {
        if let Ok(v) = val_str.parse::<usize>() {
            *val = v;
            *redraw_bodies = true;
        }
    }
}

/// This function draws the (immediate-mode) GUI.
/// [UI items](https://docs.rs/egui/latest/egui/struct.Ui.html)
pub fn ui_handler(state: &mut State, ctx: &Context, scene: &mut Scene) -> EngineUpdates {
    let mut engine_updates = EngineUpdates::default();

    // This variable prevents mutliple borrow errors.
    let mut reset_snapshot = false;
    let mut refresh_bodies = false;

    TopBottomPanel::top("0").show(ctx, |ui| {
        let snapshot = if state.snapshots.len() < state.ui.snapshot_selected {
            state.ui.snapshot_selected = 0;
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

            let mut prev_model = state.ui.galaxy_model;
            ComboBox::from_id_salt(0)
                .width(120.)
                .selected_text(state.ui.galaxy_model.to_str())
                .show_ui(ui, |ui| {
                    for model in [
                        GalaxyModel::Ngc1560,
                        GalaxyModel::Ngc2685,
                        GalaxyModel::Ngc2824,
                        GalaxyModel::Ngc3626,
                        GalaxyModel::Ugc6176,
                        GalaxyModel::M31,
                    ] {
                        ui.selectable_value(&mut state.ui.galaxy_model, model, model.to_str());
                    }
                });
            if prev_model != state.ui.galaxy_model {
                state.ui.galaxy_descrip = state.ui.galaxy_model.descrip();
                refresh_bodies = true;
            }

            ui.add_space(COL_SPACING);

            ui.checkbox(&mut state.ui.add_halo, "Add halo");
        });

        ui.add_space(ROW_SPACING);

        ui.horizontal(|ui| {
            ui.label("dt:");
            ui.add_sized(
                [60., Ui::available_height(ui)],
                egui::TextEdit::singleline(&mut state.ui.dt_input),
            );
            if ui.button("Save dt").clicked() {
                if let Ok(v) = state.ui.dt_input.parse() {
                    state.config.dt = v;
                }
            }

            ui.label("Steps (x1000):");
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

            ui.label("θ:");
            ui.add_sized(
                [40., Ui::available_height(ui)],
                egui::TextEdit::singleline(&mut state.ui.θ_input),
            );

            if ui.button("Save θ").clicked() {
                if let Ok(v) = state.ui.θ_input.parse() {
                    state.config.bh_config.θ = v;
                }
            }

            ui.label("v scaler:");
            ui.add_sized(
                [36., Ui::available_height(ui)],
                egui::TextEdit::singleline(&mut state.ui.v_scaler_input),
            );
            if ui.button("Save v scaler").clicked() {
                if let Ok(v) = state.ui.v_scaler_input.parse() {
                    state.config.v_scaler = v;
                    refresh_bodies = true; // For updated v.
                }
            }

            int_field(
                &mut state.config.num_bodies_disk,
                "bodies disk",
                &mut refresh_bodies,
                ui,
            );

            int_field(
                &mut state.config.num_bodies_bulge,
                "bodies bulge",
                &mut refresh_bodies,
                ui,
            );

            // todo: Remove A/R now that cube is in snapshots.
            if ui.button("Tree").clicked() {
                // todo: Of current snapshot.
                let bb = Cube::from_bodies(&state.bodies, BOUNDING_BOX_PAD, true).unwrap();
                //
                // for body in &state.bodies {
                //     println!("Body: {:.4?}", body);
                // }

                let tree = Tree::new(&state.bodies, &bb, &state.config.bh_config);

                // todo: Subdivide the tree based on a target body here A/R.

                scene.entities = scene
                    .entities
                    .clone()
                    .into_iter()
                    .filter(|s| s.mesh != 1)
                    .collect();

                let leaves =
                    tree.leaves(lin_alg::f64::Vec3::new(2., 2., 0.), &state.config.bh_config);
                println!(
                    "Leaf count: {}. Full tree len: {}",
                    leaves.len(),
                    tree.nodes.len()
                );

                // for leaf in &leaves { // todo te p
                //     println!("Leaf: {:?}", leaf);
                // }

                for leaf in leaves {
                    let c = leaf.bounding_box.center;
                    let posit =
                        // Vec3::new(c.x as f32, c.y as f32, c.z as f32) + Vec3::new(0., 0., 1.5);
                        Vec3::new(c.x as f32, c.y as f32, c.z as f32);

                    scene.entities.push(Entity::new(
                        1,
                        posit,
                        Quaternion::new_identity(),
                        leaf.bounding_box.width as f32 * TREE_CUBE_SCALE_FACTOR,
                        TREE_COLOR,
                        TREE_SHINYNESS,
                    ));
                }

                engine_updates.entities = true;
            }

            ui.add_space(COL_SPACING);

            ui.checkbox(&mut state.config.skip_tree, "Skip tree");

            ui.checkbox(&mut state.ui.draw_tree, "Draw tree");

            ui.add_space(COL_SPACING * 2.);

            if ui.button("Field properties").clicked() {
                let dx = 0.5;
                let mut properties = HashMap::new();
                for r in linspace(0., 2.5, 12) {
                    let point = Vec3F64::new(r, 0., 0.);

                    // todo: After running, this will be the final config.
                    // let bodies = &state.snapshots[state.ui.snapshot_selected].b
                    let bodies = &state.bodies;

                    let stats = FieldProperties::new(bodies, point, dx);
                    println!("\nStats at R={r}: {stats}");
                    properties.insert(r, stats);
                }
                plot_field_properties(&properties);
            }

            if ui
                .button(RichText::new("Save").color(Color32::GOLD))
                .clicked()
            {
                if util::save(&PathBuf::from_str(SAVE_FILE).unwrap(), &state.config).is_err() {
                    println!("Error saving config.")
                }
            }
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

    if refresh_bodies {
        reset_snapshot = true;
        engine_updates.entities = true;

        state.refresh_bodies()
    }

    if reset_snapshot {
        change_snapshot(&mut scene.entities, &state.snapshots[0], &state.body_masses);
    }

    engine_updates
}
