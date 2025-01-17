use egui::{Color32, Context, RichText, Slider, TopBottomPanel, Ui};
use graphics::{EngineUpdates, Entity, Scene};
use lin_alg::f32::{Quaternion, Vec3};

use crate::{
    accel::MondFn,
    barnes_hut::{Cube, Tree},
    body_creation::GalaxyModel,
    build,
    playback::{change_snapshot, SnapShot},
    ForceModel, State, BOUNDING_BOX_PAD,
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
    let mut redraw_bodies = false;

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

            // todo: Combo box for Galaxy model;
            // todo: Populate this as you add data.
            for model in [
                GalaxyModel::Ngc1560,
                // GalaxyModel::Ngc3198,
                // GalaxyModel::Ngc3115,
                // GalaxyModel::Ngc3031,
                // GalaxyModel::Ngc7331,
                GalaxyModel::Ngc2685,
                GalaxyModel::Ngc2824,
            ] {
                if ui
                    .radio_value(&mut state.ui.galaxy_model, model, model.to_str())
                    .changed()
                {
                    state.ui.galaxy_descrip = model.descrip();
                    redraw_bodies = true;
                };
            }

            ui.add_space(COL_SPACING);

            ui.checkbox(&mut state.ui.add_halo, "Add halo");
        });
        ui.horizontal(|ui| {
            ui.add_space(ROW_SPACING);

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

            int_field(
                &mut state.config.num_bodies_disk,
                "bodies disk",
                &mut redraw_bodies,
                ui,
            );
            // int_field(
            //     &mut state.config.num_rings_disk,
            //     "rings disk",
            //     &mut redraw_bodies,
            //     ui,
            // );
            int_field(
                &mut state.config.num_bodies_bulge,
                "bodies bulge",
                &mut redraw_bodies,
                ui,
            );
            // int_field(
            //     &mut state.config.num_rings_bulge,
            //     "rings bulge",
            //     &mut redraw_bodies,
            //     ui,
            // );

            if ui.button("Tree").clicked() {
                // todo: Of current snapshot.
                let bb = Cube::from_bodies(&state.bodies, BOUNDING_BOX_PAD, true).unwrap();

                for body in &state.bodies {
                    println!("Body: {:.4?}", body);
                }

                let tree = Tree::new(&state.bodies, &bb, &state.config.bh_config);

                // todo: Subdivide the tree based on a target body here A/R.

                scene.entities = scene
                    .entities
                    .clone()
                    .into_iter()
                    .filter(|s| s.mesh != 1)
                    .collect();

                for leaf in &tree.nodes {
                    // println!("Node. {}", leaf);
                }

                let leaves = tree.leaves(
                    lin_alg::f64::Vec3::new(2., 2., 0.),
                    99999,
                    &state.config.bh_config,
                );
                println!("Leaf count: {:?}", leaves.len());

                for leaf in leaves {
                    let c = leaf.bounding_box.center;
                    let posit =
                        Vec3::new(c.x as f32, c.y as f32, c.z as f32) + Vec3::new(0., 0., 1.5);

                    scene.entities.push(Entity::new(
                        1,
                        posit,
                        Quaternion::new_identity(),
                        leaf.bounding_box.width as f32 * 0.85,
                        (50., 100., 255.),
                        1.,
                    ));
                }

                engine_updates.entities = true;
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

    if redraw_bodies {
        // todo: We don't need to change rings for all cases of `redraw_bodies`, but this is harmless
        state.config.num_rings_disk = state.config.num_bodies_disk / 10; // todo: Delegate this.
        state.config.num_rings_bulge = state.config.num_bodies_bulge / 10; // todo: Delegate this.

        state.bodies = state.ui.galaxy_descrip.make_bodies(
            state.config.num_bodies_disk,
            state.config.num_rings_disk,
            state.config.num_bodies_bulge,
            state.config.num_rings_bulge,
        );
        state.body_masses = state.bodies.iter().map(|b| b.mass as f32).collect();

        state.snapshots = Vec::new();
        state.take_snapshot(0., 0.); // Initial snapshot; t=0.

        reset_snapshot = true;
        engine_updates.entities = true;
    }

    if reset_snapshot {
        change_snapshot(&mut scene.entities, &state.snapshots[0], &state.body_masses);
    }

    engine_updates
}
