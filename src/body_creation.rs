#![allow(non_ascii_idents)]

//! This model creates distributions of bodies, e.g. ones that coarsely represent galaxies.

use std::f64::consts::TAU;

use lin_alg::{f64::Vec3, linspace};
use rand::Rng;

use crate::{
    galaxy_data,
    util::{interpolate, volume_sphere},
    Body,
};

/// Create n bodies, in circular orbit, at equal distances from each other.
pub fn make_bodies_balanced(num: usize, r: f64, mass_body: f64, mass_central: f64) -> Vec<Body> {
    let mut result = Vec::with_capacity(num);

    for i in 0..num {
        let θ = TAU / num as f64 * i as f64;
        let posit = Vec3::new(r * θ.cos(), r * θ.sin(), 0.0);

        // Velocity magnitude for circular orbit
        let v_mag = (mass_central / r).sqrt();

        // Velocity direction: perpendicular to the radius vector
        let v_x = -v_mag * θ.sin(); // Tangential velocity in x-direction
        let v_y = v_mag * θ.cos(); // Tangential velocity in y-direction

        let vel = Vec3::new(v_x, v_y, 0.0);

        result.push(Body {
            posit,
            vel,
            accel: Vec3::new_zero(),
            mass: mass_body,
        });
    }

    result
}

#[derive(Clone, Copy, PartialEq)]
pub enum GalaxyShape {
    GrandDesignSpiral,
    FlocculentSpiral,
    MultiArmSpiral,
    BarredSpiral,
    Lenticular,
    Elliptical,
    LenticularRingSeyfertType2,
}

/// todo: We assume a spiral galaxy for now
pub struct GalaxyDescrip {
    pub shape: GalaxyShape,
    /// X: r (kpc). Y:  M☉ / kpc^2.
    pub mass_density_disk: Vec<(f64, f64)>,
    // /// X: r (kpc). Y: km/s. Note: This isn't in our standard units; convert when making bodies.
    /// X: r (kpc). Y: kpc/MYR.
    pub rotation_curve_disk: Vec<(f64, f64)>,
    /// Luminosity brightness profile. r (kpc), mu (mac arcsec^-2) -
    pub luminosity_disk: Vec<(f64, f64)>,
    /// X: r (kpc). Y:  M☉ / kpc^2. Note that this is only valid in the plane of the bulge; you must
    /// map to a 3D structure using data not included here.
    pub mass_density_bulge: Vec<(f64, f64)>,
    /// X: r (kpc). Y: kpc/MYR.
    pub rotation_curve_bulge: Vec<(f64, f64)>,
    pub luminosity_bulge: Vec<(f64, f64)>,
    // todo: More A/R
    /// 0 means a circle. 1 is fully elongated.
    pub eccentricity: f64,
    pub arm_count: usize,
    /// For generating a dark matter halo. (core radius, central density)
    pub burkert_params: (f64, f64),
    /// Not a fundamental property; used to normalize mass density etc?
    /// todo: I'm not sure what this is
    pub r_s: f64,
    /// M☉
    pub mass_bulge: f64,
    /// M☉
    pub mass_disk: f64,
    /// M/L_B Used to convert luminosity to mass density. Solar masses / Mu ?
    pub mass_to_light_ratio: f64,
    /// Kpc
    pub dist_from_earth: f64,
}

fn ring_area(r: f64, dr: f64) -> f64 {
    let r_outer = r + dr / 2.;
    let r_inner = r - dr / 2.;

    let area_outer = r_outer.powi(2) * TAU / 2.;
    let area_inner = r_inner.powi(2) * TAU / 2.;

    area_outer - area_inner
}

fn ring_volume(r: f64, dr: f64) -> f64 {
    let r_outer = r + dr / 2.;
    let r_inner = r - dr / 2.;

    let vol_outer = volume_sphere(r_outer);
    let vol_inner = volume_sphere(r_inner);

    vol_outer - vol_inner
}

impl GalaxyDescrip {
    /// See the `properties` module for info on distributions
    /// todo: Luminosity A/R
    pub fn make_bodies(
        &self,
        num_bodies_disk: usize,
        num_rings_disk: usize,
        num_bodies_bulge: usize,
        num_rings_bulge: usize,
    ) -> Vec<Body> {
        let mut result = Vec::with_capacity(num_bodies_disk + num_bodies_bulge);

        // result.append(&mut self.make_disk(num_bodies_disk, num_rings_disk));
        println!("\nMaking disk bodies...");
        result.append(&mut make_distrib(
            &self.mass_density_disk,
            &self.rotation_curve_disk,
            self.mass_disk,
            self.eccentricity,
            num_bodies_disk,
            num_rings_disk,
            false,
        ));

        // println!("Bodies: {:.4?}", &result);

        // result.append(&mut self.make_distributation(num_bodies_bulge, num_rings_bulge));

        println!("\nMaking bulge bodies...");
        if num_bodies_bulge > 0 && !self.mass_density_bulge.is_empty() {
            result.append(&mut make_distrib(
                &self.mass_density_bulge,
                &self.rotation_curve_bulge,
                self.mass_bulge,
                self.eccentricity,
                num_bodies_bulge,
                num_rings_bulge,
                // todo: Temp flat.
                true,
            ));
        }

        result
    }
}

/// todo: Move specific galaxy creation to its own module A/R
#[derive(Clone, Copy, PartialEq)]
pub enum GalaxyModel {
    Ngc1560,
    Ngc3198,
    Ngc3115,
    Ngc3031,
    Ngc7331,
    Ngc2685,
    Ngc2824,
}

impl Default for GalaxyModel {
    fn default() -> Self {
        Self::Ngc1560
    }
}

impl GalaxyModel {
    pub fn to_str(&self) -> String {
        match self {
            Self::Ngc1560 => "NGC 1560",
            Self::Ngc3198 => "NGC 3198",
            Self::Ngc3115 => "NGC 3115",
            Self::Ngc3031 => "NGC 3031",
            Self::Ngc7331 => "NGC 7331",
            Self::Ngc2685 => "NGC 2685",
            Self::Ngc2824 => "NGC 2824",
        }
        .to_owned()
    }

    pub fn descrip(&self) -> GalaxyDescrip {
        match self {
            /// Ludwig, Figures 3 and 5. todo: Partial/rough
            Self::Ngc1560 => galaxy_data::ngc_1560(),
            Self::Ngc3198 => galaxy_data::ngc_3198(),
            Self::Ngc3115 => galaxy_data::ngc_3115(),
            Self::Ngc2685 => galaxy_data::ngc_2685(),
            Self::Ngc2824 => galaxy_data::ngc_2824(),
            _ => unimplemented!(), // todo
        }
    }

    // pub fn make_bodies(&self, num_bodies_disk: usize, num_rings_disk: usize, num_bodies_disk: usize, num_rings_disk: usize,) -> Vec<Body> {
    //     self.descrip().make_bodies(num_bodies, num_rings)
    // }
}

/// Create mass density from luminosity. X axis for both is r (distance from the galactic center).
pub fn mass_density_from_lum(
    luminosity: &[(f64, f64)],
    mass_total: f64,
    luminosity_arcsec: &[(f64, f64)],
) -> Vec<(f64, f64)> {
    // todo: Come back to this; re-examine converting surface brightness to solar luminosity etc.
    // todo: You are likely missing a step.
    let mut mass_density = Vec::with_capacity(luminosity.len());
    for (i, (r, lum)) in luminosity.iter().enumerate() {
        // surface brightness profile (μ)
        // μ = mag / arcsec^2
        // mag = μ * arcsec^2
        // todo. Hmm: Why are we then, dividing?
        let mut arcsec_sq = luminosity_arcsec[i].0.powi(2);
        if arcsec_sq < 1e-14 {
            arcsec_sq = 1.; // avoid div by 0. // todo: Kludge
        }
        // mass_density.push((*r, mass_to_light_ratio * r / arcsec_sq))
    }

    // Let's try a simple technique where, to find mass density, we scale luminosity by the total mass.
    // todo: You really need to QC this!
    let mut lum_total = 0.;
    for (r, lum) in luminosity {
        if *r < 1e-14 {
            continue; // todo kludege
        }
        lum_total += lum / r.powi(2);
    }

    let lum_scaler = mass_total / lum_total;
    for (r, lum) in luminosity {
        if *r < 1e-14 {
            continue; // todo kludege
        }
        mass_density.push((*r, lum_scaler * lum / r.powi(2)));
    }

    // todo temp
    // for (r, lum) in &luminosity {
    //     println!("R: {r} lum: {lum}");
    // }

    // for (r, mass) in &mass_density{
    //     println!("R: {r} mass: {mass}");
    // }

    mass_density
}

fn make_distrib(
    mass_density: &[(f64, f64)],
    vel: &[(f64, f64)],
    mass_total: f64,
    eccentricity: f64,
    num_bodies: usize,
    num_rings: usize,
    three_d: bool,
) -> Vec<Body> {
    // fn make_bulge(&self, num_bodies: usize, num_rings: usize) -> Vec<Body> {
    // todo: DRY between this and disk; consolidate.
    // todo: maybe combine these into the same fn. Or, combine common code into a helper.

    let mut result = Vec::with_capacity(num_bodies);
    let mut rng = rand::thread_rng();

    // todo: Experiment with this.
    let min_bodies_per_ring_area = 0.4;

    // Note: We are choosing the r range here based on mass density; not velocity.
    let r_last = mass_density.last().unwrap().0;
    let r_all = linspace(mass_density[0].0, r_last, num_rings);
    let dr = r_all[1] - r_all[0];

    let r_len = r_all.len();

    let total_area = r_last.powi(2) * TAU / 2.;
    let total_volume = volume_sphere(r_last);

    // let mut mass_sample_total = 0.;
    // for r in &r_all {
    //     mass_sample_total += interpolate(mass_density, *r).unwrap();
    // }


    // Compute the minimum number of bodies at each r
    // let bodies_per_r_init: Vec<f64> = r_all
    //     .iter()
    //     .map(|r|
    //     .collect();

    // // Set mass proportionally to initial body numbers.
    // let masses_per_r_init: Vec<f64> = mass_density.iter().enumerate().map(|(i, (r, mass))| {
    //     bodies_per_r_init
    // }).collect();

    let mut body_num_by_r_init = Vec::with_capacity(r_len);
    let mut mass_per_body_by_r = Vec::with_capacity(r_len);
    let mut num_bodies_init = 0.;


    // for (r, mass) in &mass_density {
    for r in &r_all {
        // Compute the minimum number of bodies at each r
        let body_count_min = r.powi(2)* TAU / 2. * min_bodies_per_ring_area;

        body_num_by_r_init.push(body_count_min);

        // Set mass proportionally to initial body numbers.
        mass_per_body_by_r.push(interpolate(mass_density, *r).unwrap() / body_count_min);

        num_bodies_init += body_count_min;
    }

    let body_n_ratio = num_bodies as f64 / num_bodies_init;

    let mut bodies_by_r = Vec::with_capacity(r_len);

    for (i, r) in r_all.iter().enumerate() {
        // Scale body count to match target.
        bodies_by_r.push((body_num_by_r_init[i] / body_n_ratio) as usize);

        // Scale mass proportionally.
        // todo: QC this.
        // Note that we we use the final body number integer to scale mass.
        mass_per_body_by_r[i] /= body_num_by_r_init[i] / bodies_by_r[i] as f64
    }

    // todo: Experimenting with fixed mass.
    // let mass_per_body = mass_total / num_bodies as f64;
    // println!("Mass per body: {:.6}", mass_per_body);

    // let mut n_bodies_by_r = Vec::with_capacity(r_len);
    // for r in &r_all {
    //     // let mass_this_r = interpolate(mass_density, *r).unwrap();
    //
    //
    //     // todo: we go below and above for this. Check we're not below dr/2?
    //     // todo: we go below and above for this. Check we're not below dr/2?
    //     // let bodies_this_r = if *r < 1.0e-8 {
    //
    //     let bodies_this_r = if three_d {
    //         if *r < dr / 2. {
    //             1
    //         } else {
    //             let volume_this_r = ring_volume(*r, dr);
    //
    //             let volume_portion = volume_this_r / total_volume;
    //
    //             // Trying this: Choose number of bodies based on volume;
    //             (volume_portion * mass_this_r) as usize
    //
    //             // (mass_this_r / volume_this_r) as usize // todo:
    //
    //             // todo: Experimenting.
    //             // (100. * volume_portion * num_bodies as f64) as usize
    //         }
    //     } else {
    //         if *r < dr / 2. {
    //             1
    //         } else {
    //             let area_this_r = ring_area(*r, dr);
    //             let area_portion = area_this_r / total_area;
    //             // let mass_portion = mass_this_r / mass_sample_total;
    //
    //             // Trying this: Choose number of bodies based on area;
    //             (area_portion * num_bodies as f64) as usize
    //         }
    //     };
    //
    //     // println!("Bodies at r{:.4}: {}", r, bodies_this_r);
    //
    //     // let bodies_this_r = (mass_portion * num_bodies as f64) as usize;
    //
    //     // M☉
    //     // let mass_per_body_ = if bodies_this_r == 0 {
    //     //     0.
    //     // } else {
    //     //     mass_this_r / bodies_this_r as f64
    //     // };
    //
    //     n_bodies_by_r.push(bodies_this_r);
    // }

    // Normalize body count.
    // let count: usize = n_bodies_by_r.iter().sum();

    for (i, r) in r_all.iter().enumerate() {
        // let bodies_this_r = n_bodies_by_r[i] * num_bodies / count;
        // let bodies_this_r = n_bodies_by_r[i];
        // todo: Dup with above!
        // let mass_this_r = interpolate(mass_density, *r).unwrap();

        println!(
            "Body data. r: {r} N bodies: {:?} mass-per-body: {:.4?}, mass-this-r: {:.4?}",
            bodies_by_r[i], mass_per_body_by_r[i], interpolate(mass_density, *r).unwrap()
        );

        let v_mag = interpolate(vel, *r).unwrap();

        for _ in 0..bodies_by_r[i] {
            let θ = rng.gen_range(0.0..TAU);

            let (posit, vel) = if three_d {
                let ϕ = if three_d {
                    // Random phi for polar angle with area weighting
                    let u: f64 = rng.gen_range(-1.0..1.0); // Uniform random variable
                    u.acos() // Inverse cosine for area-preserving sampling
                } else {
                    0.
                };

                // Convert spherical coordinates to Cartesian coordinates
                let x = r * ϕ.sin() * θ.cos();
                let y = r * ϕ.sin() * θ.sin();
                let z = if three_d { r * ϕ.cos() } else { 0. };

                let scale_x = 1.0 - eccentricity; // Eccentricity factor for x-axis
                let posit = Vec3::new(x * scale_x, y, z);

                // Velocity direction: perpendicular to the radius vector
                let v_x = θ.sin() * ϕ.cos(); // Velocity component in x-direction
                let v_y = θ.sin() * ϕ.sin();
                let v_z = θ.cos();
                let mut vel = Vec3::new(v_x / scale_x, v_y, v_z);

                // Normalize velocity vector and scale to v_mag
                let vel = vel.to_normalized() * v_mag;

                (posit, vel)
            } else {
                let scale_x = 1.0 - eccentricity; // Eccentricity factor for x-axis
                let posit = Vec3::new(r * θ.cos() * scale_x, r * θ.sin(), 0.0);

                // Velocity direction: perpendicular to the radius vector
                let v_x = -v_mag * θ.sin(); // Tangential velocity in x-direction
                let v_y = v_mag * θ.cos(); // Tangential velocity in y-direction

                let vel = Vec3::new(v_x, v_y, 0.0);

                (posit, vel)
            };

            result.push(Body {
                posit,
                vel,
                accel: Vec3::new_zero(),
                mass: mass_per_body_by_r[i]
            })
        }
    }

    let mut mass_sum = 0.;
    for body in &result {
        mass_sum += body.mass;
    }

    let mass_scaler = mass_total / mass_sum;
    for body in &mut result {
        body.mass *= mass_scaler;
    }

    // This loop is just diagnostic.
    let mut mass_count = 0.;
    for body in &mut result {
        mass_count += body.mass;
    }

    println!("Total bodies {:?}", result.len());
    println!("Total mass: {:.0?} e9", mass_count / 1e9);

    result
}
