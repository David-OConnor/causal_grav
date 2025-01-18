#![allow(non_ascii_idents)]

//! This model creates distributions of bodies, e.g. ones that coarsely represent galaxies.

use std::f64::consts::TAU;

use lin_alg::{f64::Vec3, linspace};
use rand::Rng;

use crate::{galaxy_data, units::KPC_MYR_PER_KM_S, util::interpolate, Body};

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
    /// X: r (kpc). Y:  M☉ / kpc^2. (todo: Why not / kpc^3?)
    pub mass_density_disk: Vec<(f64, f64)>,
    /// X: r (kpc). Y: km/s. Note: This isn't in our standard units; convert when making bodies.
    pub rotation_curve_disk: Vec<(f64, f64)>,
    /// Luminosity brightness profile. r (kpc), mu (mac arcsec^-2) -
    pub luminosity_disk: Vec<(f64, f64)>,
    pub mass_density_bulge: Vec<(f64, f64)>,
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

    let vol_outer = r_outer.powi(3) * 2. / 3. * TAU;
    let vol_inner = r_inner.powi(3) * 2. / 3. * TAU;

    vol_outer - vol_inner
}

impl GalaxyDescrip {
    fn make_disk(&self, num_bodies: usize, num_rings: usize) -> Vec<Body> {
        if self.mass_density_disk.is_empty() {
            return Vec::new();
        }

        let mut result = Vec::new();
        let mut rng = rand::thread_rng();

        let r_last = self.mass_density_disk.last().unwrap().0;
        let dr = r_last / num_rings as f64;
        let r_all = linspace(0., r_last, num_rings);

        let total_area = r_last.powi(2) * TAU / 2.;

        let mut mass_sample_total = 0.;
        for r in &r_all {
            mass_sample_total += interpolate(&self.mass_density_disk, *r).unwrap();
        }

        // todo: Split into disc + bulge, among other things.
        for r in &r_all {
            // todo: Most of this is a hack.

            let mass_this_r = interpolate(&self.mass_density_disk, *r).unwrap();

            let bodies_this_r = if *r < 1.0e-8 {
                1
            } else {
                let area_this_r = ring_area(*r, dr);

                let area_portion = area_this_r / total_area;
                // let mass_portion = mass_this_r / mass_sample_total;

                // Trying this: Choose number of bodies based on area;
                (area_portion * num_bodies as f64) as usize
            };

            // let bodies_this_r = (mass_portion * num_bodies as f64) as usize;

            println!(
                "N bodies disk. r={r}: {:?} rho: {:?}",
                bodies_this_r, mass_this_r
            );

            // M☉
            let mass_per_body = if bodies_this_r == 0 {
                0.
            } else {
                mass_this_r / bodies_this_r as f64
            };
            println!("MPB disk: {:?}", mass_per_body);

            // Convert from km/s to kpc/myr
            let v_mag = interpolate(&self.rotation_curve_disk, *r).unwrap() * KPC_MYR_PER_KM_S;

            // todo experimenting. It seems the initial v_mag needs to be multiplied by ~2 to
            // prevent the outer bodies from collapsing in. Note that we expect the opposite result
            // from a naive Newton rep, without dark matter.
            // let v_mag = v_mag * 2.;

            for i in 0..bodies_this_r {
                // todo: Random, or even? Even is more uniform, which may be nice, but
                // todo it may cause resonances. Maybe even, but with a random offset per ring?
                let θ = rng.gen_range(0.0..TAU);
                // let θ = TAU / bodies_this_r as f64 * i as f64;

                // Apply eccentricity: Scale radius in x-direction
                let scale_x = 1.0 - self.eccentricity; // Eccentricity factor for x-axis

                let posit = Vec3::new(r * θ.cos() * scale_x, r * θ.sin(), 0.0);

                // todo: Does v need to be a function of theta due to eccentricity?

                // Velocity direction: perpendicular to the radius vector
                let v_x = -v_mag * θ.sin(); // Tangential velocity in x-direction
                let v_y = v_mag * θ.cos(); // Tangential velocity in y-direction

                let vel = Vec3::new(v_x, v_y, 0.0);

                result.push(Body {
                    posit,
                    vel,
                    accel: Vec3::new_zero(),
                    mass: mass_per_body,
                })
            }
        }

        // Scale mass to equal our total
        let mut mass_sum = 0.;
        for body in &result {
            mass_sum += body.mass;
        }

        let mass_scaler = self.mass_disk / mass_sum;
        for body in &mut result {
            body.mass *= mass_scaler;
        }

        // This loop is just diagnostic.
        let mut mass_count = 0.;
        for body in &mut result {
            mass_count += body.mass;
        }

        println!("Total disk body count: {:?}", result.len());
        println!("Total disk mass: {:?}", mass_count);

        result
    }

    fn make_bulge(&self, num_bodies: usize, num_rings: usize) -> Vec<Body> {
        // todo: DRY between this and disk; consolidate.
        // todo: maybe combine these into the same fn. Or, combine common code into a helper.
        if self.mass_density_bulge.is_empty() {
            return Vec::new();
        }

        let mut result = Vec::with_capacity(num_bodies);
        let mut rng = rand::thread_rng();

        let r_last = self.mass_density_bulge.last().unwrap().0;
        let dr = r_last / num_rings as f64;
        let r_all = linspace(0., r_last, num_rings);

        let total_volume = r_last.powi(3) * 2. / 3. * TAU;

        let mut mass_sample_total = 0.;
        for r in &r_all {
            mass_sample_total += interpolate(&self.mass_density_bulge, *r).unwrap();
        }

        println!("\n\nBulge total vol: {:?}\n", total_volume); // todo temp

        // todo: Split into disc + bulge, among other things.
        for r in &r_all {
            // todo: Most of this is a hack.

            let mass_this_r = interpolate(&self.mass_density_bulge, *r).unwrap();

            let bodies_this_r = if *r < 1.0e-8 {
                1
            } else {
                let volume_this_r = ring_volume(*r, dr);

                let volume_portion = volume_this_r / total_volume;

                // Trying this: Choose number of bodies based on area;
                (volume_portion * num_bodies as f64) as usize
            };

            // let bodies_this_r = (mass_portion * num_bodies as f64) as usize;

            // M☉
            let mass_per_body = if bodies_this_r == 0 {
                0.
            } else {
                mass_this_r / bodies_this_r as f64
            };

            println!(
                "Bulge. r: {r} N bodies: {:?} rho: {:.4?}, mass per: {:.4?}",
                bodies_this_r, mass_this_r, mass_per_body
            );

            // Convert from km/s to kpc/myr
            let v_mag = interpolate(&self.rotation_curve_bulge, *r).unwrap() * KPC_MYR_PER_KM_S;

            for i in 0..bodies_this_r {
                // todo: Random, or even? Even is more uniform, which may be nice, but
                // todo it may cause resonances. Maybe even, but with a random offset per ring?
                let θ = rng.gen_range(0.0..TAU);
                let ϕ = rng.gen_range(0.0..TAU / 2.);
                // let θ = TAU / bodies_this_r as f64 * i as f64;

                // Convert spherical coordinates to Cartesian coordinates
                let x = r * ϕ.sin() * θ.cos();
                let y = r * ϕ.sin() * θ.sin();
                let z = r * ϕ.cos();

                // No bulge eccentricity for now.
                let posit = Vec3::new(x, y, z);

                // Velocity direction: perpendicular to the radius vector
                let v_x = v_mag * θ.sin() * ϕ.cos(); // Velocity component in x-direction
                let v_y = v_mag * θ.sin() * ϕ.sin();
                let v_z = v_mag * θ.cos();

                let vel = Vec3::new(v_x, v_y, v_z);

                result.push(Body {
                    posit,
                    vel,
                    accel: Vec3::new_zero(),
                    mass: mass_per_body,
                })
            }
        }

        let mut mass_sum = 0.;
        for body in &result {
            println!("Body mass: {:?}", body.mass);
            mass_sum += body.mass;
        }

        let mass_scaler = self.mass_bulge / mass_sum;
        for body in &mut result {
            body.mass *= mass_scaler;
        }

        // This loop is just diagnostic.
        let mut mass_count = 0.;
        for body in &mut result {
            mass_count += body.mass;
        }

        println!("Total bulge body count: {:?}", result.len());
        println!("Total bulge mass: {:?}", mass_count);

        result
    }

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

        result.append(&mut self.make_disk(num_bodies_disk, num_rings_disk));
        result.append(&mut self.make_bulge(num_bodies_bulge, num_rings_bulge));

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
