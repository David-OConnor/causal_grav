#![allow(non_ascii_idents)]

//! This model creates distributions of bodies, e.g. ones that coarsely represent galaxies.

use std::f64::consts::TAU;

use lin_alg::{f64::Vec3, linspace};
use lin_alg::f64::Quaternion;
use rand::Rng;
use rand::rngs::ThreadRng;
use crate::{Body, DISK_RING_PORTION, util::{interpolate, volume_sphere}};

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
        num_bodies_bulge: usize,
        v_scaler: f64,
    ) -> Vec<Body> {
        let mut result = Vec::with_capacity(num_bodies_disk + num_bodies_bulge);

        // result.append(&mut self.make_disk(num_bodies_disk, num_rings_disk));
        println!("\nMaking disk bodies...");
        result.append(&mut make_distrib_data_area(
            &self.mass_density_disk,
            &self.rotation_curve_disk,
            self.mass_disk,
            self.eccentricity,
            num_bodies_disk,
            false,
            v_scaler,
        ));

        // println!("Bodies: {:.4?}", &result);

        // result.append(&mut self.make_distributation(num_bodies_bulge, num_rings_bulge));

        println!("\nMaking bulge bodies...");
        if num_bodies_bulge > 0 && !self.mass_density_bulge.is_empty() {
            result.append(&mut make_distrib_data_area(
                &self.mass_density_bulge,
                &self.rotation_curve_bulge,
                self.mass_bulge,
                self.eccentricity,
                num_bodies_bulge,
                true,
                v_scaler,
            ));
        }

        result
    }
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

/// This (older, for us) approach interpolates using a radius scale decoupled from the data. It generates
/// a set of rings, and finds mass density and velocity by interpolating the data at those radii. Bodies
/// are only placed on these rings.
///
/// Angular positions are randomized.
fn make_distrib_along_rings(
    mass_density: &[(f64, f64)],
    vel: &[(f64, f64)],
    mass_total: f64,
    eccentricity: f64,
    num_bodies: usize,
    three_d: bool,
    v_scaler: f64,
) -> Vec<Body> {
    let mut result = Vec::with_capacity(num_bodies);
    let mut rng = rand::thread_rng();

    let num_rings = num_bodies / DISK_RING_PORTION;

    // Note: We are choosing the r range here based on mass density; not velocity.
    let r_last = mass_density.last().unwrap().0;
    let r_all = linspace(mass_density[0].0, r_last, num_rings);
    let dr = r_all[1] - r_all[0];

    // let total_area = r_last.powi(2) * TAU / 2.;
    // let total_volume = volume_sphere(r_last);

    let(bodies_by_r, mass_per_body_by_r) = select_num_bodies(&r_all, dr, mass_density, num_bodies);

    // Todo: If it works for your other approach, add in the center single body.

    for (i, r) in r_all.iter().enumerate() {
        println!(
            "Body data. r: {r} N bodies: {:?} mass-per-body: {:.4?}, mass-this-r: {:.4?}",
            bodies_by_r[i], mass_per_body_by_r[i], interpolate(mass_density, *r).unwrap()
        );

        let v_mag = interpolate(vel, *r).unwrap() * v_scaler;

        for _ in 0..bodies_by_r[i] {
            result.push(create_body(*r, mass_per_body_by_r[i], v_mag, eccentricity, three_d, &mut rng));
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

/// Select the number of bodies to create in a given radius-driven region.
fn select_num_bodies(r_all: &[f64], dr: f64, mass_density: &[(f64, f64)], num_bodies: usize) -> (Vec<usize>, Vec<f64>) {
    // todo: Experiment with this.
    let min_bodies_per_ring_area = 0.4;

    let num_r = r_all.len();

    let mut body_num_by_r_init = Vec::with_capacity(num_r);
    let mut mass_per_body_by_r = Vec::with_capacity(num_r);
    let mut num_bodies_init = 0.;

    // todo: If the core is too unstable, consider grouping all inner rings into a single body at the center.

    for r in r_all {
        // Compute the minimum number of bodies at each r
        let body_count_min = r.powi(2)* TAU / 2. * min_bodies_per_ring_area;

        body_num_by_r_init.push(body_count_min);

        // Set mass proportionally to initial body numbers. This line multiplies mass/area x area.
        let mass_this_area = interpolate(mass_density, *r).unwrap() * ring_area(*r, dr);
        mass_per_body_by_r.push(mass_this_area / body_count_min);

        num_bodies_init += body_count_min;
    }

    println!("Body num init: {:.4?}", body_num_by_r_init);

    let body_n_ratio = num_bodies as f64 / num_bodies_init;

    let mut bodies_by_r = Vec::with_capacity(num_r);

    for (i, r) in r_all.iter().enumerate() {
        // Scale body count to match target.
        bodies_by_r.push((body_num_by_r_init[i] / body_n_ratio) as usize);
    }

    (bodies_by_r, mass_per_body_by_r)

}

/// Add a body at a given distance from the center,  and random non-distance posit.
fn create_body(r: f64, mass: f64, v_mag: f64, eccentricity: f64, three_d: bool, rng: &mut ThreadRng) -> Body {
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
        let z = r * ϕ.cos();

        // todo: Update this logic to make the mass density sufficiently 3d?

        let scale_x = 1.0 - eccentricity; // Eccentricity factor for x-axis
        // let posit = Vec3::new(x * scale_x, y, z);
        // todo: Put eccentricity back.
        let posit = Vec3::new(x, y, z);

        // Velocity direction: perpendicular to the radius vector

        // Rotate along an arbitrary axis perpendicular to the center.
        let rot_axis = posit.to_normalized();
        let rotator = Quaternion::from_axis_angle(rot_axis, rng.gen_range(0.0..TAU));

        // Create a unit vector perpendicular to the position.
        let starting_pt = Vec3::new(1., 0., 0.);
        let perp_vec = starting_pt.cross(rot_axis);
        let vel = rotator.rotate_vec(perp_vec).to_normalized() * v_mag;

        // todo: Put eccentricity in velocity back.

        (posit, vel)
    } else {
        let x = r * θ.cos();
        let y = r * θ.sin();

        // todo: Temp mesasure for a finite-thick disk.
        let disk_thickness_div2 = 0.2;
        let z = rng.gen_range(-disk_thickness_div2..disk_thickness_div2);

        let scale_x = 1.0 - eccentricity; // Eccentricity factor for x-axis
        let posit = Vec3::new(x * scale_x, y, z);

        // Velocity direction: perpendicular to the radius vector
        let v_x = -v_mag * θ.sin(); // Tangential velocity in x-direction
        let v_y = v_mag * θ.cos(); // Tangential velocity in y-direction

        let vel = Vec3::new(v_x, v_y, 0.0);

        (posit, vel)
    };

    Body {
        posit,
        vel,
        accel: Vec3::new_zero(),
        mass,
    }
}

/// This (newer, for us) approach  maps out an area for each data piece, and fills it with bodies at random
/// positions. Position, both angular, and distance-within-ring, are randomized.
fn make_distrib_data_area(
    mass_density: &[(f64, f64)],
    vel: &[(f64, f64)],
    mass_total: f64,
    eccentricity: f64,
    num_bodies: usize,
    three_d: bool,
    v_scaler: f64,
) -> Vec<Body> {
    let mut result = Vec::with_capacity(num_bodies);
    let mut rng = rand::thread_rng();

    // let r_all: Vec<f64> = mass_density.iter().map(|(r, _mass)| r).collect();
    // let dr = r_all[1] - r_all[0];

    // let total_area = r_last.powi(2) * TAU / 2.;
    // let total_volume = volume_sphere(r_last);

    // let (bodies_by_r, mass_per_body_by_r) = select_num_bodies(&r_all, dr, mass_density, num_bodies);

    // We will model the inner most area by a single body, to minimize chaotic effects.
    // let rings_in_center = 30; // todo Crude metric that depends on the data. Starting point.
    let rings_in_center = 1; // todo Crude metric that depends on the data. Starting point.

    // todo: This is fuzzy, since we don't have a total-mass-in-region estimate.

    {
        let center_section_r = mass_density[rings_in_center].0;
        let area_center = center_section_r.powi(2) * TAU / 2.;
        // let volume_center = volume_sphere(center_section_r);


        let mass_density_center: f64 = mass_density[..rings_in_center].iter().map(|m| m.1).sum();
        // todo: This averages the mass density of the inner rings. Is this what we want?
        let mass_center =  mass_density_center / rings_in_center as f64 * area_center;

        // todo: Temp rm.

        // result.push(Body {
        //     posit: Vec3::new_zero(),
        //     vel: Vec3::new_zero(),
        //     accel: Vec3::new_zero(),
        //     mass: mass_center,
        // });
    }

    // Create bands of masses centered on each r.
    for (i, (r, mass)) in mass_density[0..mass_density.len() - 1].iter().enumerate() {
        if i < rings_in_center { // instead of indexing 1.., this keeps i in sync.
            continue
        }
        let r_prev = mass_density[i - 1].0;
        let r_this = mass_density[i].0;
        let r_next = mass_density[i + 1].0;

        let dr_prev = r_this - r_prev;
        let dr_next = r_next - r_this;

        let r_inner = r_this - dr_prev/2.;
        let r_outer = r_this + dr_next/2.;

        let area = {
            let area_outer = r_outer.powi(2) * TAU / 2.;
            let area_inner = r_inner.powi(2) * TAU / 2.;

            area_outer - area_inner
        };

        // todo: Handle the outer edge case too.
        // Set mass proportionally to initial body numbers. This line multiplies mass/area x area.
        let mass_this_area = mass * area;

        // todo temp: Even distribution.
        let body_num_this_area = num_bodies / (mass_density.len() - rings_in_center);

        let mass_per_body = mass_this_area / body_num_this_area as f64;

        println!(
            "Body data. r: {r} N bodies: {:?} mass-per-body: {:.0?}k, mass-this-r: {:.4?}",
            body_num_this_area, mass_per_body/1000., mass_this_area
        );

        for _ in 0..body_num_this_area {
            let r_body = rng.gen_range(r_inner..r_outer);
            let v_mag = interpolate(vel, r_body).unwrap() * v_scaler;

            result.push(create_body(r_body, mass_per_body, v_mag, eccentricity, three_d, &mut rng));
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
