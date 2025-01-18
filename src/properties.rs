//! Get the properties of collections of bodies, e.g. galaxies

// todo: Output types A/R. Currently x/y tuples.

// todo: You're mixing kpc (mass) with km/s (rotation velocity)

const N_SAMPLE_PTS: usize = 40;

use lin_alg::f64::Vec3;
use plotters::{
    element::PathElement,
    prelude::{BitMapBackend, ChartBuilder, Color, IntoDrawingArea, BLACK, BLUE, WHITE},
    series::LineSeries,
};

use crate::Body;
use crate::units::KPC_MYR_PER_KM_S;

fn get_nearby_pts(bodies: &[Body], center: Vec3, r: f64, dr: f64) -> Vec<&Body> {
    // Todo: Consider a fuzzy, weighted dropoff instead of these hard boundaries. Or not;
    // todo: Maybe this is fine.
    bodies
        .iter()
        .filter(|b| ((b.posit - center).magnitude() - r).abs() <= dr / 2.)
        // .map(|b2| b2.posit)
        .collect()
}

/// Gravitational potential (𝚽). X: r (kpc) Y: 𝚽 (J/kg)
pub fn gravity_potential(bodies: &[Body], center: Vec3, r_max: f64) -> Vec<(f64, f64)> {
    Vec::new()
}

fn find_r_max(bodies: &[Body], center: Vec3) -> f64 {
    let mut result = 0.;
    for body in bodies {
        let r = (body.posit - center).magnitude();
        if r > result {
            result = r;
        }
    }
    result
}

/// Normalized mass density. X: r (kpc). Y: ρ/ρ_0
pub fn mass_density(bodies: &[Body], center: Vec3) -> Vec<(f64, f64)> {
    let mut result = Vec::with_capacity(N_SAMPLE_PTS);

    let r_max = find_r_max(bodies, center);

    let dr = r_max / N_SAMPLE_PTS as f64;

    // todo: Fix this. I believe it depends on the model. Determines rho_0.
    let r0_i = 0;
    let mut rho_0 = 1.;

    for i in 0..N_SAMPLE_PTS {
        let r = i as f64 * dr;

        let nearby_masses: Vec<f64> = get_nearby_pts(bodies, center, r, dr)
            .into_iter()
            .map(|b2| b2.mass)
            .collect();

        if nearby_masses.is_empty() {
            result.push((r, 0.));
        } else {
            let mut mass = 0.;
            for nearby_mass in &nearby_masses {
                mass += nearby_mass;
            }

            if i == r0_i {
                rho_0 = mass
            }

            // todo: Density; not pure mass. Although this may be fine for now as it's normalized?
            result.push((r, mass));
        }
    }

    // Normalize.
    for r in &mut result {
        r.1 /= rho_0
    }

    result
}

/// Luminosity profile. X: r (kpc). Y: μ (mag arcsec^-2) - Surface brightness profile.
pub fn luminosity(bodies: &[Body]) -> Vec<(f64, f64)> {
    let mut result = Vec::with_capacity(N_SAMPLE_PTS);

    result
}

/// Normalized rotation curve. X: r (kpc). Y: V/c, or km/s, or kpc/MLY?
/// We specify r_max, to avoid calculations involving outliers. But, perhaps should calculate anyway.
/// todo: In km/s for now, not V/C.
pub fn rotation_curve(bodies: &[Body], center: Vec3, c: f64) -> Vec<(f64, f64)> {
    let mut result = Vec::with_capacity(N_SAMPLE_PTS);

    let r_max = find_r_max(bodies, center);

    let dr = r_max / N_SAMPLE_PTS as f64;

    for i in 0..N_SAMPLE_PTS {
        let r = i as f64 * dr;

        let nearby_pts: Vec<Vec3> = get_nearby_pts(&bodies, center, r, dr)
            .into_iter()
            .map(|b2| b2.vel)
            .collect();

        if nearby_pts.is_empty() {
            result.push((r, 0.));
        } else {
            let mut v = 0.;
            for pt in &nearby_pts {
                v += pt.magnitude();

                // println!("RCD: {:?}", self.rotation_curve_disk[10])
                // println!("V Mag: {:?} scaled: ", pt.magnitude(), pt.magnitude());

            }

            // Convert to km/s
            // result.push((r, v / (KPC_MYR_PER_KM_S * nearby_pts.len() as f64 * c)));
            result.push((r, v / (KPC_MYR_PER_KM_S * nearby_pts.len() as f64)));
        }
    }

    result
}

/// Sersic index. X: α. Y: s.
pub fn sersic(bodies: &[Body]) -> Vec<(f64, f64)> {
    let mut result = Vec::with_capacity(N_SAMPLE_PTS);

    result
}

/// Display a 2d plot of properties, e.g. rotation curve, luminosity etc.
pub fn plot(data: &[(f64, f64)], x_label: &str, y_label: &str, plot_title: &str, filename: &str) {
    // Find the x and y ranges using PartialOrd
    let x_range = data
        .iter()
        .map(|(x, _)| *x)
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(min, max), x| {
            (min.min(x), max.max(x))
        });

    // println!("Data: {:?}", data);
    // println!("X range: {:?}", x_range);

    // let x_range = (0., 10.); // todo temp

    let y_range = data
        .iter()
        .map(|(_, y)| *y)
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(min, max), y| {
            (min.min(y), max.max(y))
        });

    // Create a drawing area for the plot
    let fname = format!("{filename}.png");
    let root = BitMapBackend::new(&fname, (800, 600)).into_drawing_area();
    root.fill(&WHITE).unwrap();

    // Create a chart builder
    let mut chart = ChartBuilder::on(&root)
        .caption(plot_title, ("sans-serif", 20))
        .margin(10)
        .x_label_area_size(30)
        .y_label_area_size(30)
        .build_cartesian_2d(x_range.0..x_range.1, y_range.0..y_range.1)
        .unwrap();

    // Set labels
    chart
        .configure_mesh()
        .x_desc(x_label)
        .y_desc(y_label)
        .draw()
        .unwrap();

    // Plot the data points
    chart
        .draw_series(LineSeries::new(data.iter().cloned(), BLUE)) // Use `.cloned()` here
        .unwrap()
        .label("Data")
        .legend(|(x, y)| PathElement::new([(x, y), (x + 20, y)], BLUE));

    // Draw the legend
    chart
        .configure_series_labels()
        .background_style(WHITE.mix(0.8))
        .border_style(BLACK)
        .draw()
        .unwrap();
}

pub fn plot_rotation_curve(data: &[(f64, f64)], desc: &str) {
    plot(
        data,
        "r (kpc)",
        // "v / c",
        "km/s",
        // &format!("Normalized rotation curve of {desc}"),
        &format!("Rotation curve of {desc}"),
        &format!("rot_plot_{desc}"),
    );
}

pub fn plot_mass_density(data: &[(f64, f64)], desc: &str) {
    plot(
        data,
        "r (kpc)",
        "ρ / ρ₀",
        &format!("Normalized mass density of {desc}"),
        &format!("mass_plot_{desc}"),
    );
}
