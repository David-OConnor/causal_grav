pub fn linspace(start: f64, stop: f64, num_points: usize) -> Vec<f64> {
    if num_points < 2 {
        return vec![start];
    }

    let step = (stop - start) / (num_points - 1) as f64;
    (0..num_points).map(|i| start + i as f64 * step).collect()
}

/// This function generates an interpolated value for the given `val` based on the
/// provided `data`. The `data` is a set of (x, y) pairs, where `x` is the input
/// and `y` is the corresponding output value.
pub fn interpolate(data: &[(f64, f64)], val: f64) -> Option<f64> {
    // Ensure there are at least two points for interpolation
    if data.len() < 2 {
        panic!("At least two data points are required for interpolation.");
    }

    // Find the interval in which `val` falls
    for i in 0..data.len() - 1 {
        let (x0, y0) = data[i];
        let (x1, y1) = data[i + 1];

        if val >= x0 && val <= x1 {
            // Perform linear interpolation
            let t = (val - x0) / (x1 - x0);
            return Some(y0 + t * (y1 - y0));
        }
    }

    // If `val` is outside the range of data, extrapolate
    if val < data[0].0 {
        let (x0, y0) = data[0];
        let (x1, y1) = data[1];
        return Some(y0 + (val - x0) / (x1 - x0) * (y1 - y0));
    } else if val > data[data.len() - 1].0 {
        let (x0, y0) = data[data.len() - 2];
        let (x1, y1) = data[data.len() - 1];
        return Some(y0 + (val - x0) / (x1 - x0) * (y1 - y0));
    }

    None
}
