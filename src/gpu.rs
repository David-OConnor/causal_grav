//! GPU computation, via CUDA.

use std::{sync::Arc, time::Instant};

use cudarc::driver::{CudaDevice, CudaSlice, LaunchAsync, LaunchConfig};
use lin_alg::f64::Vec3;

// The floating point type used on the device. This makes switching between f32 and f64 easier.
type FDev = f32;

/// Convert a collection of `Vec3`s into Cuda arrays of their components.
fn alloc_vec3s(dev: &Arc<CudaDevice>, data: &[Vec3]) -> CudaSlice<FDev> {
    let mut result = Vec::new();
    // todo: Ref etcs A/R; you are making a double copy here.
    for v in data {
        result.push(v.x as FDev);
        result.push(v.y as FDev);
        result.push(v.z as FDev);
    }
    dev.htod_copy(result).unwrap()
}

/// Run coulomb attraction via the GPU. Computes per-sample potentials in paralle on the GPU; runs
/// the per-charge logic serial in the same kernel. This prevents needing to compute the sum on the CPU
/// afterwards. Returns a potential-per-sample Vec. Same API as the parallel+CPU approach above.
pub fn run_newton(
    dev: &Arc<CudaDevice>,
    posits_charge: &[Vec3],
    posits_sample: &[Vec3],
    charges: &[f64], // Corresponds 1:1 with `posit_charges`.
) -> Vec<f64> {
    let start = Instant::now();

    // allocate buffers
    let n_charges = posits_charge.len();
    let n_samples = posits_sample.len();

    let posit_charges_gpus = alloc_vec3s(&dev, posits_charge);
    let posits_sample_gpu = alloc_vec3s(&dev, posits_sample);

    // Note: This step is not required when using f64ss.
    let charges: Vec<FDev> = charges.iter().map(|c| *c as FDev).collect();

    let mut charges_gpu = dev.alloc_zeros::<FDev>(n_charges).unwrap();
    dev.htod_sync_copy_into(&charges, &mut charges_gpu).unwrap();

    let mut V_per_sample = dev.alloc_zeros::<FDev>(n_samples).unwrap();

    let kernel = dev.get_func("cuda", "coulomb_kernel").unwrap();

    let cfg = LaunchConfig::for_num_elems(n_samples as u32);

    let cfg = {
        const NUM_THREADS: u32 = 1024;
        let num_blocks = (n_samples as u32 + NUM_THREADS - 1) / NUM_THREADS;

        // Custom launch config for 2-dimensional data (?)
        LaunchConfig {
            grid_dim: (num_blocks, 1, 1),
            block_dim: (NUM_THREADS, 1, 1),
            shared_mem_bytes: 0,
        }
    };

    unsafe {
        kernel.launch(
            cfg,
            (
                &mut V_per_sample,
                &posit_charges_gpus,
                &posits_sample_gpu,
                &charges_gpu,
                n_charges,
                n_samples,
            ),
        )
    }
    .unwrap();

    let result = dev.dtoh_sync_copy(&V_per_sample).unwrap();

    // Some profiling numbers for certain grid sizes.
    // 2D, f32: 99.144 ms
    // 3D, f32: 400.06 ms
    // 2D, f64: 1_658 ms
    // 3D, f64: 1_643 ms
    // 300 ms for both large and small sizes on f32 with std::sqrt???

    let time_diff = Instant::now() - start;
    println!("GPU coulomb data collected. Time: {:?}", time_diff);

    // This step is not required when using f64.
    result.iter().map(|v| *v as f64).collect()
    // result
}
