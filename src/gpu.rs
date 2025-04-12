//! GPU computation, via CUDA.

use std::{sync::Arc, time::Instant};

use cudarc::driver::{CudaStream, CudaModule, LaunchConfig};
use lin_alg::f64::{alloc_vec3s, Vec3};

// The floating point type used on the device. This makes switching between f32 and f64 easier.
type FDev = f32;


/// Run coulomb attraction via the GPU. Computes per-sample potentials in paralle on the GPU; runs
/// the per-charge logic serial in the same kernel. This prevents needing to compute the sum on the CPU
/// afterwards. Returns a potential-per-sample Vec. Same API as the parallel+CPU approach above.
pub fn run_newton(
    stream: &Arc<CudaStream>,
    module: &Arc<CudaModule>,
    posits_src: &[Vec3],
    posits_tgt: &[Vec3],
    charges: &[f64], // Corresponds 1:1 with `posit_charges`.
) -> Vec<f64> {
    let start = Instant::now();

    // allocate buffers
    let n_sources = posits_src.len();
    let n_targets = posits_tgt.len();

    let posit_charges_gpus = alloc_vec3s(stream, posits_src);
    let posits_sample_gpu = alloc_vec3s(stream, posits_tgt);

    // Note: This step is not required when using f64ss.
    let charges: Vec<f32> = charges.iter().map(|c| *c as f32).collect();

    let mut charges_gpu = stream.alloc_zeros::<f32>(n_sources).unwrap();
    stream.memcpy_htod(&charges, &mut charges_gpu).unwrap();

    let mut V_per_sample = stream.alloc_zeros::<f32>(n_targets).unwrap();

    // let func_coulomb = module.load_function("coulomb_kernel").unwrap();
    let func_lj_V = module.load_function("lj_V_kernel").unwrap();
    // let func_lj_force = module.load_function("lj_force_kernel").unwrap();

    // let cfg = LaunchConfig::for_num_elems(n_targets as u32);

    let cfg = {
        const NUM_THREADS: u32 = 1024;
        let num_blocks = (n_targets as u32).div_ceil(NUM_THREADS);

        // Custom launch config for 2-dimensional data (?)
        LaunchConfig {
            grid_dim: (num_blocks, 1, 1),
            block_dim: (NUM_THREADS, 1, 1),
            shared_mem_bytes: 0,
        }
    };

    let mut launch_args = stream.launch_builder(&func_lj_V);

    launch_args.arg(&mut V_per_sample);
    launch_args.arg(&posit_charges_gpus);
    launch_args.arg(&posits_sample_gpu);
    launch_args.arg(&charges_gpu);
    launch_args.arg(n_sources,);
    launch_args.arg(n_targets);

    unsafe { launch_args.launch(cfg)}.unwrap();

    let result = stream.memcpy_dtov(&V_per_sample).unwrap();

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
