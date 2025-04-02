# Robust-ImageSCC  
**Robust Mean Signal Estimation and Inference for Imaging Data**

## Overview  
Robust-ImageSCC implements a robust nonparametric framework for analyzing imaging data. The methods—centered on bivariate penalized splines over triangulation and a smoothed M-estimator—address the challenges of noisy and contaminated imaging data by enabling accurate mean signal estimation and robust simultaneous confidence corridors (SCC). Detailed theoretical properties and simulation studies demonstrating the estimator’s consistency, asymptotic normality, and robustness are provided in the accompanying manuscript.

## Files  
- **`robust_bpst_fda_functions.R`**  
  Contains the helper functions for robust mean function estimation via `mean_func_m_est()` and robust SCC construction via `robust_imagescc()`, along with routines for simulating contaminated or clean imaging data.

- **`cpp_func_cls.cpp`**  
  Implements the M-estimation routine using conjugate gradient descent to efficiently update proximal solutions, integrated with R through Rcpp.

- **`m-bpst-est-sim.R`**  
  Reproduces the simulation results for robust mean function estimation as described in the paper.

- **`m-bpst-est-sim.slurm`**  
  A SLURM scheduler script for running the mean function estimation simulations. It is configured to utilize job arrays (via the `--array` option) to efficiently manage multiple simulation replications.

- **`m-bpst-scc-sim.R`**  
  Reproduces the simulation results for robust SCC construction.

- **`m-bpst-scc-sim.slurm`**  
  A SLURM scheduler script for running the SCC simulations, also configured to support job arrays.

- **`load_cpp_cache.R`**  
  Provides a helper function to load the precompiled C++ code on environments where a suitable compiler is not available.

## Usage  
The simulation studies are designed to run on cluster systems using SLURM. The provided SLURM scripts use the `--array` option to submit multiple simulation jobs efficiently. For example:

- **OPTIONAL: precompiled C++ code on the head node**
  
  Launch R on the head node and run `load_cpp_cache.R`.

- **Mean Function Estimation Simulation (replace `<num_jobs>` with the desired number of array tasks.):**  
  Submit the job array by running:  
  ```{bash}
    sbatch --array=1-<num_jobs> m-bpst-est-sim.slurm
  ```
  
- **Robust SCC Simulation:**  
  Submit the job array by running:  
  ```{bash}
    sbatch --array=1-<num_jobs> m-bpst-scc-sim.slurm
  ```
