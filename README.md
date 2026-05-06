# Minimal action shortcut to adiabaticity in a driven Kitaev chain: competing gaps in a topological transition at finite-time

This repository contains data generation scripts and plotting routines for the figures presented in "[Minimal action shortcut to adiabaticity in a driven Kitaev chain: competing gaps in a topological transition at finite-time
]([COLOQUE_O_LINK_DO_ARXIV_AQUI](https://arxiv.org/abs/2604.26146))". The data was generated using the Python library `qutip` and the Julia library `QuantumOptics`, and plotted using the Python's `Matplotlib` and Julia's `CairoMakie`.

### 🔹 Energy Gaps & Levels
* **`Fig1_c/`**: Notebooks and scripts for calculating Kitaev chain energy levels and identifying the competing gaps across the topological transition.

### 🔹 Fidelity Dynamics (MA-STA vs Linear Ramp)
* **`Fig3/`**: Bulk fidelity simulations for different system sizes ($N=20, 50, 80$) and long-time limits.
* **`Fig4_a_b/`**: Comparative analysis of fidelity and chemical potential $\mu(t)$ for MA, 2nMA, and Linear protocols.
* **`Fig5_a_b/`**: Comparison between control over first and second momentum sectors for even and odd ground states.
* **`Fig6_a_b/`**: Comparison between bulk approximations and full numerical evolution for even and odd ground states.

### 🔹 Work Distribution & Statistics
* **`Fig7_a_b/`**: Data and scripts for the work probability distribution $P(W)$. Includes datasets for both Linear Ramps (LR) and MA-STA across different $\tau$ values.
* **`Fig7_c_d_e/`**: Calculation and plotting of the moments of work (mean, variance, skewness) for $N=14$.

### 🔹 Robustness & Appendices
* **`Fig8_a_b/`**: Study of the protocol robustness against variations in the pairing potential $\Delta$, including 3D distribution plots.

### 🔹 Core Simulation Codes
* The **`codes/`** folder contains the backbone of the project:
    * `Fidel_MA_Kitaev_even_odd_gs.ipynb`, `Fidelity_n2MA_Kitaev_Qutip.ipynb`: QuTiP-based fidelity calculation notebooks with bulk approximation.
    * `masta_kitaev_deg_gs_evol.jl`: QuantumOptics-based fidelity calculation script with exact numerics.
    * `calc_eigenstate.jl`: Script for calculating full eigenspectrum for multiples values of $\tau$. Necessary for calculating the work distribution.
    * `calculate_moments_exact_diagonalization.jl`: Script for calculating arrays with evolved states for moments and work distribution.
    * `load_evol_calculate_moments_exact_diagonalization.jl`: Script for loading arrays and calculating exact moments.
    * `num_work_distro_load.jl`: Script for loading arrays and calculating the work distribution data.

---

## Installation & Requirements

**Python 3.12.3** was used, along with the following packages:
* `qutip`: v4.7.6
* `numpy`: v1.26.4
* `Matplotlib`: v3.6.3

**Julia 1.12.5** was used, along with the following libraries:
* `QuantumOptics.jl`: v1.2.4
* `CairoMakie.jl`: v0.15.8
* `GeometryBasics.jl`: v0.5.10
* `JLD2.jl`: v0.6.4
* `ProgressMeter.jl`: v1.11.0
* `LaTeXStrings.jl`: v1.4.0

---

## Citation

Feel free to use any codes or data made available. If you do, please cite us :)

> Rafael Bentes de Sales and Krissia Zawadzki. *[Minimal action shortcut to adiabaticity in a driven kitaev chain: competing gaps in a topological transition at finite-time](https://arxiv.org/abs/2604.26146)*, 2026.
