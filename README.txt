# Minimal action shortcut to adiabaticity in a driven Kitaev chain: competing gaps in a topological transition at finite-time

This repository contains data generation scripts and plotting routines for the figures presented in "[Minimal action shortcut to adiabaticity in a driven Kitaev chain: competing gaps in a topological transition at finite-time
]([COLOQUE_O_LINK_DO_ARXIV_AQUI](https://arxiv.org/abs/2604.26146))". The data was generated using the Python library `qutip` and the Julia library `QuantumOptics`, and plotted using the Python's `Matplotlib` and Julia's `CairoMakie`.

## Repository Structure

The code is organized by figure to allow for easy reproduction of the results found in the paper. Each folder contains the specific simulation parameters and plotting scripts used for that dataset.

### 🔹 Single-Qubit Dynamics
* **`Fig1_c_d/`**: Overlap of a single-qubit state with its Liouvillian left eigenvectors, before and after the application of a Controlled-Ry gate.

### 🔹 Two-Qubit Markovian Case
* **`Fig2_a_b/`**: Overlap of a two-qubit state with its Liouvillian left eigenvectors, before and after the application of a Controlled-Ry gate.
* **`Fig2_c/`**: Asymptotic speedup as a function of the $T_1$ and $T_2$ relaxation times.
* **`Fig2_d/`**: Trace distance between the ground state and 1000 two-qubit states with $q_1$ initialized as a random Haar states and $q_2$ in its excited state, along with the speedup distributions.

### 🔹 Two-Qubit Non-Markovian Case
* **`Fig4_a/`**: Spectra of the Markovian embedding Liouvillian $\mathcal{L}_{emb}$ and the reduced Liouvillian $\mathcal{L}_{red}$ as a function of time.
* **`Fig4_b/`**: Non-Markovian case: Trace distance between the ground state and 1000 two-qubit states with $q_1$ initialized as a random Haar states and $q_2$ in its excited state, along with the speedup distributions.
* **`Fig5_a_b/`**: Comparison of the two-qubit state purity and coherence time evolutions between the Markovian embedding and the reduced model, for different values of the qubit/TLS frequency.
* **`Fig6/`**: Robustness of the C-Ry–based reset protocol under imperfect control, in the Markovian and non-Markovian case.

### 🔹 Experimental Data
* **`Fig7/`**: Experimental demonstration of the Mpemba-enhanced qubit reset.

### 🔹 Appendices
* **Figure 8**: Can be reproduced using the `Fig2_d/` code by changing the ancilla state to `qt.basis(2,0)`.
* **`Fig9/`**: Speedup distributions for different ancilla excited state populations and coherences.

---

## 🛠 Installation & Requirements

To run these simulations, **Python 3.12** was used, along with the following packages:
* `scipy`: v1.14.1
* `qutip`: v5.1.0
* `numpy`: v1.26.4
* `h5py`: v3.15.1

To plot the data, **Julia 1.10** was used, along with the following libraries:
* `CairoMakie.jl`: v0.15.6
* `ColorSchemes.jl`: v3.31.0
* `GeometryBasics.jl`: v0.5.10
* `HDF5.jl`: v0.17.2
* `KernelDensity.jl`: v0.6.11
* `LaTeXStrings.jl`: v1.4.0

> **💡 Tip:**
> Most of the core Python functions used for simulations are stored in the `sim_functions.py` file.

---

## 📝 Citation

If you use this code in your research, please cite:

> Théo Lejeune, Miha Papič, John Goold, Felix C. Binder, François Damanet, Mattia Moroder, *"[Accelerating qubit reset through the Mpemba effect](COLOQUE_O_LINK_DO_ARXIV_AQUI)"*.
