# MartiniPolymerTools

# Radial Charge Density & Zeta Potential Analysis Scripts

This repository provides Python scripts for analyzing and visualizing radial charge density profiles and zeta potential from molecular dynamics (MD) trajectories using [MDAnalysis](https://www.mdanalysis.org/). The scripts are designed for systems such as nanoparticles with polymers and ions, and are intended to help researchers extract physical insights from simulation data.

---

## Files

### 1. `radial_charge_density.py`

**Purpose:**  
Calculates and visualizes the absolute radial charge density around a nanoparticle or similar core structure, separated by charged species.

**Key Features:**
- Loads a trajectory and topology using MDAnalysis.
- Selects nanoparticle core and charged species.
- Computes the radial charge density histogram for each species.
- Outputs a bar chart of the radial (absolute) charge density profile.

**Usage:**
- Edit the file paths (`/path/to/topology`, `/path/to/trajectory`) to your own data.
- Adjust `charged_species` to include all relevant ions/charged groups in your system.
- Run the script:  
  ```bash
  python radial_charge_density.py
  ```
- The output plot (`charge_histograms.png`) will be saved in the working directory.

---

### 2. `pot_rho_profile.py`

**Purpose:**  
Calculates and visualizes both the net radial charge density and the electrostatic potential profile (ϕ(r)), including the zeta potential at a specified slipping plane.

**Key Features:**
- Loads and processes the MD trajectory.
- Computes the net charge density as a function of distance from the core.
- Solves Poisson’s equation in spherical symmetry to obtain ϕ(r).
- Plots both ρ(r) and ϕ(r) on the same graph with dual y-axes.
- Calculates and prints the zeta potential at a user-specified radius.

## Equations Used

The script numerically solves the Poisson equation in spherical symmetry to obtain the electrostatic potential profile from the charge density:

The Poisson equation in spherical symmetry is:
$
\frac{1}{r^2} \frac{d}{dr} \left( r^2 \frac{d\phi}{dr} \right) = -\frac{\rho(r)}{\epsilon_0}
$

where  
- $\phi(r)$ is the electrostatic potential,  
- $\rho(r)$ is the charge density,  
- $\epsilon_0$ is the vacuum permittivity.

The radial electric field:
$
E(r) = -\frac{d\phi}{dr}
$

Numerical solution (integration steps):

First integration (charge to field):
$
E(r) = \frac{1}{\epsilon_0 r^2} \int_0^r \rho(r') r'^2 dr'
$

Second integration (field to potential):
$
\phi(r) = -\int_0^r E(r') dr'
$

The zeta potential at the slip plane:
$
\zeta = \phi(r_\mathrm{slip})
$

**Usage:**
- Edit the file paths as in the previous script.
- Adjust `slip_radius_nm` to the desired location of the slipping plane (default: 6 nm).
- Run the script:  
  ```bash
  python pot_rho_profile.py
  ```
- The output plot (`pot_rho.png`) will be saved in the specified directory.
- The script will print the zeta potential at the selected slipping plane radius.

---

## Requirements

- Python 3.7+
- [MDAnalysis](https://www.mdanalysis.org/)
- [NumPy](https://numpy.org/)
- [Matplotlib](https://matplotlib.org/)
- [SciPy](https://scipy.org/)
- [tqdm](https://tqdm.github.io/)

Install dependencies (if needed) with:
```bash
pip install mdanalysis numpy matplotlib scipy tqdm
```

---

## Customization & Notes

- **File Paths:**  
  Update the paths to your topology and trajectory files in both scripts.
- **Atom Selections:**  
  Adjust the `select_atoms` queries to match your system's residue/atom names.
- **Bins & Ranges:**  
  Change `n_bins` and `r_max_nm` if your system's size differs.

---

## Citation

If you use these scripts in your work, please cite [MDAnalysis](https://www.mdanalysis.org/pages/citations/).

---

## License

MIT License
