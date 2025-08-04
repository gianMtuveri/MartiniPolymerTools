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
$\frac{1}{r^2} \frac{d}{dr} \left( r^2 \frac{d\phi}{dr} \right) = -\frac{\rho(r)}{\epsilon_0}$

where  
- $\phi(r)$ is the electrostatic potential,  
- $\rho(r)$ is the charge density,  
- $\epsilon_0$ is the vacuum permittivity.

The radial electric field:
$E(r) = -\frac{d\phi}{dr}$

Numerical solution (integration steps):

First integration (charge to field):
$E(r) = \frac{1}{\epsilon_0 r^2} \int_0^r \rho(r') r'^2 dr'$

Second integration (field to potential):
$\phi(r) = -\int_0^r E(r') dr'$

The zeta potential at the slip plane:
$\zeta = \phi(r_\mathrm{slip})$

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




# Radial Density and Potential Analysis Scripts

This repository provides Python scripts for analyzing radial mass/charge densities and electrostatic potential profiles from molecular dynamics (MD) simulation trajectories using [MDAnalysis](https://www.mdanalysis.org/). These tools help you extract and visualize the local structure and electrostatics around nanoparticles, micelles, or similar aggregates.

---

## Files

### 1. `radial_mass_density.py`

**Purpose:**  
Calculates and visualizes the **radial mass density profiles** for different species (PEO, PLA, water) around the center of mass of a nanoparticle or micelle, and computes the radius of gyration distribution.

**Main Features:**
- Loads a trajectory and topology file via MDAnalysis.
- Selects PEO, PLA, and water atoms (modify for your system).
- Computes and averages the mass density profile in spherical shells over the last 2000 frames.
- Plots the mass density profiles as a function of distance from the aggregate center.
- Saves numerical data of the profiles to `pull/mass_density_profiles.txt`.
- Computes and saves a histogram of the radius of gyration (Rg) for the last 1000 frames.

**Usage:**
1. Edit the file paths and atom/residue names as appropriate for your system.
2. Run:
   ```bash
   python radial_mass_density.py
   ```
3. Output files:
   - `pull/radial_mass_density_profile.png`: Plot of radial mass densities.
   - `pull/mass_density_profiles.txt`: Numerical density profile data.
   - `radiusGyr_hist.png`: Histogram plot of Rg.
   - `rg_values.txt`: Numerical values of Rg.

---

### 2. `charge_potential.py`

**Purpose:**  
Calculates the **radial net charge density** and the **electrostatic potential profile** (ϕ(r)) from charge distributions, and computes the zeta potential at a specified slip plane.

**Main Features:**
- Loads simulation data and selects system core and charged/ion species.
- Computes net charge density profile as a function of radius.
- Numerically solves the spherically symmetric Poisson equation to obtain the electrostatic potential profile.
- Plots both the radial charge density and potential on dual y-axes.
- Prints the zeta potential at a user-specified radius.
- Saves profiles and plots.

## Equations Used

The script numerically solves the Poisson equation in spherical symmetry to obtain the electrostatic potential profile from the charge density:

The Poisson equation in spherical symmetry is:
$\frac{1}{r^2} \frac{d}{dr} \left( r^2 \frac{d\phi}{dr} \right) = -\frac{\rho(r)}{\epsilon_0}$

where  
- $\phi(r)$ is the electrostatic potential,  
- $\rho(r)$ is the charge density,  
- $\epsilon_0$ is the vacuum permittivity.

The radial electric field:
$E(r) = -\frac{d\phi}{dr}$

Numerical solution (integration steps):

First integration (charge to field):
$E(r) = \frac{1}{\epsilon_0 r^2} \int_0^r \rho(r') r'^2 dr'$

Second integration (field to potential):
$\phi(r) = -\int_0^r E(r') dr'$

The zeta potential at the slip plane:
$\zeta = \phi(r_\mathrm{slip})$

**Usage:**
1. Edit file paths and atom selections for your system.
2. Run:
   ```bash
   python charge_potential.py
   ```
3. Output files:
   - `pot_rho.png`: Combined charge density and potential plot.
   - Numerical results for further analysis.

---

## Requirements

- Python 3.7+
- [MDAnalysis](https://www.mdanalysis.org/)
- [NumPy](https://numpy.org/)
- [Matplotlib](https://matplotlib.org/)
- [SciPy](https://scipy.org/)
- [tqdm](https://tqdm.github.io/)

Install dependencies with:
```bash
pip install mdanalysis numpy matplotlib scipy tqdm
```

---

## Customization & Notes

- **File Paths:**  
  Set correct paths to your topology and trajectory files in both scripts.
- **Atom Selections:**  
  Adjust `select_atoms` queries to match your system's residue/atom names.
- **Bins & Ranges:**  
  Modify `n_bins`, `r_max`, etc. for your system size or resolution.

---

## Citation

If you use these scripts, please cite [MDAnalysis](https://www.mdanalysis.org/pages/citations/).

---

## License

MIT License
