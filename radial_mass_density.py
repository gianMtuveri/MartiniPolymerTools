from MDAnalysis import mda
import numpy as np
import matplotlib.pyplot as plt


# Load the simulation
u = mda.Universe("path/to/file/topology.tpr",
                 "path/to/file/trajectory.xtc")

# Select atoms from nanoparticle/micelle (example: PEO-PLA)
peo_atoms = u.select_atoms("resname PEO")
pla_atoms = u.select_atoms("resname PLA")
wat_atoms = u.select_atoms("resname PW")
all_atoms = u.select_atoms("resname PLA PEO")  # used to define system COM

# Parameters
n_bins = 500
r_max = 200.0  # Å, adjust based on system size
bin_edges = np.linspace(0, r_max, n_bins + 1)
bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

peo_counts = np.zeros(n_bins)
pla_counts = np.zeros(n_bins)
wat_counts = np.zeros(n_bins)

# Residue masses in uma
peo_mass = peo_atoms.masses
pla_mass = pla_atoms.masses
wat_mass = wat_atoms.masses

# Loop over frames
for ts in u.trajectory[-2000:]:
    com = all_atoms.center_of_mass()  # or use hydrophobic core group if needed

    # Compute distances from COM
    peo_radii = np.linalg.norm(peo_atoms.positions - com, axis=1)
    pla_radii = np.linalg.norm(pla_atoms.positions - com, axis=1)
    wat_radii = np.linalg.norm(wat_atoms.positions - com, axis=1)

    # Bin them
    peo_hist, _ = np.histogram(peo_radii, bins=bin_edges) #, weights=peo_mass)
    pla_hist, _ = np.histogram(pla_radii, bins=bin_edges) #, weights=pla_mass)
    wat_hist, _ = np.histogram(wat_radii, bins=bin_edges) #, weights=pla_mass)

    peo_counts += peo_hist
    pla_counts += pla_hist
    wat_counts += wat_hist

# Normalize by number of frames
n_frames = len(u.trajectory)
peo_counts /= n_frames
pla_counts /= n_frames
wat_counts /= n_frames

# Compute shell volumes for normalization
shell_volumes = (4/3) * np.pi * (bin_edges[1:]**3 - bin_edges[:-1]**3)  # Å³

# Compute number density
peo_density = peo_counts / shell_volumes  # number/Å³
pla_density = pla_counts / shell_volumes
wat_density = wat_counts / shell_volumes

# Set a clean style
plt.style.use("default")  # modern seaborn-compatible theme
plt.rcParams.update({
    "figure.figsize": (16, 9),
    "font.size": 16,
    "axes.titlesize": 30,
    "axes.labelsize": 33,
    "legend.fontsize": 30,
    "xtick.labelsize": 28,
    "ytick.labelsize": 28,
    "lines.linewidth": 2,
    "lines.markersize": 6,
    "axes.spines.top": False,
    "axes.spines.right": False,
})

# multiply peo_density and pla_density by 1.66054 to convert amu/Å³ → g/cm³.

# Plot with labels and colors
plt.plot(bin_centers, peo_density, label='PEG', color='#b7b7ea', linewidth=4) 
plt.plot(bin_centers, pla_density, label='PLA', color='darkorange', linewidth=4) 
plt.plot(bin_centers, wat_density, label='PW', color='black', linewidth=4) 

# Axis labels and title
plt.xlabel("Distance from Center of Mass (Å)")
plt.ylabel("Density ($\AA$⁻³)")
#plt.title("Radial Density Profile of PEO and PLA")

# axis limit
plt.xlim(0,60)

# Optional: add grid, legend, and tight layout
plt.legend(loc='upper left', frameon=False)
plt.grid(True, linestyle='--', linewidth=0.5, alpha=0.7)
plt.tight_layout()

# Save and show
plt.savefig("pull/radial_density_profile_styled.png", dpi=300)

# Stack data column-wise: radius | PEO | PLA | water
data = np.column_stack((bin_centers, peo_density, pla_density, wat_density))

# Save to file
np.savetxt("pull/density_profiles.txt", data,
           header="Distance(nm)    PEO_density(g/cm^3)    PLA_density(g/cm^3)    PW_density(g/cm^3)",
           fmt="%.5f")
plt.show()

#############  Radius of Gyration  #############
# Select your atom group (e.g., full polymersome)
selection = u.select_atoms("resname PEO PLA")

# Get total number of frames
n_total = len(u.trajectory)
start_frame = max(0, n_total - 1000)

# Store Rg values
rg_vals = []

# Loop over the last 1000 frames
for i, ts in enumerate(u.trajectory[start_frame:]):
    rg = selection.radius_of_gyration()
    rg_vals.append(rg)

# Convert to numpy array
rg_vals = np.array(rg_vals)

# Plot histogram
plt.figure(figsize=(16, 9))
plt.hist(rg_vals, bins=20, color="#2ca02c", edgecolor="white", alpha=0.9)
plt.xlabel("Radius of Gyration (Å)")
plt.ylabel("Frequency")
plt.title("Histogram of Radius of Gyration (Last 1000 Frames)")
plt.grid(True, linestyle="--", alpha=0.5)
plt.tight_layout()
plt.savefig("radiusGyr_hist.png", dpi=300)

# Save the Rg values array
np.savetxt("rg_values.txt", rg_vals, header="Radius of Gyration (Å)", fmt="%.5f")
plt.show()
