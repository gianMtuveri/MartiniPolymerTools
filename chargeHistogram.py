import numpy as np
import matplotlib.pyplot as plt
from MDAnalysis import Universe
from scipy.constants import elementary_charge as e_charge
from tqdm import tqdm
import matplotlib.colors as mcolors
import matplotlib.cm as cm

# === Load universe ===
print("Loading trajectory ...")
u = Universe("/path/to/topology", "/path/to/trajectory")
print("Universe loaded")
np_core = u.select_atoms("resname PLA PLN")       # nanoparticle core
charged_species = ["PLN", "NA+", "CL-"]           # Add all your charged atom names here

# === Parameters ===
n_bins = 300
r_max_nm = 30.0
r_max = r_max_nm * 1e-9  # m
r_edges = np.linspace(0, r_max, n_bins + 1)
r_centers = 0.5 * (r_edges[:-1] + r_edges[1:])
shell_volumes = (4/3) * np.pi * (r_edges[1:]**3 - r_edges[:-1]**3)

# === Prepare data containers ===
hist_data = {name: np.zeros(n_bins) for name in charged_species}

# === Loop over frames ===
for ts in tqdm(u.trajectory[-2500:], desc="Processing frames"):
    center = np_core.center_of_mass() * 1e-10  # Å → m

    for name in charged_species:
        ions = u.select_atoms(f"name {name}")
        positions = ions.positions * 1e-10  # Å → m
        charges = ions.charges * e_charge   # e → C

        r = np.linalg.norm(positions - center, axis=1)
        shell_q, _ = np.histogram(r, bins=r_edges, weights=np.abs(charges))
        hist_data[name] += shell_q

# === Normalize (optional) ===
n_frames = len(u.trajectory[-2500:])
for name in hist_data:
    hist_data[name] /= n_frames  # average over time
    hist_data[name] /= shell_volumes  # to get charge density (C/m³)

# === Plotting ===
fig, ax1 = plt.subplots(figsize=(30, 18))
plt.rcParams.update({
    "axes.labelweight": "bold",
    "axes.labelsize": 45,
    "xtick.labelsize": 60,
    "ytick.labelsize": 60,
    "axes.titlesize": 50,
    #"axes.titleweight": "bold"
})

colors = cm.get_cmap("tab10").colors
width = (r_edges[1] - r_edges[0]) * 1e9  # bin width in nm

for i, name in enumerate(charged_species):
    color = colors[i % len(colors)]
    plt.bar(
        r_centers * 1e9,
        hist_data[name] * 1e3,  # Convert C/m³ to mC/m³
        width=width,
        align='center',
        edgecolor=color,
        facecolor=mcolors.to_rgba(color, alpha=0.15),  # light transparent fill
        linewidth=2,
        label=name
    )

plt.xlabel("Radius (nm)", fontsize=60)
plt.ylabel("Absolute Charge Density (mC/m³)", fontsize=60)
plt.title("Radial Charge Density", fontsize=65)
plt.grid(True)
plt.legend(fontsize=50)
plt.tight_layout()
plt.savefig("charge_histograms.png", dpi=300)
plt.show()
