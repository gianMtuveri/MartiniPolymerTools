# === Load trajectory ===
print("Loding trajectory ...")
u = mda.Universe("path/to/file", "path/to/file")
print("Universe %s loaded"%ch)
# === Atom selections (adjust as needed) ===
np_core = u.select_atoms("resname PLA PLN")       # nanoparticle core
system = u.select_atoms("resname PEO PLA PLN Ion")  # nanoparticle + polymer + ions

# === Parameters ===
n_bins = 700
r_max_nm = 35.0  # max radius in nanometers
r_max = r_max_nm * 1e-9  # convert to meters
r_edges = np.linspace(0, r_max, n_bins + 1)
r_centers = 0.5 * (r_edges[:-1] + r_edges[1:])
shell_volumes = (4/3) * np.pi * (r_edges[1:]**3 - r_edges[:-1]**3)
slip_radius_nm = 5.0

# === Accumulate charge density and Z potential ===
zeta_list = []
phi_list = []
frames_used = []
rho_total = np.zeros(n_bins)
for ts in tqdm(u.trajectory[-1000:], desc="Processing frames"):
    center = np_core.center_of_mass() * 1e-10  # convert Å → m
    positions = system.positions * 1e-10       # Å → m
    charges = system.charges * e_charge       # e → C

    r = np.linalg.norm(positions - center, axis=1)
    shell_q, _ = np.histogram(r, bins=r_edges, weights=charges)
    rho_r = shell_q / shell_volumes  # C/m³
    rho_total += rho_r

    # Solve Poisson’s equation
    integrand = rho_r * r_centers**2
    E_r = cumtrapz(integrand, r_centers, initial=0) / (epsilon * r_centers**2)
    phi_r = -cumtrapz(E_r, r_centers, initial=0)
    phi_r -= np.mean(phi_r[int(0.9 * n_bins):])  # Normalize to bulk

    slip_idx = np.argmin(np.abs(r_centers - slip_radius_nm * 1e-9))
    zeta = phi_r[slip_idx] * 1e3  # V → mV

    zeta_list.append(zeta)
    phi_list.append(phi_r)
    frames_used.append(ts.frame)

# === Average over frames ===
rho_avg = rho_total / len(u.trajectory[-1000:])

# === Plot results ===
fig, ax1 = plt.subplots(figsize=(30, 18))

plt.rcParams.update({
    "axes.labelweight": "bold",
    "axes.labelsize": 45,
    "xtick.labelsize": 60,
    "ytick.labelsize": 60,
    "axes.titlesize": 50,
    #"axes.titleweight": "bold"
})


color1 = "tab:blue"
ax1.set_xlabel("Radius (nm)")
ax1.set_ylabel("Charge Density ρ(r) [mC/m³]", color=color1)
ax1.hist(r_centers * 1e9, bins=r_edges[::10]*1e9, weights=rho_avg * 1e3, edgecolor=color1, 
         histtype='step', facecolor='w', lw=6)
#1.plot(r_centers * 1e9, rho_avg * 1e3, label="ρ(r)", color=color1, lw=6)
ax1.tick_params(axis='y', labelcolor=color1)

# Add φ(r) on second axis
ax2 = ax1.twinx()
color2 = "tab:red"
ax2.set_ylabel("Electrostatic Potential φ(r) [mV]", color=color2)
ax2.plot(r_centers * 1e9, phi_r * 1e3, label="φ(r)", color=color2, lw=6)
ax2.tick_params(axis='y', labelcolor=color2)

fig.tight_layout()
plt.title("Radial Charge Density and Zeta Potential Profile")
plt.grid(True)
plt.savefig(""path/to/image/pot_rho.png", bbox_inches='tight')
plt.show()

# === Optional: Print zeta potential at ~slipping plane radius (e.g., 6 nm)
slip_radius_nm = 6
idx_slip = np.argmin(np.abs(r_centers - slip_radius_nm * 1e-9))
zeta_potential_mv = phi_r[idx_slip] * 1e3
