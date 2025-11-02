import numpy as np
import matplotlib.pyplot as plt
from ase.io import read


def extract_density_npz(traj_file="mc-traj_bead_00.xyz",
                        npz_out="density_data.npz",
                        n_molecules=125):
    """
    Extract densities and box lengths from an ASE trajectory.

    Parameters
    ----------
    traj_file : str, optional
        Path to the trajectory file (default: 'mc-traj_bead_00.xyz').
    npz_out : str, optional
        Output filename for the saved npz data.
    n_molecules : int, optional
        Number of water molecules in the system.

    Saves
    -----
    density_data.npz : npz file containing
        - box_lengths : ndarray
        - densities : ndarray
    """
    # Physical constants
    mass_water_gmol = 18.01528       # g/mol
    NA = 6.02214076e23               # Avogadro number
    mass_water_g = mass_water_gmol / NA  # g per molecule

    # Load trajectory
    traj = read(traj_file, index=":")
    print(f"Read {len(traj)} frames from {traj_file}")

    # Compute densities and box lengths
    box_lengths, densities = [], []
    for atoms in traj:
        L = atoms.cell.lengths()[0]  # assume cubic box
        V_cm3 = (L ** 3) * 1e-24     # Å³ → cm³
        rho = n_molecules * mass_water_g / V_cm3
        box_lengths.append(L)
        densities.append(rho)

    np.savez(npz_out, box_lengths=np.array(box_lengths),
             densities=np.array(densities))
    print(f"Saved densities to {npz_out}")


def compute_avg_density_and_box(npz_file="density_data.npz", plot=True):
    """
    Compute average density, its uncertainty, and mean box length.

    Parameters
    ----------
    npz_file : str
        Input npz file with box_lengths and densities.
    plot : bool
        Whether to plot density vs MC step.

    Returns
    -------
    dict
        {'avg_density': float, 'err_density': float, 'avg_box': float}
    """
    data = np.load(npz_file)
    densities = data["densities"]
    box_lengths = data["box_lengths"]

    avg_density = np.mean(densities)
    err_density = np.std(densities) / np.sqrt(len(densities))
    avg_box = np.mean(box_lengths)

    print(f"Average density     : {avg_density:.4f} g/cm³ ± {err_density:.4f}")
    print(f"Average box length  : {avg_box:.4f} Å")

    if plot:
        plt.figure(figsize=(6, 4))
        plt.plot(densities, label="Density")
        plt.xlabel("MC Step")
        plt.ylabel("Density (g/cm³)")
        plt.title("Density vs MC Step")
        plt.grid(True, ls="--", alpha=0.6)
        plt.legend()
        plt.tight_layout()
        plt.show()

    return {
        "avg_density": avg_density,
        "err_density": err_density,
        "avg_box": avg_box
    }


def compute_compressibility_from_npz(npz_file,temperature,n_intervals):
    """
    Compute and visualize isothermal compressibility from volume fluctuations.

    Parameters
    ----------
    npz_file : str, optional
        Path to npz file containing 'box_lengths'.
    temperature : float, optional
        Simulation temperature (K).
    n_intervals : int, optional
        Number of progressive intervals for convergence plotting.

    Notes
    -----
    - Computes kappa_T from volume fluctuations:
        kappa_T = (⟨V²⟩ - ⟨V⟩²) / (kB T ⟨V⟩)
    - Plots convergence of kappa_T with MC step.
    """
    kB = 1.380649e-23  # J/K
    atm_conv = 101325  # Pa → atm

    data = np.load(npz_file)
    L = data["box_lengths"]
    V = (L ** 3) * 1e-30  # Å³ → m³

    # Full data compressibility
    V_mean = np.mean(V)
    kappa_T = (np.mean(V ** 2) - V_mean ** 2) / (kB * temperature * V_mean) * atm_conv
    print(f"Isothermal compressibility: {kappa_T:.3e} (1/atm)")

    # Progressive (convergence) analysis
    step_size = len(V) // n_intervals
    indices = np.arange(step_size, len(V) + 1, step_size)
    kappa_T_values = []

    for end in indices:
        V_subset = V[:end]
        V_mean = np.mean(V_subset)
        fluct = np.mean(V_subset ** 2) - V_mean ** 2
        kappa_T_values.append(fluct / (kB * temperature * V_mean) * atm_conv)

    plt.figure(figsize=(7, 5))
    plt.plot(indices, kappa_T_values, "o-", color="darkred", lw=1.8, ms=6)
    plt.xlabel("MC Steps", fontsize=12)
    plt.ylabel(r"Isothermal Compressibility $\kappa_T$ (1/atm)", fontsize=12)
    plt.title("Convergence of $\kappa_T$", fontsize=13)
    plt.grid(True, ls="--", alpha=0.6)
    plt.tight_layout()
    plt.show()
