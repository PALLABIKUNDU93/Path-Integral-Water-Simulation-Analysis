import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib import cm
from ase.io import read
from ase.geometry.analysis import Analysis

def generate_rdf_master_npz(m_beads, nbins):
    """Compute bead-resolved O–O RDF and save to npz."""
    
    elem="O O"
    
    rdf_all_beads = []
    r_vals = None

    for bead in range(m_beads):
        filename = f"mc-traj_bead_{bead:02d}.xyz"
        traj = read(filename, index=':')
        analysis = Analysis(traj)
        print(f"Bead {bead}: {len(traj)} frames")
        analysis = Analysis(traj)
        bead_rdfs = []

        for i, frame in enumerate(traj):
            box_lengths = frame.cell.lengths()
            rmax = np.floor(np.min(box_lengths) * 0.5)

            gr = analysis.get_rdf(rmax=rmax, nbins=nbins, imageIdx=i, elements=elem.split())[0]
            if r_vals is None:
                r_vals = np.linspace(0, rmax, nbins)
            bead_rdfs.append(gr)

        bead_rdfs = np.array(bead_rdfs)  # shape: (n_frames, nbins)
        rdf_all_beads.append(bead_rdfs)

    rdf_all_beads = np.array(rdf_all_beads)  # shape: (m_beads, n_frames, nbins)

    print(f"RDF shape: {rdf_all_beads.shape}")

    np.savez('rdf_master.npz', rdf=rdf_all_beads, r=r_vals)
    print("Saved RDF master to rdf_master.npz")
    

def plot_rdf_segments_from_npz(npz_file, last_n_frames, segment_size):
    """Plot RDF with segment-wise convergence."""
    data = np.load(npz_file)
    rdf_all = data['rdf']  # shape: (m_beads, n_frames, nbins)
    r = data['r']

    if rdf_all.ndim != 3:
        raise ValueError(f"Expected RDF shape (m_beads, n_frames, nbins), got {rdf_all.shape}")

    m_beads, n_frames, nbins = rdf_all.shape
    print(f"Loaded RDF: {m_beads} beads, {n_frames} frames, {nbins} bins")

    if last_n_frames > n_frames:
        raise ValueError(f"last_n_frames ({last_n_frames}) exceeds total frames ({n_frames})")

    # Extract last frames
    rdf_recent = rdf_all[:, -last_n_frames:, :]  # (m_beads, last_n_frames, nbins)
    rdf_avg = np.mean(rdf_recent, axis=(0, 1))   # average over beads and frames

    # Plot
    plt.figure()
    plt.title("RDF OO")

    # Main line
    plt.plot(r, rdf_avg, label='Avg RDF (target)', color='black', linewidth=2.5)

    # Add reference RDF curves
    ref_paths = {
        #'298K Experimental': 'specify path',
        #'PIMD Paesani': 'specify path',
        #'RDF without target': 'specify path'
    }

    for label, path in ref_paths.items():
        try:
            ref = np.loadtxt(path, usecols=(0, 1))
            plt.plot(ref[:, 0], ref[:, 1], label=label)
        except Exception as e:
            print(f"Could not load reference '{label}': {e}")

    # Segment plots
    num_segments = last_n_frames // segment_size
    colors = cm.viridis(np.linspace(0.2, 0.95, num_segments))

    for i in range(num_segments):
        start = i * segment_size
        end = start + segment_size
        seg = rdf_recent[:, start:end, :]  # (m_beads, segment_size, nbins)
        seg_avg = np.mean(seg, axis=(0, 1))  # avg over beads and frames in segment
        plt.plot(r, seg_avg, linestyle='--', color=colors[i], label=f"Seg {i + 1}")

    plt.xlabel("r (Å)")
    plt.ylabel("g(r)")
    plt.grid(True)
    plt.legend(fontsize='x-small')
    plt.tight_layout()
    plt.show()
    
    # Save RDF data to file
    script_dir = os.path.dirname(os.path.abspath(__file__)) # Get the current script's directory
    subfolder_path = os.path.join(script_dir, 'prop_with_target') # Define the subfolder path
    os.makedirs(subfolder_path, exist_ok=True) # Create the subfolder if it doesn't exist
    fname = os.path.join(subfolder_path, 'rdf_without_target.dat') # Construct the full file path
    print(f"Saving RDF to: {fname}")
    np.savetxt(fname, np.column_stack((r, rdf_avg)), fmt="%.6f")

