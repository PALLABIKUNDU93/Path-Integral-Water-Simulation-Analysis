import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib import cm

def prop_compare_from_npz_with_segments(npz_file, total_frames, segment_size):
    """
    Compare property histograms from simulation data.

    Parameters
    ----------
    npz_file : str
        Path to .npz file containing property arrays (LSI, ASY, etc.)
    total_frames : int
        Number of frames from the end to include in analysis
    segment_size : int
        Number of frames per segment for progressive histograms
    """
    props = ['LSI', 'ASY', 'G_5', 'TET', 'ZET']
    bins_dict = {
        'LSI': (51, (0.0, 0.47)),
        'ASY': (51, (0.0, 1.95)),
        'G_5': (51, (2.70, 4.46)),
        'TET': (51, (-0.52, 1.00)),
        'ZET': (51, (-1.12, 1.59))
    }

    data = np.load(npz_file)
    m_beads, n_frames, npart = data['LSI'].shape
    print(f"Loaded {npz_file}: {m_beads} beads, {n_frames} frames, {npart} particles")

    for prop in props:
        num_bins, range_vals = bins_dict[prop]
        bin_edges = np.linspace(*range_vals, num_bins + 1)
        bin_mid = 0.5 * (bin_edges[1:] + bin_edges[:-1])

        # Compute histogram
        full_data = data[prop][:, -total_frames:, :].reshape(-1)
        hist_total, _ = np.histogram(full_data, bins=bin_edges, density=True)

        os.makedirs('prop_with_target', exist_ok=True)
        fname = os.path.join('prop_with_target', f'{prop}_histogram.dat')
        np.savetxt(fname, np.column_stack((bin_mid, hist_total)), fmt="%.6f")
        print(f"Saved histogram: {fname}")

        # Plot
        plt.figure()
        plt.title(f"{prop} Histogram")
        plt.plot(bin_mid, hist_total, color='black', linewidth=2.0, label='Current Simulation')

        # Segment histograms
        num_segments = total_frames // segment_size
        colors = cm.viridis(np.linspace(0.2, 0.95, num_segments))
        for i in range(num_segments):
            s, e = (n_frames - total_frames) + i * segment_size, (n_frames - total_frames) + (i + 1) * segment_size
            seg_data = data[prop][:, s:e, :].reshape(-1)
            seg_hist, _ = np.histogram(seg_data, bins=bin_edges, density=True)
            plt.plot(bin_mid, seg_hist, linestyle='--', color=colors[i], alpha=0.8, label=f"Segment {i+1}")

        plt.xlabel(prop)
        plt.ylabel("Probability Density")
        plt.legend(fontsize='x-small')
        plt.grid(True)
        plt.tight_layout()
        plt.show()
