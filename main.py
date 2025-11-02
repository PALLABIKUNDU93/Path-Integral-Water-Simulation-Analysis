"""
Water Simulation Analysis Toolkit
=================================

This module provides post-processing routines for
path-integral Monte Carlo (PIMC) simulations of water.

Features:
- Property histogram comparison (LSI, ASY, etc.)
- RDF computation (OO, OH, HH)
- Density and compressibility estimation
- Bead spread analysis for PIMC convergence

Author: Pallabi Kundu
Affiliation: Stockholm University
"""

import os
from prop_analysis import prop_compare_from_npz_with_segments
from rdf_analysis import (
    generate_rdf_master_npz, plot_rdf_segments_from_npz
)
from density_analysis import (
    extract_density_npz, compute_avg_density_and_box, compute_compressibility_from_npz
)
from bead_spread import avg_bead_spread_traj


# Global parameters
m_beads = 4         # number of beads
npart = 125         # number of water molecules
cycle = 500         # after how manu cycles the trajectory is printed


if __name__ == "__main__":
    """
    Example execution order:
    1. Generate and compare histograms
    2. Compute RDFs
    3. Extract densities and compressibility
    4. Analyze bead spread
    """
# ----- Property comparison ----- 
#prop_compare_from_npz_with_segments(npz_file='prop_master.npz', total_frames=100, segment_size=25)
### total_frames -> total no of snapshots you want to calculate the property histogram from
### segment_size -> how you want to break up total_frames to identify any trends. 
### Example -> total_frames=100, segment_size=25 means you calculate hist from last 100 frames and also break that up into 4 subblocks of 25 frames each.


# ----- RDF generation and plotting ----- 
#generate_rdf_master_npz(m_beads=m_beads, nbins=200)
#plot_rdf_segments_from_npz(npz_file='rdf_master.npz', last_n_frames=100, segment_size=25)
### Example -> last_n_frames=100, segment_size=25 means you calculate rdf from last 100 frames and also break that up into 4 subblocks of 25 frames each.


# ----- Density analysis ----- 
### Extract and save density data
#extract_density_npz("mc-traj_bead_00.xyz", "density_data.npz", npart)

# ----- Compute averages ----- 
#stats = compute_avg_density_and_box("density_data.npz")
#print(stats)

# ----- Compute compressibility and plot convergence ----- 
#compute_compressibility_from_npz("density_data.npz", temperature=298, n_intervals=cycle)


# ----- Bead spread analysis ----- 
###rangeini and rangefin are the initial and final points of your trajectory file where you want the bead spread analysis
### intv is the increment (or interval) between successive values
#avg_bead_spread_traj(m_beads=m_beads, npart=npart, cycle=cycle, rangeini=0, rangefin=100, intv = 5)
