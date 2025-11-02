import numpy as np
import matplotlib.pyplot as plt
import os
from ase import io

def avg_bead_spread_traj(m_beads, npart, cycle,rangeini,rangefin,intv):
    """
    Compute average bead spread for O, H1, H2 across PIMC trajectory.

    Parameters
    ----------
    m_beads : int
        Number of path integral beads.
    npart : int
        Number of molecules.
    cycle : int
        Step multiplier for plotting.
    """
    cwd = os.getcwd()
    coord_OHH = np.zeros(shape=(3*npart,3,m_beads),dtype=np.float64)
    bead_spread_O = []
    bead_spread_H1 = []
    bead_spread_H2 = []
    x = []
    
    for i in range (rangeini,rangefin,intv):
        print('Reading traj no - ',i)
        x.append(i)
        
        for bead in range(m_beads):
            filename = os.path.join(cwd, f"mc-traj_bead_{bead:02d}.xyz")        #reads the traj files kept in the folder containing this code
            # Read the trajectory
            traj = io.read(filename, index=i)
            coord_OHH[:,:,bead] = traj.positions

        dist_O = np.zeros(shape=(npart,m_beads),dtype=np.float64)
        dist_H1 = np.zeros(shape=(npart,m_beads),dtype=np.float64)
        dist_H2 = np.zeros(shape=(npart,m_beads),dtype=np.float64)
        for bead in range(m_beads):
            left = bead
            right = bead+1
            if(right == m_beads):
                right = 0
                # For all oxygen atoms
            dist_O[:, bead] = np.sqrt(np.sum((coord_OHH[0::3, :, left] - coord_OHH[0::3, :, right])**2, axis=1))
            
            # For the first hydrogens
            dist_H1[:, bead] = np.sqrt(np.sum((coord_OHH[1::3, :, left] - coord_OHH[1::3, :, right])**2, axis=1))
            
            # For the second hydrogens
            dist_H2[:, bead] = np.sqrt(np.sum((coord_OHH[2::3, :, left] - coord_OHH[2::3, :, right])**2, axis=1))

        
        bead_spread_O.append(np.mean(dist_O,axis=0))
        bead_spread_H1.append(np.mean(dist_H1,axis=0))
        bead_spread_H2.append(np.mean(dist_H2,axis=0))
        
        
    bead_spread_O = np.array(bead_spread_O)
    bead_spread_H1 = np.array(bead_spread_H1)
    bead_spread_H2 = np.array(bead_spread_H2)
        
    # X-axis values
    x = np.array(x)

    # Function to plot data
    def plot_distances(x, distances, title, ylabel):
        colors = ['r', 'g', 'b', 'c'] # Colors for the beads
        plt.figure(figsize=(8, 6))
        for bead in range(distances.shape[1]): # Iterate over beads
            plt.plot(x*cycle, distances[:, bead], label=f'Bead {bead + 1}')
        plt.title(title)
        plt.xlabel("MC Steps")
        plt.ylabel(ylabel)
        plt.legend()
        plt.grid(True)
        plt.show()

    # Plot for dist_O
    plot_distances(x, bead_spread_O, "Avg bead spread O", "Distance (Angstrom)")

    # Plot for dist_hyd1
    plot_distances(x, bead_spread_H1, "Avg bead spread H1", "Distance (Angstrom)")

    # Plot for dist_hyd2
    plot_distances(x, bead_spread_H2, "Avg bead spread H2", "Distance (Angstrom)")
            
        
