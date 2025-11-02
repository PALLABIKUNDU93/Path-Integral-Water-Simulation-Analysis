# Water Simulation Analysis Toolkit

**Author:** Pallabi Kundu  
**Affiliation:** Stockholm University

This repository provides a collection of Python scripts for post-processing
path-integral Monte Carlo simulations of water. The toolkit analyses
structural and thermodynamic properties such as:

- Property histograms (LSI, ASY, TET, ZET)
- RDFs
- Density and compressibility
- Bead-spread convergence diagnostics
(Add more as and when required...)

## Dependencies
- Python ≥ 3.8
- ASE
- NumPy
- Matplotlib
- SciPy

## Usage
```bash
python main.py

## Remember
main.py is the MAIN func!
Edit it to activate specific analyses.
Keep the traj file and the property .npz file in the same folder as your .py codes to avoid hassle.
the traj files are sample files, replace with your own traj files.
