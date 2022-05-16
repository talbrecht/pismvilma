#pismvilma
Coupling of ice sheet model PISM with solid Earth model VILMA


## coupling scripts

Submit the heterogeneous slurm batch job `sbatch slurm_vilma3d-pism.cmd`, with selected time parameters and SBATCH options. Default is 


- timeint=.1 # kyr, meaning 100 years
- ys=-246.0  # kyr before present
- ye=0       # kyr before present
- it=1       # number of iterations to decrease the misfit of preent-day bed topography

The slurm job calls the python coupler script `python run_coupler.py`, which prepares the input data and the calls the model steps `bash run_pism_step.sh` and `tcsh run_vilma3d_step.cli` and remaps the data from a global Gauss-Legendre grid to the projected PISM grid, using [CDO](https://code.mpimet.mpg.de/projects/cdo/) and [NCO](http://nco.sourceforge.net/) tools, and vise versa.
