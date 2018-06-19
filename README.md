# QPI_project
"qpi-disp-interp.f90" is used to calculate QPI intensity.  The code is parallelized using MPI.  Compile using mpif90 or equivalent compiler.  The results are published in https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.119.107003

Before running QPI calculation, you must have obtained self-energies.  Self-energies can be obtained by running "spectral_function.m" found in "Self-energy" folder.  Examples can be found in http://iopscience.iop.org/article/10.1088/0953-2048/29/5/054009/meta
