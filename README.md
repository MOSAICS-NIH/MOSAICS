## About MOSAICS

MOSAICS is a collection of tools for characterizing membrane structure and dynamics within simulated trajectories. These tools are programmed in C++ and use the MOSAICS Analysis Template (MosAT) for reading trajectory data. 

## Getting started
Installation instructions and documentation for using the tools can be found in the user manual, i.e., "mosaics_user_manual.pdf". Note that we provide a sample trajectory in the "examples/membrane_thickness/" folder so the user can reproduce much of the analysis shown in the user manual. Example input files, such as selection cards and parameter files, are also included in the "examples/" folder. We also include scripts, mainly used for ploting data, in the "scripts/" folder.

## Code layout
The MOSAICS source code is organized as follows: 
- Source code for MosAT-based tools is located in "src/MosAT/" 
- Source code for non-MosAT-based tools is located in "src/other_tools/" 
- Header files encompassing the MOSAICS library are located in "src/headers/" 
- Files sourced from MDTraj that are used for reading GROMACS trajectories are located in "src/xdr/"
- Routines sourced from the GROMACS MD simulation package are located in "src/gmx_lib/"   

## License
MOSAICS is licensed under the BSD-3-clause license, and additional details may be found in the "LICENSE" file; see also for licensing info for code sourced from GROMACS and MDTraj.  

## Citing MOSAICS
Please cite Bernhardt N, Faraldo-Gómez JD, MOSAICS: A Software Suite for Analysis of Membrane Structure and Dynamics in Simulated Trajectories, Biophysical Journal (2022), doi: https:// doi.org/10.1016/j.bpj.2022.11.005.

#### BibTeX
```
@article{BERNHARDT2022,
    title = {MOSAICS: A Software Suite for Analysis of Membrane Structure and Dynamics in Simulated Trajectories},
    journal = {Biophysical Journal},
    year = {2022},
    issn = {0006-3495},
    doi = {https://doi.org/10.1016/j.bpj.2022.11.005},
    url = {https://www.sciencedirect.com/science/article/pii/S0006349522009031},
    author = {Nathan Bernhardt and José D. Faraldo-Gómez},
}
```
