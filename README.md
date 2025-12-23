## About MOSAICS

MOSAICS is a collection of tools for characterizing membrane **structure** and **dynamics** within simulated trajectories of molecular systems. These tools focus on the efficient computation of **spatially resolved time-averaged** lipid observables like the membrane thickness, lipid tilt angle, residence time, etc. However, other types of analysis are possible; see "mosaics_user_manual.pdf" for details. To ensure peak performance, we have programmed MOSAICS in C++ and implemented a trajectory reading framework called the MOSAICS Analysis Template, or MosAT for short. MosAT reads **individual** trajectory snapshots into memory, thus reducing memory requirements and enabling the analysis of **very large** trajectory files. MosAT is also programmed to distribute the workload across available computing cores. As such, each MOSAICS tool is **fully parallelized**.    

## Software/hardware requirements
- A computing cluster
- A recent version of Open MPI
- A c++ compiler like g++
- A c++ compiler with MPI support like mpic++ 

## Installation instructions
Here we provide a brief overview of the installation process. We note that a more detailed description is provided in the user manual.
```
 cd MOSAICS
 mkdir build
 export MPI_CPP_COMP="mpic++"
 export CPP_COMP="g++"   
 sh install_commands
```

Alternatively, the user can try the following CMake recipe:
```
 cd MOSAICS
 mkdir build
 cmake src/MosAT/programs/ -S . -B build/
 cmake --build build/
```

## Getting started
Detailed instructions for using each tool can be found in the user manual, i.e., "mosaics_user_manual.pdf". Note that we provide a sample trajectory in the "examples/membrane_thickness/" folder so the user can reproduce much of the analysis shown in the user manual. Example input files, such as selection cards and parameter files, are also included in the "examples/" folder. We also include scripts, mainly used for ploting grid data, in the "scripts/" folder.

## Example Workflow
An example workflow is available in the examples/mosaics_work_flow.zip file. This item contains input files for performing a variety of MOSAICS analyses, including preparing the trajectory with Traj Prep and more. Detailed input files are included for each case, and the user is encouraged to reference them as needed. A helper script (do_analysis.sh) is also included to automate the workflow, enabling a full analysis with minimal effort. See mosaics_workflow.pdf after unpackaging mosaics_work_flow.zip for a more detailed description of the analysis covered.

## Code layout
The MOSAICS source code is organized as follows: 
- Source code for MosAT-based tools is located in "src/MosAT/" 
- Source code for non-MosAT-based tools is located in "src/other_tools/" 
- Header files encompassing the MOSAICS library are located in "src/headers/" 
- Files sourced from MDTraj that are used for reading GROMACS trajectories are located in "src/xdr/"
- Routines sourced from the GROMACS MD simulation package are located in "src/gmx_lib/"   

## License
MOSAICS is licensed under the BSD-3-clause license, and additional details may be found in the "LICENSE" file; see also for licensing info regarding code sourced from GROMACS and MDTraj.  

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

## Feedback
Please report bugs identified within the source code and other errors like typos found in the user manual using the [issues](https://github.com/MOSAICS-NIH/MOSAICS/issues) feature. Other comments and suggestions are always welcome. Thank you for your help. 

## Troubleshooting
Are you having difficulty using MOSAICS? Not certain how to prepare your simulation for analysis? Need help adding an unrecognized lipid type or other basic inputs parameters? No problem! Please post your questions in the [issues](https://github.com/MOSAICS-NIH/MOSAICS/issues) section. 
