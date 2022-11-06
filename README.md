# MOSAICS

MOSAICS is a collection of tools for characterizing membrane structure and dynamics within simulated trajectories. These tools are programmed in C++ and use the MOSAICS Analysis Template (MosAT) for reading trajectory data. 

The source code is organized as follows: First, the source code for MosAT-based tools is located in "src/MosAT/." Likewise, the source code for non-MosAT-based tools is located in" src/other_tools/," and header files are located in "src/headers/." 

Installation instructions and documentation for using the tools can be found in mosaics_user_manual.pdf. Note that we provide a sample trajectory in the "examples/membrane_thickness/" folder so the user can reproduce much of the analysis shown in the user manual. Examples of the input files, such as selection cards and parameter files, are also included in the "examples/" folder. We also include scripts, mainly used for ploting data, in the "scripts/" folder.

MOSAICS is licensed under the BSD-3-clause license, and additional details may be found in the "LICENSE" file. 

Please cite Bernhardt N, Faraldo-Gómez JD, MOSAICS: A Software Suite for Analysis of Membrane Structure and Dynamics in Simulated Trajectories, Biophysical Journal (2022), doi: https:// doi.org/10.1016/j.bpj.2022.11.005.

