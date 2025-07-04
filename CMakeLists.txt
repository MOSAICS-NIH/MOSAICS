
#try the following install commands
#cmake src/MosAT/programs/ -S . -B build/
#cmake --build build/

cmake_minimum_required(VERSION 3.4.3)

project(MOSAICS CXX)

#use c++11 
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#Set the base directory for locating header files
include_directories(src/)

#add library for reading xdr files
add_library(xdrfile_xtc src/xdr/src/xdrfile_xtc.cpp)
add_library(xdrfile_trr src/xdr/src/xdrfile_trr.cpp)
add_library(xdrfile     src/xdr/src/xdrfile.cpp)
add_library(xdr_seek    src/xdr/src/xdr_seek.cpp)

#Add targets for MosAT based programs
foreach(CMD
    2d_enrichment
    3d_enrichment
    2d_kinetics
    apl
    apl_3d
    atomic_density_3d
    atomic_rmsf
    atoms_in_2_planes
    atom_select
    b_stamp
    b_stamp_grid
    bilayer_z
    check_broken_mols
    contact_analysis
    contact_kinetics
    contact_rmsf
    dihedrals
    distances
    h_bond_kinetics
    interdigitation
    inter_leaflet_contacts
    inter_leaflet_contacts_3d
    lipid_contacts
    lipid_contacts_3d
    lipid_density
    lipid_density_3d
    lipid_distances
    lipid_distances_3d
    lipid_flip
    lipid_gyration
    lipid_h_bonds
    lipid_immobilization
    lipid_mixing
    lipid_msd
    lipid_orientation
    lipid_orientation_3d
    lipid_protein_min_dist
    lipid_salt_bridges
    mean_coords
    mean_lipid_coords
    mean_lipid_coords_3d
    mean_protein_coords
    membrane_thickness
    mosat
    nearest_neighbors
    nearest_neighbors_3d
    p2
    p2_3d
    pbc_gen
    pbc_xy
    pbc_z
    pdb_editor
    protein_contacts
    protein_lipid_contacts
    protein_orientation
    protein_residue_enrichment
    protein_translator
    surface_residue_finder
    symmetry_enforcer
    system_translator
    traj_prep
    traj_time
    zcoord
    )

  add_executable(${CMD}_mpi src/MosAT/programs/${CMD}.cpp)

  target_link_libraries(${CMD}_mpi PUBLIC xdrfile_xtc)
  target_link_libraries(${CMD}_mpi PUBLIC xdrfile_trr)
  target_link_libraries(${CMD}_mpi PUBLIC xdrfile)
  target_link_libraries(${CMD}_mpi PUBLIC xdr_seek)

  target_include_directories(${CMD}_mpi PUBLIC
                            "xdrfile_xtc"
                            )

  target_include_directories(${CMD}_mpi PUBLIC
                            "xdrfile_trr"
                            )

  target_include_directories(${CMD}_mpi PUBLIC
                            "xdrfile"
                            )

  target_include_directories(${CMD}_mpi PUBLIC
                            "xdr_seek"
                            )
endforeach()

#Add targers for other programs serial
foreach(CMD
    bilayer_free_energy     
    bonds_generator        
    data_averager           
    delta_plot              
    grid_addition           
    grid_data_excluder     
    grid_editor          
    grid_region_integrator  
    histogram              
    leaflet_averager        
    mask_maker              
    mean_coords_row_selector
    midplane                
    nan_selector            
    orientation_histogram  
    protein_mask            
    protein_mask_grower     
    )
  add_executable(${CMD} src/other_tools/${CMD}.cpp)
endforeach()

#Add targers for other programs mpi
foreach(CMD
    2d_enrichment_distance_projection
    2d_kinetics_percent_visited
    2d_kinetics_distance_projection_global
    2d_kinetics_distance_projection_window
    binding_contributors
    binding_events_analyzer
    binding_events_analyzer_single
    binding_events_merger
    binding_events_video
    binding_extractor
    binding_lipids
    binding_list
    binding_occupancy
    binding_time_line
    grid_distance_projection
    lipid_exchange
    lipid_exchange_distances
    single_frame_distributions
    single_frame_error
    solvation_shells
    )
  add_executable(${CMD}_mpi src/other_tools/${CMD}.cpp)
endforeach()

