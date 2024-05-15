
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
    atoms_in_2_planes
    atom_select
    b_stamp
    b_stamp_grid
    bilayer_z
    check_broken_mols
    contact_kinetics
    contact_rmsf
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
    lipid_gyration
    lipid_h_bonds
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
    lipid_flip
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
    delta_plot
    grid_data_excluder
    grid_addition
    grid_region_integrator
    leaflet_averager
    protein_mask
    grid_editor
    midplane
    mask_maker
    protein_mask_grower
    nan_selector
    orientation_histogram
    mean_coords_row_selector
    bilayer_free_energy
    histogram
    data_averager
    bonds_generator
    )
  add_executable(${CMD} src/other_tools/${CMD}.cpp)
endforeach()

#Add targers for other programs mpi
foreach(CMD
    grid_distance_projection
    single_frame_distributions
    2d_enrichment_distance_projection
    binding_list
    binding_events_analyzer_single
    lipid_exchange
    binding_time_line
    2d_kinetics_distance_projection_global
    2d_kinetics_distance_projection_window
    binding_events_merger
    binding_extractor
    binding_contributors
    binding_occupancy
    binding_events_analyzer
    binding_events_video
    solvation_shells
    2d_kinetics_percent_visited
    binding_lipids
    lipid_exchange_distances
    single_frame_error
    )
  add_executable(${CMD}_mpi src/other_tools/${CMD}.cpp)
endforeach()

