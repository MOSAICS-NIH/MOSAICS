#include <iostream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <fstream>
#include <vector>
#include <algorithm>
using namespace std;

#include "xdr/include/xdrfile_xtc.h"                           //Native XTC input/output
#include "xdr/include/xdr_seek.h"                              //Native trajectory seeking
#include "xdr/include/xdrfile_trr.h"                           //Native TRR input/output
#include "xdr/include/xdrfile.h"                               //Native XDR support
#include "xdr/include/trr_header.h"                            //Native TRR header support
#include "headers/multi_dim_vec.h"                             //Native explicit multidimensional vector types
#include "headers/switch.h"                                    //Native on/off switches
#include "headers/file_reader.h"                               //Native text-file routines
#include "headers/vector_mpi.h"                                //Native collection and summation of vector rows
#include "headers/mosat_routines.h"                            //Native system variables and reporting
#include "headers/file_naming.h"                               //Native output-name tagging
#include "headers/file_naming_mpi.h"                           //Native MPI extension checks
#include "headers/command_line_args_mpi.h"                     //Native MPI command-line interface
#include "MosAT/program_variables/pv_lipid_neighbor_content.h" //Analysis settings and help description
#include "headers/array.h"                                     //Native array support
#include "headers/performance.h"                               //Native performance bookkeeping
#include "headers/index.h"                                     //Native selection-card support
#include "headers/traj.h"                                      //Streaming trajectory and geometric centers
#include "headers/parallel.h"                                  //Alternative frame ownership with full trajectory access
#include "headers/leaflet_finder.h"                            //Native lipid selection and molecule traversal
#include "headers/protein_finder.h"                            //Scaffold finder definitions
#include "headers/sol_finder.h"                                //Scaffold finder definitions
#include "headers/force_serial.h"                              //Native parallel-capability check
#include "headers/param.h"                                     //Networked selection cards

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function filters pair contacts, builds lipid_comp/type_comp and stores the owned central frame.      //
//                                                                                                           //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void lipid_neighbor_content(Trajectory &traj,program_variables &p,Param &param,sv2d &center_atoms,
                            iv3d &filter,dv3d &type_comp_frames,iv2d &neighbors_stats)
{

    int i = 0;        //standard variable used in loops
    int j = 0;        //standard variable used in loops
    int k = 0;        //standard variable used in loops
    int l = 0;        //standard variable used in loops
    int m = 0;        //standard variable used in loops
    int n = 0;        //standard variable used in loops
    int o = 0;        //standard variable used in loops
    int q = 0;        //standard variable used in loops

    dv2d lipid_comp(p.num_lipids,dv1d(param.main_size_y(),0.0));   //Neighbor-type fractions for each selected central lipid
    iv1d lipid_type(p.num_lipids,0);                               //Type of each selected lipid, identified in this frame
    iv1d res_id(p.num_lipids,0);                                   //Store the residue id for each lipid
    int lipid_i_count = 0;                                         //Selected central lipid's row in lipid_comp
    int pos           = traj.current_frame%p.window;               //Circular-history position for the current observation

    for(i=0; i<traj.target_leaflet.size(); i++) //Outer target-leaflet loop: central lipid i
    {
        int min = traj.t_lip_start(i);                 //First atom of the central lipid
        int max = traj.t_lip_end(i);                   //Last atom of the central lipid
        i       = traj.next_target_lipid(i);

        for(j=0; j<param.main_size_y(); j++) //loop over lipid types i
        {
            if(traj.res_name[min] == param.param_main_s[j][0]) //lipid type is a match
            {
                dv1d center_i         = traj.center(center_atoms[j],min,max); //Central lipid's geometric center at the read frame
                double cutoff         = param.param_main_d[j][1];             //Central-type cutoff distance in nm
                int lipid_j_count     = 0;                                    //Selected neighbor's index in the pair history

                for(int k=0; k<traj.target_leaflet.size(); k++) //Inner target-leaflet loop: possible neighbor lipid j
                {
                    int min_j = traj.t_lip_start(k);   //First atom of the possible neighbor lipid
                    int max_j = traj.t_lip_end(k);     //Last atom of the possible neighbor lipid
                    k         = traj.next_target_lipid(k);

                    for(l=0; l<param.main_size_y(); l++) //loop over lipid types j
                    {
                        if(traj.res_name[min_j] == param.param_main_s[l][0])
                        {
                            lipid_type[lipid_j_count] = l;
                            res_id[lipid_j_count]     = traj.res_nr[min_j];

                            int n_ij = 0;             //Raw contact at the read frame
                            if(min_j != min)          //A lipid is never its own neighbor
                            {
                                dv1d center_j = traj.center(center_atoms[l],min_j,max_j); //Neighbor geometric center
                                for(m=-1; m<=1; m++) //loop over x-dimension
                                {
                                    for(n=-1; n<=1; n++) //loop over y-dimension
                                    {
                                        for(o=-1; o<=1; o++) //loop over x-dimension
                                        {
                                            for(q=-1; q<=1; q++) //loop over y-dimension
                                            {
                                                double dx   = (center_i[0] + m*traj.box[XX][XX]) - (center_j[0] + o*traj.box[XX][XX]);                  //Center displacement along x in nm
                                                double dy   = (center_i[1] + n*traj.box[YY][YY]) - (center_j[1] + q*traj.box[YY][YY]);                  //Center displacement along y in nm
                                                double dist = sqrt(dx*dx + dy*dy);                                                                      //Direct three-dimensional center distance
                                                if(dist < cutoff)
                                                {
                                                    n_ij = 1;
                                                }
                                            }
                                        }
                                    }
                                }
                            }

                            //add state to filter
                            filter[lipid_i_count][lipid_j_count][pos] = n_ij;

                            lipid_j_count++;
                            break;
                        }
                    }
                }
                lipid_i_count++;
                break;
            }
        }
    }

    //The inner loop is complete: divide lipid_i_count by the accumulated filtered-neighbor count.
    if((traj.current_frame-p.read_first) >= 2*p.range) //filter is full
    {
        for(i=0; i<p.num_lipids; i++) //loop over lipids i 
        {
            int neighbor_count    = 0;                                    //How many neighbors does lipid i have

            iv1d neighbors(p.num_lipids,0);                               //stores the status for each lipid j

            //filter data
            for(j=0; j<p.num_lipids; j++) //loop over lipids_j 
            {
                double hits = 0;
                for(k=0; k<(2*p.range + 1); k++) //loop over window
                {
                    hits += filter[i][j][k]; 
                }
                double freq = hits/(double)(2*p.range + 1);

                if(freq >= 0.5)
                {
                    neighbors[j] = 1; 
                }
                else
                {
                    neighbors[j] = 0;
                }
                neighbor_count += neighbors[j];
            }
            //update neighbor count stats
            neighbors_stats[lipid_type[i]].push_back(neighbor_count);           
 
            //sum over lipid type for j
            for(j=0; j<p.num_lipids; j++) //loop over lipids_j 
            {
                for(k=0; k<param.main_size_y(); k++) //loop over lipid types
                { 
                    if(lipid_type[j] == k)
                    {
                        lipid_comp[i][k] += neighbors[j];
                    }
                }
            }

            //normalize lipid_i
            int b_empty = 0;
            for(k=0; k<param.main_size_y(); k++) //loop over lipid types
            {
                if(neighbor_count > 0)
                {
                    lipid_comp[i][k] = lipid_comp[i][k]/(double)neighbor_count;
                }
                else 
                {
                    b_empty          = 1;
  		    lipid_comp[i][k] = 0;
                }
            }
            if(b_empty == 1)
            {
                p.empty_count++;
                //printf("Lipid %d res_id %d observed with no neighbors on frame %d! \n",i,res_id[i],traj.get_frame_global());
            }
        }

        //average over lipid type i 
        dv2d type_comp(param.main_size_y(),dv1d(param.main_size_y(),0.0)); //Mean neighbor composition for each central lipid type
        iv1d lipid_type_count(param.main_size_y(),0);                      //Number of central lipids contributing to each type's mean

        for(i=0; i<p.num_lipids; i++) //loop over lipids i 
        {
            for(j=0; j<param.main_size_y(); j++) //loop over lipid types 
            {
                type_comp[lipid_type[i]][j] += lipid_comp[i][j];
            }
            lipid_type_count[lipid_type[i]]++;
        }
        for(i=0; i<param.main_size_y(); i++) //loop over lipid types i
        {
            for(j=0; j<param.main_size_y(); j++) //loop over lipid types j 
            {
                type_comp[i][j] /= (double)lipid_type_count[i];
            }
        }

        //store current frame data for type comp
        for(i=0; i<param.main_size_y(); i++) //loop over lipid types i
        {
            for(j=0; j<param.main_size_y(); j++) //loop over lipid types j
            {
        	type_comp_frames[i][j][traj.current_frame-p.read_first-2*p.range] = type_comp[i][j];
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function collects type_comp_frames, averages global frame blocks and writes the three data files.    //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
double finalize_analysis(Trajectory &traj,system_variables &s,program_variables &p,Param &param,
                         dv3d &type_comp_frames,iv2d &neighbors_stats)
{
    int i = 0;                                     //standard variable used in loops
    int j = 0;                                     //standard variable used in loops
    int k = 0;                                     //standard variable used in loops

    clock_t t_i = clock();                         //Start time of collection and output
    MPI_Barrier(MPI_COMM_WORLD);

    if(s.world_rank == 0)
    {
        printf("\nCollecting frame data.\n");
    }

    //collect type_comp_frames
    for(i=0; i<param.main_size_y(); i++) //loop over lipid type i
    {
        for(j=0; j<param.main_size_y(); j++) //loop over lipid type j
        {
            collect_dv1d(s.world_size,s.world_rank,type_comp_frames[i][j]);
        }
    }

    //get stats for lipids with no neighbors
    int global_empty_count = 0;
    MPI_Reduce(&p.empty_count, &global_empty_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if(s.world_rank == 0)
    {
        int frames_analyzed = traj.get_num_frames()-2*p.range;
        double empty_frac = 100.0*(double)global_empty_count/(double)(frames_analyzed*p.num_lipids);
        printf("\nThere were %d of %d lipids detected with no neighbors (%f %s). \n",global_empty_count,frames_analyzed*p.num_lipids,empty_frac,"%");
        printf("This corresponds to about %f lipids with no neighbors per trajectory frame. \n",(double)global_empty_count/(double)frames_analyzed);
    }
    
    //compute neighbor stats (how many neighbors did each type have on average).
    for(i=0; i<param.main_size_y(); i++) //loop over lipid types
    {
        int global_neighbors_stats = 0;
        collect_iv1d(s.world_size,s.world_rank,neighbors_stats[i]);

        if(s.world_rank == 0)
        {
            double avg = 0.0;
            for(j=0; j<neighbors_stats[i].size(); j++) //loop over data
            {
                avg = avg + neighbors_stats[i][j];
            }
            avg = avg/(double)neighbors_stats[i].size();

            double stdev = 0.0;
            for(j=0; j<neighbors_stats[i].size(); j++) //loop over data
            {
                stdev = stdev + pow((avg-neighbors_stats[i][j]),2.0);
            }
            stdev = sqrt(stdev/(double)neighbors_stats[i].size());

            printf("Average neighbors (%s): %f +/- %f \n",param.param_main_s[i][0].c_str(),avg,stdev);
        } 
    }

    //report data
    if(s.world_rank == 0)
    {
        string frame_name = add_tag(p.lnc_file_name,"_frames");  //Per-frame output derived from the requested base
        string block_name = add_tag(p.lnc_file_name,"_blocks");  //Per-block output derived from the requested base
	string delta_name = add_tag(p.lnc_file_name,"_delta");   //Adjacent-block changes derived from the base

	FILE *frame_file  = fopen64(frame_name.c_str(),"w");     //Per-frame contact-fraction output stream
	FILE *block_file  = fopen64(block_name.c_str(),"w");     //Block-average output stream
        FILE *delta_file  = fopen64(delta_name.c_str(),"w");     //Per-type block-change output stream

        if(block_file == NULL || frame_file == NULL || delta_file == NULL)
        {
            printf("Failure opening Lipid Neighbor Content output files. Check the -lnc path.\n");
            MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
        }

        //print header info (type pairs POPC-POPC   POPC-POPG)
        fprintf(frame_file," %10s","#frame");
        fprintf(block_file," %10s","#block");
        fprintf(delta_file," %10s","#block_i");
        for(i=0; i<param.main_size_y(); i++) //loop over lipid types i 
        {
            for(j=0; j<param.main_size_y(); j++) //loop over lipid types j 
            {
		string pair = "#" + param.param_main_s[i][0] + "-" + param.param_main_s[j][0].c_str(); 
                fprintf(frame_file," %10s",pair.c_str());
                fprintf(block_file," %10s",pair.c_str());
            }
            string this_delta = "#"; 
	    this_delta = this_delta + param.param_main_s[i][0].c_str();
            fprintf(delta_file," %10s",this_delta.c_str());
        }
        fprintf(frame_file,"\n");
        fprintf(block_file,"\n");
        fprintf(delta_file,"\n");

        //print frame data
        for(k=0; k<type_comp_frames[0][0].size(); k++) //loop over the frames
        {
            fprintf(frame_file," %10d",p.range + k);
            for(i=0; i<param.main_size_y(); i++)    //loop over lipid types i
            {
                for(j=0; j<param.main_size_y(); j++) //loop over lipid types j
                {
                    fprintf(frame_file," %10f",type_comp_frames[i][j][k]);
                }
            }
            fprintf(frame_file,"\n");
        }
        fclose(frame_file);

        //compute block average data
        dv2d sum(param.main_size_y(),dv1d(param.main_size_y(),0.0));         //stores a sum over the blcok for each type pair
        dv3d avg(0,dv2d(param.main_size_y(),dv1d(param.main_size_y(),0.0))); //stores a avg over the blcok for each type pair 

        for(k=0; k<type_comp_frames[0][0].size(); k++) //loop over frames
        {
            for(i=0; i<param.main_size_y(); i++) //loop over lipid type i
            {
                for(j=0; j<param.main_size_y(); j++) //loop over lipid type j
                {   
                    sum[i][j] = sum[i][j] + type_comp_frames[i][j][k];
                }
            }

            if((k+1)%p.block == 0) //block completed
            {
                dv2d this_avg(param.main_size_y(),dv1d(param.main_size_y(),0.0));
                for(i=0; i<param.main_size_y(); i++) //loop over lipid type i
                {
                    for(j=0; j<param.main_size_y(); j++) //loop over lipid type j
                    { 
                        this_avg[i][j] = sum[i][j]/(double)p.block;
                        sum[i][j]      = 0.0;
                    }
                }
                avg.push_back(this_avg);
            }
        }

        //print block average data
        for(k=0; k<avg.size(); k++) //loop over the blocks
        {
            fprintf(block_file," %10d",k);
            for(i=0; i<param.main_size_y(); i++) //loop over lipid type i
            {
                for(j=0; j<param.main_size_y(); j++) //loop over lipid type j
                {
                    fprintf(block_file," %10f",avg[k][i][j]);
                }
            } 
            fprintf(block_file,"\n");
        }
        fclose(block_file);

        //compute deltas between blocks
        for(k=1; k<avg.size(); k++) //loop over block neighbors
        {
            fprintf(delta_file," %10d",k-1);

            for(i=0; i<param.main_size_y(); i++) //loop over lipid types i
            {
                dv1d sum(param.main_size_y(),0.0) ;

                for(j=0; j<param.main_size_y(); j++) //loop over lipid types j
                {
                    //average data
                    double delta = fabs(avg[k][i][j] - avg[k-1][i][j]);
 
                    sum[i] = sum[i] + delta;
                }
                sum[i] = sum[i]/param.main_size_y();

                //print delta
                fprintf(delta_file," %10f",sum[i]);
            }
            fprintf(delta_file,"\n");
        }
        fclose(delta_file);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    return (clock() - t_i)/(double)CLOCKS_PER_SEC;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                           //
// This function configures the native scaffold, streams assigned frame blocks and finalizes the analysis.   //
//                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, const char * argv[])
{
    //Initialize the MPI environment and retain the MosAT scaffold's system bookkeeping.
    MPI_Init(NULL,NULL);

    system_variables s;                           
    program_variables p;                          
    initialize_system_variables(&s);
    initialize_program_variables(&p);

    Performance perf;                            
    s.program_name = "Lipid Neighbor Content";
    enum Switch serial = off;                     
    check_serial(s.world_rank,s.world_size,serial);
    name_and_record(s.world_rank,argc,argv,s.program_name);

    //Register the native trajectory controls and the observable's short command-line tags.
    start_input_arguments_mpi(argc,argv,s.world_rank,p.program_description);
    add_argument_mpi_s(argc,argv,"-traj",   p.in_file_name,                "Input trajectory file (xtc, trr, pdb, gro)",                s.world_rank,s.cl_tags,nullptr,       1);
    add_argument_mpi_s(argc,argv,"-ref",    p.ref_file_name,               "Reference file (pdb, gro)",                                 s.world_rank,s.cl_tags,nullptr,       1);
    add_argument_mpi_i(argc,argv,"-stride", &p.stride,                     "Read every stride-th original frame (default 1)",           s.world_rank,s.cl_tags,nullptr,       0);
    add_argument_mpi_i(argc,argv,"-b",      &p.start_frame,                "Skip original frames before this number (default 0)",       s.world_rank,s.cl_tags,nullptr,       0);
    add_argument_mpi_i(argc,argv,"-e",      &p.end_frame,                  "Skip original frames after this number (default -1)",       s.world_rank,s.cl_tags,nullptr,       0);
    add_argument_mpi_s(argc,argv,"-crd",    p.param_file_name,             "Networked lipid type/cutoff/atom-card selection (crd)",     s.world_rank,s.cl_tags,nullptr,       1);
    add_argument_mpi_s(argc,argv,"-lnc",    p.lnc_file_name,               "Block averages and output filename base (dat)",             s.world_rank,s.cl_tags,nullptr,       1);
    add_argument_mpi_i(argc,argv,"-leaf",   &p.leaflet,                    "Target leaflet (0:full membrane, 1:upper, 2:lower)",        s.world_rank,s.cl_tags,nullptr,       1);
    add_argument_mpi_i(argc,argv,"-block",  &p.block,                      "Filtered central frames per averaging block (positive)",    s.world_rank,s.cl_tags,nullptr,       1);
    add_argument_mpi_i(argc,argv,"-range",  &p.range,                      "Centered noise-filter half-width (frames; nonnegative)",    s.world_rank,s.cl_tags,nullptr,       1);
    add_argument_mpi_s(argc,argv,"-lf_pdb", p.lf_pdb_file_name,            "Optional PDB with marked leaflets (pdb)",                   s.world_rank,s.cl_tags,&p.b_lf_pdb,   0);
    add_argument_mpi_s(argc,argv,"-lf_prm", p.leaflet_finder_param_name,   "Optional additional leaflet-finder parameters (prm)",       s.world_rank,s.cl_tags,&p.b_lf_param, 0);
    conclude_input_arguments_mpi(argc,argv,s.world_rank,s.program_name,s.cl_tags);

    Trajectory traj;                              
    traj.set_block_parallel(off);
    traj.set_traj(p.in_file_name);
    traj.set_ref(p.ref_file_name);
    traj.set_lsq("",0,3,0);                        
    traj.set_res(p.stride,p.start_frame,p.end_frame);

    perf.log_time(traj.build(),"Analyze Trajectory");

    //Distribute contributing frames without restricting this object's full effective-frame access.
    int num_frames = traj.get_ef_frames();        
    traj.parallelize_by_frames(num_frames);

    //Check specified extensions through the native interface and read the native networked card.
    check_extension_mpi(s.world_rank,"-crd",p.param_file_name,".crd");
    check_extension_mpi(s.world_rank,"-lnc",p.lnc_file_name,".dat");
    if(p.b_lf_pdb == 1)
    {
        check_extension_mpi(s.world_rank,"-lf_pdb",p.lf_pdb_file_name,".pdb");
    }
    if(p.b_lf_param == 1)
    {
        check_extension_mpi(s.world_rank,"-lf_prm",p.leaflet_finder_param_name,".prm");
    }

    Param param;                                  //Native networked selection: three main columns, filenames in column 2
    param.get_param(p.param_file_name,3,2,1);
    if(param.check_file() == 0)
    {
        MPI_Finalize();
        return EXIT_FAILURE;
    }

    traj.get_leaflets(p.leaflet,p.leaflet_finder_param_name,p.b_lf_param);
    traj.write_leaflets(p.lf_pdb_file_name,p.b_lf_pdb);
    traj.get_leaflet_stats();
    traj.get_lipid_selection_stats(param.get_column_s(0),"-crd");
    p.num_lipids = traj.count_target_lipids_type(param.get_column_s(0));

    //Require a full centered window and a positive statistical block length.
    if(p.block <= 0 || param.main_size_y() == 0 || traj.get_num_frames() == 0 || p.range < 0 || p.range > (traj.get_num_frames() - 1)/2)
    {
        if(s.world_rank == 0)
        {
            printf("Require -block > 0, nonempty -crd, -range >= 0, and at least 2*range + 1 effective frames.\n");
        }
        MPI_Finalize();
        return EXIT_FAILURE;
    }
    p.window     = 2*p.range + 1;

    //Extend the assigned frame block by the filter range, clipping to the full trajectory.
    if(traj.my_frames_alt > 0)
    {
        p.read_first   = traj.frame_start_alt - p.range < 0
                       ? 0 : traj.frame_start_alt - p.range;
        p.read_last    = traj.frame_end_alt + p.range > traj.get_num_frames() - 1
                       ? traj.get_num_frames() - 1 : traj.frame_end_alt + p.range;
    }
    else //Set an empty interval so this rank skips every frame in the trajectory-reading loop
    {
        p.read_first = -1;
        p.read_last  = -2;
    }

    //Allocate the frame series, linked atom selections and pair histories.
    int frame_count = 0;
    if(p.read_last - (p.read_first-1) > 2*p.range)
    {
        frame_count = p.read_last-(p.read_first-1)-(2*p.range);
    }

    dv3d type_comp_frames(param.main_size_y(),dv2d(param.main_size_y(),dv1d(frame_count,0)));        //Stored type_comp indexed by central type, local frame, neighbor type
    sv2d center_atoms(param.main_size_y(),sv1d(0));                                                  //Native atom-name selection for each lipid type
    iv3d filter(p.num_lipids,iv2d(p.num_lipids,iv1d(2*p.range + 1,0)));                              //Raw pair history indexed by central lipid, neighbor lipid, window slot
    iv2d neighbors_stats(param.main_size_y(),iv1d(0,0));                                             //Holds statistics for how many neighbors each lipid type had 

    for(int i=0; i<param.main_size_y(); i++)
    {
        center_atoms[i] = param.get_column_sec_s(i,0);
    }

    traj.workload_frames();
    traj.report_progress();

    s.t = clock();

    for(traj.current_frame=0; traj.current_frame<traj.get_num_frames(); traj.current_frame++) //Full effective-frame index range
    {
        if(traj.current_frame >= p.read_first && traj.current_frame <= p.read_last)
        {
            traj.read_traj_frame();

            lipid_neighbor_content(traj,p,param,center_atoms,filter,type_comp_frames,neighbors_stats);

            time_stats(s.t,&s.counter,traj.current_frame - p.read_first,p.read_last - p.read_first + 1,s.world_rank);
        }
    }

    perf.log_time((clock() - s.t)/CLOCKS_PER_SEC,"Main Loop");
    perf.log_time(finalize_analysis(traj,s,p,param,type_comp_frames,neighbors_stats),"Finalize Analysis");

    perf.print_stats();
    print_closing(s.world_rank);
    MPI_Finalize();
    return 0;
}
