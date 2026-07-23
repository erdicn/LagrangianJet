#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include "simu.hpp"

void initCleanLagSimu(int argc, char** args, LagSimuParams_t* simu) {
    assert(sizeof(myfloat) == sizeof(double));

    // 1. Set default values first
    simu->out_folder_name = "None/None";
    simu->µs_to_drop = -1;
    simu->dt = -1;
    simu->objective_x = -1; // mm then it is converted to reality 
    simu->nb_to_save = -1;

    // 2. Loop through arguments and parse key=value directly
    for (int i = 1; i < argc; i++) {
        if (sscanf(args[i], "injection_rate=%lf", &simu->µs_to_drop) == 1) continue;
        if (sscanf(args[i], "dt=%lf", &simu->dt) == 1) continue;
        if (sscanf(args[i], "recording_distance_in_mm=%lf", &simu->objective_x) == 1) continue;
        if (sscanf(args[i], "nb_droplets_saved=%llu", &simu->nb_to_save) == 1) continue;

        // If argument has no '=' sign, assume it's the folder path
        if (strchr(args[i], '=') == NULL) {
            simu->out_folder_name = args[i];
        }
    }

    myfloat dist_mm = simu->objective_x;
    char base_folder[256];
    snprintf(base_folder, sizeof(base_folder), "%s", simu->out_folder_name.c_str());

    assert(simu->objective_x >= 11 && "The droplets stats start from 11mm from experiments so need to put something bigger"); 

    simu->objective_x = (simu->objective_x-11)*1e3; // from mm to µm

    char folder_buf[512];
    snprintf(folder_buf, sizeof(folder_buf),
             "%s_inj%.4f_dt%.4f_dist%.0f_n%llu",
             base_folder,
             simu->µs_to_drop,
             simu->dt,
             dist_mm,
             (unsigned long long)simu->nb_to_save);

    // Save formatted path to struct
    simu->out_folder_name = folder_buf;

    createOutFolder(simu->out_folder_name.c_str());

    // Sanity check
    assert(simu->dt < 1);

    // Initialize remaining simulation state
    simu->current_saved = 0;
    simu->t = 0;
    simu->t_ind = 0; 
    simu->save_to_file_every_n_droplet = simu->nb_to_save / 100;
    assert(simu->save_to_file_every_n_droplet >= 1 && "Please enter at least 100 here");
    simu->save_file_i = 0;

    // Feedback printout for students
    printf("\n--- Simulation Config ---\n");
    printf(" Folder     : %s\n", simu->out_folder_name.c_str());
    printf(" Inj Rate   : %.4f µs\n", simu->µs_to_drop);
    printf(" Time Step  : %.4f\n", simu->dt);
    printf(" Rec Dist   : %.2f mm\n", simu->objective_x);
    printf(" To Save    : %llu droplets\n", (unsigned long long)simu->nb_to_save);
    printf("-------------------------\n\n");

    simu->next_drop_att_ind = 0;
    simu->x_vec         = allocateFVec(NULL, 1000);
    simu->r_vec         = allocateFVec(NULL, 1000);
    simu->v_vec         = allocateFVec(NULL, 1000);
    simu->r_save_buffer = allocateFVec(NULL, simu->save_to_file_every_n_droplet);
    simu->v_save_buffer = allocateFVec(NULL, simu->save_to_file_every_n_droplet);
    simu->r_save_i = 0;
    simu->current_max_index = 0;
    simu->batch_start_time = std::chrono::high_resolution_clock::now();
}