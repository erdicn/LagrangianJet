#include <cstdlib> 
#include <cassert>
#include <cstdint>

#include <iostream>
#include <string>
#include <chrono>
#include <sstream>

#ifdef _OPENMP
#   include <omp.h>
#endif

#include "matrix.h"
#include "droplet_stats.hpp"
#include "file_manips.hpp"
#include "OpenMP/test_openmp.h"
#include "simu.hpp"

#include <sys/stat.h>  // For mkdir()
#include <sys/types.h> // For mode_t
#include <errno.h>     // For errno and EEXIST
#include <stdio.h>     // For printf and perror
// #define OUT_FOLDER_NAME "Out"

#define RHO_WATER 998
#define GRAVITY 9.81f
#define RHO_AIR 1.27
#define NU_AIR 1.5e-5
#define NU_WATER 1e-3
#define FRAC_4_3 1.333333333333333333333333333333333333
#ifndef M_PI
    #define M_PI 3.141592653589793238462643383279502984
#endif

#define ENERGY_DISIPATION 0
#define ENERGY_DISSIPATION_COEF 0
#define WE_CRIT_BOUNCE 0
#define SIGMA_WATER 0.072  
#define COLISION_DIST_MULTIPLIER 5
// TODO parametric study witd tcol dia


void createOutFolder(const char* folder_name){
    // const char* folder_name = OUT_FOLDER_NAME;
    // Try to create the directory with read/write/execute permissions (0777)
    int result = mkdir(folder_name, 0777);
    
    if (result == -1) {
        // check WHY it failed.
        if (errno == EEXIST) {
            printf("The folder '%s' already exists. That's fine!\n", folder_name);
        } else {
            // It failed for some other reason 
            perror("Error creating directory");
            return; 
        }
    } else {
        printf("Successfully created the folder '%s'!\n", folder_name);
    }
}


// inline double diameterToVolume(double dia){
//     return 4./3.*M_PI*pow(dia/2., 3);
// }

// inline double radiusToVolume(double r){
//     return 4./3.*M_PI*MYPOW(r, 3);
// }

// ·
inline bool collisionTest(myfloat x1, myfloat x2, 
                          myfloat r1, myfloat r2){
    return (MYABS(x1-x2) < COLISION_DIST_MULTIPLIER*(r1+r2));
}

// inline myfloat momentumConservation(myfloat m1, myfloat m2, myfloat s1, myfloat s2){
//     return (m1*s1 + m2*s2)/(m1+m2); 
// }

inline myfloat funcReynoldsNbSphereInAir(myfloat r, myfloat vel){
    return RHO_AIR*vel*2*r/NU_AIR;
}

inline myfloat funcReynoldsNbSphereInAirµvars(myfloat r, myfloat vel){
    return RHO_AIR*vel*2*r*1e-6/NU_AIR;
}

#define DELTA0 9.4
#define C0     24./(DELTA0*DELTA0)
inline myfloat funcDragCOefIncompressible(myfloat Re){
    return C0* MYPOW(1+DELTA0/MYSQRT(Re) , 2);
}

typedef double cfloat;  // c stands for custom
typedef uint32_t cint; // c stands for custom


void checkArguments(int argc, char** args){
    puts("Arguments entered:");
    for(int i = 0; i < argc; i++){
        printf("%s ", args[i]);
    }puts("");
    assert(argc >= 6); // TODO need to do this even in the non debug versions
    puts("Atention every dimension should be in µ");
    return;
}

extern "C" int getSimuMaxIndex(LagSimuParams_t* simu) {
    return simu->current_max_index;
}

extern "C" fvec_t* getSimuRVec(LagSimuParams_t* simu) {
    return simu->r_vec;
}

extern "C" fvec_t* getSimuXVec(LagSimuParams_t* simu) {
    return simu->x_vec;
}

extern "C" fvec_t* getSimuVVec(LagSimuParams_t* simu) {
    return simu->v_vec;
}

extern "C" myfloat getSimuObjX(LagSimuParams_t* simu){
    return simu->objective_x;
}


// typedef struct LagSimuParams{
//     std::string out_folder_name;
//     myfloat µs_to_drop;
//     myfloat dt;
//     myfloat objective_x;
//     uint64_t nb_to_save;
//     uint64_t current_saved ;
//     double t ;
//     uint64_t t_ind ; 
//     uint32_t save_to_file_every_n_droplet;
//     uint32_t save_file_i;
//     uint64_t next_drop_att_ind = 0;

//     fvec_t* x_vec      ;//= allocateFVec(NULL, 1000);
//     fvec_t* r_vec      ;//= allocateFVec(NULL, 1000);
//     fvec_t* v_vec      ;//= allocateFVec(NULL, 1000);

//     fvec_t* r_save_buffer ;//= allocateFVec(NULL, simu->save_to_file_every_n_droplet);
//     fvec_t* v_save_buffer ;//= allocateFVec(NULL, simu->save_to_file_every_n_droplet);
//     int r_save_i;// = 0;

//     int current_max_index;// = 0;

//     std::chrono::time_point<std::chrono::high_resolution_clock> batch_start_time;//= std::chrono::high_resolution_clock::now();
// } LagSimuParams_t;

#include "simu.hpp"


extern "C" LagSimuParams_t* allocateSimuParams() {
    return new LagSimuParams_t(); 
}

void initCleanLagSimu(int argc, char** args, LagSimuParams_t* simu){
    assert(sizeof(myfloat) == sizeof(double));
    checkArguments(argc, args);
    // testOpeznMP(16);
    int a = 1;
    simu->out_folder_name = args[a++]; 
    createOutFolder(simu->out_folder_name.c_str()); // TODO add varable name 
    
    simu->µs_to_drop = atof(args[a++]);
    simu->dt = atof(args[a++]);
    assert(simu->dt < 1);
    simu->objective_x = atof(args[a++]); 
    simu->nb_to_save = std::stoll(args[a++]);
    simu->current_saved = 0;
    simu->t = 0;
    simu->t_ind = 0; 
    simu->save_to_file_every_n_droplet = simu->nb_to_save/100;
    assert(simu->save_to_file_every_n_droplet > 1);
    simu->save_file_i = 0;

    printf("dt = %.16lf, saving_x = %.16lf, to save = %ld\n", simu->dt, simu->objective_x, simu->nb_to_save);

    #ifdef ENERGY_DISIPATION
        puts("Energy is getting disipated in colisions");
    #endif // ENERGY_DISIPATION

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

// TODO add batch start time 
void mainLoopCleanLagSimu(LagSimuParams_t* simu, uint32_t nb_frame_simul){
    auto nb_to_save = simu->nb_to_save;
    auto x_vec = simu->x_vec;
    auto r_vec = simu->r_vec;
    auto v_vec = simu->v_vec;
    auto dt = simu->dt;
    auto µs_to_drop = simu->µs_to_drop;
    uint32_t while_i = 0;
    while (simu->current_saved < nb_to_save && while_i++ < nb_frame_simul){
        if(simu->t_ind == simu->next_drop_att_ind || simu->t_ind > simu->next_drop_att_ind){
            double dia    = distribution11mm(gen);
            // double dia    = uniform_dist_dia_100_165um(gen);
            // double volume = diameterToVolume(dia);
            // TODO 
            // TODO remove for exp
            // double vel  = uniform_dist_vel_9_11(gen);
            double vel  = dia < 75 ? getLaplaceRandom(sat_mu, sat_b, gen) : normal_dist11mm(gen);
            //                         second     to        µs    
            // double µs_to_drop = dia < 75 ? getLaplaceRandom(46000, sat_mu, gen) : normal_distfreq(gen);
            // µs_to_drop = 1./µs_to_drop * 1e6;
            // µs_to_drop = getRandomFloat(15, 30);
            simu->next_drop_att_ind +=  round(µs_to_drop/dt);
            x_vec     ->vals[simu->current_max_index] = 0;
            // volume_vec->vals[current_max_index] = volume;
            r_vec     ->vals[simu->current_max_index] = dia/2.;
            v_vec     ->vals[simu->current_max_index] = vel;
            simu->current_max_index++;

            if(simu->current_max_index == x_vec->len){
                reallocateFVec(x_vec     , simu->current_max_index*2);
                // reallocateFVec(volume_vec, current_max_index*2);
                reallocateFVec(r_vec     , simu->current_max_index*2);
                reallocateFVec(v_vec     , simu->current_max_index*2);
                std::cout << "Realocating vector to size " << simu->current_max_index*2 << std::endl<<std::flush;
            }
            // std::cout << t_ind << " " << next_drop_att_ind << std::endl << std::flush;
        } 

        for(int i = 0; i < simu->current_max_index; i++){
            x_vec->vals[i] += v_vec->vals[i]*dt; 
        }

        bool merges_occurred;
        do{ 
            /*
            merges_occurred = false;
            for(int i = 0; i < simu->current_max_index - 1; i++){
                myfloat x1 = x_vec->vals[i];
                myfloat x2 = x_vec->vals[i+1];
                myfloat r1 = r_vec->vals[i];
                myfloat r2 = r_vec->vals[i+1];
                
                if(r1 > 0 && r2 > 0 && collisionTest(x1, x2, r1, r2)){
                    myfloat vol1 = radiusToVolume(r_vec->vals[i]);
                    myfloat vol2 = radiusToVolume(r_vec->vals[i+1]);
                    myfloat v1 = v_vec->vals[i];
                    myfloat v2 = v_vec->vals[i+1];
                    
                    x_vec->vals[i]      = momentumConservation(vol1, vol2, x1, x2);
                    v_vec->vals[i]      = momentumConservation(vol1, vol2, v1, v2);
                    #ifdef ENERGY_DISIPATION
                        myfloat S = 4*vol1*vol2/MYPOW(vol1+vol2, 2); // scaling factor
                        v_vec->vals[i] *= (1-ENERGY_DISSIPATION_COEF*S);
                    #endif // ENERGY_DISIPATION
                    r_vec->vals[i]      = volumeToRadius(vol1 + vol2);
                    // volume_vec->vals[i] = vol1 + vol2;
                    
                    // Mark the consumed droplet as dead by setting radius to 0
                    r_vec->vals[i+1] = 0; 
                    merges_occurred = true; // this is in case of if we have multiple droplets coliding at once 
                    
                }
            }
            */

            merges_occurred = false;
            for(int i = 0; i < simu->current_max_index - 1; i++){
                myfloat x1 = x_vec->vals[i];
                myfloat x2 = x_vec->vals[i+1];
                myfloat r1 = r_vec->vals[i];
                myfloat r2 = r_vec->vals[i+1];
                
                if(r1 > 0 && r2 > 0 && collisionTest(x1, x2, r1, r2)){
                    myfloat vol1 = radiusToVolume(r1);
                    myfloat vol2 = radiusToVolume(r2);
                    myfloat v1 = v_vec->vals[i];
                    myfloat v2 = v_vec->vals[i+1];
                    
                    // 1. Calculate the Weber Number
                    // Radii must be converted from µm to meters for the We calculation
                    myfloat d_small = 2.0 * std::min(r1, r2) * 1e-6; 
                    myfloat u_rel = MYABS(v1 - v2);
                    myfloat We = RHO_WATER * u_rel * u_rel * d_small / SIGMA_WATER;

                    if (We < WE_CRIT_BOUNCE) {
                        // BOUNCE REGIME 
                        myfloat CR = 1-ENERGY_DISSIPATION_COEF; // Coefficient of restitution (1.0 = perfect elastic, < 1.0 = energy dissipated)
                        
                        // Calculate new velocities based on 1D inelastic collision equations
                        // Mass is directly proportional to volume, so we can use vol ratios directly
                        myfloat v1_new = (vol1 * v1 + vol2 * v2 + vol2 * CR * (v2 - v1)) / (vol1 + vol2);
                        myfloat v2_new = (vol1 * v1 + vol2 * v2 + vol1 * CR * (v1 - v2)) / (vol1 + vol2);
                        
                        v_vec->vals[i]   = v1_new;
                        v_vec->vals[i+1] = v2_new;

                        // Anti-Overlap Fix: Push droplets apart slightly to prevent infinite collision loops next time step.
                        // We push them out relative to their center of mass
                        myfloat center_mass = (vol1 * x1 + vol2 * x2) / (vol1 + vol2);
                        myfloat safe_dist = r1 + r2 + 1e-6; // Exact radius + tiny epsilon spacing

                        if (x1 < x2) {
                            x_vec->vals[i]   = center_mass - safe_dist * (vol2 / (vol1 + vol2));
                            x_vec->vals[i+1] = center_mass + safe_dist * (vol1 / (vol1 + vol2));
                        } else {
                            x_vec->vals[i]   = center_mass + safe_dist * (vol2 / (vol1 + vol2));
                            x_vec->vals[i+1] = center_mass - safe_dist * (vol1 / (vol1 + vol2));
                        }
                        
                        // Notice we DO NOT set merges_occurred = true, and we DO NOT kill r2.

                    } else {
                        // COALESCENCE REGIME 
                        
                        x_vec->vals[i]      = momentumConservation(vol1, vol2, x1, x2);
                        v_vec->vals[i]      = momentumConservation(vol1, vol2, v1, v2);
                        
                        #ifdef ENERGY_DISIPATION
                            myfloat S = 4*vol1*vol2/MYPOW(vol1+vol2, 2); 
                            v_vec->vals[i] *= (1-ENERGY_DISSIPATION_COEF*S);
                        #endif 
                        
                        r_vec->vals[i]      = volumeToRadius(vol1 + vol2);
                        
                        // Mark the consumed droplet as dead by setting radius to 0
                        r_vec->vals[i+1] = 0; 
                        merges_occurred = true; 
                    }
                }
            }

            int write_idx = 0;
            for (int read_idx = 0; read_idx < simu->current_max_index; read_idx++) {
                if (r_vec->vals[read_idx] > 0) { // If the droplet is alive
                    x_vec->vals[write_idx]      = x_vec->vals[read_idx];
                    v_vec->vals[write_idx]      = v_vec->vals[read_idx];
                    // volume_vec->vals[write_idx] = volume_vec->vals[read_idx];
                    r_vec->vals[write_idx]      = r_vec->vals[read_idx];
                    write_idx++;
                }
            }
            // Update the new size of our active arrays
            simu->current_max_index = write_idx;
        } while(merges_occurred);


        while (simu->current_max_index > 0 && x_vec->vals[0] > simu->objective_x){
            simu->v_save_buffer->vals[simu->r_save_i  ] = v_vec->vals[0];
            simu->r_save_buffer->vals[simu->r_save_i++] = r_vec->vals[0];
            for(int i = 0; i < simu->current_max_index-1; i++){
                x_vec     ->vals[i]      = x_vec->vals[i+1];
                v_vec     ->vals[i]      = v_vec->vals[i+1];
                r_vec     ->vals[i]      = r_vec->vals[i+1];
            } 
            simu->current_max_index--;
        }

        if(simu->r_save_i >= simu->r_save_buffer->len) {
            std::stringstream ss;
            ss << simu->out_folder_name << "/radius" << simu->save_file_i << ".dat";
            saveArrayToFile(simu->r_save_buffer->vals, simu->r_save_buffer->len, ss.str().c_str());

            ss.str(std::string());
            ss << simu->out_folder_name << "/velocity" << simu->save_file_i << ".dat";
            saveArrayToFile(simu->v_save_buffer->vals, simu->v_save_buffer->len, ss.str().c_str());
            simu->save_file_i++;
            simu->current_saved += simu->r_save_i;
            
            auto batch_end_time = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = batch_end_time - simu->batch_start_time;
            printf("Calculated and saved %u elements in %.4f seconds (%.2f elements/sec)\n", 
                    simu->r_save_i, elapsed.count(), simu->r_save_i / elapsed.count());
            simu->batch_start_time = std::chrono::high_resolution_clock::now();

            simu->r_save_i = 0;
        }
        simu->t_ind++;
    }
}

/*
//     folder F   dt  obj_x  nb_save
// g++ Out11  30 0.01 11000  1000000
int mainCleanLagSimu(int argc, char** args){
    assert(sizeof(myfloat) == sizeof(double));
    checkArguments(argc, args);
    testOpenMP(16);
    // #pragma omp 
    // std::cout << omp_get_num_threads() << std::endl << std::flush;
    // #pragma omp parallel for
    // for(int i = 0; i < omp_get_num_threads(); i++){
    //     printf("Thread %d of %d is processing index %d\n", omp_get_thread_num(), omp_get_num_threads(), i);
    // }
    int a = 1;
    std::string out_folder_name = args[a++]; 
    createOutFolder(out_folder_name.c_str()); // TODO add varable name 
    
    double µs_to_drop = atof(args[a++]);
    double dt = atof(args[a++]);
    assert(dt < 1);
    double objective_x = atof(args[a++]); 
    uint64_t nb_to_save = std::stoll(args[a++]);
    uint64_t current_saved = 0;
    double t = 0;
    uint64_t t_ind = 0; 
    uint32_t save_to_file_every_n_droplet = nb_to_save/100;
    assert(save_to_file_every_n_droplet > 1);
    int save_file_i = 0;

    printf("dt = %.16lf, saving_x = %.16lf, to save = %ld\n", dt, objective_x, nb_to_save);

    uint64_t next_drop_att_ind = 0;

    fvec_t* x_vec      = allocateFVec(NULL, 1000);
    fvec_t* r_vec      = allocateFVec(NULL, 1000);
    fvec_t* v_vec      = allocateFVec(NULL, 1000);

    fvec_t* r_save_buffer = allocateFVec(NULL, save_to_file_every_n_droplet);
    fvec_t* v_save_buffer = allocateFVec(NULL, save_to_file_every_n_droplet);
    int r_save_i = 0;

    int current_max_index = 0;


    std::cout << "Entering main loop\n" << std::flush;
    auto batch_start_time = std::chrono::high_resolution_clock::now();

    while (current_saved < nb_to_save){
        if(t_ind == next_drop_att_ind || t_ind > next_drop_att_ind){
            double dia    = distribution11mm(gen);
            double volume = diameterToVolume(dia);
            double vel    = dia < 75 ? getLaplaceRandom(sat_mu, sat_b, gen) : normal_dist11mm(gen);
            //                         second     to        µs    
            // double µs_to_drop = dia < 75 ? getLaplaceRandom(46000, sat_mu, gen) : normal_distfreq(gen);
            // µs_to_drop = 1./µs_to_drop * 1e6;
            // µs_to_drop = getRandomFloat(15, 30);
            next_drop_att_ind +=  round(µs_to_drop/dt);
            x_vec     ->vals[current_max_index] = 0;
            // volume_vec->vals[current_max_index] = volume;
            r_vec     ->vals[current_max_index] = dia/2.;
            v_vec     ->vals[current_max_index] = vel;
            current_max_index++;

            if(current_max_index == x_vec->len){
                reallocateFVec(x_vec     , current_max_index*2);
                // reallocateFVec(volume_vec, current_max_index*2);
                reallocateFVec(r_vec     , current_max_index*2);
                reallocateFVec(v_vec     , current_max_index*2);
                std::cout << "Realocating vector to size " << current_max_index*2 << std::endl<<std::flush;
            }
            // std::cout << t_ind << " " << next_drop_att_ind << std::endl << std::flush;
        } 

        // faster without paralelisation the time took to dispatch takes more time to calculate it// #pragma omp parallel for 
        for(int i = 0; i < current_max_index; i++){
            // myfloat acceleration = GRAVITY*1e-6 - 0.5*RHO_AIR*MYPOW(v_vec->vals[i], 2)*
            //                                      funcDragCOefIncompressible(funcReynoldsNbSphereInAirµvars(r_vec->vals[i], v_vec->vals[i]))*
            //                                      M_PI*MYPOW(r_vec->vals[i], 2)*1e-12;
            // v_vec->vals[i] += acceleration * dt;
            // // std::cout <<  v_vec->vals[i] << std::endl;
            x_vec->vals[i] += v_vec->vals[i]*dt; 
        }

        bool merges_occurred;
        do{ 
            merges_occurred = false;
            for(int i = 0; i < current_max_index - 1; i++){
                myfloat x1 = x_vec->vals[i];
                myfloat x2 = x_vec->vals[i+1];
                myfloat r1 = r_vec->vals[i];
                myfloat r2 = r_vec->vals[i+1];
                
                if(r1 > 0 && r2 > 0 && collisionTest(x1, x2, r1, r2)){
                    myfloat vol1 = radiusToVolume(r_vec->vals[i]);
                    myfloat vol2 = radiusToVolume(r_vec->vals[i+1]);
                    myfloat v1 = v_vec->vals[i];
                    myfloat v2 = v_vec->vals[i+1];
                    
                    x_vec->vals[i]      = momentumConservation(vol1, vol2, x1, x2);
                    v_vec->vals[i]      = momentumConservation(vol1, vol2, v1, v2);
                    r_vec->vals[i]      = volumeToRadius(vol1 + vol2);
                    // volume_vec->vals[i] = vol1 + vol2;
                    
                    // Mark the consumed droplet as dead by setting radius to 0
                    r_vec->vals[i+1] = 0; 
                    merges_occurred = true;
                }
            }

            int write_idx = 0;
            for (int read_idx = 0; read_idx < current_max_index; read_idx++) {
                if (r_vec->vals[read_idx] > 0) { // If the droplet is alive
                    x_vec->vals[write_idx]      = x_vec->vals[read_idx];
                    v_vec->vals[write_idx]      = v_vec->vals[read_idx];
                    // volume_vec->vals[write_idx] = volume_vec->vals[read_idx];
                    r_vec->vals[write_idx]      = r_vec->vals[read_idx];
                    write_idx++;
                }
            }
            // Update the new size of our active arrays
            current_max_index = write_idx;
        }while(merges_occurred);


        while (current_max_index > 0 && x_vec->vals[0] > objective_x){
            v_save_buffer->vals[r_save_i  ] = v_vec->vals[0];
            r_save_buffer->vals[r_save_i++] = r_vec->vals[0];
            for(int i = 0; i < current_max_index-1; i++){
                x_vec     ->vals[i]      = x_vec->vals[i+1];
                v_vec     ->vals[i]      = v_vec->vals[i+1];
                r_vec     ->vals[i]      = r_vec->vals[i+1];
            } 
            current_max_index--;
        }

        if(r_save_i >= r_save_buffer->len) {
            std::stringstream ss;
            ss << out_folder_name << "/radius" << save_file_i << ".dat";
            saveArrayToFile(r_save_buffer->vals, r_save_buffer->len, ss.str().c_str());

            ss.str(std::string());
            ss << out_folder_name << "/velocity" << save_file_i << ".dat";
            saveArrayToFile(v_save_buffer->vals, v_save_buffer->len, ss.str().c_str());
            save_file_i++;
            current_saved += r_save_i;
            
            auto batch_end_time = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = batch_end_time - batch_start_time;
            printf("Calculated and saved %u elements in %.4f seconds (%.2f elements/sec)\n", 
                    r_save_i, elapsed.count(), r_save_i / elapsed.count());
            batch_start_time = std::chrono::high_resolution_clock::now();

            r_save_i = 0;
        }
        t_ind++;
    }

    return 0;
}
*/

