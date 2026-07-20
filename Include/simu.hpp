#ifndef SIMU_HPP
#define SIMU_HPP

#include "matrix.h"
#include "myfloat.h"
#include "stdint.h"


#ifdef __cplusplus
#include <string>
#include <chrono>
// Full definition ONLY visible to C++
typedef struct LagSimuParams {
    std::string out_folder_name;
    myfloat µs_to_drop;
    myfloat dt;
    myfloat objective_x;
    uint64_t nb_to_save;
    uint64_t current_saved ;
    double t ;
    uint64_t t_ind ; 
    uint32_t save_to_file_every_n_droplet;
    uint32_t save_file_i;
    uint64_t next_drop_att_ind = 0;

    fvec_t* x_vec      ;//= allocateFVec(NULL, 1000);
    fvec_t* r_vec      ;//= allocateFVec(NULL, 1000);
    fvec_t* v_vec      ;//= allocateFVec(NULL, 1000);

    fvec_t* r_save_buffer ;//= allocateFVec(NULL, simu->save_to_file_every_n_droplet);
    fvec_t* v_save_buffer ;//= allocateFVec(NULL, simu->save_to_file_every_n_droplet);
    int r_save_i;// = 0;

    int current_max_index;// = 0;
    std::chrono::time_point<std::chrono::high_resolution_clock> batch_start_time;
} LagSimuParams_t;

extern "C" {
#else
// Forward declaration ONLY visible to C
typedef struct LagSimuParams LagSimuParams_t;
#endif


int getSimuMaxIndex(LagSimuParams_t* simu);
fvec_t* getSimuRVec(LagSimuParams_t* simu);
fvec_t* getSimuXVec(LagSimuParams_t* simu);
fvec_t* getSimuVVec(LagSimuParams_t* simu);
myfloat getSimuObjX(LagSimuParams_t* simu);


/* * This tells the C++ compiler to expose this function with C linkage.
 * The C compiler just sees a standard function declaration.
 */
void initCleanLagSimu(int argc, char** args, LagSimuParams_t* simu);
void mainLoopCleanLagSimu(LagSimuParams_t* simu, uint32_t nb_frame_simul);
LagSimuParams_t* allocateSimuParams();

void checkArguments(int argc, char** args);
void createOutFolder(const char* folder_name);

#ifndef M_PI
    #define M_PI 3.141592653589793238462643383279502984
#endif

#ifndef M_1_PI
    #define M_1_PI           0.31830988618379067154  /* 1/pi */
#endif

// Define them as static inline so the compiler knows exactly what to inject
static inline double diameterToVolume(double dia){
    return 4./3.*M_PI*pow(dia/2., 3);
}

static inline double radiusToVolume(double r){
    return 4./3.*M_PI*MYPOW(r, 3);
}

static inline double volumeToRadius(double V){
    return MYCBRT(V*M_1_PI*0.75);
}

static inline myfloat momentumConservation(myfloat m1, myfloat m2, myfloat s1, myfloat s2){
    return (m1*s1 + m2*s2)/(m1+m2); 
}
// inline myfloat momentumConservation(myfloat m1, myfloat m2, myfloat s1, myfloat s2);
// inline double diameterToVolume(double dia);
// inline double radiusToVolume(double r);
// inline double volumeToRadius(double V);
// LagSimuParams_t* allocateSimuParams() {
//     return new LagSimuParams_t(); 
// }

#ifdef __cplusplus
}
#endif


#endif /* SIMU_HPP */
