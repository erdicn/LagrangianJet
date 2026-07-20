#include <iostream>

#include "simu.hpp"

int main(int argc, char** args){
    LagSimuParams_t simu = {};
    initCleanLagSimu(argc, args, &simu);
    // return mainCleanLagSimu(argc, args);
    std::cout << "Entering main loop\n" << std::flush;
    simu.batch_start_time = std::chrono::high_resolution_clock::now();
    while(simu.current_saved < simu.nb_to_save)
        mainLoopCleanLagSimu(&simu, 1);

    return 0;
}