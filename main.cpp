// This file is part of necsim project which is released under MIT license.
// See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.
/**
 * @file main.cpp
 * @brief A generic simulator for spatially explicit coalescence models suitable for HPC applications.
 * It contains all functions for running large-scale simulations backwards in time using coalescence techniques.
 * Outputs include an SQLite database containing spatial and temporal information about tracked lineages, and allow for
 * rebuilding of the coalescence tree.
 *
 * Currently, a fat-tailed dispersal kernel or normal distribution can be used for dispersal processes.
 *
 * Run with -h to see full input options.
 *
 * Outputs include
 * - habitat map file(s)
 * - species richness and species abundances for the supplied minimum speciation rate.
 * - SQL database containing full spatial data. This can be later analysed by the Speciation_Counter program for
 *   applying higher speciation rates.
 *
 * Contact: samuel.thompson14@imperial.ac.uk or thompsonsed@gmail.com
 *
 * Based heavily on code written by James Rosindell
 *
 * Contact: j.rosindell@imperial.ac.uk
 *
 *
 * @author Samuel Thompson
 *
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 *
*/



#ifdef DEBUG
#ifndef verbose
#define verbose
#endif
#endif

#include "Logger.h"
#include "Logging.h"
#include "ProtractedSpatialTree.h"
#include "SimulationTemplates.h"
#include "GenericTree.h"
// #define historical_mode // not required unless you experience problems.
// This performs a more thorough check after each move operation.
// Currently, it will also check that the historical state value is greater than the returned value within every map cell.
// Note that this may cause problems if the historical state is not the state with the highest number of individuals.



/************************************************************
        MAIN ROUTINE AND COMMAND LINE ARG ROUTINES

 ************************************************************/



/**
 * @brief Main function containing program structure
 * @param argc the number of command-line arguments provided
 * @param argv a pointer to the arguments
 * @return a program exit code, 0 if successful, -1 (generally) indicates an error.
 */
int main(int argc, char* argv[])
{
//        logger = new Logger();
//        vector<string> comargs;
//        importArgs(static_cast<const unsigned int &>(argc), argv, comargs);
//        const string &config_file = getConfigFileFromCmdArgs(comargs);
//        runMain<SpatialTree>(config_file);
//        delete logger;
//        return 0;
//    necsim::ProtractedSpatialTree test2;
//    test3 = necsim::ProtractedSpatialTree();
//    necsim::ProtractedSpatialTree test = necsim::ProtractedSpatialTree();
//    test2 = test3;
//    test2 = std::move(test);
//    test1 = test3;
    necsim::ProtractedSpatialTree tree();
    necsim::logger = new necsim::Logger();
    int __pyx_lineno = 0;
    const char *__pyx_filename = NULL;
    int __pyx_clineno = 0;
//    try{
//        tree = necsim::ProtractedSpatialTree();
//    }
//    catch(...){
//        return 1;
//    }

    return 0;
}