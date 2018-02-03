//Copyright: 2015, Samuel Thompson, thompsonsed@gmail.com
// License: BSD-3
//All rights reserved.
//
//Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
//
//1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
//
//2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
//
//3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

//THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

//
/**
 * @file main.cpp
 * @brief A generic simulator for spatially explicit coalescence models suitable for HPC applications.
 * It contains all functions for running large-scale simulations backwards in time using coalescence techniques.
 * Outputs include an SQLite database containing spatial and temporal information about tracked lineages, and allow for rebuilding of the coalescence tree.
 * Currently, a fat-tailed dispersal kernel or normal distribution can be used for dispersal processes.
 *
 * Run with -h to see full input options.
 *
 * Outputs include
 * - habitat map file(s)
 * - species richness and species abundances for the supplied minimum speciation rate.
 * - SQL database containing full spatial data. This can be later analysed by the Speciation_Counter program for applying higher speciation rates.
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
 * @copyright <a href="https://opensource.org/licenses/BSD-3-Clause">BSD-3 Licence.</a>
 *
*/



#ifdef DEBUG
#ifndef verbose
#define verbose
#endif
#endif

#include "SpatialTree.h"
#include "SimulationTemplates.h"

// #define pristine_mode // not required unless you experience problems.
// This performs a more thorough check after each move operation.
// Currently, it will also check that the pristine state value is greater than the returned value within every map cell.
// Note that this may cause problems if the pristine state is not the state with the highest number of individuals.



/************************************************************
		MAIN ROUTINE AND COMMAND LINE ARG ROUTINES

 ************************************************************/



/**
 * @brief Main function containing program structure
 * @param argc the number of command-line arguments provided
 * @param argv a pointer to the arguments
 * @return a program exit code, 0 if successful, -1 (generally) indicates an error.
 */
int main(int argc, char *argv[])
{
	vector<string> comargs;
	importArgs(static_cast<const unsigned int &>(argc), argv, comargs);
	const string &config_file = getConfigFileFromCmdArgs(comargs);
	runMain<SpatialTree>(config_file);
	return 0;
}
