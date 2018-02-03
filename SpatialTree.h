// This file is part of NECSim project which is released under BSD-3 license.
// See file **LICENSE.txt** or visit https://opensource.org/licenses/BSD-3-Clause) for full license details.
//
/**
 * @author Sam Thompson
 * @date 31/08/2016
 * @file SpatialTree.h
 * @brief Contains the SpatialTree class for running simulations and outputting the phylogenetic tree.
 *
 * Contact: samuel.thompson14@imperial.ac.uk or thompsonsed@gmail.com
 * @copyright <a href="https://opensource.org/licenses/BSD-3-Clause">BSD-3 Licence.</a>
 *
 */
#ifndef SPATIALTREE_H
#define SPATIALTREE_H

#include <cstdio>
#include <fstream>
#include <vector>
#include <iostream>
#include <string>
#include <cstring>
#include <cmath>
#include <iomanip>
#include <cmath>
#include <ctime>
#include <ctime>
#include <sqlite3.h>
#include <unistd.h>
#include <algorithm>
#include <stdexcept>
//#define with_gdal
// extra boost include - this requires the installation of boost on the system
// note that this requires compilation with the -lboost_filesystem and -lboost_system linkers.
#include <boost/filesystem.hpp>

// include fast-csv-parser by Ben Strasser (available from https://github.com/ben-strasser/fast-cpp-csv-parser)
// for fast file reading
#ifdef use_csv
#include "fast-cpp-csv-parser/csv.h"
#endif
/**
* @defgroup DEFINES Preprocessor Definitions
*/
/**
*@ingroup DEFINES
*
* Macro for using the fast-cpp-csv-parser from Ben Strasser (available from
*https://github.com/ben-strasser/fast-cpp-csv-parser).
* This enables much faster csv reading, but can cause problems on systems where this module is not fully tested.
*/

/**
*@ingroup DEFINES
* Macro to compile using RAM for storage of the active SQL database.
* Without this, the database will be written directly to disc (which is slower, but an option if RAM requirements get
*too huge).
* For HPC systems, it is recommended to use this option as write speeds are generally fast and large simulations don't
* have a linear increase in the SQL database size (at least in RAM).
*/
#ifndef sql_ram
#define sql_ram
#endif

// other includes for required files
#include "Tree.h"
#include "Matrix.h"
#include "NRrand.h"
#include "SimParameters.h"
#include "DataPoint.h"
#include "TreeNode.h"
#include "SpeciesList.h"
#include "Map.h"
#include "Community.h"
#include "Setup.h"
#include "DispersalCoordinator.h"
#include "ReproductionMap.h"
#include "Logging.h"

using namespace std;

/**
* @class SpatialTree
* @brief Represents the output phylogenetic tree, when run on a spatially-explicit landscape.
*
* Contains all functions for running simulations, outputting data and calculating coalescence tree structure.
*/
class SpatialTree : public virtual Tree
{
protected:
	// Our dispersal coordinator for getting dispersal distances and managing calls from the habitat_map
	DispersalCoordinator dispersal_coordinator;
	// The reproduction map object
	ReproductionMap rep_map;
	// A list of new variables which will contain the relevant information for maps and grids.
	//  strings containing the file names to be imported.
	string fine_map_input, coarse_map_input;
	string pristine_fine_map_input, pristine_coarse_map_input;
	// the time since pristine forest and the rate of change of the rainforest.
	double gen_since_pristine, habitat_change_rate;
	// the variables for the grid containing the initial individuals.
	unsigned long grid_x_size, grid_y_size;
	// The fine map variables at the same resolution as the grid.
	// the coarse map variables at a scaled resolution of the fine map.
	long fine_map_x_size, fine_map_y_size, fine_map_x_offset, fine_map_y_offset;
	long coarse_map_x_size, coarse_map_y_size, coarse_map_x_offset, coarse_map_y_offset, coarse_map_scale;
	// Map object containing both the coarse and fine maps for checking whether or not there is forest at a particular
	// location.
	Map habitat_map;
	// An indexing spatial positioning of the lineages
	Matrix<SpeciesList> grid;
	// dispersal and sigma references
	double sigma, tau;
	// the cost for moving through non-forest. 1.0 means there is no cost. 10 means that movement is 10x
	// slower through forest.
	double dispersal_relative_cost;
	// the desired number of species we are aiming for. If it is 0, we will carry on forever.
	unsigned long desired_specnum;
	// contains the DataMask for where we should start lineages from.
	DataMask samplegrid;
public:
	/**
	 * @brief The constructor for the tree object.
	 *
	 * Sets all uninitiated variables to false, except log_all.
	 * log_all should be changed to false if minimal text output during simulations is desired.
	 */
	SpatialTree() : Tree()
	{
		sigma = 0.0;
		tau = 0.0;
		outdatabase = nullptr;
		gen_since_pristine = 0.0;
		habitat_change_rate = 0.0;
		grid_x_size = 0;
		grid_y_size = 0;
		fine_map_x_size = 0;
		fine_map_y_size = 0;
		fine_map_x_offset = 0;
		fine_map_y_offset = 0;
		coarse_map_x_size = 0;
		coarse_map_y_size = 0;
		coarse_map_x_offset = 0;
		coarse_map_y_offset = 0;
		coarse_map_scale = 1;
		dispersal_relative_cost = 1.0;
		desired_specnum = 1;
	}

	~SpatialTree() = default;
	/**
	 * @brief Imports the simulation variables from the config file.
	 *
	 * Also checks for paused simulations and file existence.
	 * @param configfile a vector containing the configfile to import from
	 */
	void importSimulationVariables(const string &configfile) override;

	/**
	 * @brief Parses the command line arguments and saves the flags in Tree for fullmode,
	 * resuming and other important variables.
	 * @deprecated command line arguments are no longer supported. Function will be removed at next major update.
	 * @param comargs
	 */
	void parseArgs(vector<string> &comargs);

	/**
	 * @brief Checks that the folders exist and the files required for the simulation also exist.
	 */
	void checkFolders();

	/**
	 * @brief Sets the map object with the correct variables, taking the SimParameters structure defined elsewhere for the
	 * parameters.
	 * 
	 * Requires that parameters have already been imported into the SimParameters
	 *
	 * This function can only be run once, otherwise a Main_Exception will be thrown
	 *	 */
	void setParameters() override;


	// Imports the maps using the variables stored in the class. This function must be run after the set_mapvars() in
	// order to function correctly.
	/**
	 * @brief Imports the maps into the forestmap object.
	 *
	 * The simulation variables should have already been imported by setParameters(), otherwise a Fatal_Exception will be
	 * thrown.
	 */
	void importMaps();

	/**
	 * @brief Imports the reproduction map from file.
	 */
	void importReproductionMap();

	/**
	 * @brief Counts the number of individuals that exist on the spatial grid.
	 * @return the number of individuals that will be initially simulated
	 */
	unsigned long getInitialCount() override;

	/**
	 * @brief Sets up the dispersal coordinator by linking to the correct functions and choosing the appropriate dispersal
	 * method.
	 */
	void setupDispersalCoordinator();

	/**
	 * @brief Contains the setup routines for a spatial landscape.
	 * It also checks for paused simulations and imports data if necessary from paused files.
	 * importMaps() is called for importing the map files
	 *
	 * @bug For values of dispersal, forest transform rate and time since pristine (and any other double values),
	 * they will not be correctly outputted to the SIMULATION_PARAMETERS table if the value is smaller than 10e-6.
	 * The solution is to implement string output mechanisms using boost::lexical_cast(),
	 * but this has so far only been deemed necessary for the speciation rate (which is intrinsically very small).
	 *
	 */
	void setup() override;

	/**
	 * @brief Fill the active, data and grid objects with the starting lineages.
	 * @param initial_count the number of individuals expected to exist
	 * @return the number of lineages added (for validation purposes)
	 */
	unsigned long fillObjects(const unsigned long &initial_count) override;

	/**
	 * @brief Gets the number of individuals to be sampled at the particular point and time.
	 * Round the number down to the nearest whole number for numbers of individuals.
	 * @param x the x location for individuals to be sampled
	 * @param y the y location for individuals to be sampled
	 * @param x_wrap the number of x wraps for the cell
	 * @param y_wrap the number of y wraps for the cell
	 * @param current_gen the current generation timer
	 * @return the number of individuals to sample at this location.
	 */
	unsigned long getIndividualsSampled(const long &x, const long &y,
								 const long &x_wrap, const long &y_wrap, const double &current_gen);
	/**
	 * @brief Removes the old position within active by checking any wrapping and removing connections.
	 *
	 * The function also corrects the linked list to identify the correct nwrap for every wrapped lineage in that space.
	 *
	 * @param chosen the desired active reference to remove from the grid.
	 */
	void removeOldPosition(const unsigned long &chosen) override;

	/**
	 * @brief Calculate the move, given a start x,y coordinates and wrapping.
	 *
	 * The provided parameters will be altered to contain the new values so no record of the old variables remains after
	 * function running.
	 * Current dispersal methods use a fattailed dispersal.
	 */
	void calcMove();

	/**
	 * @brief  Calculates the minmax for a given branch.
	 *
	 * Calculates the speciation rate required for speciation to have occured on this branch.
	 *
	 * @param current the current active reference to perform calculations over.
	 */
	long double calcMinMax(const unsigned long &current);

	/**
	 * @brief Calculates the new position, checking whether coalescence has occured and with which lineage.
	 *
	 * This involves correct handling of checking wrapped lineages (outside the original grid). The probability of
	 * coalescence is also calculated.
	 *
	 * @param coal boolean for whether coalescence occured or not
	 * @param chosen the chosen lineage
	 * @param coalchosen the lineaged that is coalescing (if required)
	 * @param  oldx the old x position
	 * @param oldy the old y position
	 * @param oldxwrap the old x wrapping
	 * @param oldywrap the old y wrapping
	 *
	 */
	void calcNewPos(bool &coal,
					const unsigned long &chosen,
					unsigned long &coalchosen,
					const long &oldx,
					const long &oldy,
					const long &oldxwrap,
					const long &oldywrap);

	/**
	 * @brief Switches the chosen position with the endactive position.
	 * @param chosen the chosen lineage to switch with endactive.
	 */
	void switchPositions(const unsigned long &chosen) override;

	/**
	 * @brief Calculates the next step for the simulation.
	 */
	void calcNextStep() override;

	/**
	 * @brief Estimates the species number from the second largest minimum speciation rate remaining in active.
	 *
	 *  This allows for halting of the simulation once this threshold has been reached.
	 * However, the function is not currently in use as calculating the coalescence tree is very computionally
	 * intensive.
	 * @deprecated this function is currently obselete and not implemented, but was still functional as of version 3.2
	 *
	 */
	unsigned long estSpecnum();

#ifdef pristine_mode
	/**
	 * @brief Checks that pristine maps make sense. Only relevant for pristine mode.
	 */
	void pristineStepChecks();
#endif


	/**
	 * @brief Increments the generation counter and step references, then updates the map for any changes to habitat
	 * cover.
	 */
	void incrementGeneration() override ;

	/**
	 * @brief Updates the coalescence variables in the step object.
	 */
	void updateStepCoalescenceVariables() override;

	/**
	 * @brief Expands the map, generating the new lineages where necessary.
	 *
	 * The samplemask provided is used for expansion. Any empty spaces are filled with a new lineage.
	 * Lineages which have not moved are changed to tips, with a new data entry so that original and new generations are
	 * recorded.
	 *
	 * @param generation_in the generation that the expansion is occuring at. This is used in recording the new tips
	 */
	void addLineages(double generation_in) override;

	/**
	 * @brief Creates a string containing the SQL insertion statement for the simulation parameters.
	 * @return string containing the SQL insertion statement
	 */
	string simulationParametersSqlInsertion() override;


	/**
	 * @brief Pause the simulation and dump data from memory.
	 */
	void simPause() override;

	/**
	 * @brief Saves the map object to file.
	 * @param pause_folder the folder to save files into
	 */
	void dumpMap(string pause_folder);

	/**
	 * @brief Resumes the simulation from a previous state.
	 *
	 * Reads in the parameters and objects from file and re-starts the simulation.
	 */
	void simResume() override;

	/**
	 * @brief Loads the grid from the save file into memory.
	 */
	void loadGridSave();

	/**
	 * @brief Loads the map from the save file into memory.
	 */
	void loadMapSave();

	/**
	 * @brief Checks that the reproduction map makes sense with the fine density map.
	 */
	void verifyReproductionMap();

	/**
	 * @brief Adds the lineage to the correct point in the linked list of active lineages
	 * @param numstart the active position to add
	 * @param x the x position of the lineage
	 * @param y the y position of the lineage
	 */
	void addWrappedLineage(unsigned long numstart, long x, long y);

	/**
	 * @brief Counts the number of lineages at a particular location that need to be added, after making the correct
	 * number of those that already exist into tips.
	 *
	 * @param x the x coordinate of the location of interest
	 * @param y the y coordinate of the location of interest
	 * @param xwrap the x wrapping of the location
	 * @param ywrap the y wrapping of the location
	 * @param generationin the generation to assign to new tips
	 * @param make_tips if true, stores the tips as well as counting them
	 * @return the number of new lineages that need to be added.
	 */
	unsigned long countCellExpansion(const long &x, const long &y, const long &xwrap, const long &ywrap,
									 const double &generationin, const bool &make_tips);

	/**
	 * @brief Expands the cell at the desired location by adding the supplied number of lineages
	 *
	 * This takes into account wrapping to correctly add the right number
	 *
	 * @param x the x coordinate to add at
	 * @param y the y coordinate to add at
	 * @param x_wrap the x wrapping to add at
	 * @param y_wrap the y wrapping to add at
	 * @param generation_in the generation to set the new lineages to
	 * @param add the total number of lineages to add at this location
	 */
	void expandCell(long x, long y, long x_wrap, long y_wrap, double generation_in, unsigned long add);




#ifdef DEBUG
	/**
	 * @brief Validates all lineages have been set up correctly. This may take considerable time for larger simulations.
	 *
	 * This function is only relevant in debug mode.
	 */
	void validateLineages() override;

	/**
	 * @brief Runs the debug checks upon a dispersal event.
	 *
	 */
	void debugDispersal();
	/**
	 * @brief Checks that adding a lineage has resulting in the correct structures being created.
	 *
	 * This function is only relevant in debug mode.
	 *
	 * @param numstart the active lineage to check for
	 * @param x the x coordinate
	 * @param y the y coordinate
	 */
	void debugAddingLineage(unsigned long numstart, long x, long y);

	/**
	 * @brief Run checks at the end of each cycle which make certain the move has been successful.
	 * @param chosen the chosen lineage to check
	 * @param coalchosen the lineage which is coalescing with the chosen lineage which we are also required to check
	 */
	void runChecks(const unsigned long &chosen, const unsigned long &coalchosen) override;
#endif
};

#endif  // SPATIALTREE_H
