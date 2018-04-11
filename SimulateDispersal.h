// This file is part of NECSim project which is released under BSD-3 license.
// See file **LICENSE.txt** or visit https://opensource.org/licenses/BSD-3-Clause) for full license details.

/**
 * @author Samuel Thompson
 * @file SimulateDispersal.h
 * @brief Contains the ability to simulate a given dispersal kernel on a specified density map, outputting the effect 
 * dispersal distance distribution to an SQL file after n number of dispersal events (specified by the user).
 * @copyright <a href="https://opensource.org/licenses/BSD-3-Clause">BSD-3 Licence.</a>
 */

#ifndef DISPERSAL_TEST
#define DISPERSAL_TEST
#ifndef PYTHON_COMPILE
#define PYTHON_COMPILE
#endif
#include<string>
#include <cstdio>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <stdexcept>
#include <sqlite3.h>
#include "Landscape.h"
#include "DispersalCoordinator.h"
#include "NRrand.h"
#include "Cell.h"

/**
 * @brief Calculates the distance between two cells
 *
 * @param c1 Cell containing one point
 * @param c2 Cell containing second point
 * @return the distance between the two points
 */
double distanceBetweenCells(Cell &c1, Cell &c2);

/**
 * @class SimulateDispersal
 * @brief Contains routines for importing a density map file, running a dispersal kernel n times on a landscape and record the 
 * dispersal distances.
 */
class SimulateDispersal
{
protected:
	// The density map object
	Landscape density_landscape;
	// Dispersal coordinator
	DispersalCoordinator dispersal_coordinator{};
	// Stores all key simulation parameters for the Landscape object
	SimParameters  * simParameters;
	// The random number generator object
	NRrand random;
	// The random number seed
	unsigned long seed;
	// The sqlite3 database object for storing outputs
	sqlite3 * database;
	// Vector for storing successful dispersal distances
	vector<double> distances;
	// Vector for storing the cells (for randomly choosing from)
	vector<Cell> cells;
	// The number of repeats to run the dispersal loop for
	unsigned long num_repeats;
	// The number of num_steps within each dispersal loop for the average distance travelled/
	unsigned long num_steps;
	// generation counter
	double generation;
	// If true, sequentially selects dispersal probabilities, default is true
	bool is_sequential;
	// Reference number for this set of parameters in the database output
	unsigned long parameter_reference;
public:
	SimulateDispersal()
	{
		simParameters = nullptr;
		num_repeats = 0;
		num_steps = 0;
		database = nullptr;
		seed = 0;
		is_sequential = false;
		parameter_reference = 0;
		generation = 0.0;
	}
	
	~SimulateDispersal()
	{
		sqlite3_close(database);
	}
	
	/**
	 * @brief Sets the is_sequential flag
	 * @param bSequential if true, dispersal events are selected using the end point of the last dispersal 
	 * distance for the start of the next move event
	 */
	void setSequential(bool bSequential);

	/**
	 * @brief Sets the pointer to the simulation parameters object
	 * @param sim_parameters pointer to the simulation parameters to use
	 */
	void setSimulationParameters(SimParameters * sim_parameters);

	/**
	 * @brief Import the maps from the prescribed files.
	 * @param fine_map_file path to the fine map file
	 * @param fine_map_x fine map x dimension
	 * @param fine_map_y fine map y dimension
	 * @param fine_map_x_offset fine map x offset
	 * @param fine_map_y_offset fine map y offset
	 * @param coarse_map_file file path to the coarse map file
	 * @param coarse_map_x coarse map x dimension
	 * @param coarse_map_y coarse map y dimension
	 * @param coarse_map_x_offset coarse map x offset
	 * @param coarse_map_y_offset coarse map y offset
	 */
	void importMaps();
	
	/**
	 * @brief Sets the seed for the random number generator
	 * @param s the seed 
	 */
	void setSeed(unsigned long s)
	{
		seed = s;
		random.setSeed(s);
	}

	/**
	 * @brief Sets the output database for writing results to
	 * @param out_database path to the output database
	 */
	void setOutputDatabase(string out_database);
	
	/**
	 * @brief Sets the number of repeats to run the dispersal kernel for
	 * @param n the number of repeats
	 */
	void setNumberRepeats(unsigned long n);

	/**
	 * @brief Sets the number of steps to run each repeat of the dispersal kernel for when recording mean distance
	 * travelled
	 * @param s the number of steps
	 */
	void setNumberSteps(unsigned long s);
	/**
	 * @brief Calculates the list of cells to choose randomly from 
	 */
	void storeCellList();
	
	/**
	 * @brief Gets a random cell from the list of cells
	 * @return a Cell object reference containing the x and y positions to choose from
	 */
	const Cell& getRandomCell();

	/**
	 * @brief Checks the density a given distance from the start point, calling the relevant landscape function.
	 * 
	 * This also takes into account the rejection sampling of density based on the maximal density value from the map.
	 * 
	 * @param this_cell Cell containing the x and y coordinates of the starting position
	 */
	void getEndPoint(Cell &this_cell);
	
	/**
	 * @brief Simulates the dispersal kernel for the set parameters, storing the mean dispersal distance
	 */
	void runMeanDispersalDistance();

	/**
	 * @brief Simulates the dispersal kernel for the set parameters, storing the mean distance travelled.
	 */
	void runMeanDistanceTravelled();
	
	/**
	 * @brief Writes out the distances to the SQL database.
	 * @param table_name the name of the table to output to, either 'DISPERSAL_DISTANCE' or 'DISTANCES_TRAVELLED'
	 */
	void writeDatabase(string table_name);

	/**
	 * @brief Writes the simulation parameters to the output SQL database.
	 * @param table_name the name of the table to output to, either 'DISPERSAL_DISTANCE' or 'DISTANCES_TRAVELLED'
	 */
	void writeParameters(string table_name);

	/**
	 * @brief Gets the maximum parameter reference from the output SQL database and saves val + 1 to parameter_reference
	 * Assumes that the database exists.
	 *
	 */
	void checkMaxParameterReference();

	/**
	 * @brief Gets the maximum id number from the output SQL database and returns val + 1
	 * Assumes that the database exists.
	 * @note this function does not check for SQL injection attacks and should not be used with variable function names.
	 * @param table_name: the name of the table to check for max(id) in
	 * @returns the maximum id + 1 from the given table
	 */
	unsigned long checkMaxIdNumber(string table_name);
};

#endif