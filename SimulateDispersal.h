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
#include <string>
#include <stdio.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <stdexcept>
#include <sqlite3.h>
#include "Matrix.h"
#include "NRrand.h"
/**
 * @class Cell
 * @brief Simple structure containing the x and y positions of a cell
 */
struct Cell
{
	long x;
	long y;
	/**
	 * @brief Overloading equality operator
	 * @param c the Cell containing the values to overload
	 * @return the cell with the new values
	 */
	Cell &operator=(Cell const& c)
	= default;
};

/**
 * @brief Calculates the distance between two cells
 *
 * @param c1 Cell containing one point
 * @param c2 Cell containing second point
 * @return the distance between the two points
 */
double distanceBetween(Cell &c1, Cell &c2);

/**
 * @class SimulateDispersal
 * @brief Contains routines for importing a density map file, running a dispersal kernel n times on a landscape and record the 
 * dispersal distances.
 */
class SimulateDispersal
{
protected:
	// The density map object
	Matrix<uint32_t> density_map;
	// Set to true when the size of the density map has been set
	bool has_set_size;
	// The random number generator object
	NRrand random;
	// The map file path
	string map_name;
	// The random number seed
	unsigned long seed;
	// The dispersal method
	string dispersal_method;
	// The dispersal sigma value
	double sigma;
	// The dispersal nu value (for fat-tailed dispersal kernels)
	double tau;
	// The dispersal m_probability - chance of picking from a uniform distribution (for norm-uniform dispersal kernels)
	double m_prob;
	// The maximum dispersal distance for the norm-uniform dispersal distance
	double cutoff;
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
	// The maximal density value
	unsigned long max_density;
	// If true, sequentially selects dispersal probabilities, default is true
	bool is_sequential;
	// Reference number for this set of parameters in the database output
	unsigned long parameter_reference;
	// Function pointer for the landscape function
	typedef bool (SimulateDispersal::*landscape_fptr)(const double &dist, const double &angle,
													  const Cell &this_cell, Cell &end_cell);
	landscape_fptr getValFptr;
public:
	SimulateDispersal()
	{
		has_set_size = false;
		sigma = 0.0;
		tau = 0.0;
		m_prob = 0.0;
		cutoff = 0.0;
		num_repeats = 0;
		num_steps = 0;
		database = nullptr;
		max_density = 0;
		seed = 0;
		is_sequential = false;
		parameter_reference = 0;
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
	 * @brief Sets the sizes of the density map
	 * @param x the x dimension (number of columns) in the density map
	 * @param y the y dimension (number of rows) in the density map
	 */
	void setSizes(unsigned long x, unsigned long y);

	/**
	 * @brief Import the map from the prescribed file.
	 * @param map_file the map file to import from
	 */
	void importMaps(string map_file);
	
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
	 * @brief Sets the dispersal parameters
	 * @param dispersal_method_in the dispersal method (e.g. "normal")
	 * @param sigma_in the sigma value for normal and fat-tailed dispersals
	 * @param tau_in the nu value for fat-tailed dispersals
	 * @param m_prob_in the m_prob for norm-uniform dispersals
	 * @param cutoff_in the maximum dispersal distance for norm-uniform dispersal
	 * @param landscape_type string containing the landscape type (one of "closed", "tiled" or "infinite").
	 */
	void setDispersalParameters(string dispersal_method_in, double sigma_in, double tau_in, double m_prob_in,
								 double cutoff_in, string landscape_type);

	void setLandscapeType(string landscape_type);

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
	 * @brief Calculates the new position from the start cell based on the distance and angle moved.
	 * Stores the new x and y location in the end cell.
	 * @param dist the distance to move
	 * @param angle the direction from the start cell to move
	 * @param start_cell the cell containing the start x and y position
	 * @param end_cell the cell to contain the end x and y position
	 */
	void calculateNewPosition(const double &dist, const double &angle, const Cell &start_cell, Cell &end_cell);

	/**
	 * @brief Checks the density is greater than 0 a given distance from the start point on an infinite null landscape.
	 *
	 * This also takes into account the rejection sampling of density based on the maximal density value from the map.
	 *
	 * @param dist the distance of dispersal
	 * @param angle the angle of dispersal
	 * @param this_cell Cell containing the x and y coordinates of the starting position
	 * @param end_cell Cell to store the end point in
	 * @return true if the point has a density > 0, otherwise generates a random number between 0 and the max size and
	 * if it is > the density, return true, false otherwise.
	 */
	bool getEndPointInfinite(const double &dist, const double &angle, const Cell &this_cell, Cell &end_cell);

	/**
	 * @brief Checks the density is greater than 0 a given distance from the start point on an infinite tiled landscape.
	 *
	 * This also takes into account the rejection sampling of density based on the maximal density value from the map.
	 *
	 * @param dist the distance of dispersal
	 * @param angle the angle of dispersal
	 * @param this_cell Cell containing the x and y coordinates of the starting position
	 * @param end_cell Cell to store the end point in

	 * @return true if the point has a density > 0, otherwise generates a random number between 0 and the max size and
	 * if it is > the density, return true, false otherwise.
	 */
	bool getEndPointTiled(const double &dist, const double &angle, const Cell &this_cell, Cell &end_cell);

	/**
	 * @brief Checks the density a given distance from the start point on an closed landscape.
	 *
	 * This also takes into account the rejection sampling of density based on the maximal density value from the map.
	 *
	 * @param dist the distance of dispersal
	 * @param angle the angle of dispersal
	 * @param this_cell Cell containing the x and y coordinates of the starting position
	 * @param end_cell Cell to store the end point in

	 * @return true if the point has a density > 0, otherwise generates a random number between 0 and the max size and
	 * if it is > the density, return true, false otherwise.
	 */
	bool getEndPointClosed(const double &dist, const double &angle, const Cell &this_cell, Cell &end_cell);

	/**
	 * @brief Checks the density a given distance from the start point, calling the relevant landscape function.
	 * 
	 * This also takes into account the rejection sampling of density based on the maximal density value from the map.
	 * 
	 * @param dist the distance of dispersal
	 * @param angle the angle of dispersal
	 * @param this_cell Cell containing the x and y coordinates of the starting position
	 * @param end_cell Cell to store the end point in

	 * @return true if the point has a density > 0, otherwise generates a random number between 0 and the max size and
	 * if it is > the density, return true, false otherwise.
	 */
	bool getEndPoint(const double &dist, const double &angle, const Cell &this_cell, Cell &end_cell);
	
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