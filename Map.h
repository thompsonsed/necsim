//This file is part of NECSim project which is released under BSD-3 license.
//See file **LICENSE.txt** or visit https://opensource.org/licenses/BSD-3-Clause) for full license details
/**
 * @author Samuel Thompson
 * @date 31/08/16
 * @file Map.h
 *
 * @brief Contains the Map object for easy referencing of the respective coarse and fine map within the same coordinate system.
 * @copyright <a href="https://opensource.org/licenses/BSD-3-Clause">BSD-3 Licence.</a>
 */

#ifndef MAP
#define MAP

# include <string>
# include <stdio.h>
#include <vector>
# include <iostream>
# include <fstream>
# include <math.h>
# include <stdexcept>
# include <boost/filesystem.hpp>

#include "Matrix.h"
#include "DataMask.h"
#include "SimParameters.h"

using namespace std;


/************************************************************
					MAP OBJECT
 ************************************************************/
/**
 * @class Map
 * @brief Contains all maps and provides the functions for accessing a grid cell in the correct temporal and spacial location.
 * The function runDispersal() also provides the move routine, provided two alternative methods for moving individuals.
 */
// Object containing both the maps (the coarse and fine version) and routines for easy setting up and switching between the different coordinate systems.
class Map
{
protected:
	// The map files which are read in (or generated if running with "null" as the map file".
	// Pristine maps are meant for before any deforestation occured, whereas the other maps are intended for modern day maps.
	// A linear transformation from modern to pristine maps is used, approaching the habitat_change_rate variable times the difference between the pristine and modern maps.
	// Once the gen_since_pristine number of generations has been reached, the map will jump to the pristine condition.
	// the finer grid for the area around the sample area.
	Matrix<uint32_t> fine_map;
	// the pristine finer map.
	Matrix<uint32_t> pristine_fine_map;
	// the coarser grid for the wider zone.
	Matrix<uint32_t> coarse_map;
	// the pristine coarser map.
	Matrix<uint32_t> pristine_coarse_map;
	// for importing and storing the simulation set-up options.
	SimParameters mapvars;
	// the minimum values for each dimension for offsetting.
	long fine_x_min, fine_y_min, coarse_x_min, coarse_y_min;
	// the maximum values for each dimension for offsetting.
	long fine_x_max, fine_y_max, coarse_x_max, coarse_y_max;
	// the offsetting of the map in FINE map units.
	long fine_x_offset, fine_y_offset, coarse_x_offset, coarse_y_offset;
	// the scale of the coarse map compared with the smaller map.
	unsigned long scale{};
	// the length of the grid where the species start.
	long x_dim;
	// the height of the grid where the species start.
	long y_dim;
	unsigned long deme{};
	// for checking that the dimensions have been set before attempting to importSpatialParameters the maps.
	bool check_set_dim;
	// for setting the movement cost through forest.
	double dispersal_relative_cost{};
	// the last time the map was updated, in generations.
	double update_time;
	// the rate at which the habitat transforms from the modern forest map to the pristine habitat map.
	// A value of 1 will give a smooth curve from the present day to pristine habitat.
	double habitat_change_rate;
	// the number of generations at which point the habitat becomes entirely pristine.
	double gen_since_pristine;
	// the time the current map was updated.
	double current_map_time;
	// checks whether the simulation has already been set to the pristine state.
	bool is_pristine;
	// flag of whether the simulation has a pristine state or not.
	bool has_pristine;
	// the maximum value for habitat
	unsigned long habitat_max;
	// the maximum value on the fine map file
	unsigned long fine_max;
	// the maximum value on the coarse map file
	unsigned long coarse_max;
	// the maximum value on the pristine fine map file
	unsigned long pristine_fine_max;
	// the maximum value on the pristine coarse map file
	unsigned long pristine_coarse_max;
	// if true, dispersal is possible from anywhere, only the fine map spatial structure is preserved
	string landscape_type;
	string NextMap;
	// If this is false, there is no coarse map defined, so ignore the boundaries.
	bool bCoarse;
	// the number of updates to have occured.
	unsigned int nUpdate{};
	// Typedef for single application of the infinite landscape verses bounded landscape.
	typedef unsigned long (Map::*fptr)(const double &x, const double &y, const long &xwrap, const long &ywrap,
									   const double &dCurrentGen);

	fptr getValFunc;
public:
	/**
	 * @brief The default constructor.
	 */
	Map()
	{
		check_set_dim = false; // sets the check to false.
		is_pristine = false;
		current_map_time = 0;
		habitat_max = 1;
		getValFunc = nullptr;
		bCoarse = false;
		has_pristine = false;
		landscape_type = "closed";
		fine_max = 0;
		coarse_max = 0;
		pristine_fine_max = 0;
		pristine_coarse_max = 0;
	}

	/**
	 * @brief Gets the maximum habitat value from any map
	 * @return the maximum habitat value
	 */
	unsigned long getHabitatMax();

	/**
	 * @brief Sets the dimensions of the grid, the area where the species are initially sampled from.
	 * This function must be run before any of the calc map functions to allow for the correct deme allocation.
	 * 
	 * @param mapvarsin the SimParameters object containing the map variables to import
	 */
	void setDims(SimParameters mapvarsin);

	bool checkMapExists();

	/********************************************
	 * CALC MAP FUNCTIONS
	 ********************************************/

	/**
	 * @brief Imports the fine map object from file and calculates the correct values at each point.
	 * Without a map to input, the fine map will simply be a matrix of 1s.
	 */
	void calcFineMap();

	/**
	  * @brief Imports the pristine fine map object from file and calculates the correct values at each point.
	  * Without a map to input, the pristine fine map will simply be a matrix of 1s.
	  * This has the potential to be changed easily in future versions.
	  */
	void calcPristineFineMap();

	/**
	  * @brief Imports the coarse map object from file and calculates the correct values at each point.
	  * Without a map to input, the coarse map will simply be a matrix of 1s.
	  * This has the potential to be changed easily in future versions.
	  */
	void calcCoarseMap();

	/**
	  * @brief Imports the pristine coarse map object from file and calculates the correct values at each point.
	  * Without a map to input, the pristine coarse map will simply be a matrix of 1s.
	  * This has the potential to be changed easily in future versions.
	  */
	void calcPristineCoarseMap();

	/**
	 * @brief Sets the time variables.
	 * @param gen_since_pristine_in the time (in generations) since a pristine habitat state was achieved.
	 * @param habitat_change_rate_in the rate of transform of the habitat up until the pristine time.
	 * A value of 0.2 would mean 20% of the change occurs linearlly up until the pristine time and the remaining 80%
	 * occurs in a jump to the pristine state.
	 */
	void setTimeVars(double gen_since_pristine_in, double habitat_change_rate_in);


	/**
	 * @brief Calculates the offset and extremeties of the fine map.
	 * 
	 * Note that setting dispersal_relative_cost to a value other than 1 can massively increase simulation time.
	 * 
	 */
	void calcOffset();
	/********************************************
	 * VALIDATE MAPS
	 ********************************************/

	/**
	 * @brief Checks that the map file sizes are correct and that each value on the fragmented maps is less than the pristine maps.
	 * This should be disabled in simulations where habitat sizes are expected to shrink as well as grow.
	 */
	void validateMaps();

	/********************************************
	 * CHANGE MAP FUNCTIONS
	 ********************************************/

	/**
	 * @brief Updates the maps to the newer map.
	 */
	void updateMap(double generation);

	/**
	 * @brief Gets the pristine boolean.
	 * @return the pristine map state.
	 */
	bool isPristine()
	{
		if(has_pristine)
		{
			return is_pristine;
		}
		return true;
	}

	/**
	 * @brief Sets the pristine state of the system.
	 * @param bPristinein the pristine state.
	 */
	void setPristine(const bool &bPristinein)
	{
		is_pristine = bPristinein;
	}

	/**
	 * @brief Get the pristine map time
	 * @return double the pristine map time
	 */
	double getPristine()
	{
		return gen_since_pristine;
	}

	string getLandscapeType()
	{
		return landscape_type;
	}

	/**
	 * @brief Checks if the pristine state has been reached.
	 * 
	 * If there are no pristine maps, this function will do nothing.
	 * @param generation the time to check at.
	 */
	void checkPristine(double generation)
	{
		if(has_pristine)
		{
			if(generation >= gen_since_pristine)
			{
				is_pristine = true;
			}
		}
	}
	// 
	/********************************************
	 * GET VAL FUNCTIONS
	 ********************************************/


	// Function for getting the val at a particular coordinate from either the coarse or fine map
	// altered to use the current generation as well to determine the value.


	/**
	 * @brief Sets the landscape functions to either infinite or finite
	 * @param is_infinite a string of either closed, infinite, tiled_fine or tiled_coarse,
	 * corresponding to the relevant landscape type.
	 * 
	 * 
	 */
	void setLandscape(string is_infinite);

	/**
	 * @brief Gets the value at a particular coordinate from the correct map.
	 * Takes in to account temporal and spatial referencing.
	 * This version involves a call to the function pointer, *getValFunc, so that the correct call to either
	 * getValFinite() or getValInfinite is made.
	 * @param x the x position on the grid.
	 * @param y the y position on the grid.
	 * @param xwrap the number of wraps in the x dimension..
	 * @param ywrap the number of wraps in the y dimension..
	 * @param current_generation the current generation time.
	 * @return the value on the correct map at the correct space.
	 */
	unsigned long getVal(const double &x, const double &y,
						 const long &xwrap, const long &ywrap, const double &current_generation);

	/**
	 * @brief Gets the value from the coarse maps, including linear interpolating between the pristine and present maps
	 * @param xval the x coordinate
	 * @param yval the y coordinate
	 * @param current_generation the current generation timer
	 * @return the value of the map at the given coordinates and time
	 */
	unsigned long getValCoarse(const double &xval, const double &yval, const double &current_generation);

	/**
	 * @brief Gets the value from the fine maps, including linear interpolating between the pristine and present maps
	 * @param xval the x coordinate
	 * @param yval the y coordinate
	 * @param current_generation the current generation timer
	 * @return the value of the map at the given coordinates and time
	 */
	unsigned long getValFine(const double &xval, const double &yval, const double &current_generation);

	/**
	 * @brief Gets the value at a particular coordinate from the correct map.
	 * Takes in to account temporal and spatial referencing. This version assumes finite landscape.
	 * @param x the x position on the grid.
	 * @param y the y position on the grid.
	 * @param xwrap the number of wraps in the x dimension..
	 * @param ywrap the number of wraps in the y dimension..
	 * @param current_generation the current generation time.
	 * @return the value on the correct map at the correct space.
	 */
	unsigned long
	getValFinite(const double &x, const double &y, const long &xwrap, const long &ywrap, const double &current_generation);


	/**
	 * @brief Gets the value at a particular coordinate from the correct map.
	 * Takes in to account temporal and spatial referencing. This version assumes an infinite landscape.
	 * @param x the x position on the grid.
	 * @param y the y position on the grid.
	 * @param xwrap the number of wraps in the x dimension..
	 * @param ywrap the number of wraps in the y dimension..
	 * @param current_generation the current generation time.
	 * @return the value on the correct map at the correct space.
	 */
	unsigned long
	getValInfinite(const double &x, const double &y, const long &xwrap, const long &ywrap, const double &current_generation);


	/**
	 * @brief Gets the value at a particular coordinate from the correct map.
	 * Takes in to account temporal and spatial referencing.
	 * This version assumes an infinite landscape of tiled coarse maps.
	 * @param x the x position on the grid.
	 * @param y the y position on the grid.
	 * @param xwrap the number of wraps in the x dimension..
	 * @param ywrap the number of wraps in the y dimension..
	 * @param current_generation the current generation time.
	 * @return the value on the correct map at the correct space.
	 */
	unsigned long getValCoarseTiled(const double &x, const double &y, const long &xwrap, const long &ywrap,
									const double &current_generation);


	/**
	 * @brief Gets the value at a particular coordinate from the correct map.
	 * Takes in to account temporal and spatial referencing. 
	 * This version assumes an infinite landscape of tiled fine maps.
	 * @param x the x position on the grid.
	 * @param y the y position on the grid.
	 * @param xwrap the number of wraps in the x dimension..
	 * @param ywrap the number of wraps in the y dimension..
	 * @param current_generation the current generation time.
	 * @return the value on the correct map at the correct space.
	 */
	unsigned long getValFineTiled(const double &x, const double &y, const long &xwrap, const long &ywrap,
								  const double &current_generation);

	/**
	 * @brief Gets the x position on the fine map, given an x and x wrapping.
	 * 
	 * Note that this function will not check if the value is actually within bounds of the fine map, 
	 * and an error will likely be thrown by the matrix referencing if this is the case.
	 * @param x the x coordinate on the sample mask
	 * @param xwrap the x wrapping of the sample mask.
	 * @return the x location on the fine map
	 */
	unsigned long convertSampleXToFineX(const unsigned long &x, const long &xwrap);

	/**
	 * @brief Gets the y position on the fine map, given a y and y wrapping.
	 * 
	 * Note that this function will not check if the value is actually within bounds of the fine map, 
	 * and an error will likely be thrown by the matrix referencing if this is the case.
	 * @param y the y coordinate on the sample mask
	 * @param ywrap the y wrapping of the sample mask.
	 * @return the y location on the fine map
	 */
	unsigned long convertSampleYToFineY(const unsigned long &y, const long &ywrap);

	/**
	 * @brief Converts the fine map coordinates to the sample grid coordinates.
	 * Main conversion is in a call to convertCoordinates, but also makes sure the returned types are long integers.
	 * @param x the x coordinate to modify
	 * @param xwrap the x wrapping to modify
	 * @param y the y coordinate to modify
	 * @param ywrap the y wrapping to modify
	 */
	void convertFineToSample(long &x, long &xwrap, long &y, long &ywrap);


	/**
	 * @brief Counts the number of spaces available in the initial species space. Requires the samplemask to check the sampling area.
	 * @param dSample the sample proportion (from 0 to 1).
	 * @param samplemask the DataMask object to sample from.
	 * @return the total number of individuals predicted to initially exist on the map.
	 */
	unsigned long getInitialCount(double dSample, DataMask &samplemask);

	/**
	 * @brief Gets the mapvars object for referencing simulation parameters.
	 * @return 
	 */
	SimParameters getSimParameters();

	/********************************************
	 * CHECK MAP FUNCTIONS
	 ********************************************/
	/**
	 * @brief Checks whether the point is habitat or non-habitat.
	 *@param x the x position on the grid.
	 * @param y the y position on the grid.
	 * @param xwrap the number of wraps in the x dimension.
	 * @param ywrap the number of wraps in the y dimension.
	 * @param generation the current generation time.
	 * @return a boolean of whether the map is habitat or non-habitat.
	 */
	bool checkMap(const double &x, const double &y, const long &xwrap, const long &ywrap, const double generation);


	/**
	 * @brief  Checks whether the point comes from the fine grid.
	 * @param x the x position.
	 * @param y the y position.
	 * @param xwrap the number of wraps in the x dimension.
	 * @param ywrap the number of wraps in the y dimension.
	 * @return a boolean of whether the location is on the fine map.
	 */
	bool checkFine(const double &x, const double &y, const long &xwrap, const long &ywrap);


	/**
	 * @brief Converts the coordinates to within the original grid, altering the xwrap and ywrap consequently.
	 * @param x the x position.
	 * @param y the y position.
	 * @param xwrap the number of wraps in the x dimension.
	 * @param ywrap the number of wraps in the y dimension.
	 */
	void convertCoordinates(double &x, double &y, long &xwrap, long &ywrap);

	/********************************************
	 * MAIN DISPERSAL FUNCTION
	 ********************************************/
	/**
	 * @brief The function that actually performs the dispersal. 
	 * It is included here for easier  programming and efficiency as the function doesn't need to perform all the checks until the edge of the fine grid.
	 * @param dist the distance travelled (or "distance energy" if dispersal_relative_cost is not 1).
	 * @param angle the angle of movement.
	 * @param startx the start x position.
	 * @param starty the start y position.
	 * @param startxwrap the start number of wraps in the x dimension.
	 * @param startywrap the start number of wraps in the y dimension.
	 * @param disp_comp a boolean of whether the dispersal was complete or not. This value is returned true if dispersal is to habitat, false otherwise.
	 * @param generation the time in generations since the start of the simulation.
	 * @return the density value at the end dispersal point
	 */
	unsigned long runDispersal(const double &dist, const double &angle, long &startx, long &starty, long &startxwrap,
							   long &startywrap, bool &disp_comp, const double &generation);


	/**
	 * @brief Operator for outputting the Map object variables to an output stream.
	 * This is used for storing the Map object to file.
	 * @param os the output stream.
	 * @param r the Map object to output.
	 * @return the output stream.
	 */
	friend ostream &operator<<(ostream &os, const Map &r)
	{
		os << r.mapvars << "\n" << r.fine_x_min << "\n" << r.fine_x_max << "\n" << r.coarse_x_min << "\n" << r.coarse_x_max;
		os << "\n" << r.fine_y_min << "\n" << r.fine_y_max << "\n" << r.coarse_y_min << "\n" << r.coarse_y_max << "\n";
		os << r.fine_x_offset << "\n" << r.fine_y_offset << "\n" << r.coarse_x_offset << "\n" << r.coarse_y_offset << "\n";
		os << r.scale << "\n" << r.x_dim << "\n" << r.y_dim << "\n" << r.deme << "\n" << r.check_set_dim << "\n"
		   << r.dispersal_relative_cost << "\n";
		os << r.update_time << "\n" << r.habitat_change_rate << "\n" << r.gen_since_pristine << "\n" << r.current_map_time << "\n"
		   << r.is_pristine << "\n";
		os << r.NextMap << "\n" << r.nUpdate << "\n" << r.landscape_type << "\n" << r.fine_max << "\n"
		   << r.coarse_max << "\n";
		os << r.pristine_fine_max << "\n" << r.pristine_coarse_max << "\n" << r.habitat_max << "\n"
		   << r.bCoarse << "\n" << r.has_pristine << "\n";
		return os;
	}

	/**
	 * @brief Operator for inputting the Map object variables from an input stream.
	 * This is used for reading the Map object from file.
	 * @param is the input stream.
	 * @param r the Map object to input to.
	 * @return the input stream.
	 */
	friend istream &operator>>(istream &is, Map &r)
	{
		is >> r.mapvars >> r.fine_x_min;
		is >> r.fine_x_max >> r.coarse_x_min;
		is >> r.coarse_x_max >> r.fine_y_min >> r.fine_y_max;
		is >> r.coarse_y_min >> r.coarse_y_max;
		is >> r.fine_x_offset >> r.fine_y_offset >> r.coarse_x_offset >> r.coarse_y_offset >> r.scale >> r.x_dim >> r.y_dim
		   >> r.deme >> r.check_set_dim >> r.dispersal_relative_cost;
		is >> r.update_time >> r.habitat_change_rate >> r.gen_since_pristine >> r.current_map_time >> r.is_pristine;
		getline(is, r.NextMap);
		is >> r.nUpdate;
		is >> r.landscape_type;
		is >> r.fine_max >> r.coarse_max;
		is >> r.pristine_fine_max >> r.pristine_coarse_max;
		is >> r.habitat_max >> r.bCoarse >> r.has_pristine;
		r.setLandscape(r.mapvars.landscape_type);
		r.calcFineMap();
		r.calcCoarseMap();
		r.calcPristineFineMap();
		r.calcPristineCoarseMap();
		r.recalculateHabitatMax();
		return is;
	}

	/**
	 * @brief Prints some selected Map variables to the terminal.
	 * @return the string containing the map variables to print
	 */
	string printVars();

	/**
	 * @brief Wipes the map of all variables. Only really useful for testing purposes.
	 */
	void clearMap();

	/**
	 * @brief Recalculates the habitat map maximum by checking the maximums for each of the relevant map files
	 * (fine, coarse and pristines).
	 */
	void recalculateHabitatMax();

};

#endif // MAP
 