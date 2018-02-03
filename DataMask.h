/**
 * @author Samuel Thompson
 * @file DataMask.h
 * @brief  Contains the DataMask class for describing the spatial sampling pattern on a landscape.
 *
 * @copyright <a href="https://opensource.org/licenses/BSD-3-Clause">BSD-3 Licence.</a>
 */
#ifndef SPECIATIONCOUNTER_DataMask_H
#define SPECIATIONCOUNTER_DataMask_H

/**
 * @brief Forward declaration for the Map class
 */
class Map;

#include <string>

#include "Matrix.h"
#include "SimParameters.h"


// Class which contains the DataMask object, telling us where to sample from within the habitat map.
/**
 * @class DataMask
 * @brief Contains the DataMask object, a Matrix of booleans describing the spatial sampling pattern.
 */
class DataMask
{
protected:
	// the file to read in from
	string inputfile;
	bool bDefault;
	unsigned long x_offset, y_offset;
	// Stores the size of the grid which is stored as a full species list
	unsigned long x_dim, y_dim;
	// Stores the size of the samplemask from which spatially sampling is read
	unsigned long mask_x_dim, mask_y_dim;
	// Function pointer for obtaining the proportional sampling from the sample mask.
	typedef double (DataMask::*fptr)(const long &x, const long &y, const long &xwrap, const long &ywrap);
	fptr getProportionfptr;
public:
	Matrix<bool> sample_mask; /** A binary grid telling whether or not the cell should be sampled.*/
	// Stores the exact values from the input tif file.
	Matrix<double> sample_mask_exact; /** Exact grid for determining sampling proportion.*/

	/**
	 * @brief The DataMask constructor.
	 */
	DataMask();

	~DataMask() = default;

	/**
	 * @brief Returns if the simulation is using the a null samplemask, and therefore does not need to store the full
	 * sample grid in memory.
	 * @return true if using a null samplemask
	 */
	bool getDefault();

	/**
	 * @brief Sets the parameters for the datamask, including the dimensions of the map, the offsets from the grid and
	 * the dimensions of the grid itself for recalculating coordinates.
	 * @param sample_mask_file the file to import the samplemask from (or "null")
	 * @param x_in x dimension of the grid
	 * @param y_in y dimension of the grid
	 * @param mask_x_in x dimension of the sample mask
	 * @param mask_y_in y dimension of the sample mask
	 * @param x_offset_in x offset of the sample mask from the grid
	 * @param y_offset_in y offset of the sample mask from the grid
	 * @return true if using a "null" samplemask
	 */
	bool setup(const string &sample_mask_file, const unsigned long &x_in, const unsigned long &y_in,
			   const unsigned long &mask_x_in, const unsigned long &mask_y_in,
			   const unsigned long &x_offset_in, const unsigned long &y_offset_in);

	/**
	 * @brief Imports the sample mask as a boolean mask and sets the relevant sample mask dimensions.
	 * @param xdim the x dimension of the grid area
	 * @param ydim the y dimension of the grid area
	 * @param mask_xdim the x dimension of the sample map file
	 * @param mask_ydim the y dimension of the sample map file
	 * @param xoffset the x offset of the grid area from the sample map file
	 * @param yoffset the y offset of the grid area from the sample map file
	 * @param inputfile the path to the sample map file
	 */
	void importBooleanMask(unsigned long xdim, unsigned long ydim, unsigned long mask_xdim, unsigned long mask_ydim,
						   unsigned long xoffset, unsigned long yoffset, string inputfile);

	/**
	 * @brief Imports the boolean map object.
	 */
	void doImport();

	/**
	 * @brief Imports the specified file for the sampling percentage within each cell.
	 *
	 * The map should consist of floating points representing the relative sampling rates in each cell. Note that the
	 * actual sampling proportion is equal to the cell value multiplied by global deme sampling proportion.
	 * @param mapvarin the SimParameters object containing the samplemask file location and dimensions
	 */
	void importSampleMask(SimParameters &mapvarin);


	/**
	 * @brief Calculates the matrix value at the provided x, y location.
	 * If everywhere is sampled, simply returns true, as no sample_mask will be stored in memory.
	 * This is to save RAM where possible.
	 *
	 * @param x the x position on the grid
	 * @param y the y position on the grid
	 * @param xval the number of x wraps
	 * @param yval the number of y wraps
	 * @return the sample_mask value at x,y (or true if the file was "null")
	 */
	bool getVal(const long &x, const long &y, const long &xwrap, const long &ywrap);

	/**
	 * @brief Separate return function for always returning 1.0 as density value.
	 * @param x the x position on the grid
	 * @param y the y position on the grid
	 * @param xwrap the number of x wraps around the map
	 * @param ywrap the number of y wraps around the map
	 * @return the sample_mask_exact value at (x, y)
	 */
	double getNullProportion(const long &x, const long &y, const long &xwrap, const long &ywrap);

	/**
	 * @brief Returns the exact value from the spatial sampling map, for calculating the proportion of individuals
	 * to be sampled in each cell.
	 * @note this function assumes that the file is not "null" and the exact sampling mask has been imported. No error
	 * checks on these conditions are performed except in debugging mode.
	 * @param x the x position on the grid
	 * @param y the y position on the grid
	 * @param xwrap the number of x wraps around the map
	 * @param ywrap the number of y wraps around the map
	 * @return the sample_mask_exact value at (x, y)
	 */
	double getBoolProportion(const long &x, const long &y, const long &xwrap, const long &ywrap);

	/**
	 * @brief Returns the exact value from the spatial sampling map, for calculating the proportion of individuals
	 * to be sampled in each cell.
	 * @note this function assumes that the file is not "null" and the exact sampling mask has been imported. No error
	 * checks on these conditions are performed except in debugging mode.
	 * @param x the x position on the grid
	 * @param y the y position on the grid
	 * @param xwrap the number of x wraps around the map
	 * @param ywrap the number of y wraps around the map
	 * @return the sample_mask_exact value at (x, y)
	 */
	double getSampleProportion(const long &x, const long &y, const long &xwrap, const long &ywrap);

	/**
	 * @brief Returns the exact value from the spatial sampling map, as returned by the pointer function.
	 * @param x the x position on the grid
	 * @param y the y position on the grid
	 * @param xwrap the number of x wraps around the map
	 * @param ywrap the number of y wraps around the map
	 * @return the sample_mask_exact value at (x, y)
	 */
	double getExactValue(const long &x, const long &y, const long &xwrap, const long &ywrap);

	/**
	 * @brief Converts the spatial map into the boolean grid required for continued simulation.
	 * This is done so that the faster boolean accesses are possible.
	 * @param map1 the map object to obtain density values from
	 * @param deme_sampling the proportion of individuals to sample
	 * @param generation the generation individuals are added at
	 */
	void convertBoolean(Map &map1, const double &deme_sampling, const double &generation);

	/**
	 * @brief Removes the spatial mask from memory. This should be performed if no more map expansions are required.
	 */
	void clearSpatialMask();

	/**
	 * @brief Converts the coordinates back into the grid format. Changes the values in the provided variables to
	 * be correct.
	 *
	 * @param x the x value to convert
	 * @param y the y value to convert
	 * @param x_wrap the xwrap variable to place the value into
	 * @param y_wrap the ywrap variable to place the value into
	 */
	void recalculate_coordinates(long &x, long &y, long &x_wrap, long &y_wrap);
};


#endif //SPECIATIONCOUNTER_DataMask_H
