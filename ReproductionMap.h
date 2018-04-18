//This file is part of NECSim project which is released under MIT license.
//See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.

/**
 * @author Samuel Thompson
 * @date 16/08/2017
 * @file ReproductionMap.h
 * @brief Contains the ReproductionMap, which inherits from Matrix and adds a few extra parameters.
 *
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 */

#ifndef REPRODUCTIONMAP_H
#define REPRODUCTIONMAP_H

#include <cstring>
#include <string>
#include <cstdio>
#include <iostream>
#include <fstream>

#include "Map.h"
#include "NRrand.h"


/**
 * @brief Contains the routines for importing the reproduction map and getting a cell value from the map.
 */
class ReproductionMap
{
protected:
	// Matrix containing the relative reproduction probabilities
	Map<double> reproduction_map;
	// Path to the map file
	string map_file;
	// Maximum value across the map
	double max_val;
	// If true,
	bool null_map;
	// The fine map offsets and the sample map dimensions
	unsigned long offset_x, offset_y, x_dim, y_dim;
	// Function pointer for our reproduction map checker
	typedef bool (ReproductionMap::*rep_ptr)(NRrand &random_no,
											 const unsigned long &x, const unsigned long &y,
											 const long &xwrap, const long &ywrap);
	// once setup will contain the end check function to use for this simulation.
	rep_ptr reproductionMapChecker_fptr;
public:
	ReproductionMap()
	{
		map_file = "none";
		max_val = 0;
		null_map = true;
	}

	/**
	 * @brief Imports the map file from the given path
	 * Requires the dimensions to be identical as the fine map file dimensions
	 * @param file_name the path to the reproduction map to import
	 * @param size_x the x dimensions of the map file
	 * @param size_y the y dimensions of the map file
	 */
	void import(string file_name, unsigned long size_x, unsigned long size_y);

	/**
	 * @brief Correctly sets the reproduction function to either rejectionSampleNull or rejectionSample
	 * depending on if a reproduction map is used or not.
	 */
	void setReproductionFunction();

	/**
	 * @brief Sets the offsets for the reproduction map from the sample grid.
	 * @param x_offset the x offset from the sample grid
	 * @param y_offset the y offset from the sample grid
	 * @param xdim the x dimension of the sample grid
	 * @param ydim the y dimension of the sample grid
	 */
	void setOffsets(const unsigned long &x_offset, const unsigned long &y_offset,
					const unsigned long &xdim, const unsigned long &ydim);

	/**
	 * @brief Returns true for all cell values
	 * Function to be pointed to in cases where there is no reproduction map
	 * @param random_number random number object to pass forward
	 * @param x x coordinate of the lineage on the sample grid
	 * @param y y coordinate of the lineage on the sample grid
	 * @param xwrap x wrapping of the lineage
	 * @param ywrap y wrapping of the lineage
	 * @return true always
	 */
	bool rejectionSampleNull(NRrand &random_number, const unsigned long &x, const unsigned long &y, const long &xwrap,
							 const long &ywrap);

	/**
	 * @brief Returns true for all cell values
	 * Function to be pointed to in cases where there is no reproduction map
	 * @param random_number random number object to pass forward
	 * @param x x coordinate of the lineage on the sample grid
	 * @param y y coordinate of the lineage on the sample grid
	 * @param xwrap x wrapping of the lineage
	 * @param ywrap y wrapping of the lineage
	 * @return true always
	 */
	bool rejectionSample(NRrand &random_number, const unsigned long &x, const unsigned long &y,
					 const long &xwrap, const long &ywrap);

	/**
	 * @brief Gets the value of the reproduction map at that location
	 * @param x x coordinate of the lineage on the sample grid
	 * @param y y coordinate of the lineage on the sample grid
	 * @param xwrap x wrapping of the lineage
	 * @param ywrap y wrapping of the lineage
	 * @return value of the reproduction map at the required location
	 */
	double getVal(const unsigned long &x, const unsigned long &y, const long &xwrap, const long &ywrap);

	/**
	 * @brief
	 * @param random_number random number object to draw from
	 * @param x x coordinate of the lineage on the sample grid
	 * @param y y coordinate of the lineage on the sample grid
	 * @param xwrap x wrapping of the lineage
	 * @param ywrap y wrapping of the lineage
	 * @return value of the reproduction map at the required location
	 */
	bool hasReproduced(NRrand &random_number, const unsigned long &x, const unsigned long &y,
					   const long &xwrap, const long &ywrap);

	/**
	 * @brief Operator [] for getting values directly from the reproduction map.
	 * @param index the index to get the row of
	 * @return the row present at that index
	 */
	Row<double> operator[](long index)
	{
		return reproduction_map[index];
	}

	/**
	 * @brief Operator for outputting to an ostream.
	 * @param os the ostream to output to
	 * @param r the ReproductionMap to read from
	 * @return the os object
	 */
	friend ostream& operator<<(ostream& os, ReproductionMap&r)
	{
		os << r.map_file << "\n";
		os << r.reproduction_map.getCols() << "\n";
		os << r.reproduction_map.getRows() << "\n";
		os << r.offset_x << "\n";
		os << r.offset_y << "\n";
		os << r.x_dim << "\n";
		os << r.y_dim << "\n";
		return os;
	}

	/**
	 * @brief Operator for inputting from an istream.
	 * @param is the istream to input from
	 * @param r the ReproductionMap to input to
	 * @return the is object
	 */
	friend istream& operator>>(istream &is, ReproductionMap &r)
	{
		is.ignore();
		getline(is, r.map_file);
		unsigned long col, row;
		is >> col >> row;
		is >> r.offset_x >> r.offset_y >> r.x_dim >> r.y_dim;
		r.import(r.map_file, col, row);
		return is;
	}


};


#endif //REPRODUCTIONMAP_H
