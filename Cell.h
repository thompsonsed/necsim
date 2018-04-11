// This file is part of NECSim project which is released under BSD-3 license.
// See file **LICENSE.txt** or visit https://opensource.org/licenses/BSD-3-Clause) for full license details.

/**
 * @author Samuel Thompson
 * @file Cell.h
 * @brief Basic container for location data.
 * @copyright <a href="https://opensource.org/licenses/BSD-3-Clause">BSD-3 Licence.</a>
 */

#ifndef CELL_H
#define CELL_H

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

#endif // CELL_H