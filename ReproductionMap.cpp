//This file is part of NECSim project which is released under MIT license.
//See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.

/**
 * @author Samuel Thompson
 * @date 16/08/2017
 * @file ReproductionMap.cpp
 * @brief Contains the ReproductionMap, which inherits from Matrix and adds a few extra parameters.
 *
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 */

#include "ReproductionMap.h"

void ReproductionMap::import(string file_name, unsigned long size_x, unsigned long size_y)
{
	map_file = file_name;
	if(file_name == "null" || file_name == "none")
	{
		null_map = true;
	}
	else
	{
		null_map = false;
		reproduction_map.setSize(size_y, size_x);
		reproduction_map.import(file_name);
		for(unsigned long y = 0; y < reproduction_map.getRows(); y++)
		{
			for(unsigned long x = 0; x < reproduction_map.getCols(); x++)
			{
				if(reproduction_map[y][x] > max_val)
				{
					max_val = reproduction_map[y][x];
				}
			}
		}
		reproduction_map.close();
	}
	setReproductionFunction();
}

void ReproductionMap::setReproductionFunction()
{
	if(null_map)
	{
		reproductionMapChecker_fptr = &ReproductionMap::rejectionSampleNull;
	}
	else
	{
		reproductionMapChecker_fptr = &ReproductionMap::rejectionSample;
	}
}

void ReproductionMap::setOffsets(const unsigned long &x_offset, const unsigned long &y_offset,
								 const unsigned long &xdim, const unsigned long &ydim)
{
	offset_x = x_offset;
	offset_y = y_offset;
	x_dim = xdim;
	y_dim = ydim;
}


bool ReproductionMap::rejectionSampleNull(NRrand &random_number, const unsigned long &x, const unsigned long &y,
									  const long &xwrap, const long &ywrap)
{
	return true;
}

bool ReproductionMap::rejectionSample(NRrand &random_number, const unsigned long &x, const unsigned long &y,
									  const long &xwrap, const long &ywrap)
{
	return random_number.d01() <= getVal(x, y, xwrap, ywrap);
}

double ReproductionMap::getVal(const unsigned long &x, const unsigned long &y, const long &xwrap, const long &ywrap)
{
	unsigned long x_ref = x + (xwrap * x_dim) + offset_x;
	unsigned long y_ref = y + (ywrap * y_dim) + offset_y;
	return reproduction_map[y_ref][x_ref];
}

bool ReproductionMap::hasReproduced(NRrand &random_number, const unsigned long &x, const unsigned long &y, const long &xwrap,
							   const long &ywrap)
{
	return (this->*reproductionMapChecker_fptr)(random_number, x, y, xwrap, ywrap);
}






