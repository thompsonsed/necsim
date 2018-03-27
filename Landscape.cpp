// This file is part of NECSim project which is released under BSD-3 license.
// See file **LICENSE.txt** or visit https://opensource.org/licenses/BSD-3-Clause) for full license details
/**
 * @author Samuel Thompson
 * @file Landscape.cpp
 * @brief Contains the Landscape class implementation for easy referencing of the respective coarse and fine map within the
 * same coordinate system.
 * @copyright <a href="https://opensource.org/licenses/BSD-3-Clause">BSD-3 Licence.</a>
 */

#include "Landscape.h"
#include "Filesystem.h"

uint32_t importToMapAndRound(string map_file, Map<uint32_t> &matrix_in, unsigned long matrix_x,
							 unsigned long matrix_y,
							 unsigned long scalar)
{
#ifndef SIZE_LIMIT
	if(matrix_x > 1000000 || matrix_y > 1000000)
	{
		throw runtime_error("Extremely large map sizes set for " + map_file + ": " + to_string(matrix_x) + ", " +
									to_string(matrix_y) + "\n");
	}
#endif

	Map<float> temp_matrix;
//	temp_matrix.open(map_file);
	temp_matrix.setSize(matrix_y, matrix_x);
#ifdef DEBUG
	writeInfo("Calculating fine map");
#endif
	if(map_file == "null")
	{
		for(unsigned long i = 0; i < matrix_y; i++)
		{
			for(unsigned long j = 0; j < matrix_x; j++)
			{
				temp_matrix[i][j] = 1.0;
			}
		}
	}
	else  // There is a map to read in.
	{
		temp_matrix.import(map_file);
	}
#ifdef DEBUG
	writeInfo("import complete");
#endif
	uint32_t max_value = 0;
	matrix_in.setSize(temp_matrix.getRows(), temp_matrix.getCols());
	for(unsigned long i = 0; i < matrix_y; i++)
	{
		for(unsigned long j = 0; j < matrix_x; j++)
		{
			matrix_in[i][j] = (uint32_t)(max(round(temp_matrix[i][j] * (double)scalar), 0.0));
			if(matrix_in[i][j] > max_value)
			{
				max_value = matrix_in[i][j];
			}
		}
	}
	temp_matrix.close();
	return max_value;
}


void Landscape::setDims(SimParameters mapvarsin)
{
	if(!check_set_dim)  // checks to make sure it hasn't been run already.
	{
		mapvars = mapvarsin;
		deme = mapvarsin.deme;
		x_dim = mapvarsin.grid_x_size;
		y_dim = mapvarsin.grid_y_size;
		scale = mapvarsin.coarse_map_scale;
		nUpdate = 0;
		check_set_dim = true;
		update_time = 0;
		gen_since_pristine = mapvarsin.gen_since_pristine;
		if(gen_since_pristine == 0)
		{
			gen_since_pristine = 0.000000000000000001;
		}
		habitat_change_rate = mapvarsin.habitat_change_rate;
		landscape_type = mapvarsin.landscape_type;
	}
	else
	{
		writeError("ERROR_MAP_001: Dimensions have already been set");
	}
}

bool Landscape::checkMapExists()
{
	for(unsigned int i = 0; i < mapvars.configs.getSectionOptionsSize(); i++)
	{
		string tmppath = mapvars.configs[i].getOption("path");
		if(!doesExistNull(tmppath))
		{
			return false;
		}
	}
	return true;
}

void Landscape::calcFineMap()
{
	string fileinput = mapvars.fine_map_file;
	unsigned long mapxsize = mapvars.fine_map_x_size;
	unsigned long mapysize = mapvars.fine_map_y_size;
	if(!check_set_dim)  // checks that the dimensions have been set.
	{
		throw FatalException("ERROR_MAP_002: dimensions not set.");
	}
	// Note that the default "null" type is to have 100% forest cover in every cell.
	fine_max = importToMapAndRound(fileinput, fine_map, mapxsize, mapysize, deme);
}

void Landscape::calcPristineFineMap()
{
	string file_input = mapvars.pristine_fine_map_file;
	unsigned long map_x_size = mapvars.fine_map_x_size;
	unsigned long map_y_size = mapvars.fine_map_y_size;
	if(!check_set_dim)  // checks that the dimensions have been set.
	{
		throw FatalException("ERROR_MAP_002: dimensions not set.");
	}
	has_pristine = file_input != "none";
	pristine_fine_max = 0;
	if(has_pristine)
	{
		pristine_fine_max = importToMapAndRound(file_input, pristine_fine_map, map_x_size, map_y_size, deme);
	}
}

void Landscape::calcCoarseMap()
{
	string file_input = mapvars.coarse_map_file;
	unsigned long map_x_size = mapvars.coarse_map_x_size;
	unsigned long map_y_size = mapvars.coarse_map_y_size;
	if(!check_set_dim)  // checks that the dimensions have been set.
	{
		throw FatalException("ERROR_MAP_003: dimensions not set.");
	}
	has_coarse = file_input != "none";
	coarse_max = 0;
	if(has_coarse)
	{
		coarse_max = importToMapAndRound(file_input, coarse_map, map_x_size, map_y_size, deme);
	}
}

void Landscape::calcPristineCoarseMap()
{
	string file_input = mapvars.pristine_coarse_map_file;
	unsigned long map_x_size = mapvars.coarse_map_x_size;
	unsigned long map_y_size = mapvars.coarse_map_y_size;
	if(!check_set_dim)  // checks that the dimensions have been set.
	{
		throw FatalException("ERROR_MAP_003: dimensions not set.");
	}
	pristine_coarse_max = 0;
	if(has_coarse)
	{
		has_pristine = file_input != "none";
		if(has_pristine)
		{
			pristine_coarse_max = importToMapAndRound(file_input, pristine_coarse_map, map_x_size, map_y_size, deme);
		}
	}
}

void Landscape::setTimeVars(double gen_since_pristine_in, double habitat_change_rate_in)
{
	update_time = 0;
	gen_since_pristine = gen_since_pristine_in;
	habitat_change_rate = habitat_change_rate_in;
}

void Landscape::calcOffset()
{
	if(mapvars.times_file != "null")
	{
		mapvars.setPristine(0);
	}
	//	os << mapvars.times_file << endl;
	if(fine_map.getCols() == 0 || fine_map.getRows() == 0)
	{
		throw FatalException("ERROR_MAP_004: fine map not set.");
	}
	if(coarse_map.getCols() == 0 || coarse_map.getRows() == 0)
	{
		if(has_coarse)
		{
			coarse_map.setSize(fine_map.getRows(), fine_map.getCols());
		}
		//		throw FatalException("ERROR_MAP_004: coarse map not set.");
	}
	if(checkAllDimensionsZero())
	{
		calculateOffsetsFromMaps();
	}
	else
	{
		calculateOffsetsFromParameters();
	}
	dispersal_relative_cost = mapvars.dispersal_relative_cost;
#ifdef DEBUG
	stringstream os;
	os << "\nfinex: " << fine_x_min << "," << fine_x_max << endl;
	os << "finey: " << fine_y_min << "," << fine_y_max << endl;
	os << "coarsex: " << coarse_x_min << "," << coarse_x_max << endl;
	os << "coarsey: " << coarse_y_min << "," << coarse_y_max << endl;
	os << "offsets: "
		 << "(" << fine_x_offset << "," << fine_y_offset << ")(" << coarse_x_offset << "," << coarse_y_offset << ")" << endl;
	os << "pristine fine file: " << pristine_fine_map << endl;
	os << "pristine coarse file: " << pristine_coarse_map << endl;
	writeInfo(os.str());
#endif
	//		os << "fine variables: " << finexmin << "," << fine_x_max << endl;
	//		os << "coarse variabes: " << coarse_x_min << "," << coarse_x_max << endl;
	if(fine_x_min < coarse_x_min || fine_x_max > coarse_x_max || (fine_x_max - fine_x_min) < x_dim || (fine_y_max - fine_y_min) < y_dim)
	{
		throw FatalException(
			"ERROR_MAP_006: FATAL - fine map extremes outside coarse map or sample grid larger than fine map");
	}
}

bool Landscape::checkAllDimensionsZero()
{
	return mapvars.fine_map_x_offset == 0 && mapvars.fine_map_y_offset == 0 && mapvars.coarse_map_x_offset == 0 &&
			mapvars.coarse_map_y_offset == 0 && mapvars.sample_x_offset == 0 && mapvars.sample_y_offset == 0 &&
			mapvars.fine_map_x_size == 0 && mapvars.fine_map_y_size == 0 && mapvars.coarse_map_x_size == 0 &&
			mapvars.coarse_map_y_size == 0;
}

void Landscape::calculateOffsetsFromMaps()
{
	// TODO complete this functionality - required for rcoalescence

}

void Landscape::calculateOffsetsFromParameters()
{
	fine_x_offset = mapvars.fine_map_x_offset + mapvars.sample_x_offset;
	fine_y_offset = mapvars.fine_map_y_offset + mapvars.sample_y_offset;
	coarse_x_offset = mapvars.coarse_map_x_offset;
	coarse_y_offset = mapvars.coarse_map_y_offset;
	scale = mapvars.coarse_map_scale;
	// this is the location of the top left (or north west) corner of the respective map
	// and the x and y distance from the top left of the grid object that contains the initial lineages.
	fine_x_min = -fine_x_offset;
	fine_y_min = -fine_y_offset;
	fine_x_max = fine_x_min + (fine_map.getCols());
	fine_y_max = fine_y_min + (fine_map.getRows());
	if(has_coarse) // Check if there is a coarse map
	{
		coarse_x_min = -coarse_x_offset - fine_x_offset;
		coarse_y_min = -coarse_y_offset - fine_y_offset;
		coarse_x_max = coarse_x_min + scale * (coarse_map.getCols());
		coarse_y_max = coarse_y_min + scale * (coarse_map.getRows());
	}
	else // Just set the offsets to the same as the fine map
	{
		coarse_x_min = fine_x_min;
		coarse_y_min = fine_y_min;
		coarse_x_max = fine_x_max;
		coarse_y_max = fine_y_max;
		scale = 1;
	}
}



void Landscape::validateMaps()
{
	stringstream os;
	os << "\rValidating maps..." << flush;
	double dTotal = fine_map.getCols() + coarse_map.getCols();
	unsigned long iCounter = 0;
	if(has_pristine)
	{
		if(fine_map.getCols() == pristine_fine_map.getCols() && fine_map.getRows() == pristine_fine_map.getRows() &&
				coarse_map.getCols() == pristine_coarse_map.getCols() && coarse_map.getRows() ==
																			pristine_coarse_map.getRows())
		{
			os << "\rValidating maps...map sizes okay" << flush;
			writeInfo(os.str());
		}
		else
		{
			throw FatalException(
					"ERROR_MAP_009: Landscape validation failed - modern and pristine maps are not the same dimensions.");
		}
		for(unsigned long i = 0; i < fine_map.getCols(); i++)
		{
			for(unsigned long j = 0; j < fine_map.getRows(); j++)
			{
				if(fine_map[j][i] > pristine_fine_map[j][i])
				{
#ifdef DEBUG
					stringstream ss;
					ss << "fine map: " << fine_map[j][i] << " pristine map: " << pristine_fine_map[j][i];
					ss << " x,y: " << i << "," << j << endl;
					writeLog(50, ss);
#endif //DEBUG
					throw FatalException("ERROR_MAP_007: Landscape validation failed - fine map value larger "
											  "than pristine fine map value.");
				}
			}
			double dPercentComplete = 100 * ((double)(i + iCounter) / dTotal);
			if(i % 1000 == 0)
			{
				os.str("");
				os << "\rValidating maps..." << dPercentComplete << "%                " << flush;
				writeInfo(os.str());
			}
		}
	}
	iCounter = fine_map.getCols();
	if(has_pristine)
	{
		for(unsigned long i = 0; i < coarse_map.getCols(); i++)
		{
			for(unsigned long j = 0; j < coarse_map.getRows(); j++)
			{
				if(coarse_map[j][i] > pristine_coarse_map[j][i])
				{
#ifdef DEBUG
					stringstream ss;
					ss << "coarse map: " << coarse_map[j][i] << " pristine map: " << pristine_coarse_map[j][i];
					ss << " coarse map x+1: " << coarse_map[j][i + 1]
						 << " pristine map: " << pristine_coarse_map[j][i + 1];
					ss << " x,y: " << i << "," << j;
					writeLog(50, ss);
#endif // DEBUG
					throw FatalException("ERROR_MAP_008: Landscape validation failed - coarse map value larger "
											  "than pristine coarse map value.");
				}
			}
			double dPercentComplete = 100 * ((double)(i + iCounter) / dTotal);
			if(i % 1000 == 0)
			{
				os.str("");
				os << "\rValidating maps..." << dPercentComplete << "%                " << flush;
				writeInfo(os.str());
			}
		}
		
	}
	os.str("");
	os << "\rValidating maps complete                                       " << endl;
	writeInfo(os.str());
}

void Landscape::updateMap(double generation)
{
	// only update the map if the pristine state has not been reached.
	if(!mapvars.is_pristine && has_pristine)
	{
		if(mapvars.gen_since_pristine < generation)
		{
			// Only update the map if the maps have actually changed
			if(mapvars.setPristine(nUpdate+1))
			{
				nUpdate++;
				// pristine_fine_map = mapvars.pristine_fine_map_file;
				// pristine_coarse_map = mapvars.pristine_coarse_map_file;
				current_map_time = gen_since_pristine;
				gen_since_pristine = mapvars.gen_since_pristine;
				if(gen_since_pristine == 0)
				{
					gen_since_pristine = 0.000000000000000001;
				}
				habitat_change_rate = mapvars.habitat_change_rate;
				fine_max = pristine_fine_max;
				fine_map = pristine_fine_map;
				coarse_max = pristine_coarse_max;
				coarse_map = pristine_coarse_map;
				calcPristineCoarseMap();
				calcPristineFineMap();
				if(has_pristine)
				{
					is_pristine = mapvars.is_pristine;
				}
				recalculateHabitatMax();
				
				
			}
		}
	}
}

void Landscape::setLandscape(string landscape_type)
{
	if(landscape_type == "infinite")
	{
		writeInfo("Setting infinite landscape.\n");
		getValFunc = &Landscape::getValInfinite;
	}
	else if(landscape_type == "tiled_coarse")
	{
		writeInfo("Setting tiled coarse infinite landscape.\n");
		getValFunc = &Landscape::getValCoarseTiled;
	}
	else if(landscape_type == "tiled_fine")
	{
		writeInfo("Setting tiled fine infinite landscape.\n");
		getValFunc = &Landscape::getValFineTiled;
	}
	else if(landscape_type == "closed")
	{
		getValFunc = &Landscape::getValFinite;
	}
	else
	{
		throw FatalException("Provided landscape type is not a valid option: " + landscape_type);
	}
}

unsigned long Landscape::getVal(const double& x, const double& y,
						  const long& xwrap, const long& ywrap, const double& current_generation)
{
	return (this->*getValFunc)(x, y, xwrap, ywrap, current_generation);
}

unsigned long Landscape::getValInfinite(
	const double& x, const double& y, const long& xwrap, const long& ywrap, const double& current_generation)
{
	double xval, yval;
	xval = x + (x_dim * xwrap);  //
	yval = y + (y_dim * ywrap);
	//		// return 0 if the requested coordinate is completely outside the map
	if(xval < coarse_x_min || xval >= coarse_x_max || yval < coarse_y_min || yval >= coarse_y_max)
	{
		return deme;
	}
	return getValFinite(x, y, xwrap, ywrap, current_generation);
}

unsigned long Landscape::getValCoarseTiled(
	const double& x, const double& y, const long& xwrap, const long& ywrap, const double& current_generation)
{
	double newx = fmod(x + (xwrap * x_dim) + fine_x_offset + coarse_x_offset, coarse_map.getCols());
	double newy = fmod(y + (ywrap * y_dim) + fine_x_offset + coarse_x_offset, coarse_map.getRows());
	if(newx < 0)
	{
		newx += coarse_map.getCols();
	}
	if(newy < 0)
	{
		newy += coarse_map.getRows();
	}
	return getValCoarse(newx, newy, current_generation);
}

unsigned long Landscape::getValFineTiled(
	const double& x, const double& y, const long& xwrap, const long& ywrap, const double& current_generation)
{

	double newx = fmod(x + (xwrap * x_dim) + fine_x_offset, fine_map.getCols());
	double newy = fmod(y + (ywrap * y_dim) + fine_y_offset, fine_map.getRows());
	// Now adjust for incorrect wrapping behaviour of fmod
	if(newx < 0)
	{
		newx += fine_map.getCols();
	}
	if(newy < 0)
	{
		newy += fine_map.getRows();
	}
#ifdef DEBUG
	if(newx >= fine_map.getCols() || newx < 0 || newy >= fine_map.getRows() || newy < 0)
	{
		stringstream ss;
		ss << "Fine map indexing out of range of fine map." << endl;
		ss << "x, y: " << newx << ", " << newy << endl;
		ss << "cols, rows: " << fine_map.getCols() << ", " << fine_map.getRows() << endl;
		throw out_of_range(ss.str());
	}
#endif
	return getValFine(newx, newy, current_generation);
}

unsigned long Landscape::getValCoarse(const double &xval, const double &yval, const double &current_generation)
{
	unsigned long retval = 0;
	if(has_pristine)
	{
		if(is_pristine || pristine_coarse_map[yval][xval] == coarse_map[yval][xval])
		{
			return pristine_coarse_map[yval][xval];
		}
		else
		{
			double currentTime = current_generation - current_map_time;
			retval = (unsigned long)floor(coarse_map[yval][xval] +
										   (habitat_change_rate *
											((pristine_coarse_map[yval][xval] - coarse_map[yval][xval]) /
													(gen_since_pristine-current_map_time)) * currentTime));
		}
	}
	else
	{
		return coarse_map[yval][xval];
	}
#ifdef pristine_mode
	if(retval > pristine_coarse_map[yval][xval])
		{
			string ec =
				"Returned value greater than pristine value. Check file input. (or disable this error before "
				"compilation.\n";
			ec += "pristine value: " + to_string((long long)pristine_coarse_map[yval][xval]) +
				  " returned value: " + to_string((long long)retval);
			throw FatalException(ec);
		}
// Note that debug mode will throw an exception if the returned value is less than the pristine state

#endif
	return retval;
}

unsigned long Landscape::getValFine(const double&xval, const double &yval, const double& current_generation)
{
	unsigned long retval = 0;
	if(has_pristine)
	{
		if(is_pristine || pristine_fine_map[yval][xval] == fine_map[yval][xval])
		{
			retval = pristine_fine_map[yval][xval];
		}
		else
		{
			double currentTime = current_generation - current_map_time;
			retval = (unsigned long)floor(fine_map[yval][xval] +
										   (habitat_change_rate * ((pristine_fine_map[yval][xval] - fine_map[yval][xval]) /
												   (gen_since_pristine-current_map_time)) * currentTime));
		}
	}
	else
	{
		return fine_map[yval][xval];
	}
// os <<fine_map[yval][xval] << "-"<< retval << endl;
// Note that debug mode will throw an exception if the returned value is less than the pristine state
#ifdef pristine_mode
	if(has_pristine)
	{
		if(retval > pristine_fine_map[yval][xval])
		{
			throw FatalException("Returned value greater than pristine value. Check file input. (or disable this "
									  "error before compilation.");
		}
	}
#endif
	return retval;
}

unsigned long Landscape::getValFinite(
	const double& x, const double& y, const long& xwrap, const long& ywrap, const double& current_generation)
{

	double xval, yval;
	xval = x + (x_dim * xwrap);  //
	yval = y + (y_dim * ywrap);
	//		// return 0 if the requested coordinate is completely outside the map
	if(xval < coarse_x_min || xval >= coarse_x_max || yval < coarse_y_min || yval >= coarse_y_max)
	{
		return 0;
	}
	if((xval < fine_x_min || xval >= fine_x_max || yval < fine_y_min ||
	   yval >= fine_y_max) && has_coarse)  // check if the coordinate comes from the coarse resolution map.
	{
		// take in to account the fine map offsetting
		xval += fine_x_offset;
		yval += fine_y_offset;
		// take in to account the coarse map offsetting and the increased scale of the larger map.
		xval = floor((xval + coarse_x_offset) / scale);
		yval = floor((yval + coarse_y_offset) / scale);
		return getValCoarse(xval, yval, current_generation);
	}
	// take in to account the fine map offsetting
	// this is done twice to avoid having all the comparisons involve additions.
	xval += fine_x_offset;
	yval += fine_y_offset;
	return getValFine(xval, yval, current_generation);

}

unsigned long Landscape::convertSampleXToFineX(const unsigned long &x, const long &xwrap)
{
	return x + fine_x_offset + (xwrap * x_dim);
}

unsigned long Landscape::convertSampleYToFineY(const unsigned long &y, const long &ywrap)
{
	return y + fine_y_offset + (ywrap * y_dim);
}

void Landscape::convertFineToSample(long & x, long & xwrap, long &y, long &ywrap)
{
	auto tmpx = double(x);
	auto tmpy = double(y);
	convertCoordinates(tmpx, tmpy, xwrap, ywrap);
	x = static_cast<long>(floor(tmpx));
	y = static_cast<long>(floor(tmpy));
}


unsigned long Landscape::getInitialCount(double dSample, DataMask& samplemask)
{
	unsigned long toret;
	toret = 0;
	long x, y;
	long xwrap, ywrap;
	unsigned long max_x, max_y;
	if(samplemask.getDefault())
	{
		max_x = fine_map.getCols();
		max_y = fine_map.getRows();
	}
	else
	{
		max_x = samplemask.sample_mask.getCols();
		max_y = samplemask.sample_mask.getRows();
	}
	for(unsigned long i = 0; i < max_x; i++)
	{
		for(unsigned long j = 0; j < max_y; j++)
		{
			x = i;
			y = j;
			xwrap = 0;
			ywrap = 0;
			samplemask.recalculate_coordinates(x, y, xwrap, ywrap);
			toret += (unsigned long) (max(floor(dSample * (getVal(x, y, xwrap, ywrap, 0)) *
												samplemask.getExactValue(x, y, xwrap, ywrap)), 0.0));
		}
	}
	return toret;
}

SimParameters Landscape::getSimParameters()
{
	return mapvars;
}

bool Landscape::checkMap(const double& x, const double& y, const long& xwrap, const long& ywrap, const double generation)
{
	return getVal(x, y, xwrap, ywrap, generation) != 0;
}

bool Landscape::checkFine(const double& x, const double& y, const long& xwrap, const long& ywrap)
{
	double tmpx, tmpy;
	tmpx = x + xwrap * x_dim;
	tmpy = y + ywrap * y_dim;
	return !(tmpx < fine_x_min || tmpx >= fine_x_max || tmpy < fine_y_min || tmpy >= fine_y_max);
}

void Landscape::convertCoordinates(double& x, double& y, long& xwrap, long& ywrap)
{
	xwrap += floor(x / x_dim);
	ywrap += floor(y / y_dim);
	x = x - xwrap * x_dim;
	y = y - ywrap * y_dim;
}

unsigned long Landscape::runDispersal(const double& dist,
					   const double& angle,
					   long& startx,
					   long& starty,
					   long& startxwrap,
					   long& startywrap,
					   bool& disp_comp,
					   const double& generation)
{
// Checks that the start point is not out of matrix - this might have to be disabled to ensure that when updating the
// map, it doesn't cause problems.
#ifdef pristine_mode
	if(!checkMap(startx, starty, startxwrap, startywrap, generation))
	{
		disp_comp = true;
		return;
	}
#endif

	// Different calculations for each quadrant to ensure that the dispersal reads the probabilities correctly.
	double newx, newy;
	newx = startx + (x_dim * startxwrap) + 0.5;
	newy = starty + (y_dim * startywrap) + 0.5;
	if(dispersal_relative_cost ==1)
	{
		// then nothing complicated is required and we can jump straight to the final point.
		newx += dist * cos(angle);
		newy += dist * sin(angle);
	}
	else  // we need to see which deforested patches we pass over
	{
		long boost;
		boost = 1;
		double cur_dist, tot_dist, l;
		cur_dist = 0;
		tot_dist = 0;
		// Four different calculations for the different quadrants.
		if(angle > 7 * M_PI_4 || angle <= M_PI_4)
		{
			// Continue while the dist travelled is less than the dist energy
			while(cur_dist < dist)
			{
				// Check if the starting position of the loop is in the fine map or not.
				if(checkFine(newx, newy, 0, 0))
				{
					// Keep the standard movement rate
					boost = 1;
				}
				else
				{
					// Accellerate the travel speed if the point is outside the fine grid.
					// Note this means that lineages travelling from outside the fine grid to within the
					// fine grid may
					// see 1 grid's worth of approximation, rather than exact values.
					// This is an acceptable approximation!
					boost = deme;
				}

				// Add the value to the new x and y values.
				newx = newx + boost;
				newy = newy + boost * tan(angle);
				// Check if the new point is within forest.
				if(checkMap(newx, newy, 0, 0, generation))
				{
					l = 1;
				}
				else
				{
					l = dispersal_relative_cost;
				}
				// Move forward different dists based on the difficulty of moving through forest.
				cur_dist = cur_dist + l * boost * (1 / cos(angle));
				tot_dist = tot_dist + boost * (1 / cos(angle));
			}
		}
		else if(angle > 3 * M_PI_4 && angle <= 5 * M_PI_4)
		{
			while(cur_dist < dist)
			{
				if(checkFine(newx, newy, 0, 0))
				{
					boost = 1;
				}
				else
				{
					boost = deme;
				}
				// Add the change to the new x and y values.
				newx = newx - boost;
				newy = newy + boost * tan(M_PI - angle);
				if(checkMap(newx, newy, 0, 0, generation))
				{
					l = 1;
				}
				else
				{
					l = dispersal_relative_cost;
				}
				cur_dist = cur_dist + boost * l * (1 / cos(M_PI - angle));
				tot_dist = tot_dist + boost * (1 / cos(M_PI - angle));
			}
		}
		else if(angle > M_PI_4 && angle <= 3 * M_PI_4)
		{
			while(cur_dist < dist)
			{
				if(checkFine(newx, newy, 0, 0))
				{
					boost = 1;
				}
				else
				{
					boost = deme;
				}
				// Add the change to the new x and y values.
				newx = newx + boost * tan(angle - M_PI_2);
				newy = newy + boost;
				if(checkMap(newx, newy, 0, 0, generation))
				{
					l = 1;
				}
				else
				{
					l = dispersal_relative_cost;
				}
				cur_dist = cur_dist + l * boost / cos(angle - M_PI_2);
				tot_dist = tot_dist + boost / cos(angle - M_PI_2);
			}
		}
		else if(angle > 5 * M_PI_4 && angle <= 7 * M_PI_4)
		{
			//				os << "...ang4..." <<  flush;
			while(cur_dist < dist)
			{
				if(checkFine(newx, newy, 0, 0))
				{
					boost = 1;
				}
				else
				{
					boost = deme;
				}
				newx = newx + boost * tan(3 * M_PI_2 - angle);
				newy = newy - boost;
				if(checkMap(newx, newy, 0, 0, generation))
				{
					l = 1;
				}
				else
				{
					l = dispersal_relative_cost;
				}
				cur_dist = cur_dist + l * boost / cos(3 * M_PI_2 - angle);
				tot_dist = tot_dist + boost / cos(3 * M_PI_2 - angle);
			}
		}
		// Move the point back to get the exact placement
		if(checkMap(newx, newy, 0, 0, generation))
		{
			tot_dist = tot_dist - min(cur_dist - dist, (double(boost) - 0.001));
		}
		else
		{
			disp_comp = true;
		}
		newx = startx + 0.5 + tot_dist * cos(angle);
		newy = starty + 0.5 + tot_dist * sin(angle);
	}
	unsigned long ret = getVal(newx, newy, 0, 0, generation);
	if(ret >0)
	{
		long newxwrap, newywrap;
		newxwrap = 0;
		newywrap = 0;
		convertCoordinates(newx, newy, newxwrap, newywrap);
#ifdef DEBUG
		if(!checkMap(newx, newy, newxwrap, newywrap, generation))
		{
			throw FatalException(string(
				"ERROR_MOVE_007: Dispersal attempted to non-forest. Check dispersal function. Forest cover: " +
				to_string((long long)getVal(newx, newy, newxwrap, newywrap, generation))));
		}
#endif
		startx = newx;
		starty = newy;
		startxwrap = newxwrap;
		startywrap = newywrap;
		disp_comp = false;
	}
	return ret;
};

void Landscape::clearMap()
{
	current_map_time = 0;
	check_set_dim = false;
	is_pristine = false;
}

string Landscape::printVars()
{
	stringstream os;
	os << "fine x limits: " << fine_x_min << " , " << fine_x_max << endl;
	os << "fine y limits: " << fine_y_min << " , " << fine_y_max << endl;
	os << "fine map offset: " << fine_x_offset << " , " << fine_y_offset << endl;
	os << "coarse x limits: " << coarse_x_min << " , " << coarse_x_max << endl;
	os << "coarse y limits: " << coarse_y_min << " , " << coarse_y_max << endl;
	os << "x,y dims: " << x_dim << " , " << y_dim << endl;
	return os.str();
}

unsigned long Landscape::getHabitatMax()
{
	return habitat_max;
}

void Landscape::recalculateHabitatMax()
{
	habitat_max = 0;
	if(is_pristine && has_pristine)
	{
		if(habitat_max < pristine_fine_max)
		{
			habitat_max = pristine_fine_max;
		}
		if(habitat_max < pristine_coarse_max)
		{
			habitat_max = pristine_coarse_max;
		}
	}
	else
	{
		if(habitat_max < fine_max)
		{
			habitat_max = fine_max;
		}
		if(habitat_max < coarse_max)
		{
			habitat_max = coarse_max;
		}
		if(habitat_max < pristine_fine_max)
		{
			habitat_max = pristine_fine_max;
		}
		if(habitat_max < pristine_coarse_max)
		{
			habitat_max = pristine_coarse_max;
		}
	}
#ifdef DEBUG
	if(habitat_max > 10000)
	{
		stringstream ss;
		writeLog(10, "habitat_max may be unreasonably large: " + to_string(habitat_max));
		ss << "fine, coarse, pfine, pcoarse: " << fine_max << ", " << coarse_max;
		ss << ", " << pristine_fine_max << ", " << pristine_coarse_max << endl;
	}
#endif
}

