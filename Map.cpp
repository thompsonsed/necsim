// This file is part of NECSim project which is released under BSD-3 license.
// See file **LICENSE.txt** or visit https://opensource.org/licenses/BSD-3-Clause) for full license details
/**
 * @author Samuel Thompson
 * @file Map.cpp
 * @brief Contains the Map class implementation for easy referencing of the respective coarse and fine map within the
 * same coordinate system.
 * @copyright <a href="https://opensource.org/licenses/BSD-3-Clause">BSD-3 Licence.</a>
 */

#include "Map.h"
#include "Filesystem.h"
#include "CustomExceptions.h"

void Map::setDims(SimParameters mapvarsin)
{
	if(!check_set_dim)  // checks to make sure it hasn't been run already.
	{
		mapvars = mapvarsin;
		deme = mapvarsin.deme;
		x_dim = mapvarsin.vargridxsize;
		y_dim = mapvarsin.vargridysize;
		scale = mapvarsin.varcoarsemapscale;
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
		cerr << "ERROR_MAP_001: Dimensions have already been set" << endl;
		return;
	}
	return;
}

bool Map::checkMapExists()
{
	for(unsigned int i = 0; i < mapvars.configs.getSectionOptionsSize(); i++)
	{
		string tmppath = mapvars.configs[i].getOption("path");
		if(!doesExistNull(tmppath))
		{
			return (false);
		}
	}
	return (true);
}

void Map::calcFineMap()
{
	string fileinput = mapvars.finemapfile;
	unsigned long mapxsize = mapvars.varfinemapxsize;
	unsigned long mapysize = mapvars.varfinemapysize;
	if(!check_set_dim)  // checks that the dimensions have been set.
	{
		throw FatalException("ERROR_MAP_002: dimensions not set.");
	}
	// Note that the default "null" type is to have 100% forest cover in every cell.
	#ifndef SIZE_LIMIT
	if(mapxsize > 1000000 || mapysize > 1000000)
	{
		throw runtime_error("Extremely large map sizes set for " + fileinput + ": " + to_string(mapxsize) + ", " + to_string(mapysize) + "\n");
	}
	#endif
	Matrix<float> toret;
	toret.SetSize(mapysize, mapxsize);
	fine_map.SetSize(mapysize, mapxsize);
#ifdef DEBUG
	writeInfo("Calculating fine map");
#endif
	if(fileinput == "null")
	{
		for(unsigned long i = 0; i < mapysize; i++)
		{
			for(unsigned long j = 0; j < mapxsize; j++)
			{
				toret[i][j] = 1.0;
			}
		}
	}
	else  // There is a map to read in.
	{
		toret.import(fileinput);
	}
#ifdef DEBUG
	writeInfo("import complete");
#endif
	fine_max = 0;

	for(unsigned long i = 0; i < mapysize; i++)
	{
		for(unsigned long j = 0; j < mapxsize; j++)
		{
			fine_map[i][j] = (unsigned long)(max((double)round(toret[i][j] * deme), 0.0));
			if(fine_map[i][j] > fine_max)
			{
				fine_max = fine_map[i][j];
			}
		}
	}
}

void Map::calcPristineFineMap()
{
	string fileinput = mapvars.pristinefinemapfile;
	unsigned long mapxsize = mapvars.varfinemapxsize;
	unsigned long mapysize = mapvars.varfinemapysize;
	if(!check_set_dim)  // checks that the dimensions have been set.
	{
		throw FatalException("ERROR_MAP_002: dimensions not set.");
	}
	#ifndef SIZE_LIMIT
	if(mapxsize > 1000000 || mapysize > 1000000)
	{
		throw runtime_error("Extremely large map sizes set for " + fileinput + ": " + to_string(mapxsize) + ", " + to_string(mapysize) + "\n");
	}
	#endif
	// Note that the default "null" type is to have 100% forest cover in every cell.
	Matrix<float> toret;
	has_pristine = true;
	pristine_fine_max = 0;
	if(fileinput == "null")
	{
		toret.SetSize(mapysize, mapxsize);
		pristine_fine_map.SetSize(mapysize, mapxsize);
		for(unsigned long i = 0; i < mapxsize; i++)
		{
			for(unsigned long j = 0; j < mapysize; j++)
			{
				toret[j][i] = 1.0;
			}
		}
	}
	else if(fileinput == "none")
	{
		has_pristine = false;
	}
	else  // There is a map to read in.
	{
		toret.SetSize(mapysize, mapxsize);
		pristine_fine_map.SetSize(mapysize, mapxsize);
		toret.import(fileinput);
	}
	// os << toret << endl;
	if(has_pristine)
	{
		for(unsigned long i = 0; i < mapysize; i++)
		{
			for(unsigned long j = 0; j < mapxsize; j++)
			{
				pristine_fine_map[i][j] = (unsigned long)(max((double)round(toret[i][j] * deme), 0.0));
				if(pristine_fine_map[i][j] > pristine_fine_max)
				{
					pristine_fine_max = pristine_fine_map[i][j];
				}
			}
		}
	}
	return;
}

void Map::calcCoarseMap()
{
	string fileinput = mapvars.coarsemapfile;
	unsigned long mapxsize = mapvars.varcoarsemapxsize;
	unsigned long mapysize = mapvars.varcoarsemapysize;
	if(!check_set_dim)  // checks that the dimensions have been set.
	{
		throw FatalException("ERROR_MAP_003: dimensions not set.");
	}
	#ifndef SIZE_LIMIT
	if(mapxsize > 1000000 || mapysize > 1000000)
	{
		throw runtime_error("Extremely large map sizes set for " + fileinput + ": " + to_string(mapxsize) + ", " + to_string(mapysize) + "\n");
	}
	#endif
	// Note that the default "null" type for the coarse type is to have a density of 1 in every cell. "none" defaults to no
	// pristine map
	Matrix<float> toret;
	bCoarse = true;
	coarse_max = 0;
	if(fileinput == "null")
	{
		toret.SetSize(mapysize, mapxsize);
		coarse_map.SetSize(mapysize, mapxsize);
		for(unsigned long i = 0; i < mapxsize; i++)
		{
			for(unsigned long j = 0; j < mapysize; j++)
			{
				toret[j][i] = 1.0;
			}
		}
	}
	else if(fileinput == "none")
	{
		bCoarse = false;
	}
	else  // There is a map to read in.
	{
		toret.SetSize(mapysize, mapxsize);
		coarse_map.SetSize(mapysize, mapxsize);
		toret.import(fileinput);
	}
	if(bCoarse)
	{
		for(unsigned long i = 0; i < mapysize; i++)
		{
			for(unsigned long j = 0; j < mapxsize; j++)
			{
				coarse_map[i][j] = (unsigned long)(max((double)round(toret[i][j] * deme), 0.0));
				if(coarse_map[i][j] > coarse_max)
				{
					coarse_max = coarse_map[i][j];
				}
			}
		}
	}
	return;
}

void Map::calcPristineCoarseMap()
{
	//	os << "pristine coarse map file: " << mapvars.pristinecoarsemapfile << endl;
	string fileinput = mapvars.pristinecoarsemapfile;
	unsigned long mapxsize = mapvars.varcoarsemapxsize;
	unsigned long mapysize = mapvars.varcoarsemapysize;
	if(!check_set_dim)  // checks that the dimensions have been set.
	{
		throw FatalException("ERROR_MAP_003: dimensions not set.");
	}
	#ifndef SIZE_LIMIT
	if(mapxsize > 1000000 || mapysize > 1000000)
	{
		throw runtime_error("Extremely large map sizes set for " + fileinput + ": " + to_string(mapxsize) + ", " + to_string(mapysize) + "\n");
	}
	#endif
	// Note that the default "null" type for the coarse type is to have non-forest in every cell.
	Matrix<float> toret;
	pristine_coarse_max = 0;
	if(bCoarse)
	{
		if(fileinput == "null")
		{
			toret.SetSize(mapysize, mapxsize);
			pristine_coarse_map.SetSize(mapysize, mapxsize);
			for(unsigned long i = 0; i < mapxsize; i++)
			{
				for(unsigned long j = 0; j < mapysize; j++)
				{
					toret[j][i] = 1.0;
				}
			}
		}
		else if(fileinput == "none")
		{
			has_pristine = false;
		}
		else  // There is a map to read in.
		{
			toret.SetSize(mapysize, mapxsize);
			pristine_coarse_map.SetSize(mapysize, mapxsize);
			toret.import(fileinput);
		}
	}
	if(bCoarse && has_pristine)
	{
		for(unsigned long i = 0; i < mapysize; i++)
		{
			for(unsigned long j = 0; j < mapxsize; j++)
			{
				pristine_coarse_map[i][j] = (unsigned long)(max((double)round(toret[i][j] * deme), 0.0));
				if(pristine_coarse_map[i][j] > pristine_coarse_max)
				{
					pristine_coarse_max = pristine_coarse_map[i][j];
				}
			}
		}
	}
	return;
}

void Map::setTimeVars(double gen_since_pristine_in, double habitat_change_rate_in)
{
	update_time = 0;
	gen_since_pristine = gen_since_pristine_in;
	habitat_change_rate = habitat_change_rate_in;
}

void Map::calcOffset()
{
	if(mapvars.times_file != "null")
	{
		mapvars.setPristine(0);
	}
	//	os << mapvars.times_file << endl;
	if(fine_map.GetCols() == 0 || fine_map.GetRows() == 0)
	{
		throw FatalException("ERROR_MAP_004: fine map not set.");
	}
	if(coarse_map.GetCols() == 0 || coarse_map.GetRows() == 0)
	{
		if(bCoarse)
		{
			coarse_map.SetSize(fine_map.GetRows(), fine_map.GetCols());
		}
		//		throw FatalException("ERROR_MAP_004: coarse map not set.");
	}
	fine_x_offset = mapvars.varfinemapxoffset + mapvars.varsamplexoffset;
	fine_y_offset = mapvars.varfinemapyoffset + mapvars.varsampleyoffset;
	coarse_x_offset = mapvars.varcoarsemapxoffset;
	coarse_y_offset = mapvars.varcoarsemapyoffset;
	scale = mapvars.varcoarsemapscale;
	// this is the location of the top left (or north west) corner of the respective map
	// and the x and y distance from the top left of the grid object that contains the initial lineages.
	fine_x_min = -fine_x_offset;
	fine_y_min = -fine_y_offset;
	fine_x_max = fine_x_min + (fine_map.GetCols());
	fine_y_max = fine_y_min + (fine_map.GetRows());
	if(bCoarse) // Check if there is a coarse map
	{
		coarse_x_min = -coarse_x_offset - fine_x_offset;
		coarse_y_min = -coarse_y_offset - fine_y_offset;
		coarse_x_max = coarse_x_min + scale * (coarse_map.GetCols());
		coarse_y_max = coarse_y_min + scale * (coarse_map.GetRows());
	}
	else // Just set the offsets to the same as the fine map
	{
		coarse_x_min = fine_x_min;
		coarse_y_min = fine_y_min;
		coarse_x_max = fine_x_max;
		coarse_y_max = fine_y_max;
		scale = 1;
	}
	dispersal_relative_cost = mapvars.dispersal_relative_cost;
#ifdef DEBUG
	stringstream os;
	os << "finex: " << fine_x_min << "," << fine_x_max << endl;
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
	return;
}

void Map::validateMaps()
{
	stringstream os;
	os << "\rValidating maps..." << flush;
	double dTotal = fine_map.GetCols() + coarse_map.GetCols();
	unsigned long iCounter = 0;
	if(has_pristine)
	{
		if(fine_map.GetCols() == pristine_fine_map.GetCols() && fine_map.GetRows() == pristine_fine_map.GetRows() &&
		   coarse_map.GetCols() == pristine_coarse_map.GetCols() && coarse_map.GetRows() == pristine_coarse_map.GetRows())
		{
			os << "\rValidating maps...map sizes okay" << flush;
			writeInfo(os.str());
		}
		else
		{
			throw FatalException(
					"ERROR_MAP_009: Map validation failed - modern and pristine maps are not the same dimensions.");
		}
		for(unsigned long i = 0; i < fine_map.GetCols(); i++)
		{
			for(unsigned long j = 0; j < fine_map.GetRows(); j++)
			{
				if(fine_map[j][i] > pristine_fine_map[j][i])
				{
#ifdef DEBUG
					stringstream ss;
					ss << "fine map: " << fine_map[j][i] << " pristine map: " << pristine_fine_map[j][i];
					ss << " x,y: " << i << "," << j << endl;
					writeLog(50, ss);
#endif //DEBUG
					throw FatalException("ERROR_MAP_007: Map validation failed - fine map value larger "
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
	iCounter = fine_map.GetCols();
	if(has_pristine)
	{
		for(unsigned long i = 0; i < coarse_map.GetCols(); i++)
		{
			for(unsigned long j = 0; j < coarse_map.GetRows(); j++)
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
					throw FatalException("ERROR_MAP_008: Map validation failed - coarse map value larger "
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

void Map::updateMap(double generation)
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
				// pristine_fine_map = mapvars.pristinefinemapfile;
				// pristine_coarse_map = mapvars.pristinecoarsemapfile;
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

void Map::setLandscape(string landscape_type)
{
	if(landscape_type == "infinite")
	{
		writeInfo("Setting infinite landscape.\n");
		getValFunc = &Map::getValInfinite;
	}
	else if(landscape_type == "tiled_coarse")
	{
		writeInfo("Setting tiled coarse infinite landscape.\n");
		getValFunc = &Map::getValCoarseTiled;
	}
	else if(landscape_type == "tiled_fine")
	{
		writeInfo("Setting tiled fine infinite landscape.\n");
		getValFunc = &Map::getValFineTiled;
	}
	else if(landscape_type == "closed")
	{
		getValFunc = &Map::getValFinite;
	}
	else
	{
		throw FatalException("Provided landscape type is not a valid option: " + landscape_type);
	}
}

unsigned long Map::getVal(const double& x, const double& y,
						  const long& xwrap, const long& ywrap, const double& current_generation)
{
	return (this->*getValFunc)(x, y, xwrap, ywrap, current_generation);
}

unsigned long Map::getValInfinite(
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

unsigned long Map::getValCoarseTiled(
	const double& x, const double& y, const long& xwrap, const long& ywrap, const double& current_generation)
{
	double newx = fmod(x + (xwrap * x_dim) + fine_x_offset + coarse_x_offset, coarse_map.GetCols());
	double newy = fmod(y + (ywrap * y_dim) + fine_x_offset + coarse_x_offset, coarse_map.GetRows());
	if(newx < 0)
	{
		newx += coarse_map.GetCols();
	}
	if(newy < 0)
	{
		newy += coarse_map.GetRows();
	}
	return getValCoarse(newx, newy, current_generation);
}

unsigned long Map::getValFineTiled(
	const double& x, const double& y, const long& xwrap, const long& ywrap, const double& current_generation)
{

	double newx = fmod(x + (xwrap * x_dim) + fine_x_offset, fine_map.GetCols());
	double newy = fmod(y + (ywrap * y_dim) + fine_y_offset, fine_map.GetRows());
	// Now adjust for incorrect wrapping behaviour of fmod
	if(newx < 0)
	{
		newx += fine_map.GetCols();
	}
	if(newy < 0)
	{
		newy += fine_map.GetRows();
	}
#ifdef DEBUG
	if(newx >= fine_map.GetCols() || newx < 0 || newy >= fine_map.GetRows() || newy < 0)
	{
		stringstream ss;
		ss << "Fine map indexing out of range of fine map." << endl;
		ss << "x, y: " << newx << ", " << newy << endl;
		ss << "cols, rows: " << fine_map.GetCols() << ", " << fine_map.GetRows() << endl;
		throw out_of_range(ss.str());
	}
#endif
	return getValFine(newx, newy, current_generation);
}

unsigned long Map::getValCoarse(const double &xval, const double &yval, const double &current_generation)
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

unsigned long Map::getValFine(const double&xval, const double &yval, const double& current_generation)
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

unsigned long Map::getValFinite(
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
	   yval >= fine_y_max) && bCoarse)  // check if the coordinate comes from the coarse resolution map.
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

unsigned long Map::convertSampleXToFineX(const unsigned long &x, const long &xwrap)
{
	return x + fine_x_offset + (xwrap * x_dim);
}

unsigned long Map::convertSampleYToFineY(const unsigned long &y, const long &ywrap)
{
	return y + fine_y_offset + (ywrap * y_dim);
}

void Map::convertFineToSample(long & x, long & xwrap, long &y, long &ywrap)
{
	auto tmpx = double(x);
	auto tmpy = double(y);
	convertCoordinates(tmpx, tmpy, xwrap, ywrap);
	x = floor(tmpx);
	y = floor(tmpy);
}


unsigned long Map::getInitialCount(double dSample, DataMask& samplemask)
{
	unsigned long toret;
	toret = 0;
	long x, y;
	long xwrap, ywrap;
	unsigned long max_x, max_y;
	if(samplemask.getDefault())
	{
		max_x = fine_map.GetCols();
		max_y = fine_map.GetRows();
	}
	else
	{
		max_x = samplemask.sample_mask.GetCols();
		max_y = samplemask.sample_mask.GetRows();
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

SimParameters Map::getSimParameters()
{
	return mapvars;
}

bool Map::checkMap(const double& x, const double& y, const long& xwrap, const long& ywrap, const double generation)
{
	return getVal(x, y, xwrap, ywrap, generation) != 0;
}

bool Map::checkFine(const double& x, const double& y, const long& xwrap, const long& ywrap)
{
	double tmpx, tmpy;
	tmpx = x + xwrap * x_dim;
	tmpy = y + ywrap * y_dim;
	return !(tmpx < fine_x_min || tmpx >= fine_x_max || tmpy < fine_y_min || tmpy >= fine_y_max);
}

void Map::convertCoordinates(double& x, double& y, long& xwrap, long& ywrap)
{
	xwrap += floor(x / x_dim);
	ywrap += floor(y / y_dim);
	x = x - xwrap * x_dim;
	y = y - ywrap * y_dim;
}

unsigned long Map::runDispersal(const double& dist,
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

void Map::clearMap()
{
	current_map_time = 0;
	check_set_dim = false;
	is_pristine = false;
}

string Map::printVars()
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

unsigned long Map::getHabitatMax()
{
	return habitat_max;
}

void Map::recalculateHabitatMax()
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
