#include "DataMask.h"
#include "Map.h"
#include "Logging.h"

/**
 * @author Samuel Thompson
 * @file DataMask.cpp
 * @brief  Contains the DataMask class for describing the spatial sampling pattern on a landscape.
 *
 * @copyright <a href="https://opensource.org/licenses/BSD-3-Clause">BSD-3 Licence.</a>
 */

DataMask::DataMask()
{
	bDefault = true;
	x_dim = 0;
	y_dim = 0;
	x_offset = 0;
	y_offset = 0;
	mask_x_dim = 0;
	mask_y_dim = 0;
	getProportionfptr = &DataMask::getBoolProportion;
}

bool DataMask::getDefault()
{
	return bDefault;
}

bool DataMask::setup(const string &sample_mask_file, const unsigned long &x_in, const unsigned long &y_in,
					 const unsigned long &mask_x_in, const unsigned long &mask_y_in,
					 const unsigned long &x_offset_in, const unsigned long &y_offset_in)
{
#ifdef DEBUG
	if((x_in > mask_x_in || y_in > mask_y_in) && !bDefault)
	{
		writeLog(50, "Grid size: " + to_string(x_in) + ", " + to_string(y_in));
		writeLog(50, "Sample mask size: " + to_string(mask_x_in) + ", " + to_string(mask_y_in));
		throw FatalException("Datamask dimensions do not make sense");
	}
#endif // DEBUG
	inputfile = sample_mask_file;
	x_dim = x_in;
	y_dim = y_in;
	mask_x_dim = mask_x_in;
	mask_y_dim  = mask_y_in;
	x_offset = x_offset_in;
	y_offset = y_offset_in;
	if(inputfile == "null" || inputfile == "none")
	{
		bDefault = true;
	}
	else
	{
		bDefault = false;
	}
	return bDefault;
}

void DataMask::importBooleanMask(unsigned long xdim, unsigned long ydim, unsigned long mask_xdim,
								 unsigned long mask_ydim,
								 unsigned long xoffset, unsigned long yoffset, string inputfile)
{
	if(!setup(inputfile, xdim, ydim, mask_xdim, mask_ydim, xoffset, yoffset))
	{
		doImport();
	}
}
void DataMask::doImport()
{
	sample_mask.SetSize(mask_y_dim, mask_x_dim);
	sample_mask.import(inputfile);
	getProportionfptr = &DataMask::getBoolProportion;

}

void DataMask::importSampleMask(SimParameters &mapvarin)
{
	if(!setup(mapvarin.samplemaskfile, mapvarin.vargridxsize, mapvarin.vargridysize,
			 mapvarin.varsamplexsize, mapvarin.varsampleysize, mapvarin.varsamplexoffset, mapvarin.varsampleyoffset))
	{
		if(mapvarin.uses_spatial_sampling)
		{
#ifdef DEBUG
			writeLog(10, "Using spatial sampling.");
			writeLog(10, "Mask dimensions: " + to_string(mask_x_dim) + ", " + to_string(mask_y_dim));
#endif // DEBUG
			sample_mask_exact.SetSize(mask_y_dim, mask_x_dim);
			sample_mask_exact.import(inputfile);
			getProportionfptr = &DataMask::getSampleProportion;
		}
		else
		{
			doImport();
		}
	}
	else
	{
		if(mapvarin.uses_spatial_sampling)
		{
			// This could perhaps be a warning, but I'd prefer to have the warning/prohibit potential in python
			// and throw a full exception here.
			throw FatalException("Cannot use a spatial sampling routine when the map file is null.");
		}
		getProportionfptr = &DataMask::getNullProportion;
	}
}

bool DataMask::getVal(const long &x, const long &y, const long &xwrap, const long &ywrap)
{
	long xval = x + (xwrap * x_dim) + x_offset;
	long yval = y + (ywrap * y_dim) + y_offset;
	if(bDefault)
	{
		return true;
	}
#ifdef DEBUG
	if(xval < 0 || xval >= (long) mask_x_dim || yval < 0 || yval >= (long) mask_y_dim)
	{
		stringstream ss;
		ss << "Get value on samplemask requested for non index." << endl;
		ss << "x, y: " << x << ", " << y << endl;
		ss << "dimensions x,y: " << mask_x_dim << ", " << mask_y_dim << endl;
		ss << "x, y wrap: " << xwrap << ", " << ywrap << endl;
		ss << "xval, yval: " << xval << ", " << yval << endl;
		ss << "offsets x, y: " << x_offset << ", " << y_offset << endl;
		writeLog(50, ss);
		ss.str("Get value on samplemask requested for non index.");
		throw out_of_range(ss.str());
	}
#endif
	return sample_mask[yval][xval];
}


double DataMask::getNullProportion(const long &x, const long &y, const long &xwrap, const long &ywrap)
{
	return 1.0;
}

double DataMask::getBoolProportion(const long &x, const long &y, const long &xwrap, const long &ywrap)
{

	if(getVal(x, y, xwrap, ywrap))
	{
		return 1.0;
	}
	else
	{
		return 0.0;
	}
}

double DataMask::getSampleProportion(const long &x, const long &y, const long &xwrap, const long &ywrap)
{
#ifdef DEBUG
	if(bDefault || sample_mask_exact.GetCols() == 0)
	{
		throw out_of_range("Cannot get the exact value from a samplemask if we are using a null mask, or the "
								   "exact samplemask has not been properly imported.");
	}
#endif // DEBUG
	long xval = x + (xwrap * x_dim) + x_offset;
	long yval = y + (ywrap * y_dim) + y_offset;
	return sample_mask_exact[yval][xval];
}

double DataMask::getExactValue(const long &x, const long &y, const long &xwrap, const long &ywrap)
{
	return (this->*getProportionfptr)(x, y, xwrap, ywrap);
}

void DataMask::convertBoolean(Map &map1, const double &deme_sampling, const double &generation)
{
	// Clear the old boolean object and set the new size
	sample_mask.SetSize(y_dim, x_dim);
	for(unsigned long y = 0; y < y_dim; y++)
	{
		for(unsigned long x = 0; x < x_dim; x++)
		{
			long tmp_x = x;
			long tmp_y = y;
			long tmp_xwrap = 0;
			long tmp_ywrap = 0;
			recalculate_coordinates(tmp_x, tmp_y, tmp_xwrap, tmp_ywrap);
			double density = map1.getVal(tmp_x, tmp_y, tmp_xwrap, tmp_ywrap, generation) * deme_sampling;
			sample_mask[y][x] = density >= 1.0;
		}
	}
}

void DataMask::clearSpatialMask()
{
	sample_mask_exact.SetSize(0, 0);
}

void DataMask::recalculate_coordinates(long &x, long &y, long &x_wrap, long &y_wrap)
{
	if(!bDefault)
	{
		x_wrap = (long)((floor((x - (double) x_offset) / (double) x_dim)));
		y_wrap = (long)((floor((y - (double) y_offset) / (double) y_dim)));
		x += -x_offset - (x_wrap * x_dim);
		y += -y_offset - (y_wrap * y_dim);
	}
}
