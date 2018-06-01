//This file is part of NECSim project which is released under MIT license.
//See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.

/**
 * @author Samuel Thompson
 * @date 07/08/2017
 * @file DispersalCoordinator.cpp
 * @brief Contains the DispersalCoordinator, which contains all routines related to dispersal
 * including utilisation of density maps and dispersal probability maps.
 * 
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 */
#include "DispersalCoordinator.h"
#include "CustomExceptions.h"
#include "Logger.h"

DispersalCoordinator::DispersalCoordinator()
{
}

DispersalCoordinator::~DispersalCoordinator()
{

}

void DispersalCoordinator::setRandomNumber(NRrand *NR_ptr)
{
	NR = NR_ptr;
}

void DispersalCoordinator::setHabitatMap(Landscape *map_ptr)
{
	landscape = map_ptr;
	if(!map_ptr)
	{
		throw FatalException("Attempting to set map pointer to null pointer in DispersalCoordinator.");
	}
	xdim = landscape->getSimParameters()->fine_map_x_size;
}

void DispersalCoordinator::setGenerationPtr(double *generation_ptr)
{
	generation = generation_ptr;
}

void DispersalCoordinator::setDispersal(const string &dispersal_method, const string &dispersal_file,
										const unsigned long &dispersal_x, const unsigned long &dispersal_y,
										const double &m_probin, const double &cutoffin,
										const double &sigmain, const double &tauin, const bool &restrict_self)
{
	// Open our file connection
	if(dispersal_file == "none")
	{
		if(!NR)
		{
			throw FatalException("Random number generator pointer has not been set in DispersalCoordinator.");
		}
		writeInfo("Using dispersal kernel.\n");
		setEndPointFptr(restrict_self);
		NR->setDispersalParams(sigmain, tauin);
		NR->setDispersalMethod(dispersal_method, m_probin, cutoffin);
		doDispersal = &DispersalCoordinator::disperseDensityMap;
	}
	else if(dispersal_file == "null")
	{
		writeInfo("Using null dispersal file.\n");
		doDispersal = &DispersalCoordinator::disperseNullDispersalMap;
	}
	else
	{
		writeInfo("Using dispersal file.\n");
		doDispersal = &DispersalCoordinator::disperseDispersalMap;
		// Check file existance
		ifstream infile(dispersal_file);
		if(!infile.good())
		{
			string msg =
					"Could not access dispersal map file " + dispersal_file + ". Check file exists and is readable.";
			throw FatalException(msg);
		}
		infile.close();
		dispersal_prob_map.setSize(dispersal_x * dispersal_y, dispersal_x * dispersal_y);
		dispersal_prob_map.import(dispersal_file);
		fixDispersal();
		dispersal_prob_map.close();
	}
}

void DispersalCoordinator::setDispersal(SimParameters *simParameters)
{
	if(!simParameters)
	{
		throw FatalException("Simulation parameters pointer has not been set for DispersalCoordinator.");
	}
	setDispersal(simParameters->dispersal_method, simParameters->dispersal_file,
				 simParameters->fine_map_x_size, simParameters->fine_map_y_size, simParameters->m_prob,
				 simParameters->cutoff, simParameters->sigma, simParameters->tau, simParameters->restrict_self);
}

void DispersalCoordinator::fixDispersal()
{
	for(unsigned long row = 0; row < dispersal_prob_map.getRows(); row++)
	{
		fixDispersalRow(row);
	}
}

void DispersalCoordinator::fixDispersalRow(unsigned long row)
{
	// check if the row needs to be fixed
	if(checkDispersalRow(row))
	{
		double total_value = 0.0;
		for(unsigned long i = 0; i < dispersal_prob_map.getCols(); i++)
		{
			total_value += dispersal_prob_map[row][i];
		}
		if(total_value == 0.0)
		{
			return;
		}
		dispersal_prob_map[row][0] = dispersal_prob_map[row][0] / total_value;
		for(unsigned long i = 1; i < dispersal_prob_map.getCols(); i++)
		{
			dispersal_prob_map[row][i] = dispersal_prob_map[row][i - 1] + (dispersal_prob_map[row][i] / total_value);
		}
#ifdef DEBUG
		if(checkDispersalRow(row))
		{
			throw FatalException("Dispersal probability map not correctly fixed to sum to 1.0.");
		}
#endif // DEBUG
	}
}

bool DispersalCoordinator::checkDispersalRow(unsigned long row)
{
	if(abs(dispersal_prob_map[row][dispersal_prob_map.getCols() - 1] - 1.0) > 0.00000001)
	{
		return true;
	}
	for(unsigned long i = 0; i < dispersal_prob_map.getCols() - 1; i++)
	{
		if(dispersal_prob_map[row][i] > dispersal_prob_map[row][i + 1])
		{
			return true;
		}
	}
	return false;
}

void DispersalCoordinator::verifyDispersalMap()
{
	if(dispersal_prob_map.getCols() > 0)
	{
		writeInfo("Verifying dispersal...\n");
		if(dispersal_prob_map.getCols() != dispersal_prob_map.getRows())
		{
			throw FatalException("Dispersal probability map dimensions do not match.");
		}
		bool has_printed = false;
		for(unsigned long y = 0; y < dispersal_prob_map.getRows(); y++)
		{
			Step origin_step;
			calculateCellCoordinates(origin_step, y);
#ifdef DEBUG
			assertReferenceMatches(y);
#endif // DEBUG
			bool origin_value = landscape->getVal(origin_step.oldx, origin_step.oldy,
												  origin_step.oldxwrap, origin_step.oldywrap, 0.0) > 0;
			double dispersal_total = 0.0;
			for(unsigned long x = 0; x < dispersal_prob_map.getCols(); x++)
			{
				Step destination_step;
				calculateCellCoordinates(destination_step, x);
				assertReferenceMatches(x);
				bool destination_value = landscape->getVal(destination_step.oldx, destination_step.oldy,
														   destination_step.oldxwrap,
														   destination_step.oldywrap, 0.0) > 0;
				double dispersal_prob = 0.0;
				if(x == 0)
				{
					dispersal_prob = dispersal_prob_map[y][0];
				}
				else
				{
					dispersal_prob = dispersal_prob_map[y][x] - dispersal_prob_map[y][x - 1];
				}
				dispersal_total += dispersal_prob;
				if(dispersal_prob > 0.0)
				{
					if(!destination_value && origin_value)
					{
						stringstream ss;
						ss << "Dispersal from " << origin_step.oldx << ", " << origin_step.oldy << " (";
						ss << origin_step.oldxwrap << ", " << origin_step.oldywrap << ") to ";
						ss << destination_step.oldx << ", " << destination_step.oldy << " ("
						   << destination_step.oldxwrap;
						ss << ", " << destination_step.oldywrap << ")" << endl;
						ss << "Source row: " << y << " destination row: " << x << endl;
						ss << "Dispersal map value: " << dispersal_prob << endl;
						ss << "Origin density: " << landscape->getVal(origin_step.oldx, origin_step.oldy,
																	  origin_step.oldxwrap, origin_step.oldywrap,
																	  0.0) << endl;
						ss << "Destination density: " << landscape->getVal(destination_step.oldx, destination_step.oldy,
																		   destination_step.oldxwrap,
																		   destination_step.oldywrap, 0.0) << endl;
						writeError(ss.str());
						throw FatalException("Dispersal map is non zero where density is 0.");
					}
					if(!origin_value && !has_printed)
					{
						has_printed = true;
						writeWarning("Dispersal values exist for non-zero density values.");
					}
				}
			}
			if(dispersal_total == 0.0 && origin_value)
			{
				stringstream ss;
				ss << "No dispersal probabilities from cell at " << origin_step.oldx << ", " << origin_step.oldy;
				ss << " (" << origin_step.oldxwrap << ", " << origin_step.oldywrap;
				ss << ") to any other cell, despite non-zero density." << endl;
				throw FatalException(ss.str());
			}
		}
	}
}

void DispersalCoordinator::assertReferenceMatches(unsigned long expected)
{
	unsigned long row_ref = expected;
	Step step;
	calculateCellCoordinates(step, row_ref);
	auto actual = calculateCellReference(step);
	if(actual != expected)
	{
		stringstream ss;
		ss << "Expected reference " << expected << endl;
		ss << "Actual reference " << actual << endl;
		ss << "Coordinates: " << step.oldx << ", " << step.oldy << "(" << step.oldxwrap << ", ";
		ss << step.oldywrap << ")" << endl;
		ss << "Converted values: " << landscape->convertSampleXToFineX(step.oldx, step.oldxwrap);
		ss << ", " << landscape->convertSampleYToFineY(step.oldy, step.oldywrap) << endl;
		throw FatalException(ss.str());
	}
}

void DispersalCoordinator::disperseNullDispersalMap(Step &this_step)
{
	// Pick a random cell - that's all we need
	this_step.oldx = floor(NR->d01() * (xdim - 1));
	this_step.oldy = floor(NR->d01() * (xdim - 1));
}

void DispersalCoordinator::disperseDispersalMap(Step &this_step)
{
	// Generate random number 0-1
	double random_no = NR->d01();
	// Now find the cell with that value
	// Now we get the cell reference
	unsigned long row_ref = calculateCellReference(this_step);
	// Interval bisection on the cells to get the dispersal value
	unsigned long min_col = 0;
	unsigned long max_col = dispersal_prob_map.getCols() - 1;
	while(max_col - min_col > 1)
	{
		auto to_check = static_cast<unsigned long>(floor(double(max_col - min_col) / 2.0) + min_col);
		if(dispersal_prob_map[row_ref][to_check] < random_no)
		{
			min_col = to_check;
		}
		else
		{
			max_col = to_check;
		}
	}
	// Now get the coordinates of our cell reference
	calculateCellCoordinates(this_step, max_col);
#ifdef DEBUG
	if(landscape->getVal(this_step.oldx, this_step.oldy, this_step.oldxwrap, this_step.oldywrap, *generation) < 1.0)
	{
		stringstream ss;
		ss << "Dispersal attempted to cell of zero density " << this_step.oldx << ", " << this_step.oldy;
		ss << " (" << this_step.oldxwrap << ", " << this_step.oldywrap << ")" << endl;
		throw FatalException(ss.str());
	}
#endif // DEBUG
}

void DispersalCoordinator::calculateCellCoordinates(Step &this_step, const unsigned long &col_ref)
{
	this_step.oldx = long(floor(fmod(double(col_ref), xdim)));
	this_step.oldy = long(floor(double(col_ref) / xdim));
	this_step.oldxwrap = 0;
	this_step.oldywrap = 0;
	// Convert back to sample map
	landscape->convertFineToSample(this_step.oldx, this_step.oldxwrap, this_step.oldy, this_step.oldywrap);

}

unsigned long DispersalCoordinator::calculateCellReference(Step &this_step)
{
	unsigned long x = landscape->convertSampleXToFineX(this_step.oldx, this_step.oldxwrap);
	unsigned long y = landscape->convertSampleYToFineY(this_step.oldy, this_step.oldywrap);
	return x + (y * xdim);
}

void DispersalCoordinator::disperseDensityMap(Step &this_step)
{
	bool fail;
	fail = true;
	// Store the starting positions
	long startx, starty, startxwrap, startywrap;
	startx = this_step.oldx;
	starty = this_step.oldy;
	startxwrap = this_step.oldxwrap;
	startywrap = this_step.oldywrap;
	// keep looping until we reach a viable place to move from.
	// Store the density in the end location.
	unsigned long density;
	double dist, angle;
	while(fail)
	{
		angle = NR->direction();
		dist = NR->dispersal();
		density = landscape->runDispersal(dist, angle, this_step.oldx,
										  this_step.oldy, this_step.oldxwrap, this_step.oldywrap, fail, *generation);
		if(!fail)
		{
			fail = !checkEndPoint(density, this_step.oldx, this_step.oldy, this_step.oldxwrap, this_step.oldywrap,
								  startx, starty, startxwrap, startywrap);
		}
		// Discard the dispersal event a percentage of the time, based on the maximum value of the habitat map.
		// This is to correctly mimic less-dense cells having a lower likelihood of being the parent to the cell.

#ifdef DEBUG
		if(landscape->getVal(this_step.oldx, this_step.oldy, this_step.oldxwrap, this_step.oldywrap, *generation) ==
		   0 &&
		   !fail)
		{
			stringstream ss;
			ss << "x,y: " << this_step.oldx << "," << this_step.oldy;
			ss << " x,y wrap: " << this_step.oldxwrap << "," << this_step.oldywrap << "Habitat cover: ";
			ss << landscape->getVal(this_step.oldx, this_step.oldy, this_step.oldxwrap,
									this_step.oldywrap, *generation) << endl;
			writeLog(50, ss);
			throw FatalException("ERROR_MOVE_007: Dispersal attempted to non-habitat. Check dispersal function.");
		}
#endif
	}
}

void DispersalCoordinator::setEndPointFptr(const bool &restrict_self)
{
	if(restrict_self)
	{
		checkEndPointFptr = &DispersalCoordinator::checkEndPointRestricted;
	}
	else
	{
		checkEndPointFptr = &DispersalCoordinator::checkEndPointDensity;
	}
}

bool DispersalCoordinator::checkEndPoint(const unsigned long &density, long &oldx, long &oldy,
										 long &oldxwrap, long &oldywrap, const long &startx, const long &starty,
										 const long &startxwrap, const long &startywrap)
{
	return (this->*checkEndPointFptr)(density, oldx, oldy, oldxwrap, oldywrap, startx, starty, startxwrap, startywrap);
}

bool DispersalCoordinator::checkEndPointDensity(const unsigned long &density, long &oldx, long &oldy,
												long &oldxwrap, long &oldywrap, const long &startx, const long &starty,
												const long &startxwrap, const long &startywrap)
{
	if((double(density) / double(landscape->getHabitatMax())) <
	   NR->d01())
	{
		oldx = startx;
		oldy = starty;
		oldxwrap = startxwrap;
		oldywrap = startywrap;
		return false;
	}
	return true;
}

bool DispersalCoordinator::checkEndPointRestricted(const unsigned long &density, long &oldx, long &oldy, long &oldxwrap,
												   long &oldywrap, const long &startx, const long &starty,
												   const long &startxwrap, const long &startywrap)
{
	if(startx == oldx && starty == oldy && startxwrap == oldxwrap && startywrap == oldywrap)
	{
		return false;
	}
	return checkEndPointDensity(density, oldx, oldy, oldxwrap, oldywrap, startx, starty, startxwrap, startywrap);
}

void DispersalCoordinator::disperse(Step &this_step)
{
	(this->*doDispersal)(this_step);
}



