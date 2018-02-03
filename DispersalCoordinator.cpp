//This file is part of NECSim project which is released under BSD-3 license.
//See file **LICENSE.txt** or visit https://opensource.org/licenses/BSD-3-Clause) for full license details.

/**
 * @author Samuel Thompson
 * @date 07/08/2017
 * @file DispersalCoordinator.cpp
 * @brief Contains the DispersalCoordinator, which contains all routines related to dispersal
 * including utilisation of density maps and dispersal probability maps.
 * 
 * @copyright <a href="https://opensource.org/licenses/BSD-3-Clause">BSD-3 Licence.</a>
 */
#include "DispersalCoordinator.h"
#include "CustomExceptions.h"

DispersalCoordinator::DispersalCoordinator()
{
}

DispersalCoordinator::~DispersalCoordinator()
{
	
}

void DispersalCoordinator::setRandomNumber(NRrand * NR_ptr)
{
	NR = NR_ptr;
}


void DispersalCoordinator::setHabitatMap(Map *map_ptr)
{
	habitat_map = map_ptr;
	xdim = habitat_map->getSimParameters().varfinemapxsize;
}

void DispersalCoordinator::setGenerationPtr(double * generation_ptr)
{
	generation = generation_ptr;
}


void DispersalCoordinator::setDispersal(const string &dispersal_method, const string &dispersal_file,
							  const unsigned long dispersal_x, const unsigned long dispersal_y,
							  const double &m_probin, const double &cutoffin,
							  const double &sigmain, const double &tauin, const bool &restrict_self)
{
	// Open our file connection
	if(dispersal_file == "none")
	{
		setEndPointFptr(restrict_self);
		NR->setDispersalParams(sigmain, tauin);
		NR->setDispersalMethod(dispersal_method, m_probin, cutoffin);
		doDispersal = &DispersalCoordinator::disperseDensityMap;
	}
	else if(dispersal_file == "null")
	{
		doDispersal = &DispersalCoordinator::disperseNullDispersalMap;
	}
	else
	{
		doDispersal = &DispersalCoordinator::disperseDispersalMap;
		// Check file existance
		ifstream infile(dispersal_file);
		if(!infile.good())
		{
			string msg = "Could not access dispersal map file " + dispersal_file + ". Check file exists and is readable.";
			throw FatalException(msg);
		}
		infile.close();
		dispersal_prob_map.SetSize(dispersal_x * dispersal_y, dispersal_x * dispersal_y);
		dispersal_prob_map.import(dispersal_file);
	}
}

void DispersalCoordinator::disperseNullDispersalMap(Step &this_step)
{
	// Pick a random cell - that's all we need
	this_step.oldx = floor(NR->d01()*(xdim-1));
	this_step.oldy = floor(NR->d01()*(xdim - 1));
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
	unsigned long max_col = dispersal_prob_map.GetCols() - 1;
	while(max_col - min_col > 1)
	{
		unsigned long to_check = floor(double(max_col-min_col)/2.0) + min_col;
		if(dispersal_prob_map[row_ref][to_check] > random_no)
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
}

void DispersalCoordinator::calculateCellCoordinates(Step & this_step, const unsigned long &col_ref)
{
	this_step.oldx = long(floor(fmod(double(col_ref), xdim)));
	this_step.oldy = long(floor(double(col_ref)/xdim));
	this_step.oldxwrap = 0;
	this_step.oldywrap = 0;
	// Convert back to sample map
	habitat_map->convertFineToSample(this_step.oldx, this_step.oldxwrap, this_step.oldy, this_step.oldywrap);
	
}

unsigned long DispersalCoordinator::calculateCellReference(Step &this_step)
{
	unsigned long x = habitat_map->convertSampleXToFineX(this_step.oldx, this_step.oldxwrap);
	unsigned long y = habitat_map->convertSampleYToFineY(this_step.oldy, this_step.oldywrap);
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
		density = habitat_map->runDispersal(dist, angle, this_step.oldx,
										  this_step.oldy, this_step.oldxwrap, this_step.oldywrap, fail, *generation);
		if(!fail)
		{
			fail = !checkEndPoint(density, this_step.oldx, this_step.oldy, this_step.oldxwrap, this_step.oldywrap,
								  startx, starty, startxwrap, startywrap);
		}
		// Discard the dispersal event a percentage of the time, based on the maximum value of the habitat map.
		// This is to correctly mimic less-dense cells having a lower likelihood of being the parent to the cell.
		
#ifdef DEBUG
		if(habitat_map->getVal(this_step.oldx, this_step.oldy, this_step.oldxwrap, this_step.oldywrap, *generation) == 0 &&
			!fail)
		{
			stringstream ss;
			ss << "x,y: " << this_step.oldx << "," << this_step.oldy;
			ss << " x,y wrap: " << this_step.oldxwrap << "," << this_step.oldywrap << "Habitat cover: ";
			ss << habitat_map->getVal(this_step.oldx, this_step.oldy, this_step.oldxwrap,
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


bool DispersalCoordinator::checkEndPoint(const unsigned long & density, long &oldx, long &oldy,
										 long &oldxwrap, long &oldywrap, const long &startx, const long &starty,
										 const long &startxwrap, const long &startywrap)
{
	return (this->*checkEndPointFptr)(density, oldx, oldy, oldxwrap, oldywrap, startx, starty, startxwrap, startywrap);
}

bool DispersalCoordinator::checkEndPointDensity(const unsigned long &density, long &oldx, long &oldy,
												long &oldxwrap, long &oldywrap, const long &startx, const long &starty,
												const long &startxwrap, const long &startywrap)
{
	if((double(density) / double(habitat_map->getHabitatMax())) <
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

bool DispersalCoordinator::checkEndPointRestricted(const unsigned long &density, long &oldx, long &oldy, long &oldxwrap, long &oldywrap, const long &startx, const long &starty,
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



