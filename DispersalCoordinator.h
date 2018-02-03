//This file is part of NECSim project which is released under BSD-3 license.
//See file **LICENSE.txt** or visit https://opensource.org/licenses/BSD-3-Clause) for full license details.

/**
 * @author Samuel Thompson
 * @date 07/08/2017
 * @file DispersalCoordinator.h
 * @brief Contains the DispersalCoordinator, which contains all routines related to dispersal
 * including utilisation of density maps and dispersal probability maps.
 * 
 * @copyright <a href="https://opensource.org/licenses/BSD-3-Clause">BSD-3 Licence.</a>
 */
#ifndef DISPERSALCOORDINATOR_H
#define DISPERSALCOORDINATOR_H

#include <cstring>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <cmath>

#include "NRrand.h"
#include "Matrix.h"
#include "Step.h"
#include "Map.h"

/**
 * @class DispersalCoordinator
 * @brief Class for generating dispersal distances and provide routines for reading dispersal distance maps
 * as a unwound map-of-maps. This class also handles reading density maps for rejection sampling.
 * 
 * It requires linking to a density map, random number generator and a generation counter from the Tree class.
 * 
 * Note that no element of this object is recorded during a paused simulation, as all objects pointed to are stored 
 * elsewhere and behaviours are recalculated upon simulation resume.
 */
class DispersalCoordinator
{
protected:
	
	// Our map of dispersal probabilities (if required)
	// This will contain cummulative probabilities across rows
	Matrix<double> dispersal_prob_map;
	// Our random number generator for dispersal distances
	// This is a pointer so that the random number generator is the same
	// across the program.
	
	NRrand * NR;
	// Pointer to the habitat map object for getting density values.
	Map * habitat_map;
	// Pointer to the generation counter for the simulation
	double * generation;
	
	// function ptr for our getDispersal function
	typedef void (DispersalCoordinator::*dispersal_fptr)(Step &this_step);
	dispersal_fptr doDispersal;
	
	// Function pointer for end checks
	typedef bool (DispersalCoordinator::*end_fptr)(const unsigned long &density, long &oldx, long &oldy, long &oldxwrap,
										 long &oldywrap, const long &startx, const long &starty, 
										 const long &startxwrap, const long &startywrap); 
	// once setup will contain the end check function to use for this simulation.
	end_fptr checkEndPointFptr;
	unsigned long xdim;
	
public:
	DispersalCoordinator();
	
	~DispersalCoordinator();
	
	/**
	 * @brief Sets the random number pointer to an NRrand instance.
	 * @param NR_ptr the random number object to set to
	 */
	void setRandomNumber(NRrand * NR_ptr);
	
	/**
	 * @brief Sets the pointer to the Map object
	 * @param map_ptr pointer to a Map object
	 */
	void setHabitatMap(Map *map_ptr);
	
	/**
	 * @brief Sets the generation pointer to the provided double
	 * @param generation_ptr pointer to the generation double
	 */
	void setGenerationPtr(double * generation_ptr);


	/**
	 * @brief Sets the dispersal method and parameters
	 * 
	 * @param dispersal_method string containing the dispersal type. Can be one of [normal, fat-tail, norm-uniform]
	 * @param dispersal_file string containing the dispersal file, or "none" if using dispersal kernel
	 * @param m_probin the probability of drawing from the uniform distribution. Only relevant for uniform dispersals
	 * @param cutoffin the maximum value to be drawn from the uniform dispersal. Only relevant for uniform dispersals
	 * @param sigmain the fatness of the fat-tailed dispersal kernel
	 * @param tauin the width of the fat-tailed dispersal kernel
	 * @param restrict_self if true, denies possibility that dispersal comes from the same cell as the parent
	 */
	void setDispersal(const string &dispersal_method, const string &dispersal_file,
					  const unsigned long dispersal_x, const unsigned long dispersal_y,
					  const double &m_probin, const double &cutoffin,
					  const double &sigmain, const double &tauin, const bool &restrict_self);
	
	/**
	 * @brief Picks a random cell from the whole map and stores the value in the step object
	 * @param this_step the step object to store end points in
	 */
	void disperseNullDispersalMap(Step &this_step);
	
	/**
	 * @brief Picks a random dispersal distance from the dispersal map
	 * @param this_step the step object to store end points in
	 */
	void disperseDispersalMap(Step &this_step);
	
	
	/**
	 * @brief Calculates the new coordinates for a column reference.
	 * This includes converting between the fine map and sample map.
	 * New coordinates are saved in this_step.
	 * @param this_step the step to save new coordinates in.
	 * @param col_ref the column reference for 
	 */
	void calculateCellCoordinates(Step & this_step, const unsigned long &col_ref);
	
	/**
	 * @brief Calculates the cell reference for a particular coordinate
	 * 
	 * The formula for this calculation is x + (y * xdim) where xdim is the dimensions of the fine map, 
	 * and x and y are the coordinates for the fine map
	 * 
	 * @param this_step the step object containing the x, y location, and x,y wrapping
	 * @return the cell reference from the dispersal_prob_map which corresponds to the required cell
	 */
	unsigned long calculateCellReference(Step &this_step);
	
	/**
	 * @brief Calls the dispersal kernel from the supplied dispersal distribution.
	 * @param this_step the step object to store end points in
	 */
	void disperseDensityMap(Step &this_step);
	
	/**
	 * @brief Sets the end point function pointer correctly, based on whether it is restricted or not.
	 * @param restrict_self if true, denies possibility that dispersal comes from the same cell as the parent
	 */
	void setEndPointFptr(const bool &restrict_self);
	
	/**
	 * @brief Check the end point for the given coordinates and density
	 * @param density the density at the end point - avoids an extra call to Map::getVal()
	 * @param oldx the old x position
	 * @param oldy the old y position
	 * @param oldxwrap the old x wrap
	 * @param oldywrap the old y wrap
	 * @param startx the starting x position
	 * @param starty the starting y position
	 * @param startxwrap the starting x wrap
	 * @param startywrap the ending y wrap
	 * @return true if the end point passes the density and restricted checks
	 */
	bool checkEndPoint(const unsigned long & density, long &oldx, long &oldy, long &oldxwrap, long &oldywrap,
					    const long &startx, const long &starty, const long &startxwrap, const long &startywrap);
	
	
	/**
	 * @brief Check the end point for the given coordinates and density
	 * @param density the density at the end point - avoids an extra call to Map::getVal()
	 * @param oldx the old x position
	 * @param oldy the old y position
	 * @param oldxwrap the old x wrap
	 * @param oldywrap the old y wrap
	 * @param startx the starting x position
	 * @param starty the starting y position
	 * @param startxwrap the starting x wrap
	 * @param startywrap the ending y wrap
	 * @return true if the end point passes the density and restricted checks
	 */
	bool checkEndPointDensity(const unsigned long &density, long &oldx, long &oldy, long &oldxwrap, long &oldywrap,
							   const long &startx, const long &starty, const long &startxwrap, const long &startywrap);
	
	
	/**
	 * @brief Check the end point for the given coordinates and density
	 * @param density the density at the end point - avoids an extra call to Map::getVal()
	 * @param oldx the old x position
	 * @param oldy the old y position
	 * @param oldxwrap the old x wrap
	 * @param oldywrap the old y wrap
	 * @param startx the starting x position
	 * @param starty the starting y position
	 * @param startxwrap the starting x wrap
	 * @param startywrap the ending y wrap
	 * @return true if the end point passes the density and restricted checks
	 */
	bool checkEndPointRestricted(const unsigned long &density, long &oldx, long &oldy, long &oldxwrap, long &oldywrap,
								  const long &startx, const long &starty, const long &startxwrap, const long &startywrap);
	
	
	/**
	 * @brief Performs the dispersal routine using the Step object to read starting positions and record the end positions.
	 * @param this_step the Step object for reading starting position and storing output distances and angles
	 */
	void disperse(Step &this_step);
	
};

#endif // DISPERSAL_H
