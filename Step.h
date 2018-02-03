//This file is part of NECSim project which is released under BSD-3 license.
//See file **LICENSE.txt** or visit https://opensource.org/licenses/BSD-3-Clause) for full license details.

/**
 * @class Step
 * @author Sam Thompson
 * @date 09/08/2017
 * @file Step.h
 * @brief Contains the Step class for storing required data during a single step of a coalescence simulation.
 * @copyright <a href="https://opensource.org/licenses/BSD-3-Clause">BSD-3 Licence.</a>
 */
#ifndef STEP_H
#define STEP_H
/**
 * @class Step
 * @brief Stores the elements associated with a single step in a coalescence simulation.
 * 
 * This object should only contain transient variables that are used within a single simulation step and therefore
 * should not be important for pausing/resuming simulations.
 */
struct Step
{
	unsigned long chosen, coalchosen;
	long oldx, oldy, oldxwrap, oldywrap;
	bool coal, bContinueSim;
	unsigned int time_reference;
	double distance;
	double angle;
#ifdef verbose
	long number_printed;
#endif

	/**
	 * @brief Step constructor
	 * @return 
	 */
	Step()
	{
		chosen = 0;
		coalchosen = 0;
		oldx = 0;
		oldy = 0;
		oldxwrap = 0;
		oldywrap = 0;
		coal = false;
		bContinueSim = true;
		time_reference = 0;
		distance = 0.0;
		angle = 0.0;
#ifdef verbose
		number_printed =0;
#endif
	}
	
	
	/**
	 * @brief Removes all stored data from the step.
	 * This should be run at the start of a single coalescence step.
	 */
	void wipeData()
	{
		chosen = 0;
		coalchosen = 0;
		oldx = 0;
		oldy = 0;
		oldxwrap = 0;
		oldywrap = 0;
		coal = false;
		distance = 0.0;
		angle = 0.0;
	}
	
	
};

#endif