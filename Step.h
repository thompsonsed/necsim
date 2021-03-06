//This file is part of necsim project which is released under MIT license.
//See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.

/**
 * @author Sam Thompson
 * @date 09/08/2017
 * @file Step.h
 * @brief Contains the Step class for storing required data during a single step of a coalescence simulation.
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 */
#ifndef STEP_H
#define STEP_H

#include "Cell.h"
#include "MapLocation.h"
namespace necsim
{
    /**
     * @class Step
     * @brief Stores the elements associated with a single step in a coalescence simulation.
     *
     * This object should only contain transient variables that are used within a single simulation step and therefore
     * should not be important for pausing/resuming simulations.
     */
    struct Step : virtual public MapLocation
    {
        unsigned long chosen, coalchosen;
        bool coal, bContinueSim;
        unsigned int time_reference;
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
            coal = false;
            bContinueSim = true;
            time_reference = 0;
#ifdef verbose
            number_printed = 0;
#endif
        }

        Step(const Cell &cell)
        {
            x = cell.x;
            y = cell.y;
            xwrap = 0;
            ywrap = 0;
            coal = false;
            bContinueSim = true;
        }

        /**
         * @brief Removes all stored data from the step.
         * This should be run at the start of a single coalescence step.
         */
        void wipeData()
        {
            chosen = 0;
            coalchosen = 0;
            x = 0;
            y = 0;
            xwrap = 0;
            ywrap = 0;
            coal = false;
        }

    };
}
#endif