// This file is part of necsim project which is released under MIT license.
// See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.

/**
 * @author Samuel Thompson
 * @date 24/03/17
 * @file SpatialTree.cpp
 *
 * @brief  Contains SpatialTree implementation as the main simulation object for spatially explicit
 * coalescence simulations.
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 */

#include <algorithm>
#include "SpatialTree.h"

#ifdef WIN_INSTALL
#include <windows.h>
#include <io.h>
#define dup2 _dup2
#endif

#include "eastl/heap.h"

namespace necsim
{
    void SpatialTree::runFileChecks()
    {
        // Now check that our folders exist
        checkFolders();
        // Now check for paused simulations
        checkSims();
    }

    void SpatialTree::checkFolders()
    {

        stringstream os;
        os << "Checking folder existance..." << flush;
        bool bFineMap, bCoarseMap, bFineMapHistorical, bCoarseMapHistorical, bSampleMask, bOutputFolder;
        try
        {
            bFineMap = doesExistNull(sim_parameters->fine_map_file);
        }
        catch(FatalException &fe)
        {
            writeError(fe.what());
            bFineMap = false;
        }
        try
        {
            bCoarseMap = doesExistNull(sim_parameters->coarse_map_file);
        }
        catch(FatalException &fe)
        {
            writeError(fe.what());
            bCoarseMap = false;
        }
        try
        {
            bFineMapHistorical = doesExistNull(sim_parameters->historical_fine_map_file);
        }
        catch(FatalException &fe)
        {
            writeError(fe.what());
            bFineMapHistorical = false;
        }
        try
        {
            bCoarseMapHistorical = doesExistNull(sim_parameters->historical_coarse_map_file);
        }
        catch(FatalException &fe)
        {
            writeError(fe.what());
            bCoarseMapHistorical = false;
        }
        bOutputFolder = checkOutputDirectory();
        try
        {
            bSampleMask = doesExistNull(sim_parameters->sample_mask_file);
        }
        catch(FatalException &fe)
        {
            writeError(fe.what());
            bSampleMask = false;
        }
        if(bFineMap && bCoarseMap && bFineMapHistorical && bCoarseMapHistorical && bOutputFolder && bSampleMask)
        {
            os << "\rChecking folder existance...done.                                                                "
               << endl;
            writeInfo(os.str());
            return;
        }
        else
        {
            throw FatalException("Required files do not all exist. Check program inputs.");
        }
    }

    void SpatialTree::setParameters()
    {
        if(!has_imported_vars)
        {
            Tree::setParameters();
            // Set the variables equal to the value from the Mapvars object.
            fine_map_input = sim_parameters->fine_map_file;
            coarse_map_input = sim_parameters->coarse_map_file;
            // historical map information
            historical_fine_map_input = sim_parameters->historical_fine_map_file;
            historical_coarse_map_input = sim_parameters->historical_coarse_map_file;
            desired_specnum = sim_parameters->desired_specnum;
            if(sim_parameters->landscape_type == "none")
            {
                sim_parameters->landscape_type = "closed";
            }
            if(sim_parameters->dispersal_method == "none")
            {
                sim_parameters->dispersal_method = "normal";
            }
        }
        else
        {
            throw FatalException("ERROR_MAIN_001: Variables already imported.");
        }
    }

    void SpatialTree::importMaps()
    {
        if(has_imported_vars)
        {
            // Set the dimensions
            landscape->setDims(sim_parameters);
            try
            {
                // Set the time variables
                landscape->checkMapExists();
                // landscape->setTimeVars(gen_since_historical,habitat_change_rate);
                // Import the fine map
                landscape->calcFineMap();
                // Import the coarse map
                landscape->calcCoarseMap();
                // Calculate the offset for the extremeties of each map
                landscape->calcOffset();
                // Import the historical maps;
                landscape->calcHistoricalFineMap();
                landscape->calcHistoricalCoarseMap();
                // Calculate the maximum values
                landscape->recalculateHabitatMax();
                importActivityMaps();
                samplegrid.importSampleMask(sim_parameters);
            }
            catch(FatalException &fe)
            {
                stringstream ss;
                ss << "Problem setting up map files: " << fe.what() << endl;
                throw FatalException(ss.str());
            }
        }
        else
        {
            throw FatalException("ERROR_MAIN_002: Variables not imported.");
        }
    }

    void SpatialTree::importActivityMaps()
    {
        death_map->import(sim_parameters->death_file,
                          sim_parameters->fine_map_x_size,
                          sim_parameters->fine_map_y_size,
                          NR);
        death_map->setOffsets(sim_parameters->fine_map_x_offset,
                              sim_parameters->fine_map_y_offset,
                              sim_parameters->grid_x_size,
                              sim_parameters->grid_y_size);
        if(sim_parameters->death_file == sim_parameters->reproduction_file)
        {
            reproduction_map = death_map;
        }
        else
        {

            reproduction_map->import(sim_parameters->reproduction_file,
                                     sim_parameters->fine_map_x_size,
                                     sim_parameters->fine_map_y_size,
                                     NR);
            reproduction_map->setOffsets(sim_parameters->fine_map_x_offset,
                                         sim_parameters->fine_map_y_offset,
                                         sim_parameters->grid_x_size,
                                         sim_parameters->grid_y_size);
        }
        // Now verify that the reproduction map is always non-zero when the density is non-zero.
        verifyActivityMaps();
    }

    unsigned long SpatialTree::getInitialCount()
    {
        unsigned long initcount = 0;
        // Get a count of the number of individuals on the grid.
        try
        {
            long max_x, max_y;
            if(samplegrid.isNull())
            {
                max_x = sim_parameters->fine_map_x_size;
                max_y = sim_parameters->fine_map_y_size;
            }
            else
            {
                if(sim_parameters->uses_spatial_sampling)
                {
                    max_x = samplegrid.sample_mask_exact.getCols();
                    max_y = samplegrid.sample_mask_exact.getRows();
                }
                else
                {
                    max_x = samplegrid.sample_mask.getCols();
                    max_y = samplegrid.sample_mask.getRows();
                }
            }
            long x, y, xwrap, ywrap;
            for(long i = 0; i < max_y; i++)
            {
                for(long j = 0; j < max_x; j++)
                {
                    x = j;
                    y = i;
                    xwrap = 0;
                    ywrap = 0;
                    samplegrid.recalculateCoordinates(x, y, xwrap, ywrap);
                    initcount += getIndividualsSampled(x, y, xwrap, ywrap, 0.0);
                }
            }
        }
        catch(exception &e)
        {
            throw FatalException(e.what());
        }
        // Set active and data at the correct sizes.
        if(initcount == 0)
        {
            throw runtime_error("Initial count is 0. No individuals to simulate. Exiting program.");
        }
        else
        {
            writeInfo("Initial count is " + to_string(initcount) + "\n");
        }
        if(initcount > 10000000000)
        {
            writeWarning("Initial count extremely large, RAM issues likely: " + to_string(initcount));
        }
        return initcount;
    }

    void SpatialTree::setupDispersalCoordinator()
    {
        dispersal_coordinator.setMaps(landscape, reproduction_map);
        dispersal_coordinator.setRandomNumber(NR);
        dispersal_coordinator.setGenerationPtr(&generation);
        dispersal_coordinator.setDispersal(sim_parameters->dispersal_method,
                                           sim_parameters->dispersal_file,
                                           sim_parameters->fine_map_x_size,
                                           sim_parameters->fine_map_y_size,
                                           sim_parameters->m_prob,
                                           sim_parameters->cutoff,
                                           sim_parameters->sigma,
                                           sim_parameters->tau,
                                           sim_parameters->restrict_self);
    }

    void SpatialTree::setup()
    {
        printSetup();
        if(has_paused)
        {
            if(!has_imported_pause)
            {
                setResumeParameters();
            }
            simResume();
            setupDispersalCoordinator();
        }
        else
        {
            setParameters();
            setInitialValues();
            importMaps();
            landscape->setLandscape(sim_parameters->landscape_type);
            setupDispersalCoordinator();
#ifdef DEBUG
            landscape->validateMaps();
#endif
            generateObjects();
        }
    }

    unsigned long SpatialTree::fillObjects(const unsigned long &initial_count)
    {
        active[0].setup(0, 0, 0, 0, 0, 0, 0);
        grid.setSize(sim_parameters->grid_y_size, sim_parameters->grid_x_size);
        unsigned long number_start = 0;
        stringstream os;
        os << "\rSetting up simulation...filling grid                           " << flush;
        writeInfo(os.str());
        // Add the individuals to the grid, and add wrapped individuals to their correct locations.
        // This loop adds individuals to data and active (for storing the coalescence tree and active lineage tracking)
        try
        {
            long x, y;
            long x_wrap, y_wrap;
            for(unsigned long i = 0; i < sim_parameters->sample_x_size; i++)
            {
                for(unsigned long j = 0; j < sim_parameters->sample_y_size; j++)
                {

                    x = i;
                    y = j;
                    x_wrap = 0;
                    y_wrap = 0;
                    samplegrid.recalculateCoordinates(x, y, x_wrap, y_wrap);
                    if(grid.get(y, x).getListSize() == 0)
                    {
                        unsigned long stored_next = grid.get(y, x).getNext();
                        unsigned long stored_nwrap = grid.get(y, x).getNwrap();
                        grid.get(y, x).initialise(landscape->getVal(x, y, 0, 0, 0));
                        grid.get(y, x).setNwrap(stored_nwrap);
                        grid.get(y, x).setNext(stored_next);
                    }
                    if(x_wrap == 0 && y_wrap == 0)
                    {
                        unsigned long sample_amount = getIndividualsSampled(x, y, 0, 0, 0.0);
                        if(sample_amount >= 1)
                        {
                            for(unsigned long k = 0; k < sample_amount; k++)
                            {
                                if(k >= grid.get(y, x).getMaxSize() && deme_sample <= 1.0)
                                {
                                    break;
                                }
                                if(number_start + 1 > initial_count)
                                {
                                    stringstream msg;
                                    msg << "Number start greater than initial count. Please report this error!" << endl;
                                    msg << "Number start: " << number_start << ". Initial count: " << initial_count
                                        << endl;
                                    throw out_of_range(msg.str());
                                }
                                else
                                {
                                    number_start++;
                                    unsigned long list_position_in = grid.get(y, x).addSpecies(number_start);
                                    // Add the species to active
                                    active[number_start].setup(x, y, 0, 0, number_start, list_position_in, 1);
                                    // Add a tip in the TreeNode for calculation of the coalescence tree at the
                                    // end of the simulation.
                                    // This also contains the start x and y position of the species.
                                    (*data)[number_start].setup(true, x, y, 0, 0);
                                    (*data)[number_start].setSpec(NR->d01());
                                    endactive++;
                                    enddata++;
                                }
                            }
                        }
                    }
                    else
                    {
                        unsigned long sample_amount = getIndividualsSampled(x, y, x_wrap, y_wrap, 0.0);
                        if(sample_amount >= 1)
                        {
                            for(unsigned long k = 0; k < sample_amount; k++)
                            {
                                if(number_start + 1 > initial_count)
                                {
                                    stringstream msg;
                                    msg << "Number start greater than initial count. Please report this error!";
                                    msg << "Number start: " << number_start << ". Initial count: " << initial_count
                                        << endl;
                                    throw out_of_range(msg.str());
                                }
                                else
                                {
                                    number_start++;
                                    // Add the lineage to the wrapped lineages
                                    active[number_start].setup((unsigned long) x,
                                                               (unsigned long) y,
                                                               x_wrap,
                                                               y_wrap,
                                                               number_start,
                                                               0,
                                                               1);
                                    addWrappedLineage(number_start, x, y);
                                    // Add a tip in the TreeNode for calculation of the coalescence tree at the
                                    // end of the simulation.
                                    // This also contains the start x and y position of the species.
                                    (*data)[number_start].setup(true, x, y, x_wrap, y_wrap);
                                    (*data)[number_start].setSpec(NR->d01());
                                    endactive++;
                                    enddata++;
                                }
                            }
                        }
                    }

                }
            }
            if(sim_parameters->uses_spatial_sampling)
            {

                samplegrid.convertBoolean(landscape, deme_sample, generation);
                // if there are no additional time points to sample at, we can remove the sample mask from memory.
                if(!(uses_temporal_sampling && this_step.time_reference < reference_times.size()))
                {
                    samplegrid.clearSpatialMask();
                }
            }
        }
        catch(out_of_range &out_of_range1)
        {
            stringstream ss;
            ss << "Fatal exception thrown when filling grid (out_of_range): " << out_of_range1.what() << endl;
            throw FatalException(ss.str());
        }
        catch(exception &fe)
        {
            throw FatalException("Fatal exception thrown when filling grid (other) \n");
        }

        if(number_start == initial_count)  // Check that the two counting methods match up.
        {
        }
        else
        {
            if(initial_count > 1.1 * number_start)
            {
                writeCritical("Data usage higher than neccessary - check allocation of individuals to the grid.");
                stringstream ss;
                ss << "Initial count: " << initial_count << "  Number counted: " << number_start << endl;
                writeWarning(ss.str());
            }
        }
#ifdef DEBUG
        validateLineages();
#endif
        return number_start;
    }

    unsigned long SpatialTree::getIndividualsSampled(const long &x, const long &y, const long &x_wrap,
                                                     const long &y_wrap, const double &current_gen)
    {
        //	if(sim_parameters->uses_spatial_sampling)
        //	{
        return static_cast<unsigned long>(max(floor(deme_sample * landscape->getVal(x, y, x_wrap, y_wrap, current_gen)
                                                    * samplegrid.getExactValue(x, y, x_wrap, y_wrap)), 0.0));
        //	}
        //	else
        //	{
        //		return static_cast<unsigned long>(max(floor(deme_sample * landscape->getVal(x, y, x_wrap, y_wrap, 0.0)), 0.0));
        //	}
    }

    unsigned long SpatialTree::getNumberLineagesAtLocation(const MapLocation &location)
    {
        if(location.isOnGrid())
        {
            return grid.get(location.y, location.x).getListSize();
        }
        unsigned long next = grid.get(location.y, location.x).getNext();
        unsigned long total = 0;
        while(next != 0)
        {
            if(static_cast<MapLocation>(active[next]) == location)
            {
                total++;
            }
            next = active[next].getNext();
        }
        return total;
    }

    unsigned long SpatialTree::getNumberIndividualsAtLocation(const MapLocation &location)
    {
        return landscape->getVal(location.x, location.y, location.xwrap, location.ywrap, generation);
    }

    void SpatialTree::removeOldPosition(const unsigned long &chosen)
    {
        long nwrap = active[chosen].getNwrap();
        long oldx = active[chosen].getXpos();
        long oldy = active[chosen].getYpos();
        if(nwrap == 0)
        {
#ifdef DEBUG

            if(active[chosen].getXwrap() != 0 || active[chosen].getYwrap() != 0)
            {
                active[chosen].logActive(50);
                throw FatalException("ERROR_MOVE_015: Nwrap not set correctly. Nwrap 0, but x and y wrap not 0. ");
            }
#endif // DEBUG
            // Then the lineage exists in the main lineage_indices;
            // debug (can be removed later)
#ifdef historical_mode
            if(grid.get(y, x).getMaxsize() < active[chosen].getListpos())
            {
                stringstream ss;
                ss << "grid maxsize: " << grid.get(y, x).getMaxsize() << endl;
                writeCritical(ss.str());
                throw FatalException("ERROR_MOVE_001: Listpos outside maxsize. Check move programming function.");
            }
#endif
            // delete the species from the lineage_indices
            grid.get(oldy, oldx).deleteSpecies(active[chosen].getListpos());
            // clear out the variables.
            active[chosen].setNext(0);
            active[chosen].setNwrap(0);
            active[chosen].setListPosition(0);
        }
        else  // need to loop over the nwrap to check nexts
        {
            if(nwrap == 1)
            {
                grid.get(oldy, oldx).setNext(active[chosen].getNext());
                // Now reduce the nwrap of the lineages that have been effected.
                long nextpos = active[chosen].getNext();
                // loop over the rest of the lineage_indices, reducing the nwrap
                while(nextpos != 0)
                {
                    active[nextpos].decreaseNwrap();
                    nextpos = active[nextpos].getNext();
                }
                // decrease the nwrap
                grid.get(oldy, oldx).decreaseNwrap();
                active[chosen].setNwrap(0);
                active[chosen].setNext(0);
                active[chosen].setListPosition(0);
            }
            else
            {
                long lastpos = grid.get(oldy, oldx).getNext();
                while(active[lastpos].getNext()
                      != chosen)  // loop until we reach the next, then set the next correctly.
                {
                    lastpos = active[lastpos].getNext();
                }
                if(lastpos != 0)
                {
                    active[lastpos].setNext(active[chosen].getNext());
#ifdef DEBUG
                    if(active[lastpos].getNwrap() != (active[chosen].getNwrap() - 1))
                    {
                        writeLog(50, "Logging last position: ");
                        active[lastpos].logActive(50);
                        writeLog(50, "Logging chosen position: ");
                        active[chosen].logActive(50);
                        throw FatalException("ERROR_MOVE_022: nwrap setting of either chosen or the "
                                             "lineage wrapped before chosen. Check move function.");
                    }
#endif // DEBUG
                    lastpos = active[lastpos].getNext();
                    while(lastpos != 0)
                    {
                        active[lastpos].decreaseNwrap();
                        lastpos = active[lastpos].getNext();
                    }
                }
                else
                {
#ifdef DEBUG
                    writeLog(50, "Logging chosen");
                    active[chosen].logActive(50);
#endif // DEBUG
                    throw FatalException("ERROR_MOVE_024: Last position before chosen is 0 - this is impossible.");
                }
                grid.get(oldy, oldx).decreaseNwrap();
                active[chosen].setNwrap(0);
                active[chosen].setNext(0);
                active[chosen].setListPosition(0);
            }
#ifdef DEBUG
            unsigned long iCount = 1;
            long pos = grid.get(oldy, oldx).getNext();
            if(pos == 0)
            {
                iCount = 0;
            }
            else
            {
                unsigned long c = 0;
                while(active[pos].getNext() != 0)
                {
                    c++;
                    iCount++;
                    pos = active[pos].getNext();
                    if(c > std::numeric_limits<unsigned long>::max())
                    {
                        throw FatalException("ERROR_MOVE_014: Wrapping exceeds numeric limits.");
                    }
                }
            }

            if(iCount != grid.get(oldy, oldx).getNwrap())
            {
                stringstream ss;
                ss << "Nwrap: " << grid.get(oldy, oldx).getNwrap() << " Counted lineages: " << iCount << endl;
                writeLog(50, ss);
                throw FatalException("Nwrap not set correctly after move for grid cell");
            }
#endif // DEBUG
        }
    }

    void SpatialTree::calcMove()
    {
        dispersal_coordinator.disperse(this_step);
    }

    long double SpatialTree::calcMinMax(const unsigned long &current)
    {
        // this formula calculates the speciation rate required for speciation to have occured on this branch.
        // need to allow for the case that the number of gens was 0
        long double newminmax = 1;
        long double oldminmax = active[current].getMinmax();
        if((*data)[active[current].getReference()].getGenRate() == 0)
        {
            newminmax = (*data)[active[current].getReference()].getSpecRate();
        }
        else
        {
            // variables need to be defined separately for the decimal division to function properly.
            long double tmpdSpec = (*data)[active[current].getReference()].getSpecRate();
            long double tmpiGen = (*data)[active[current].getReference()].getGenRate();
            newminmax = 1 - (pow(1 - tmpdSpec, (1 / tmpiGen)));
        }
        long double toret = min(newminmax, oldminmax);
        return toret;
    }

    void SpatialTree::calcNewPos()
    {

        // Calculate the new position of the move, whilst also calculating the probability of coalescence.
        unsigned long nwrap = active[this_step.chosen].getNwrap();
        if(this_step.isOnGrid())
        {
            // Debug check (to remove later)
            if(nwrap != 0)
            {
                throw FatalException("Nwrap not set correctly. Check move programming function.");
            }
            // then the procedure is relatively simple.
            // check for coalescence
            // check if the grid needs to be updated.
            if(grid.get(this_step.y, this_step.x).getMaxSize()
               != landscape->getVal(this_step.x, this_step.y, 0, 0, generation))
            {
                grid.get(this_step.y, this_step.x)
                    .setMaxsize(landscape->getVal(this_step.x, this_step.y, 0, 0, generation));
            }
            this_step.coalchosen = grid.get(this_step.y, this_step.x).getRandLineage(NR);
#ifdef DEBUG
            if(this_step.coalchosen != 0)
            {
                if(active[this_step.coalchosen].getXpos() != (unsigned long) this_step.x ||
                   active[this_step.coalchosen].getYpos() != (unsigned long) this_step.y ||
                   active[this_step.coalchosen].getXwrap() != this_step.xwrap|| active[this_step.coalchosen].getYwrap() != this_step.ywrap)
                {
                    writeLog(50, "Logging this_step.chosen:");
                    active[this_step.chosen].logActive(50);
                    writeLog(50, "Logging this_step.coalchosen: ");
                    active[this_step.coalchosen].logActive(50);
                    throw FatalException("Nwrap not set correctly. Please report this bug.");
                }
            }
#endif
            if(this_step.coalchosen == 0)  // then the lineage can be placed in the empty space.
            {
                long tmplistindex = grid.get(this_step.y, this_step.x).addSpecies(this_step.chosen);
                // check
                if(grid.get(this_step.y, this_step.x).getLineageIndex(tmplistindex) != this_step.chosen)
                {
                    throw FatalException("Grid index not set correctly for species. Check move programming function.");
                }
#ifdef historical_mode
                if(grid.get(y, x).getListsize() > grid.get(y, x).getMaxsize())
                {
                    throw FatalException(
                        "ERROR_MOVE_001: Listpos outside maxsize. Check move programming function.");
                }
#endif
                active[this_step.chosen].setNwrap(0);
                active[this_step.chosen].setListPosition(tmplistindex);
                this_step.coal = false;
            }
            else  // then coalescence has occured
            {
                active[this_step.chosen].setNwrap(0);
                active[this_step.chosen].setListPosition(0);
                // DO THE COALESCENCE STUFF
                this_step.coal = true;
            }
        }
        else  // need to check all the possible places the lineage could be.
        {
            if(nwrap != 0)
            {
                throw FatalException("Nwrap not set correctly in move.");
            }
            nwrap = grid.get(this_step.y, this_step.x).getNwrap();
            if(nwrap != 0)  // then coalescence is possible and we need to loop over the nexts to check those that are
                // in the same position
            {
                calcWrappedCoalescence(nwrap);
            }
            else  // just add the lineage to next.
            {
                if(grid.get(this_step.y, this_step.x).getNext() != 0)
                {
                    throw FatalException("No nwrap recorded, but next is non-zero.");
                }
                this_step.coalchosen = 0;
                this_step.coal = false;
                grid.get(this_step.y, this_step.x).setNext(this_step.chosen);
                active[this_step.chosen].setNwrap(1);
                active[this_step.chosen].setNext(0);
                grid.get(this_step.y, this_step.x).increaseNwrap();
#ifdef DEBUG
                if(grid.get(this_step.y, this_step.x).getNwrap() != 1)
                {
                    throw FatalException("Nwrap not set correctly in move.");
                }
#endif
            }
            if(this_step.coalchosen != 0)
            {
                if(active[this_step.coalchosen].getXpos() != (unsigned long) this_step.x
                   || active[this_step.coalchosen].getYpos() != (unsigned long) this_step.y
                   || active[this_step.coalchosen].getXwrap() != this_step.xwrap
                   || active[this_step.coalchosen].getYwrap() != this_step.ywrap)
                {
#ifdef DEBUG
                    writeLog(50, "Logging this_step.chosen:");
                    active[this_step.chosen].logActive(50);
                    writeLog(50, "Logging this_step.coalchosen: ");
                    active[this_step.coalchosen].logActive(50);
#endif // DEBUG
                    throw FatalException("Nwrap not set correctly. Check move "
                                         "programming function.");
                }
            }
        }
    }

    void SpatialTree::calcWrappedCoalescence(const unsigned long &nwrap)
    {

        // Count the possible matches of the position.
        unsigned long matches = 0;
        // Create an array containing the lineage_indices of active references for those that match as
        // this stops us having to loop twice over the same lineage_indices.
        vector<unsigned long> match_list(nwrap);
        unsigned long next_active;
        next_active = grid.get(this_step.y, this_step.x).getNext();
        // Count if the first "next" matches
        if(active[next_active].getXwrap() == this_step.xwrap && active[next_active].getYwrap() == this_step.ywrap)
        {
#ifdef DEBUG
            if(active[next_active].getNwrap() != 1)
                {
                    throw FatalException("ERROR_MOVE_022a: Nwrap not set correctly in move.");
                }
#endif
            match_list[matches] = next_active;  // add the match to the lineage_indices of matches.
            matches++;
        }
        // Now loop over the remaining nexts counting matches
        //#ifdef DEBUG
        unsigned long ncount = 1;
        //#endif
        while(active[next_active].getNext() != 0)
        {
            next_active = active[next_active].getNext();
            if(active[next_active].getXwrap() == this_step.xwrap && active[next_active].getYwrap() == this_step.ywrap)
            {
                match_list[matches] = next_active;
                matches++;
            }
            // check
            //#ifdef DEBUG
            ncount++;
#ifdef DEBUG
            if(active[next_active].getNwrap() != ncount)
                {
                    throw FatalException("ERROR_MOVE_022d: Nwrap not set correctly in move.");
                }
#endif
        }
        if(nwrap != ncount)
        {
            throw FatalException("Nwrap not set correctly in move.");
        }
        // Matches now contains the number of lineages at the exact x,y, xwrap and ywrap position.
        // Check if there were no matches at all
        if(matches == 0)
        {
            this_step.coalchosen = 0;
            this_step.coal = false;
            active[next_active].setNext(this_step.chosen);
            grid.get(this_step.y, this_step.x).increaseNwrap();
            active[this_step.chosen].setNwrap(grid.get(this_step.y, this_step.x).getNwrap()

                                             );
            active[this_step.chosen].setListPosition(0);
        }
        else  // if there were matches, generate a random number to see if coalescence occured or not
        {
            unsigned long randwrap = floor(NR->d01() * (landscape
                    ->getVal(this_step.x, this_step.y, this_step.xwrap, this_step.ywrap, generation)) + 1);
            // Get the random reference from the match lineage_indices.
            // If the movement is to an empty space, then we can update the chain to include the new
            // lineage.
            if(randwrap > matches)  // coalescence has not occured
            {
                // os << "This shouldn't happen" << endl;
                this_step.coalchosen = 0;
                this_step.coal = false;
                active[next_active].setNext(this_step.chosen);
                grid.get(this_step.y, this_step.x).increaseNwrap();

                active[this_step.chosen].setNwrap(grid.get(this_step.y, this_step.x).getNwrap());
                active[this_step.chosen].setListPosition(0);
            }
            else  // coalescence has occured
            {
                this_step.coal = true;
                this_step.coalchosen = match_list[randwrap - 1];
                active[this_step.chosen].setEndpoint<Step>(this_step);
                if(this_step.coalchosen == 0)
                {
                    throw FatalException("Coalescence attempted with lineage of 0.");
                }
            }
        }
#ifdef historical_mode
        if(grid.get(this_step.y, this_step.x).getMaxsize() < active[this_step.chosen].getListpos())
            {
                throw FatalException("Listpos outside maxsize. Check move programming function.");
            }
#endif
    }

    void SpatialTree::switchPositions(const unsigned long &chosen)
    {
#ifdef DEBUG
        if(chosen > endactive)
        {
            stringstream ss;
            ss << "chosen: " << chosen << " endactive: " << endactive << endl;
            writeLog(50, ss);
            throw FatalException("ERROR_MOVE_023: Chosen is greater than endactive. Check move function.");
        }
#endif // DEBUG
        if(chosen != endactive)
        {
            // This routine assumes that the previous chosen position has already been deleted.
            DataPoint tmpdatactive;
            tmpdatactive.setup(active[chosen]);
            // now need to remove the chosen lineage from memory, by replacing it with the lineage that lies in the last
            // place.
            if(active[endactive].getXwrap() == 0
               && active[endactive].getYwrap() == 0)  // if the end lineage is simple, we can just copy it across.
            {
                // check endactive
                if(active[endactive].getNwrap() != 0)
                {
                    stringstream ss;
                    ss << "Nwrap is not set correctly for endactive (nwrap should be 0, but is ";
                    ss << active[endactive].getNwrap() << " ). Identified during switch of positions." << endl;
                    writeError(ss.str());
                }
                grid.get(active[endactive].getYpos(), active[endactive].getXpos())
                    .setSpecies(active[endactive].getListpos(), chosen);
                active[chosen].setup(active[endactive]);
                active[endactive].setup(tmpdatactive);
                active[endactive].setNwrap(0);
                active[endactive].setNext(0);
            }
            else  // else the end lineage is wrapped, and needs to be processed including the wrapping routines.
            {
                if(active[endactive].getNwrap() == 0)
                {
                    stringstream ss;
                    ss << "Nwrap is not set correctly for endactive (nwrap incorrectly 0).";
                    ss << "Identified during switch of positions." << endl;
                    writeError(ss.str());
                }
                //				os << "wrap"<<endl;
                long tmpactive = grid.get(active[endactive].getYpos(), active[endactive].getXpos()).getNext();
                unsigned long tmpnwrap = active[endactive].getNwrap();

                // if the wrapping is just once, we need to set the grid next to the chosen variable.
                if(tmpnwrap == 1)
                {
                    // check
                    if(grid.get(active[endactive].getYpos(), active[endactive].getXpos()).getNext() != endactive)
                    {
                        throw FatalException(string("Nwrap for endactive not set correctly. Nwrap is 1, but "
                                                    "lineage at 1st position is " + to_string((long long) grid
                                .get(active[endactive].getYpos(), active[endactive].getXpos()).getNext())
                                                    + ". Identified during the move."));
                    }
                    grid.get(active[endactive].getYpos(), active[endactive].getXpos()).setNext(chosen);
                }
                else  // otherwise, we just set the next to chosen instead of endactive.
                {
                    unsigned long tmpcount = 0;
                    // loop over nexts until we reach the right lineage.
                    while(active[tmpactive].getNext() != endactive)
                    {
                        tmpactive = active[tmpactive].getNext();
                        tmpcount++;
#ifdef DEBUG
                        if(tmpcount > tmpnwrap)
                        {
                            writeLog(30, "ERROR_MOVE_013: NON FATAL. Looping has not encountered a match, "
                                         "despite going further than required. Check nwrap counting.");
                            if(tmpactive == 0)
                            {
                                stringstream ss;
                                ss << "gridnext: "
                                   << grid.get(active[endactive].getYpos(), active[endactive].getXpos()).getNext()
                                   << endl;
                                ss << "endactive: " << endactive << endl;
                                ss << "tmpactive: " << tmpactive << endl;
                                ss << "tmpnwrap: " << tmpnwrap << " tmpcount: " << tmpcount
                                   << endl;
                                writeLog(50, ss);
                                writeLog(50, "Logging chosen:");
                                active[chosen].logActive(50);
                                throw FatalException("No match found, please report this bug.");
                            }
                        }
#endif // DEBUG
                    }
                    active[tmpactive].setNext(chosen);
                }
                active[chosen].setup(active[endactive]);
                active[endactive].setup(tmpdatactive);

                // check - debugging
                unsigned long testwrap = active[chosen].getNwrap();
                unsigned long testnext = grid.get(active[chosen].getYpos(), active[chosen].getXpos()).getNext();
                for(unsigned long i = 1; i < testwrap; i++)
                {
                    testnext = active[testnext].getNext();
                }

                if(testnext != chosen)
                {
                    throw FatalException("ERROR_MOVE_009: Nwrap position not set correctly after coalescence. "
                                         "Check move process.");
                }
            }
        }
        endactive--;
    }

    void SpatialTree::calcNextStep()
    {
        calcMove();
        // Calculate the new position, perform the move if coalescence doesn't occur or
        // return the variables for the coalescence event if coalescence does occur.
        active[this_step.chosen].setEndpoint<Step>(this_step);
        calcNewPos();
    }

    unsigned long SpatialTree::estSpecnum()
    {
        // This bit has been removed as it has a very significant performance hit and is not required for most simulations.
        // As of version 3.2 it was fully compatible with the rest of the simulation, however. See estSpecnum for commented
        // code
        // (removed from here to make things tidier).
        // This bit was moved from runSimulation() to make things tidier there.
        /*
        if(steps%1000000==0)
    {
                time(&now);
                if(now - time_taken>200&&dPercentComplete>95)
                {
                                time(&time_taken);
                                unsigned long specnum = est_specnum();
                                os << "Estimated number of species: " << specnum <<
                                flush;
                                if(specnum<desired_specnum)
                                {
                                                os << " - desired
                                                number of species reached." << endl << "Halting
                                                simulations..." << endl;
                                                bContinueSim = false;
                                }
                                else
                                {
                                                os << endl;
                                }
                }
    }
    //*/
        long double dMinmax = 0;
        // first loop to find the maximum speciation rate required
        for(unsigned int i = 1; i <= endactive; i++)
        {
            long double tmpminmax = calcMinMax(i);
            active[i].setMinmax(tmpminmax);
            dMinmax = (long double) max(dMinmax, tmpminmax);
        }
        for(unsigned long i = 0; i <= enddata; i++)
        {
            if((*data)[i].isTip())
            {
                (*data)[i].setExistence(true);
            }
            double maxret = 1;
            if((*data)[i].getGenRate() == 0)
            {
                maxret = 1;
            }
            else
            {
                maxret = (*data)[i].getGenRate();
            }
            // This is the line that compares the individual random numbers against the speciation rate.
            if((*data)[i].getSpecRate() < (1 - pow(double(1 - dMinmax), maxret)))
            {
                (*data)[i].speciate();
            }
        }
        bool loop = true;
        while(loop)
        {
            loop = false;
            for(unsigned int i = 0; i <= enddata; i++)
            {
                if((*data)[i].exists() && !(*data)[(*data)[i].getParent()].exists() && !(*data)[i].hasSpeciated())
                {
                    loop = true;
                    (*data)[(*data)[i].getParent()].setExistence(true);
                }
            }
        }
        unsigned long iSpecies = 0;
        for(unsigned int i = 0; i <= enddata; i++)
        {
            if((*data)[i].exists() && (*data)[i].hasSpeciated())
            {
                iSpecies++;
            }
        }
        for(unsigned int i = 0; i <= enddata; i++)
        {
            (*data)[i].qReset();
        }
        //		os << "Estimated species number is: " << iSpecies << endl;
        return iSpecies;
    }

#ifdef historical_mode
    void SpatialTree::historicalStepChecks()
    {
        if(landscape->getVal(this_step.x, this_step.y, this_step.xwrap, this_step.ywrap, generation) == 0)
        {
            throw FatalException(
                string("ERROR_MOVE_008: Dispersal attempted from non-forest. Check dispersal function. Forest "
                       "cover: " +
                       to_string((long long)landscape->getVal(this_step.x, this_step.y, this_step.xwrap,
                                                             this_step.ywrap, generation))));
        }
    }
#endif

    void SpatialTree::incrementGeneration()
    {
        Tree::incrementGeneration();
        if(landscape->updateMap(generation))
        {
            dispersal_coordinator.updateDispersalMap();
        }

        checkTimeUpdate();
        // check if the map is historical yet
        landscape->checkHistorical(generation);

    }

#ifdef DEBUG

    void SpatialTree::debugDispersal()
    {
        if(landscape->getVal(this_step.x, this_step.y, this_step.xwrap, this_step.ywrap, generation) == 0)
        {
            throw FatalException(
                    string("ERROR_MOVE_007: Dispersal attempted to non-forest. "
                           "Check dispersal function. Forest cover: " +
                           to_string((long long) landscape->getVal(this_step.x, this_step.y, this_step.xwrap,
                                                                   this_step.ywrap, generation))));
        }
    }

#endif

    void SpatialTree::updateStepCoalescenceVariables()
    {
        while(!death_map->actionOccurs(active[this_step.chosen].getXpos(),
                                       active[this_step.chosen].getYpos(),
                                       active[this_step.chosen].getXwrap(),
                                       active[this_step.chosen].getYwrap()))
        {
            this_step.chosen = NR->i0(endactive - 1) + 1;  // cannot be 0
        }
        recordLineagePosition();
#ifdef historical_mode
        historicalStepChecks();
#endif
    }

    void SpatialTree::recordLineagePosition()
    {
        Tree::updateStepCoalescenceVariables();
        // record old position of lineage
        this_step.x = active[this_step.chosen].getXpos();
        this_step.y = active[this_step.chosen].getYpos();
        this_step.xwrap = active[this_step.chosen].getXwrap();
        this_step.ywrap = active[this_step.chosen].getYwrap();
    }

    void SpatialTree::addLineages(double generation_in)
    {
        // Store all tree nodes to add in a vector
        vector<TreeNode> data_added{};
        // Store all added active lineages in a vector
        vector<DataPoint> active_added{};
        // Update the sample grid boolean mask, if required.
        if(sim_parameters->uses_spatial_sampling)
        {
            samplegrid.convertBoolean(landscape, deme_sample, generation_in);
        }
        for(unsigned long i = 0; i < sim_parameters->sample_x_size; i++)
        {
            for(unsigned long j = 0; j < sim_parameters->sample_y_size; j++)
            {
                long x, y;
                x = i;
                y = j;
                long xwrap, ywrap;
                xwrap = 0;
                ywrap = 0;
                samplegrid.recalculateCoordinates(x, y, xwrap, ywrap);
                if(samplegrid.getVal(x, y, xwrap, ywrap))
                {
                    unsigned long num_to_add = countCellExpansion(x, y, xwrap, ywrap, generation_in, data_added);
                    expandCell(x, y, xwrap, ywrap, generation_in, num_to_add, data_added, active_added);
                }
            }
        }
        // now resize data and active if necessary
        checkSimSize(data_added.size(), active_added.size());
        // Add all the data_added and active_added lineages
        for(auto &item : data_added)
        {
            enddata++;
            (*data)[enddata] = item;
        }
        for(auto &item : active_added)
        {
            endactive++;
            active[endactive] = item;
            if(item.getXwrap() != 0 || item.getYwrap() != 0)
            {
                addWrappedLineage(endactive, item.getXpos(), item.getYpos());
            }
        }
        // double check sizes
        if(enddata >= data->size() || endactive >= active.size())
        {
            throw FatalException("ERROR_MAIN_012: FATAL. Enddata or endactive is greater than the size of the "
                                 "relevant object. Programming error likely.");
        }
        if(endactive > startendactive)
        {
            startendactive = endactive;
        }
#ifdef DEBUG
        validateLineages();
#endif
    }

    string SpatialTree::simulationParametersSqlInsertion()
    {
        string to_execute;
        to_execute = "INSERT INTO SIMULATION_PARAMETERS VALUES(" + to_string((long long) seed) + ","
                     + to_string((long long) job_type);
        to_execute += ",'" + out_directory + "'," + boost::lexical_cast<std::string>((long double) spec) + ","
                      + to_string((long double) sim_parameters->sigma) + ",";
        to_execute += to_string((long double) sim_parameters->tau) + "," + to_string((long double) sim_parameters->deme)
                      + ",";
        to_execute += to_string((long double) sim_parameters->deme_sample) + "," + to_string((long long) maxtime) + ",";
        to_execute += to_string((long double) sim_parameters->dispersal_relative_cost) + ","
                      + to_string((long long) desired_specnum) + ",";
        to_execute += to_string((long double) sim_parameters->habitat_change_rate) + ",";
        to_execute += to_string((long double) sim_parameters->gen_since_historical) + ",'" + sim_parameters->times_file
                      + "','";
        to_execute += coarse_map_input + "'," + to_string((long long) sim_parameters->coarse_map_x_size) + ",";
        to_execute += to_string((long long) sim_parameters->coarse_map_y_size) + ","
                      + to_string((long long) sim_parameters->coarse_map_x_offset) + ",";
        to_execute += to_string((long long) sim_parameters->coarse_map_y_offset) + ","
                      + to_string((long long) sim_parameters->coarse_map_scale) + ",'";
        to_execute += fine_map_input + "'," + to_string((long long) sim_parameters->fine_map_x_size) + ","
                      + to_string((long long) sim_parameters->fine_map_y_size);
        to_execute += "," + to_string((long long) sim_parameters->fine_map_x_offset) + ","
                      + to_string((long long) sim_parameters->fine_map_y_offset) + ",'";
        to_execute += sim_parameters->sample_mask_file + "'," + to_string((long long) sim_parameters->grid_x_size) + ","
                      + to_string((long long) sim_parameters->grid_y_size) + ","
                      + to_string((long long) sim_parameters->sample_x_size) + ", ";
        to_execute += to_string((long long) sim_parameters->sample_y_size) + ", ";
        to_execute += to_string((long long) sim_parameters->sample_x_offset) + ", ";
        to_execute += to_string((long long) sim_parameters->sample_y_offset) + ", '";
        to_execute += historical_coarse_map_input + "','" + historical_fine_map_input + "'," + to_string(sim_complete);
        to_execute += ", '" + sim_parameters->dispersal_method + "', ";
        to_execute += boost::lexical_cast<std::string>(sim_parameters->m_prob) + ", ";
        to_execute += to_string((long double) sim_parameters->cutoff) + ", ";
        to_execute += to_string(sim_parameters->restrict_self) + ", '";
        to_execute += sim_parameters->landscape_type + "', ";
        // Now save the protracted speciation variables (not relevant in this simulation scenario)
        to_execute += protractedVarsToString();
        to_execute += ", '" + sim_parameters->dispersal_file + "'";
        to_execute += ");";
        return to_execute;
    }

    void SpatialTree::simPause()
    {
        // This function dumps all simulation data to a file.
        auto out1 = initiatePause();
        dumpMain(out1);
        dumpMap(out1);
        dumpActive(out1);
        dumpGrid(out1);
        dumpData(out1);
        completePause(out1);
    }

    void SpatialTree::dumpMap(shared_ptr<ofstream> out)
    {
        try
        {
            // Output the data object
            *out << *landscape;
        }
        catch(exception &e)
        {
            stringstream ss;
            ss << "Failed to perform dump of map: " << e.what() << endl;
            writeCritical(ss.str());
        }
    }

    void SpatialTree::dumpGrid(shared_ptr<ofstream> out)
    {
        try
        {
            // Output the data object
            *out << grid;
        }
        catch(exception &e)
        {
            stringstream ss;
            ss << "Failed to perform dump of grid: " << e.what() << endl;
            writeCritical(ss.str());
        }
    }

    void SpatialTree::simResume()
    {
        initiateResume();
        auto is = openSaveFile();
        // now load the objects
        loadMainSave(is);
        loadMapSave(is);
        setObjectSizes();
        loadActiveSave(is);
        loadGridSave(is);
        loadDataSave(is);
        time(&sim_start);
        writeInfo("\rLoading data from temp file...done.\n");
        sim_parameters->printVars();
    }

    void SpatialTree::loadGridSave(shared_ptr<ifstream> in1)
    {
        grid.setSize(sim_parameters->grid_y_size, sim_parameters->grid_x_size);
        *in1 >> grid;
        try
        {
            stringstream os;
            os << "\rLoading data from temp file...grid..." << flush;
            // New method for re-creating grid data from active lineages
            // First initialise the empty grid object
            writeInfo(os.str());
            for(unsigned long i = 0; i < sim_parameters->grid_y_size; i++)
            {
                for(unsigned long j = 0; j < sim_parameters->grid_x_size; j++)
                {
                    grid.get(i, j).initialise(landscape->getVal(j, i, 0, 0, generation));
                }
            }
            // Now fill the grid object with lineages from active. Only need to loop once.
            for(unsigned long i = 1; i <= endactive; i++)
            {
                if(active[i].getXwrap() == 0 && active[i].getYwrap() == 0)
                {
                    grid.get(active[i].getYpos(), active[i].getXpos()).setSpeciesEmpty(active[i].getListpos(), i);
                }
                else
                {
                    if(active[i].getNwrap() == 0)
                    {
                        throw runtime_error("Nwrap should not be 0 if x and y wrap are not 0. Programming error likely.");
                    }
                    if(active[i].getNwrap() == 1)
                    {
                        grid.get(active[i].getYpos(), active[i].getXpos()).setNext(i);
                    }
                    grid.get(active[i].getYpos(), active[i].getXpos()).increaseNwrap();
                }
            }
        }
        catch(exception &e)
        {
            string msg;
            msg = "Failure to import grid from temp grid: " + string(e.what());
            throw FatalException(msg);
        }
    }

    void SpatialTree::loadMapSave(shared_ptr<ifstream> in1)
    {
        // Input the map object
        try
        {
            stringstream os;
            os << "\rLoading data from temp file...map..." << flush;
            writeInfo(os.str());
            landscape->setDims(sim_parameters);
            *in1 >> *landscape;
            samplegrid.importSampleMask(sim_parameters);
            importActivityMaps();
        }
        catch(exception &e)
        {
            string msg;
            msg = "Failure to import data from temp map: " + string(e.what());
            throw FatalException(msg);
        }
    }

    void SpatialTree::verifyActivityMaps()
    {
        bool has_printed = false;
        if(!(sim_parameters->death_file == "none" || sim_parameters->death_file == "null") && !death_map->isNull())
        {
            for(unsigned long i = 0; i < sim_parameters->fine_map_y_size; i++)
            {
                for(unsigned long j = 0; j < sim_parameters->fine_map_x_size; j++)
                {
                    if(death_map->get(i, j) == 0.0 && landscape->getValFine(j, i, 0.0) != 0)
                    {
                        stringstream ss;
                        ss << "Location: " << j << ", " << i << endl;
                        ss << "Death value: " << death_map->get(i, j) << endl;
                        ss << "Density: " << landscape->getValFine(j, i, 0.0) << endl;
                        writeInfo(ss.str());
                        throw FatalException("Death map is zero where density is non-zero. "
                                             "This will cause an infinite loop.");
                    }

#ifdef DEBUG
                    if(!has_printed)
                    {
                        if(landscape->getValFine(j, i, 0.0) == 0 && death_map->get(i, j) != 0.0)
                        {
                            stringstream ss;
                            ss << "Density is zero where death map is non-zero for " << j << ", " << i << endl;
                            ss << "Density: " << landscape->getValFine(j, i, 0.0) << endl;
                            ss << "Death map: " << death_map->get(i, j) << endl;
                            ss << "This is likely incorrect." << endl;
                            writeCritical(ss.str());
                        }
                    }
#else // NDEBUG
                    if(!has_printed)
                    {
                        if(landscape->getValFine(j, i, 0.0) == 0 && death_map->get(i, j) != 0.0)
                        {
                            has_printed = true;
                            writeCritical("Density is zero where death map is non-zero. This is likely incorrect.");
                        }
                    }
#endif // DEBUG
                }
            }
#ifdef DEBUG
            writeLog(10, "\nActivity map validation complete.");
#endif // DEBUG
        }
        if(!(sim_parameters->reproduction_file == "none" || sim_parameters->reproduction_file == "null")
           && !reproduction_map->isNull())
        {
            has_printed = false;
            for(unsigned long i = 0; i < sim_parameters->fine_map_y_size; i++)
            {
                for(unsigned long j = 0; j < sim_parameters->fine_map_x_size; j++)
                {
                    if(reproduction_map->get(i, j) == 0.0 && landscape->getValFine(j, i, 0.0) != 0)
                    {
                        stringstream ss;
                        ss << "Location: " << j << ", " << i << endl;
                        ss << "Reproduction value: " << reproduction_map->get(i, j) << endl;
                        ss << "Density: " << landscape->getValFine(j, i, 0.0) << endl;
                        writeInfo(ss.str());
                        throw FatalException("Reproduction map is zero where density is non-zero. "
                                             "This will cause an infinite loop.");
                    }
#ifdef DEBUG
                    if(landscape->getValFine(j, i, 0.0) == 0 && reproduction_map->get(i, j) != 0.0)
                    {
                        stringstream ss;
                        ss << "Density is zero where reproduction map is non-zero for " << j << ", " << i << endl;
                        ss << "Density: " << landscape->getValFine(j, i, 0.0) << endl;
                        ss << "Reproduction map: " << reproduction_map->get(i, j) << endl;
                        ss << "This is likely incorrect." << endl;
                        writeCritical(ss.str());
                    }
#else // NDEBUG
                    if(!has_printed)
                    {
                        if(landscape->getValFine(j, i, 0.0) == 0 && reproduction_map->get(i, j) != 0.0)
                        {
                            has_printed = true;
                            writeCritical(
                                    "Density is zero where reproduction map is non-zero. This is likely incorrect.");
                        }
                    }
#endif // NDEBUG
                }
            }
        }

    }

    void SpatialTree::addWrappedLineage(unsigned long numstart, long x, long y)
    {

        if(grid.get(y, x).getNwrap() == 0)
        {
            grid.get(y, x).setNext(numstart);
            grid.get(y, x).setNwrap(1);
            active[numstart].setNwrap(1);
        }
        else
        {
            unsigned long tmp_next = grid.get(y, x).getNext();
            unsigned long tmp_last = tmp_next;
            unsigned long tmp_nwrap = 0;
            while(tmp_next != 0)
            {
                tmp_nwrap++;
                tmp_last = tmp_next;
                tmp_next = active[tmp_next].getNext();
            }
            grid.get(y, x).increaseNwrap();
            active[tmp_last].setNext(numstart);
            active[numstart].setNwrap(tmp_nwrap + 1);
        }
#ifdef DEBUG
        debugAddingLineage(numstart, x, y);
#endif
    }

    unsigned long SpatialTree::countCellExpansion(const long &x,
                                                  const long &y,
                                                  const long &xwrap,
                                                  const long &ywrap,
                                                  const double &generation_in,
                                                  vector<TreeNode> &data_added)
    {
        unsigned long map_cover = landscape->getVal(x, y, xwrap, ywrap, generation_in);
        unsigned long num_to_add = getIndividualsSampled(x, y, xwrap, ywrap, generation_in);
        double proportion_added = double(num_to_add) / double(map_cover);
        if(xwrap == 0 && ywrap == 0)
        {
            // Check that the species lineage_indices sizings make sense
            unsigned long ref = 0;
            if(map_cover != grid.get(y, x).getMaxSize())
            {
                if(map_cover > grid.get(y, x).getMaxSize())
                {
                    grid.get(y, x).changePercentCover(map_cover);
                }
                else
                {
                    grid.get(y, x).setMaxsize(map_cover);
                }
            }
            if(map_cover > grid.get(y, x).getListLength())
            {
                grid.get(y, x).changePercentCover(map_cover);
            }
            // Add the lineages
            while(ref < grid.get(y, x).getListLength() && num_to_add > 0)
            {
                unsigned long tmp_active = grid.get(y, x).getLineageIndex(ref);
                if(tmp_active != 0)
                {
                    if(checkProportionAdded(proportion_added))
                    {
                        makeTip(tmp_active, generation_in, data_added);
                        num_to_add--;
                    }
                }
                ref++;
            }
        }
        else
        {
            unsigned long next = grid.get(y, x).getNext();
            while(next != 0 && num_to_add > 0)
            {
                if(active[next].getXwrap() == xwrap && active[next].getYwrap() == ywrap)
                {
                    if(checkProportionAdded(proportion_added))
                    {
                        num_to_add--;
                        makeTip(next, generation_in, data_added);
                    }
                }
                next = active[next].getNext();
            }
        }
        return num_to_add;
    }

    void SpatialTree::expandCell(long x,
                                 long y,
                                 long x_wrap,
                                 long y_wrap,
                                 double generation_in,
                                 unsigned long num_to_add,
                                 vector<TreeNode> &data_added,
                                 vector<DataPoint> &active_added)
    {
        if(num_to_add > 0)
        {
            for(unsigned long k = 0; k < num_to_add; k++)
            {
                TreeNode tmp_tree_node{};
                DataPoint tmp_data_point{};
                unsigned long listpos = 0;
                // Add the species to active
                if(x_wrap == 0 && y_wrap == 0)
                {
                    listpos = grid.get(y, x).addSpecies(endactive + active_added.size() + 1);
                }
                tmp_data_point.setup(x, y, x_wrap, y_wrap, enddata + data_added.size() + 1, listpos, 1);
                if(enddata >= data->size())
                {
                    throw FatalException("Cannot add lineage - no space in data-> "
                                         "Check size calculations.");
                }
                if(endactive >= active.size())
                {
                    throw FatalException("Cannot add lineage - no space in active. "
                                         "Check size calculations.");
                }

                // Add a tip in the TreeNode for calculation of the coalescence tree at the end of the simulation.
                // This also contains the start x and y position of the species.
                tmp_tree_node.setup(true, x, y, x_wrap, y_wrap, generation_in);
                tmp_tree_node.setSpec(NR->d01());
                active_added.emplace_back(tmp_data_point);
                data_added.emplace_back(tmp_tree_node);

            }
        }
    }

    void SpatialTree::addGillespie(const double &g_threshold)
    {
        stringstream ss;
        ss << "Using gillespie algorithm in simulation from generation " << g_threshold << "." << endl;
        writeInfo(ss.str());
        gillespie_threshold = g_threshold;
        using_gillespie = true;
    }

    bool SpatialTree::runSimulationGillespie()
    {
        do
        {
            runSingleLoop();
        }
        while((endactive < gillespie_threshold) && (endactive > 1) &&
              ((steps < 100) || difftime(sim_end, start) < maxtime) && this_step.bContinueSim);
        // Switch to gillespie
        writeInfo("Switching to Gillespie algorithm.\n");
        setupGillespie();
        writeInfo("Starting Gillespie event loop...\n");
        do
        {
            runGillespieLoop();
        }
        while(endactive > 1);
        return stopSimulation();

    }

    void SpatialTree::runGillespieLoop()
    {

        // Decide what event and execute
        EventType next_event = heap.front().event_type;
        // Update the event timer
        steps += (heap.front().time_of_event - generation) * double(endactive);
        generation = heap.front().time_of_event;
        stringstream ss; //  TODO remove
        ss << "\tEvent at " << heap.front().time_of_event << "has " << active.size() - 1 << " lineages remaining...\n"
           << endl;
        writeInfo(ss.str());
        // Estimate the number of steps that have occurred.
        switch(next_event)
        {
        case EventType::cell_event:
        {
            GillespieProbability &origin = probabilities.get(heap.front().cell.y, heap.front().cell.x);
            gillespieCellEvent(origin);
            break;
        }

        case EventType::map_event:
            gillespieUpdateMap();
            break;

        case EventType::sample_event:
            gillespieSampleIndividuals();
            break;

        case EventType::undefined:
            throw FatalException("Undefined event in Gillespie algorithm. Please report this bug.");
        }

    }

    void SpatialTree::setupGillespie()
    {
        setupGillespieMaps();
        findLocations();
        updateAllProbabilities();
        createEventList();
        checkMapEvents();
        checkSampleEvents();
        sortEvents();
        // TODO add making dispersal map cumulative without the self-dispersal events
    }

    void SpatialTree::setupGillespieMaps()
    {
        if(dispersal_coordinator.isFullDispersalMap())
        {
            writeInfo("\tCreating cumulative dispersal map, excluding self-dispersal events...\n");
            self_dispersal_probabilities.setSize(sim_parameters->fine_map_y_size, sim_parameters->fine_map_x_size);
            dispersal_coordinator.reimportRawDispersalMap();
            for(unsigned long i = 0; i < sim_parameters->fine_map_y_size; i++)
            {
                for(unsigned long j = 0; j < sim_parameters->fine_map_x_size; j++)
                {
                    self_dispersal_probabilities.get(i, j) = dispersal_coordinator
                            .getSelfDispersalProbability(Cell(j, i));
                }
            }
            dispersal_coordinator.removeSelfDispersal();
        }
        probabilities.setSize(sim_parameters->fine_map_y_size, sim_parameters->fine_map_x_size);
    }

    Cell SpatialTree::getCellOfMapLocation(const MapLocation &location)
    {
        Cell cell{};
        cell.x = landscape->convertSampleXToFineX(location.x, location.xwrap);
        cell.y = landscape->convertSampleYToFineY(location.y, location.ywrap);
        return cell;
    }

    void SpatialTree::findLocations()
    {
        writeInfo("\tFinding all locations in the simulated world...\n");
        for(unsigned long y = 0; y < sim_parameters->fine_map_y_size; y++)
        {
            for(unsigned long x = 0; x < sim_parameters->fine_map_x_size; x++)
            {
                long x_pos = x;
                long y_pos = y;
                long x_wrap = 0;
                long y_wrap = 0;
                landscape->convertFineToSample(x_pos, x_wrap, y_pos, y_wrap);
                MapLocation location(x_pos, y_pos, x_wrap, y_wrap);
                addLocation(location);
            }
        }

    }

    void SpatialTree::checkMapEvents()
    {
        double next_map_update = 0.0;
        if(landscape->requiresUpdate())
        {
            next_map_update = sim_parameters->all_historical_map_parameters.front().generation;
            if(next_map_update > 0.0 && next_map_update < generation)
            {
                heap.emplace_back(GillespieHeapNode(generation, EventType::map_event));
            }
        }
    }

    void SpatialTree::checkSampleEvents()
    {
        for(const auto &item: reference_times)
        {
            // Find the first time that's after this point in time.
            if(item > generation)
            {
                heap.emplace_back(GillespieHeapNode(generation, EventType::sample_event));
            }
        }
    }

    void SpatialTree::gillespieCellEvent(GillespieProbability &origin)
    {
        CellEventType cell_event = origin.generateRandomEvent(NR);
        origin.setRandomNumber(NR->d01());
        switch(cell_event)
        {
        case CellEventType::coalescence_event:
            // implement coalescence
            writeInfo("Coalescence event."); // TODO remove
            gillespieCoalescenceEvent(origin);
            break;

        case CellEventType::dispersal_event:
            // choose dispersal
            writeInfo("Dispersal event."); // TODO remove
            gillespieDispersalEvent(origin);
            break;

        case CellEventType::speciation_event:
            writeInfo("Speciation event."); // TODO remove
            gillespieSpeciationEvent(origin);
            break;

        case CellEventType::undefined:
            throw FatalException("Undefined cell event type. Please report this bug.");
            break;
        }

    }

    void SpatialTree::gillespieUpdateMap()
    {
        // First delete all existing objects
        clearGillespieObjects();

        // Update the existing landscape structure
        if(landscape->updateMap(generation))
        {
            dispersal_coordinator.updateDispersalMap();
            // Update all the heap variables
            findLocations();
            updateAllProbabilities();
            createEventList();

            // Now need to get the next update map event.
            checkMapEvents();
            checkSampleEvents();
            sortEvents();
        }
        else
        {
            stringstream ss;
            ss << "Didn't update map at generation " << generation
               << ". Incorrect placement of map_event on events queue. Please report this bug." << endl;
            throw FatalException(ss.str());
        }

    }

    void SpatialTree::gillespieSampleIndividuals()
    {
        clearGillespieObjects();
        addLineages(generation);
        findLocations();
        updateAllProbabilities();
        createEventList();
        checkMapEvents();
        checkSampleEvents();
        sortEvents();
    }

    void SpatialTree::gillespieCoalescenceEvent(GillespieProbability &origin)
    {
        auto lineages = selectTwoRandomLineages(origin.getMapLocation());
        gillespieUpdateGeneration(lineages.first);

        //        stringstream ss; // TODO remove
        //        ss << "Location at " << origin.getMapLocation().x << ", " << origin.getMapLocation().y << endl;
        //        ss << "Coalescing lineages at " << lineages.first << " and " << lineages.second << endl;
        //        writeInfo(ss.str());

        if(lineages.first > active.size() || lineages.second > active.size())
        {
            throw FatalException("Lineage indexing incorrect. Please report this bug."); // TODO remove
        }

        coalescenceEvent(lineages.first, lineages.second);

        const MapLocation &location = origin.getMapLocation();
        updateCellCoalescenceProbability(origin, getNumberIndividualsAtLocation(location));
        updateInhabitedCellOnHeap(convertMapLocationToCell(location));
    }

    void SpatialTree::gillespieDispersalEvent(GillespieProbability &origin)
    {
        // Select a lineage in the target cell.
        unsigned long chosen = selectRandomLineage(origin.getMapLocation());
        // Sets the old location for the lineage and zeros out the coalescence stuff
        recordLineagePosition();
        // Remove the chosen lineage from the cell
        removeOldPosition(chosen);
        // Performs the move and calculates any coalescence events
        calcNextStep();
        // Get the destination cell and update the probabilities.
        Cell destination_cell(convertMapLocationToCell(active[chosen]));
        auto x = destination_cell.x;
        auto y = destination_cell.y;
        unsigned long n = getNumberLineagesAtLocation(origin.getMapLocation());
        if(n > 0)
        {
            updateCellCoalescenceProbability(origin, getNumberIndividualsAtLocation(origin.getMapLocation()));
            updateInhabitedCellOnHeap(destination_cell);
        }
        else
        {
            removeHeapTop();
            cellToHeapPositions.get(y, x) = SpatialTree::UNUSED;
        }

        GillespieProbability &destination = probabilities.get(y, x);
        if(cellToHeapPositions.get(y, x) == SpatialTree::UNUSED)
        {
            addNewEvent(x, y);
        }
        else if(!this_step.coal)
        {
            // Needs to update destination
            destination.setCoalescenceProbability(calculateCoalescenceProbability(destination.getMapLocation()));
            const double local_death_rate = getLocalDeathRate(active[chosen]);
            destination.setRandomNumber(NR->d01());
            const double t = destination.calcTimeToNextEvent(local_death_rate,
                                                             summed_death_rate,
                                                             getNumberIndividualsAtLocation(destination
                                                                                                    .getMapLocation()));
            heap[cellToHeapPositions.get(y, x)].time_of_event = t;
            updateInhabitedCellOnHeap(destination_cell);
        }

    }

    void SpatialTree::gillespieSpeciationEvent(GillespieProbability &origin)
    {
        const MapLocation &location = origin.getMapLocation();
        unsigned long chosen = selectRandomLineage(location);
        gillespieUpdateGeneration(chosen);
        speciateLineage(active[chosen].getReference());
        updateCellCoalescenceProbability(origin, getNumberIndividualsAtLocation(location));
        updateInhabitedCellOnHeap(convertMapLocationToCell(origin.getMapLocation()));
    }

    template<typename T> const double SpatialTree::getLocalDeathRate(const T &location)
    {
        const Cell cell = convertMapLocationToCell(location);
        if(death_map->isNull())
        {
            return 1.0;
        }
        else
        {
            return death_map->get(cell.y, cell.x);
        }
    }

    template<typename T> const double SpatialTree::getLocalSelfDispersalRate(const T &location)
    {
        const Cell cell = convertMapLocationToCell(location);
        if(dispersal_coordinator.isFullDispersalMap())
        {
            return 1.0;
        }
        return self_dispersal_probabilities.get(cell.y, cell.x);
    }

    void SpatialTree::clearGillespieObjects()
    {
        cellToHeapPositions.fill(0);
        heap.clear();
        for(auto &item : probabilities)
        {
            item.reset();
        }
    }

    void SpatialTree::gillespieUpdateGeneration(const unsigned long &lineage)
    {// TODO remove
        if(lineage == 0 || lineage > active.size())
        {
            stringstream ss;
            ss << "Lineage " << lineage << " out of range of active." << endl;
            throw FatalException(ss.str());
        }
        TreeNode &tree_node = (*data)[active[lineage].getReference()];
        unsigned long generations_existed = round(generation) - tree_node.getGeneration();
        tree_node.setGeneration(generations_existed);
    }

    void SpatialTree::updateCellCoalescenceProbability(GillespieProbability &origin, const unsigned long &n)
    {
        const MapLocation &location = origin.getMapLocation();
        origin.setCoalescenceProbability(calculateCoalescenceProbability(location));
        origin.setRandomNumber(NR->d01());
        heap.front().time_of_event = origin.calcTimeToNextEvent(getLocalDeathRate(location), summed_death_rate, n);
    }

    void SpatialTree::updateInhabitedCellOnHeap(const Cell &pos)
    {
        eastl::change_heap(heap.begin(), heap.size(), cellToHeapPositions.get(pos.y, pos.x));
        
        gillespieValidateHeap(); // TODO remove
    }

    void SpatialTree::gillespieValidateHeap() // TODO remove
    {
        if(!eastl::is_heap(heap.begin(), heap.end()))
        {
            throw FatalException("The heap property has been broken. Please report this bug."); // TODO remove
        }
        
        for (size_t i = 0; i < heap.size(); i++)
        {
            if(*(heap[i].locator) != i)
            {
                throw FatalException("The heap locator has been broken. Please report this bug."); // TODO remove
            }
        }
    }

    void SpatialTree::updateAllProbabilities()
    {
        writeInfo("\tCalculating global mean death rate and total number of individuals...\n");

        if(!death_map->isNull())
        {
            summed_death_rate = 0.0;
            for(unsigned long y = 0; y < sim_parameters->fine_map_y_size; y++)
            {
                for(unsigned long x = 0; x < sim_parameters->fine_map_y_size; x++)
                {
                    summed_death_rate += death_map->get(y, x) * landscape->getValFine(x, y, generation);
                }
            }
        }
        else
        {
            // calculate global death rate mean if death rates are equal
            summed_death_rate = std::accumulate(landscape->getFineMap().begin(),
                                                landscape->getFineMap().end(),
                                                (unsigned long) 0);
        }
    }

    void SpatialTree::removeHeapTop()
    {
        eastl::pop_heap(heap.begin(), heap.end());
        heap.pop_back();

        gillespieValidateHeap(); // TODO remove
    }

    void SpatialTree::createEventList()
    {
        writeInfo("\tAdding events to event list...\n");
        cellToHeapPositions.setSize(sim_parameters->fine_map_y_size, sim_parameters->fine_map_x_size);
        cellToHeapPositions.fill(SpatialTree::UNUSED);

        for(unsigned long y = 0; y < sim_parameters->fine_map_y_size; y++)
        {
            for(unsigned long x = 0; x < sim_parameters->fine_map_x_size; x++)
            {
                addNewEvent < false > (x, y);
            }
        }
    }

    void SpatialTree::sortEvents()
    {
        eastl::make_heap(heap.begin(), heap.end());
        
        gillespieValidateHeap(); // TODO remove
    }

    template<bool restoreHeap = true> void SpatialTree::addNewEvent(const unsigned long &x, const unsigned long &y)
    {
        const MapLocation &location = probabilities.get(y, x).getMapLocation();
        if(getNumberLineagesAtLocation(location) > 0)
        {
            cellToHeapPositions.get(y, x) = heap.size();
            
            heap.emplace_back(GillespieHeapNode(Cell(x, y), (probabilities.get(y, x).calcTimeToNextEvent(
                    getLocalDeathRate(location), summed_death_rate, getNumberIndividualsAtLocation(location)) +
                    generation), EventType::cell_event, &heap, &cellToHeapPositions.get(y, x)));
            
            if(restoreHeap)
            {
                eastl::push_heap(heap.begin(), heap.end());
            }
        }
    }

    void SpatialTree::addLocation(const MapLocation &location)
    {
        Cell cell = getCellOfMapLocation(location);
        GillespieProbability gp(location);
        // check if any lineages exist there

        gp.setDispersalOutsideCellProbability(1.0 - getLocalSelfDispersalRate(location));
        gp.setCoalescenceProbability(calculateCoalescenceProbability(location));
        gp.setSpeciationProbability(spec);
        gp.setRandomNumber(NR->d01());

        probabilities.get(static_cast<const unsigned long &>(cell.y), static_cast<const unsigned long &>(cell.x)) = gp;
    }

    double SpatialTree::calculateCoalescenceProbability(const MapLocation &location)
    {
        unsigned long max_number_individuals = landscape
                ->getVal(location.x, location.y, location.xwrap, location.ywrap, generation);
        unsigned long current_number = getNumberLineagesAtLocation(location);
        if(current_number == 1)
        {
            return 0.0;
        }
        return min((double(current_number) - 1.0) / double(max_number_individuals), 1.0);

    }

    unsigned long SpatialTree::selectRandomLineage(const MapLocation &location) const
    {
        vector<unsigned long> lineage_ids = detectLineages(location);
        unsigned long random_index = NR->i0(lineage_ids.size() - 1);
        return lineage_ids[random_index];
    }

    pair<unsigned long, unsigned long> SpatialTree::selectTwoRandomLineages(const MapLocation &location) const
    {
        vector<unsigned long> lineage_ids = detectLineages(location);
        if(lineage_ids.size() < 2)
        {
            throw FatalException("Cannot select two lineages when fewer than two exist at location.");
        }
        pair<unsigned long, unsigned long> selected_lineages;
        selected_lineages.first = lineage_ids[NR->i0(lineage_ids.size() - 1)];

        do
        {
            selected_lineages.second = lineage_ids[NR->i0(lineage_ids.size() - 1)];
        }
        while(selected_lineages.second == selected_lineages.first);
        return selected_lineages;
    }

    vector<unsigned long> SpatialTree::detectLineages(const MapLocation &location) const
    {
        const SpeciesList &species_list = grid.get(location.y, location.x);
        vector<unsigned long> lineage_ids;
        if(location.isOnGrid())
        {
            lineage_ids.reserve(species_list.getListSize());
            for(unsigned long i = 0; i < species_list.getListLength(); i++)
            {
                unsigned long lineage_index = species_list.getLineageIndex(i);
                if(lineage_index != 0)
                {
                    lineage_ids.push_back(lineage_index);
                    if(lineage_ids.size() > species_list.getListSize())
                    {
                        break;
                    }
                }
            }
        }
        else
        {
            lineage_ids.reserve(species_list.getNwrap());
            unsigned long next = species_list.getNext();
            do
            {
                const DataPoint &datapoint = active[next];
                if(datapoint == location)
                {
                    lineage_ids.push_back(next);
                }
                next = active[next].getNext();
            }
            while(next != 0);
        }
        // TODO remove or move to DEBUG
        for(const auto &item: lineage_ids)
        {
            if(item == 0)
            {
                stringstream ss;
                ss << "Lineages not correctly calculated for location " << location.x << ", " << location.y << "("
                   << location.xwrap << ", " << location.ywrap << ")" << endl;
                throw FatalException(ss.str());
            }
        }
        return lineage_ids;

    }

#ifdef DEBUG

    void SpatialTree::validateLineages()
    {
        bool fail = false;
        writeInfo("\nStarting lineage validation...");
#ifdef historical_mode
        unsigned long printed = 0;
#endif // historical_mode
        // Basic checks
        if(endactive >= active.size() || enddata >= data->size())
        {
            stringstream ss;
            ss << "Endactive (size):" << endactive << "(" << active.size() << ")" << endl;
            ss << "Enddata (size):" << enddata << "(" << data->size() << ")" << endl;
            writeCritical(ss.str());
            throw FatalException("Endactive out of range of active or enddata out of range of data-> "
                                 "Please report this bug.");
        }
        for(unsigned long i = 1; i < endactive; i++)
        {
            stringstream ss;
            DataPoint tmp_datapoint = active[i];
            // Validate the location exists
#ifdef historical_mode
            if(landscape->getVal(tmp_datapoint.getXpos(), tmp_datapoint.getYpos(),
                                tmp_datapoint.getXwrap(), tmp_datapoint.getYwrap(), generation) == 0)
            {
                if(printed < 100)
                {
                    printed++;
                    ss << "Map value: " << landscape->getVal(tmp_datapoint.getXpos(), tmp_datapoint.getYpos(),
                                                            tmp_datapoint.getXwrap(), tmp_datapoint.getYwrap(),
                                                            generation) << endl;
                }
                fail = true;
            }
#endif // historical mode
            if(tmp_datapoint.getXwrap() == 0 && tmp_datapoint.getYwrap() == 0)
            {
                if(tmp_datapoint.getNwrap() != 0)
                {
                    fail = true;
                }
                else
                {
                    if(i !=
                         grid.get(tmp_datapoint.getYpos(), tmp_datapoint.getXpos()).getLineageIndex(tmp_datapoint.getListpos()))
                    {
                        fail = true;
                    }
                }
            }
            else
            {
                if(tmp_datapoint.getNwrap() == 0)
                {
                    fail = true;
                }
                else
                {
                    unsigned long tmp_next = grid.get(tmp_datapoint.getYpos(), tmp_datapoint.getXpos()).getNext();
                    unsigned long count = 0;
                    while(tmp_next != 0)
                    {
                        count++;
                        if(count != active[tmp_next].getNwrap())
                        {
                            ss << "problem in wrap: " << count << " != " << active[tmp_next].getNwrap() << endl;
                            fail = true;
                        }
                        tmp_next = active[tmp_next].getNext();
                    }
                    if(count == 0 && count != grid.get(tmp_datapoint.getYpos(), tmp_datapoint.getXpos()).getNwrap())
                    {
                        fail = true;
                    }
                    if(count != grid.get(tmp_datapoint.getYpos(), tmp_datapoint.getXpos()).getNwrap())
                    {
                        fail = true;
                    }
                }
            }
            if(fail)
            {
                stringstream ss;
                ss << "Active reference: " << i << endl;
                ss << "Grid wrapping: " << grid.get(tmp_datapoint.getYpos(), tmp_datapoint.getXpos()).getNwrap() << endl;
                ss << "Endactive: " << endactive << endl;
                ss << "Active size: " << active.size() << endl;
                ss << "Enddata: " << enddata << endl;
                ss << "Data size: " << data->size() << endl;
                writeLog(50, ss);
                tmp_datapoint.logActive(50);
                (*data)[tmp_datapoint.getReference()].logLineageInformation(50);
                throw FatalException("Failure in lineage validation. Please report this bug.");
            }
        }
        writeInfo("done.\n");
    }

    void SpatialTree::debugAddingLineage(unsigned long numstart, long x, long y)
    {
        unsigned long tmp_next = grid.get(y, x).getNext();
        unsigned long tmp_nwrap = 0;
        while(tmp_next != 0)
        {
            tmp_nwrap++;
            if(active[tmp_next].getNwrap() != tmp_nwrap)
            {
                stringstream ss;
                ss << "tmp_nwrap: " << tmp_nwrap << endl;
                ss << "next = " << tmp_next << endl;
                ss << "numstart: " << numstart << endl;
                writeLog(50, ss);
                active[tmp_nwrap].logActive(50);
                throw FatalException("Incorrect setting of nwrap in wrapped lineage, please report this bug.");
            }
            tmp_next = active[tmp_next].getNext();
        }
        if(tmp_nwrap != grid.get(y, x).getNwrap())
        {
            stringstream ss;
            ss << "Grid nwrap: " << grid.get(y, x).getNwrap() << endl;
            ss << "Counted wrapping: " << tmp_nwrap << endl;
            ss << "active: " << numstart << endl;
            tmp_next = grid.get(y, x).getNext();
            tmp_nwrap = 0;
            while(tmp_next != 0 && tmp_nwrap < grid.get(y, x).getNwrap())
            {
                tmp_nwrap++;
                ss << "tmp_next: " << tmp_next << endl;
                ss << "tmp_nwrap: " << tmp_nwrap << endl;
                tmp_next = active[tmp_next].getNext();
            }
            writeLog(50, ss);
            throw FatalException("Grid wrapping value not set correctly");
        }
    }

    void SpatialTree::runChecks(const unsigned long &chosen, const unsigned long &coalchosen)
    {
        // final checks
#ifdef historical_mode
        if(active[chosen].getListpos() > grid.get(active, chosen).getYpos()][active[chosen].getXpos()].getMaxsize() &&
           active[chosen].getNwrap() == 0)
        {
            throw FatalException("ERROR_MOVE_001: Listpos outside maxsize.");
        }

        if(active[coalchosen].getListpos() >
               grid.get(active, coalchosen).getYpos()][active[coalchosen].getXpos()].getMaxsize() &&
           active[coalchosen].getNwrap() == 0 && coalchosen != 0)
        {
            throw FatalException("Coalchosen list_position outside maxsize. Please report this bug.");
        }
#endif
        Tree::runChecks(chosen, coalchosen);
        if(active[chosen].getNwrap() != 0)
        {
            unsigned long tmpactive = grid.get(active[chosen].getYpos(), active[chosen].getXpos()).getNext();
            for(unsigned long i = 1; i < active[chosen].getNwrap(); i++)
            {
                tmpactive = active[tmpactive].getNext();
            }

            if(tmpactive != chosen)
            {
                active[chosen].logActive(50);
                throw FatalException("ERROR_MOVE_003: Nwrap not set correctly.");
            }
        }

        if(active[chosen].getNwrap() != 0)
        {
            if(active[chosen].getXwrap() == 0 && active[chosen].getYwrap() == 0)
            {
                throw FatalException("ERROR_MOVE_10: Nwrap set to non-zero, but x and y wrap 0.");
            }
        }
        if(active[endactive].getNwrap() != 0)
        {
            unsigned long nwrap = active[endactive].getNwrap();
            if(nwrap == 1)
            {
                if(grid.get(active[endactive].getYpos(), active[endactive].getXpos()).getNext() != endactive)
                {
                    stringstream ss;
                    ss << "Lineage at 1st position: "
                       << grid.get(active[endactive].getYpos(), active[endactive].getXpos()).getNext() << endl;
                    ss << "endactive: " << endactive << endl
                       << "nwrap: " << nwrap << endl;
                    ss << "chosen: " << chosen << endl;
                    writeLog(10, ss);
                    throw FatalException("ERROR_MOVE_016: Nwrap for endactive not set correctly. Nwrap is 1, "
                                         "but the lineage at 1st position is not endactive.");
                }
            }
            else
            {
                unsigned long tmpcheck = grid.get(active[endactive].getYpos(), active[endactive].getXpos()).getNext();
                unsigned long tmpnwrap = 1;
                while(tmpcheck != endactive)
                {
                    tmpnwrap++;
                    tmpcheck = active[tmpcheck].getNext();
                    if(tmpnwrap > nwrap + 1)
                    {
                        stringstream ss;
                        ss << "ERROR_MOVE_017: NON FATAL. Nrap for endactive not set correctly; looped "
                              "beyond nwrap and not yet found enactive."
                           << endl;
                        ss << "endactive: " << endactive << endl
                           << "nwrap: " << nwrap << endl
                           << "x,y: " << active[endactive].getXpos() << "," << active[endactive].getYpos()
                           << endl;
                        ss << "chosen: " << chosen << endl;
                        writeLog(10, ss);
                    }
                }
                if(tmpnwrap != nwrap)
                {
                    stringstream ss;
                    ss << "ERROR_MOVE_018: NON FATAL. Nwrap for endactive not set correctly. Nwrap is "
                       << nwrap << " but endactive is at position " << tmpnwrap << endl;
                    ss << "endactive: " << endactive << endl
                       << "nwrap: " << nwrap << endl
                       << "x,y: " << active[endactive].getXpos() << "," << active[endactive].getYpos()
                       << endl;
                    ss << "chosen: " << chosen << endl;
                    writeLog(10, ss);
                }
            }
        }
    }



#endif
};
