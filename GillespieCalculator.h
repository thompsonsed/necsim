//
// Created by sam on 07/09/19.
//

#ifndef NECSIM_GILLESPIECALCULATOR_H
#define NECSIM_GILLESPIECALCULATOR_H

#include <map>
#include <algorithm>

#include "heap.hpp"

#include "ActivityMap.h"
#include "Cell.h"
#include "DispersalCoordinator.h"
#include "RNGController.h"
#include "DataPoint.h"
#include "SpatialTree.h"

/**
 * @brief Container for the different event types that can occur during the Gillespie Algorithm.
 */
enum class EventType
{
    undefined, cell_event, map_update, sample_time
};

enum class CellEventType
{
    dispersal_event, coalescence_event, speciation_event
};

class GillespieProbability
{
protected:
    double dispersal_outside_cell_probability;
    double coalescence_probability;
    double speciation_probability;
    double cell_turnover_probability;
    double random_number;
    EventType event_type;
    MapLocation location;
    unsigned long *pos;

public:

    GillespieProbability(MapLocation c, unsigned long *pos) : dispersal_outside_cell_probability(0.0),
                                                              coalescence_probability(0.0),
                                                              speciation_probability(0.0),
                                                              random_number(0.0), location(c), time_to_event(0.0),
                                                              event_type(undefined), pos(pos)
    {

    }

    void setDispersalOutsideCellProbability(const double &d)
    {
        this->dispersal_outside_cell_probability = d;
    }

    void setCoalescenceProbability(const double &c)
    {
        this->coalescence_probability = c;
    }

    void setSpeciationProbability(const double &s)
    {
        this->speciation_probability = speciation_probability;
    }

    void setRandomNumber(const double &r)
    {
        random_number = r;
    }

    double getInCellProbability() const
    {
        return speciation_probability + (1 - speciation_probability) *
                                        ((1 - dispersal_outside_cell_probability) * coalescence_probability +
                                         dispersal_outside_cell_probability);
    }

    EventType &getEventType() const
    {
        return event_type;
    }

    CellEventType &generateRandomEvent(shared_ptr<RNGCoordinator> rng) const
    {
        double p = rng.d01() * getInCellProbability();
        if(p < speciation_probability)
        {
            return speciation_event;
        }
        else
        {
            if(p < speciation_probability + (1 - speciation_probability) * dispersal_probability)
            {
                return dispersal_event;
            }
            else
            {
                return coalescence_event;
            }
        }
    }

    MapLocation &getMapLocation() const
    {
        return location;
    }

    /**
     * @brief Gets the parameter for the exponential function
     *
     * The rate is per birth-death event on the whole landscape.
     *
     * @param mean_death_rate
     * @param n
     * @param total_n
     * @return lambda
     */
    double getLambda(const double &mean_death_rate, const unsigned long &n, const unsigned long &total_n) const
    {
        return getInCellProbability() * mean_death_rate * double(n) / double(total_n);
    }

    void setEvent(EventType event)
    {
        event_type = event;
    }

    double calcTimeToNextEvent(const double &mean_death_rate, const unsigned long &n, const unsigned long &total_n)
    {
        return RNGController::exponentialDistribution(getLambda(mean_death_rate, n, total_n), random_number);
    }

    /*bool operator<(const GillespieProbability &gp) const
    {
        return (time_to_event < gp.time_to_event);
    }*/
};

class GillespieHeapNode
{
public:
    Cell cell;
    double timeOfEvent;
    // Pointer to index in heap
    unsigned long *pos;

    GillespieHeapNode(Cell cell, double timeOfEvent, unsigned long *pos)
            : cell(cell), timeOfEvent(timeOfEvent), pos(pos){}

    GillespieHeapNode() : GillespieHeapNode(Cell(), 0.0, nullptr){}
};

bool operator<(const GillespieHeapNode &n1, const GillespieHeapNode &n2)
{
    return n1.timeOfEvent > n2.timeOfEvent;
}

namespace std
{
    template<>
    void swap<GillespieHeapNode>(GillespieHeapNode &lhs, GillespieHeapNode &rhs)
    {
        std::swap(lhs.cell, rhs.cell);
        std::swap(lhs.timeOfEvent, rhs.timeOfEvent);
        std::swap(*lhs.pos, *rhs.pos);
    }
}

class GillespieCalculator : public SpatialTree
{
protected:
    // Objects from the simulation that are required by the Gillespie algorithm
    // Map of dispersal probabilities is stored in the dispersal coordinator
    shared_ptr<DispersalCoordinator> dispersal_coordinator;
    // Map of reproduction rates
    shared_ptr<ActivityMap> reproduction_map;
    // Map of death rates
    shared_ptr<ActivityMap> death_map;
    // The random number generator
    shared_ptr<RNGController> random;
    // Landscape containing the spatial referencing between maps
    shared_ptr<Landscape> landscape;
    // Contains the coalescence tree, will be some combination of spatial or non-spatial tree and protracted or point
    // speciation

    // The active lineages in the simulation
    shared_ptr<vector<DataPoint>> active;

    Matrix<GillespieProbability> probabilities;
    vector<GillespieHeapNode> heap;
    // Index to heap position, or UNUSED if cell is not used.
    Matrix<unsigned long> cellToHeapPositions;

    unsigned long global_individuals;
    double mean_death_rate;

    static const unsigned long UNUSED = static_cast<long>(-1);

public:
    GillespieCalculator() : locations(), probabilities(){}

    /**
     * @brief
     */
    void runGillespie()
    {
        setupGillespie();
        do
        {
            runSingleLoop(); // TODO update this to the gillespie loop
        }
        while(endactive > 1);

    }

    void runGillespieLoop()
    {
        GillespieProbability &origin = probabilities[heap.front().gillespieProbabilityIndex];

        // Decide what event and execute
        EventType next_event = origin.getEventType();
        switch(next_event)
        {
            case cell_event:
                gillespieCellEvent(origin);
                break;

            case map_event:
                gillespieUpdateMap();

//                if(landscape->updateMap(generation))
//                {
//                    dispersal_coordinator.updateDispersalMap();
//                }
                break;

            case sample_event:
                gillespieSampleIndividuals();
                break;

            case undefined:
                throw FatalException("Undefined event in Gillespie algorithm. Please report this bug.");
        }

    }

    void gillespieCellEvent(GillespieProbability &origin)
    {
        CellEventType cell_event = origin.generateRandomEvent(NR);
        switch(cell_event)
        {
            case coalescence_event:
                // implement coalescence
                gillespieCoalescenceEvent(origin);
                break;

            case dispersal_event:
                // choose dispersal
                gillespieDispersalEvent();
                break;

            case speciation_event:
                break;
        }

    }

    void gillespieUpdateMap()
    {

    }

    void gillespieSampleIndividuals()
    {

    }

    void gillespieCoalescenceEvent(GillespieProbability &origin)
    {

    }

    void gillespieDispersalEvent(GillespieProbability &origin)
    {
        // Select a lineage in the target cell.
        unsigned long chosen = selectLineage(origin);
        // Sets the old location for the lineage and zeros out the coalescence stuff
        recordLineagePosition();
        // Remove the chosen lineage from the cell
        removeOldPosition(chosen);
        // Performs the move and calculates any coalescence events
        calcNextStep();
        checkOriginCellInhabited();
        // Get the destination cell and update the probabilities.
        unsigned long x = active[chosen].getXpos();
        unsigned long y = active[chosen].getYpos();
        Cell destination_cell(x, y);
        unsigned long n = getNumberIndividualsAtLocation(origin.getMapLocation());
        if(n > 0)
        {
            origin.setCoalescenceProbability(calculateCoalescenceProbability(origin.getMapLocation()));
            origin.setRandomNumber(NR.d01());
            heap.front().timeOfEvent = origin.calcTimeToNextEvent(mean_death_rate, n, global_individuals);
            updateInhabitedCellOnHeap(destination_cell);
        }
        else
        {
            std::pop_heap(heap.begin(), heap.end());
            heap.pop_back();
            cellToHeapPositions.get(y, x) = GillespieCalculator::UNUSED;
        }

        GillespieProbability &destination = probabilities.get(y, x);
        if(cellToHeapPositions.get(y, x) == GillespieCalculator::UNUSED)
        {
            addNewEvent(x, y);
        }
        else if(!this_step.coal)
        {
            // Needs to update destination
            destination.setCoalescenceProbability(calculateCoalescenceProbability(destination.getMapLocation()));
            // TODO check if the new random number definitely doesn't need to be generated
            heap[cellToHeapPositions.get(y, x)].timeOfEvent =
                    destination.calcTimeToNextEvent(mean_death_rate,
                                                    getNumberIndividualsAtLocation(
                                                            destination.getMapLocation()),
                                                    global_individuals);
            updateInhabitedCellOnHeap(destination_cell);
        }

    }

    void updateHeapOrigin(GillespieProbability &origin)
    {
        // If cell is still inhabited
        {
            Cell pos = convertMapLocationToXY(origin.getMapLocation());

            // Update probability and times of origin
            updateInhabitedCellOnHeap(pos);

        }
        // Else
        {
            std::pop_heap(heap.begin(), heap.end());
            heap.pop_back();
        }
    }

    void updateHeapDestination(GillespieProbability &destination)
    {
        // If dispersal to other cell

        Cell pos = convertMapLocationToXY(destination.getMapLocation());

        // If cell is already inhabited
        {
            // Update probability and times of destination
            updateInhabitedCellOnHeap(pos);
        }
        // Else
        {
            heap.emplace_back(GillespieHeapNode(pos,
                                                destination.calcTimeToNextEvent(mean_death_rate,
                                                                                getNumberIndividualsAtLocation(
                                                                                        destination.getMapLocation()),
                                                                                global_individuals) +
                                                generation),
                              &cellToHeapPositions.get(pos.y, pos.x));
            std::push_heap(heap.begin(), heap.end());
        }
    }

    void updateInhabitedCellOnHeap(const Cell &pos)
    {
        std::update_heap(heap.begin(), heap.end(),
                         heap.begin() + cellToHeapPositions.get(pos.y, pos.x));
    }

    /**
     * @brief Sets up the Gillespie algorithm.
     */
    void setupGillespie()
    {
        findLocations();
        updateAllProbabilities();
        createEventList();
        // TODO add making dispersal map cumulative without the self-dispersal events
    }

    /**
     * @brief Calculates the x, y position on the fine map of the lineage
     * @param location the map location
     * @return cell object containing the x, y location
     */
    Cell getCellOfMapLocation(const MapLocation &location)
    {
        Cell cell{};
        cell.x = landscape->convertSampleXToFineX(location.x, location.xwrap);
        cell.y = landscape->convertSampleYToFineY(location.y, location.ywrap);
        return cell;
    }

    /**
     * @brief Finds the locations that lineages are at and adds them to the list of locations.
     *
     * This also involves calculating the event probabilities for each cell.
     */
    void findLocations()
    {
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

    void updateAllProbabilities()
    {
        // calculate global death rate mean
        calculateGlobalDeathRate();
        // calculate total number of individuals
        calculateGlobalIndividuals();
    }

    Cell convertMapLocationToXY(const MapLocation &location)
    {
        unsigned long x = landscape->convertSampleXToFineX(location.x, location.xwrap);
        unsigned long y = landscape->convertSampleXToFineX(location.y, location.ywrap);

        return Cell(x, y);
    }

    /**
     * @brief Calculates the times for each event and sorts the event list.
     */
    void createEventList()
    {
        cellToHeapPositions.setSize(sim_parameters->fine_map_y_size, sim_parameters->fine_map_x_size);
        cellToHeapPositions.fill(GillespieCalculator::UNUSED);

        for(unsigned long y = 0; y < sim_parameters->fine_map_y_size; y++)
        {
            for(unsigned long x = 0; x < sim_parameters->fine_map_x_size; x++)
            {
                addNewEvent(x, y);
            }
        }

        std::make_heap(heap.begin(), heap.end());
    }

    void addNewEvent(unsigned long x, unsigned long y)
    {
        local_individuals = getNumberIndividualsAtLocation(
                probabilities.get(y, x).getMapLocation());
        if(local_individuals > 0)
        {
            cellToHeapPositions.get(y, x) = heap.size();
            heap.emplace_back(
                    GillespieHeapNode(Cell(x, y),
                                      probabilities.get(y, x).calcTimeToNextEvent(
                                              mean_death_rate,
                                              local_individuals,
                                              global_individuals) +
                                      generation),
                    &cellToHeapPositions.get(y, x));
        }
    }

    /**
     * @brief Adds the given location.
     *
     * Calculates the probabilities of coalescence, dispersal and speciation.
     * @param location the location to add and calculate values for
     */
    void addLocation(const MapLocation &location)
    {
        Cell cell = getCellOfMapLocation(location);

        GillespieProbability gp = GillespieProbability(location);
        // check if any lineages exist there

        gp.setDispersalProbability(1.0 - dispersal_coordinator->getSelfDispersalProbability(cell));
        gp.setCoalescenceProbability(calculateCoalescenceProbability(location));
        gp.setSpeciationProbability(spec);
        gp.setRandomNumber(NR->d01());

        probabilities.get(cell.x, cell.y) = gp;
    }

    double calculateCoalescenceProbability(const MapLocation &location)
    {
        unsigned long max_number_individuals = landscape->getVal(location.x, location.y, location.xwrap, location.ywrap,
                                                                 generation);
        unsigned long current_number = getNumberIndividualsAtLocation(location);
        return min(double(current_number) / double(max_number_individuals), 1.0);

    }

};

#endif //NECSIM_GILLESPIECALCULATOR_H
