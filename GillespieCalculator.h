//
// Created by sam on 07/09/19.
//

#ifndef NECSIM_GILLESPIECALCULATOR_H
#define NECSIM_GILLESPIECALCULATOR_H

#include "ActivityMap.h"
#include "Cell.h"
#include "DispersalCoordinator.h"
#include "RNGController.h"
#include "DataPoint.h"
#include "SpatialTree.h"
#include <map>

/**
 * @brief Container for the different event types that can occur during the Gillespie Algorithm.
 */
enum EventType {undefined, cell_event, map_update, sample_time};

class GillespieProbabilities
{
protected:
    long double dispersal_probability;
    long double coalescence_probability;
    long double speciation_probability;
    double cell_turnover_probability;
    long double random_number;
    MapLocation location;
    double time_to_event;
    EventType event_type;


public:

    GillespieProbabilities(MapLocation c) : dispersal_probability(0.0), coalescence_probability(0.0),
                                            speciation_probability(0.0),
                                            random_number(0.0), location(c), time_to_event(0.0), event_type(undefined)
    {


    }

    void setDispersalProbability(const long double &d)
    {
        this->dispersal_probability = d;
    }

    void setCoalescenceProbability(const long double &c)
    {
        this->coalescence_probability = c;
    }

    void setSpeciationProbability(const long double &s)
    {
        this->speciation_probability = speciation_probability;
    }

    void setRandomNumber(const long double &r)
    {
        random_number = r;
    }

    double getTotalProbability() const
    {
        return cell_turnover_probability * (dispersal_probability + coalescence_probability + speciation_probability);
    }

    void setEvent(double t, EventType event)
    {
        time_to_event = t;
        event_type = event;
    }

    void calclTimeToNextEvent()
    {
        setEvent(RNGController::exponentialDistribution(getTotalProbability(), random_number), cell_event);
    }

    double getTimeToNextEvent()
    {
        return time_to_event;
    }

    bool operator < (const GillespieProbabilities & gp) const
    {
        return (time_to_event < gp.time_to_event);
    }
};

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
    // Store the list of all locations which contain individuals and the respective probabilities.
    set<MapLocation> locations;
    vector<GillespieProbabilities> probabilities;


public:
    GillespieCalculator() : locations(), probabilities()
    { }

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

    /**
     * @brief Sets up the Gillespie algorithm.
     */
    void setupGillespie()
    {
        findLocations();
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
        cell.y = landscape->convertSampleYToFineY(location.y, location.ywrap));
        return cell;
    }



    /**
     * @brief Finds the locations that lineages are at and adds them to the list of locations.
     *
     * This also involves calculating the event probabilities for each cell.
     */
    void findLocations()
    {
        for(const auto &item: *active)
        {
            MapLocation location = static_cast<MapLocation>(item);
            if(!locations.count(location))
            {
                addLocation(location);
            }
        }

    }

    /**
     * @brief Calculates the times for each event and sorts the event list.
     */
    void createEventList()
    {
        for(auto &item: probabilities)
        {
            item.calclTimeToNextEvent();
        }
        // TODO add in events for updating the map files, extra sampling processes, or changes to the protracted
        //  speciation process (maybe we just don't support protracted speciation for now?)
        sort(probabilities.begin(), probabilities.end());
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

        probabilities.emplace_back(GillespieProbabilities(location));
        GillespieProbabilities gp = *probabilities.rbegin();
        gp.setDispersalProbability(1.0-dispersal_coordinator->getSelfDispersalProbability(cell));
        gp.setCoalescenceProbability(calculateCoalescenceProbability(location));
        gp.setSpeciationProbability(spec);
        gp.setRandomNumber(NR->d01());
        locations.insert(location);

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
