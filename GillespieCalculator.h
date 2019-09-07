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


public:

    GillespieProbabilities(MapLocation c) : dispersal_probability(0.0), coalescence_probability(0.0),
                                            speciation_probability(0.0),
                                            random_number(0.0), location(c), time_to_event(0.0)
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

    void setRandomNumber(const long double &randomNumber)
    {
        random_number = randomNumber;
    }

    double getTotalProbability() const
    {
        return cell_turnover_probability * (dispersal_probability + coalescence_probability + speciation_probability);
    }

    void setTimeToNextEvent()
    {
        time_to_event = RNGController::exponentialDistribution(getTotalProbability(), random_number);
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

    void runGillespie()
    {
        setupGillespie();
        do
        {
            runSingleLoop();
        }
        while(endactive > 1);


    }

    void setupGillespie()
    {
        findLocations();
        createEventList();
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

    void createEventList()
    {
        for(auto &item: probabilities)
        {
            item.setTimeToNextEvent();
        }
        // TODO add in events for updating the map files, extra sampling processes, or changes to the protracted
        //  speciation process (maybe we just don't support protracted speciation for now?)
        sort(probabilities.begin(), probabilities.end());

    }

    void addLocation(const MapLocation &location)
    {
        Cell cell = getCellOfMapLocation(location);

        probabilities.emplace_back(GillespieProbabilities(location));
        GillespieProbabilities gp = *probabilities.rbegin();
        gp.setDispersalProbability(dispersal_coordinator->getSelfDispersalProbability(cell));
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
