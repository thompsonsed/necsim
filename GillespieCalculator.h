//
// Created by sam on 07/09/19.
//

#ifndef NECSIM_GILLESPIECALCULATOR_H
#define NECSIM_GILLESPIECALCULATOR_H

#include <map>
#include <algorithm>
#include <numeric>

#include "heap.hpp"

#include "ActivityMap.h"
#include "Cell.h"
#include "DispersalCoordinator.h"
#include "RNGController.h"
#include "DataPoint.h"
#include "SpatialTree.h"

namespace necsim
{
    using namespace random_numbers;

    /**
     * @brief Container for the different event types that can occur during the Gillespie Algorithm.
     */
    enum class EventType
    {
        undefined, cell_event, map_event, sample_event
    };

    enum class CellEventType
    {
        undefined, dispersal_event, coalescence_event, speciation_event
    };

    class GillespieProbability
    {
    protected:
        double dispersal_outside_cell_probability;
        double coalescence_probability;
        double speciation_probability;
        double cell_turnover_probability;
        double random_number;
        MapLocation location;

    public:

        GillespieProbability() : GillespieProbability(MapLocation())
        { }

        GillespieProbability(const MapLocation &c) : dispersal_outside_cell_probability(0.0),
                                                     coalescence_probability(0.0), speciation_probability(0.0),
                                                     random_number(0.0), location(c)
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
            return speciation_probability + (1 - speciation_probability)
                                            * ((1 - dispersal_outside_cell_probability) * coalescence_probability
                                               + dispersal_outside_cell_probability);
        }

        CellEventType generateRandomEvent(const shared_ptr<RNGController> &rng) const
        {
            double p = rng->d01() * getInCellProbability();
            if(p < speciation_probability)
            {
                return CellEventType::speciation_event;
            }
            else
            {
                if(p < speciation_probability + (1 - speciation_probability) * dispersal_outside_cell_probability)
                {
                    return CellEventType::dispersal_event;
                }
                else
                {
                    return CellEventType::coalescence_event;
                }
            }
        }

        MapLocation &getMapLocation()
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

        double calcTimeToNextEvent(const double &mean_death_rate, const unsigned long &n, const unsigned long &total_n)
        {
            return RNGController::exponentialDistribution(getLambda(mean_death_rate, n, total_n), random_number);
        }

        void reset()
        {

            dispersal_outside_cell_probability = 0.0;
            coalescence_probability = 0.0;
            speciation_probability = 0.0;
            cell_turnover_probability = 0.0;
            random_number = 0.0;
        }

        friend ostream &operator<<(ostream &os, const GillespieProbability &gp)
        {
            os << gp.random_number << "," << gp.speciation_probability << "," << gp.coalescence_probability << ","
               << gp.cell_turnover_probability << "," << gp.dispersal_outside_cell_probability << "," << gp.location
               << endl;
            return os;
        }

        friend std::istream &operator>>(std::istream &is, GillespieProbability &gp)
        {
            char delim;
            is >> gp.random_number >> delim >> gp.speciation_probability >> delim >> gp.coalescence_probability >> delim
               >> gp.cell_turnover_probability >> delim >> gp.dispersal_outside_cell_probability >> delim
               >> gp.location;
            return is;
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
        double time_of_event;
        // Pointer to index in heap
        unsigned long* pos;
        EventType event_type;

        GillespieHeapNode(Cell cell, double time_of_event, unsigned long* pos, const EventType &e) : cell(cell),
                                                                                                     time_of_event(
                                                                                                             time_of_event),
                                                                                                     pos(pos),
                                                                                                     event_type(e)
        { }

        GillespieHeapNode() : GillespieHeapNode(Cell(), 0.0, nullptr, EventType::undefined)
        { }

        GillespieHeapNode(const double time_of_event, const EventType &e) : GillespieHeapNode(Cell(),
                                                                                              time_of_event,
                                                                                              nullptr,
                                                                                              e)
        { }
    };

    bool operator<(const GillespieHeapNode &n1, const GillespieHeapNode &n2)
    {
        return n1.time_of_event > n2.time_of_event;
    }

    class GillespieCalculator : public SpatialTree
    {
    protected:
        // Matrix of all the probabilities at every location in the map.
        Matrix<GillespieProbability> probabilities;
        // Vector used for holding the priority queue as a binary heap
        vector<GillespieHeapNode> heap;
        // Index to heap position, or UNUSED if cell is not used.
        Matrix<unsigned long> cellToHeapPositions;

        unsigned long global_individuals;
        double mean_death_rate;

        static const unsigned long UNUSED = static_cast<long>(-1);

    public:
        GillespieCalculator() : probabilities()
        { }

        bool runSimulationGillespie()
        {
            do
            {
                runSingleLoop();
            }
            while((endactive < gillespie_threshold) && (endactive > 1) && (steps < 100)
                  || difftime(sim_end, start) < maxtime && this_step.bContinueSim);
            // Switch to gillespie
            setupGillespie();
            do
            {
                runGillespieLoop();
            }
            while(endactive > 1);
            return stopSimulation();

        }

        void runGillespieLoop()
        {

            // Decide what event and execute
            EventType next_event = heap.front().event_type;
            // Update the event timer
            steps += (heap.front().time_of_event - generation) * double(endactive);
            generation = heap.front().time_of_event;

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

        /**
         * @brief Sets up the Gillespie algorithm.
         */
        void setupGillespie()
        {
            findLocations();
            updateAllProbabilities();
            createEventList();
            checkMapEvents();
            checkSampleEvents();
            sortEvents();
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

        void checkMapEvents()
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

        void checkSampleEvents()
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

        void gillespieCellEvent(GillespieProbability &origin)
        {
            CellEventType cell_event = origin.generateRandomEvent(NR);
            origin.setRandomNumber(NR->d01());
            switch(cell_event)
            {
            case CellEventType::coalescence_event:
                // implement coalescence
                gillespieCoalescenceEvent(origin);
                break;

            case CellEventType::dispersal_event:
                // choose dispersal
                gillespieDispersalEvent(origin);
                break;

            case CellEventType::speciation_event:
                gillespieSpeciationEvent(origin);
                break;

            case CellEventType::undefined:
                throw FatalException("Undefined cell event type. Please report this bug.");
                break;
            }

        }

        void gillespieUpdateMap()
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

        void gillespieSampleIndividuals()
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

        void gillespieCoalescenceEvent(GillespieProbability &origin)
        {
            auto lineages = selectTwoRandomLineages(origin.getMapLocation());
            gillespieUpdateGeneration(lineages.first);
            coalescenceEvent(lineages.first, lineages.second);
            const MapLocation &location = origin.getMapLocation();
            updateCellCoalescenceProbability(origin, getNumberIndividualsAtLocation(location));
            updateInhabitedCellOnHeap(convertMapLocationToCell(location));

        }

        void gillespieDispersalEvent(GillespieProbability &origin)
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
            unsigned long x = active[chosen].getXpos();
            unsigned long y = active[chosen].getYpos();
            Cell destination_cell(x, y);
            unsigned long n = getNumberIndividualsAtLocation(origin.getMapLocation());
            if(n > 0)
            {
                updateCellCoalescenceProbability(origin, n);
                updateInhabitedCellOnHeap(destination_cell);
            }
            else
            {
                removeHeapTop();
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
                const double t = destination.calcTimeToNextEvent(mean_death_rate,
                                                                 getNumberIndividualsAtLocation(destination
                                                                                                        .getMapLocation()),
                                                                 global_individuals);
                heap[cellToHeapPositions.get(y, x)].time_of_event = t;
                updateInhabitedCellOnHeap(destination_cell);
            }

        }

        void gillespieSpeciationEvent(GillespieProbability &origin)
        {
            const MapLocation &location = origin.getMapLocation();
            unsigned long chosen = selectRandomLineage(location);
            gillespieUpdateGeneration(chosen);
            speciateLineage(active[chosen].getReference());
            updateCellCoalescenceProbability(origin, getNumberIndividualsAtLocation(location));
            updateInhabitedCellOnHeap(convertMapLocationToCell(origin.getMapLocation()));
        }

        void clearGillespieObjects()
        {
            cellToHeapPositions.fill(0);
            heap.clear();
            for(auto &item : probabilities)
            {
                item.reset();
            }
        }

        //        void updateHeapOrigin(GillespieProbability &origin)
        //        {
        //            // If cell is still inhabited
        //            {
        //                Cell pos = convertMapLocationToCell(origin.getMapLocation());
        //
        //                // Update probability and times of origin
        //                updateInhabitedCellOnHeap(pos);
        //
        //            }
        //            // Else
        //            {
        //                std::pop_heap(heap.begin(), heap.end());
        //                heap.pop_back();
        //            }
        //        }

        void gillespieUpdateGeneration(const unsigned long &lineage)
        {
            TreeNode &tree_node = (*data)[active[lineage].getReference()];
            unsigned long generations_existed = round(generation) - tree_node.getGeneration();
            tree_node.setGeneration(generations_existed);
        }

        //        void updateHeapDestination(GillespieProbability &destination)
        //        {
        //            // If dispersal to other cell
        //
        //            Cell pos = convertMapLocationToCell(destination.getMapLocation());
        //
        //            // If cell is already inhabited
        //            {
        //                // Update probability and times of destination
        //                updateInhabitedCellOnHeap(pos);
        //            }
        //            // Else
        //            {
        //                heap.emplace_back(GillespieHeapNode(pos,
        //                                                    (destination.calcTimeToNextEvent(mean_death_rate,
        //                                                                                     getNumberIndividualsAtLocation(
        //                                                                                             destination
        //                                                                                                     .getMapLocation()),
        //                                                                                     global_individuals) + generation),
        //                                                    &cellToHeapPositions.get(pos.y, pos.x)));
        //                std::push_heap(heap.begin(), heap.end());
        //            }
        //        }

        void updateCellCoalescenceProbability(GillespieProbability &origin, const unsigned long &n)
        {
            origin.setCoalescenceProbability(calculateCoalescenceProbability(origin.getMapLocation()));
            origin.setRandomNumber(NR->d01());
            heap.front().time_of_event = origin.calcTimeToNextEvent(mean_death_rate, n, global_individuals);
        }

        void updateInhabitedCellOnHeap(const Cell &pos)
        {
            std::update_heap(heap.begin(), heap.end(), heap.begin() + cellToHeapPositions.get(pos.y, pos.x));
        }

        void updateAllProbabilities()
        {
            // calculate global death rate mean
            mean_death_rate = death_map->getMean();
            // calculate total number of individuals
            global_individuals = std::accumulate(landscape->getFineMap().begin(),
                                                 landscape->getFineMap().end(),
                                                 (unsigned long) 0);
        }

        void removeHeapTop()
        {
            std::pop_heap(heap.begin(), heap.end());
            heap.pop_back();
        }

        Cell convertMapLocationToCell(const MapLocation &location)
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
        }

        void sortEvents()
        {
            std::make_heap(heap.begin(), heap.end());
        }

        void addNewEvent(const unsigned long &x, const unsigned long &y)
        {
            unsigned long local_individuals = getNumberIndividualsAtLocation(probabilities.get(y, x).getMapLocation());
            if(local_individuals > 0)
            {
                cellToHeapPositions.get(y, x) = heap.size();
                heap.emplace_back(GillespieHeapNode(Cell(x, y),
                                                    (probabilities.get(y, x).calcTimeToNextEvent(mean_death_rate,
                                                                                                 local_individuals,
                                                                                                 global_individuals)
                                                     + generation),
                                                    &cellToHeapPositions.get(y, x),
                                                    EventType::cell_event));
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

            GillespieProbability gp(location);
            // check if any lineages exist there

            gp.setDispersalOutsideCellProbability(1.0 - dispersal_coordinator.getSelfDispersalProbability(cell));
            gp.setCoalescenceProbability(calculateCoalescenceProbability(location));
            gp.setSpeciationProbability(spec);
            gp.setRandomNumber(NR->d01());

            probabilities.get(cell.y, cell.x) = gp;
        }

        double calculateCoalescenceProbability(const MapLocation &location)
        {
            unsigned long max_number_individuals = landscape
                    ->getVal(location.x, location.y, location.xwrap, location.ywrap, generation);
            unsigned long current_number = getNumberIndividualsAtLocation(location);
            return min(double(current_number) / double(max_number_individuals), 1.0);

        }

        unsigned long selectRandomLineage(const MapLocation &location) const
        {
            vector<unsigned long> lineage_ids = detectLineages(location);
            unsigned long random_index = NR->i0(lineage_ids.size());
            return lineage_ids[random_index];
        }

        pair<unsigned long, unsigned long> selectTwoRandomLineages(const MapLocation &location) const
        {
            vector<unsigned long> lineage_ids = detectLineages(location);
            if(lineage_ids.size() < 2)
            {
                throw FatalException("Cannot select two lineages when fewer than two exist at location.");
            }
            pair<unsigned long, unsigned long> selected_lineages;
            selected_lineages.first = lineage_ids[NR->i0(lineage_ids.size())];

            do
            {
                selected_lineages.second = lineage_ids[NR->i0(lineage_ids.size())];
            }
            while(selected_lineages.second == selected_lineages.first);
            return selected_lineages;
        }

        vector<unsigned long> detectLineages(const MapLocation &location) const
        {
            SpeciesList species_list = grid.getCopy(location.y, location.x);
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
            return lineage_ids;

        }

    };
}
namespace std
{
    template<>
    void swap<necsim::GillespieHeapNode>(necsim::GillespieHeapNode &lhs, necsim::GillespieHeapNode &rhs) noexcept
    {
        std::swap(lhs.cell, rhs.cell);
        std::swap(lhs.time_of_event, rhs.time_of_event);
        std::swap(lhs.event_type, rhs.event_type);
        std::swap(*lhs.pos, *rhs.pos);
    }
}

#endif //NECSIM_GILLESPIECALCULATOR_H
