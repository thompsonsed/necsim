//
// Created by sam on 07/09/19.
//


#ifndef NECSIM_GILLESPIECALCULATOR_H
#define NECSIM_GILLESPIECALCULATOR_H

#include <map>
#include <algorithm>
#include <numeric>
#include <memory>

#include "heap.hpp"

#include "Cell.h"
#include "RNGController.h"
#include "MapLocation.h"

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

        GillespieProbability() : GillespieProbability(MapLocation()){}

        GillespieProbability(const MapLocation &c) : dispersal_outside_cell_probability(0.0),
                                                     coalescence_probability(0.0), speciation_probability(0.0),
                                                     random_number(0.0), location(c)
        {

        }

        void setDispersalOutsideCellProbability(const double &d);

        void setCoalescenceProbability(const double &c);

        void setSpeciationProbability(const double &s);

        void setRandomNumber(const double &r);

        double getInCellProbability() const;

        CellEventType generateRandomEvent(const shared_ptr<RNGController> &rng) const;

        MapLocation &getMapLocation();

        /**
         * @brief Gets the parameter for the exponential function
         *
         * The rate is per birth-death event on the whole landscape.
         *
         * @param summed_death_rate
         * @param n
         * @param total_n
         * @return lambda
         */
        double getLambda(const double &local_death_rate, const double &summed_death_rate, const unsigned long &n) const;

        double calcTimeToNextEvent(const double &local_death_rate,
                                   const double &mean_death_rate,
                                   const unsigned long &n) const;

        void reset();

        friend ostream &operator<<(ostream &os, const GillespieProbability &gp);

        friend std::istream &operator>>(std::istream &is, GillespieProbability &gp);

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
        unsigned long *pos;
        EventType event_type;
        vector<GillespieHeapNode> *heap;

        GillespieHeapNode(const Cell cell,
                          const double time_of_event,
                          unsigned long *pos,
                          const EventType &e,
                          vector<GillespieHeapNode> *heap) : cell(cell), time_of_event(time_of_event), pos(pos),
                                                             event_type(e), heap(heap){}

        GillespieHeapNode() : GillespieHeapNode(Cell(), 0.0, nullptr, EventType::undefined, nullptr){}

        GillespieHeapNode(const double time_of_event, const EventType &e) : GillespieHeapNode(Cell(), time_of_event,
                                                                                              nullptr, e, nullptr){}

        GillespieHeapNode(GillespieHeapNode &&other) noexcept
        {
            cell = other.cell;
            time_of_event = other.time_of_event;
            event_type = other.event_type;
            pos = other.pos;
            heap = other.heap;
            
            if (this->pos != nullptr && this->heap != nullptr)
            {
                if (this >= heap->data() && this < (heap->data() + heap->size()))
                {
                    /*stringstream ss;
                    ss << "Custom move constructor to " << (this - heap->data()) << endl;
                    necsim::writeInfo(ss.str());*/
                    
                    *pos = this - heap->data();
                }
            }
        }

        GillespieHeapNode &operator=(GillespieHeapNode &&other)
        {
            cell = other.cell;
            time_of_event = other.time_of_event;
            event_type = other.event_type;
            pos = other.pos;
            heap = other.heap;
            
            if (this->pos != nullptr && this->heap != nullptr)
            {
                if (this >= heap->data() && this < (heap->data() + heap->size()))
                {
                    /*stringstream ss;
                    ss << "Custom move assignment to " << (this - heap->data()) << endl;
                    necsim::writeInfo(ss.str());*/
                    
                    *pos = this - heap->data();
                }
            }

            return *this;
        }

        GillespieHeapNode(const GillespieHeapNode &other) = default;

        bool operator<(const GillespieHeapNode &n) const
        {
            return time_of_event > n.time_of_event;
        }
    };

}

#endif //NECSIM_GILLESPIECALCULATOR_H
