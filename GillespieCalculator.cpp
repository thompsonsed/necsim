//
// Created by sam on 07/09/19.
//

#include "GillespieCalculator.h"

namespace necsim
{

    void GillespieProbability::setDispersalOutsideCellProbability(const double &d)
    {
        this->dispersal_outside_cell_probability = d;
    }

    void GillespieProbability::setCoalescenceProbability(const double &c)
    {
        this->coalescence_probability = c;
    }

    void GillespieProbability::setSpeciationProbability(const double &s)
    {
        this->speciation_probability = s;
    }

    void GillespieProbability::setRandomNumber(const double &r)
    {
        random_number = r;
    }

    double GillespieProbability::getInCellProbability() const
    {
        return speciation_probability + (1.0 - speciation_probability)
                                        * ((1.0 - dispersal_outside_cell_probability) * coalescence_probability
                                           + dispersal_outside_cell_probability);
    }

    CellEventType GillespieProbability::generateRandomEvent(const shared_ptr<RNGController> &rng) const
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
//                stringstream ss; // TODO remove
//                ss << "Selecting coalescence probability with change " << coalescence_probability << "/"
//                   << getInCellProbability() << endl;
//                writeInfo(ss.str());
                return CellEventType::coalescence_event;
            }
        }
    }

    MapLocation &GillespieProbability::getMapLocation()
    {
        return location;
    }

    double GillespieProbability::getLambda(const double &local_death_rate, const double &summed_death_rate,
                                           const unsigned long &n) const
    {
        return getInCellProbability() * local_death_rate * double(n) / summed_death_rate;
    }

    double GillespieProbability::calcTimeToNextEvent(const double &local_death_rate, const double &summed_death_rate,
                                                     const unsigned long &n) const
    {
        //        stringstream ss;
        //        ss << "Generating event with lambda " << getLambda(local_death_rate, summed_death_rate, n)
        //           << " and random number " << random_number << endl;
        //        ss << "Rates: " << local_death_rate << ", "  << summed_death_rate << ", " << n << endl;
        //        writeInfo(ss.str()); // TODO remove
        return RNGController::exponentialDistribution(getLambda(local_death_rate, summed_death_rate, n), random_number);
    }

    void GillespieProbability::reset()
    {

        dispersal_outside_cell_probability = 0.0;
        coalescence_probability = 0.0;
        speciation_probability = 0.0;
        cell_turnover_probability = 0.0;
        random_number = 0.0;
    }

    ostream &operator<<(ostream &os, const GillespieProbability &gp)
    {
        os << gp.random_number << "," << gp.speciation_probability << "," << gp.coalescence_probability << ","
           << gp.cell_turnover_probability << "," << gp.dispersal_outside_cell_probability << "," << gp.location
           << endl;
        return os;
    }

    std::istream &operator>>(std::istream &is, GillespieProbability &gp)
    {
        char delim;
        is >> gp.random_number >> delim >> gp.speciation_probability >> delim >> gp.coalescence_probability >> delim
           >> gp.cell_turnover_probability >> delim >> gp.dispersal_outside_cell_probability >> delim >> gp.location;
        return is;
    }

    /*void swap(GillespieHeapNode &lhs, GillespieHeapNode &rhs) noexcept
    {
        writeInfo("Custom swap\n");
        using std::swap;
        swap(lhs.cell, rhs.cell);
        swap(lhs.time_of_event, rhs.time_of_event);
        swap(lhs.event_type, rhs.event_type);
        swap(*lhs.pos, *rhs.pos);
    }*/

}