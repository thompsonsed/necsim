//This file is part of necsim project which is released under MIT license.
//See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.
/**
 * @author Samuel Thompson
 * @date 30/08/2016
 * @file Datapoint.cpp
 *
 * @brief Contains the Datapoint class for storing objects during simulation run time.
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 * This class is only used during simulation runs and is not outputted to a database.
 * A Row of Datapoint objects is utilised by the main Tree objects.
 */
#include "Logging.h"
#include "DataPoint.h"
namespace necsim
{
    void DataPoint::setup(unsigned long x, unsigned long y, long xwrap_in, long ywrap_in, unsigned long reference_in,
                          unsigned long list_position_in, double min_max_in)
    {
        this->x = x;
        this->y = y;
        xwrap = xwrap_in;
        ywrap = ywrap_in;
        next_lineage = 0;
        reference = reference_in;
        list_position = list_position_in;
        nwrap = 0;
        min_max = min_max_in;
    }

    void DataPoint::setup(unsigned long reference_in, unsigned long list_position_in, double min_max_in)
    {
        setup(0, 0, 0, 0, reference_in, list_position_in, min_max_in);
    }

    void DataPoint::setup(const DataPoint &datin)
    {
        x = datin.getXpos();
        y = datin.getYpos();
        xwrap = datin.getXwrap();
        ywrap = datin.getYwrap();
        next_lineage = datin.getNext();
        //		last = datin.get_last(); // removed as of version 3.1
        reference = datin.getReference();
        list_position = datin.getListpos();
        nwrap = datin.getNwrap();
        min_max = datin.getMinmax();
    }

    void DataPoint::setReference(unsigned long z)
    {
        reference = z;
    }

    void DataPoint::setNext(unsigned long x)
    {
        next_lineage = x;
    }

    void DataPoint::setListPosition(unsigned long l)
    {
        list_position = l;
    }

    void DataPoint::setNwrap(unsigned long n)
    {
        nwrap = n;
    }

    void DataPoint::setMinmax(double d)
    {
        min_max = d;
    }

    unsigned long DataPoint::getXpos() const
    {
        return x;
    }

    unsigned long DataPoint::getYpos() const
    {
        return y;
    }

    long DataPoint::getXwrap() const
    {
        return xwrap;
    }

    long DataPoint::getYwrap() const
    {
        return ywrap;
    }

    unsigned long DataPoint::getReference() const
    {
        return reference;
    }

    unsigned long DataPoint::getNext() const
    {
        return next_lineage;
    }

    unsigned long DataPoint::getListpos() const
    {
        return list_position;
    }

    unsigned long DataPoint::getNwrap() const
    {
        return nwrap;
    }

    double DataPoint::getMinmax() const
    {
        return min_max;
    }

    void DataPoint::decreaseNwrap()
    {
        if(nwrap == 0)
        {
            throw std::out_of_range("ERROR_DATA_001: Trying to decrease  nwrap less than 0.");
        }
        else
        {
            nwrap--;
        }
    }


    std::ostream &operator<<(std::ostream &os, const DataPoint &d)
    {
        os << d.x << "," << d.y << "," << d.xwrap << "," << d.ywrap << "," << d.next_lineage << "," << d.reference
           << "," << d.list_position << "," << d.nwrap << ",";
        os << d.min_max << "\n";
        return os;
    }

    std::istream &operator>>(std::istream &is, DataPoint &d)
    {
        char delim;
        is >> d.x >> delim >> d.y >> delim >> d.xwrap >> delim >> d.ywrap >> delim >> d.next_lineage >> delim
           >> d.reference >> delim >> d.list_position >> delim >> d.nwrap >> delim;
        is >> d.min_max;
        return is;
    }

#ifdef DEBUG
    void DataPoint::logActive(const int &level)
    {
        writeLog(50, "x, y, (x wrap, y wrap): " + std::to_string(x) + ", " + std::to_string(y) + ", (" +
                     std::to_string(xwrap) + ", " + std::to_string(ywrap) + ")");
        writeLog(50, "Lineage next: " + std::to_string(next_lineage));
        writeLog(50, "Reference: " + std::to_string(reference));
        writeLog(50, "List position: " + std::to_string(list_position));
        writeLog(50, "Number in wrapped lineages: " + std::to_string(nwrap));
        writeLog(50, "Minimum maximum: " + std::to_string(min_max));
    }
#endif // DEBUG
}