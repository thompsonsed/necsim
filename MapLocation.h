//
// Created by sam on 07/09/19.
//

#ifndef NECSIM_MAPLOCATION_H
#define NECSIM_MAPLOCATION_H

struct MapLocation
{
    long x;
    long y;
    long xwrap;
    long ywrap;

    MapLocation():x(0), y(0), xwrap(0), ywrap(0)
    {}

    /**
     * @brief Checks if the location is directly on the grid without wrapping (i.e. xwrap and ywrap are 0).
     * @return true if on the grid
     */
    bool isOnGrid() const
    {
        return xwrap == 0 && ywrap == 0;
    }

    /**
     * @brief Equality operator for MapLocation
     * @param m the MapLocation object to compare against
     * @return true if the x, y, xwrap and ywrap are identical
     */
    bool operator==(MapLocation const &m)
    {
        return x == m.x && y == m.y && xwrap == m.xwrap && ywrap == m.ywrap;
    }

    /**
     * @brief Inequality operator for MapLocation
     * @param m the MapLocation object to compare against
     * @return true if locations are not identical
     */
    bool operator!=(MapLocation const &m)
    {
        return !(this->operator==(m));
    }
};

#endif //NECSIM_MAPLOCATION_H
