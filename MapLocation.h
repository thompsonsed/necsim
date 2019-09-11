//
// Created by sam on 07/09/19.
//

#ifndef NECSIM_MAPLOCATION_H
#define NECSIM_MAPLOCATION_H

namespace necsim
{
    struct MapLocation
    {
        long x;
        long y;
        long xwrap;
        long ywrap;

        MapLocation() : x(0), y(0), xwrap(0), ywrap(0)
        { }

        MapLocation(long x, long y, long xwrap, long ywrap) : x(x), y(y), xwrap(xwrap), ywrap(ywrap)
        { }

        /**
         * @brief Checks if the location is directly on the grid without wrapping (i.e. xwrap and ywrap are 0).
         * @return true if on the grid
         */
        bool isOnGrid() const;

        /**
         * @brief Equality operator for MapLocation
         * @param m the MapLocation object to compare against
         * @return true if the x, y, xwrap and ywrap are identical
         */
        bool operator==(MapLocation const &m);

        /**
         * @brief Inequality operator for MapLocation
         * @param m the MapLocation object to compare against
         * @return true if locations are not identical
         */
        bool operator!=(MapLocation const &m);
    };
}
#endif //NECSIM_MAPLOCATION_H
