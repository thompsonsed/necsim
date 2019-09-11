//
// Created by sam on 07/09/19.
//

#include "MapLocation.h"
namespace necsim
{
    bool MapLocation::isOnGrid() const
    {
        return xwrap == 0 && ywrap == 0;
    }

    bool MapLocation::operator==(MapLocation const &m)
    {
        return x == m.x && y == m.y && xwrap == m.xwrap && ywrap == m.ywrap;
    }

    bool MapLocation::operator!=(MapLocation const &m)
    {
        return !(this->operator==(m));
    }
}
