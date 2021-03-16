//
// Created by sam on 15/03/2021.
//
#include "ProtractedSpatialTree.h"

namespace necsim
{
    ProtractedSpatialTree::ProtractedSpatialTree() : Tree(), SpatialTree(), ProtractedTree()
    {

    }

    ProtractedSpatialTree::ProtractedSpatialTree(const ProtractedSpatialTree &other)
     : SpatialTree(other) {
        speciation_generation_min = other.speciation_generation_min;
        speciation_generation_max = other.speciation_generation_max;
    }

    ProtractedSpatialTree &ProtractedSpatialTree::operator=(const ProtractedSpatialTree &other) noexcept
    {
        static_cast<SpatialTree &>(*this) = static_cast<const SpatialTree &>(other);
        speciation_generation_min = other.speciation_generation_min;
        speciation_generation_max = other.speciation_generation_max;
        return *this;
    }

    ProtractedSpatialTree &ProtractedSpatialTree::operator=(ProtractedSpatialTree &&other) noexcept
    {
        // Also move-assigns `Tree`
        static_cast<SpatialTree &>(*this) = std::move(static_cast<SpatialTree &>(other));
        // Copy over other member variables
        speciation_generation_min = other.speciation_generation_min;
        speciation_generation_max = other.speciation_generation_max;
        return *this;
    }
}
