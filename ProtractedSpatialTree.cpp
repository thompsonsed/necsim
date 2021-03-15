//
// Created by sam on 15/03/2021.
//
#include "ProtractedSpatialTree.h"

namespace necsim
{
    ProtractedSpatialTree::ProtractedSpatialTree() : Tree(), SpatialTree(), ProtractedTree()
    {

    }

    ProtractedSpatialTree &ProtractedSpatialTree::operator=(ProtractedSpatialTree &&other) noexcept
    {
        // Also move-assigns `Tree`
        static_cast<SpatialTree &>(*this) = std::move(static_cast<SpatialTree &>(other));
        // Copy over other member variables
        speciation_generation_min = other.speciation_generation_min;
        // The number of generations a lineage can exist before speciating.
        // All remaining lineages are speciated at this time.
        // If this value is 0, it has no effect.
        speciation_generation_max = other.speciation_generation_max;
    }
}
