//
// Created by sam on 15/03/2021.
//
#include "ProtractedSpatialTree.h"

namespace necsim
{
    ProtractedSpatialTree::ProtractedSpatialTree() : SpatialTree(), ProtractedTree()
    {
        std::cout << "Creating protracted spatial tree..." << std::endl; // TODO remove
    }

    ProtractedSpatialTree::ProtractedSpatialTree(ProtractedSpatialTree &&other) noexcept
    {
        *this = std::move(other);
    }

    ProtractedSpatialTree::ProtractedSpatialTree(const ProtractedSpatialTree &other) noexcept : Tree(other), SpatialTree(other),
                                                                                       ProtractedTree(other)
    {
        std::cout << "Copy constructor protracted spatial tree..." << std::endl; // TODO remove
    }

    ProtractedSpatialTree &ProtractedSpatialTree::operator=(const ProtractedSpatialTree &other) noexcept
    {
        std::cout << "Copy operator spatial tree..." << std::endl; // TODO remove
        static_cast<SpatialTree &>(*this) = static_cast<const SpatialTree &>(other);
        speciation_generation_min = other.speciation_generation_min;
        speciation_generation_max = other.speciation_generation_max;
        return *this;
    }

    ProtractedSpatialTree &ProtractedSpatialTree::operator=(ProtractedSpatialTree &&other) noexcept
    {
        std::cout << "Move operator protracted spatial tree..." << std::endl; // TODO remove
        // Also move-assigns `Tree`
        static_cast<SpatialTree &>(*this) = std::move(static_cast<SpatialTree &&>(other));
        // Copy over other member variables
        speciation_generation_min = other.speciation_generation_min;
        speciation_generation_max = other.speciation_generation_max;
        return *this;
    }
}
