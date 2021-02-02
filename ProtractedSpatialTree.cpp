// This file is part of necsim project which is released under MIT license.
// See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.
//
/**
 * @file ProtractedSpatialTree.cpp
 * @brief Contains the ProtractedSpatialTree class for running simulations and outputting the phylogenetic trees using
 * protracted speciation.
 *
 * Contact: samuel.thompson14@imperial.ac.uk or thompsonsed@gmail.com
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 *
 */
#include "ProtractedSpatialTree.h"

namespace necsim
{
    void ProtractedSpatialTree::setParameters()
    {
        SpatialTree::setParameters();
        setProtractedVariables(sim_parameters->max_speciation_gen, sim_parameters->max_speciation_gen);
    }
}