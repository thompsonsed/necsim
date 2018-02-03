// This file is part of NECSim project which is released under BSD-3 license.
// See file **LICENSE.txt** or visit https://opensource.org/licenses/BSD-3-Clause) for full license details.
//
/**
 * @author Sam Thompson
 * @file ProtractedTree.h
 * @brief Contains the ProtractedTree class for running simulations and outputting the phylogenetic trees using
 * protracted speciation.
 *
 * Contact: samuel.thompson14@imperial.ac.uk or thompsonsed@gmail.com
 * @copyright <a href="https://opensource.org/licenses/BSD-3-Clause">BSD-3 Licence.</a>
 *
 */

#include "SpatialTree.h"
#include "ProtractedTree.h"

/**
 * @class ProtractedTree
 * @author Sam Thompson
 * @date 10/07/2017
 * @file ProtractedTree.h
 * @brief Contains the protracted tree class, for running simulations with procated speciation.
 */
#ifndef SPECIATIONCOUNTER_PROTRACTEDSPATIALTREE_H
#define SPECIATIONCOUNTER_PROTRACTEDSPATIALTREE_H


class ProtractedSpatialTree : public SpatialTree, public ProtractedTree
{

};


#endif //SPECIATIONCOUNTER_PROTRACTEDSPATIALTREE_H
