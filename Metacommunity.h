// This file is part of NECSim project which is released under BSD-3 license.
// See file **LICENSE.txt** or visit https://opensource.org/licenses/BSD-3-Clause) for full license details.

/**
 * @author Samuel Thompson
 * @file Metacommunity.cpp
 *
 * @copyright <a href="https://opensource.org/licenses/BSD-3-Clause">BSD-3 Licence.</a>
 * @brief Contains the Metacommunity class for generating a neutral metacommunity.
 *
 * For use with completed simulations from NECSim, using the SpeciationCounter program. Individuals will be drawn from
 * the metacommunity for each speciation event, instead of creating a new species each time.
 *
 * Contact: samuel.thompson14@imperial.ac.uk or thompsonsed@gmail.com
 */

#ifndef SPECIATIONCOUNTER_METACOMMUNITY_H
#define SPECIATIONCOUNTER_METACOMMUNITY_H

#include <string>
#include <sqlite3.h>
#include <set>
#include "Tree.h"
#include "Community.h"
#include "NRrand.h"
#include "SpecSimParameters.h"

using namespace std;
/**
 * @brief Contains routines for importing a metacommunity from a database, saving to a database and
 */
class Metacommunity : public virtual Community
{
protected:
	// The number of individuals in the metacommunity
	unsigned long community_size;
	// The speciation rate used for creation of the metacommunity
	long double speciation_rate;
	// Simulation seed and task (read from the output database or set to 1)
	unsigned long seed;
	unsigned long task;
	bool parameters_checked;
	Row<unsigned long> * metacommunity_cumulative_abundances;
	NRrand random;
	Tree metacommunity_tree;
public:

	Metacommunity();

	~Metacommunity() = default;

	/**
	 * @brief Sets the parameters for the metacommunity
	 * @param community_size_in the number of individuals in the metacommunity
	 * @param speciation_rate_in the speciation rate to use for metacommunity creation
	 * @param database_name_in the path to the database to store the metacommunity in
	 */
	void setCommunityParameters(unsigned long community_size_in, long double speciation_rate_in, string database_name_in);

	/**
	 * @brief Gets the seed and the task from the SIMULATION_PARAMETERS database and stores them in the relevant
	 * variables.
	 *
	 * @note Should only be called once, and will have no effect if called multiple times.
	 */
	void checkSimulationParameters();

	/**
	 * @brief Speciates TreeNode and updates the species count, including adding to the set of species ids, if a new
	 * species has been selected from the metacommunity
	 *
	 * @param species_count the total number of species currently in the community
	 * @param tree_node pointer to the TreeNode object for this lineage
	 * @param species_list the set of all species ids.
	 */
	void addSpecies(unsigned long &species_count, TreeNode *tree_node, set<unsigned long> &species_list) override;

	/**
	 * @brief Creates the metacommunity in memory using a non-spatially_explicit neutral model, which is run using the
	 * Tree class
	 */
	void createMetacommunityNSENeutralModel();

	/**
	 * @brief Selects a random lineage from the metacommunity (rOut), which should be a cumulative sum of species
	 * abundances.
	 * Performs a binary search on rOut
	 * @return the species id for the lineage
	 */
	unsigned long selectLineageFromMetacommunity();

	/**
	 * @brief Applies the speciation parameters to the completed simulation, including running the spatially-implicit
	 * for the metacommunity structure
	 * @param sp the speciation parameters to apply
	 */
	void apply(SpecSimParameters *sp) override ;

};


#endif //SPECIATIONCOUNTER_METACOMMUNITY_H
