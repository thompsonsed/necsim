// This file is part of necsim project which is released under MIT license.
// See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.

/**
 * @author Samuel Thompson
 * @file Metacommunity.h
 *
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 * @brief Generates a neutral metacommunity.
 *
 * Individuals will be drawn from the metacommunity for each speciation event, instead of creating a new species each
 * time. The metacommunity itself is generated using spatially implicit neutral simulations.
 *
 * Contact: samuel.thompson14@imperial.ac.uk or thompsonsed@gmail.com
 */

#ifndef SPECIATIONCOUNTER_METACOMMUNITY_H
#define SPECIATIONCOUNTER_METACOMMUNITY_H

#include <string>
#include <sqlite3.h>
#include <set>
#include <memory>
#include "Community.h"
#include "Tree.h"
#include "RNGController.h"
#include "SpecSimParameters.h"
#include "SpeciesAbundancesHandler.h"
#include "SimulatedSpeciesAbundancesHandler.h"

using std::unique_ptr;
using std::make_shared;
using std::shared_ptr;
using std::make_unique;
namespace necsim
{
    /**
     * @brief Generates a metacommunity using spatially implicit neutral simulations, which is used to draw individuals
     * from a community.
     */
    class Metacommunity : public virtual Community
    {
    protected:
        // Simulation seed and task (read from the output database or set to 1)
        unsigned long seed;
        unsigned long task;
        bool parameters_checked;
        unique_ptr<SpeciesAbundancesHandler> species_abundances_handler;
        shared_ptr<RNGController> random;
        Tree metacommunity_tree;
    public:

        /**
         * @brief Default constructor
         */
        Metacommunity() : seed(0), task(0), parameters_checked(false),
                          species_abundances_handler(make_unique<SimulatedSpeciesAbundancesHandler>()),
                          random(make_shared<RNGController>()), metacommunity_tree()
        {
        }

        ~Metacommunity() override = default;

        Metacommunity(Metacommunity &&other) noexcept : Metacommunity()
        {
            *this = std::move(other);
        }

        Metacommunity(const Metacommunity &other) : Metacommunity()
        {
            *this = other;
        };

        Metacommunity &operator=(Metacommunity other) noexcept
        {
            other.swap(*this);
            return *this;
        }

        void swap(Metacommunity &other) noexcept
        {
            if(this != &other)
            {
                Community::swap(other);
                std::swap(seed, other.seed);
                std::swap(task, other.task);
                std::swap(parameters_checked, other.parameters_checked);
                std::swap(species_abundances_handler, other.species_abundances_handler);
                std::swap(random, other.random);
                std::swap(metacommunity_tree, other.metacommunity_tree);
            }
        }

        /**
         * @brief Sets the parameters for the metacommunity
         * @param community_size_in the number of individuals in the metacommunity
         * @param speciation_rate_in the speciation rate to use for metacommunity creation
         * @param database_name_in the path to the database to store the metacommunity in
         * @param metacommunity_option_in the metacommunity option, either "simulated", "analytical" or a path to a file
         * @param metacommunity_reference_in the metacommunity reference in the input metacommunity database
         */
        void setCommunityParameters(shared_ptr<MetacommunityParameters> metacommunity_parameters);

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
        void addSpecies(unsigned long &species_count, TreeNode* tree_node, std::set<unsigned long> &species_list) override;

        /**
         * @brief Creates the metacommunity in memory using a non-spatially_explicit neutral model, which is run using the
         * Tree class
         */
        void createMetacommunityNSENeutralModel();

        /**
         * @brief Applies the speciation parameters to the completed simulation, including running the spatially implicit
         * for the metacommunity structure, but doesn't write the output
         * @param sp speciation parameters to apply, including speciation rate, times and spatial sampling procedure.
         */
        void applyNoOutput(shared_ptr<SpecSimParameters> sp) override;

        /**
         * @brief Applies the speciation parameters to the completed simulation, including running the spatially implicit
         * for the metacommunity structure, but doesn't write the output
          * @param sp speciation parameters to apply, including speciation rate, times and spatial sampling procedure.
          */
        void applyNoOutput(shared_ptr<SpecSimParameters> sp, shared_ptr<vector<TreeNode>> tree_data) override;

        /**
         * @brief Approximates the SAD from a NSE neutral model, based on Chisholm and Pacala (2010).
         */
        void approximateSAD();

        /**
         * @brief Reads the SAD from the database provided in metacommunity_option.
         * The database must exist and have a table called SPECIES_ABUNDANCES with the relevant community parameter.
         */
        void readSAD();

        /**
         * @brief Prints the metacommunity parameters to the logger.
         */
        void printMetacommunityParameters();

    };
}
#endif //SPECIATIONCOUNTER_METACOMMUNITY_H
