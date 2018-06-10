// This file is part of NECSim project which is released under MIT license.
// See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.

/**
 * @author Samuel Thompson
 * @file Metacommunity.cpp
 *
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 * @brief Contains the Metacommunity class for generating a neutral metacommunity.
 *
 * For use with completed simulations from NECSim, using the SpeciationCounter program. Individuals will be drawn from
 * the metacommunity for each speciation event, instead of creating a new species each time.
 *
 * Contact: samuel.thompson14@imperial.ac.uk or thompsonsed@gmail.com
 */

#include "Metacommunity.h"
#include "LogFile.h"

Metacommunity::Metacommunity(): community_size(0), speciation_rate(0.0), seed(0), task(0), parameters_checked(false),
								metacommunity_cumulative_abundances(nullptr), random(), metacommunity_tree()
{

}

void Metacommunity::setCommunityParameters(unsigned long community_size_in, long double speciation_rate_in)
{
	community_size = community_size_in;
	speciation_rate = speciation_rate_in;
}

void Metacommunity::checkSimulationParameters()
{
	if(!parameters_checked)
	{
		if(database == nullptr)
		{
			throw FatalException("Cannot read simulation parameters as database is null pointer.");
		}
		// Now do the same for times
		sqlite3_stmt *stmt = nullptr;
		string sql_call = "SELECT seed, task from SIMULATION_PARAMETERS";
		int rc = sqlite3_prepare_v2(database, sql_call.c_str(), static_cast<int>(strlen(sql_call.c_str())), &stmt,
									nullptr);
		if(rc != SQLITE_DONE && rc != SQLITE_OK)
		{
			sqlite3_close(database);
		}
		sqlite3_step(stmt);
		seed = static_cast<unsigned long>(sqlite3_column_int(stmt, 1));
		random.setSeed(seed);
		task = static_cast<unsigned long>(sqlite3_column_int(stmt, 1));
		sqlite3_step(stmt);
		sqlite3_finalize(stmt);
		parameters_checked = true;
	}
}

void Metacommunity::addSpecies(unsigned long &species_count, TreeNode *tree_node, set<unsigned long> &species_list)
{

	auto species_id = selectLineageFromMetacommunity();
	if(species_list.empty() || species_list.find(species_id) != species_list.end())
	{
		species_list.insert(species_id);
		species_count ++;
	}
	tree_node->burnSpecies(species_id);
}

void Metacommunity::createMetacommunityNSENeutralModel()
{
#ifdef DEBUG
	writeLog(10, "Running spatially-implicit model for metacommunity generation.");
#endif //DEBUG
	// First set up a non-spatial coalescence simulation to generate our metacommunity
	SimParameters temp_parameters;
	temp_parameters.setMetacommunityParameters(community_size, speciation_rate, seed, task);
	metacommunity_tree.internalSetup(temp_parameters);
	// Run our simulation and calculate the species abundance distribution (as this is all that needs to be stored).
	if(!metacommunity_tree.runSimulation())
	{
		throw FatalException("Completion of the non-spatial coalescence simulation "
									 "to create the metacommunity did not finish in time.");
	}
	metacommunity_tree.applySpecRateInternal(speciation_rate, 0.0);
	// row_out now contains the number of individuals per species
	// Make it cumulative to increase the speed of indexing using binary search.
	metacommunity_cumulative_abundances = metacommunity_tree.getCumulativeAbundances();
#ifdef DEBUG
	writeLog(10, "Spatially-implicit simulation completed.");
#endif //DEBUG

}

unsigned long Metacommunity::selectLineageFromMetacommunity()
{
	auto max_indices = metacommunity_cumulative_abundances->size() - 1;
	auto random_value = random.i0(community_size - 1);
#ifdef DEBUG
	// binary search
	if(random_value > (*metacommunity_cumulative_abundances)[max_indices])
	{
		throw FatalException("Random number generation out of range of the community size in lineage selection.");
	}
#endif //DEBUG
	unsigned long mid_point;
	unsigned long min_indices = 0;
	while(min_indices < max_indices-1)
	{
		mid_point = static_cast<unsigned long>(floor(((max_indices - min_indices) / 2) + min_indices));
		if(random_value == (*metacommunity_cumulative_abundances)[mid_point])
		{
			min_indices = mid_point;
			max_indices = mid_point;
		}
		if(random_value <= (*metacommunity_cumulative_abundances)[mid_point])
		{
			max_indices = mid_point;
		}
		else
		{
			min_indices = mid_point;
		}
	}
	if(min_indices == max_indices - 1)
	{
		return max_indices;
	}
#ifdef DEBUG
	if(min_indices != max_indices)
	{
		throw FatalException("Error in binary search algorithm for lineage selection. Please report this bug.");
	}
#endif // DEBUG
	return min_indices;
}

void Metacommunity::applyNoOutput(SpecSimParameters *sp)
{
#ifdef DEBUG
	writeLog(10, "********************");
	writeLog(10, "Metacommunity application");
#endif //DEBUG
	setCommunityParameters(sp->metacommunity_size, sp->metacommunity_speciation_rate);
	// Make sure that the connection is opened to file.
	openSqlConnection(sp->filename);
	checkSimulationParameters();
	closeSqlConnection();
	createMetacommunityNSENeutralModel();
#ifdef DEBUG
	writeLog(10, "Creating coalescence tree from metacommunity...");
#endif //DEBUG
	Community::applyNoOutput(sp);
}





