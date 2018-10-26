// This file is part of NECSim project which is released under MIT license.
// See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.

/**
 * @author Samuel Thompson
 * @file SimulatedSpeciesAbundances.cpp
 *
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 * @brief Class for repeatedly selecting random species from a distribution of species abundances.
 *
 * Contact: samuel.thompson14@imperial.ac.uk or thompsonsed@gmail.com
 */

#include "SimulatedSpeciesAbundancesHandler.h"

SimulatedSpeciesAbundancesHandler::SimulatedSpeciesAbundancesHandler() : abundance_list(
		make_shared<vector<unsigned long>>()), total_species_number(0), number_of_individuals(0){}

void SimulatedSpeciesAbundancesHandler::setAbundanceList(
		const shared_ptr<map<unsigned long, unsigned long>> &abundance_list_in)
{
	max_species_id = 0;
	number_of_individuals = 0;
	for(const auto &item : *abundance_list_in)
	{
		species_richness_per_abundance[item.second] ++;
		species_abundances[item.second].push_back(item.first);
		abundance_list->push_back(item.second);
		max_species_id = max(max_species_id, item.first);
	}
}

void SimulatedSpeciesAbundancesHandler::setAbundanceList(shared_ptr<vector<unsigned long>> abundance_list_in)
{
	abundance_list = abundance_list_in;
	max_species_id = 0;
	number_of_individuals = 0;
	generateAbundanceTable();
}

void SimulatedSpeciesAbundancesHandler::generateAbundanceTable()
{
	max_species_id = 0;
	for(const auto &item: (*abundance_list))
	{
		max_species_id++;
		species_abundances[item].push_back(max_species_id);
	}
	for(const auto &item: species_abundances)
	{
		species_richness_per_abundance[item.first] = item.second.size();
	}
}

unsigned long SimulatedSpeciesAbundancesHandler::getRandomAbundance()
{
	unsigned long result = 0;
	if(abundance_list->size() == 1)
	{
		if((*abundance_list)[0] == 0)
		{
			throw FatalException("Only one abundance found in abundance list, and it is 0. Please report this bug.");
		}
		return (*abundance_list)[0];
	}
	do
	{
		result = (*abundance_list)[random->i0(abundance_list->size() - 1)];
	}
	while(result == 0);
	return result;
}

unsigned long SimulatedSpeciesAbundancesHandler::getSpeciesRichnessOfAbundance(const unsigned long &abundance)
{
	if(species_richness_per_abundance.count(abundance) == 0)
	{
		stringstream ss;
		ss << "Cannot find species of abundance " << abundance << " - please report this bug." << endl;
		throw FatalException(ss.str());
	}
	return species_richness_per_abundance[abundance];
}
