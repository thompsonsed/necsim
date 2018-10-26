// This file is part of NECSim project which is released under MIT license.
// See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.

/**
 * @author Samuel Thompson
 * @file SpeciesAbundancesHandler.cpp
 *
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 * @brief Base class for storing and generating species abundances
 *
 * Contact: samuel.thompson14@imperial.ac.uk or thompsonsed@gmail.com
 */

#include "SpeciesAbundancesHandler.h"

SpeciesAbundancesHandler::SpeciesAbundancesHandler() : species_abundances(), species_richness_per_abundance(),
													   random(make_shared<NRrand>()), max_species_id(0),
													   community_size(0), speciation_rate(0.0){}

void SpeciesAbundancesHandler::setup(shared_ptr<NRrand> random, const unsigned long &community_size,
									 const long double &speciation_rate)
{
	SpeciesAbundancesHandler::random = std::move(random);
	SpeciesAbundancesHandler::community_size = community_size;
	SpeciesAbundancesHandler::speciation_rate = speciation_rate;
}

unsigned long SpeciesAbundancesHandler::getRandomSpeciesID()
{
	unsigned long random_abundance = getRandomAbundance();
	if(species_abundances.count(random_abundance) == 0)
	{
		max_species_id++;
		species_abundances[random_abundance].push_back(max_species_id);
		species_richness_per_abundance[random_abundance] = getSpeciesRichnessOfAbundance(random_abundance);
		return max_species_id;
	}
	unsigned long random_species_index = random->i0(species_richness_per_abundance[random_abundance]-1);
	if(random_species_index >= species_abundances[random_abundance].size())
	{
		max_species_id++;
		species_abundances[random_abundance].push_back(max_species_id);
		return max_species_id;
	}
	return species_abundances[random_abundance][random_species_index];
}

void SpeciesAbundancesHandler::setAbundanceList(const shared_ptr<map<unsigned long, unsigned long>> &abundance_list_in)
{

}

void SpeciesAbundancesHandler::setAbundanceList(shared_ptr<vector<unsigned long>> abundance_list_in)
{

}
