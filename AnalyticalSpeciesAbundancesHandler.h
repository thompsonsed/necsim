//
// Created by Sam Thompson on 16/10/2018.
//

#ifndef ANALYICAL_SPECIES_ABUNDANCES_H
#define ANALYICAL_SPECIES_ABUNDANCES_H

#include "SpeciesAbundancesHandler.h"
#include "neutral_analytical.h"
#include "NRrand.h"

namespace na = neutral_analytical;
using namespace std;

class AnalyticalSpeciesAbundancesHandler : public virtual SpeciesAbundancesHandler
{
public:



	/**
	 * @brief Gets a random species abundance by sampling from the logarithmic distribution.
	 * @return the randomly generated abundance
	 */
	unsigned long getRandomAbundance() override
	{
		return static_cast<unsigned long>(max(static_cast<unsigned long>(
													  min(random->randomLogarithmic(1.0 - speciation_rate),
														  community_size)), (unsigned long) 1));
	}


	/**
	 * @brief Gets the species richness of a particular abundance class.
	 *
	 * This method is slightly incorrect, as it does not account for the variance in the species abundance, but
	 * will approximate the mean behaviour.
	 *
	 * @param abundance the abundance class of the species
	 * @return the number of species with that abundance
	 */
	unsigned long getSpeciesRichnessOfAbundance(const unsigned long &abundance) override
	{
		auto random_number_species = na::nseMetacommunitySpeciesWithAbundance(abundance, community_size,
																			  speciation_rate);
		auto min_species_number = static_cast<unsigned long>(floor(random_number_species));
		if(random->d01() < random_number_species - min_species_number)
		{
			min_species_number ++;
		}
		return min_species_number;

	}
};

#endif //ANALYICAL_SPECIES_ABUNDANCES_H
