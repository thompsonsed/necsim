// This file is part of necsim project which is released under MIT license.
// See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.

/**
 * @author Samuel Thompson
 * @file AnalyticalSpeciesAbundancesHandler.cpp
 *
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 * @brief Class for repeatedly selecting random species from a distribution of species abundances using analytical
 * solutions from Vallade and Houchmandzadeh (2003) and Alonso and McKane (2004).
 *
 * Contact: samuel.thompson14@imperial.ac.uk or thompsonsed@gmail.com
 */

#include "AnalyticalSpeciesAbundancesHandler.h"
#include "custom_exceptions.h"

AnalyticalSpeciesAbundancesHandler::AnalyticalSpeciesAbundancesHandler() : seen_no_individuals(0), ind_to_species()
{

}

void
AnalyticalSpeciesAbundancesHandler::setup(shared_ptr<RNGController> random, const unsigned long &metacommunity_size,
                                          const long double &speciation_rate,
                                          const unsigned long &local_community_size)
{
    SpeciesAbundancesHandler::setup(random, metacommunity_size, speciation_rate, local_community_size);
    generateSpeciesAbundances();
}

void AnalyticalSpeciesAbundancesHandler::generateSpeciesAbundances()
{
    writeInfo("burning in species abundance..."); // TODO replace or remove
    auto expected_richness = static_cast<unsigned long>(neutral_analytical::nseSpeciesRichness(metacommunity_size, speciation_rate));
    for(unsigned long i = 0; i < expected_richness; i ++)
    {
        addNewSpecies();
    }
//    for(unsigned long i = 0; i < local_community_size; i++)
//    {
//        addNewSpecies();
//        if(seen_no_individuals >= metacommunity_size)
//        {
//            break;
//        }
//    }
    //    while(seen_no_individuals < metacommunity_size)
    //    {
    //        addNewSpecies();
    //    }
    // Make sure that we've seen at least as many individuals as in the local community.
    if(seen_no_individuals < local_community_size)
    {
        stringstream ss;
        ss << "Seen number of individuals (" << seen_no_individuals << ") is not more than local community size (";
        ss << local_community_size << ") - please report this bug" << endl;
        throw FatalException(ss.str());
    }
    writeInfo("done.\n");

}

unsigned long AnalyticalSpeciesAbundancesHandler::getRandomSpeciesID()
{
    // Select a random individual from the seen number of individuals
    auto individual_id = random->i0(metacommunity_size- 1);

    // Pick out a new individual
    if(individual_id >= seen_no_individuals) // TODO clean up or remove
    {
        addNewSpecies();
        stringstream ss; // TODO remove
        ss << "Added new species for id of " << individual_id << " with species id of " << ind_to_species.rbegin()->second;
        ss << " and abundance of " << ind_to_species.rbegin()->first <<endl;
        writeInfo(ss.str());
        return ind_to_species.rbegin()->second;
    }
    stringstream ss; // TODO remove
    ss << "Individual id was " << individual_id << " with species id " << pickPreviousIndividual(individual_id) << endl;
    writeInfo(ss.str());

#ifdef DEBUG
    if(individual_id > seen_no_individuals)
    {
        stringstream ss;
        ss << "Random individual ID (" << individual_id << ") is greater than the number of seen indiviudals (";
        ss << seen_no_individuals << ")" << endl;
        throw FatalException(ss.str());
    }
    if(ind_to_species.empty())
    {
        throw FatalException(
                "No individuals have been seen yet, but an individual ID was generated. Please report this bug.");
    }
#endif // DEBUG
    return pickPreviousIndividual(individual_id);
}

unsigned long AnalyticalSpeciesAbundancesHandler::pickPreviousIndividual(const unsigned long &individual_id)
{
    return ind_to_species.upper_bound(individual_id)->second;
}

void AnalyticalSpeciesAbundancesHandler::addNewSpecies()
{
    max_species_id++;
    unsigned long new_abundance = 0;
    do
    {
        new_abundance = getRandomAbundanceOfSpecies();
    }while(new_abundance > metacommunity_size - seen_no_individuals);
    unsigned long cumulative_abundance;
    if(seen_no_individuals == 0)
    {
        cumulative_abundance = new_abundance;
    }
    else
    {
        cumulative_abundance = new_abundance + ind_to_species.rbegin()->first;
    }
    ind_to_species[cumulative_abundance] = max_species_id;
    seen_no_individuals += new_abundance;
    if(ind_to_species.rbegin()->first != seen_no_individuals)
    {
        stringstream ss;
        ss << "ind_to_species end does not equal seen no inds: " << ind_to_species.rbegin()->first << "!=";
        ss << seen_no_individuals << endl;
        throw FatalException(ss.str());
    }
    if(ind_to_species.rbegin()->second != max_species_id)
    {
        stringstream ss;
        ss << "Last species id has not been set correctly: " << ind_to_species.rbegin()->second << "!=";
        ss << max_species_id << endl;
        throw FatalException(ss.str());
    }
}

unsigned long AnalyticalSpeciesAbundancesHandler::getRandomAbundanceOfSpecies()
{
    // First generate a random abundance class
    return static_cast<unsigned long>(max(static_cast<double>(
                                                  min(random->randomLogarithmic(1.0 - speciation_rate),
                                                      metacommunity_size)), 1.0));
}

unsigned long AnalyticalSpeciesAbundancesHandler::getRandomAbundanceOfIndividual()
{
    long double p = 1 - speciation_rate;
    unsigned long unseen_individuals = metacommunity_size - seen_no_individuals;
    long double val = round(log(((random->d01() * (pow(p, unseen_individuals) - 1)) + 1)) / log(p));
    unsigned long res = max(min(static_cast<unsigned long>(val), unseen_individuals),
                            (unsigned long) 1);
//    stringstream ss;
//    ss << "speciation rate is " << speciation_rate << endl;
//    ss << "metacommunity size is " << metacommunity_size << endl;
//    ss << "seen individuals: " << seen_no_individuals << endl;
////    ss << "unseen individuals: " << unseen_individuals << endl;
//    ss << "Results was " << res << " for val " << val << endl; // TODO remove
//    writeInfo(ss.str());
    return res;
}


