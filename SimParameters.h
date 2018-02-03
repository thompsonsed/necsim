//
// Created by Sam Thompson on 10/09/2017.
//

#ifndef SPECIATIONCOUNTER_SIMPARAMETERS_H
#define SPECIATIONCOUNTER_SIMPARAMETERS_H
#include <string>
#include <vector>
#include "ConfigFileParser.h"
#include "Logging.h"

/************************************************************
					MAPVARS STRUCTURE
 ************************************************************/
/**
 * @struct SimParameters
 * @brief Stores and imports the variables required by the Map object.
 * Used to setting the Map variables in a more elegant way.
 */
struct SimParameters
{
	string finemapfile, coarsemapfile,outdirectory,pristinefinemapfile,pristinecoarsemapfile, samplemaskfile;
	 // for file naming purposes.
	long long the_task{}, the_seed{};
	// the variables for the grid containing the initial individuals.
	unsigned long vargridxsize{}, vargridysize{};
	// The variables for the sample grid, which may or may not be the same as the main simulation grid
	unsigned long varsamplexsize{}, varsampleysize{};
	unsigned long varsamplexoffset{}, varsampleyoffset{};
	// The fine map variables at the same resolution as the grid.
	unsigned long varfinemapxsize{}, varfinemapysize{}, varfinemapxoffset{}, varfinemapyoffset{};
	// the coarse map variables at a scaled resolution of the fine map.
	unsigned long varcoarsemapxsize{}, varcoarsemapysize{}, varcoarsemapxoffset{}, varcoarsemapyoffset{};
	unsigned long varcoarsemapscale{};
	unsigned long desired_specnum{};
	// the relative cost of moving through non-forest
	double dispersal_relative_cost{};
	// the size of each square of habitat in numbers of individuals
	unsigned long deme{};
	// the sample proportion,
	 double deme_sample{};
	// the speciation rate.
	long double  spec{};
	// the variance of the dispersal kernel.
	double sigma{};
	// max time to run for
	unsigned long maxtime{};
	// the number of generations since a pristine landscape was encountered.
	double gen_since_pristine{};
	// the transform rate of the forest from pristine to modern forest.
	double habitat_change_rate{};
	// the fatness of the dispersal kernel
	double tau{};
	// dispersal method - should be one of [normal, fat-tail, norm-uniform]
	string dispersal_method;
	// the probability of selecting from a uniform dispersal kernel (for uniformally-modified dispersals)
	double m_prob{};
	// the cutoff for the normal dispersal in cells.
	double cutoff{};
	// if true, restricts dispersal from the same cell.
	bool restrict_self{};
	// file containing the points to record data from
	string times_file;
	// Stores the full list of configs imported from file
	ConfigOption configs;
	// Set to true if the completely pristine state has been reached.
	bool is_pristine{};
	// if the sample file is not null, this variable tells us whether different points in space require different
	// numbers of individuals to be sampled. If this is the case, the actual values are read from the sample mask as a
	// proportion of individuals sampled, from 0-1. Otherwise, it is treated as a boolean mask, with values > 0.5
	// representing sampling in the cell.
	bool uses_spatial_sampling{};
	// This can be closed, infinite and tiled (which is also infinite)
	string landscape_type;
	// The protracted speciation parameters - these DON'T need to be stored upon pausing simulations
	bool is_protracted{};
	double min_speciation_gen{}, max_speciation_gen{};

	// a map of dispersal values, where each row corresponds to the probability of moving from one cell
	// to any other.
	string dispersal_file;

	// a map of relative reproduction probabilities.
	string reproduction_file;

	/**
	 * @brief Default constructor
	 */
	SimParameters()
	{
		finemapfile = "none";
		coarsemapfile = "none";
		outdirectory = "none";
		pristinefinemapfile = "none";
		pristinecoarsemapfile = "none";
		samplemaskfile = "none";
		times_file = "null";
		dispersal_method = "none";
		landscape_type = "none";
		reproduction_file = "none";
		dispersal_file = "none";
		min_speciation_gen = 0.0;
		max_speciation_gen = 0.0;
		is_protracted = false;
		restrict_self = false;
		m_prob = 0;
		cutoff = 0;
		tau =0;
	}

	/**
	 * @brief Imports the spatial variables from a vector of command line arguments.
	 * If importing from a config file, the comargs should just contain the path to the config file to import from.
	 * @param comargs vector of strings containing the command line arguments
	 */
	void importParameters(const string &conf_in)
	{
		// do the importSpatialParameters of the values from combination of command-line arguments and file.
		configs.setConfig(conf_in, false, true);
		configs.parseConfig();
		varsamplexsize = stoul(configs.getSectionOptions("sample_grid", "x", "0"));
		varsampleysize = stoul(configs.getSectionOptions("sample_grid", "y", "0"));
		varsamplexoffset = stoul(configs.getSectionOptions("sample_grid", "x_off", "0"));
		varsampleyoffset = stoul(configs.getSectionOptions("sample_grid", "y_off", "0"));
		uses_spatial_sampling = static_cast<bool>(stoi(configs.getSectionOptions("sample_grid",
																				 "uses_spatial_sampling", "0")));
		if(configs.hasSection("grid_map"))
		{
			vargridxsize = stoul(configs.getSectionOptions("grid_map", "x"));
			vargridysize = stoul(configs.getSectionOptions("grid_map", "y"));
		}
		else
		{
			vargridxsize = varsamplexsize;
			vargridysize = varsampleysize;
		}
		samplemaskfile = configs.getSectionOptions("sample_grid","mask", "null");
		finemapfile = configs.getSectionOptions("fine_map", "path", "none");
		varfinemapxsize = stoul(configs.getSectionOptions("fine_map", "x", "0"));
		varfinemapysize = stoul(configs.getSectionOptions("fine_map", "y", "0"));
		varfinemapxoffset = stoul(configs.getSectionOptions("fine_map", "x_off", "0"));
		varfinemapyoffset = stoul(configs.getSectionOptions("fine_map", "y_off", "0"));
		coarsemapfile = configs.getSectionOptions("coarse_map", "path", "none");
		varcoarsemapxsize = stoul(configs.getSectionOptions("coarse_map", "x", "0"));
		varcoarsemapysize = stoul(configs.getSectionOptions("coarse_map", "y", "0"));
		varcoarsemapxoffset = stoul(configs.getSectionOptions("coarse_map", "x_off", "0"));
		varcoarsemapyoffset = stoul(configs.getSectionOptions("coarse_map", "y_off", "0"));
		varcoarsemapscale = stoul(configs.getSectionOptions("coarse_map", "scale", "0"));
		pristinefinemapfile = configs.getSectionOptions("pristine_fine0", "path", "none");
		pristinecoarsemapfile = configs.getSectionOptions("pristine_coarse0", "path", "none");
		dispersal_method = configs.getSectionOptions("dispersal", "method", "none");
		m_prob = stod(configs.getSectionOptions("dispersal", "m_probability", "0"));
		cutoff = stod(configs.getSectionOptions("dispersal", "cutoff", "0.0"));
		// quick and dirty conversion for string to bool
		restrict_self = static_cast<bool>(stoi(configs.getSectionOptions("dispersal", "restrict_self", "0")));
		landscape_type = configs.getSectionOptions("dispersal", "infinite_landscape", "none");
		dispersal_file = configs.getSectionOptions("dispersal", "dispersal_file", "none");
		reproduction_file = configs.getSectionOptions("reproduction", "map", "none");
		outdirectory = configs.getSectionOptions("main", "output_directory", "Default");
		the_seed = stol(configs.getSectionOptions("main", "seed", "0"));
		the_task = stol(configs.getSectionOptions("main", "job_type", "0"));
		tau = stod(configs.getSectionOptions("main", "tau", "0.0"));
		sigma = stod(configs.getSectionOptions("main", "sigma", "0.0"));
		deme = stoul(configs.getSectionOptions("main", "deme"));
		deme_sample = stod(configs.getSectionOptions("main", "sample_size"));
		maxtime = stoul(configs.getSectionOptions("main", "max_time"));
		dispersal_relative_cost = stod(configs.getSectionOptions("main", "dispersal_relative_cost", "0"));
		times_file = configs.getSectionOptions("main", "time_config");
		spec = stod(configs.getSectionOptions("main", "min_spec_rate"));
		desired_specnum = stoul(configs.getSectionOptions("main", "min_species", "1"));
		if(configs.hasSection("protracted"))
		{
			is_protracted = static_cast<bool>(stoi(
					configs.getSectionOptions("protracted", "has_protracted", "0")));
			min_speciation_gen = stod(configs.getSectionOptions("protracted", "min_speciation_gen", "0.0"));
			max_speciation_gen = stod(configs.getSectionOptions("protracted", "max_speciation_gen"));
		}
		setPristine(0);
	}

	/**
	 * @brief Alters the pristine parameters to the configuration matching the input number. If no configuration option
	 * exists for this number, bPristine will be set to true.
	 * @param n the pristine map number to check.
	 * @return bool true if we need to re-import the maps (i.e. the pristine maps have changed between updates)
	 */

	bool setPristine(unsigned int n)
	{
		is_pristine = true;
		bool finemapcheck = false;
		bool coarsemapcheck = false;
		// Loop over each element in the config file (each line) and check if it is pristine fine or pristine coarse.
		for(unsigned int i = 0; i < configs.getSectionOptionsSize(); i ++ )
		{

			if(configs[i].section.find("pristine_fine") == 0)
			{
				// Then loop over each element to find the number, and check if it is equal to our input number.
				is_pristine = false;
				if(stol(configs[i].getOption("number")) == n)
				{
					string tmpmapfile;
					tmpmapfile = configs[i].getOption("path");
					if(pristinefinemapfile != tmpmapfile)
					{
						finemapcheck = true;
						pristinefinemapfile = tmpmapfile;
					}
					habitat_change_rate = stod(configs[i].getOption("rate"));
					gen_since_pristine = stod(configs[i].getOption("time"));
				}
			}
			else if(configs[i].section.find("pristine_coarse") == 0)
			{
				if(stol(configs[i].getOption("number")) == n)
				{
					string tmpmapfile;
					tmpmapfile = configs[i].getOption("path");
					is_pristine = false;
					if(tmpmapfile != pristinecoarsemapfile)
					{
						coarsemapcheck=true;
						pristinecoarsemapfile = tmpmapfile;
						// check matches
						if(habitat_change_rate != stod(configs[i].getOption("rate")) || gen_since_pristine != stod(configs[i].getOption("time")))
						{
							cerr << "Forest transform values do not match between fine and coarse maps. Using fine values." << endl;
						}
					}
				}
			}
		}
		// if one of the maps has changed, we need to update, so return true.
		if(finemapcheck != coarsemapcheck)
		{
			return true;
		}
		else
		{
			// finemapcheck should therefore be the same as coarsemapcheck
			return finemapcheck;
		}
	}
	/**
	 * @brief Prints selected important variables to the terminal.
	 */
	void printVars()
	{
		stringstream os;
		os << "Seed: " << the_seed << endl;
		os << "Speciation rate: " << spec << endl;
		os << "Dispersal (tau, sigma): " << tau << ", " << sigma << endl;
		os << "Dispersal method: " << dispersal_method << endl;
		if(dispersal_method == "norm-uniform")
		{
			os << "Dispersal (m, cutoff): " << m_prob << ", " << cutoff << endl;
		}
		if(is_protracted)
		{
			os << "Protracted variables: " << min_speciation_gen << ", " << max_speciation_gen << endl;
		}
		os << "Job Type: " << the_task << endl;
		os << "Max time: " << maxtime << endl;
		os << "Fine input file: " << finemapfile  << endl;
		os << "-dimensions: (" << varfinemapxsize << ", " << varfinemapysize <<")"<< endl;
		os << "-offset: (" << varfinemapxoffset << ", " << varfinemapyoffset << ")" << endl;
		os << "Coarse input file: " << coarsemapfile  << endl;
		os << "-dimensions: (" << varcoarsemapxsize << ", " << varcoarsemapysize <<")"<< endl;
		os << "-offset: (" << varcoarsemapxoffset << ", " << varcoarsemapyoffset << ")" << endl;
		os << "-scale: " << varcoarsemapscale << endl;
		os << "Sample grid" << endl;
		os << "-dimensions: (" << varsamplexsize << ", " << varsampleysize << ")" << endl;
		os << "-optimised area: (" << vargridxsize << ", " << vargridysize << ")" << endl;
		os << "-optimised offsets: (" << varsamplexoffset << ", " << varsampleyoffset << ")" << endl;
		os << "-deme: " << deme << endl;
		os << "-deme sample: " << deme_sample << endl;
		os << "Output directory: " << outdirectory << endl;
		os << "Disp Rel Cost: " << dispersal_relative_cost << endl;
		writeInfo(os.str());
	}

	void setMetacommunityParameters(const unsigned long &metacommunity_size,
									const double &speciation_rate,
									const unsigned long &seed,
									const unsigned long &job)
	{
		outdirectory = "Default";
		// randomise the seed slightly so that we get a different starting number to the initial simulation
		the_seed = static_cast<long long int>(seed * job);
		the_task = (long long int) job;
		deme = metacommunity_size;
		deme_sample = 1.0;
		spec = speciation_rate;
		// Default to 1000 seconds - should be enough for most simulation sizes, but can be changed later if needed.
		maxtime = 1000;
		times_file = "null";
		min_speciation_gen = 0.0;
		max_speciation_gen = 0.0;
	}

	/**
	 * @brief Overloading the << operator for outputting to the output stream
	 * @param os the output stream.
	 * @param m the SimParameters object.
	 * @return os the output stream.
	 */
	friend ostream& operator<<(ostream& os,const SimParameters& m)
	{
		os << m.finemapfile << "\n" << m.coarsemapfile << "\n" << m.pristinefinemapfile << "\n";
		os << m.pristinecoarsemapfile << "\n" << m.samplemaskfile << "\n";
		os << m.the_seed << "\n" <<  m.the_task << "\n" <<  m.vargridxsize << "\n" << m.vargridysize << "\n";
		os << m.varsamplexsize << "\n" << m.varsampleysize << "\n" << m.varsamplexoffset << "\n" << m.varsampleyoffset << "\n";
		os << m.varfinemapxsize << "\n" << m.varfinemapysize << "\n";
		os << m.varfinemapxoffset << "\n" << m.varfinemapyoffset << "\n" << m.varcoarsemapxsize << "\n" << m.varcoarsemapysize << "\n" << m.varcoarsemapxoffset << "\n";
		os << m.varcoarsemapyoffset << "\n" << m.varcoarsemapscale << "\n" << m.desired_specnum << "\n";
		os << m.dispersal_relative_cost << "\n" << m.deme << "\n" << m.deme_sample<< "\n";
		os << m.spec << "\n" << m.sigma << "\n" << m.maxtime << "\n" << m.gen_since_pristine << "\n" << m. habitat_change_rate << "\n" << m.tau;
		os << "\n" << m.dispersal_method << "\n";
		os << m.m_prob << "\n" << m.cutoff << "\n" << m.restrict_self <<"\n" << m.landscape_type << "\n" << m.times_file << "\n";
		os << m.dispersal_file << "\n" << m.uses_spatial_sampling << "\n";
		os << m.configs;
		return os;
	}

	/**
	 * @brief Overloading the >> operator for inputting from an input stream
	 * @param is the input stream
	 * @param m the mapvars object
	 * @return is the input stream
	 */
	friend istream& operator>>(istream& is, SimParameters& m)
	{
		getline(is, m.finemapfile);
		getline(is, m.coarsemapfile);
		getline(is, m.pristinefinemapfile);
		getline(is, m.pristinecoarsemapfile);
		getline(is, m.samplemaskfile);
		is >> m.the_seed >> m.the_task >>  m.vargridxsize >> m.vargridysize;
		is >> m.varsamplexsize >> m.varsampleysize >> m.varsamplexoffset >> m.varsampleyoffset;
		is >> m.varfinemapxsize >> m.varfinemapysize;
		is >> m.varfinemapxoffset >> m.varfinemapyoffset >> m.varcoarsemapxsize >> m.varcoarsemapysize >> m.varcoarsemapxoffset ;
		is >> m.varcoarsemapyoffset >> m.varcoarsemapscale >> m.desired_specnum >> m.dispersal_relative_cost >> m.deme >> m.deme_sample;
		is >> m.spec >> m.sigma >> m.maxtime >> m.gen_since_pristine >> m.habitat_change_rate >> m.tau;
		is.ignore();
		getline(is, m.dispersal_method);
		is >> m.m_prob >> m.cutoff >> m.restrict_self >> m.landscape_type;
		is.ignore();
		getline(is, m.times_file);
		getline(is, m.dispersal_file);
		is >> m.uses_spatial_sampling;
		is >> m.configs;
		return is;
	}
};


#endif //SPECIATIONCOUNTER_SIMPARAMETERS_H
