//
// Created by Sam Thompson on 10/09/2017.
//

#ifndef SPECIATIONCOUNTER_SIMPARAMETERS_H
#define SPECIATIONCOUNTER_SIMPARAMETERS_H
#include <string>
#include <vector>
#include "ConfigFileParser.h"
#include "Logging.h"
using namespace std;
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
	string fine_map_file, coarse_map_file, output_directory;
	string pristine_fine_map_file, pristine_coarse_map_file, sample_mask_file;
	 // for file naming purposes.
	long long the_task{}, the_seed{};
	// the variables for the grid containing the initial individuals.
	unsigned long grid_x_size{}, grid_y_size{};
	// The variables for the sample grid, which may or may not be the same as the main simulation grid
	unsigned long sample_x_size{}, sample_y_size{};
	unsigned long sample_x_offset{}, sample_y_offset{};
	// The fine map variables at the same resolution as the grid.
	unsigned long fine_map_x_size{}, fine_map_y_size{}, fine_map_x_offset{}, fine_map_y_offset{};
	// the coarse map variables at a scaled resolution of the fine map.
	unsigned long coarse_map_x_size{}, coarse_map_y_size{}, coarse_map_x_offset{}, coarse_map_y_offset{};
	unsigned long coarse_map_scale{};
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
	unsigned long max_time{};
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
		fine_map_file = "none";
		coarse_map_file = "none";
		output_directory = "none";
		pristine_fine_map_file = "none";
		pristine_coarse_map_file = "none";
		sample_mask_file = "none";
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
	 * @brief Links to the provided ConfigOption.
	 * Assumes that the parameters have already been parsed from the config file.
	 * @param configOption the pointer to the parsed ConfigOption object
	 */
	void importParameters(ConfigOption *configOption)
	{
		configs = *configOption;
		importParameters();
	}

	/**
	 * @brief Imports the spatial variables from a path to the config file..
	 * @param config_in string of the path to the config file
	 */
	void importParameters(const string &conf_in)
	{
		// do the importSpatialParameters of the values from combination of command-line arguments and file.
		configs.setConfig(conf_in, false, true);
		configs.parseConfig();
		importParameters();
	}

	/**
	 * @brief Main import of parameters from the config file option.
	 */
	void importParameters()
	{
		sample_x_size = stoul(configs.getSectionOptions("sample_grid", "x", "0"));
		sample_y_size = stoul(configs.getSectionOptions("sample_grid", "y", "0"));
		sample_x_offset = stoul(configs.getSectionOptions("sample_grid", "x_off", "0"));
		sample_y_offset = stoul(configs.getSectionOptions("sample_grid", "y_off", "0"));
		uses_spatial_sampling = static_cast<bool>(stoi(configs.getSectionOptions("sample_grid",
																				 "uses_spatial_sampling", "0")));
		if(configs.hasSection("grid_map"))
		{
			grid_x_size = stoul(configs.getSectionOptions("grid_map", "x"));
			grid_y_size = stoul(configs.getSectionOptions("grid_map", "y"));
		}
		else
		{
			grid_x_size = sample_x_size;
			grid_y_size = sample_y_size;
		}
		sample_mask_file = configs.getSectionOptions("sample_grid","mask", "null");
		fine_map_file = configs.getSectionOptions("fine_map", "path", "none");
		fine_map_x_size = stoul(configs.getSectionOptions("fine_map", "x", "0"));
		fine_map_y_size = stoul(configs.getSectionOptions("fine_map", "y", "0"));
		fine_map_x_offset = stoul(configs.getSectionOptions("fine_map", "x_off", "0"));
		fine_map_y_offset = stoul(configs.getSectionOptions("fine_map", "y_off", "0"));
		coarse_map_file = configs.getSectionOptions("coarse_map", "path", "none");
		coarse_map_x_size = stoul(configs.getSectionOptions("coarse_map", "x", "0"));
		coarse_map_y_size = stoul(configs.getSectionOptions("coarse_map", "y", "0"));
		coarse_map_x_offset = stoul(configs.getSectionOptions("coarse_map", "x_off", "0"));
		coarse_map_y_offset = stoul(configs.getSectionOptions("coarse_map", "y_off", "0"));
		coarse_map_scale = stoul(configs.getSectionOptions("coarse_map", "scale", "0"));
		pristine_fine_map_file = configs.getSectionOptions("pristine_fine0", "path", "none");
		pristine_coarse_map_file = configs.getSectionOptions("pristine_coarse0", "path", "none");
		dispersal_method = configs.getSectionOptions("dispersal", "method", "none");
		m_prob = stod(configs.getSectionOptions("dispersal", "m_probability", "0"));
		cutoff = stod(configs.getSectionOptions("dispersal", "cutoff", "0.0"));
		// quick and dirty conversion for string to bool
		restrict_self = static_cast<bool>(stoi(configs.getSectionOptions("dispersal", "restrict_self", "0")));
		landscape_type = configs.getSectionOptions("dispersal", "infinite_landscape", "none");
		dispersal_file = configs.getSectionOptions("dispersal", "dispersal_file", "none");
		reproduction_file = configs.getSectionOptions("reproduction", "map", "none");
		output_directory = configs.getSectionOptions("main", "output_directory", "Default");
		the_seed = stol(configs.getSectionOptions("main", "seed", "0"));
		the_task = stol(configs.getSectionOptions("main", "job_type", "0"));
		tau = stod(configs.getSectionOptions("main", "tau", "0.0"));
		sigma = stod(configs.getSectionOptions("main", "sigma", "0.0"));
		deme = stoul(configs.getSectionOptions("main", "deme"));
		deme_sample = stod(configs.getSectionOptions("main", "sample_size"));
		max_time = stoul(configs.getSectionOptions("main", "max_time"));
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
	 * @brief Sets the main simulation parameters
	 * @param task_in the task reference number, used for file referencing
	 * @param seed_in the seed to set random number generation
	 * @param output_directory_in the output directory
	 * @param max_time_in the maximum time to simulate for
	 * @param desired_specnum_in the desired number of species to aim towards (currently not functional)
	 * @param times_file_in the file containing a list of temporal sampling points
	 */
	void setKeyParameters(const long long &task_in, const long long &seed_in, const string &output_directory_in,
						  const unsigned long &max_time_in, unsigned long desired_specnum_in, const string &times_file_in)
	{
		the_task = task_in;
		the_seed = seed_in;
		output_directory = output_directory_in;
		max_time = max_time_in;
		desired_specnum = desired_specnum_in;
		times_file = times_file_in;

	}

	/**
	 * @brief Sets the speciation parameters for the simulation.
	 * @param spec_in the speciation rate to use
	 * @param is_protracted_in if true, simulates as a protracted simulation
	 * @param min_speciation_gen_in the minimum speciation generation for protracted simulations
	 * @param max_speciation_gen_in the maximum speciation generation for protracted simulations
	 */
	void setSpeciationParameters(const long double &spec_in, bool is_protracted_in, const double &min_speciation_gen_in,
								 const double &max_speciation_gen_in)
	{
		spec = spec_in;
		is_protracted = is_protracted_in;
		min_speciation_gen = min_speciation_gen_in;
		max_speciation_gen = max_speciation_gen_in;
	}

	/**
	 * @brief Sets the dispersal parameters for the simulation.
	 * @param dispersal_method_in the method of individuals dispersing (normal, fat-tailed or norm-uniform)
	 * @param sigma_in the sigma value for a normal distribution
	 * @param tau_in the tau value for the fat-tailed distribution
	 * @param m_prob_in the probability of uniform dispersal for the norm-uniform distribution
	 * @param cutoff_in the maximum dispersal distance for the uniform distribution
	 * @param dispersal_relative_cost_in the relative cost of dispersing through non-forest
	 * @param restrict_self_in if true, prevents dispersal from the same cell
	 * @param landscape_type_in the landscape type (infinite, tiled or closed)
	 * @param dispersal_file_in a map of dispersal probabilities
	 * @param reproduction_file_in a map of reproduction probabilities
	 */
	void setDispersalParameters(const string &dispersal_method_in, const double &sigma_in, const double &tau_in,
								const double &m_prob_in, const double &cutoff_in,
								const double &dispersal_relative_cost_in, bool restrict_self_in,
								const string &landscape_type_in, const string &dispersal_file_in,
								const string &reproduction_file_in)
	{
		dispersal_method = dispersal_method_in;
		sigma = sigma_in;
		tau = tau_in;
		m_prob = m_prob_in;
		cutoff = cutoff_in;
		dispersal_relative_cost = dispersal_relative_cost_in;
		restrict_self = restrict_self_in;
		landscape_type = landscape_type_in;
		dispersal_file = dispersal_file_in;
		reproduction_file = reproduction_file_in;
	}

	/**
	 * @brief Sets the pristine map parameters for the simulation.
	 * @param pristine_fine_file_map_in the fine resolution pristine file
	 * @param pristine_coarse_map_file_in the coarse resolution pristine file
	 * @param gen_since_pristine_in the number of generations since the pristine state was achieved
	 * @param habitat_change_rate_in the rate of habitat change towards the pristine state
	 */
	void setPristineMapParameters(const string &pristine_fine_file_map_in, const string &pristine_coarse_map_file_in,
								  const double &gen_since_pristine_in, const double &habitat_change_rate_in)
	{
		pristine_fine_map_file = pristine_fine_file_map_in;
		pristine_coarse_map_file = pristine_coarse_map_file_in;
		gen_since_pristine = gen_since_pristine_in;
		habitat_change_rate = habitat_change_rate_in;
	}

	/**
	 * @brief Sets the map parameters for the simulation.
	 * @param fine_map_file_in the fine resolution density map
	 * @param coarse_map_file_in the coarse resolution density map
	 * @param sample_mask_file_in the spatial sampling mask
	 * @param grid_x_size_in the x dimension of the grid
	 * @param grid_y_size_in the y dimension of the grid
	 * @param sample_x_size_in the x dimension of the sample mask
	 * @param sample_y_size_in the y dimension of the sample mask
	 * @param sample_x_offset_in the x offset of the sample mask from the grid
	 * @param sample_y_offset_in the y offset of the sample mask from the grid
	 * @param fine_map_x_size_in the x dimension of the fine map
	 * @param fine_map_y_size_in the y dimension of the fine map
	 * @param fine_map_x_offset_in the x offset of the fine map from the sample mask
	 * @param fine_map_y_offset_in the y offset of the fine map from the sample mask
	 * @param coarse_map_x_size_in the x dimension of the coarse map
	 * @param coarse_map_y_size_in the y dimension of the coarse map
	 * @param coarse_map_x_offset_in the x offset of the coarse map from the fine map
	 * @param coarse_map_y_offset_in the y offset of the coarse map from the fine map
	 * @param coarse_map_scale_in the scale of the coarse map compared to the fine map
	 * @param deme_in the number of individuals per cell
	 * @param deme_sample_in the proportion of individuals to sample from each cell
	 * @param uses_spatial_sampling_in if the sample mask denotes differing spatial sampling proportions
	 */
	void setMapParameters(const string &fine_map_file_in, const string &coarse_map_file_in,
						  const string &sample_mask_file_in, const unsigned long &grid_x_size_in,
						  const unsigned long &grid_y_size_in, const unsigned long &sample_x_size_in,
						  const unsigned long &sample_y_size_in, const unsigned long &sample_x_offset_in,
						  const unsigned long &sample_y_offset_in, const unsigned long &fine_map_x_size_in,
						  const unsigned long &fine_map_y_size_in, const unsigned long &fine_map_x_offset_in,
						  const unsigned long &fine_map_y_offset_in, const unsigned long &coarse_map_x_size_in,
						  const unsigned long &coarse_map_y_size_in, const unsigned long &coarse_map_x_offset_in,
						  const unsigned long &coarse_map_y_offset_in, const unsigned long &coarse_map_scale_in,
						  const unsigned long &deme_in, const double &deme_sample_in, bool uses_spatial_sampling_in)
	{
		fine_map_file = fine_map_file_in;
		coarse_map_file = coarse_map_file_in;
		sample_mask_file = sample_mask_file_in;
		grid_x_size = grid_x_size_in;
		grid_y_size = grid_y_size_in;
		sample_x_size = sample_x_size_in;
		sample_y_size = sample_y_size_in;
		sample_x_offset = sample_x_offset_in;
		sample_y_offset = sample_y_offset_in;
		fine_map_x_size = fine_map_x_size_in;
		fine_map_y_size = fine_map_y_size_in;
		fine_map_x_offset = fine_map_x_offset_in;
		fine_map_y_offset = fine_map_y_offset_in;
		coarse_map_x_size = coarse_map_x_size_in;
		coarse_map_y_size = coarse_map_y_size_in;
		coarse_map_x_offset = coarse_map_x_offset_in;
		coarse_map_y_offset = coarse_map_y_offset_in;
		coarse_map_scale = coarse_map_scale_in;
		deme = deme_in;
		deme_sample = deme_sample_in;
		uses_spatial_sampling = uses_spatial_sampling_in;
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
					if(pristine_fine_map_file != tmpmapfile)
					{
						finemapcheck = true;
						pristine_fine_map_file = tmpmapfile;
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
					if(tmpmapfile != pristine_coarse_map_file)
					{
						coarsemapcheck=true;
						pristine_coarse_map_file = tmpmapfile;
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
		os << "Max time: " << max_time << endl;
		os << "Fine input file: " << fine_map_file  << endl;
		os << "-dimensions: (" << fine_map_x_size << ", " << fine_map_y_size <<")"<< endl;
		os << "-offset: (" << fine_map_x_offset << ", " << fine_map_y_offset << ")" << endl;
		os << "Coarse input file: " << coarse_map_file  << endl;
		os << "-dimensions: (" << coarse_map_x_size << ", " << coarse_map_y_size <<")"<< endl;
		os << "-offset: (" << coarse_map_x_offset << ", " << coarse_map_y_offset << ")" << endl;
		os << "-scale: " << coarse_map_scale << endl;
		os << "Sample grid" << endl;
		os << "-dimensions: (" << sample_x_size << ", " << sample_y_size << ")" << endl;
		os << "-optimised area: (" << grid_x_size << ", " << grid_y_size << ")" << endl;
		os << "-optimised offsets: (" << sample_x_offset << ", " << sample_y_offset << ")" << endl;
		os << "-deme: " << deme << endl;
		os << "-deme sample: " << deme_sample << endl;
		os << "Output directory: " << output_directory << endl;
		os << "Disp Rel Cost: " << dispersal_relative_cost << endl;
		writeInfo(os.str());
	}

	void setMetacommunityParameters(const unsigned long &metacommunity_size,
									const double &speciation_rate,
									const unsigned long &seed,
									const unsigned long &job)
	{
		output_directory = "Default";
		// randomise the seed slightly so that we get a different starting number to the initial simulation
		the_seed = static_cast<long long int>(seed * job);
		the_task = (long long int) job;
		deme = metacommunity_size;
		deme_sample = 1.0;
		spec = speciation_rate;
		// Default to 1000 seconds - should be enough for most simulation sizes, but can be changed later if needed.
		max_time = 1000;
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
		os << m.fine_map_file << "\n" << m.coarse_map_file << "\n" << m.pristine_fine_map_file << "\n";
		os << m.pristine_coarse_map_file << "\n" << m.sample_mask_file << "\n";
		os << m.the_seed << "\n" <<  m.the_task << "\n" <<  m.grid_x_size << "\n" << m.grid_y_size << "\n";
		os << m.sample_x_size << "\n" << m.sample_y_size << "\n" << m.sample_x_offset << "\n" << m.sample_y_offset << "\n";
		os << m.fine_map_x_size << "\n" << m.fine_map_y_size << "\n";
		os << m.fine_map_x_offset << "\n" << m.fine_map_y_offset << "\n" << m.coarse_map_x_size << "\n" << m.coarse_map_y_size << "\n" << m.coarse_map_x_offset << "\n";
		os << m.coarse_map_y_offset << "\n" << m.coarse_map_scale << "\n" << m.desired_specnum << "\n";
		os << m.dispersal_relative_cost << "\n" << m.deme << "\n" << m.deme_sample<< "\n";
		os << m.spec << "\n" << m.sigma << "\n" << m.max_time << "\n" << m.gen_since_pristine << "\n" << m. habitat_change_rate << "\n" << m.tau;
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
		getline(is, m.fine_map_file);
		getline(is, m.coarse_map_file);
		getline(is, m.pristine_fine_map_file);
		getline(is, m.pristine_coarse_map_file);
		getline(is, m.sample_mask_file);
		is >> m.the_seed >> m.the_task >>  m.grid_x_size >> m.grid_y_size;
		is >> m.sample_x_size >> m.sample_y_size >> m.sample_x_offset >> m.sample_y_offset;
		is >> m.fine_map_x_size >> m.fine_map_y_size;
		is >> m.fine_map_x_offset >> m.fine_map_y_offset >> m.coarse_map_x_size >> m.coarse_map_y_size >> m.coarse_map_x_offset ;
		is >> m.coarse_map_y_offset >> m.coarse_map_scale >> m.desired_specnum >> m.dispersal_relative_cost >> m.deme >> m.deme_sample;
		is >> m.spec >> m.sigma >> m.max_time >> m.gen_since_pristine >> m.habitat_change_rate >> m.tau;
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
