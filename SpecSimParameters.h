//
// Created by Sam Thompson on 08/01/2018.
//

#ifndef SPECIATIONCOUNTER_SPECSIMPARAMETERS_H
#define SPECIATIONCOUNTER_SPECSIMPARAMETERS_H


/**
 * @class SpecSimParameters
 * @brief Contains the simulation parameters that are read from the command line.
 *
 */
struct SpecSimParameters
{
	bool use_spatial;
	bool bMultiRun;
	bool use_fragments;
	string filename;
	vector<double> all_speciation_rates;
	string samplemask;
	string times_file;
	vector<double> all_times;
	string fragment_config_file;
	double min_speciation_gen, max_speciation_gen;
	unsigned long metacommunity_size;
	double metacommunity_speciation_rate;

	/**
	 * @brief Sets the application arguments for the inputs. Intended for use with the applyspecmodule for
	 * integration with python.
	 *
	 * @param file_in the database to apply speciation rates to
	 * @param use_spatial_in if true, record full spatial data
	 * @param sample_file the sample file to select lineages from the map
	 * @param time_config the time config file to use
	 * @param use_fragments_in fragment file, or "T"/"F" for automatic detection/no detection
	 * @param speciation_rates the speciation rates to apply
	 * @param min_speciation_gen_in the minimum generation rate for speciation in protracted simulations
	 * @param max_speciation_gen_in the maximum generation rate for speciation in protracted simulations
	 */
	void setup(string file_in, bool use_spatial_in, string sample_file, string time_config, string use_fragments_in,
			   vector<double> speciation_rates, double min_speciation_gen_in, double max_speciation_gen_in)
	{
		setup(file_in, use_spatial_in, sample_file, time_config, use_fragments_in, speciation_rates,
			  min_speciation_gen_in, max_speciation_gen_in, 0, 0.0);
	}

	/**
	 * @brief Sets the application arguments for the inputs. Intended for use with the applyspecmodule for
	 * integration with python.
	 *
	 * @param file_in the database to apply speciation rates to
	 * @param use_spatial_in if true, record full spatial data
	 * @param sample_file the sample file to select lineages from the map
	 * @param time_config the time config file to use
	 * @param use_fragments_in fragment file, or "T"/"F" for automatic detection/no detection
	 * @param speciation_rates the speciation rates to apply
	 * @param min_speciation_gen_in the minimum generation rate for speciation in protracted simulations
	 * @param max_speciation_gen_in the maximum generation rate for speciation in protracted simulations
	 * @param metacommunity_size_in
	 * @param metacommunity_speciation_rate_in
	 */
	void setup(string file_in, bool use_spatial_in, string sample_file, string time_config, string use_fragments_in,
			   vector<double> speciation_rates, double min_speciation_gen_in, double max_speciation_gen_in,
			   unsigned long metacommunity_size_in, double metacommunity_speciation_rate_in)
	{
		filename = std::move(file_in);
		use_spatial = use_spatial_in;
		samplemask = std::move(sample_file);
		times_file = std::move(time_config);
		min_speciation_gen = std::move(min_speciation_gen_in);
		max_speciation_gen = std::move(max_speciation_gen_in);
		importTimeConfig();
		use_fragments = !(use_fragments_in == "F");
		fragment_config_file = use_fragments_in;
		bMultiRun = speciation_rates.size() > 1;
		for(double speciation_rate : speciation_rates)
		{
			all_speciation_rates.push_back(speciation_rate);
		}
		metacommunity_size = metacommunity_size_in;
		metacommunity_speciation_rate = metacommunity_speciation_rate_in;
	}

	/**
	 * @brief Import the time config file, if there is one
	 */
	void importTimeConfig()
	{
		if(times_file == "null")
		{
			all_times.push_back(0.0);
		}
		else
		{
			vector<string> tmpimport;
			ConfigOption tmpconfig;
			tmpconfig.setConfig(times_file, false);
			tmpconfig.importConfig(tmpimport);
			for(unsigned int i = 0; i < tmpimport.size(); i++)
			{
				all_times.push_back(stod(tmpimport[i]));
				//					os << "t_i: " << sp.reference_times[i] << endl;
			}
		}
	}
};


#endif //SPECIATIONCOUNTER_SPECSIMPARAMETERS_H
