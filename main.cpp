// This file is part of necsim project which is released under MIT license.
// See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.
/**
 * @file main.cpp
 * @brief A generic simulator for spatially explicit coalescence models suitable for HPC applications.
 * It contains all functions for running large-scale simulations backwards in time using coalescence techniques.
 * Outputs include an SQLite database containing spatial and temporal information about tracked lineages, and allow for
 * rebuilding of the coalescence tree.
 *
 * Currently, a fat-tailed dispersal kernel or normal distribution can be used for dispersal processes.
 *
 * Run with -h to see full input options.
 *
 * Outputs include
 * - habitat map file(s)
 * - species richness and species abundances for the supplied minimum speciation rate.
 * - SQL database containing full spatial data. This can be later analysed by the Speciation_Counter program for
 *   applying higher speciation rates.
 *
 * Contact: samuel.thompson14@imperial.ac.uk or thompsonsed@gmail.com
 *
 * Based heavily on code written by James Rosindell
 *
 * Contact: j.rosindell@imperial.ac.uk
 *
 *
 * @author Samuel Thompson
 *
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 *
*/



#ifdef DEBUG
#ifndef verbose
#define verbose
#endif
#endif

#include "Logger.h"
#include "Logging.h"
#include "ProtractedSpatialTree.h"
#include "SimulationTemplates.h"
#include "GenericTree.h"
// #define historical_mode // not required unless you experience problems.
// This performs a more thorough check after each move operation.
// Currently, it will also check that the historical state value is greater than the returned value within every map cell.
// Note that this may cause problems if the historical state is not the state with the highest number of individuals.



/************************************************************
        MAIN ROUTINE AND COMMAND LINE ARG ROUTINES

 ************************************************************/
namespace necsim
{
    class A
    {
    protected:
        // storing the coalescence tree itself
        shared_ptr<vector<TreeNode>> data;
        // a reference for the last written point in data.
        unsigned long enddata{};
        // Stores the command line current_metacommunity_parameters and parses the required information.
        shared_ptr<SimParameters> sim_parameters{};
        // random number generator
        shared_ptr<RNGController> NR{};
        // Storing the speciation rates for later reference.
        vector<long double> speciation_rates{};
        // flag for having set the simulation seed.
        bool seeded{};
        // random seed
        long long seed{};
        // for file naming - good to know which task in a series is being executed here
        long long task{};
        // The map file containing the times that we want to expand the model and record all lineages again.
        // If this is null, uses_temporal_sampling will be false and the vector will be empty.
        string times_file{};
        vector<double> reference_times{};
        // Set to true if we are recording at times other than the present day.
        bool uses_temporal_sampling{};
        // The time variables (for timing the simulation in real time)
        time_t start{}, sim_start{}, sim_end{}, now{}, sim_finish{}, out_finish{};
        time_t time_taken{};
        // Active lineages stored as a row of datapoints
        vector<DataPoint> active{};
        // Stores the point of the end of the active vector. 0 is reserved as null
        unsigned long endactive{};
        // the maximum size of endactive
        unsigned long startendactive{};
        // the maximum simulated number of individuals in the present day.
        unsigned long maxsimsize{};
        // for create the link to the speciationcounter object which handles everything.
        Community community{};
        // This might need to be updated for simulations that have large changes in maximum population size over time.
        // number of simulation num_steps
        long steps{};
        // Maximum time to run for (in seconds)
        unsigned long maxtime{};
        // number of generations passed,
        double generation{};
        // The number of individuals per cell
        double deme{};
        // The proportion of individuals to sample
        double deme_sample{};
        // the speciation rate
        long double spec{};
        // Path to output directory
        string out_directory{};
        // sqlite3 object that stores all the data
        shared_ptr<SQLiteHandler> database{};
        // only set to true if the simulation has finished, otherwise will be false.
        bool sim_complete{};
        // set to true when variables are imported
        bool has_imported_vars{};
        // If sql database is written first to memory, then need another object to contain the in-memory database.
        // Create the step object that will be retained for the whole simulation.
        // Does not need saving on simulation pause.
        Step this_step{};
        string sql_output_database{};
        // If true, means the command-line imports were under the (deprecated) fullmode.
        bool bFullMode{};
        // If true, the simulation is to be resumed.
        bool bResume{};
        // If true, a config file contains the simulation variables.
        bool bConfig{};
        // If true, simulation can be resumed.
        bool has_paused{}, has_imported_pause{};
        // Should always be false in the base class
        bool bIsProtracted{};
        // variable for storing the paused sim location if files have been moved during paused/resumed simulations!
        string pause_sim_directory{};
        // Set to true to use the gillespie method - this is currently only supported for spatial simulations using a
        // dispersal map and point speciation (i.e. the method is unsupported for non-spatial simulations,
        // spatial simulations not using a dispersal map and those that use protracted speciation).
        bool using_gillespie{};
#ifdef sql_ram
        SQLiteHandler outdatabase{};
#endif

    public:
        A() : data(make_shared<vector<TreeNode >>()), enddata(0), sim_parameters(make_shared<SimParameters>()),
              NR(make_shared<RNGController>()), speciation_rates(), seeded(false), seed(-1), task(-1),
              times_file("null"), reference_times(), uses_temporal_sampling(false), start(0), sim_start(0), sim_end(0),
              now(0), sim_finish(0), out_finish(0), time_taken(0), active(), endactive(0), startendactive(0),
              maxsimsize(0), community(data), steps(0), maxtime(0), generation(0.0), deme(0.0), deme_sample(0.0),
              spec(0.0), out_directory(""), database(make_shared<SQLiteHandler>()), sim_complete(false),
              has_imported_vars(false), this_step(), sql_output_database("null"), bFullMode(false), bResume(false),
              bConfig(true), has_paused(false),
#ifdef sql_ram
              outdatabase(),
#endif
              has_imported_pause(false), bIsProtracted(false), pause_sim_directory("null"), using_gillespie(false)
        {
            sql_output_database = "null"; // TODO remove?
            std::cout << "Making A at " << this << std::endl;
        }

        virtual ~A() = default;

        A(A &&other) noexcept : A()
        {
            *this = std::move(other);
        }

        A(const A &other) : A()
        {
            *this = other;
        };

        A &operator=(A other) noexcept
        {
            other.swap(*this);
            return *this;
        };

        void swap(A &other) noexcept
        {
            if(this != &other)
            {
                std::cout << "this: " << this << std::endl; // TODO remove
                std::cout << "other: " << &other << std::endl;
                sql_output_database = "test";
                std::cout << "this output: " << this->sql_output_database << std::endl;
                std::cout << "other output: " << other.sql_output_database << std::endl;

                std::swap(data, other.data);
                std::swap(enddata, other.enddata);
                std::swap(sim_parameters, other.sim_parameters);
                std::swap(NR, other.NR);
                std::swap(speciation_rates, other.speciation_rates);
                std::swap(seeded, other.seeded);
                std::swap(seed, other.seed);
                std::swap(task, other.task);
                std::swap(times_file, other.times_file);
                std::swap(reference_times, other.reference_times);
                std::swap(uses_temporal_sampling, other.uses_temporal_sampling);
                std::swap(start, other.start);
                std::swap(sim_start, other.sim_start);
                std::swap(sim_end, other.sim_end);
                std::swap(now, other.now);
                std::swap(sim_finish, other.sim_finish);
                std::swap(out_finish, other.out_finish);
                std::swap(time_taken, other.time_taken);
                std::swap(active, other.active);
                std::swap(endactive, other.endactive);
                std::swap(startendactive, other.startendactive);
                std::swap(maxsimsize, other.maxsimsize);
                std::swap(community, other.community);
                std::swap(steps, other.steps);
#ifdef sql_ram
                std::swap(outdatabase, other.outdatabase);
#endif
                std::swap(maxtime, other.maxtime);
                std::swap(generation, other.generation);
                std::swap(deme, other.deme);
                std::swap(deme_sample, other.deme_sample);
                std::swap(spec, other.spec);
                std::swap(out_directory, other.out_directory);
                std::swap(database, other.database);
                std::swap(sim_complete, other.sim_complete);
                std::swap(has_imported_vars, other.has_imported_vars);
                std::swap(this_step, other.this_step);
                std::swap(sql_output_database, other.sql_output_database);
                std::swap(bFullMode, other.bFullMode);
                std::swap(bResume, other.bResume);
                std::swap(bConfig, other.bConfig);
                std::swap(has_paused, other.has_paused);
                std::swap(has_imported_pause, other.has_imported_pause);
                std::swap(bIsProtracted, other.bIsProtracted);
                std::swap(pause_sim_directory, other.pause_sim_directory);
                std::swap(using_gillespie, other.using_gillespie);
            }
        }

    };

    class B : public virtual A
    {
    protected:
        std::string another_test{};
    public:
        B() : A(), another_test()
        { }
    };

}

/**
 * @brief Main function containing program structure
 * @param argc the number of command-line arguments provided
 * @param argv a pointer to the arguments
 * @return a program exit code, 0 if successful, -1 (generally) indicates an error.
 */
int main(int argc, char* argv[])
{
//        logger = new Logger();
//        vector<string> comargs;
//        importArgs(static_cast<const unsigned int &>(argc), argv, comargs);
//        const string &config_file = getConfigFileFromCmdArgs(comargs);
//        runMain<SpatialTree>(config_file);
//        delete logger;
//        return 0;
//    necsim::ProtractedSpatialTree test2;
//    test3 = necsim::ProtractedSpatialTree();
//    necsim::ProtractedSpatialTree test = necsim::ProtractedSpatialTree();
//    test2 = test3;
//    test2 = std::move(test);
//    test1 = test3;
    necsim::Tree tree;
    necsim::logger = new necsim::Logger();
    int __pyx_lineno = 0;
    const char* __pyx_filename = NULL;
    int __pyx_clineno = 0;
    try
    {
        tree = necsim::Tree();
    }
    catch(...)
    {
        return 1;
    }

    return 0;
}