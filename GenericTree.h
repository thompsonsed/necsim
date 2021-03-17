//
// Created by sam on 12/03/2021.
//
#include <string>
#include <vector>
#include <memory>

#ifndef NECSIM_GENERICTREE_H
#define NECSIM_GENERICTREE_H

namespace necsim
{

    template<class T>
    class GenericTree
    {
    private:
        std::unique_ptr<T> tree;
    public:

        GenericTree<T>() : tree(std::make_unique<T>())
        { }

//        GenericTree<T> &operator=(const GenericTree<T> &other) noexcept = default;
//
//        GenericTree<T> &operator=(GenericTree<T> &&other) noexcept = default;
//
        void wipeSimulationVariables()
        {
            tree->wipeSimulationVariables();
        }

        void importSimulationVariables(std::string config_file)
        {
            tree->importSimulationVariables(config_file);
        }

        void importSimulationVariablesFromString(std::string config_string)
        {
            tree->importSimulationVariablesFromString(config_string);
        }

        void setup()
        {
            tree->setup();
        }

        bool runSimulation()
        {
            return tree->runSimulation();
        }

        void setResumeParameters(std::string pause_directory_str,
                                 std::string out_directory_str,
                                 long seed,
                                 long task,
                                 long max_time)
        {
            tree->setResumeParameters(pause_directory_str, out_directory_str, seed, task, max_time);
        }

        void checkSims(std::string pause_directory, long seed, long task)
        {
            tree->checkSims(pause_directory, seed, task);
        }

        bool hasPaused()
        {
            return tree->hasPaused();
        }

        void addGillespie(const double &g_threshold)
        {
            tree->addGillespie(g_threshold);
        }

        void addSpeciationRates(std::vector<long double> spec_rates_long)
        {
            tree->addSpeciationRates(spec_rates_long);
        }

        void applyMultipleRates()
        {
            tree->applyMultipleRates();
        }
    };
}
#endif //NECSIM_GENERICTREE_H
