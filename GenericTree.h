//
// Created by sam on 12/03/2021.
//
#include <string>
#include <vector>
#include <memory>
#include <iostream>

#ifndef NECSIM_GENERICTREE_H
#define NECSIM_GENERICTREE_H

namespace necsim
{

    template<class T>
    class GenericTree
    {
    private:
        std::shared_ptr<T> tree_ptr{};
    public:

        explicit GenericTree() noexcept : tree_ptr(std::make_shared<T>())
        {

        }

        virtual ~GenericTree() = default;

        GenericTree(GenericTree &&other) noexcept : GenericTree()
        {
            *this = std::move(other);
        }

        GenericTree(const GenericTree &other) noexcept: GenericTree()
        {
            this->tree_ptr = other.tree_ptr;
        };

        GenericTree &operator=(GenericTree other) noexcept
        {
            other.swap(*this);
            return *this;
        }

        void swap(GenericTree &other) noexcept
        {
            if(this != &other)
            {
                this->tree_ptr.swap(other.tree_ptr);

            }
        }

        void wipeSimulationVariables()
        {
            tree_ptr->wipeSimulationVariables();
        }

        void importSimulationVariables(std::string config_file)
        {
            tree_ptr->importSimulationVariables(config_file);
        }

        void importSimulationVariablesFromString(std::string config_string)
        {
            tree_ptr->importSimulationVariablesFromString(config_string);
        }

        void setup()
        {
            tree_ptr->setup();
        }

        bool runSimulation()
        {
            return tree_ptr->runSimulation();
        }

        void setResumeParameters(std::string pause_directory_str,
                                 std::string out_directory_str,
                                 long seed,
                                 long task,
                                 long max_time)
        {
            tree_ptr->setResumeParameters(pause_directory_str, out_directory_str, seed, task, max_time);
        }

        void checkSims(std::string pause_directory, long seed, long task)
        {
            tree_ptr->checkSims(pause_directory, seed, task);
        }

        bool hasPaused()
        {
            return tree_ptr->hasPaused();
        }

        void addGillespie(const double &g_threshold)
        {
            tree_ptr->addGillespie(g_threshold);
        }

        void addSpeciationRates(std::vector<long double> spec_rates_long)
        {
            tree_ptr->addSpeciationRates(spec_rates_long);
        }

        void applyMultipleRates()
        {
            tree_ptr->applyMultipleRates();
        }
    };
}
#endif //NECSIM_GENERICTREE_H
