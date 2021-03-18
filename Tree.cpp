// This file is part of necsim project which is released under MIT license.
// See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.

/**
 * @file Tree.cpp
 * @brief Contains the main simulation object for spatially implicit coalescence simulations.
 * Provides the basis for spatially explicit versions in SpatialTree, and protracted speciation versions in
 * ProtractedTree and ProtractedSpatialTree.
 *
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 */

#include <algorithm>

#include "cpp17_includes.h"
#include "Tree.h"
#include "Logging.h"
#include "LogFile.h"

namespace necsim
{


    void Tree::importSimulationVariables(string configfile)
    {
        sim_parameters->importParameters(configfile);
        runFileChecks();
    }

    void Tree::importSimulationVariables(ConfigParser config)
    {
        sim_parameters->importParameters(config);
        runFileChecks();
    }

    void Tree::importSimulationVariablesFromString(string config_string)
    {
        std::stringstream ss;
        ss << config_string;
        std::istream &istream1 = ss;
        ConfigParser config;
        config.parseConfig(istream1);
        sim_parameters->importParameters(config);
        runFileChecks();
    }

    void Tree::runFileChecks()
    {
        checkOutputDirectory();
        // Now check for paused simulations
        checkSims();
    }

    void Tree::wipeSimulationVariables()
    {
        sim_parameters = make_shared<SimParameters>();
    }

    void Tree::internalSetup(shared_ptr<SimParameters> sim_parameters_in)
    {
        sim_parameters = std::move(sim_parameters_in);
        setup();
    }

    bool Tree::checkOutputDirectory()
    {
        if(sim_parameters->output_directory != "null")
        {
            try
            {
                doesExist(sim_parameters->output_directory);
            }
            catch(std::runtime_error &re)
            {
                writeInfo("Output folder does not exist... creating...");
                bool bOutputFolder = fs::create_directory(sim_parameters->output_directory);
                if(bOutputFolder)
                {
                    writeInfo("done.\n");
                }
                else
                {
                    writeError(re.what());
                }
            }
        }
        else
        {
            throw FatalException("ERROR_MAIN_009: FATAL. Output folder cannot be null.");
        }
        return true;
    }

    void Tree::checkSims()
    {
        checkSims(sim_parameters->output_directory, sim_parameters->seed, sim_parameters->task);
    }

    void Tree::checkSims(string output_dir, long seed_in, long task)
    {

        std::stringstream os;
        os << "Checking for unfinished simulations..." << std::flush;
        std::ifstream out;
        string file_to_open;
        //	char file_to_open[100];
        //	sprintf (file_to_open, "%s/Pause/Data_%i.csv",outdirect,int(task));
        file_to_open = output_dir + string("/Pause/Dump_main_") + std::to_string((unsigned long long) task) + "_"
                       + std::to_string((unsigned long long) seed_in) + string(".csv");
        out.open(file_to_open);
        if(out.good())
        {
            os << "done." << std::endl << "File found containing unfinished simulations." << std::endl;
            writeInfo(os.str());
            if(!has_imported_pause)
            {
                setResumeParameters(sim_parameters->output_directory,
                                    sim_parameters->output_directory,
                                    static_cast<unsigned long>(sim_parameters->seed),
                                    static_cast<unsigned long>(sim_parameters->task),
                                    sim_parameters->max_time);
            }
            has_paused = true;
        }
        else
        {
            os << "done." << std::endl << "No files found containing unfinished simulations." << std::endl;
            writeInfo(os.str());
            has_paused = false;
        }
    }

    void Tree::setParameters()
    {
        if(!has_imported_vars)
        {
            out_directory = sim_parameters->output_directory;

            task = sim_parameters->task;
            seed = sim_parameters->seed;

            deme = sim_parameters->deme;
            deme_sample = sim_parameters->deme_sample;
            spec = sim_parameters->spec;
            maxtime = sim_parameters->max_time;
            times_file = sim_parameters->times_file;
            setProtractedVariables(sim_parameters->min_speciation_gen, sim_parameters->max_speciation_gen);
            has_imported_vars = true;
        }
        else
        {
            throw FatalException("Setup has already been started: variables already imported.");
        }
    }

    void Tree::setProtractedVariables(double speciation_gen_min, double speciation_gen_max)
    { }

    bool Tree::hasPaused()
    {
        return has_paused;
    }

    vector<double> Tree::getTemporalSampling()
    {
        if(uses_temporal_sampling)
        {
            return reference_times;
        }
        else
        {
            vector<double> tmp;
            tmp.push_back(0.0);
            return tmp;
        }
    }

    long long Tree::getSeed()
    {
        return seed;
    }

    long long Tree::getJobType()
    {
        return task;
    }

    void Tree::setSeed(long long seed_in)
    {
        if(!seeded)
        {
            // There are problems if the random seed is ever 0 as it will produce an identical output to if the seed is 1
            // therefore the user should be informed (but an error is not thrown).
            if(seed_in == 0)
            {
                std::stringstream ss;
                ss << "Seed is set as 0 - this will produce identical behaviour to if the seed is 1." << std::endl;
                writeCritical(ss.str());
            }
            NR->setSeed(seed_in);
            seed = seed_in;
            seeded = true;
        }
    }

    double Tree::getGeneration() const
    {
        return generation;
    }

    unsigned long Tree::getInitialCount()
    {
        return static_cast<unsigned long>(floor(deme * deme_sample));
    }

    unsigned long Tree::setObjectSizes()
    {
        unsigned long initial_count = getInitialCount();
        active.resize(initial_count + 1);
        data->resize(2 * initial_count + 1);
        return initial_count;
    }

    void Tree::setup()
    {
        printSetup();
        if(has_imported_pause)
        {
            setResumeParameters();
            simResume();
        }
        else
        {
            // Start the timer
            time(&start);
            setParameters();
            setInitialValues();
            generateObjects();
        }
    }

    void Tree::setInitialValues()
    {
        // other variables
        steps = 0;
        generation = 0;
        // Set the seed
        setSeed(seed);
        setTimes();
        sim_parameters->printVars();
        // Determine the speciation rates which will be applied after the simulation completes.
        determineSpeciationRates();
    }

    void Tree::setSimStartVariables()
    {
        this_step.bContinueSim = true;
        this_step.time_reference = 0;
        if(uses_temporal_sampling && generation > 0.0)
        {
            for(unsigned int i = 0; i < reference_times.size(); i++)
            {
                if(reference_times[i] > generation)
                {
                    this_step.time_reference = i + 1;
                    break;
                }
            }
        }
    }

    void Tree::printSetup()
    {
        std::stringstream os;
        os << "*************************************************" << std::endl;
        os << "Setting up simulation..." << std::endl;
        writeInfo(os.str());
        os.str("");
        time(&start);
    }

    void Tree::setTimes()
    {
        // Import the time sample points
        if(!reference_times.empty())
        {
            throw FatalException("Reference times have already been set.");
        }
        if(times_file == "set")
        {
            uses_temporal_sampling = true;
            reference_times = sim_parameters->times;
            sort(reference_times.begin(), reference_times.end());
        }
        if(reference_times.size() <= 1)
        {
            times_file = "null";
            reference_times.clear();
            reference_times.push_back(0.0);
        }
    }

    void Tree::determineSpeciationRates()
    {
        if(bConfig)
        {
            if(sim_parameters->configs.hasSection("spec_rates"))
            {
                vector<string> spec_rates = sim_parameters->configs.getSectionValues("spec_rates");
                for(const auto &spec_rate : spec_rates)
                {
                    speciation_rates.push_back(stod(spec_rate));
                }
            }
        }
        else
        {
            speciation_rates.push_back(spec);
        }
        sort(speciation_rates.begin(), speciation_rates.end());
    }

    void Tree::addSpeciationRates(vector<long double> spec_rates_in)
    {
        if(speciation_rates.empty())
        {
            speciation_rates.push_back(spec);
        }
        for(const auto &item : spec_rates_in)
        {
            if(item > spec)
            {
                speciation_rates.push_back(item);
            }
            else if(doubleCompare(spec, item, item * 0.000001))
            {
                speciation_rates.push_back(spec);
            }
            else
            {
                std::stringstream ss;
                ss << "Speciation rate of " << item << " is less than the minimum possible (" << spec;
                ss << ") - skipping." << std::endl;
                throw FatalException(ss.str());
            }
        }
        // Sort the speciation rates remove duplicates
        sort(speciation_rates.begin(), speciation_rates.end());
        speciation_rates.erase(unique(speciation_rates.begin(), speciation_rates.end()), speciation_rates.end());
    }

    void Tree::generateObjects()
    {
        unsigned long initial_count = setObjectSizes();
        endactive = 0;
        unsigned long number_start = fillObjects(initial_count);
        std::stringstream os;
        os << "\rSetting up simulation...done.                           " << std::endl;
        os << "Number of individuals simulating: " << endactive << std::endl;
        writeInfo(os.str());
        maxsimsize = enddata;
        if(active.size() < endactive || endactive == 0)
        {

            if(endactive == 0)
            {
                throw std::runtime_error("No individuals to simulate! Check set up. Exiting...");
            }
            else
            {
                std::stringstream ss;
                ss << "ERROR_MAIN_007: FATAL. Sizing error - endactive is greater than the size of active. ";
                ss << "Please report this bug" << std::endl;
                ss << "endactive: " << endactive << std::endl;
                ss << "active.size: " << active.size() << std::endl;
                ss << "initial_count: " << initial_count << std::endl;
                ss << "number_start: " << number_start << std::endl;
                throw FatalException(ss.str());
            }
        }
        startendactive = endactive;
    }

    unsigned long Tree::fillObjects(const unsigned long &initial_count)
    {
        active[0].setup(0, 0, 0, 0, 0, 0, 0);
        unsigned long number_start = 0;
        std::stringstream os;
        os << "\rSetting up simulation...filling grid                           " << std::flush;
        writeInfo(os.str());
        // This loop adds individuals to data and active (for storing the coalescence tree and active lineage tracking)
        double sample_number = floor(deme_sample * deme);
        for(unsigned long i = 0; i < sample_number; i++)
        {
            number_start++;
            // Add the species to active
            active[number_start].setup(number_start, i, 1);
            // Add a tip in the TreeNode for calculation of the coalescence tree at the
            // end of the simulation.
            // This also contains the start x and y position of the species.
            (*data)[number_start].setup(true);
            (*data)[number_start].setSpec(NR->d01());
            endactive++;
            enddata++;
        }
        if(number_start == initial_count)  // Check that the two counting methods match up.
        {
        }
        else
        {
            if(initial_count > 1.1 * number_start)
            {
                writeWarning("Data usage higher than neccessary - check allocation of individuals to the grid.");
                std::stringstream ss;
                ss << "Initial count: " << initial_count << "  Number counted: " << number_start << std::endl;
                writeWarning(ss.str());
            }
        }
#ifdef DEBUG
        validateLineages();
#endif
        return number_start;
    }

    void Tree::runSingleLoop()
    {
        chooseRandomLineage();
        writeStepToConsole();
        // See estSpecnum for removed code.
        // Check that we still want to continue the simulation.
        if(this_step.bContinueSim)
        {
            auto chosen_reference = active[this_step.chosen].getReference();
            // increase the counter of the number of moves (or generations) the lineage has undergone.
            (*data)[chosen_reference].increaseGen();
            // Check if speciation happens
            if(calcSpeciation((*data)[chosen_reference].getSpecRate(),
                              0.99999 * spec,
                              (*data)[chosen_reference].getGenerationRate()))
            {
                speciation(this_step.chosen);
            }
            else
            {
                // remove the species data from the species lineage_indices to be placed somewhere new.
                removeOldPosition(this_step.chosen);
                calcNextStep();
#ifdef DEBUG
                debugCoalescence();
#endif
                if(this_step.coal)
                {
                    coalescenceEvent(this_step.chosen, this_step.coalchosen);
                }
            }
        }

#ifdef DEBUG
        debugEndStep();
#endif
        if(uses_temporal_sampling && endactive == 1)
        {
            // Check whether we need to continue simulating at a later time.
            if(reference_times[this_step.time_reference] > generation)
            {
                // Then we need to expand the map
                // This is a hack, I know it's a hack and is wrong, and I aint gonna change it :)
                (*data)[active[endactive].getReference()].setSpec(0.0);
                // First speciate the remaining lineage
                speciation(endactive);
                generation = reference_times[this_step.time_reference] + 0.000000000001;
                checkTimeUpdate();
                if(endactive < 2)
                {
                    this_step.bContinueSim = false;
                }
            }
            // TODO fix this to account for potential speciation of the remaining lineage!
        }
    }

    bool Tree::runSimulation()
    {

        writeSimStartToConsole();
        // Main while loop to process while there is still time left and the simulation is not complete.
        // Ensure that the step object contains no data->
        this_step.wipeData();
        setSimStartVariables();
        if(endactive < 2)
        {
            return stopSimulation();
        }
        if(using_gillespie)
        {
            return runSimulationGillespie();
        }
        return runSimulationNoGillespie();
    }

    bool Tree::runSimulationNoGillespie()
    {
        // Create the move object
        do
        {
            runSingleLoop();
        }
        while((endactive > 1) && (steps < 100 || difftime(sim_end, start) < maxtime) && this_step.bContinueSim);
        // If the simulations finish correctly, output the completed data->
        // Otherwise, pause the simulation and save objects to file.
        return stopSimulation();
    }

    bool Tree::stopSimulation()
    {
        if(endactive > 1)
        {
            std::stringstream os;
            time(&sim_finish);
            time_taken += sim_finish - start;
            os.str("");
            os << "........out of time!" << std::endl;
            os << "Pausing simulation: add extra time or re-run to ensure simulation completion." << "\n";
            os << "Lineages remaining: " << endactive << "\n";
            writeInfo(os.str());
            simPause();
            return false;
        }
        else
        {
            for(unsigned int i = 0; i <= endactive; i++)
            {
                speciateLineage(active[i].getReference());
                (*data)[active[i].getReference()].setSpec(0.0);
            }
            sim_complete = true;
            time(&sim_finish);
            time_taken += sim_finish - start;
            if(!this_step.bContinueSim)
            {
                writeInfo("done - desired number of species achieved!\n");
                return true;
            }
            else
            {
                writeInfo("done.\n");
                return true;
            }
        }
    }

    void Tree::writeSimStartToConsole()
    {
        // now do the calculations required to build the tree
        std::stringstream os;
        os << "*************************************************" << std::endl;
        os << "Beginning simulations..." << std::flush;
        writeInfo(os.str());
        os.str("");

        //		double current_gen =0;
        // check time
        time(&sim_start);
        time(&sim_end);
        time(&now);
    }

    void Tree::writeStepToConsole()
    {
        if(steps % 10000 == 0)
        {
            time(&sim_end);
#ifdef verbose
            if(sim_end - now > 0.2)  // output every 0.2 seconds
            {
                double dPercentComplete = 20 * (1 - (double(endactive) / double(startendactive)));
                time(&now);
                if(this_step.number_printed < dPercentComplete)
                {
                    std::stringstream os;
                    os << "\rBeginning simulations...";
                    this_step.number_printed = 0;
                    while(this_step.number_printed < dPercentComplete)
                    {
                        os << ".";

                        this_step.number_printed++;
                    }
                    os << std::flush;
                    writeInfo(os.str());
                }
            }
#endif // verbose
        }
    }

    void Tree::incrementGeneration()
    {
        steps++;
        // increment generation counter
        generation += 2.0 / (double(endactive));
    }

    void Tree::chooseRandomLineage()
    {
        incrementGeneration();
        // choose a random lineage to die and be reborn out of those currently active
        this_step.chosen = NR->i0(endactive - 1) + 1;  // cannot be 0
        // Rejection sample based on reproductive potential
        updateStepCoalescenceVariables();
    }

    void Tree::updateStepCoalescenceVariables()
    {
        this_step.coalchosen = 0;
        this_step.coal = false;
    }

    void Tree::speciation(const unsigned long &chosen)
    {
        // alter the data such that it reflects the speciation event.
        const unsigned long data_position = active[chosen].getReference();
#ifdef DEBUG
        // Store debug information in DEBUG mode
        if((*data)[data_position].hasSpeciated())
        {
            std::stringstream ss;
            ss << "Chosen: " << chosen << std::endl;
            writeLog(50, ss);
            ss.str("");
            ss << "Endactive: " << endactive << std::endl;
            writeLog(50, ss);
            (*data)[data_position].logLineageInformation(50);
            active[chosen].logActive(50);
            throw FatalException("ERROR_MOVE_028: Attempting to speciate a speciated species.");
        }
#endif
        speciateLineage(data_position);
        // Now remove the old chosen lineage from the active directory.
        removeOldPosition(chosen);
        switchPositions(chosen);
    }

    void Tree::speciateLineage(const unsigned long &data_position)
    {
        (*data)[data_position].speciate();
    }

    void Tree::removeOldPosition(const unsigned long &chosen)
    {
        // This may seem a bit stupid, but this function is overwridden with more complex routines in child classes.
        active[chosen].setListPosition(0);
    }

    void Tree::switchPositions(const unsigned long &chosen)
    {
#ifdef DEBUG

        if(chosen > endactive)
        {
            std::stringstream ss;
            ss << "chosen: " << chosen << " endactive: " << endactive << std::endl;
            writeLog(50, ss);
            throw FatalException("ERROR_MOVE_023: Chosen is greater than endactive. Check move function.");
        }
#endif // DEBUG
        if(chosen != endactive)
        {
            std::swap(active[chosen], active[endactive]);
        }
        endactive--;

    }

    void Tree::calcNextStep()
    {
        unsigned long random_lineage = NR->i0(static_cast<unsigned long>(deme)) + 1;
        if(random_lineage != this_step.chosen && random_lineage <= endactive)
        {
            // then we have a coalescence event
            this_step.coal = true;
            this_step.coalchosen = random_lineage;
        }
    }

    bool Tree::calcSpeciation(const long double &random_number,
                              const long double &speciation_rate,
                              const unsigned long &no_generations)
    {
        return checkSpeciation(random_number, speciation_rate, no_generations);
    }

    void Tree::coalescenceEvent(const unsigned long &chosen, unsigned long &coalchosen)
    {
        // coalescence occured, so we need to adjust the data appropriatedly
        // our chosen lineage has merged with the coalchosen lineage, so we need to sync up the data->
        enddata++;
        (*data)[enddata].setup(0,
                               active[chosen].getXpos(),
                               active[chosen].getYpos(),
                               active[chosen].getXwrap(),
                               active[chosen].getYwrap(),
                               generation);

        // First perform the move
        (*data)[active[chosen].getReference()].setParent(enddata);
        (*data)[active[coalchosen].getReference()].setParent(enddata);
        active[coalchosen].setMinmax(std::max(active[coalchosen].getMinmax(),
                                         active[chosen].getMinmax()));  // set the new minmax to the maximum of the two minimums.
        active[chosen].setMinmax(active[coalchosen].getMinmax());
        (*data)[enddata].setGenerationRate(0);
        (*data)[enddata].setSpec(NR->d01());
        active[chosen].setReference(enddata);
        active[coalchosen].setReference(enddata);
        //		removeOldPosition(chosen);
        switchPositions(chosen);
    }

    void Tree::checkTimeUpdate()
    {
        if(uses_temporal_sampling && this_step.time_reference < reference_times.size())
        {
            // check if we need to update
            if(reference_times[this_step.time_reference] <= generation)
            {
                //					os << "check2" << std::endl;
                if(reference_times[this_step.time_reference] > 0.0)
                {
                    std::stringstream os;
                    os << "\n" << "expanding map at generation " << generation << std::endl;
                    addLineages(reference_times[this_step.time_reference]);
                    writeInfo(os.str());
                }
                this_step.time_reference++;
            }
        }
    }

    void Tree::addLineages(double generation_in)
    {
        auto number_added = static_cast<unsigned long>(floor(deme_sample * deme));
        // Store all the data lineages to add in a vector
        vector<TreeNode> data_to_add{};
        // change those that already exist to tips
        for(unsigned long i = 0; i < endactive; i++)
        {
            // With probability deme_sample, just change the active lineage to a tip.
            if(checkProportionAdded(deme_sample) && number_added > 0)
            {
                number_added--;
                makeTip(endactive, generation_in, data_to_add);
            }
        }
        checkSimSize(data_to_add.size() + number_added, number_added);
        for(auto &item : data_to_add)
        {
            enddata++;
            (*data)[enddata] = item;
        }
        for(unsigned long i = 0; i < number_added; i++)
        {
            enddata++;
            endactive++;
            active[endactive].setup(enddata, endactive, 1.0);
            (*data)[enddata].setup(true, 0, 0, 0, 0, generation_in);
            (*data)[enddata].setSpec(NR->d01());
        }
    }

    bool Tree::checkProportionAdded(const double &proportion_added)
    {
        return NR->d01() < proportion_added;
    }

    void Tree::checkSimSize(unsigned long req_data, unsigned long req_active)
    {
        unsigned long min_active = endactive + req_active + 2;
        unsigned long min_data = enddata + req_data + 2;
        // Take into account future coalescence events
        min_data += min_active * 2;
        if(data->size() < min_data)
        {
            // change the size of data
            data->resize(min_data);
        }

        if(active.size() < min_active)
        {
            // change the size of active.
            active.resize(min_active);
        }
    }

    void Tree::makeTip(const unsigned long &tmp_active, const double &generationin, vector<TreeNode> &data_added)
    {
        auto cur_active = &active[tmp_active];
        auto cur_data = &(*data)[cur_active->getReference()];
        if(cur_data->isTip())
        {
            createNewTip(tmp_active, generationin, data_added);
        }
        else
        {
            cur_data->setGeneration(generationin);
            cur_data->setTip(true);
            cur_data->setPosition(cur_active->getXpos(),
                                  cur_active->getYpos(),
                                  cur_active->getXwrap(),
                                  cur_active->getYwrap());
        }
    }

    void Tree::createNewTip(unsigned long i, double generationin, vector<TreeNode> &data_added)
    {
        auto cur_active = &active[i];
        TreeNode tmp_tree_node;
        tmp_tree_node.setup(true,
                            cur_active->getXpos(),
                            cur_active->getYpos(),
                            cur_active->getXwrap(),
                            cur_active->getYwrap(),
                            generationin);
        // Now link the old tip to the new tip
        auto data_pos = enddata + data_added.size() + 1;
        (*data)[cur_active->getReference()].setParent(data_pos);
        tmp_tree_node.setGenerationRate(0);
        tmp_tree_node.setSpec(NR->d01());
        cur_active->setReference(data_pos);
        data_added.emplace_back(tmp_tree_node);
    }

    void Tree::applySpecRate(long double sr, double t)
    {
        setupCommunityCalculation(sr, t);
        community.createDatabase();
#ifdef record_space
        community.recordSpatial();
#endif
    }

    void Tree::applySpecRateInternal(long double sr, double t)
    {
        setupCommunityCalculation(sr, t);
        community.calculateCoalescenceTree();
        community.calcSpeciesAbundance();
    }

    shared_ptr<vector<unsigned long>> Tree::getCumulativeAbundances()
    {
        return community.getCumulativeAbundances();
    }

    shared_ptr<std::map<unsigned long, unsigned long>> Tree::getSpeciesAbundances(const unsigned long &community_reference)
    {
        return community.getSpeciesAbundances(community_reference);
    }

    shared_ptr<vector<unsigned long>> Tree::getSpeciesAbundances()
    {
        return community.getSpeciesAbundances();
    }

    ProtractedSpeciationParameters Tree::setupCommunity()
    {
        return community.setupInternal(sim_parameters, database);
    }

    void Tree::setupCommunityCalculation(long double sr, double t)
    {
        auto tmp = setupCommunity();
        MetacommunityParameters null_parameters;
        community.addCalculationPerformed(sr, t, false, null_parameters, tmp);
    }

    void Tree::applySpecRate(long double sr)
    {
        applySpecRate(sr, 0.0);
    }

    void Tree::applyMultipleRates()
    {
        if(!sim_complete)
        {
            throw FatalException("Simulation is not complete - cannot apply speciation rates.");
        }
        std::stringstream os;
        if(speciation_rates.empty())
        {
            os << "No additional speciation rates to apply." << std::endl;
        }
        speciation_rates.push_back(spec);
        // Get only unique speciation rates
        vector<long double> unique_speciation_rates;
        for(const long double &s : speciation_rates)
        {
            bool add = true;
            for(const long double &u : unique_speciation_rates)
            {
                if(doubleCompare(u, s, s * 0.00001))
                {
                    add = false;
                }
            }
            if(add)
            {
                unique_speciation_rates.push_back(s);
            }
        }
        speciation_rates = unique_speciation_rates;
        os << "Speciation rate";
        if(speciation_rates.size() > 1)
        {
            os << "s are: ";
        }
        else
        {
            os << " is: ";
        }
        for(unsigned long i = 0; i < speciation_rates.size(); i++)
        {
            os << speciation_rates[i];
            if(i + 1 == speciation_rates.size())
            {
                os << "." << std::endl;
            }
            else
            {
                os << ", ";
            }
        }
        // Now check to make sure repeat speciation rates aren't done twice (this is done to avoid the huge number of errors
        // SQL throws if you try to add identical data
        sortData();
        sqlCreate();
        vector<double> temp_sampling = getTemporalSampling();
        os << "Time";
        if(temp_sampling.size() > 1)
        {
            os << "s are: ";
        }
        else
        {
            os << " is: ";
        }
        for(unsigned long i = 0; i < temp_sampling.size(); i++)
        {
            os << temp_sampling[i];
            if(i + 1 == temp_sampling.size())
            {
                os << "." << std::endl;
            }
            else
            {
                os << ", ";
            }
        }
        writeInfo(os.str());
        for(const long double &i: speciation_rates)
        {
            for(double k : temp_sampling)
            {
                if(i > spec)
                {
                    applySpecRate(i, k);
                }
                else if(i == spec)
                {
                    // Use the run spec if the rates are very close to equal
                    applySpecRate(spec, k);
                }
            }
        }
        community.writeNewCommunityParameters();
        outputData();
    }

    bool Tree::getProtracted()
    {
        return false;
    }

    string Tree::getProtractedVariables()
    {
        std::stringstream ss;
        ss << "0.0\n0.0\n";
        return ss.str();
    }

    double Tree::getProtractedGenerationMin()
    {
        return 0.0;
    }

    double Tree::getProtractedGenerationMax()
    {
        return 0.0;
    }

    void Tree::sqlOutput()
    {
#ifdef sql_ram
        // open connection to the database file
        remove(sql_output_database.c_str());
        std::stringstream os;
        os << "\tWriting to " << sql_output_database << "..." << std::endl;
        writeInfo(os.str());
        outdatabase.open(sql_output_database);
        // create the backup object to write data to the file from memory.
        outdatabase.backupFrom(*database);
#endif
    }

    void Tree::createAndOutputData()
    {
        sortData();
        sqlCreate();
        // Run the data sorting functions and output the data into the correct format.
        outputData();
    }

    void Tree::outputData()
    {
        time(&out_finish);
#ifdef sql_ram
        sqlOutput();
#endif
        time(&sim_end);
        writeTimes();
    }

    void Tree::sortData()
    {
        // Sort and process the species lineage_indices so that the useful information can be extracted from it.
        std::stringstream os;
        os << "Finalising data..." << std::flush;
        writeInfo(os.str());
        os.str("");
        // coalescence finished - process speciation
        // check the data structure
        if(enddata > data->size())
        {
#ifdef DEBUG
            std::stringstream ss;
            ss << "enddata: " << enddata << std::endl;
            ss << "data->size(): " << data->size() << std::endl;
            writeLog(50, ss);
#endif // DEBUG
            throw FatalException("Enddata greater than data size. Programming error likely.");
        }
        // Now make sure those left in endactive will definitely speciate.
        for(unsigned long i = 1; i <= endactive; i++)
        {
            (*data)[active[i].getReference()].setSpec(0.0);
        }
        // Double check speciation events have been counted.
        unsigned long spec_up_to = 0;
        for(unsigned int i = 1; i <= enddata; i++)
        {
            if(calcSpeciation((*data)[i].getSpecRate(), spec, (*data)[i].getGenerationRate()))
            {
                spec_up_to++;
                (*data)[i].speciate();
            }
        }
        try
        {
            for(unsigned long i = 1; i <= enddata; i++)
            {
                if((!((*data)[i].hasSpeciated())) && ((*data)[i].getParent() == 0 && (*data)[i].exists()))
                {
                    throw FatalException(string(std::to_string((long long) i) + " has not speciated and parent is 0."));
                }
            }
            // here we check the data is valid - alternative validity check.
            for(unsigned long i = 1; i <= enddata; i++)
            {
                if(!((*data)[i].hasSpeciated()) && (*data)[i].exists())
                {
                    long j = i;
                    while(!((*data)[j].hasSpeciated()))
                    {
                        j = (*data)[j].getParent();
                        if(j == 0)
                        {
                            throw FatalException("0 found in parent while following speciation trail.");
                        }
                    }
                }
            }
        }
        catch(FatalException &me)
        {
#ifdef DEBUG
            writeLog(30, me.what());
            writeLog(30, "Returning max possible size (may cause RAM issues).");
#endif // DEBUG
            std::stringstream ss;
            ss << "\nError found when validating coalescence tree post-simulation: " << me.what() << std::endl;
            writeCritical(ss.str());
        }
        writeInfo("done.\n");
    }

    void Tree::writeTimes()
    {
        std::stringstream os;
        os << "Total generations simulated (steps): " << generation << " (" << steps << ")" << std::endl;
        os << "Setup time was " << floor((sim_start - start) / 60) << " minutes " << (sim_start - start) % 60
           << " seconds" << std::endl;
        os << "Simulation time was " << floor((sim_finish - sim_start) / 3600) << " hours "
           << (floor((sim_finish - sim_start) / 60) - 60 * floor((sim_finish - sim_start) / 3600)) << " minutes "
           << (sim_finish - sim_start) % 60 << " seconds" << std::endl;
        os << "File output and species calculation time was " << floor((out_finish - sim_finish) / 60) << " minutes "
           << (out_finish - sim_finish) % 60 << " seconds" << std::endl;
        os << "SQL output time was " << floor((sim_end - out_finish) / 60) << " minutes " << (sim_end - out_finish) % 60
           << " seconds" << std::endl;
        time_taken += (sim_end - sim_finish);
        os << "Total simulation and output time was " << floor((time_taken) / 3600) << " hours " << std::flush;
        os << (floor((time_taken) / 60) - 60 * floor((time_taken) / 3600)) << std::flush;
        os << " minutes " << (time_taken) % 60 << " seconds" << std::endl;
        writeInfo(os.str());
    }

    void Tree::openSQLDatabase()
    {
        if(!database->isOpen())
        {
#ifdef sql_ram
            database->open(":memory:");
#endif
#ifndef sql_ram
            database->open(sql_output_database);
#endif
        }
    }

    void Tree::sqlCreate()
    {
        time(&out_finish);
        std::stringstream os;
        os << "Creating SQL database file..." << std::endl;
        os << "\tChecking for existing folders...." << std::endl;
        writeInfo(os.str());
        os.str("");
        // Create the folder if it doesn't exist
        setupOutputDirectory();
        os.str("");
        os << "\tGenerating species list...." << std::endl;
        writeInfo(os.str());
        // Open a SQL database in memory. This will be written to disk later.
        // A check here can be done to write to disc directly instead to massively reduce RAM consumption
        openSQLDatabase();
        setupCommunity();
        // Create the command to be executed by adding to the string.
        community.createSpeciesList();
        community.writeSpeciesList(enddata);
        // Vacuum the file so that the file size is reduced (reduces by around 3%)
        try
        {
            database->execute("VACUUM;");
        }
        catch(FatalException &fe)
        {
            std::stringstream ss;
            ss << "Error thrown whilst vacuuming the database: " << fe.what() << std::endl;
            ss << "Continuing..." << std::endl;
            writeCritical(ss.str());
        }
        sqlCreateSimulationParameters();
    }

    void Tree::setupOutputDirectory()
    {
        if(sql_output_database == "null" || !fs::exists(fs::path(sql_output_database)))
        {
            sql_output_database = out_directory;
            string sqlfolder = out_directory;
            try
            {
                if(!fs::exists(fs::path(sqlfolder)))
                {
                    fs::create_directory(fs::path(sqlfolder));
                }
                sql_output_database += string("/data_") + std::to_string(task) + "_" + std::to_string(seed) + ".db";
            }
            catch(FatalException &fe)
            {
                writeWarning(fe.what());
                sql_output_database = string("data_") + std::to_string(task) + "_" + std::to_string(seed) + ".db";
            }
            fs::remove(fs::path(sql_output_database));
        }
    }

    void Tree::sqlCreateSimulationParameters()
    {
        // Now additionally store the simulation current_metacommunity_parameters (extremely useful data)
        string to_execute = "CREATE TABLE SIMULATION_PARAMETERS (seed INT PRIMARY KEY not null, task INT NOT NULL,";
        to_execute += "output_dir TEXT NOT NULL, speciation_rate DOUBLE NOT NULL, sigma DOUBLE NOT NULL,tau DOUBLE NOT "
                      "NULL, deme DOUBLE NOT NULL, ";
        to_execute += "sample_size DOUBLE NOT NULL, max_time INT NOT NULL, dispersal_relative_cost DOUBLE NOT NULL, "
                      "min_num_species ";
        to_execute += "INT NOT NULL, habitat_change_rate DOUBLE NOT NULL, gen_since_historical DOUBLE NOT NULL, ";
        to_execute += "time_config_file TEXT NOT NULL, coarse_map_file TEXT NOT NULL, coarse_map_x INT NOT NULL, "
                      "coarse_map_y INT NOT NULL,";
        to_execute += "coarse_map_x_offset INT NOT NULL, coarse_map_y_offset INT NOT NULL, coarse_map_scale DOUBLE NOT "
                      "NULL, fine_map_file TEXT NOT NULL, fine_map_x INT NOT NULL,";
        to_execute += "fine_map_y INT NOT NULL, fine_map_x_offset INT NOT NULL, fine_map_y_offset INT NOT NULL, ";
        to_execute += "sample_file TEXT NOT NULL, grid_x INT NOT NULL, grid_y INT NOT NULL, sample_x INT NOT NULL, ";
        to_execute += "sample_y INT NOT NULL, sample_x_offset INT NOT NULL, sample_y_offset INT NOT NULL, ";
        to_execute += "historical_coarse_map TEXT NOT NULL, historical_fine_map TEXT NOT NULL, sim_complete INT NOT NULL, ";
        to_execute += "dispersal_method TEXT NOT NULL, m_probability DOUBLE NOT NULL, cutoff DOUBLE NOT NULL, ";
        to_execute += "restrict_self INT NOT NULL, landscape_type TEXT NOT NULL, protracted INT NOT NULL, ";
        to_execute += "min_speciation_gen DOUBLE NOT NULL, max_speciation_gen DOUBLE NOT NULL, dispersal_map TEXT NOT NULL);";
        database->execute(to_execute);
        to_execute = simulationParametersSqlInsertion();
        database->execute(to_execute);
    }

    string Tree::simulationParametersSqlInsertion()
    {
        string to_execute;
        std::stringstream ss1, ss2;
        ss1 << std::setprecision(64);
        ss1 << spec;
        to_execute = "INSERT INTO SIMULATION_PARAMETERS VALUES(" + std::to_string((long long) seed) + ","
                     + std::to_string((long long) task);
        to_execute += ",'" + out_directory + "'," + ss1.str() + "," + std::to_string(0.0) + ",";
        to_execute += std::to_string(0.0) + "," + std::to_string((long double) deme) + ",";
        to_execute += std::to_string((long double) deme_sample) + "," + std::to_string((long long) maxtime) + ",";
        to_execute += std::to_string(0.0) + "," + std::to_string(0.0) + ",";
        to_execute += std::to_string((long double) sim_parameters->habitat_change_rate) + ",";
        to_execute += std::to_string((long double) sim_parameters->gen_since_historical) + ",'" + sim_parameters->times_file
                      + "','";
        to_execute += "none', 0, 0, 0, 0, 0, 'null', 0, 0, 0, 0, 'none', 1, 1, 1, 1, 0, 0, 'none', 'none',";
        to_execute += std::to_string(sim_complete);
        to_execute += ", 'none', 0.0, 0, 0, 'none', ";
        // Now save the protracted speciation variables (not relevant in this simulation scenario)
        to_execute += protractedVarsToString();
        to_execute += ", 'none');";
        return to_execute;
    }

    string Tree::protractedVarsToString()
    {
        string tmp = std::to_string(false) + ", " + std::to_string(0.0) + ", " + std::to_string(0.0);
        return tmp;
    }

    void Tree::simPause()
    {
        auto out1 = initiatePause();
        dumpMain(out1);
        dumpActive(out1);
        dumpData(out1);
        completePause(out1);
    }

    shared_ptr<std::ofstream> Tree::initiatePause()
    {
        std::stringstream os;
        os << "Pausing simulation..." << std::endl << "Saving data to temp file in " << out_directory << "/Pause/ ..."
           << std::flush;
        writeInfo(os.str());
        os.str("");
        // Create the pause directory
        string pause_folder = out_directory + "/Pause/";
        fs::path pause_dir(pause_folder);
        if(!fs::exists(pause_dir))
        {
            try
            {
                fs::create_directory(pause_dir);
            }
            catch(std::exception &e)
            {
                std::stringstream ss;
                ss << "Failure to create " << out_directory << "/Pause/" << "." << std::endl;
                ss << e.what() << std::endl;
                ss << "Writing directly to output directory." << std::endl;
                writeError(ss.str());
                pause_folder = out_directory;
            }
        }
        string file_to_open = pause_folder + "Dump_main_" + std::to_string(task) + "_" + std::to_string(seed) + ".csv";
        shared_ptr<std::ofstream> out = make_shared<std::ofstream>();
        out->open(file_to_open.c_str());
        *out << std::setprecision(64);
        return out;
    }

    void Tree::completePause(shared_ptr<std::ofstream> out)
    {
        out->close();
        std::stringstream os;
        os << "done." << std::endl;
        os << "SQL dump started" << std::endl;
        writeInfo(os.str());
        os.str("");
        time(&out_finish);
        sqlCreate();
        sqlOutput();
        os << "Data dump complete" << std::endl;
        writeInfo(os.str());
        time(&sim_end);
        writeTimes();
    }

    void Tree::dumpMain(shared_ptr<std::ofstream> out)
    {
        try
        {
            // Save that this simulation was not a protracted speciation sim
            *out << bIsProtracted << "\n";
            // Saving the initial data to one file.
            *out << enddata << "\n" << seeded << "\n" << seed << "\n" << task << "\n" << times_file << "\n"
                 << uses_temporal_sampling << "\n";
            *out << out_directory << "\n";
            *out << has_imported_vars << "\n" << start << "\n" << sim_start << "\n";
            *out << sim_end << "\n" << now << "\n" << time_taken << "\n" << sim_finish << "\n" << out_finish << "\n";
            *out << endactive << "\n" << startendactive << "\n" << maxsimsize << "\n" << steps << "\n";
            *out << generation << "\n" << "\n" << maxtime << "\n";
            *out << deme_sample << "\n" << spec << "\n" << deme << "\n";
            *out << sql_output_database << "\n" << *NR << "\n" << *sim_parameters << "\n";
            // now output the protracted speciation variables (there should be two of these).
            *out << getProtractedVariables();
        }
        catch(std::exception &e)
        {
            std::stringstream ss;
            ss << "Failed to perform dump of main: " << e.what() << std::endl;
            writeError(ss.str());
        }
    }

    void Tree::dumpActive(shared_ptr<std::ofstream> out)
    {
        try
        {
            // Output the active object
            *out << active;
        }
        catch(std::exception &e)
        {
            std::stringstream ss;
            ss << "Failed to perform dump of active: " << e.what() << std::endl;
            writeError(ss.str());
        }
    }

    void Tree::dumpData(shared_ptr<std::ofstream> out)
    {
        try
        {
            // Output the data object
            *out << *data;
        }
        catch(std::exception &e)
        {
            std::stringstream ss;
            ss << "Failed to perform dump of data: " << e.what() << std::endl;
            writeError(ss.str());
        }
    }

    void Tree::setResumeParameters()
    {
        if(!has_imported_pause)
        {
            pause_sim_directory = out_directory;
            has_imported_pause = true;
        }
    }

    shared_ptr<std::ifstream> Tree::openSaveFile()
    {
        shared_ptr<std::ifstream> in1 = make_shared<std::ifstream>();
        string file_to_open =
                pause_sim_directory + string("/Pause/Dump_main_") + std::to_string(task) + "_" + std::to_string(seed)
                + string(".csv");
        in1->open(file_to_open);
        if(!*in1)
        {
            std::stringstream es;
            es << "Cannot open file at " << file_to_open << std::endl;
            throw FatalException(es.str());
        }
        return in1;
    }

    void Tree::setResumeParameters(string pausedir,
                                   string outdir,
                                   unsigned long seed,
                                   unsigned long task,
                                   unsigned long new_max_time)
    {
        if(!has_imported_pause)
        {
            pause_sim_directory = move(pausedir);
            out_directory = move(outdir);
            this->seed = static_cast<long long int>(seed);
            this->task = static_cast<long long int>(task);
            maxtime = new_max_time;
            has_imported_pause = true;
        }
    }

    void Tree::loadMainSave(shared_ptr<std::ifstream> in1)
    {
        try
        {
            std::stringstream os;
            os << "\rLoading data from temp file...main..." << std::flush;
            writeInfo(os.str());
            os.str("");
            // Reading the initial data
            string string1;
            // First read our boolean which just determines whether the simulation is a protracted simulation or not.
            // For these simulations, it should not be.
            bool tmp;
            *in1 >> tmp;
            if(tmp != getProtracted())
            {
                if(getProtracted())
                {
                    throw FatalException("Paused simulation is not a protracted speciation simulation. "
                                         "Cannot be resumed by this program. Please report this bug");
                }
                else
                {
                    throw FatalException("Paused simulation is a protracted speciation simulation. "
                                         "Cannot be resumed by this program. Please report this bug");
                }
            }
            *in1 >> enddata >> seeded >> seed >> task;
            in1->ignore(); // Ignore the endline character
            getline(*in1, times_file);
            *in1 >> uses_temporal_sampling;
            in1->ignore();
            getline(*in1, string1);
            time_t tmp_time;
            *in1 >> has_imported_vars >> tmp_time;
            *in1 >> sim_start >> sim_end >> now;
            *in1 >> time_taken >> sim_finish >> out_finish >> endactive >> startendactive >> maxsimsize >> steps;
            unsigned long tempmaxtime = maxtime;
            *in1 >> generation >> maxtime;
            has_imported_vars = false;
            *in1 >> deme_sample >> spec >> deme;
            in1->ignore();
            getline(*in1, sql_output_database);
            *in1 >> *NR;
            in1->ignore();
            *in1 >> *sim_parameters;
            if(maxtime == 0)
            {
                sim_parameters->max_time = tempmaxtime;
            }
#ifdef DEBUG
            if(maxtime == 0 && tempmaxtime == 0)
            {
                throw FatalException("Time set to 0 on resume!");
            }
#endif
            NR->setDispersalMethod(sim_parameters->dispersal_method, sim_parameters->m_prob, sim_parameters->cutoff);
            if(has_imported_pause)
            {
                sim_parameters->output_directory = out_directory;
            }
            setParameters();
            double tmp1, tmp2;
            *in1 >> tmp1 >> tmp2;
            setProtractedVariables(tmp1, tmp2);
            if(times_file == "null")
            {
                if(uses_temporal_sampling)
                {
                    throw std::runtime_error("uses_temporal_sampling should not be true");
                }
            }
            else
            {
                if(!uses_temporal_sampling)
                {
                    throw std::runtime_error("uses_temporal_sampling should not be false");
                }
                vector<string> tmpimport;
                ConfigParser tmpconfig;
                tmpconfig.setConfig(times_file, false);
                tmpconfig.importConfig(tmpimport);
                for(const auto &i : tmpimport)
                {
                    reference_times.push_back(stod(i));
                    //					os << "t_i: " << reference_times[i] << std::endl;
                }
            }
        }
        catch(std::exception &e)
        {
            string msg;
            msg = "Failure to import current_metacommunity_parameters from temp main: " + string(e.what());
            throw FatalException(msg);
        }
    }

    void Tree::loadDataSave(shared_ptr<std::ifstream> in1)
    {
        try
        {
            std::stringstream os;
            os << "\rLoading data from temp file...data..." << std::flush;
            writeInfo(os.str());
            *in1 >> *data;
        }
        catch(std::exception &e)
        {
            string msg;
            msg = "Failure to import data from temp data: " + string(e.what());
            throw FatalException(msg);
        }
    }

    void Tree::loadActiveSave(shared_ptr<std::ifstream> in1)
    {
        string file_to_open;
        try
        {
            std::stringstream os;
            os << "\rLoading data from temp file...active..." << std::flush;
            writeInfo(os.str());
            // Input the active object
            *in1 >> active;
        }
        catch(std::exception &e)
        {
            string msg;
            msg = "Failure to import data from temp active: " + string(e.what());
            throw FatalException(msg);
        }
    }

    void Tree::initiateResume()
    {
        // Start the timer
        // Only resume the simulation if there is a simulation to resume from.
        if(!has_paused)
        {
            return;
        }
        time(&start);
        // Loads the data from the files into the relevant objects.
        std::stringstream os;
#ifdef DEBUG
        writeLog(10, "Paused directory: " + pause_sim_directory);
        writeLog(10, "Output directory: " + out_directory);
        writeLog(10, "Seed: " + std::to_string(seed));
        writeLog(10, "Task: " + std::to_string(task));
        writeLog(10, "Max time: " + std::to_string(maxtime));
#endif // DEBUG
        os << "Resuming simulation..." << std::endl << "Loading data from temp file..." << std::flush;
        writeInfo(os.str());
        os.str("");

    }

    void Tree::simResume()
    {
        initiateResume();
        // open the save file
        auto is = openSaveFile();
        // now load the objects
        loadMainSave(is);
        setObjectSizes();
        loadActiveSave(is);
        loadDataSave(is);
        time(&sim_start);
        writeInfo("\rLoading data from temp file...done.\n");
    }

    void Tree::addGillespie(const double &g_threshold)
    {
        std::stringstream ss;
        ss << "The gillespie algorithm is not supported for non-spatial coalescence trees yet. "
           << " Cannot run with Gillespie threshold of " << g_threshold << ". "
           << "Please contact the project maintainer if this is a feature you would like to see.";
        throw FatalException(ss.str());
    }

    bool Tree::runSimulationGillespie()
    {
        throw FatalException("The gillespie algorithm is not supported for non-spatial coalescence trees yet. "
                             "Please contact the project maintainer if this is a feature you would like to see.");
    }

#ifdef DEBUG

    void Tree::validateLineages()
    {
        bool fail = false;
        writeInfo("\nStarting lineage validation...");
        for(unsigned long i = 1; i < endactive; i++)
        {
            std::stringstream ss;
            DataPoint tmp_datapoint = active[i];
            if(tmp_datapoint.getXwrap() == 0 && tmp_datapoint.getYwrap() == 0)
            {
                if(tmp_datapoint.getNwrap() != 0)
                {
                    fail = true;
                }
            }
            else
            {
                fail = true;
            }
            if(fail)
            {
                ss << "\nFailure in map expansion. Please report this bug." << std::endl;
                ss << "active reference: " << i << std::endl;
                (*data)[active[i].getReference()].logLineageInformation(50);
                throw FatalException(ss.str());
            }
        }
        writeInfo("done.\n");
        validateCoalescenceTree();
    }

    void Tree::validateCoalescenceTree()
    {
        writeInfo("Validating coalescence tree...");
        // Get the active lineages
        std::set<unsigned long> active_lineage_refs;
        for(unsigned long i = 1; i < endactive + 1; i++)
        {
            active_lineage_refs.insert(active[i].getReference());
        }
        try
        {
            for(unsigned long i = 1; i < enddata - 1; i++)
            {
                if(active_lineage_refs.count(i) == 0)
                {
                    const auto &tree_node = (*data)[i];
                    if(checkSpeciation(tree_node.getSpecRate(), spec, tree_node.getGenerationRate()))
                    {
                        if(tree_node.getParent() != 0)
                        {
                            std::stringstream ss;
                            ss << "Tree node at " << i << " can speciate, but parent is not 0. Please report this bug."
                               << std::endl;
                            throw FatalException(ss.str());
                        }
                    }
                    else
                    {
                        if(tree_node.getParent() == 0)
                        {

                            std::stringstream ss;
                            ss << "Tree node at " << i << " has not speciated, but parent is 0. Please report this bug."
                               << std::endl;
                            unsigned long j = 0;
                            for(unsigned long z = 0; z < endactive; z++)
                            {
                                if(active[z].getReference() == i)
                                {
                                    j = z;
                                    break;
                                }
                            }
                            ss << "Location in active is: " << j << std::endl;
                            throw FatalException(ss.str());
                        }
                    }
                }
            }
            writeInfo("done.\n");
        }
        catch(FatalException &fe)
        {
            std::stringstream ss;
            ss << "Error validating coalescence tree: " << fe.what();
            throw FatalException(ss.str());
        }
    }

    void Tree::debugEndStep()
    {
        try
        {
            runChecks(this_step.chosen, this_step.coalchosen);
            // runs the debug every 10,000 time steps
            if(steps % 10000 == 0)
            {
                for(unsigned long i = 0; i <= endactive; i++)
                {
                    runChecks(i, i);
                }
            }
        }
        catch(FatalException &fe)
        {
            writeLog(50, "Logging chosen:");
            active[this_step.chosen].logActive(50);
            writeLog(50, "Logging coalchosen");
            active[this_step.coalchosen].logActive(50);
            std::stringstream ss;
            ss << "dumping data file..." << std::endl;
            sqlCreate();
#ifdef sql_ram
            sqlOutput();
#endif
            ss << "done." << std::endl;
            writeWarning(ss.str());
            throw fe;
        }

    }

    void Tree::debugCoalescence()
    {
        if(this_step.coalchosen == 0)
        {
            return;
        }
        std::stringstream ss;
        if(active[this_step.coalchosen].getXpos() != active[this_step.chosen].getXpos()
           || active[this_step.coalchosen].getYpos() != active[this_step.chosen].getYpos()
           || active[this_step.coalchosen].getXwrap() != active[this_step.chosen].getXwrap()
           || active[this_step.coalchosen].getYwrap() != active[this_step.chosen].getYwrap())
        {
            writeLog(50, "Logging chosen: " + std::to_string(this_step.chosen));
            (*data)[active[this_step.chosen].getReference()].logLineageInformation(50);
            writeLog(50, "Logging coalchosen: " + std::to_string(this_step.coalchosen));
            (*data)[active[this_step.coalchosen].getReference()].logLineageInformation(50);
            ss << "Nwrap not set correctly. Check move programming function." << std::endl;
            throw FatalException(ss.str());
        }
        if(active[this_step.coalchosen].getXpos() != (unsigned long) this_step.x
           || active[this_step.coalchosen].getYpos() != (unsigned long) this_step.y
           || active[this_step.coalchosen].getXwrap() != this_step.xwrap
           || active[this_step.coalchosen].getYwrap() != this_step.ywrap)
        {
            writeLog(50, "Logging chosen: " + std::to_string(this_step.chosen));
            (*data)[active[this_step.chosen].getReference()].logLineageInformation(50);
            writeLog(50, "Logging coalchosen: " + std::to_string(this_step.coalchosen));
            (*data)[active[this_step.coalchosen].getReference()].logLineageInformation(50);
            ss << "Nwrap not set correctly. Check move programming function." << std::endl;
            throw FatalException(ss.str());
        }
    }

    void Tree::runChecks(const unsigned long &chosen, const unsigned long &coalchosen)
    {
        miniCheck(chosen);
        miniCheck(coalchosen);
    }

    void Tree::miniCheck(const unsigned long &chosen)
    {
        if(chosen == 0)
        {
            return;
        }
        if(active[chosen].getReference() == 0)
        {
            throw FatalException("Active reference should not be 0.");
        }
        if((*data)[active[chosen].getReference()].getParent() != 0)
        {
            writeLog(50, "Active: " + std::to_string(chosen));
            (*data)[active[chosen].getReference()].logLineageInformation(50);
            throw FatalException("Parent not set to 0 for active lineage.");
        }
    }

#endif // DEBUG
}