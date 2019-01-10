// This file is part of necsim project which is released under MIT license.
// See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.

/**
 * @file Community.cpp
 * @brief Contains various objects used for reconstructing the coalescence tree after simulations are complete.
 *
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 */
//#define use_csv
#include <algorithm>
#include <set>
#include <unordered_map>
#include <utility>
#include <numeric>
#include "Community.h"
#include "double_comparison.h"
#include "file_system.h"

bool checkSpeciation(const long double &random_number, const long double &speciation_rate,
					 const unsigned long &no_generations)
{
	// bool result1, result2, result3, result4;
	long double res = 1.0 - pow(double(1.0 - speciation_rate), double(no_generations));
	return random_number <= res;
}

Samplematrix::Samplematrix()
{
	bIsFragment = false;
	bIsNull = false;
}

bool Samplematrix::getTestVal(unsigned long xval, unsigned long yval, long xwrap, long ywrap)
{
	return getVal(xval, yval, xwrap, ywrap);
}

bool Samplematrix::getMaskVal(unsigned long x1, unsigned long y1, long x_wrap, long y_wrap)
{
	if(bIsFragment)
	{
		long x, y;
		x = x1 + (x_wrap * x_dim) + x_offset;
		y = y1 + (y_wrap * y_dim) + y_offset;
		return fragment.x_west <= x && x <= fragment.x_east && fragment.y_north <= y &&
			   y <= fragment.y_south;
	}
	if(bIsNull)
	{
		return true;
	}
	return getVal(x1, y1, x_wrap, y_wrap);
}

void Samplematrix::setFragment(Fragment &fragment_in)
{
	fragment = fragment_in;
	bIsFragment = true;
}

void Samplematrix::removeFragment()
{
	bIsFragment = false;
}

void Community::setList(shared_ptr<vector<TreeNode>> l)
{
	nodes = std::move(l);
}

void Community::setDatabase(sqlite3 *dbin)
{
	if(!database_set)
	{
		database = dbin;
	}
	else
	{
		throw SpeciesException("ERROR_SPEC_002: Attempt to set database - database link has already been set");
	}
	database_set = true;  // this just specifies that the database has been created in memory.
}

bool Community::hasImportedData()
{
	return has_imported_data;
}

long double Community::getMinimumSpeciation()
{
	return min_spec_rate;
}

void Community::importSamplemask(string sSamplemask)
{
	// Check that the sim data has been imported.
	if(!has_imported_data)
	{
		throw SpeciesException(
				"ERROR_SPEC_003: Attempt to import samplemask object before simulation current_metacommunity_parameters: dimensions not known");
	}
	// Check that the main data has been imported already, otherwise the dimensions of the samplemask will not be correct
	if(!has_imported_samplemask)
	{
		stringstream os;
		samplemask.importBooleanMask(grid_x_size, grid_y_size, samplemask_x_size, samplemask_y_size,
									 samplemask_x_offset, samplemask_y_offset, sSamplemask);
		if(sSamplemask != "null")
		{
			unsigned long total = 0;
			for(unsigned long i = 0; i < samplemask.sample_mask.getCols(); i++)
			{
				for(unsigned long j = 0; j < samplemask.sample_mask.getRows(); j++)
				{
					if(samplemask.sample_mask[j][i])
					{
						total++;
					}
				}
			}
			os << "Sampling " << total << " cells." << endl;
		}
		else
		{

#ifdef DEBUG
			os << "Sampling all areas." << endl;
#endif
		}
		writeInfo(os.str());
		has_imported_samplemask = true;
	}
}

unsigned long Community::countSpecies()
{
	unsigned int precount = 0;
	for(unsigned long i = 1; i < nodes->size(); i++)
	{
		if((*nodes)[i].hasSpeciated())
		{
			precount++;
		}
	}
	return precount;
}

unsigned long Community::calcSpecies()
{
	resetTree();
	if(!has_imported_samplemask)
	{
#ifdef DEBUG
		writeInfo("No samplemask imported. Defaulting to null.\n");
#endif
		importSamplemask("null");
	}
	unsigned long species_count = 0;  // start at 2 because the last species has been burnt already.
	// check that tips exist within the spatial and temporal frame of interest.
#ifdef DEBUG
	writeLog(10, "Assigning tips.");
#endif // DEBUG
	for(unsigned long i = 1; i < nodes->size(); i++)
	{
		TreeNode *this_node = &(*nodes)[i];
#ifdef DEBUG
		if((*nodes)[i].getParent() >= nodes->size())
		{
			writeLog(50, "i: " + to_string(i));
			this_node->logLineageInformation(50);
			writeLog(50, "size: " + to_string(nodes->size()));
			throw SpeciesException("ERROR_SQL_017: The parent is outside the size of the the data object. Bug "
										   "in expansion of data structures or object set up likely.");
		}
#endif //DEBUG
		this_node->setExistence(this_node->isTip() && samplemask.getMaskVal(this_node->getXpos(), this_node->getYpos(),
																			this_node->getXwrap(),
																			this_node->getYwrap()) &&
								doubleCompare(this_node->getGeneration(), current_community_parameters->time, 0.0001));
		// Calculate if speciation occured at any point in the lineage's branch
		if(protracted)
		{
			long double lineage_age = this_node->getGeneration() + this_node->getGenRate();
			if(lineage_age >= applied_protracted_parameters.min_speciation_gen)
			{
				if(checkSpeciation(this_node->getSpecRate(), current_community_parameters->speciation_rate,
								   this_node->getGenRate()))
				{
					this_node->speciate();
				}
				if(lineage_age >= applied_protracted_parameters.max_speciation_gen)
				{
					this_node->speciate();
				}
			}
		}
		else
		{
			if(checkSpeciation(this_node->getSpecRate(), current_community_parameters->speciation_rate,
							   this_node->getGenRate()))
			{
				this_node->speciate();
			}
		}
	}

#ifdef DEBUG
	writeLog(10, "Calculating lineage existence.");
#endif // DEBUG
	// now continue looping to calculate species identities for lineages given the new speciation probabilities.
	bool bSorter = true;
	while(bSorter)
	{
		bSorter = false;
		for(unsigned long i = 1; i < nodes->size(); i++)
		{
			TreeNode *this_node = &(*nodes)[i];
			// check if any parents exist
			if(!(*nodes)[this_node->getParent()].exists() && this_node->exists() &&
			   !this_node->hasSpeciated())
			{
				bSorter = true;
				(*nodes)[this_node->getParent()].setExistence(true);
			}
		}
	}
#ifdef DEBUG
	writeLog(10, "Speciating lineages.");
#endif // DEBUG
	species_count = 0;
	set<unsigned long> species_list;
	// Now loop again, creating a new species for each species that actually exists.
	auto start = nodes->begin();
	start++;
	for(auto item = start; item != nodes->end(); item++)
	{
		if(item->exists() && item->hasSpeciated())
		{
			addSpecies(species_count, &(*item), species_list);
		}
	}

	// Compress the species IDs so that the we have full mapping of species_ids to integers in range 0:n
	// Only do so if the numbers do not match initially

#ifdef DEBUG
	writeLog(10, "Counting species.");
#endif // DEBUG
	if(!species_list.empty() && species_count != *species_list.rbegin())
	{
		stringstream ss;
		ss << "initial count is " << species_count << " species" << endl;
		writeInfo(ss.str());
		ss.str("");
		ss << "Rescaling species ids for " << species_list.size() << " species...";
		writeInfo(ss.str());
		unsigned long initial_species_count = species_count;
		unsigned long tmp_species_count = 0;
		// Maps old to new species ids
		map<unsigned long, unsigned long> old_ids_to_new_ids;
//		old_ids_to_new_ids.reserve(species_count);
		for(unsigned long i = 1; i < nodes->size(); i++)
		{
			TreeNode *this_node = &(*nodes)[i];
			if(this_node->hasSpeciated() && this_node->exists())
			{
				auto map_id = old_ids_to_new_ids.find(this_node->getSpeciesID());
				if(map_id == old_ids_to_new_ids.end())
				{
					tmp_species_count++;
					old_ids_to_new_ids[this_node->getSpeciesID()] = tmp_species_count;
					this_node->resetSpecies();
					this_node->burnSpecies(tmp_species_count);
				}
				else
				{
					this_node->resetSpecies();
					this_node->burnSpecies(map_id->second);
				}
			}
		}
		if(tmp_species_count > initial_species_count)
		{
			writeInfo("\n");
			stringstream ss;
			ss << "Number of generated species (" << tmp_species_count << ") is more than the initial species count ";
			ss << "(" << initial_species_count << ") - please report this bug." << endl;
			throw FatalException(ss.str());
		}
		species_count = tmp_species_count;
		writeInfo("done!\n");
	}
	else
	{
		species_count = 0;
		for(unsigned long i = 0; i < nodes->size(); i++)
		{
			TreeNode *this_node = &(*nodes)[i];
			// count all speciation events, not just the ones that exist!
			if(this_node->hasSpeciated() && this_node->exists() && this_node->getSpeciesID() != 0)
			{
				species_count++;
			}
		}
	}
	// now loop to correctly assign each species id
	bool loopon = true;
	bool error_printed = false;
#ifdef DEBUG
	writeLog(10, "Generating species IDs.");
#endif // DEBUG
	while(loopon)
	{
		loopon = false;
		// if we start at the end of the loop and work backwards, we should remove some of the repeat
		// speciation events.
		for(unsigned long i = (nodes->size()) - 1; i > 0; i--)
		{
			TreeNode *this_node = &(*nodes)[i];
			//				os << i << endl;
			if(this_node->getSpeciesID() == 0 && this_node->exists())
			{
				loopon = true;
				this_node->burnSpecies((*nodes)[this_node->getParent()].getSpeciesID());
#ifdef DEBUG
				if((*nodes)[this_node->getParent()].getSpeciesID() == 0 &&
				   doubleCompare(this_node->getGeneration(), current_community_parameters->time, 0.001))
				{

					if(!error_printed)
					{
						stringstream ss;
						ss << "Potential parent ID error in " << i << " - incomplete simulation likely." << endl;
						writeCritical(ss.str());
						writeLog(50, "Lineage information:");
						this_node->logLineageInformation(50);
						writeLog(50, "Parent information:");
						(*nodes)[this_node->getParent()].logLineageInformation(50);
						error_printed = true;
						break;
					}
				}
#endif

			}
		}
		if(error_printed)
		{
			throw FatalException("Parent ID error when calculating coalescence tree.");
		}
	}
	// count the number of species that have been created
#ifdef DEBUG
	writeLog(10, "Completed tree creation.");
#endif // DEBUG
	iSpecies = species_count;
	//		os << "iSpecies: " << iSpecies << endl;
	return species_count;
}

void Community::addSpecies(unsigned long &species_count, TreeNode *tree_node, set<unsigned long> &species_list)
{
	species_count++;
	tree_node->burnSpecies(species_count);
}

void Community::calcSpeciesAbundance()
{
	species_abundances = make_shared<vector<unsigned long>>();
	species_abundances->resize(iSpecies + 1, 0);
	for(unsigned long i = 1; i < nodes->size(); i++)
	{
		TreeNode *this_node = &(*nodes)[i];
		if(this_node->isTip() &&
		   doubleCompare(this_node->getGeneration(), current_community_parameters->time, 0.0001) &&
		   this_node->exists())
		{
#ifdef DEBUG
			if(this_node->getSpeciesID() >= species_abundances->size())
			{
				throw out_of_range("Node index out of range of abundances size. Please report this bug.");
			}
#endif // DEBUG
			// The line that counts the number of individuals
			species_abundances->operator[](this_node->getSpeciesID())++;
#ifdef DEBUG
			if(!samplemask.getMaskVal(this_node->getXpos(), this_node->getYpos(),
									  this_node->getXwrap(), this_node->getYwrap()) &&
			   doubleCompare(this_node->getGeneration(), current_community_parameters->time, 0.0001))
			{
				stringstream ss;
				ss << "x,y " << (*nodes)[i].getXpos() << ", " << (*nodes)[i].getYpos() << endl;
				ss << "tip: " << (*nodes)[i].isTip() << " Existance: " << (*nodes)[i].exists()
					 << " samplemask: " << samplemask.getMaskVal(this_node->getXpos(), this_node->getYpos(),
																 this_node->getXwrap(), this_node->getYwrap()) << endl;
				ss
						<< "ERROR_SQL_005: Tip doesn't exist. Something went wrong either in the import or "
								"main simulation running."
						<< endl;
				writeWarning(ss.str());

			}
			if(this_node->getSpeciesID() == 0 && samplemask.getMaskVal(this_node->getXpos(), this_node->getYpos(),
																	   this_node->getXwrap(), this_node->getYwrap()) &&
			   doubleCompare(this_node->getGeneration(), current_community_parameters->time, 0.0001))
			{
				stringstream ss;
				ss << "x,y " << this_node->getXpos() << ", " << this_node->getYpos() << endl;
				ss << "generation (point,required): " << this_node->getGeneration() << ", "
				   << current_community_parameters->time << endl;
				TreeNode *p_node = &(*nodes)[this_node->getParent()];
				ss << "samplemasktest: " << samplemask.getTestVal(this_node->getXpos(), this_node->getYpos(),
																  this_node->getXwrap(), this_node->getYwrap()) << endl;
				ss << "samplemask: " << samplemask.getVal(this_node->getXpos(), this_node->getYpos(),
														  this_node->getXwrap(), this_node->getYwrap()) << endl;
				ss << "parent (tip, exists, generations): " << p_node->isTip() << ", "
				   << p_node->exists() << ", " << p_node->getGeneration() << endl;
				ss << "species id zero - i: " << i << " parent: " << p_node->getParent()
				   << " speciation_probability: " << p_node->getSpecRate() << "has speciated: " << p_node->hasSpeciated()
				   << endl;
				writeCritical(ss.str());
				throw runtime_error("Fatal, exiting program.");
			}
#endif
		}
	}

}

void Community::resetTree()
{
	for(auto &item: *nodes)
	{
		item.qReset();
	}
}

void Community::detectDimensions(string db)
{
	sqlite3 *tmpdb;
	int rc = sqlite3_open_v2(db.c_str(), &tmpdb, SQLITE_OPEN_READWRITE, "unix-dotfile");
	string to_exec = "SELECT MAX(xval),MAX(yval) FROM SPECIES_LIST;";
	sqlite3_stmt *stmt = nullptr;
	rc = sqlite3_prepare_v2(tmpdb, to_exec.c_str(), static_cast<int>(strlen(to_exec.c_str())), &stmt, nullptr);
	unsigned long xvalmax, yvalmax;
	rc = sqlite3_step(stmt);
	xvalmax = static_cast<unsigned long>(sqlite3_column_int(stmt, 0) + 1);
	yvalmax = static_cast<unsigned long>(sqlite3_column_int(stmt, 1) + 1);
	samplemask.sample_mask.setSize(xvalmax, yvalmax);
	// close the old statement
	rc = sqlite3_finalize(stmt);
	if(rc != SQLITE_OK && rc != SQLITE_DONE)
	{
		stringstream ss;
		ss << "Could not detect dimensions: " << rc << " (" << sqlite3_errmsg(tmpdb) << ")" << endl;
		throw SpeciesException(ss.str());
	}
}

void Community::openSqlConnection(string inputfile)
{
	// open the database objects
	sqlite3_backup *backupdb = nullptr;
	sqlite3 *outdatabase = nullptr;
	// open one db in memory and one from the file.
	if(!boost::filesystem::exists(inputfile))
	{
		stringstream ss;
		ss << "Output database does not exist at " << inputfile << ": cannot open sql connection." << endl;
		throw FatalException(ss.str());
	}
	try
	{
		openSQLiteDatabase(":memory:", database);
		openSQLiteDatabase(inputfile, outdatabase);
		in_mem = true;
		// copy the db from file into memory.
		backupdb = sqlite3_backup_init(database, "main", outdatabase, "main");
		int rc = sqlite3_backup_step(backupdb, -1);

		if(rc != SQLITE_DONE && rc != SQLITE_OK)
		{
			sqlite3_close(outdatabase);
			sqlite3_open(inputfile.c_str(), &outdatabase);
			backupdb = sqlite3_backup_init(database, "main", outdatabase, "main");
		}
		rc = sqlite3_backup_finish(backupdb);
		//			os << "rc: " << rc << endl;
		if(rc != SQLITE_DONE && rc != SQLITE_OK)
		{
			sqlite3_close(database);
			sqlite3_close(outdatabase);
			throw SpeciesException("ERROR_SQL_002: FATAL. Source file cannot be opened.");
		}
		sqlite3_close(outdatabase);
	}
	catch(FatalException &fe)
	{
		writeWarning("Can't open in-memory database. Writing to file instead (this will be slower).\n");
		in_mem = false;
		sqlite3_close(database);
		int rc = sqlite3_open_v2(inputfile.c_str(), &database, SQLITE_OPEN_READWRITE, "unix-dotfile");
		// Revert to different VFS file opening method if the backup hasn't started properly.
		// Two different versions will be attempted before an error will be thrown.
		// A different way of assigning the VFS method and opening the file correctly could be implemented later.
		// Currently "unix-dotfile" works for HPC runs and "unix" works for PC runs.
		if(rc != SQLITE_OK)
		{
			throw SpeciesException("ERROR_SQL_002: FATAL. Source file cannot be opened. Error: " + string(fe.what()) +
								   " and " + to_string(rc));
		}
	}
	bSqlConnection = true;
}

void Community::closeSqlConnection()
{
	sqlite3_close(database);
	bSqlConnection = false;
}

void Community::setInternalDatabase()
{
	{
		if(!database_set)
		{
			openSQLiteDatabase(":memory:", database);
		}
		internalOption();
	}
}

void Community::internalOption()
{
	has_imported_data = true;
	bSqlConnection = true;
	database_set = true;
	in_mem = true;
}

void Community::importData(string inputfile)
{
	if(!bSqlConnection)
	{
		openSqlConnection(inputfile);
	}
	if(!has_imported_data)
	{
		importSimParameters(inputfile);
	}
	if(!nodes->empty())
	{
		return;
	}
	writeInfo("Beginning data import...");
	// The sql statement to store the sql statement message object
	sqlite3_stmt *stmt = nullptr;

	// Now find out the max size of the species_id_list, so we have a count to work from
	string count_command = "SELECT COUNT(*) FROM SPECIES_LIST;";
	sqlite3_prepare_v2(database, count_command.c_str(), static_cast<int>(strlen(count_command.c_str())), &stmt,
					   nullptr);
	unsigned long datasize;
	// skip first row (should be blank)
	sqlite3_step(stmt);
	datasize = static_cast<unsigned long>(sqlite3_column_int(stmt, 0));
	sqlite3_finalize(stmt);

	// Create db query
	string all_commands = "SELECT * FROM SPECIES_LIST;";
	sqlite3_prepare_v2(database, all_commands.c_str(), static_cast<int>(strlen(all_commands.c_str())), &stmt, nullptr);
	nodes->resize(datasize);
	// Check that the file opened correctly.
	sqlite3_step(stmt);
	// Copy the data across to the TreeNode data structure.
	// For storing the number of ignored lineages so this can be subtracted off the parent number.
	unsigned long ignored_lineages = 0;
#ifdef DEBUG
	bool has_printed_error = false;
#endif
	for(unsigned long i = 0; i < datasize; i++)
	{
		auto species_id = static_cast<unsigned long>(sqlite3_column_int(stmt, 1));
		//		os << species_id << endl;
		long xval = sqlite3_column_int(stmt, 2);
		long yval = sqlite3_column_int(stmt, 3);
		long xwrap = sqlite3_column_int(stmt, 4);
		long ywrap = sqlite3_column_int(stmt, 5);
		auto tip = bool(sqlite3_column_int(stmt, 6));
		auto speciation = bool(sqlite3_column_int(stmt, 7));
		auto parent = static_cast<unsigned long>(sqlite3_column_int(stmt, 8));
		auto iGen = static_cast<unsigned long>(sqlite3_column_int(stmt, 11));
		auto existence = bool(sqlite3_column_int(stmt, 9));
		double dSpec = sqlite3_column_double(stmt, 10);
		long double generationin = sqlite3_column_double(stmt, 12);
		// the -1 is to ensure that the species_id_list includes all lineages, but fills the output from the beginning
		unsigned long index = i - ignored_lineages;
		(*nodes)[index].setup(tip, xval, yval, xwrap, ywrap, generationin);
		(*nodes)[index].burnSpecies(species_id);
		(*nodes)[index].setSpec(dSpec);
		(*nodes)[index].setExistence(existence);
		(*nodes)[index].setGenerationRate(iGen);
		(*nodes)[index].setParent(parent - ignored_lineages);
		if(index == parent && parent != 0)
		{
			stringstream ss;
			ss << "ERROR_SQL_001: Import failed as parent is self. Please report this bug." << endl;
			ss << " i: " << index << " parent: " << parent << endl;
			throw SpeciesException(ss.str());
		}
		(*nodes)[index].setSpeciation(speciation);
		sqlite3_step(stmt);
#ifdef DEBUG
		if(parent < index && !speciation)
		{
			if(!has_printed_error)
			{
				stringstream ss;
				ss << "parent: " << parent << " index: " << index << endl;
				ss << "Parent before index error. Check program." << endl;
				has_printed_error = true;
				writeWarning(ss.str());
			}
		}
#endif
//		}
	}
	// Now we need to blank all objects
	sqlite3_finalize(stmt);
	// Now read the useful information from the SIMULATION_PARAMETERS table
	writeInfo("\rBeginning data import...done\n");
}

void Community::getMaxSpeciesAbundancesID()
{
	if(!bSqlConnection)
	{
		throw FatalException("Attempted to get from sql database without opening database connection.");
	}
	if(max_species_id == 0)
	{
		sqlite3_stmt *stmt = nullptr;
		// Now find out the max size of the species_id_list, so we have a count to work from
		string count_command = "SELECT MAX(ID) FROM SPECIES_ABUNDANCES;";
		sqlite3_prepare_v2(database, count_command.c_str(), static_cast<int>(strlen(count_command.c_str())), &stmt,
						   nullptr);
		sqlite3_step(stmt);
		max_species_id = static_cast<unsigned long>(sqlite3_column_int(stmt, 0)) + 1;
		// close the old statement
		sqlite3_finalize(stmt);
	}
}

shared_ptr<vector<unsigned long>> Community::getCumulativeAbundances()
{
	shared_ptr<vector<unsigned long>> out = make_shared<vector<unsigned long>>();
	out->reserve(species_abundances->size());
	partial_sum(species_abundances->begin(), species_abundances->end(), out->begin());
	return species_abundances;
}

shared_ptr<vector<unsigned long>> Community::getRowOut()
{
	return species_abundances;
}

unsigned long Community::getSpeciesNumber()
{
	return iSpecies;
}

void Community::getMaxSpeciesLocationsID()
{
	if(!bSqlConnection)
	{
		throw FatalException("Attempted to get from sql database without opening database connection.");
	}
	if(max_locations_id == 0)
	{
		sqlite3_stmt *stmt = nullptr;
		// Now find out the max size of the species_id_list, so we have a count to work from
		string count_command = "SELECT MAX(ID) FROM SPECIES_LOCATIONS;";
		sqlite3_prepare_v2(database, count_command.c_str(), static_cast<int>(strlen(count_command.c_str())), &stmt,
						   nullptr);
		sqlite3_step(stmt);
		max_locations_id = static_cast<unsigned long>(sqlite3_column_int(stmt, 0)) + 1;
		// close the old statement
		sqlite3_finalize(stmt);
	}
}

void Community::getMaxFragmentAbundancesID()
{
	if(!bSqlConnection)
	{
		throw FatalException("Attempted to get from sql database without opening database connection.");
	}
	if(max_fragment_id == 0)
	{
		sqlite3_stmt *stmt = nullptr;
		// Now find out the max size of the species_id_list, so we have a count to work from
		string count_command = "SELECT MAX(ID) FROM FRAGMENT_ABUNDANCES;";
		sqlite3_prepare_v2(database, count_command.c_str(), static_cast<int>(strlen(count_command.c_str())), &stmt,
						   nullptr);
		sqlite3_step(stmt);
		max_fragment_id = static_cast<unsigned long>(sqlite3_column_int(stmt, 0)) + 1;
		// close the old statement
		sqlite3_finalize(stmt);
	}
}

void Community::createDatabase()
{
	generateCoalescenceTree();
	stringstream os;
	os << "Generating new SQL table for speciation rate " << current_community_parameters->speciation_rate
	   << "..." << flush;
	writeInfo(os.str());
	string table_command = "CREATE TABLE IF NOT EXISTS SPECIES_ABUNDANCES (ID int PRIMARY KEY NOT NULL, "
						   "species_id INT NOT NULL, no_individuals INT NOT "
						   "NULL, community_reference INT NOT NULL);";
	int rc = sqlite3_exec(database, table_command.c_str(), nullptr, nullptr, nullptr);
	if(rc != SQLITE_OK)
	{
		throw SpeciesException("ERROR_SQL_002b: Could not create SPECIES_ABUNDANCES table.");
	}
	getMaxSpeciesAbundancesID();
	outputSpeciesAbundances();
}

void Community::generateCoalescenceTree()
{
	writeInfo("Calculating tree structure...");
	// Search through past speciation rates
	if(current_community_parameters->speciation_rate < min_spec_rate)
	{
		if(doubleCompare(current_community_parameters->speciation_rate, min_spec_rate, min_spec_rate * 0.000001))
		{
			current_community_parameters->speciation_rate = min_spec_rate;
		}
		else
		{
			stringstream ss;
			ss << "ERROR_SQL_018: Speciation rate of " << current_community_parameters->speciation_rate;
			ss << " is less than the minimum possible (" << min_spec_rate << ". Skipping." << endl;
			throw SpeciesException(ss.str());
		}
	}
	unsigned long nspec = calcSpecies();
	calcSpeciesAbundance();
	stringstream os;
	os << "done!" << endl;
	os << "Number of species: " << nspec << endl;
	writeInfo(os.str());
	os.str("");
}

void Community::outputSpeciesAbundances()
{
	// Only write to SPECIES_ABUNDANCES if the speciation rate has not already been applied
	if(!current_community_parameters->updated)
	{
//#ifdef DEBUG
		if(checkSpeciesAbundancesReference())
		{
			stringstream ss;
			ss << "Duplicate insertion of " << current_community_parameters->reference << " into SPECIES_ABUNDANCES.";
			ss << endl;
			writeWarning(ss.str());
			return;
		}
//#endif // DEBUG
		sqlite3_stmt *stmt = nullptr;
		string table_command = "INSERT INTO SPECIES_ABUNDANCES (ID, species_id, "
							   "no_individuals, community_reference) VALUES (?,?,?,?);";
		sqlite3_prepare_v2(database, table_command.c_str(), static_cast<int>(strlen(table_command.c_str())), &stmt,
						   nullptr);

		// Start the transaction
		sqlite3_exec(database, "BEGIN TRANSACTION;", nullptr, nullptr, nullptr);
		for(unsigned long i = 0; i < species_abundances->size(); i++)
		{
			// lexical cast fixes a precision problem that allows for printing of very small doubles.
			sqlite3_bind_int(stmt, 1, static_cast<int>(max_species_id++));
			sqlite3_bind_int(stmt, 2, static_cast<int>(i));
			sqlite3_bind_int(stmt, 3, static_cast<int>(species_abundances->operator[](i)));
			sqlite3_bind_int(stmt, 4, static_cast<int>(current_community_parameters->reference));
			int step = sqlite3_step(stmt);
			// makes sure the while loop doesn't go forever.
			time_t start_check, end_check;
			time(&start_check);
			time(&end_check);
			while(step != SQLITE_DONE && (end_check - start_check) < 1)
			{
				step = sqlite3_step(stmt);
				time(&end_check);
			}
			if(step != SQLITE_DONE)
			{
				stringstream os;
				os << "SQLITE error code: " << step << endl;
				os << "ERROR_SQL_004d: Could not insert into database. Check destination file has not "
					  "been moved or deleted and that an entry doesn't already exist with the same ID."
				   << endl;
				os << sqlite3_errmsg(database) << endl;
				sqlite3_clear_bindings(stmt);
				sqlite3_reset(stmt);
				throw FatalException(os.str());
			}
			sqlite3_clear_bindings(stmt);
			sqlite3_reset(stmt);

		}
		// execute the command and close the connection to the database
		int rc1 = sqlite3_exec(database, "END TRANSACTION;", nullptr, nullptr, nullptr);
		// Need to finalise the statement
		int rc2 = sqlite3_finalize(stmt);
		if(rc1 != SQLITE_OK || rc2 != SQLITE_OK)
		{
			writeError("ERROR_SQL_013: Could not complete SQL transaction. Check memory database assignment and "
					   "SQL commands. Ensure SQL statements are properly cleared and that you are not attempting "
					   "to insert repeat IDs into the database.");
		}
		else
		{
			stringstream ss;
			ss << "\rGenerating new SQL table for speciation rate " << current_community_parameters->speciation_rate
			   << "...done!" << endl;
			writeInfo(ss.str());
		}
	}
	else
	{
		stringstream ss;
		ss << "current_metacommunity_parameters already applied, not outputting SPECIES_ABUNDANCES table..." << endl;
		writeInfo(ss.str());
	}
}

bool
Community::checkCalculationsPerformed(const long double &speciation_rate, const double &time, const bool &fragments,
									  const MetacommunityParameters &metacomm_parameters,
									  const ProtractedSpeciationParameters &proc_parameters)
{
	auto metacommunity_reference = past_metacommunities.getReference(metacomm_parameters);
	if(metacommunity_reference == 0 && metacomm_parameters.isMetacommunityOption())
	{
		return false;
	}
	bool has_pair = past_communities.hasPair(speciation_rate, time, fragments, metacommunity_reference,
											 proc_parameters);
	if(fragments && past_communities.hasPair(speciation_rate, time, false, metacommunity_reference, proc_parameters))
	{
		return false;

	}
	if(!fragments && past_communities.hasPair(speciation_rate, time, true, metacommunity_reference, proc_parameters))
	{
		return true;
	}
//	if(past_communities.hasOption(speciation_rate, time, !fragments,
//								past_metacommunities.getReference(metacommunity_speciation_rate, metacommunity_size)))
//	{
//		return !fragments || has_pair;
//	}
	return has_pair;
}

void Community::createFragmentDatabase(const Fragment &f)
{
	//		os << "Generating new SQL table for speciation rate " << s << "..." << flush;
	string table_command = "CREATE TABLE IF NOT EXISTS FRAGMENT_ABUNDANCES (ID int PRIMARY KEY NOT NULL, fragment "
						   "TEXT NOT NULL, area DOUBLE NOT NULL, size INT NOT NULL,  species_id INT NOT NULL, "
						   "no_individuals INT NOT NULL, community_reference int NOT NULL);";
	sqlite3_exec(database, table_command.c_str(), nullptr, nullptr, nullptr);
	getMaxFragmentAbundancesID();
	sqlite3_stmt *stmt = nullptr;
	table_command = "INSERT INTO FRAGMENT_ABUNDANCES (ID, fragment, area, size, species_id, "
					"no_individuals, community_reference) VALUES (?,?,?,?,?,?,?);";
	sqlite3_prepare_v2(database, table_command.c_str(), static_cast<int>(strlen(table_command.c_str())), &stmt,
					   nullptr);

	// Start the transaction
	sqlite3_exec(database, "BEGIN TRANSACTION;", nullptr, nullptr, nullptr);
	for(unsigned long i = 0; i < species_abundances->size(); i++)
	{
		auto tmp_row = &species_abundances->operator[](i);
		if(*tmp_row != 0)
		{
			// fixed precision problem - lexical cast allows for printing of very small doubles.
			sqlite3_bind_int(stmt, 1, static_cast<int>(max_fragment_id++));
			sqlite3_bind_text(stmt, 2, f.name.c_str(), -1, SQLITE_STATIC);
			sqlite3_bind_double(stmt, 3, f.area);
			sqlite3_bind_int(stmt, 4, static_cast<int>(f.num));
			sqlite3_bind_int(stmt, 5, static_cast<int>(i));
			sqlite3_bind_int(stmt, 6, static_cast<int>(*tmp_row));
			sqlite3_bind_int(stmt, 7, static_cast<int>(current_community_parameters->reference));
			int step = sqlite3_step(stmt);
			// makes sure the while loop doesn't go forever.
			time_t start_check, end_check;
			time(&start_check);
			time(&end_check);
			while(step != SQLITE_DONE && (end_check - start_check) < 10)
			{
				step = sqlite3_step(stmt);
				time(&end_check);
			}
			if(step != SQLITE_DONE)
			{
				stringstream ss;
				ss << "ERROR_SQL_004e: Could not insert into database. Check destination file has not "
					  "been moved or deleted and that an entry doesn't already exist with the same ID."
				   << endl;
				ss << "SQLITE error code: " << step << endl;
				ss << sqlite3_errmsg(database) << endl;
				writeWarning(ss.str());
				sqlite3_clear_bindings(stmt);
				sqlite3_reset(stmt);
				break;
			}
			sqlite3_clear_bindings(stmt);
			sqlite3_reset(stmt);
		}
	}
	// execute the command and close the connection to the database
	int rc1 = sqlite3_exec(database, "END TRANSACTION;", nullptr, nullptr, nullptr);
	// Need to finalise the statement
	int rc2 = sqlite3_finalize(stmt);
	if(rc1 != SQLITE_OK || rc2 != SQLITE_OK)
	{
		writeError("ERROR_SQL_013: Could not complete SQL transaction. Check memory database assignment and SQL "
				   "commands. Ensure SQL statements are properly cleared and that you are not attempting to insert "
				   "repeat IDs into the database.");
	}
}

void Community::exportDatabase()
{
	if(in_mem)
	{
		stringstream ss;
		stringstream os;
		os << "Writing out to " << spec_sim_parameters->filename << "..." << flush;
		// Now write the database to the file object.
		sqlite3 *outdatabase2;
		writeInfo(os.str());
		int rc = sqlite3_open_v2(spec_sim_parameters->filename.c_str(), &outdatabase2, SQLITE_OPEN_READWRITE,
								 "unix-dotfile");
		// check that the connection to file has opened correctly
		if(rc != SQLITE_OK && rc != SQLITE_DONE)
		{
			// attempt other output method
			sqlite3_close(outdatabase2);
			rc = sqlite3_open(spec_sim_parameters->filename.c_str(), &outdatabase2);
			if(rc != SQLITE_OK && rc != SQLITE_DONE)
			{
				ss << "Connection to output database cannot be opened at " << spec_sim_parameters->filename;
				ss << ". Check write access on output folder. Error code " << rc << ": ";
				ss << sqlite3_errmsg(outdatabase2) << endl;
				throw FatalException(ss.str());
			}
		}

		// create the backup object to write data to the file from memory.

		sqlite3_backup *backupdb;
		backupdb = sqlite3_backup_init(outdatabase2, "main", database, "main");
		if(!backupdb)
		{
			ss << "ERROR_SQL_003: Could not backup to SQL database. Check destination file has not been "
				  "moved or deleted."
			   << endl;
			throw FatalException(ss.str());
		}
		// Perform the backup
		rc = sqlite3_backup_step(backupdb, -1);
		int counter = 0;
		while(counter < 10 && (rc == SQLITE_BUSY || rc == SQLITE_LOCKED))
		{
			counter++;
			rc = sqlite3_backup_step(backupdb, -1);
			sleep(1);
		}
		if(rc != SQLITE_OK && rc != SQLITE_DONE)
		{
			ss << "Database backup cannot be started to " << spec_sim_parameters->filename;
			ss << ". Check write access on output folder. Error code " << rc << ": ";
			ss << sqlite3_errmsg(outdatabase2) << endl;
			throw FatalException(ss.str());
		}
		rc = sqlite3_backup_finish(backupdb);
		if(rc != SQLITE_OK && rc != SQLITE_DONE)
		{
			ss << "Database backup cannot be completed to " << spec_sim_parameters->filename;
			ss << ". Check write access on output folder. Error code " << rc << ": ";
			ss << sqlite3_errmsg(outdatabase2) << endl;
			throw FatalException(ss.str());
		}
		sqlite3_close(outdatabase2);
		writeInfo("done!\n");
	}
	closeSqlConnection();
}

bool Community::checkSpeciesLocationsReference()
{
	if(!bSqlConnection)
	{
		throw FatalException("Attempted to get from sql database without opening database connection.");
	}

	sqlite3_stmt *stmt = nullptr;
	// Now find out the max size of the species_id_list, so we have a count to work from
	string count_command = "SELECT COUNT(*) FROM SPECIES_LOCATIONS WHERE community_reference == ";
	count_command += to_string(current_community_parameters->reference) + ";";
	sqlite3_prepare_v2(database, count_command.c_str(), static_cast<int>(strlen(count_command.c_str())), &stmt,
					   nullptr);
	sqlite3_step(stmt);
	int tmp_val = sqlite3_column_int(stmt, 0);
	// close the old statement
	sqlite3_finalize(stmt);
	return tmp_val > 0;
}

bool Community::checkSpeciesAbundancesReference()
{
	if(!bSqlConnection)
	{
		throw FatalException("Attempted to get from sql database without opening database connection.");
	}

	sqlite3_stmt *stmt = nullptr;
	// Now find out the max size of the species_id_list, so we have a count to work from
	string count_command = "SELECT COUNT(*) FROM SPECIES_ABUNDANCES WHERE community_reference = ";
	count_command += to_string(current_community_parameters->reference) + ";";
	sqlite3_prepare_v2(database, count_command.c_str(), static_cast<int>(strlen(count_command.c_str())), &stmt,
					   nullptr);
	sqlite3_step(stmt);
	int tmp_val = sqlite3_column_int(stmt, 0);
	// close the old statement
	sqlite3_finalize(stmt);
	return tmp_val > 0;
}

void Community::recordSpatial()
{
//	os << "Recording spatial data for speciation rate " << current_community_parameters->speciation_rate << "..." << flush;
	string table_command = "CREATE TABLE IF NOT EXISTS SPECIES_LOCATIONS (ID int PRIMARY KEY NOT NULL, species_id INT "
						   "NOT NULL, x INT NOT NULL, y INT NOT NULL, community_reference INT NOT NULL);";
	sqlite3_exec(database, table_command.c_str(), nullptr, nullptr, nullptr);
	getMaxSpeciesLocationsID();
	sqlite3_stmt *stmt = nullptr;
	// Checks that the SPECIES_LOCATIONS table doesn't already have a reference in matching the current reference
	if(current_community_parameters->updated)
	{
		if(checkSpeciesLocationsReference())
		{
			return;
		}
	}
	table_command = "INSERT INTO SPECIES_LOCATIONS (ID,species_id, x, y, community_reference) VALUES (?,?,?,?,?);";

	sqlite3_prepare_v2(database, table_command.c_str(), static_cast<int>(strlen(table_command.c_str())), &stmt,
					   nullptr);
	//		os << "test1" << endl;
	// Start the transaction
	sqlite3_exec(database, "BEGIN TRANSACTION;", nullptr, nullptr, nullptr);
	// Make sure only the tips which we want to check are recorded
	//		os << "nodes->size(): " << nodes->size() << endl;
	for(unsigned long i = 1; i < nodes->size(); i++)
	{
		TreeNode *this_node = &(*nodes)[i];
		//			os << nodes[i].exists() << endl;
		if(this_node->isTip() &&
		   this_node->exists() && doubleCompare(static_cast<double>(this_node->getGeneration()),
												static_cast<double>(current_community_parameters->time), 0.0001))
		{
			if(samplemask.getMaskVal(this_node->getXpos(), this_node->getYpos(),
									 this_node->getXwrap(), this_node->getYwrap()))
			{
				long x = this_node->getXpos();
				long y = this_node->getYpos();
				long xwrap = this_node->getXwrap();
				long ywrap = this_node->getYwrap();
				long xval = x + (xwrap * grid_x_size) + samplemask_x_offset;
				long yval = y + (ywrap * grid_y_size) + samplemask_y_offset;
				sqlite3_bind_int(stmt, 1, static_cast<int>(max_locations_id++));
				sqlite3_bind_int(stmt, 2, static_cast<int>(this_node->getSpeciesID()));
				sqlite3_bind_int(stmt, 3, static_cast<int>(xval));
				sqlite3_bind_int(stmt, 4, static_cast<int>(yval));
				sqlite3_bind_int(stmt, 5, static_cast<int>(current_community_parameters->reference));
				int step = sqlite3_step(stmt);
				// makes sure the while loop doesn't go forever.
				time_t start_check, end_check;
				time(&start_check);
				time(&end_check);
				while(step != SQLITE_DONE && (end_check - start_check) < 10 && step != SQLITE_OK)
				{
					step = sqlite3_step(stmt);
					time(&end_check);
				}
				if(step != SQLITE_DONE)
				{
					stringstream ss;
					ss << "SQLITE error code: " << step << endl;
					ss << "ERROR_SQL_004f: Could not insert into database. Check destination file has not "
						  "been moved or deleted and that an entry doesn't already exist with the same ID."
					   << endl;
					ss << sqlite3_errmsg(database) << endl;
					writeWarning(ss.str());
					break;
				}
				sqlite3_clear_bindings(stmt);
				sqlite3_reset(stmt);
			}
		}
	}
	// execute the command and close the connection to the database
	int rc1 = sqlite3_exec(database, "END TRANSACTION;", nullptr, nullptr, nullptr);
	// Need to finalise the statement
	int rc2 = sqlite3_finalize(stmt);
	if(rc1 != SQLITE_OK || rc2 != SQLITE_OK)
	{
		writeError("ERROR_SQL_013: Could not complete SQL transaction. Check memory database assignment and SQL "
				   "commands. Ensure SQL statements are properly cleared and that you are not attempting to insert "
				   "repeat IDs into the database.");
	}
}

void Community::calcFragments(string fragment_file)
{
	// Loop over every grid cell in the samplemask to determine if it is the start (top left corner) of a fragment.
	// Note that fragment detection only works for squares and rectangles. Adjacent squares and rectangles will be
	// treated as separate fragments if they are different sizes.
	// Downwards shapes are prioritised (i.e. a vertical rectangle on top of a horizontal rectangle will produce 3
	// fragments instead of two - this is a known bug).
	if(fragment_file == "null")
	{
		unsigned long fragment_number = 0;
		for(unsigned long i = 0; i < samplemask.sample_mask.getCols(); i++)
		{
			for(unsigned long j = 0; j < samplemask.sample_mask.getRows(); j++)
			{
				bool in_fragment = false;
				// Make sure is isn't on the top or left edge
				if(samplemask.sample_mask[j][i])
				{
					if(i > 0 && j > 0)
					{
						// Perform the check
						in_fragment = !(samplemask.sample_mask[j][i - 1] || samplemask.sample_mask[j - 1][i]);
					}
						// if it is on an edge, we need to check the fragment
					else
					{
						// if it is on the left edge we need to check above it - if there is forest
						// there, it is not a fragment.
						if(i == 0 && j > 0)
						{
							if(!samplemask.sample_mask[j - 1][i])
							{
								in_fragment = true;
							}
						}
							// if it is on the top edge, need to check to the left of it -  if there is
							// forest there, it is not a fragment.
						else if(j == 0 && i > 0)
						{
							if(!samplemask.sample_mask[j][i - 1])
							{
								in_fragment = true;
							}
						}
						else if(i == 0 && j == 0)
						{
							in_fragment = true;
						}
					}
				}
				if(in_fragment)
				{
					// Now move along the x and y axis (separately) until we hit a non-forest patch.
					// This marks the edge of the fragment and the value is recorded.
					bool x_continue = true;
					bool y_continue = true;
					unsigned long x, y;
					x = i;
					y = j;
					fragment_number++;
					// Also need to check that fragments that lie partly next to each other aren't
					// counted twice.
					// So count along the x axis until we hit non-habitat. Then count down the y axis
					// checking both extremes of the square for non-habitat.
					// Perform a check on the x axis to make sure that the square above is empty, as
					// fragments give priority in a downwards motion.
					while(x_continue)
					{
						x++;
						if(samplemask.sample_mask[j][x])
						{
							// Check we're not on top edge of the map.
							if(j > 0)
							{
								// if the cell above is non-fragment then we don't need to
								// continue (downwards fragments get priority).
								if(samplemask.sample_mask[j - 1][x])
								{
									x_continue = true;
								}
								else
								{
									x_continue = false;
								}
							}
							else
							{
								x_continue = true;
							}
						}
						else
						{
							x_continue = false;
						}
					}
					while(y_continue)
					{
						y++;
						// Make sure both extremes of the rectangle are still within patch.
						if(samplemask.sample_mask[y][i] && samplemask.sample_mask[y][i - 1])
						{
							y_continue = true;
						}
						else
						{
							y_continue = false;
						}
					}
					// Create the fragment to add.
					Fragment to_add;
					to_add.name = to_string((long long) fragment_number);
					to_add.x_west = i;
					to_add.x_east = x - 1;
					to_add.y_north = j;
					to_add.y_south = y - 1;
					// calculate the square area of the plot and record it.
					to_add.area = (x - i) * (y - j);
					// Now store the size of the fragment in the vector.
					fragments.push_back(to_add);
				}
			}
		}
	}
	else
	{
		stringstream os;
		os << "Importing fragments from " << fragment_file << endl;
		writeInfo(os.str());
#ifdef use_csv
		writeInfo("Using fast-cpp-csv-parse");
		// Then use the fast-cpp-csv-parser
		// There is a config file to import - here we use a specific piece of import code to parse the csv file.
		// first count the number of lines
		int number_of_lines = 0;
		string line;
		ifstream fragment_configs(fragment_file);
		while(getline(fragment_configs, line))
		{
			number_of_lines++;
		}
		//			os << "Number of lines in text file: " << number_of_lines << endl;
		fragment_configs.close();
		io::LineReader in(fragment_file);
		// Keep track of whether we've printed to terminal or not.
		bool bPrint = false;
		fragments.resize(number_of_lines);
//		os << "size: "  << fragments.capacity() << endl;
		for(int i = 0; i < number_of_lines; i++)
		{
			char *line = in.next_line();
			if(line == nullptr)
			{
				if(!bPrint)
				{
					writeError("Input dimensions incorrect - read past end of file.");
					bPrint = true;
				}
				break;
			}
			else
			{
				char *dToken;
				dToken = strtok(line, ",");
				for(int j = 0; j < 6; j++)
				{
					//						os << j << endl;
					if(dToken == nullptr)
					{
						if(!bPrint)
						{
							writeError("Input dimensions incorrect - read past end of file.");
							bPrint = true;
						}
						break;
					}
					else
					{
						//							os << "-" << endl;
						switch(j)
						{
							case 0:
								fragments[i].name = string(dToken);
								break;
							case 1:
								fragments[i].x_west = atoi(dToken);
								break;
							case 2:
								fragments[i].y_north = atoi(dToken);
								break;
							case 3:
								fragments[i].x_east = atoi(dToken);
								break;
							case 4:
								fragments[i].y_south = atoi(dToken);
								break;
							case 5:
								fragments[i].area = atof(dToken);
								break;
						}
						dToken = strtok(NULL, ",");
					}
				}
			}
		}
#endif
#ifndef use_csv
		ifstream fragment_configs(fragment_file);
		vector<vector<string>> tmp_raw_read;
		while(fragment_configs.good())
		{
			tmp_raw_read.push_back(getCsvLineAndSplitIntoTokens(fragment_configs));
		}
		fragments.resize(tmp_raw_read.size());
		for(unsigned long i = 0; i < tmp_raw_read.size(); i++)
		{
			vector<string> *this_fragment = &tmp_raw_read[i];
			if(this_fragment->size() != 6)
			{
				// Only throw an error if the size is not 1 (which usually means that there is an extra line
				// at the end of the file)
				if(this_fragment->size() != 1)
				{
					stringstream ss;
					ss << "Could not parse fragments file, " << this_fragment->size() << " columns detected";
					ss << ", requires 6 (name, x_west, y_north, x_east, y_south, area)." << endl;
					throw FatalException(ss.str());
				}
				fragments.resize(tmp_raw_read.size() - 1);
				break;
			}
			// Fragment name and dimensions
			try
			{
				for(auto &item : (*this_fragment))
				{
					if(item.empty())
					{
						throw invalid_argument("Cannot convert empty argument.");
					}
				}
				fragments[i].name = (*this_fragment)[0];
				fragments[i].x_west = stoi((*this_fragment)[1]);
				fragments[i].y_north = stoi((*this_fragment)[2]);
				fragments[i].x_east = stoi((*this_fragment)[3]);
				fragments[i].y_south = stoi((*this_fragment)[4]);
				fragments[i].area = stof((*this_fragment)[5]);
			}
			catch(invalid_argument &argument)
			{
				stringstream ss;
				ss << "Could not convert row arguments: ";
				for(auto &n : (*this_fragment))
				{
					ss << n << ", ";
				}
				ss << " should be str, int, int, int, int, float." << endl;
				throw FatalException(ss.str());
			}
		}
		fragment_configs.close();
#endif
	}
//	os << "Completed fragmentation analysis: " << fragments.size() << " fragments identified." << endl;
}

void Community::applyFragments()
{
	// For each fragment in the vector, perform the analysis and record the data in to a new data object, which will
	// then be outputted to an SQL file.
	for(unsigned int i = 0; i < fragments.size(); i++)
	{
		stringstream os;
		os << "\rApplying fragments... " << (i + 1) << "/" << fragments.size() << "      " << flush;
		writeInfo(os.str());
		// Set the new samplemask to the fragment
		samplemask.setFragment(fragments[i]);
		// Now filter only those lineages which exist in the fragments.
		// We also want to count the number of individuals that actually exist
		unsigned long iSpecCount = 0;
		for(unsigned long j = 0; j < nodes->size(); j++)
		{
			TreeNode *this_node = &(*nodes)[j];
			if(this_node->isTip() && samplemask.getMaskVal(this_node->getXpos(), this_node->getYpos(),
														   this_node->getXwrap(), this_node->getYwrap()) &&
			   doubleCompare(this_node->getGeneration(), current_community_parameters->time, 0.0001))
			{
				// if they exist exactly in the generation of interest.
				this_node->setExistence(true);
				iSpecCount++;
			}
			else if(this_node->isTip())
			{
				this_node->setExistence(false);
			}
		}
		fragments[i].num = iSpecCount;
		// Now calculate the species abundance. This will create a vector with lots of zeros in it. However, the
		// database creation will filter these out.
		calcSpeciesAbundance();
		createFragmentDatabase(fragments[i]);
		//			os << "done!" << endl;
	}
	samplemask.removeFragment();
	writeInfo("done!\n");
}

void Community::setSimParameters(const shared_ptr<SimParameters> sim_parameters)
{
	if(!has_imported_data)
	{
		min_spec_rate = sim_parameters->spec;
		grid_x_size = sim_parameters->grid_x_size;
		grid_y_size = sim_parameters->grid_y_size;
		protracted = sim_parameters->is_protracted;
		min_speciation_gen = sim_parameters->min_speciation_gen;
		max_speciation_gen = sim_parameters->max_speciation_gen;
		samplemask_x_offset = sim_parameters->sample_x_offset;
		samplemask_y_offset = sim_parameters->sample_y_offset;
		samplemask_x_size = sim_parameters->sample_x_size;
		samplemask_y_size = sim_parameters->sample_y_size;
		if(protracted)
		{
			if(max_speciation_gen == 0.0)
			{
				throw SpeciesException("Protracted speciation does not make sense when maximum speciation gen is 0.0.");
			}
			if(min_speciation_gen > max_speciation_gen)
			{
				throw SpeciesException("Cannot have simulation with minimum speciation generation less than maximum!");
			}
		}
	}
	has_imported_data = true;
}

void Community::importSimParameters(string file)
{
	if(has_imported_data)
	{
		return;
	}
	if(!bSqlConnection)
	{
#ifdef DEBUG
		stringstream os;
		os << "opening connection..." << flush;
#endif
		openSqlConnection(file);
#ifdef DEBUG
		os << "done!" << endl;
		writeInfo(os.str());
#endif
	}
	try
	{
#ifdef DEBUG
		stringstream os;
		os << "Reading current_metacommunity_parameters..." << flush;
#endif
		sqlite3_stmt *stmt2;
		string sql_parameters = "SELECT speciation_rate, grid_x, grid_y, protracted, min_speciation_gen, max_speciation_gen, "
								"sample_x_offset, sample_y_offset, sample_x, sample_y  FROM SIMULATION_PARAMETERS;";
		int rc = sqlite3_prepare_v2(database, sql_parameters.c_str(), static_cast<int>(strlen(sql_parameters.c_str())),
									&stmt2, nullptr);
		if(rc != SQLITE_DONE && rc != SQLITE_OK)
		{
			stringstream ss;
			ss << "Could not open simulation parameters in " << file << ": Error code: " << rc << ": ";
			ss << sqlite3_errmsg(database);
			sqlite3_close(database);
			throw SpeciesException(ss.str());
		}
		sqlite3_step(stmt2);
		min_spec_rate = sqlite3_column_double(stmt2, 0);
		grid_x_size = static_cast<unsigned long>(sqlite3_column_int(stmt2, 1));
		grid_y_size = static_cast<unsigned long>(sqlite3_column_int(stmt2, 2));
		protracted = bool(sqlite3_column_int(stmt2, 3));
		min_speciation_gen = sqlite3_column_double(stmt2, 4);
		max_speciation_gen = sqlite3_column_double(stmt2, 5);
		samplemask_x_offset = static_cast<unsigned long>(sqlite3_column_int(stmt2, 6));
		samplemask_y_offset = static_cast<unsigned long>(sqlite3_column_int(stmt2, 7));
		samplemask_x_size = static_cast<unsigned long>(sqlite3_column_int(stmt2, 8));
		samplemask_y_size = static_cast<unsigned long>(sqlite3_column_int(stmt2, 9));
		if(protracted)
		{
			if(max_speciation_gen == 0.0)
			{
				throw SpeciesException("Protracted speciation does not make sense when maximum speciation gen is 0.0.");
			}
			if(min_speciation_gen > max_speciation_gen)
			{
				throw SpeciesException("Cannot have simulation with minimum speciation generation less than maximum!");
			}
		}
		sqlite3_step(stmt2);
		sqlite3_finalize(stmt2);
#ifdef DEBUG
		os << "done!" << endl;
		writeInfo(os.str());
#endif
	}
	catch(exception &er)
	{
		throw SpeciesException(er.what());
	}
	has_imported_data = true;
}

void Community::forceSimCompleteParameter()
{
	if(!bSqlConnection)
	{
		openSqlConnection(spec_sim_parameters->filename);
	}
	string update_command = "UPDATE SIMULATION_PARAMETERS SET sim_complete=1 WHERE sim_complete=0;";
	sqlite3_exec(database, update_command.c_str(), nullptr, nullptr, nullptr);
}

bool Community::isSetDatabase()
{
	return database_set;
}

void Community::setProtractedParameters(const ProtractedSpeciationParameters &protracted_params)
{
	applied_protracted_parameters = protracted_params;
	if(min_speciation_gen > 0 && max_speciation_gen > 0 &&
	   (applied_protracted_parameters.min_speciation_gen > min_speciation_gen ||
		applied_protracted_parameters.max_speciation_gen > max_speciation_gen))
	{
#ifdef DEBUG
		writeLog(50, "Applied speciation current_metacommunity_parameters: " + to_string(applied_protracted_parameters.min_speciation_gen) + ", " +
				to_string(applied_protracted_parameters.max_speciation_gen));
		writeLog(50, "Simulated speciation current_metacommunity_parameters: " + to_string(min_speciation_gen) + ", " +
					 to_string(max_speciation_gen));
#endif // DEBUG
		stringstream ss;
		ss << "Applied protracted speciation current_metacommunity_parameters: "
		   << applied_protracted_parameters.min_speciation_gen << ", ";
		ss << applied_protracted_parameters.max_speciation_gen << endl;
		ss << "Original protracted speciation current_metacommunity_parameters: " << min_speciation_gen << ", "
		   << max_speciation_gen << endl;
		writeCritical(ss.str());
		throw SpeciesException(
				"Cannot use protracted current_metacommunity_parameters with minimum > simulated minimum or "
				"maximum > simulated maximums.");
	}
}

void Community::overrideProtractedParameters(const ProtractedSpeciationParameters &protracted_params)
{
	min_speciation_gen = protracted_params.min_speciation_gen;
	max_speciation_gen = protracted_params.max_speciation_gen;
	applied_protracted_parameters = protracted_params;

}

void Community::setProtracted(bool protracted_in)
{
	protracted = protracted_in;
}

void Community::getPreviousCalcs()
{
	writeInfo("Getting previous calculations...");
	// Read the community current_metacommunity_parameters and store them in the relevant objects
	sqlite3_stmt *stmt1;
	string call1 = "select count(type) from sqlite_master where type='table' and name='COMMUNITY_PARAMETERS'";
	int rc = sqlite3_prepare_v2(database, call1.c_str(), static_cast<int>(strlen(call1.c_str())), &stmt1, nullptr);
	if(rc != SQLITE_DONE && rc != SQLITE_OK)
	{
		stringstream ss;
		ss << "Could not detect COMMUNITY_PARAMETERS table while finding previous calculations. Error code: " << rc
		   << ": ";
		ss << sqlite3_errmsg(database) << endl;
		sqlite3_close(database);
		throw FatalException(ss.str());
	}
	sqlite3_step(stmt1);
	auto has_community_parameters = static_cast<bool>(sqlite3_column_int(stmt1, 0));
	sqlite3_step(stmt1);
	sqlite3_finalize(stmt1);
	// Read the speciation rates from the community_parameters table
	if(has_community_parameters)
	{
		sqlite3_stmt *stmt2;
		string call2 = "SELECT reference, speciation_rate, time, fragments, metacommunity_reference ";
		if(protracted)
		{
			call2 += ", min_speciation_gen, max_speciation_gen ";
		}
		call2 += " FROM COMMUNITY_PARAMETERS";
		rc = sqlite3_prepare_v2(database, call2.c_str(), static_cast<int>(strlen(call2.c_str())), &stmt2,
								nullptr);
		if(rc != SQLITE_DONE && rc != SQLITE_OK)
		{
			stringstream ss;
			ss << "Could not select from COMMUNITY_PARAMETERS table. Error code: " << rc << ": ";
			ss << sqlite3_errmsg(database) << endl;
			sqlite3_close(database);
			throw SpeciesException(ss.str());
		}
		rc = sqlite3_step(stmt2);
		if(rc == SQLITE_ROW)
		{
			writeInfo("previous calculations detected.\n");
		}
		else
		{
			writeInfo("COMMUNITY_PARAMETERS table detected, but not previous calculations found.\n");
		}
		while(rc == SQLITE_ROW)
		{
			auto row_val = sqlite3_column_int(stmt2, 0);
			if(row_val == 0)
			{
				writeWarning(
						"Reference of 0 found in community current_metacommunity_parameters in database, skipping...\n");
			}
			else
			{
				ProtractedSpeciationParameters tmp{};
				if(protracted)
				{
					tmp.min_speciation_gen = sqlite3_column_double(stmt2, 5);
					tmp.max_speciation_gen = sqlite3_column_double(stmt2, 6);
				}
				past_communities.pushBack(static_cast<unsigned long>(row_val), sqlite3_column_double(stmt2, 1),
										  sqlite3_column_double(stmt2, 2), bool(sqlite3_column_int(stmt2, 3)),
										  static_cast<unsigned long>(sqlite3_column_int(stmt2, 4)), tmp);
			}
			rc = sqlite3_step(stmt2);
		}
		if(rc != SQLITE_OK && rc != SQLITE_DONE)
		{
			stringstream ss;
			ss << "ERROR_SQL_020b: FATAL. Could not read community current_metacommunity_parameters." << endl;
			ss << "Code: " << rc << endl << "Errmsg: ";
			ss << sqlite3_errmsg(database) << endl;
			sqlite3_clear_bindings(stmt2);
			sqlite3_reset(stmt2);
			throw SpeciesException(ss.str());
		}
		sqlite3_finalize(stmt2);
	}
	else
	{
		writeInfo("no previous calculations detected.\n");
	}
	// And the same for metacommunity current_metacommunity_parameters
	sqlite3_stmt *stmt3;
	string call3 = "select count(type) from sqlite_master where type='table' and name='METACOMMUNITY_PARAMETERS'";
	rc = sqlite3_prepare_v2(database, call3.c_str(), static_cast<int>(strlen(call3.c_str())), &stmt3, nullptr);
	if(rc != SQLITE_DONE && rc != SQLITE_OK)
	{
		sqlite3_close(database);
		throw SpeciesException(
				"ERROR_SQL_020: FATAL. Could not check for METACOMMUNITY_PARAMETERS table. Error code: " +
				to_string(rc));
	}
	sqlite3_step(stmt3);
	has_community_parameters = static_cast<bool>(sqlite3_column_int(stmt3, 0));
	sqlite3_step(stmt3);
	sqlite3_finalize(stmt3);
	// Read the speciation rates from the community_parameters table
	if(has_community_parameters)
	{
		sqlite3_stmt *stmt4;
		string call4 = "SELECT reference, speciation_rate, metacommunity_size, option, external_reference FROM ";
		call4 += "METACOMMUNITY_PARAMETERS";
		rc = sqlite3_prepare_v2(database, call4.c_str(), static_cast<int>(strlen(call4.c_str())), &stmt4,
								nullptr);
		if(rc != SQLITE_DONE && rc != SQLITE_OK)
		{
			sqlite3_close(database);
			throw SpeciesException(
					"ERROR_SQL_020: FATAL. Could not detect METACOMMUNITY_PARAMETERS table. Error code: " +
					to_string(rc));
		}
		rc = sqlite3_step(stmt4);
		while(rc == SQLITE_ROW)
		{
			past_metacommunities.pushBack(static_cast<unsigned long>(sqlite3_column_int(stmt4, 0)),
										  static_cast<unsigned long>(sqlite3_column_int(stmt4, 2)),
										  sqlite3_column_double(stmt4, 1), (char *) (sqlite3_column_text(stmt4, 3)),
										  static_cast<const unsigned long &>(sqlite3_column_int(stmt4, 4)));
			rc = sqlite3_step(stmt4);
		}
		if(rc != SQLITE_OK && rc != SQLITE_DONE)
		{
			stringstream ss;
			ss << "ERROR_SQL_020: FATAL. Could not read metacommunity current_metacommunity_parameters." << endl;
			ss << "Code: " << rc << endl << "Errmsg: ";
			ss << sqlite3_errmsg(database) << endl;
			sqlite3_clear_bindings(stmt4);
			sqlite3_reset(stmt4);
			throw SpeciesException(ss.str());
		}
		sqlite3_step(stmt4);
		sqlite3_finalize(stmt4);
	}
}

void Community::addCalculationPerformed(const long double &speciation_rate, const double &time, const bool &fragments,
										const MetacommunityParameters &metacomm_parameters,
										const ProtractedSpeciationParameters &protracted_parameters)
{
	auto meta_reference = past_metacommunities.getReference(metacomm_parameters);
	if(meta_reference == 0 && metacomm_parameters.isMetacommunityOption())
	{
#ifdef DEBUG
		stringstream ss;
		ss << "Adding metacommunity (" << metacomm_parameters.metacommunity_size << ", " <<
		   metacomm_parameters.speciation_rate << ")" << endl;
		writeInfo(ss.str());
#endif
		meta_reference = past_metacommunities.addNew(metacomm_parameters);
	}
	current_community_parameters = past_communities.addNew(speciation_rate, time, fragments, meta_reference,
														   protracted_parameters);
#ifdef DEBUG
	for(const auto &i : past_communities.comm_parameters)
	{
		if(doubleCompare(i->time, current_community_parameters->time, 0.00001) &&
			doubleCompare(i->speciation_rate, current_community_parameters->speciation_rate,
						  i->speciation_rate*0.00001) &&
				i->protracted_parameters == current_community_parameters->protracted_parameters &&
				i->metacommunity_reference == current_community_parameters->metacommunity_reference &&
				i->reference != current_community_parameters->reference)
		{
			throw FatalException("Communities are identical, but references differ! Please report this bug.");
		}
	}
#endif // DEBUG
}

vector<unsigned long> Community::getUniqueCommunityRefs()
{
	vector<unsigned long> unique_community_refs;
	// Read the community parameters and store them in the relevant objects
	sqlite3_stmt *stmt1;
	string call1 = "select count(type) from sqlite_master where type='table' and name='COMMUNITY_PARAMETERS'";
	int rc = sqlite3_prepare_v2(database, call1.c_str(), static_cast<int>(strlen(call1.c_str())), &stmt1, nullptr);
	if(rc != SQLITE_DONE && rc != SQLITE_OK)
	{
		stringstream ss;
		ss << "Could not detect COMMUNITY_PARAMETERS table while getting unique community references. Error code: ";
		ss << rc << ": " << sqlite3_errmsg(database) << endl;
		sqlite3_close(database);
		throw FatalException(ss.str());
	}
	sqlite3_step(stmt1);
	auto has_community_parameters = static_cast<bool>(sqlite3_column_int(stmt1, 0));
	sqlite3_step(stmt1);
	sqlite3_finalize(stmt1);
	// Read the speciation rates from the community_parameters table
	if(has_community_parameters)
	{
		sqlite3_stmt *stmt2;
		string call2 = "SELECT DISTINCT(reference) FROM COMMUNITY_PARAMETERS";
		rc = sqlite3_prepare_v2(database, call2.c_str(), static_cast<int>(strlen(call2.c_str())), &stmt2,
								nullptr);
		if(rc != SQLITE_DONE && rc != SQLITE_OK)
		{
			stringstream ss;
			ss << "Could not get distinct references from COMMUNITY_PARAMETERS table. Error code: " << rc << ": ";
			ss << sqlite3_errmsg(database) << endl;
			sqlite3_close(database);
			throw FatalException(ss.str());
		}
		rc = sqlite3_step(stmt2);
		while(rc != SQLITE_DONE)
		{
			unique_community_refs.push_back(static_cast<unsigned long>(sqlite3_column_int(stmt2, 0)));
			rc = sqlite3_step(stmt2);
			if(rc > 10000)
			{
				throw SpeciesException("Could not read speciation rates.");
			}
		}
		sqlite3_step(stmt2);
		sqlite3_finalize(stmt2);
	}
	return unique_community_refs;
}

vector<unsigned long> Community::getUniqueMetacommunityRefs()
{
	vector<unsigned long> unique_metacommunity_refs;
	// Read the community current_metacommunity_parameters and store them in the relevant objects
	sqlite3_stmt *stmt1;
	string call1 = "select count(type) from sqlite_master where type='table' and name='METACOMMUNITY_PARAMETERS'";
	int rc = sqlite3_prepare_v2(database, call1.c_str(), static_cast<int>(strlen(call1.c_str())), &stmt1, nullptr);
	if(rc != SQLITE_DONE && rc != SQLITE_OK)
	{
		sqlite3_close(database);
		throw SpeciesException(
				"ERROR_SQL_020: FATAL. Could not check for METACOMMUNITY_PARAMETERS table. Error code: " +
				to_string(rc));
	}
	sqlite3_step(stmt1);
	auto has_metacommunity_parameters = static_cast<bool>(sqlite3_column_int(stmt1, 0));
	sqlite3_step(stmt1);
	sqlite3_finalize(stmt1);
	// Read the speciation rates from the community_parameters table
	if(has_metacommunity_parameters)
	{
		sqlite3_stmt *stmt2;
		string call2 = "SELECT DISTINCT(reference) FROM METACOMMUNITY_PARAMETERS";
		rc = sqlite3_prepare_v2(database, call2.c_str(), static_cast<int>(strlen(call2.c_str())), &stmt2,
								nullptr);
		if(rc != SQLITE_DONE && rc != SQLITE_OK)
		{
			sqlite3_close(database);
			throw SpeciesException(
					"ERROR_SQL_020: FATAL. Could not detect METACOMMUNITY_PARAMETERS table. Error code: " +
					to_string(rc));
		}
		rc = sqlite3_step(stmt2);
		while(rc != SQLITE_DONE)
		{
			unique_metacommunity_refs.push_back(static_cast<unsigned long>(sqlite3_column_int(stmt2, 0)));
			rc = sqlite3_step(stmt2);
			if(rc > 10000)
			{
				throw SpeciesException("ERROR_SQL_020: FATAL. Could not read speciation rates.");
			}
		}
		sqlite3_step(stmt2);
		sqlite3_finalize(stmt2);
	}
	return unique_metacommunity_refs;
}

void Community::writeNewCommunityParameters()
{
	// Find new community current_metacommunity_parameters to add
	auto unique_community_refs = getUniqueCommunityRefs();
	CommunitiesArray communities_to_write;
	for(auto &community_param : past_communities.comm_parameters)
	{
		if(find(unique_community_refs.begin(),
				unique_community_refs.end(), community_param->reference) == unique_community_refs.end())
		{
			communities_to_write.pushBack(community_param);
			unique_community_refs.push_back(community_param->reference);
		}
	}
	if(!communities_to_write.comm_parameters.empty())
	{
		// Create the table if it doesn't exist
		string table_command = "CREATE TABLE IF NOT EXISTS COMMUNITY_PARAMETERS (reference INT PRIMARY KEY NOT NULL,"
							   " speciation_rate DOUBLE NOT NULL, time DOUBLE NOT NULL, fragments INT NOT NULL, "
							   "metacommunity_reference INT";
		string table_command2 = "INSERT INTO COMMUNITY_PARAMETERS (reference, speciation_rate, time, fragments,"
								" metacommunity_reference";
		string table_command3 = "VALUES (?,?,?,?,?";
		if(protracted)
		{
			table_command += ", min_speciation_gen DOUBLE NOT NULL, max_speciation_gen DOUBLE NOT NULL";
			table_command2 += ", min_speciation_gen, max_speciation_gen";
			table_command3 += ", ?, ?";
		}
		table_command += ");";
		table_command2 += ") " + table_command3 + ");";

		sqlite3_exec(database, table_command.c_str(), nullptr, nullptr, nullptr);
		sqlite3_stmt *stmt = nullptr;

		sqlite3_prepare_v2(database, table_command2.c_str(), static_cast<int>(strlen(table_command.c_str())), &stmt,
						   nullptr);
		// Then add the required elements
		sqlite3_exec(database, "BEGIN TRANSACTION;", nullptr, nullptr, nullptr);
		for(auto &item : communities_to_write.comm_parameters)
		{
			if(item->reference == 0)
			{
				continue;
			}
			sqlite3_bind_int(stmt, 1, static_cast<int>(item->reference));
			sqlite3_bind_double(stmt, 2, static_cast<double>(item->speciation_rate));
			sqlite3_bind_double(stmt, 3, static_cast<double>(item->time));
			sqlite3_bind_int(stmt, 4, static_cast<int>(item->fragment));
			sqlite3_bind_int(stmt, 5, static_cast<int>(item->metacommunity_reference));
			if(protracted)
			{
				sqlite3_bind_double(stmt, 6, item->protracted_parameters.min_speciation_gen);
				sqlite3_bind_double(stmt, 7, item->protracted_parameters.max_speciation_gen);
			}
			time_t start_check, end_check;
			time(&start_check);
			time(&end_check);
			int step = sqlite3_step(stmt);
			while(step != SQLITE_DONE && (end_check - start_check) < 10 && step != SQLITE_OK)
			{
				step = sqlite3_step(stmt);
				time(&end_check);
			}
			if(step != SQLITE_DONE)
			{
				stringstream ss;
				ss << "SQLITE error code: " << step << endl;
				ss << sqlite3_errmsg(database) << endl;
				ss << "ERROR_SQL_004a: Could not insert into database. Check destination file has not "
					  "been moved or deleted and that an entry doesn't already exist with the same ID."
				   << endl;
				sqlite3_clear_bindings(stmt);
				sqlite3_reset(stmt);
				writeWarning(ss.str());
				break;
			}
			sqlite3_clear_bindings(stmt);
			sqlite3_reset(stmt);
		}
		int rc1 = sqlite3_exec(database, "END TRANSACTION;", nullptr, nullptr, nullptr);
		// Need to finalise the statement
		int rc2 = sqlite3_finalize(stmt);
		if(rc1 != SQLITE_OK || rc2 != SQLITE_OK)
		{
			stringstream ss;
			ss << "ERROR_SQL_013: Could not complete SQL transaction. Check memory database assignment and SQL "
				  "commands. Please report this bug." << endl;
			writeWarning(ss.str());
		}
	}
}

void Community::writeNewMetacommunityParameters()
{
	auto unique_metacommunity_refs = getUniqueMetacommunityRefs();
	MetacommunitiesArray metacommunities_to_write;
	if(unique_metacommunity_refs.empty())
	{
		for(auto &community_param : past_metacommunities.metacomm_parameters)
		{
			metacommunities_to_write.pushBack(community_param);
		}
	}
	else
	{
		for(auto &community_param : past_metacommunities.metacomm_parameters)
		{
			if(find(unique_metacommunity_refs.begin(),
					unique_metacommunity_refs.end(), community_param->reference) == unique_metacommunity_refs.end())
			{
				metacommunities_to_write.pushBack(community_param);
				unique_metacommunity_refs.push_back(community_param->reference);
			}
		}
	}
	if(!metacommunities_to_write.metacomm_parameters.empty())
	{
		// Create the table if it doesn't exist
		string table_command = "CREATE TABLE IF NOT EXISTS METACOMMUNITY_PARAMETERS (reference INT PRIMARY KEY NOT NULL,"
							   " speciation_rate DOUBLE NOT NULL, metacommunity_size DOUBLE NOT NULL, "
							   "option TEXT NOT NULL, external_reference INT NOT NULL);";
		sqlite3_exec(database, table_command.c_str(), nullptr, nullptr, nullptr);
		sqlite3_stmt *stmt = nullptr;
		table_command = "INSERT INTO METACOMMUNITY_PARAMETERS (reference, speciation_rate, metacommunity_size, "
						"option, external_reference) VALUES (?,?,?, ?, ?);";
		sqlite3_prepare_v2(database, table_command.c_str(), static_cast<int>(strlen(table_command.c_str())), &stmt,
						   nullptr);
		// Then add the required elements
		sqlite3_exec(database, "BEGIN TRANSACTION;", nullptr, nullptr, nullptr);
		for(auto &item : metacommunities_to_write.metacomm_parameters)
		{
			if(item->reference == 0)
			{
				continue;
			}
			sqlite3_bind_int(stmt, 1, static_cast<int>(item->reference));
			sqlite3_bind_double(stmt, 2, static_cast<double>(item->speciation_rate));
			sqlite3_bind_int(stmt, 3, static_cast<int>(item->metacommunity_size));
			sqlite3_bind_text(stmt, 4, item->option.c_str(), static_cast<int>(item->option.length()), SQLITE_TRANSIENT);
			sqlite3_bind_int(stmt, 5, static_cast<int>(item->external_reference));
			time_t start_check, end_check;
			time(&start_check);
			time(&end_check);
			int step = sqlite3_step(stmt);
			while(step != SQLITE_DONE && (end_check - start_check) < 10 && step != SQLITE_OK)
			{
				step = sqlite3_step(stmt);
				time(&end_check);
			}
			if(step != SQLITE_DONE)
			{
#ifdef DEBUG
				stringstream ss;
				ss << "SQLITE error code: " << step << endl;
				ss << "Metacommunity reference: " << item->reference << endl;
				ss << "Speciation rate: " << item->speciation_rate << ", metacommunity size: " << item->metacommunity_size << endl;
				ss << sqlite3_errmsg(database) << endl;
				writeLog(10, ss);
#endif // DEBUG
				throw SpeciesException("ERROR_SQL_004b: Could not insert into database. Check destination file has not "
									   "been moved or deleted and that an entry doesn't already exist with the"
									   " same ID.");
			}
			sqlite3_clear_bindings(stmt);
			sqlite3_reset(stmt);
		}
		int rc1 = sqlite3_exec(database, "END TRANSACTION;", nullptr, nullptr, nullptr);
		// Need to finalise the statement
		int rc2 = sqlite3_finalize(stmt);
		if(rc1 != SQLITE_OK || rc2 != SQLITE_OK)
		{
			stringstream ss;
			ss << "ERROR_SQL_013: Could not complete SQL transaction. Check memory database assignment and SQL "
				  "commands. Please report this bug." << endl;
			ss << sqlite3_errmsg(database) << endl;
			writeWarning(ss.str());
		}
	}
}

void Community::createSpeciesList()
{
	string create_species_list;
	create_species_list =
			"CREATE TABLE SPECIES_LIST (ID int PRIMARY KEY NOT NULL, unique_spec INT NOT NULL, xval INT NOT NULL,";
	create_species_list += "yval INT NOT NULL, xwrap INT NOT NULL, ywrap INT NOT NULL, tip INT NOT NULL, speciated INT NOT "
						   "NULL, parent INT NOT NULL, existence INT NOT NULL, randnum DOUBLE NOT NULL, gen_alive INT NOT "
						   "NULL, gen_added DOUBLE NOT NULL);";

	// Create the table within the SQL database
	char *sErrMsg = nullptr;
	int rc = sqlite3_exec(database, create_species_list.c_str(), nullptr, nullptr, &sErrMsg);
	if(rc != SQLITE_OK)
	{
		stringstream ss;
		ss << "Error creating SPECIES_LIST table in database: " << sErrMsg << endl;
		throw FatalException(ss.str());
	}
}

void Community::deleteSpeciesList()
{
	string wipe_species_list;
	wipe_species_list =
			"DROP TABLE IF EXISTS SPECIES_LIST;";
	// Drop the table from the SQL database
	char *sErrMsg = nullptr;
	int rc = sqlite3_exec(database, wipe_species_list.c_str(), nullptr, nullptr, &sErrMsg);
	if(rc != SQLITE_OK)
	{
		stringstream ss;
		ss << "Error dropping SPECIES_LIST table in database: " << sErrMsg << endl;
		throw FatalException(ss.str());
	}
}

void Community::writeSpeciesList(const unsigned long &enddata)
{
	sqlite3_stmt *stmt = nullptr;
	char *sErrMsg = nullptr;
	// Now create the prepared statement into which we shall insert the values from the table
	string insert_species_list = "INSERT INTO SPECIES_LIST "
								 "(ID,unique_spec,xval,yval,xwrap,ywrap,tip,speciated,parent,existence,randnum,gen_alive,gen_added) "
								 "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)";
	sqlite3_prepare_v2(database, insert_species_list.c_str(), static_cast<int>(strlen(insert_species_list.c_str())),
					   &stmt, nullptr);

	// Start the transaction
	int rc = sqlite3_exec(database, "BEGIN TRANSACTION;", nullptr, nullptr, &sErrMsg);
	if(rc != SQLITE_OK)
	{
		writeError("ERROR_SQL_008: Cannot start SQL transaction. Check memory database assignment and SQL commands.");
	}
	for(unsigned int i = 0; i <= enddata; i++)
	{
		sqlite3_bind_int(stmt, 1, i);
		sqlite3_bind_int(stmt, 2, static_cast<int>((*nodes)[i].getSpeciesID()));
		sqlite3_bind_int(stmt, 3, static_cast<int>((*nodes)[i].getXpos()));
		sqlite3_bind_int(stmt, 4, static_cast<int>((*nodes)[i].getYpos()));
		sqlite3_bind_int(stmt, 5, static_cast<int>((*nodes)[i].getXwrap()));
		sqlite3_bind_int(stmt, 6, static_cast<int>((*nodes)[i].getYwrap()));
		sqlite3_bind_int(stmt, 7, (*nodes)[i].isTip());
		sqlite3_bind_int(stmt, 8, (*nodes)[i].hasSpeciated());
		sqlite3_bind_int(stmt, 9, static_cast<int>((*nodes)[i].getParent()));
		sqlite3_bind_int(stmt, 10, (*nodes)[i].exists());
		sqlite3_bind_double(stmt, 11, static_cast<double>((*nodes)[i].getSpecRate()));
		sqlite3_bind_int(stmt, 12, static_cast<int>((*nodes)[i].getGenRate()));
		sqlite3_bind_double(stmt, 13, static_cast<double>((*nodes)[i].getGeneration()));
		sqlite3_step(stmt);
		sqlite3_clear_bindings(stmt);
		sqlite3_reset(stmt);
	}
	stringstream os;
	os << "\r    Executing SQL commands...." << flush;
	writeInfo(os.str());
	// execute the command and close the connection to the database
	rc = sqlite3_exec(database, "END TRANSACTION;", nullptr, nullptr, &sErrMsg);
	if(rc != SQLITE_OK)
	{
		stringstream ss;
		ss << "ERROR_SQL_008: Cannot complete SQL transaction. Check memory database assignment and SQL "
			  "commands. Ensure SQL statements are properly cleared."
		   << endl;
		ss << "Error code: " << rc << endl;
		// try again
		int i = 0;
		while((rc != SQLITE_OK && rc != SQLITE_DONE) && i < 10)
		{
			sleep(1);
			i++;
			rc = sqlite3_exec(database, "END TRANSACTION;", nullptr, nullptr, &sErrMsg);
			ss << "Attempt " << i << " failed..." << endl;
			ss << "ERROR_SQL_008: Cannot complete SQL transaction. Check memory database assignment and SQL "
				  "commands. Ensure SQL statements are properly cleared." << endl;
		}
		writeError(ss.str());
	}
	// Need to finalise the statement
	rc = sqlite3_finalize(stmt);
	if(rc != SQLITE_OK)
	{
		stringstream ss;
		ss << "ERROR_SQL_008: Cannot complete SQL transaction. Check memory database assignment and SQL "
			  "commands. Ensure SQL statements are properly cleared."
		   << endl;
		ss << "Error code: " << rc << endl;
	}

}

void Community::updateCommunityParameters()
{
	for(const auto &parameter : past_communities.comm_parameters)
	{
		if(parameter->updated)
		{
			if(!bSqlConnection)
			{
				throw FatalException("Attempted to update sql database without opening database connection.");
			}

			// Now find out the max size of the species_id_list, so we have a count to work from
			string count_command = "UPDATE COMMUNITY_PARAMETERS SET fragments = 1 WHERE reference = ";
			count_command += to_string(parameter->reference) + ";";
			int rc = sqlite3_exec(database, count_command.c_str(), nullptr, nullptr, nullptr);
			// Need to finalise the statement
			if(rc != SQLITE_OK && rc != SQLITE_DONE)
			{
				stringstream ss;
				ss << "ERROR_SQL_013: Could not update sql database. Check file write access. ";
				ss << "Otherwise, please report this bug." << endl;
				ss << sqlite3_errmsg(database) << endl;
				writeWarning(ss.str());
			}
		}
	}
}

void Community::writeSpeciationRates()
{
	stringstream os;
	os << "***************************" << endl;
	os << "STARTING CALCULATIONS" << endl;
	os << "Input file is " << spec_sim_parameters->filename << endl;
	if(!spec_sim_parameters->bMultiRun)
	{
		os << "Speciation rate is " << *spec_sim_parameters->all_speciation_rates.begin() << endl;
	}
	else
	{
		os << "Speciation rates are: " << flush;
		unsigned long i = 0;
		for(const auto &item :spec_sim_parameters->all_speciation_rates)
		{
			os << item << flush;
			i++;
			if(i == spec_sim_parameters->all_speciation_rates.size())
			{
				os << "." << endl;
			}
			else
			{
				os << ", " << flush;
			}
		}
	}
	writeInfo(os.str());
	if(!spec_sim_parameters->protracted_parameters.empty())
	{
		os.str("");
		os << "Protracted speciation parameters (min, max) are: " << endl;
		for(const auto i : spec_sim_parameters->protracted_parameters)
		{
			os << i.min_speciation_gen << ", " << i.max_speciation_gen << endl;
		}
		writeInfo(os.str());
	}
	os.str("");
	for(const auto &item: spec_sim_parameters->metacommunity_parameters)
	{
		if(item->metacommunity_size > 0)
		{
			os << "Metacommunity size: " << item->metacommunity_size << endl;
			os << "Metacommunity speciation rate: " << item->speciation_rate << endl;
		}
	}
	writeInfo(os.str());
}

void Community::calculateTree()
{
	stringstream os;
	for(const auto &protracted_params : spec_sim_parameters->protracted_parameters)
	{
		setProtractedParameters(protracted_params);
		for(const auto &sr : spec_sim_parameters->all_speciation_rates)
		{
			os << "Calculating speciation rate " << sr << endl;
			writeInfo(os.str());
			os.str("");
			for(auto time : spec_sim_parameters->all_times)
			{
				os.str("");
				os << "Calculating generation " << time << "\n";
				writeInfo(os.str());
				resetTree();
				if(!checkCalculationsPerformed(sr, time, spec_sim_parameters->use_fragments,
											   *current_metacommunity_parameters, applied_protracted_parameters))
				{
					addCalculationPerformed(sr, time, spec_sim_parameters->use_fragments,
											*current_metacommunity_parameters, applied_protracted_parameters);
					createDatabase();
					if(spec_sim_parameters->use_spatial)
					{
						recordSpatial();
					}
					if(spec_sim_parameters->use_fragments)
					{
						applyFragments();
					}
				}
				else
				{
					os.str("");
					os << "calculation already performed for speciation rate=" << sr << ", time=" << time;
					os << " and protracted parameters " << protracted_params.min_speciation_gen << ", ";
					os << protracted_params.max_speciation_gen << endl;
					writeInfo(os.str());
				}
			}
		}
	}
}

void Community::output()
{
	writeNewCommunityParameters();
	writeNewMetacommunityParameters();
	updateCommunityParameters();
	exportDatabase();
}

void Community::printEndTimes(time_t tStart, time_t tEnd)
{
	time(&tEnd);
	stringstream os;
	os << "Calculations complete." << endl;
	os << "Time taken was " << floor((tEnd - tStart) / 3600) << " hours "
	   << (floor((tEnd - tStart) / 60) - 60 * floor((tEnd - tStart) / 3600)) << " minutes " << (tEnd - tStart) % 60
	   << " seconds" << endl;
	writeInfo(os.str());
}

void Community::apply(shared_ptr<SpecSimParameters> sp)
{
	time_t tStart{};
	time_t tEnd{};
	// Start the clock
	time(&tStart);
	applyNoOutput(std::move(sp));
	output();
	printEndTimes(tStart, tEnd);
}

void Community::applyNoOutput(shared_ptr<SpecSimParameters> sp)
{
	doApplication(std::move(sp));
}

void Community::doApplication(shared_ptr<SpecSimParameters> sp)
{
	shared_ptr<vector<TreeNode>> data = make_shared<vector<TreeNode>>();
	doApplication(std::move(sp), data);
}

void Community::doApplication(shared_ptr<SpecSimParameters> sp, shared_ptr<vector<TreeNode>> data)
{
	spec_sim_parameters = sp;
	writeSpeciationRates();
	// Set up the objects
	setList(std::move(data));
	importSimParameters(sp->filename);
	importSamplemask(sp->samplemask);
	importData(sp->filename);
	getPreviousCalcs();
	if(sp->use_fragments)
	{
		calcFragments(sp->fragment_config_file);
		stringstream os;
		os << "Total fragments: " << fragments.size() << endl;
		writeInfo(os.str());
	}
	if(spec_sim_parameters->metacommunity_parameters.empty())
	{
		spec_sim_parameters->metacommunity_parameters.addNull();
		current_metacommunity_parameters = spec_sim_parameters->metacommunity_parameters.metacomm_parameters[0];
	}
	if(spec_sim_parameters->protracted_parameters.empty())
	{
		ProtractedSpeciationParameters tmp;
		spec_sim_parameters->protracted_parameters.emplace_back(tmp);
	}

	calculateTree();
}

void Community::doApplicationInternal(shared_ptr<SpecSimParameters> sp, shared_ptr<vector<TreeNode>> data)
{
	setInternalDatabase();
	doApplication(std::move(sp), std::move(data));
}

void Community::speciateRemainingLineages(const string &filename)
{
	importSimParameters(filename);
	importSamplemask("null");
	importData(filename);
	spec_sim_parameters->filename = filename;
	// Skip the first entry as it's always blank
	for(unsigned long i = 1; i < nodes->size(); i++)
	{
		TreeNode *this_node = &(*nodes)[i];
		if(this_node->getParent() == 0 && !checkSpeciation(this_node->getSpecRate(),
														   min_spec_rate, this_node->getGenRate()))
		{
			this_node->setSpec(0.0);
		}
	}
	deleteSpeciesList();
	createSpeciesList();
	writeSpeciesList(nodes->size());
	forceSimCompleteParameter();
	exportDatabase();

}

unsigned long Community::getSpeciesRichness(const unsigned long &community_reference)
{
	if(!bSqlConnection)
	{
		throw FatalException("Attempted to get from sql database without opening database connection.");
	}
	sqlite3_stmt *stmt = nullptr;
	// Now find out the max size of the species_id_list, so we have a count to work from
	string count_command = "SELECT COUNT(DISTINCT(species_id)) FROM SPECIES_ABUNDANCES WHERE no_individuals > 0 ";
	count_command += "AND community_reference == ";
	count_command += to_string(community_reference) + ";";
	sqlite3_prepare_v2(database, count_command.c_str(), static_cast<int>(strlen(count_command.c_str())), &stmt,
					   nullptr);
	sqlite3_step(stmt);
	int tmp_val = sqlite3_column_int(stmt, 0);
	// close the old statement
	sqlite3_finalize(stmt);
	return static_cast<unsigned long>(tmp_val);
}

shared_ptr<map<unsigned long, unsigned long>> Community::getSpeciesAbundances(const unsigned long &community_reference)
{
	if(!bSqlConnection)
	{
		throw FatalException("Attempted to get from sql database without opening database connection.");
	}
	// Check that the table exists
	sqlite3_stmt *stmt1;
	string call1 = "select count(type) from sqlite_master where type='table' and name='SPECIES_ABUNDANCES'";
	int rc = sqlite3_prepare_v2(database, call1.c_str(), static_cast<int>(strlen(call1.c_str())), &stmt1, nullptr);
	if(rc != SQLITE_DONE && rc != SQLITE_OK)
	{
		sqlite3_close(database);
		throw SpeciesException("Could not check for SPECIES_ABUNDANCES table. Error code: " +
							   to_string(rc));
	}
	sqlite3_step(stmt1);
	auto has_species_abundances = static_cast<bool>(sqlite3_column_int(stmt1, 0));
	sqlite3_step(stmt1);
	sqlite3_finalize(stmt1);
	// Get the number of species in the table matching the community reference
	if(!has_species_abundances)
	{
		throw SpeciesException("No SPECIES_ABUNDANCES table has been written yet.");
	}
	string call2 = "select count(distinct(species_id)) from SPECIES_ABUNDANCES where no_individuals>0 and ";
	call2 += "community_reference == " + to_string(community_reference);
	rc = sqlite3_prepare_v2(database, call2.c_str(), static_cast<int>(strlen(call2.c_str())), &stmt1, nullptr);
	if(rc != SQLITE_DONE && rc != SQLITE_OK)
	{
		sqlite3_close(database);
		throw SpeciesException("Could not check for SPECIES_ABUNDANCES table. Error code: " +
							   to_string(rc));
	}
	sqlite3_step(stmt1);
	auto no_species = static_cast<unsigned int>(sqlite3_column_int(stmt1, 0));
	if(no_species == 0)
	{
		stringstream ss;
		ss << "No species found in SPECIES_ABUNDANCES for reference of " << community_reference << endl;
		throw SpeciesException(ss.str());
	}
	sqlite3_step(stmt1);
	sqlite3_finalize(stmt1);
	// Now fetch the species abundances
	string all_commands = "SELECT species_id, no_individuals FROM SPECIES_ABUNDANCES WHERE community_reference ==";
	all_commands += to_string(community_reference) + ";";
	sqlite3_prepare_v2(database, all_commands.c_str(), static_cast<int>(strlen(all_commands.c_str())), &stmt1, nullptr);
	sqlite3_step(stmt1);
	// Copy the data across to the TreeNode data structure.
	// For storing the number of ignored lineages so this can be subtracted off the parent number.
	shared_ptr<map<unsigned long, unsigned long>> output_species_abundances =
			make_shared<map<unsigned long, unsigned long>>();
	unsigned long i = 0;
	while(i < no_species)
	{
		auto species_id = static_cast<unsigned long>(sqlite3_column_int(stmt1, 0));
		auto no_individuals = static_cast<unsigned long>(sqlite3_column_int(stmt1, 1));
		if(no_individuals > 0)
		{
			(*output_species_abundances)[species_id] = no_individuals;
			i++;
		}
		sqlite3_step(stmt1);
	}
	// Now we need to blank all objects
	sqlite3_finalize(stmt1);
	return output_species_abundances;
}

shared_ptr<vector<unsigned long>> Community::getSpeciesAbundances()
{
	return species_abundances;
}


