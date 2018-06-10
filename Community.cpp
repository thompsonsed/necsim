// This file is part of NECSim project which is released under MIT license.
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
#include "Community.h"
#include "DoubleComparison.h"
#include "Filesystem.h"

bool checkSpeciation(const long double &random_number, const long double &speciation_rate,
					 const unsigned long &no_generations)
{
	// bool result1, result2, result3, result4;
	long double res = 1.0 - pow(double(1.0 - speciation_rate), double(no_generations));
	return random_number <= res;
}

CommunityParameters::CommunityParameters(unsigned long reference_in, long double speciation_rate_in, long double time_in,
										 bool fragment_in, unsigned long metacommunity_reference_in,
										 const ProtractedSpeciationParameters &protracted_params)
{
	setup(reference_in, speciation_rate_in, time_in, fragment_in, metacommunity_reference_in, protracted_params);
}

void CommunityParameters::setup(unsigned long reference_in, long double speciation_rate_in, long double time_in, bool fragment_in,
								unsigned long metacommunity_reference_in, const ProtractedSpeciationParameters &protracted_params)
{
	time = time_in;
	speciation_rate = speciation_rate_in;
	fragment = fragment_in;
	reference = reference_in;
	metacommunity_reference = metacommunity_reference_in;
	protracted_parameters = protracted_params;
	updated = false;

}

bool CommunityParameters::compare(long double speciation_rate_in, long double time_in, bool fragment_in,
								  unsigned long metacommunity_reference_in,
								  const ProtractedSpeciationParameters &protracted_params)
{
	if(doubleCompare(double(time_in), double(0.0), 0.00001))
	{
#ifdef DEBUG
		stringstream os;
		os << "Detected generation at t=0.0." << endl;
		writeLog(10, os);
#endif // DEBUG
		return doubleCompare(speciation_rate, speciation_rate_in, speciation_rate * 0.000001) &&
			   fragment == fragment_in && metacommunity_reference == metacommunity_reference_in &&
			   protracted_params == protracted_parameters;
	}
	return doubleCompare(speciation_rate, speciation_rate_in, speciation_rate * 0.000001) &&
		   doubleCompare(time, time_in, time * 0.0001) && fragment == fragment_in &&
		   metacommunity_reference == metacommunity_reference_in &&
			protracted_params == protracted_parameters;
}

bool CommunityParameters::compare(long double speciation_rate_in, long double time_in,
								  unsigned long metacommunity_reference_in,
								  const ProtractedSpeciationParameters &protracted_params)
{
	return doubleCompare(speciation_rate, speciation_rate_in, speciation_rate * 0.000001) &&
		   doubleCompare(time, time_in, 0.0001) && metacommunity_reference == metacommunity_reference_in &&
			protracted_params == protracted_parameters;
}

bool CommunityParameters::compare(unsigned long reference_in)
{
	return reference == reference_in;
}

void CommunitiesArray::pushBack(unsigned long reference, long double speciation_rate, long double time, bool fragment,
								unsigned long metacommunity_reference,
								const ProtractedSpeciationParameters &protracted_params)
{
	CommunityParameters tmp_param(reference, speciation_rate, time, fragment, metacommunity_reference, protracted_params);
	communityParameters.push_back(tmp_param);
}

void CommunitiesArray::pushBack(CommunityParameters tmp_param)
{
	communityParameters.push_back(tmp_param);
}

CommunityParameters &CommunitiesArray::addNew(long double speciation_rate, long double time, bool fragment,
											  unsigned long metacommunity_reference,
											  const ProtractedSpeciationParameters &protracted_params)
{
	unsigned long max_reference = 1;
	for(auto &i : communityParameters)
	{
		if(i.compare(speciation_rate, time, metacommunity_reference, protracted_params))
		{
			if(i.fragment == fragment || !fragment)
			{
				stringstream ss;
				ss << "Non-unique parameter set: " << endl;
				ss << "-speciation rate: " << speciation_rate << endl;
				ss << "-time: " << time << endl;
				ss << "-fragment: " << fragment << endl;
				ss << "-metacommunity reference: " << metacommunity_reference << endl;
				ss << "-protracted speciation min: " << protracted_params.min_speciation_gen << endl;
				ss << "-protracted speciation max: " << protracted_params.max_speciation_gen<< endl;
				writeCritical(ss.str());
				throw FatalException("Tried to get reference for non-unique parameter set in communities. "
									 "Please report this bug.");
			}
			else
			{
				i.fragment = true;
				i.updated = true;
				return i;
			}
		}
		else
		{
			if(i.reference >= max_reference)
			{
				max_reference = i.reference + 1;
			}
		}
	}
	CommunityParameters tmp_param(max_reference, speciation_rate, time, fragment, metacommunity_reference,
								  protracted_params);
	communityParameters.push_back(tmp_param);
	return communityParameters.back();
}

bool CommunitiesArray::hasPair(long double speciation_rate, double time, bool fragment,
							   unsigned long metacommunity_reference,
							   const ProtractedSpeciationParameters &protracted_params)
{
	for(auto &i : communityParameters)
	{
		if(i.compare(speciation_rate, time, fragment, metacommunity_reference, protracted_params))
		{
			return true;
		}
	}
	return false;
}

MetacommunityParameters::MetacommunityParameters(unsigned long reference_in, long double speciation_rate_in,
												 unsigned long metacommunity_size_in)
{
	metacommunity_size = metacommunity_size_in;
	speciation_rate = speciation_rate_in;
	reference = reference_in;
}

bool MetacommunityParameters::compare(long double speciation_rate_in, unsigned long metacommunity_size_in)
{
	return doubleCompare(speciation_rate, speciation_rate_in, speciation_rate * 0.000001) &&
		   metacommunity_size == metacommunity_size_in;
}

bool MetacommunityParameters::compare(unsigned long reference_in)
{
	return reference == reference_in;
}

void MetacommunitiesArray::pushBack(unsigned long reference, long double speciation_rate,
									unsigned long metacommunity_size)
{
	MetacommunityParameters tmp_param(reference, speciation_rate, metacommunity_size);
	calc_array.push_back(tmp_param);
}

void MetacommunitiesArray::pushBack(MetacommunityParameters tmp_param)
{
	calc_array.push_back(tmp_param);
}

unsigned long MetacommunitiesArray::addNew(long double speciation_rate, unsigned long metacommunity_size)
{
	unsigned long max_reference = 1;
	for(auto &i : calc_array)
	{
		if(i.compare(speciation_rate, metacommunity_size))
		{

			throw FatalException("Tried to get reference for non-unique parameter set in metacommunities. "
								 "Please report this bug.");
		}
		else
		{
			if(i.reference >= max_reference)
			{
				max_reference = i.reference + 1;
			}
		}
	}
	MetacommunityParameters tmp_param(max_reference, speciation_rate, metacommunity_size);
	calc_array.push_back(tmp_param);
	return max_reference;
}

bool MetacommunitiesArray::hasPair(long double speciation_rate, unsigned long metacommunity_size)
{
	for(auto &i : calc_array)
	{
		if(i.compare(speciation_rate, metacommunity_size))
		{
			return true;
		}
	}
	return false;
}

bool MetacommunitiesArray::hasPair(unsigned long reference)
{
	for(auto &i : calc_array)
	{
		if(i.compare(reference))
		{
			return true;
		}
	}
	return false;
}

unsigned long MetacommunitiesArray::getReference(long double speciation_rate, unsigned long metacommunity_size)
{
	if(metacommunity_size == 0)
	{
		return 0;
	}
	for(auto &i : calc_array)
	{
		if(i.compare(speciation_rate, metacommunity_size))
		{
			return i.reference;
		}
	}
	return 0;
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
//			os << "W,E,N,S: " << fragment.x_west << ", " << fragment.x_east << ", " << fragment.y_north << ", " << fragment.y_south << endl;
	bIsFragment = true;
}

void Samplematrix::removeFragment()
{
	bIsFragment = false;
}

void Community::setList(Row<TreeNode> *l)
{
	nodes = l;
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
				"ERROR_SPEC_003: Attempt to import samplemask object before simulation parameters: dimensions not known");
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
	//		os << "listsize: " << nodes->size() << endl;
	unsigned long iSpecCount = 0;  // start at 2 because the last species has been burnt already.
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
			if(!(*nodes)[this_node->getParent()].getExistence() && this_node->getExistence() &&
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
	iSpecCount = 0;
	set<unsigned long> species_list;
	// Now loop again, creating a new species for each species that actually exists.
	for(unsigned long i = 1; i < nodes->size(); i++)
	{
		TreeNode *this_node = &(*nodes)[i];
		if(this_node->getExistence() && this_node->hasSpeciated())
		{
			addSpecies(iSpecCount, this_node, species_list);
		}
	}
	// Compress the species IDs so that the we have full mapping of species_ids to integers in range 0:n
	// Only do so if the numbers do not match initially
#ifdef DEBUG
	writeLog(10, "Counting species.");
#endif // DEBUG
	if(!species_list.empty() && iSpecCount != *species_list.rbegin())
	{
		unsigned long tmp_species_count = 0;
		unordered_map<unsigned long, unsigned long> ids_map;
		ids_map.reserve(iSpecCount);
		for(unsigned long i = 1; i < nodes->size(); i++)
		{
			TreeNode *this_node = &(*nodes)[i];
			if(this_node->hasSpeciated() && this_node->getExistence())
			{
				auto map_id = ids_map.find(this_node->getSpeciesID());
				if(map_id == ids_map.end())
				{
					tmp_species_count++;
					this_node->resetSpecies();
					this_node->burnSpecies(tmp_species_count);
					ids_map.emplace(this_node->getSpeciesID(), tmp_species_count);
				}
				else
				{
					this_node->resetSpecies();
					this_node->burnSpecies(map_id->second);
				}
			}
		}
		iSpecCount = tmp_species_count;
	}
	else
	{
		iSpecCount = 0;
		for(unsigned long i = 0; i < nodes->size(); i++)
		{
			TreeNode *this_node = &(*nodes)[i];
			// count all speciation events, not just the ones that exist!
			if(this_node->hasSpeciated() && this_node->getExistence() && this_node->getSpeciesID() != 0)
			{
				iSpecCount++;
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
			if(this_node->getSpeciesID() == 0 && this_node->getExistence())
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
	iSpecies = iSpecCount;
	//		os << "iSpecies: " << iSpecies << endl;
	return iSpecCount;
}

void Community::addSpecies(unsigned long &species_count, TreeNode *tree_node, set<unsigned long> &species_list)
{
	species_count++;
	tree_node->burnSpecies(species_count);
}

void Community::calcSpeciesAbundance()
{
	row_out.setSize(iSpecies + 1);
	//		os << "iSpecies: " << iSpecies << endl;
	for(unsigned long i = 0; i < row_out.size(); i++)
	{
		row_out[i] = 0;
	}
	for(unsigned long i = 1; i < nodes->size(); i++)
	{
		TreeNode *this_node = &(*nodes)[i];
		if(this_node->isTip() &&
		   doubleCompare(this_node->getGeneration(), current_community_parameters->time, 0.0001) &&
		   this_node->getExistence())
		{
#ifdef DEBUG
			if(this_node->getSpeciesID() >= row_out.size())
			{
				throw out_of_range("Node index out of range of abundances size. Please report this bug.");
			}
#endif // DEBUG
			// The line that counts the number of individuals
			row_out[this_node->getSpeciesID()]++;
#ifdef DEBUG
			if(!samplemask.getMaskVal(this_node->getXpos(), this_node->getYpos(),
									  this_node->getXwrap(), this_node->getYwrap()) &&
			   doubleCompare(this_node->getGeneration(), current_community_parameters->time, 0.0001))
			{
				stringstream ss;
				ss << "x,y " << (*nodes)[i].getXpos() << ", " << (*nodes)[i].getYpos() << endl;
				ss << "tip: " << (*nodes)[i].isTip() << " Existance: " << (*nodes)[i].getExistence()
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
				   << p_node->getExistence() << ", " << p_node->getGeneration() << endl;
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
	for(unsigned long i = 0; i < nodes->size(); i++)
	{
		(*nodes)[i].qReset();
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
	if(nodes->size() != 0)
	{
		return;
	}
	writeInfo("Beginning data import...");
	// The sql statement to store the sql statement message object
	sqlite3_stmt *stmt = nullptr;

	// Now find out the max size of the list, so we have a count to work from
	string count_command = "SELECT COUNT(*) FROM SPECIES_LIST;";
	sqlite3_prepare_v2(database, count_command.c_str(), static_cast<int>(strlen(count_command.c_str())), &stmt,
					   nullptr);
	unsigned long datasize;
	// skip first row (should be blank)
	sqlite3_step(stmt);
	datasize = static_cast<unsigned long>(sqlite3_column_int(stmt, 0));
	//		os << "datasize: " << datasize << endl;
	// close the old statement
	sqlite3_finalize(stmt);

	// Create db query
	string all_commands = "SELECT * FROM SPECIES_LIST;";
	sqlite3_prepare_v2(database, all_commands.c_str(), static_cast<int>(strlen(all_commands.c_str())), &stmt, nullptr);
	nodes->setSize(datasize + 1);
	// Check that the file opened correctly.
	sqlite3_step(stmt);
	// Copy the data across to the TreeNode data structure.
	// For storing the number of ignored lineages so this can be subtracted off the parent number.
	unsigned long ignored_lineages = 0;
#ifdef DEBUG
	bool has_printed_error = false;
#endif
	for(unsigned long i = 1; i <= datasize; i++)
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
		//			os << xval << ", " << yval << endl;
		// ignored lineages are now not ignored!
		// TODO fix this properly and check functionality
//		if(tip && !samplemask.getVal(xval, yval) && generationin > generation && false)
//		{
//			ignored_lineages++;
//			sqlite3_step(stmt);
//		}
//		else
//		{
		// the -1 is to ensure that the list includes all lineages, but fills the output from the beginning
		unsigned long index = i - 1 - ignored_lineages;
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
		// Now find out the max size of the list, so we have a count to work from
		string count_command = "SELECT MAX(ID) FROM SPECIES_ABUNDANCES;";
		sqlite3_prepare_v2(database, count_command.c_str(), static_cast<int>(strlen(count_command.c_str())), &stmt,
						   nullptr);
		sqlite3_step(stmt);
		max_species_id = static_cast<unsigned long>(sqlite3_column_int(stmt, 0)) + 1;
		// close the old statement
		sqlite3_finalize(stmt);
	}
}

Row<unsigned long> *Community::getCumulativeAbundances()
{
	unsigned long total = 0;
	for(unsigned long i = 0; i < row_out.size(); i++)
	{
		total += row_out[i];
		row_out[i] = total;
	}
	return &row_out;
}

Row<unsigned long> Community::getRowOut()
{
	return row_out;
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
		// Now find out the max size of the list, so we have a count to work from
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
		// Now find out the max size of the list, so we have a count to work from
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
			ss << "Duplicate insertion of " << current_community_parameters->reference << "into SPECIES_ABUNDANCES.";
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
		for(unsigned long i = 0; i < row_out.size(); i++)
		{
			// only do all the export itself if the value of i is not 0
			// if(row_out[i] != 0)
			//{

			// fixed precision problem - lexical cast allows for printing of very small doubles.
			sqlite3_bind_int(stmt, 1, static_cast<int>(max_species_id++));
			sqlite3_bind_int(stmt, 2, static_cast<int>(i));
			sqlite3_bind_int(stmt, 3, static_cast<int>(row_out[i]));
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
		ss << "parameters already applied, not outputting SPECIES_ABUNDANCES table..." << endl;
		writeInfo(ss.str());
	}
}

bool Community::checkCalculationsPerformed(long double speciation_rate, double time, bool fragments,
										   unsigned long metacommunity_size, long double metacommunity_speciation_rate,
										   ProtractedSpeciationParameters proc_parameters)
{
	auto metacommunity_reference = past_metacommunities.getReference(metacommunity_speciation_rate, metacommunity_size);
	if(metacommunity_reference == 0 && metacommunity_size != 0)
	{
		return false;
	}
	bool has_pair = past_communities.hasPair(speciation_rate, time, fragments, metacommunity_reference, proc_parameters);
#ifdef DEBUG
	stringstream os;
	os << "Checking for calculations with sr=" << speciation_rate << ", t=" << time;
	os << " and ref: " << past_metacommunities.getReference(metacommunity_speciation_rate,
															metacommunity_size) << ": " << has_pair << endl;
	writeLog(10, os);
#endif // DEBUG
	if(fragments && past_communities.hasPair(speciation_rate, time, false, metacommunity_reference, proc_parameters))
	{
		return false;

	}
	if(!fragments && past_communities.hasPair(speciation_rate, time, true, metacommunity_reference, proc_parameters))
	{
		return true;
	}
//	if(past_communities.hasPair(speciation_rate, time, !fragments,
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
	for(unsigned long i = 0; i < row_out.size(); i++)
	{
		if(row_out[i] != 0)
		{
			// fixed precision problem - lexical cast allows for printing of very small doubles.
			sqlite3_bind_int(stmt, 1, static_cast<int>(max_fragment_id++));
			sqlite3_bind_text(stmt, 2, f.name.c_str(), -1, SQLITE_STATIC);
			sqlite3_bind_double(stmt, 3, f.area);
			sqlite3_bind_int(stmt, 4, static_cast<int>(f.num));
			sqlite3_bind_int(stmt, 5, static_cast<int>(i));
			sqlite3_bind_int(stmt, 6, static_cast<int>(row_out[i]));
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
				ss << "ERROR_SQL_016: Connection to output database cannot be opened. Check write access "
					  "on output folder. Error code: "
				   << rc << "." << endl;
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
		if(rc != SQLITE_OK && rc != SQLITE_DONE)
		{
			ss << "ERROR_SQL_016: Connection to output database cannot be opened. Check write access on "
				  "output folder. Error code: "
			   << rc << "." << endl;
			throw FatalException(ss.str());
		}
		rc = sqlite3_backup_finish(backupdb);
		if(rc != SQLITE_OK && rc != SQLITE_DONE)
		{
			ss << "ERROR_SQL_016: Connection to output database cannot be opened. Check write access on "
				  "output folder. Error code: "
			   << rc << "." << endl;
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
	// Now find out the max size of the list, so we have a count to work from
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
	// Now find out the max size of the list, so we have a count to work from
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
		//			os << nodes[i].getExistence() << endl;
		if(this_node->isTip() &&
		   this_node->getExistence() && doubleCompare(static_cast<double>(this_node->getGeneration()),
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
				break;
			}
			// Fragment name and dimensions
			try
			{
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

void Community::setSimParameters(const SimParameters *sim_parameters)
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
		os << "Reading parameters..." << flush;
#endif
		sqlite3_stmt *stmt2;
		string sql_parameters = "SELECT speciation_rate, grid_x, grid_y, protracted, min_speciation_gen, max_speciation_gen, "
								"sample_x_offset, sample_y_offset, sample_x, sample_y  FROM SIMULATION_PARAMETERS;";
		int rc = sqlite3_prepare_v2(database, sql_parameters.c_str(), static_cast<int>(strlen(sql_parameters.c_str())),
									&stmt2, nullptr);
		if(rc != SQLITE_DONE && rc != SQLITE_OK)
		{
			stringstream ss;
			ss << "ERROR_SQL_020: FATAL. Could not open simulation parameters in " << file << ". Error code: ";
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
		writeLog(50, "Applied speciation parameters: " + to_string(applied_protracted_parameters.min_speciation_gen) + ", " +
				to_string(applied_protracted_parameters.max_speciation_gen));
		writeLog(50, "Simulated speciation parameters: " + to_string(min_speciation_gen) + ", " +
					 to_string(max_speciation_gen));
#endif // DEBUG
		stringstream ss;
		ss << "Applied protracted speciation parameters: " << applied_protracted_parameters.min_speciation_gen << ", ";
		ss << applied_protracted_parameters.max_speciation_gen << endl;
		ss << "Original protracted speciation parameters: " << min_speciation_gen << ", " << max_speciation_gen << endl;
		writeCritical(ss.str());
		throw SpeciesException("Cannot use protracted parameters with minimum > simulated minimum or "
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
	// Read the community parameters and store them in the relevant objects
	sqlite3_stmt *stmt1;
	string call1 = "select count(type) from sqlite_master where type='table' and name='COMMUNITY_PARAMETERS'";
	int rc = sqlite3_prepare_v2(database, call1.c_str(), static_cast<int>(strlen(call1.c_str())), &stmt1, nullptr);
	if(rc != SQLITE_DONE && rc != SQLITE_OK)
	{
		sqlite3_close(database);
		throw SpeciesException("ERROR_SQL_020: FATAL. Could not check for COMMUNITY_PARAMETERS table. Error code: " +
							   to_string(rc));
		//				exit(EXIT_FAILURE);
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
			sqlite3_close(database);
			throw SpeciesException("ERROR_SQL_020: FATAL. Could not detect COMMUNITY_PARAMETERS table. Error code: " +
								   to_string(rc));
		}
		rc = sqlite3_step(stmt2);
		while(rc == SQLITE_ROW)
		{
			auto row_val = sqlite3_column_int(stmt2, 0);
			if(row_val == 0)
			{
				writeWarning("Reference of 0 found in community parameters in database, skipping...\n");
			}
			else
			{
				ProtractedSpeciationParameters tmp{};
				if(protracted)
				{
					tmp.min_speciation_gen = sqlite3_column_double(stmt2, 5);
					tmp.max_speciation_gen= sqlite3_column_double(stmt2, 6);
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
			ss << "ERROR_SQL_020b: FATAL. Could not read community parameters." << endl;
			ss << "Code: " << rc << endl << "Errmsg: ";
			ss << sqlite3_errmsg(database) << endl;
			sqlite3_clear_bindings(stmt2);
			sqlite3_reset(stmt2);
			throw SpeciesException(ss.str());
		}
		sqlite3_finalize(stmt2);
	}
	// And the same for metacommunity parameters
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
		string call4 = "SELECT reference, speciation_rate, metacommunity_size FROM ";
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
										  sqlite3_column_double(stmt4, 1),
										  static_cast<unsigned long>(sqlite3_column_int(stmt4, 2)));
			rc = sqlite3_step(stmt4);
		}
		if(rc != SQLITE_OK && rc != SQLITE_DONE)
		{
			stringstream ss;
			ss << "ERROR_SQL_020: FATAL. Could not read metacommunity parameters." << endl;
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

void Community::addCalculationPerformed(long double speciation_rate, double time, bool fragments,
										unsigned long metacommunity_size, long double metacommunity_speciation_rate,
										const ProtractedSpeciationParameters &protracted_params)
{
	auto meta_reference = past_metacommunities.getReference(metacommunity_speciation_rate,
															metacommunity_size);
	if(meta_reference == 0 && metacommunity_size != 0)
	{
#ifdef DEBUG
		stringstream ss;
		ss << "Adding metacommunity (" << metacommunity_size << ", " << metacommunity_speciation_rate << ")" << endl;
		writeInfo(ss.str());
#endif
		meta_reference = past_metacommunities.addNew(metacommunity_speciation_rate, metacommunity_size);
	}
	current_community_parameters = &past_communities.addNew(speciation_rate, time, fragments, meta_reference,
															protracted_params);
#ifdef DEBUG
	for(auto &i : past_communities.communityParameters)
	{
		if(doubleCompare(i.time, current_community_parameters->time, 0.00001) &&
			doubleCompare(i.speciation_rate, current_community_parameters->speciation_rate,
						  i.speciation_rate*0.00001) &&
				i.protracted_parameters == current_community_parameters->protracted_parameters &&
				i.metacommunity_reference == current_community_parameters->metacommunity_reference &&
				i.reference != current_community_parameters->reference)
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
		sqlite3_close(database);
		throw SpeciesException("ERROR_SQL_020: FATAL. Could not check for COMMUNITY_PARAMETERS table. Error code: " +
							   to_string(rc));
		//				exit(EXIT_FAILURE);
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
			sqlite3_close(database);
			throw SpeciesException("ERROR_SQL_020: FATAL. Could not detect COMMUNITY_PARAMETERS table. Error code: " +
								   to_string(rc));
		}
		rc = sqlite3_step(stmt2);
		while(rc != SQLITE_DONE)
		{
			unique_community_refs.push_back(static_cast<unsigned long>(sqlite3_column_int(stmt2, 0)));
			rc = sqlite3_step(stmt2);
			if(rc > 10000)
			{
				throw SpeciesException("ERROR_SQL_020: FATAL. Could not read speciation rates.");
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
	// Read the community parameters and store them in the relevant objects
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
	// Find new community parameters to add
	auto unique_community_refs = getUniqueCommunityRefs();
	CommunitiesArray communities_to_write;
	for(auto &community_param : past_communities.communityParameters)
	{
		if(find(unique_community_refs.begin(),
				unique_community_refs.end(), community_param.reference) == unique_community_refs.end())
		{
			communities_to_write.pushBack(community_param);
			unique_community_refs.push_back(community_param.reference);
		}
	}
	if(!communities_to_write.communityParameters.empty())
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
		for(auto &item : communities_to_write.communityParameters)
		{
			if(item.reference == 0)
			{
				continue;
			}
			sqlite3_bind_int(stmt, 1, static_cast<int>(item.reference));
			sqlite3_bind_double(stmt, 2, static_cast<double>(item.speciation_rate));
			sqlite3_bind_double(stmt, 3, static_cast<double>(item.time));
			sqlite3_bind_int(stmt, 4, static_cast<int>(item.fragment));
			sqlite3_bind_int(stmt, 5, static_cast<int>(item.metacommunity_reference));
			if(protracted)
			{
				sqlite3_bind_int(stmt, 6, static_cast<int>(item.protracted_parameters.min_speciation_gen));
				sqlite3_bind_int(stmt, 7, static_cast<int>(item.protracted_parameters.max_speciation_gen));
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

void Community::writeNewMetacommuntyParameters()
{
	auto unique_metacommunity_refs = getUniqueMetacommunityRefs();
	MetacommunitiesArray metacommunities_to_write;
	if(unique_metacommunity_refs.empty())
	{
		for(auto &community_param : past_metacommunities.calc_array)
		{
			metacommunities_to_write.pushBack(community_param);
		}
	}
	else
	{
		for(auto &community_param : past_metacommunities.calc_array)
		{
			if(find(unique_metacommunity_refs.begin(),
					unique_metacommunity_refs.end(), community_param.reference) == unique_metacommunity_refs.end())
			{
				metacommunities_to_write.pushBack(community_param);
				unique_metacommunity_refs.push_back(community_param.reference);
			}
		}
	}
	if(!metacommunities_to_write.calc_array.empty())
	{
		// Create the table if it doesn't exist
		string table_command = "CREATE TABLE IF NOT EXISTS METACOMMUNITY_PARAMETERS (reference INT PRIMARY KEY NOT NULL,"
							   " speciation_rate DOUBLE NOT NULL, metacommunity_size DOUBLE NOT NULL);";
		sqlite3_exec(database, table_command.c_str(), nullptr, nullptr, nullptr);
		sqlite3_stmt *stmt = nullptr;
		table_command = "INSERT INTO METACOMMUNITY_PARAMETERS (reference, speciation_rate, metacommunity_size"
						") VALUES (?,?,?);";
		sqlite3_prepare_v2(database, table_command.c_str(), static_cast<int>(strlen(table_command.c_str())), &stmt,
						   nullptr);
		// Then add the required elements
		sqlite3_exec(database, "BEGIN TRANSACTION;", nullptr, nullptr, nullptr);
		for(auto &item : metacommunities_to_write.calc_array)
		{
			if(item.reference == 0)
			{
				continue;
			}
			sqlite3_bind_int(stmt, 1, static_cast<int>(item.reference));
			sqlite3_bind_double(stmt, 2, static_cast<double>(item.speciation_rate));
			sqlite3_bind_int(stmt, 3, static_cast<int>(item.metacommunity_size));
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
				ss << "Metacommunity reference: " << item.reference << endl;
				ss << "Speciation rate: " << item.speciation_rate << ", metacommunity size: " << item.metacommunity_size << endl;
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

void Community::updateCommunityParameters()
{
	for(auto parameter : past_communities.communityParameters)
	{
		if(parameter.updated)
		{
			if(!bSqlConnection)
			{
				throw FatalException("Attempted to update sql database without opening database connection.");
			}

			// Now find out the max size of the list, so we have a count to work from
			string count_command = "UPDATE COMMUNITY_PARAMETERS SET fragments = 1 WHERE reference = ";
			count_command += to_string(parameter.reference) + ";";
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
	sort(spec_sim_parameters->all_speciation_rates.begin(), spec_sim_parameters->all_speciation_rates.end());
	if(!spec_sim_parameters->bMultiRun)
	{
		os << "Speciation rate is " << spec_sim_parameters->all_speciation_rates[0] << endl;
	}
	else
	{
		os << "Speciation rates are: " << flush;
		for(unsigned int i = 0; i < spec_sim_parameters->all_speciation_rates.size(); i++)
		{
			os << spec_sim_parameters->all_speciation_rates[i] << flush;
			if(i + 1 == spec_sim_parameters->all_speciation_rates.size())
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
	if(spec_sim_parameters->metacommunity_size > 0)
	{
		os.str("");
		os << "Metacommunity size: " << spec_sim_parameters->metacommunity_size << endl;
		os << "Metacommunity speciation rate: " << spec_sim_parameters->metacommunity_speciation_rate << endl;
	}
}

void Community::calculateTree()
{
	stringstream os;
	for(auto protracted_params : spec_sim_parameters->protracted_parameters)
	{
		setProtractedParameters(protracted_params);
		for(auto sr : spec_sim_parameters->all_speciation_rates)
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
											   spec_sim_parameters->metacommunity_size,
											   spec_sim_parameters->metacommunity_speciation_rate,
											   applied_protracted_parameters))
				{
					addCalculationPerformed(sr, time, spec_sim_parameters->use_fragments,
											spec_sim_parameters->metacommunity_size,
											spec_sim_parameters->metacommunity_speciation_rate,
											applied_protracted_parameters);
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
					os << "calculation already performed for " << sr << " at time " << time << endl;
					writeInfo(os.str());
				}
			}
		}
	}
}

void Community::output()
{
	writeNewCommunityParameters();
	writeNewMetacommuntyParameters();
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

void Community::apply(SpecSimParameters *sp)
{
	time_t tStart{};
	time_t tEnd{};
	// Start the clock
	time(&tStart);
	applyNoOutput(sp);
	output();
	printEndTimes(tStart, tEnd);
}

void Community::applyNoOutput(SpecSimParameters *sp)
{
	doApplication(sp);
}



void Community::doApplication(SpecSimParameters *sp)
{
	Row<TreeNode> data;
	doApplication(sp, &data);
}

void Community::doApplication(SpecSimParameters *sp, Row<TreeNode> *data)
{
	spec_sim_parameters = sp;
	writeSpeciationRates();
	// Set up the objects
	setList(data);
	importSimParameters(sp->filename);
	importSamplemask(sp->samplemask);
	importData(sp->filename);
	getPreviousCalcs();
	if(sp->use_fragments)
	{
		calcFragments(sp->fragment_config_file);
	}
	calculateTree();
}

void Community::doApplicationInternal(SpecSimParameters *sp, Row<TreeNode> *data)
{
	setInternalDatabase();
	doApplication(sp, data);

}

