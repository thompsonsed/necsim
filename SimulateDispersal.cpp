// This file is part of NECSim project which is released under MIT license.
// See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.

/**
 * @author Samuel Thompson
 * @file SimulateDispersal.cpp
 * @brief Contains the ability to simulate a given dispersal kernel on a specified density map, outputting the effect
 * dispersal distance distribution to an SQL file after n number of dispersal events (specified by the user).
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 */


#include "SimulateDispersal.h"
#include "Logger.h"
#include "CustomExceptions.h"
#include "Filesystem.h"
#include "Community.h"

#include <utility>

void SimulateDispersal::setSequential(bool bSequential)
{
	is_sequential = bSequential;
}

void SimulateDispersal::setSimulationParameters(SimParameters * sim_parameters, bool print)
{
	if(print)
	{
		writeInfo("********************************\n");
		writeInfo("Setting simulation parameters...\n");
	}
	simParameters = sim_parameters;
	if(print)
	{
		simParameters->printSpatialVars();
	}
}

void SimulateDispersal::importMaps()
{
	writeInfo("Starting map import...\n");
	if(simParameters == nullptr)
	{
		throw FatalException("Simulation parameters have not been set.");
	}
	setDispersalParameters();
	density_landscape.setDims(simParameters);
	dispersal_coordinator.setHabitatMap(&density_landscape);
	density_landscape.calcFineMap();
	density_landscape.calcCoarseMap();
	density_landscape.calcOffset();
	density_landscape.calcHistoricalFineMap();
	density_landscape.calcHistoricalCoarseMap();
	density_landscape.setLandscape(simParameters->landscape_type);
	density_landscape.recalculateHabitatMax();
	dataMask.importSampleMask(*simParameters);
}

void SimulateDispersal::setSizes()
{
	checkMaxParameterReference();
	unsigned long i = max_parameter_reference;
	if(num_steps.empty())
	{
		num_steps.insert(1);
	}
	for(const auto &item : num_steps)
	{
		parameter_references.insert(make_pair(item, i));
		i ++;
	}
	distances.clear();
	distances.resize(num_repeats * num_steps.size());

}

void SimulateDispersal::setDispersalParameters()
{
	dispersal_coordinator.setRandomNumber(&random);
	dispersal_coordinator.setGenerationPtr(&generation);
	dispersal_coordinator.setDispersal(simParameters);

}

void SimulateDispersal::setOutputDatabase(string out_database)
{
	// Check the file is a database
	if(out_database.substr(out_database.length() - 3) != ".db")
	{
		throw FatalException("Output database is not a .db file, check file name.");
	}
	// Open our SQL connection to the database
	openSQLiteDatabase(out_database, database);
//	int o2 = sqlite3_open_v2(out_database.c_str(), &database, SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE, "unix-dotfile");
//	if(o2 != SQLITE_OK && o2 != SQLITE_DONE)
//	{
//		throw FatalException("Database file cannot be opened or created.");
//	}
	checkMaxParameterReference();
}

void SimulateDispersal::setNumberRepeats(unsigned long n)
{
	num_repeats = n;
}

void SimulateDispersal::setNumberSteps(const vector<unsigned long> &s)
{
	for(const auto item : s)
	{
		num_steps.insert(item);
	}
}

unsigned long SimulateDispersal::getMaxNumberSteps()
{
	unsigned long max_number_steps = 0;
	if(!num_steps.empty())
	{
		max_number_steps = *num_steps.rbegin();
	}
	else
	{
		throw FatalException("No steps have been set.");
	}
	return max_number_steps;
}

void SimulateDispersal::storeCellList()
{
	unsigned long total = 0;
	unsigned long cell_total = 0;
	// First count the number of density cells and pick a cell size
	for(unsigned long i = 0; i < simParameters->sample_y_size; i++)
	{
		for(unsigned long j = 0; j < simParameters->sample_x_size; j++)
		{
			if(dataMask.getVal(j, i, 0, 0))
			{
				cell_total++;
				total += density_landscape.getVal(j, i, 0, 0, 0.0);
			}
		}
	}
	writeInfo("Choosing from " + to_string(cell_total) + " cells of " + to_string(total) + " individuals.\n");
	cells.resize(total);
	unsigned long ref = 0;
	for(unsigned long i = 0; i < simParameters->sample_y_size; i++)
	{
		for(unsigned long j = 0; j < simParameters->sample_x_size; j++)
		{
			for(unsigned long k = 0; k < density_landscape.getVal(j, i, 0, 0, 0.0); k++)
			{
				cells[ref].x = j;
				cells[ref].y = i;
				ref ++;
			}
		}
	}
}

const Cell& SimulateDispersal::getRandomCell()
{
	auto index = static_cast<unsigned long>(floor(random.d01() * cells.size()));
	return cells[index];
}

void SimulateDispersal::getEndPoint(Cell &this_cell)
{
	Step tmp_step(this_cell);
	dispersal_coordinator.disperse(tmp_step);
	this_cell.x = tmp_step.oldx + tmp_step.oldxwrap * simParameters->sample_x_size;
	this_cell.y  = tmp_step.oldy + tmp_step.oldywrap * simParameters->sample_y_size;
//	return (this->*getValFptr)(dist, angle, this_cell, end_cell);
}

void SimulateDispersal::runMeanDispersalDistance()
{
	writeInfo("Simulating dispersal " + to_string(num_repeats) + " times.\n");
	storeCellList();
	Cell this_cell{};
	this_cell = getRandomCell();
	// Set up the parameter reference
	setSizes();
	for(unsigned long i = 0; i < num_repeats; i++)
	{
		Cell start_cell{};
		if(!is_sequential)
		{
			// This takes into account rejection sampling based on density due to
			// setup process for the cell list
			this_cell = getRandomCell();
		}
		start_cell = this_cell;
		// Check the end point
		getEndPoint(this_cell);
		// Now store the output location
		distances[i] = make_pair(1, distanceBetweenCells(this_cell, start_cell));
	}
	writeInfo("Dispersal simulation complete.\n");
}

void SimulateDispersal::runMeanDistanceTravelled()
{
	stringstream ss;
	ss << "Simulating dispersal " << num_repeats << " times for (";
	// The maximum number of steps
	setSizes();
	unsigned long max_number_steps = getMaxNumberSteps();
	for(const auto &item : num_steps)
	{
		ss << item;
		if(item != max_number_steps)
		{
			ss << ", ";
		}
	}
	ss << ") generations.\n";
	writeInfo(ss.str());
	storeCellList();
	Cell this_cell{}, start_cell{};
	// Reference for the distances vector
	unsigned long dist_i = 0;
	for(unsigned long i = 0; i < num_repeats; i ++)
	{
		writeRepeatInfo(i);
		// iterator for elements in the set.
		auto step_iterator = num_steps.begin();
		this_cell = getRandomCell();
		start_cell = this_cell;
		generation = 0.0;
		// Keep looping until we get a valid end point
		for(unsigned long j = 1; j <= max_number_steps; j ++)
		{
			getEndPoint(this_cell);
			generation += 0.5;
			if(j == *step_iterator)
			{
				distances[dist_i] = make_pair(j, distanceBetweenCells(start_cell, this_cell));
				step_iterator ++;
				dist_i ++;
			}
		}
	}
	writeRepeatInfo(num_repeats);
	writeInfo("\nDispersal simulation complete.\n");
}

void SimulateDispersal::writeRepeatInfo(unsigned long i)
{
	stringstream os;
	os << "\rSimulating dispersal " << i << "/" << num_repeats;
	writeInfo(os.str());
}

void SimulateDispersal::writeDatabase(string table_name)
{
	if(database)
	{
		if(table_name != "DISTANCES_TRAVELLED" && table_name != "DISPERSAL_DISTANCES")
		{
			string message = "Table name " + table_name;
			message += "  is not one of 'DISTANCES_TRAVELLED' or 'DISPERSAL_DISTANCES'.";
			throw FatalException(message);
		}
		// Write out the parameters
		checkMaxParameterReference();
		writeParameters(table_name);
		// Do the sql output
		// First create the table
		char* sErrMsg = nullptr;
		sqlite3_stmt *stmt = nullptr;
		string create_table = "CREATE TABLE IF NOT EXISTS " + table_name + " (id INT PRIMARY KEY not null, ";
		create_table += " distance DOUBLE not null, parameter_reference INT NOT NULL);";
		int rc = sqlite3_exec(database, create_table.c_str(), nullptr, nullptr, &sErrMsg);
		int step;
		if(rc != SQLITE_OK)
		{
			string message = "Could not create " + table_name + " table in database: ";
			throw FatalException(message.append(sErrMsg));
		}
		// Now add the objects to the database
		string insert_table = "INSERT INTO " + table_name + " (id, distance, parameter_reference) VALUES (?, ?, ?);";
		sqlite3_prepare_v2(database, insert_table.c_str(),
						   static_cast<int>(strlen(insert_table.c_str())), &stmt, nullptr);
		// Start the transaction
		rc = sqlite3_exec(database, "BEGIN TRANSACTION;", nullptr, nullptr, nullptr);
		if(rc != SQLITE_OK)
		{
			throw FatalException("Cannot start SQL transaction.");
		}
		unsigned long max_id = checkMaxIdNumber(table_name);
		for(unsigned long i = 0; i < distances.size(); i++)
		{
			sqlite3_bind_int(stmt, 1, static_cast<int>(max_id + i));
			auto iter = parameter_references.find(distances[i].first);

#ifdef DEBUG
			if(iter == parameter_references.end())
			{
				throw FatalException("Cannot find parameter reference. Please report this bug.");
			}
#endif // DEBUG
			unsigned long reference = iter->second;
			if(reference > max_parameter_reference)
			{
				max_parameter_reference = reference;
			}
			sqlite3_bind_double(stmt, 2, distances[i].second);
			sqlite3_bind_int(stmt, 3, static_cast<int>(reference));
			step = sqlite3_step(stmt);
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
				ss << "SQLITE error code: " << step << endl;
				ss << sqlite3_errmsg(database) << endl;
				ss << "Could not insert into database." << endl;
				throw  FatalException(ss.str());
			}
			sqlite3_clear_bindings(stmt);
			sqlite3_reset(stmt);
		}
		rc = sqlite3_exec(database, "END TRANSACTION;", nullptr, nullptr, &sErrMsg);
		if(rc != SQLITE_OK)
		{
			string message = "Cannot end the SQL transaction: ";
			throw FatalException(message.append(sErrMsg));
		}
		// Need to finalise the statement
		rc = sqlite3_finalize(stmt);
		if(rc != SQLITE_OK)
		{
			string message = "Cannot finalise the SQL transaction: ";
			throw FatalException(message.append(sErrMsg));
		}

	}
	else
	{
		throw FatalException("Database connection has not been opened, check programming.");
	}
	clearParameters();
}

void SimulateDispersal::writeParameters(string table_name)
{
	// Now add the parameters
	string create_table = "CREATE TABLE IF NOT EXISTS PARAMETERS (ref INT PRIMARY KEY not null,";
	create_table += "simulation_type TEXT not null, ";
	create_table += " sigma DOUBLE not null, tau DOUBLE not null, m_prob DOUBLE not null, cutoff DOUBLE NOT NULL,";
	create_table += "dispersal_method TEXT not null, map_file TEXT not null, seed INT NOT NULL, number_steps ";
	create_table += "INT NOT NULL, number_repeats INT NOT NULL);";
	char * sErrMsg;
	int rc = sqlite3_exec(database, create_table.c_str(), nullptr, nullptr, &sErrMsg);
	if(rc != SQLITE_OK)
	{
		string message = "Could not create PARAMETERS table in database: ";
		throw FatalException(message.append(sErrMsg));
	}
	for(const auto &item : parameter_references)
	{
		string insert_table =
				"INSERT INTO PARAMETERS VALUES(" + to_string(item.second) + ", '" + table_name + "',";
		insert_table += to_string((long double) simParameters->sigma) + ",";
		insert_table +=
				to_string((long double) simParameters->tau) + ", " + to_string((long double) simParameters->m_prob);
		insert_table +=
				", " + to_string((long double) simParameters->cutoff) + ", '" + simParameters->dispersal_method + "','";
		insert_table += simParameters->fine_map_file + "', " + to_string(seed) + ", " + to_string(item.first) + ", ";
		insert_table += to_string(num_repeats) + ");";
		rc = sqlite3_exec(database, insert_table.c_str(), nullptr, nullptr, &sErrMsg);
		if(rc != SQLITE_OK)
		{
			string message = "Could not insert into PARAMETERS table in database. \n";
			message += "Error: ";
			throw FatalException(message.append(sErrMsg));
		}
	}
}

void SimulateDispersal::clearParameters()
{
	distances.clear();
	parameter_references.clear();
	num_steps.clear();
}

void SimulateDispersal::checkMaxParameterReference()
{
	string to_exec = "SELECT CASE WHEN COUNT(1) > 0 THEN MAX(ref) ELSE 0 END AS [Value] FROM PARAMETERS;";
	sqlite3_stmt *stmt = nullptr;
	sqlite3_prepare_v2(database, to_exec.c_str(), static_cast<int>(strlen(to_exec.c_str())), &stmt, nullptr);
	int rc = sqlite3_step(stmt);
	max_parameter_reference = static_cast<unsigned long>(sqlite3_column_int(stmt, 0) + 1);
	// close the old statement
	rc = sqlite3_finalize(stmt);
	if(rc != SQLITE_OK && rc != SQLITE_DONE)
	{
		stringstream ss;
		ss << "Could not check max parameter reference. Error code: " << rc << "\n";
		throw SpeciesException(ss.str());
	}
}

unsigned long SimulateDispersal::checkMaxIdNumber(string table_name)
{
	string to_exec = "SELECT CASE WHEN COUNT(1) > 0 THEN MAX(id) ELSE 0 END AS [Value] FROM " + table_name +";";
	sqlite3_stmt *stmt = nullptr;
	sqlite3_prepare_v2(database, to_exec.c_str(), static_cast<int>(strlen(to_exec.c_str())), &stmt, nullptr);
	int rc = sqlite3_step(stmt);
	auto max_id = static_cast<unsigned long>(sqlite3_column_int(stmt, 0) + 1);
	// close the old statement
	rc = sqlite3_finalize(stmt);
	if(rc != SQLITE_OK && rc != SQLITE_DONE)
	{
		stringstream ss;
		ss << "Could not check max id number. Error code: " << rc << "\n";
		throw SpeciesException(ss.str());
	}
	return max_id;
}


