// This file is part of NECSim project which is released under BSD-3 license.
// See file **LICENSE.txt** or visit https://opensource.org/licenses/BSD-3-Clause) for full license details.

/**
 * @author Samuel Thompson
 * @file SimulateDispersal.cpp
 * @brief Contains the ability to simulate a given dispersal kernel on a specified density map, outputting the effect
 * dispersal distance distribution to an SQL file after n number of dispersal events (specified by the user).
 * @copyright <a href="https://opensource.org/licenses/BSD-3-Clause">BSD-3 Licence.</a>
 */


#include "SimulateDispersal.h"
#include "CustomExceptions.h"
#include "Filesystem.h"
#include "Community.h"

#include <utility>

double distanceBetween(Cell &c1, Cell &c2)
{
	return pow(pow(c1.x - c2.x, 2) + pow(c1.y - c2.y, 2), 0.5);
}
void SimulateDispersal::setSequential(bool bSequential)
{
	is_sequential = bSequential;
}

void SimulateDispersal::setSizes(unsigned long x, unsigned long y)
{
	if(!has_set_size)
	{
		density_map.SetSize(y, x);
		has_set_size = true;
	}
	else
	{
		throw FatalException("Dimensions of the density map already set.");
	}
}

void SimulateDispersal::importMaps(string map_file)
{
	if(has_set_size)
	{
		map_name = map_file;
		if(map_file != "null")
		{
			density_map.import(map_file);
			// Now loop over the density map to find the maximum value
			for(unsigned long i = 0; i < density_map.GetRows(); i ++)
			{
				for(unsigned long j = 0; j < density_map.GetCols(); j ++)
				{
					if(density_map[i][j] > max_density)
					{
						max_density = density_map[i][j];
					}
				}
			}
			if(max_density < 1)
			{
				throw FatalException("Maximum density on density map is less than 1. Please check your maps.");
			}
		}
		else
		{
			for(unsigned long i = 0; i < density_map.GetRows(); i ++)
			{
				for(unsigned long j = 0; j < density_map.GetCols(); j ++)
				{
					density_map[i][j] = 1;
				}
			}
			max_density = 1;
		}
	}
	else
	{
		throw FatalException("Dimensions of density map not set before importSpatialParameters");
	}
}

void SimulateDispersal::setDispersalParameters(
	string dispersal_method_in, double sigma_in, double tau_in, double m_prob_in, double cutoff_in,
	string landscape_type)
{
	random.setDispersalMethod(dispersal_method_in, m_prob_in, cutoff_in);
	random.setDispersalParams(sigma_in, tau_in);
	setLandscapeType(std::move(landscape_type));
	dispersal_method = dispersal_method_in;
	sigma = sigma_in;
	tau = tau_in;
	m_prob = m_prob_in;
	cutoff = cutoff_in;
}

void SimulateDispersal::setLandscapeType(string landscape_type)
{
	if(landscape_type == "infinite")
	{
		getValFptr = &SimulateDispersal::getEndPointInfinite;
	}
	else if(landscape_type == "closed")
	{
		getValFptr = &SimulateDispersal::getEndPointClosed;
	}
	else if(landscape_type == "tiled")
	{
		getValFptr = &SimulateDispersal::getEndPointTiled;
	}
	else
	{
		throw FatalException("Landscape type not compatible: " + landscape_type);
	}
}

void SimulateDispersal::setOutputDatabase(string out_database)
{
	// Check the file is a database
	if(out_database.substr(out_database.length() - 3) != ".db")
	{
		throw FatalException("Output database is not a .db file, check file name.");
	}
	// Open our SQL connection to the database
	int o2 = sqlite3_open_v2(out_database.c_str(), &database, SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE, "unix-dotfile");
	if(o2 != SQLITE_OK && o2 != SQLITE_DONE)
	{
		throw FatalException("Database file cannot be opened or created.");
	}
}

void SimulateDispersal::setNumberRepeats(unsigned long n)
{
	num_repeats = n;
	distances.resize(num_repeats);
}

void SimulateDispersal::setNumberSteps(unsigned long s)
{
	num_steps = s;
}

void SimulateDispersal::storeCellList()
{
	unsigned long total = 0;
	// First count the number of density cells and pick a cell size
	for(unsigned long i = 0; i < density_map.GetRows(); i++)
	{
		for(unsigned long j = 0; j < density_map.GetCols(); j++)
		{
			total += density_map[i][j];
		}
	}
	cells.resize(total);
	unsigned long ref = 0;
	for(unsigned long i = 0; i < density_map.GetRows(); i++)
	{
		for(unsigned long j = 0; j < density_map.GetCols(); j++)
		{
			for(unsigned long k = 0; k < density_map[i][j]; k++)
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

void SimulateDispersal::calculateNewPosition(const double &dist, const double &angle,
											 const Cell &start_cell, Cell &end_cell)
{
	end_cell.x = (long) floor(start_cell.x + 0.5 + dist * cos(angle));
	end_cell.y = (long) floor(start_cell.y + 0.5 + dist * sin(angle));
}

bool SimulateDispersal::getEndPointInfinite(const double &dist, const double &angle,
											const Cell &this_cell, Cell&end_cell)
{
	if(getEndPointTiled(dist, angle, this_cell, end_cell))
	{
		return true;
	}
	return end_cell.x >= (long) (density_map.GetCols()) || end_cell.x > 0 ||
			end_cell.y >= (long) (density_map.GetRows()) || end_cell.y < 0;
}

bool SimulateDispersal::getEndPointTiled(const double &dist, const double &angle,
										 const Cell &this_cell, Cell &end_cell)
{
	calculateNewPosition(dist, angle, this_cell, end_cell);
	return double(density_map[end_cell.y % density_map.GetCols()][end_cell.x % density_map.GetRows()]) >
			(random.d01() * double(max_density));
}

bool SimulateDispersal::getEndPointClosed(const double &dist, const double &angle,
										  const Cell &this_cell, Cell &end_cell)
{
	calculateNewPosition(dist, angle, this_cell, end_cell);
	return !(end_cell.x >= (long) density_map.GetCols() || end_cell.x > 0 ||
			end_cell.y >= (long) density_map.GetRows() || end_cell.y < 0) &&
		   getEndPointTiled(dist, angle, this_cell, end_cell);
}

bool SimulateDispersal::getEndPoint(const double &dist, const double &angle, const Cell &this_cell, Cell &end_cell)
{
	return (this->*getValFptr)(dist, angle, this_cell, end_cell);
}

void SimulateDispersal::runMeanDispersalDistance()
{
	storeCellList();
	Cell this_cell{};
	this_cell = getRandomCell();
	for(unsigned long i = 0; i < num_repeats; i++)
	{
		if(!is_sequential)
		{
			// This takes into account rejection sampling based on density due to
			// setup process for the cell list
			this_cell = getRandomCell();
		}
		Cell end_cell{};
		bool fail;
		double dist, angle;
		// Keep looping until we get a valid end point
		do
		{
			// Get a random dispersal distance
			dist = random.dispersal();
			angle = random.direction();
			// Check the end point
			fail = !getEndPoint(dist, angle, this_cell, end_cell);
		} while(fail);
		// Copy the end location into this cell
		this_cell = end_cell;
		// Now store the output location
		distances[i] = dist;
	}
}

void SimulateDispersal::runMeanDistanceTravelled()
{
	storeCellList();
	Cell this_cell{}, start_cell{}, end_cell{};
	for(unsigned long i = 0; i < num_repeats; i ++)
	{
		this_cell = getRandomCell();
		start_cell = this_cell;
		bool fail;
		double dist, angle;
		// Keep looping until we get a valid end point
		for(unsigned long j = 0; j < num_steps; j ++)
		{
			do
			{
				dist = random.dispersal();
				angle = random.direction();
				fail = !getEndPoint(dist, angle, this_cell, end_cell);
			}
			while(fail);
			this_cell = end_cell;
		}
		// Now stores the distance travelled
		distances[i] = distanceBetween(start_cell, this_cell);
	}
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
		char* sErrMsg;
		sqlite3_stmt* stmt;
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
			sqlite3_bind_double(stmt, 2, distances[i]);
			sqlite3_bind_int(stmt, 3, static_cast<int>(parameter_reference));
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
	string insert_table = "INSERT INTO PARAMETERS VALUES(" + to_string(parameter_reference) + ", '" + table_name + "',";
	insert_table += to_string((long double)sigma) + ",";
	insert_table += to_string((long double)tau) + ", " +  to_string((long double)m_prob);
	insert_table += ", " + to_string((long double)cutoff) + ", '" + dispersal_method + "','";
	insert_table += map_name + "', " + to_string(seed) + ", " + to_string(num_steps) + ", ";
	insert_table += to_string(num_repeats) + ");";
	rc = sqlite3_exec(database, insert_table.c_str(), nullptr, nullptr, &sErrMsg);
	if(rc != SQLITE_OK)
	{
		string message = "Could not insert into PARAMETERS table in database. \n";
		message += "Error: ";
		throw FatalException(message.append(sErrMsg));
	}
}

void SimulateDispersal::checkMaxParameterReference()
{
	string to_exec = "SELECT CASE WHEN COUNT(1) > 0 THEN MAX(ref) ELSE 0 END AS [Value] FROM PARAMETERS;";
	sqlite3_stmt *stmt;
	sqlite3_prepare_v2(database, to_exec.c_str(), static_cast<int>(strlen(to_exec.c_str())), &stmt, nullptr);
	int rc = sqlite3_step(stmt);
	parameter_reference = static_cast<unsigned long>(sqlite3_column_int(stmt, 0) + 1);
	// close the old statement
	rc = sqlite3_finalize(stmt);
	if(rc != SQLITE_OK && rc != SQLITE_DONE)
	{
		cerr << "rc: " << rc << endl;
		throw SpeciesException("Could not detect dimensions");
	}
}

unsigned long SimulateDispersal::checkMaxIdNumber(string table_name)
{
	string to_exec = "SELECT CASE WHEN COUNT(1) > 0 THEN MAX(id) ELSE 0 END AS [Value] FROM " + table_name +";";
	sqlite3_stmt *stmt;
	sqlite3_prepare_v2(database, to_exec.c_str(), static_cast<int>(strlen(to_exec.c_str())), &stmt, nullptr);
	int rc = sqlite3_step(stmt);
	unsigned long max_id = static_cast<unsigned long>(sqlite3_column_int(stmt, 0) + 1);
	// close the old statement
	rc = sqlite3_finalize(stmt);
	if(rc != SQLITE_OK && rc != SQLITE_DONE)
	{
		cerr << "rc: " << rc << endl;
		throw SpeciesException("Could not detect dimensions");
	}
	return max_id;
}


