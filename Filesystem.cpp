// This file is part of NECSim project which is released under BSD-3 license.
// See file **LICENSE.txt** or visit https://opensource.org/licenses/BSD-3-Clause) for full license details.

/**
 * @author Samuel Thompson
 * @date 19/07/2017
 * @file Filesystem.cpp
 *
 * @copyright <a href="https://opensource.org/licenses/BSD-3-Clause">BSD-3 Licence.</a>
 * @brief Contains routines for checking files and folder exist, opening sqlite databases safely, with support for various
 * virtual filesystems, and checking parents of a file exist.
 *
 * Contact: samuel.thompson14@imperial.ac.uk or thompsonsed@gmail.com
 */
#include <string>
#include <sstream>
#include <zconf.h>
#include <boost/filesystem.hpp>
#include "Filesystem.h"
#include "CustomExceptions.h"
#include "Logging.h"

void openSQLiteDatabase(const string &database_name, sqlite3 * &database)
{
	int rc = sqlite3_open_v2(database_name.c_str(), &database, SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE, "unix-dotfile");
	if(rc != SQLITE_OK && rc != SQLITE_DONE)
	{
		int i = 0;
		while((rc != SQLITE_OK && rc != SQLITE_DONE) && i < 10)
		{
			i++;
			sleep(1);
			rc = sqlite3_open_v2(database_name.c_str(), &database, SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE,
								 "unix-dotfile");
		}
		// Attempt different opening method if the first fails.
		int j = 0;
		while((rc != SQLITE_OK && rc != SQLITE_DONE) && j < 10)
		{
			j++;
			sleep(1);
			rc = sqlite3_open(database_name.c_str(), &database);
		}
		if(rc != SQLITE_OK && rc != SQLITE_DONE)
		{
			stringstream ss;
			ss << "ERROR_SQL_010: SQLite database file could not be opened. Check the folder exists and you "
					"have write permissions. (REF1) Error code: "
				 << rc << endl;
			ss << " Attempted call " << max(i, j) << " times" << endl;
			throw FatalException(ss.str());
		}
	}
}

void createParent(const string &file)
{
	boost::filesystem::path file_path(file);
	if(!boost::filesystem::exists(file_path.parent_path()))
	{
		if(!boost::filesystem::create_directories(file_path.parent_path()))
		{
			throw FatalException("Cannot create parent folder for " + file);
		}
	}
}



bool doesExist(string testfile)
{
	if(boost::filesystem::exists(testfile))
	{
		stringstream os;
		os << "\rChecking folder existance..." << testfile << " exists!               " << endl;
		writeInfo(os.str());
		return true;
	}
	else
	{
		throw runtime_error(string("ERROR_MAIN_008: FATAL. Input or output folder does not exist: " + testfile + "."));
	}
}

bool doesExistNull(string testfile)
{
	return testfile == "null" || testfile == "none" || doesExist(testfile);
}

unsigned long cantorPairing(unsigned long x1, unsigned long x2)
{
	return ((x1 + x2) * (x1 + x2 + 1)/2) + x2;
}