/**
 * @author Sam Thompson
 * @file LogFile.h
 * @brief Contains a class for logging to a logfile, including reporting level and timestamps.
 * @copyright <a href="https://opensource.org/licenses/BSD-3-Clause">BSD-3 Licence.</a>
 */

#ifndef LOGFILE_H
#define LOGFILE_H

#include <cstring>
#include <fstream>
#include <ctime>
#include <map>
#include "CustomExceptions.h"
#define LOGNAME_FORMAT "%d%m%Y_%H%M%S"

using namespace std;
/**
 * @brief Gets the default log file path.
 * Stored at log/DDMMYYYY_HHMMSS.log where DDMMYYYY_HHMMSS are replaced by the current date and time.
 * @return location of the default logging position.
 */
string getDefaultLogFile();

/**
 * @brief Modifies the file name so that it doesn't point to an existing file.
 * @param basic_string the file name to modify
 */
void getUniqueFileName(string &basic_string);

/**
 * @brief Contains routines for writing to log files.
 */
class LogFile
{
protected:
	// output stream to log file
	ofstream output_stream;
	// log file name
	string file_name;
	// mapping integer levels to logging level
	map<int, string> levels_map;
	// Makes the class non-copyable as we don't want to copy file streams
	LogFile(const LogFile&) = delete;
	LogFile& operator=(const LogFile&) = delete;
public:
	/**
	 * @brief Default constructor.
	 */
	LogFile();

	/**
	 * @brief Constructor taking location of a log file.
	 * @param file_name_in the path to the log file to open
	 */
	explicit LogFile(const string &file_name_in);
	/**
	 * @brief Default destructor, including writing closure out to log file.
	 */
	~LogFile();

	void init(const string &file_name_in);

	/**
	 * @brief Writes the message out to the logfile at the specified logging level.
	 * @param level the level of logging severity
	 * @param message the message to write out
	 */
	void write(const int &level, string message);

	/**
	 * @brief Writes the message out to the logfile at the specified logging level.
	 * @param level the level of logging severity
	 * @param message the message to write out
	 */
	void write(const int &level, stringstream &message);
};


#endif //LOGFILE_H
