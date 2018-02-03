/**
 * @author Sam Thompson
 * @file Logging.cpp
 * @brief Routines for writing to cout. Intended to be overloaded for pythonic versions with the logging module.
 * @copyright <a href="https://opensource.org/licenses/BSD-3-Clause">BSD-3 Licence.</a>
 */

#include <sstream>
#include "Logging.h"


using namespace std;

void writeInfo(string message)
{
#ifdef DEBUG
	writeLog(10, message);
#endif // DEBUG
	cout << message << flush;
}

void writeWarning(string message)
{
#ifdef DEBUG
	writeLog(20, message);
#endif // DEBUG
	cerr << message << flush;
}

void writeError(string message)
{
#ifdef DEBUG
	writeLog(30, message);
#endif // DEBUG
	cerr << message << flush;
}

void writeCritical(string message)
{
#ifdef DEBUG
	writeLog(40, message);
#endif // DEBUG
	cerr << message << flush;
}

#ifdef DEBUG
void writeLog(const int &level, string message)
{
	static LogFile logfile;
	logfile.write(level, std::move(message));
}

void writeLog(const int &level, stringstream &message)
{
	writeLog(level, message.str());
}

#endif // DEBUG