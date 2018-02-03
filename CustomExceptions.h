//This file is part of NECSim project which is released under BSD-3 license.
//See file **LICENSE.txt** or visit https://opensource.org/licenses/BSD-3-Clause) for full license details.
// Author: Samuel Thompson
// Contact: samuel.thompson14@imperial.ac.uk or thompsonsed@gmail.com
/**
 * @author Samuel Thompson
 * @file CustomExceptions.h
 * @brief Contains the various exceptions used by NECSim.
 * 
 * @copyright <a href="https://opensource.org/licenses/BSD-3-Clause">BSD-3 Licence.</a>
 */

#ifndef EXCEPTION
#define EXCEPTION

#include <stdexcept>
#include "Logging.h"

using namespace std;


/*!
 * @struct FatalException
 * @brief  This is called any time a fatal exception is called and the program is unwound and ended.
 */
struct FatalException : public runtime_error
{
	FatalException() : runtime_error("Fatal exception thrown at run time, quitting program. "){}

	explicit FatalException(string msg) : runtime_error(msg)
	{
#ifdef DEBUG
		writeLog(50, msg);
#endif //DEBUG
	}
};

/**
 * @struct ConfigException
 * @brief A structure for all exceptions thrown within config processes.
 */
struct ConfigException : public FatalException
{
	ConfigException() : FatalException("Exception thrown at run time in config: "){}

	explicit ConfigException(string msg) : FatalException(msg){}
};


/**
 * @struct SpeciesException
 * @brief An exception thrown whenever a non-fatal Species exception is thrown.
 */

struct SpeciesException : public FatalException
{
	/**
	 * @brief Throws a runtime_error with a custom message indicating source.
	 */
	SpeciesException() : FatalException("Exception thrown at run time in SpeciationCounter: "){}

	/**
	 * @brief Overloaded runtime_error call which provides error message parsing
	 * @param msg the message to be passed to the runtime_error
	 */
	explicit SpeciesException(string msg) : FatalException(msg){}
};

#endif