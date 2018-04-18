// This file is part of NECSim project which is released under MIT license.
// See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.
/**
 * @author Sam Thompson
 * @file Logging.h
 * @brief Routines for writing to cout. Intended to be overloaded for pythonic versions with the logging module.
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 */
#ifndef LOGGING_IMPORT_H
#define LOGGING_IMPORT_H
#include <string>
#include <iostream>
#include <cstdio>
#include <stdexcept>
#include <sstream>
#include "LogFile.h"
#include "CPLCustomHandler.h"

using namespace std;

/**
 * @brief Writes to cout, or to info in logging module if being compiled with python
 * @param message the message to write out
 */
void writeInfo(string message);

/**
 * @brief Writes to cerr, or to warning in logging module if being compiled with python
 * @param message the message to write out
 */
void writeWarning(string message);

/**
 * @brief Writes to cerr, or to error in logging module if being compiled with python
 * @param message the message to write out
 */
void writeError(string message);

/**
 * @brief Writes to cerr, or to critical in logging module if being compiled with python
 * @param message the message to write out
 */
void writeCritical(string message);


#ifdef DEBUG
/**
 * @brief Calls the static logger object for logging out
 * @param level the level of logging severity
 * @param message the message to pass out as a string
 */
void writeLog(const int &level, string message);

/**
 * @brief Calls the static logger object for logging out
 * @param level the level of logging severity
 * @param message the message to pass out as a stringstream
 */
void writeLog(const int &level, stringstream &message);
#endif // DEBUG
#endif // LOGGING_IMPORT_H
