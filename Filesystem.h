// This file is part of NECSim project which is released under MIT license.
// See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.

/**
 * @author Samuel Thompson
 * @date 19/07/2017
 * @file Filesystem.h
 *
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 * @brief Contains routines for checking files and folder exist, opening sqlite databases safely, with support for various
 * virtual filesystems, and checking parents of a file exist.
 *
 * Contact: samuel.thompson14@imperial.ac.uk or thompsonsed@gmail.com
 */

#define _USE_MATH_DEFINES
#include <cmath>
#include <sqlite3.h>
#include <string>

#ifndef SPECIATIONCOUNTER_FILESYSTEM_H
#define SPECIATIONCOUNTER_FILESYSTEM_H

using namespace std;

/**
 * @brief Safely opens a connection to the provided SQLite database.
 *
 * Adds type safety for usage on different filesystems.
 * @param database_name
 * @param database
 */
void openSQLiteDatabase(const string &database_name, sqlite3 *& database);

/**
 * @brief Checks that parent folder to the supplied file exists, and if it doesn't creates it.
 * @param file the file path to check for
 */
void createParent(string file);


/**
* @brief Checks the existance of a file on the hard drive.
* @param testfile the file to examine
* @return if true, file exists
*/
bool doesExist(string testfile);

/**
 * @brief Checks for the existance of a file, but returns true if the file name is 'null'.
 * Note: this function just calls doesExist().
 * @param testfile the file to examine
 * @return if true, file exists (or is null).
 */
bool doesExistNull(string testfile);

/**
 * @brief Generates a unique ID for the pair of provided parameters.
 *
 * 	Maps ZxZ -> N, so only relevant for positive numbers.
	For any A and B, generates C such that no D and E produce C unless D=A and B=E.

 *@deprecated Should not be used for large integers, or of unknown size, as integer overflows are likely. Cantor pairing
 * explodes in size of return value.
 *
 * @param x1 the first integer reference
 * @param x2 the second integer reference
 * @return a unique reference for the two provided integers
 */
unsigned long cantorPairing(unsigned long x1, unsigned long x2);

/**
 * @brief Gets the next line from a csv filestream and splits the row into a vector of strings, where each string is the
 * value from the csv file, delimited by a comma (i.e. each column of the row).
 * @param str the input stream from the csv file.
 * @return a vector where each element corresponds to the respective row from the csv.
 */
vector<string> getCsvLineAndSplitIntoTokens(istream& str);

#endif //SPECIATIONCOUNTER_FILESYSTEM_H
