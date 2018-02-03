//This file is part of NECSim project which is released under BSD-3 license.
//See file **LICENSE.txt** or visit https://opensource.org/licenses/BSD-3-Clause) for full license details.
// 
/**
 * @author Sam Thompson
 * @file Setup.cpp
 * @brief Contains the command line parsing and setup options for NECSim.
 * 
 * Contact: samuel.thompson14@imperial.ac.uk or thompsonsed@gmail.com
 * @copyright <a href="https://opensource.org/licenses/BSD-3-Clause">BSD-3 Licence.</a>
 */
 
#include "Setup.h"
#include "Logging.h"
// Global variables
// store the log file name for access anywhere.
string log_name = "null"; 
// the old stdout 
int saved_stdout;


#ifndef verbose

void openLogFile(bool append)
{
	// If verbose mode is not defined, we output to a log file instead of to the terminal
	// Get the current time
	saved_stdout = dup(fileno(stdout));
	//dup2(stdout, 1);
	auto t = time(nullptr);
	auto tp = *localtime(&t);
	// get the time string
	ostringstream oss;
	oss << std::put_time(&tp, "%d-%m-%Y-%H:%M:%S");
	if(log_name == "null")
	{
		log_name = "Logs/Log_"+oss.str() + ".txt";
	}
	// Check that the Log folder exists, and create if necessary.
	if(!boost::filesystem::exists("Logs"))
	{
		if(!boost::filesystem::create_directory("Logs"))
		{
			cerr << "Cannot create log directory (check write access) - defaulting to terminal." << endl;
		}
	}
	if(boost::filesystem::exists("Logs"))
	{
		// Open the log file for writing to.
		FILE * tmpfileptr;
		if(append)
		{
			tmpfileptr = freopen(log_name.c_str(), "a", stdout);
		}
		else
		{
			tmpfileptr = freopen(log_name.c_str(), "w", stdout);
		}
		if(stdout == nullptr || tmpfileptr == nullptr)
		{
			cerr << "Cannot create log file (check write access) - defaulting to terminal." << endl;
			dup2(saved_stdout, fileno(stdout));
			close(saved_stdout);
		}
	}
	//cerr<< "logfile2 : " << log_name << endl;
}
#endif /* verbose */

#ifdef PROFILE

ofstream csv_output;

void logToCsv(string place, time_t start, time_t end)
{
	if(!csv_output.good())
	{
		try
		{
			csv_output.open("csvout.csv");
		}
		catch(exception &e)
		{
			throw Fatal_Main_Exception("Csv logging output not good: " + e.what());
		}
	}
	csv_output << place << "," << start << "," << end << endl;
}

void closeCsv()
{
	if(csv_output.good())
	{
		csv_output.close();
	}
}
#endif

void runAsDefault(vector<string> &comargs)
{
	writeInfo("Setting default variables on small grid.\n");
	comargs.push_back("-f");
	comargs.push_back("1");
	comargs.push_back("10");
	comargs.push_back("10");
	comargs.push_back("null");
	comargs.push_back("150");
	comargs.push_back("150");
	comargs.push_back("25");
	comargs.push_back("25");
	comargs.push_back("null");
	comargs.push_back("2000");
	comargs.push_back("2000");
	comargs.push_back("500");
	comargs.push_back("500");
	comargs.push_back("100");
	comargs.push_back("Default/");
	comargs.push_back("0.000009");
	comargs.push_back("2");
	comargs.push_back("1");
	comargs.push_back("1");
	comargs.push_back("4");
	comargs.push_back("1");
	comargs.push_back("0");
	comargs.push_back("100");
	comargs.push_back("null");
	comargs.push_back("null");
	comargs.push_back("0.5");
	comargs.push_back("20.0");
	comargs.push_back("2.0");
	comargs.push_back("null");
	comargs.push_back("null");
	comargs.push_back("0.000009");
}

void runLarge(vector<string> &comargs)
{
	writeInfo("Setting default variables on large grid.\n");
	comargs.push_back("-f");
	comargs.push_back("1");
	comargs.push_back("500");
	comargs.push_back("500");
	comargs.push_back("null");
	comargs.push_back("500");
	comargs.push_back("500");
	comargs.push_back("0");
	comargs.push_back("0");
	comargs.push_back("null");
	comargs.push_back("100");
	comargs.push_back("100");
	comargs.push_back("2500");
	comargs.push_back("2500");
	comargs.push_back("100");
	comargs.push_back("Default/");
	comargs.push_back("0.0001");
	comargs.push_back("2");
	comargs.push_back("10");
	comargs.push_back("1");
	comargs.push_back("3600");
	comargs.push_back("1");
	comargs.push_back("1");
	comargs.push_back("50000");
	comargs.push_back("null");
	comargs.push_back("null");
	comargs.push_back("0.5");
	comargs.push_back("20.0");
	comargs.push_back("2.0");
	comargs.push_back("null");
	comargs.push_back("null");
	comargs.push_back("0.001");
}

void runXL(vector<string> &comargs)
{
	writeInfo("Setting default variables on extra large grid.\n");
	comargs.push_back("-f");
	comargs.push_back("1");
	comargs.push_back("6000");
	comargs.push_back("6400");
	comargs.push_back("null");
	comargs.push_back("34000");
	comargs.push_back("28000");
	comargs.push_back("8800");
	comargs.push_back("14800");
	comargs.push_back("null");
	comargs.push_back("24000");
	comargs.push_back("20000");
	comargs.push_back("10320");
	comargs.push_back("8080");
	comargs.push_back("10");
	comargs.push_back("Default/");
	comargs.push_back("0.0000001");
	comargs.push_back("2");
	comargs.push_back("49");
	comargs.push_back("0.2");
	comargs.push_back("21600");
	comargs.push_back("1");
	comargs.push_back("3");
	comargs.push_back("600");
	comargs.push_back("null");
	comargs.push_back("null");
	comargs.push_back("0");
	comargs.push_back("2.2");
	comargs.push_back("1.0");
	comargs.push_back("null");
	comargs.push_back("null");
	comargs.push_back("0.000009");
}

void removeComOption(unsigned long &argc, vector<string> &comargs)
{
	// stupidly long list of possible arguments, but can't think of a better way to check this.
	if(comargs[1] == "-d" || comargs[1] == "-D" ||  comargs[1] == "-dl" ||  comargs[1] == "-dL" ||  comargs[1] == " -Dl" ||  comargs[1] == "-DL" ||
		comargs[1] == "-dx" || comargs[1] == "-dX" ||  comargs[1] == "-DX" ||  comargs[1] == " -Dx" ||  comargs[1] == "-c" ||  comargs[1] == "-C" ||
		comargs[1] == "-config" ||  comargs[1] == "-Config" || comargs[1] == "-f" || comargs[1] == "-h" || comargs[1] == "-H" || comargs[1] == "-F)")
	{
		comargs.erase(comargs.begin() + 1);
		argc --;
	}
	return;
}




