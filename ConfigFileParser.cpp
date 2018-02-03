//This file is part of NECSim project which is released under BSD-3 license.
//See file **LICENSE.txt** or visit https://opensource.org/licenses/BSD-3-Clause) for full license details.
// 
/**
 * @author Sam Thompson
 * @date 31/08/2016
 * @file Config.cpp
 * @brief Contains implementation of the Config.h functions.
 * 
 * Contact: samuel.thompson14@imperial.ac.uk or thompsonsed@gmail.com
 * @copyright <a href="https://opensource.org/licenses/BSD-3-Clause">BSD-3 Licence.</a>
 *
 */

#include "ConfigFileParser.h"
#include "CustomExceptions.h"
#include "Logging.h"

void importArgs(const unsigned int &argc, char *argv[], vector<string> &comargs)
{
	for(unsigned int i = 0; i < argc; i++)
	{
		comargs.emplace_back(argv[i]);
	}
	// check size is correct
	if(comargs.size() != argc)
	{
		cerr << "ERROR_MAIN_010: Incorrect command line parsing." << endl;
	}
}

string SectionOption::getOption(string refval)
{
	for(unsigned int i = 0; i < refs.size(); i++)
	{
		if(refs[i] == refval)
		{
			return (val[i]);
		}
	}
#ifdef DEBUG
	stringstream ss;
		ss << "Reference " << refval << " not found in keyoption." << endl;
	writeInfo(ss.str());
#endif
	return ("null");
}

ostream &operator<<(ostream &os, const SectionOption &k)
{
	os << k.section << "\n" << k.val.size() << "\n" << k.refs.size() << "\n";
	for(const auto &i : k.val)
	{
		os << i << "\n";
	}
	for(const auto &ref : k.refs)
	{
		os << ref << "\n";
	}
	return os;
}

istream &operator>>(istream &is, SectionOption &k)
{
	// os << m.numRows<<" , "<<m.numCols<<" , "<<endl;
	unsigned int valsize, refsize;
	is >> k.section >> valsize >> refsize;
	is.ignore();
	string tmp;
	for(unsigned int i = 0; i < valsize; i++)
	{
		getline(is, tmp);
		k.val.push_back(tmp);
	}
	for(unsigned int i = 0; i < refsize; i++)
	{
		getline(is, tmp);
		k.refs.push_back(tmp);
	}
	return is;
}

void ConfigOption::setConfig(const string &file, bool main, bool full_parse)
{
	if(!bConfig)
	{
		bMain = main;
		configfile = file;
		bConfig = true;
		bFullParse = full_parse;
	}
	else
	{
		throw ConfigException("Attempt to set config file twice.");
	}
}

void ConfigOption::parseConfig()
{
	ifstream is_file;
	try
	{
		is_file.open(configfile);
	}
	catch(...)
	{
		throw ConfigException(
				"ERROR_CONF_004c: Could not open the config file. Check file exists and is readable.");
	}
	if(!is_file.fail() || !is_file.good())
	{
		string line;
		// Get the first line of the file.
		while(getline(is_file, line))
		{
//				os << line << endl;
			istringstream is_line(line);
			string key;
			string val;
			// Skip all whitespace
			is_line >> skipws;
			// start a new section
			if(line[0] == '[')
			{
				SectionOption tempSections;
				// get the section name
				string section;
				if(getline(is_line, section, ']'))
				{
					section = section.erase(0, 1);
					tempSections.section = section;
//						os << section << endl;
				}
				// read each line
				while(getline(is_file, line))
				{
					// end the section when a new one starts.
					if(line[0] == '[' || line.size() == 0)
					{
						break;
					}
					istringstream is_line2(line); // update the input-line stream
					if(getline(is_line2, key, '='))
					{

						key.erase(std::remove(key.begin(), key.end(), ' '), key.end());
						is_line2 >> skipws;
					}
					if(!is_line2)
					{
//							os << is_line2 << endl;
						throw ConfigException("ERROR_CONF_001: Read error in config file.");
					}
					if(getline(is_line2, val))
					{
//						 	This line has been removed to allow for white spaces in file names and paths
//							val.erase(std::remove(val.begin(), val.end(), ' '), val.end());
						while(val[0] == ' ')
						{
							val.erase(val.begin(), val.begin() + 1);
						}

					}
					if(!is_line2)
					{
						throw ConfigException("ERROR_CONF_001: Read error in config file.");
					}
					tempSections.refs.push_back(key);
					tempSections.val.push_back(val);
				}
				configs.push_back(tempSections);
			}
		}
	}
	else
	{
		throw ConfigException(
				"ERROR_CONF_004b: Could not open the config file " + configfile +
				". Check file exists and is readable.");
	}
	if(is_file.eof())
	{
		is_file.close();
	}
	else
	{
		throw ConfigException("ERROR_CONF_002: End of file not reached. Check input file formating.");
	}
}

vector<SectionOption> ConfigOption::getSectionOptions()
{
	return configs;
}

SectionOption ConfigOption::operator[](int index)
{
	return (configs[index]);
}

unsigned int ConfigOption::getSectionOptionsSize()
{
	return (configs.size());
}

vector<string> ConfigOption::getSections()
{
	vector<string> toret;
	for(unsigned i = 0; i < configs.size(); i++)
	{
		toret.push_back(configs[i].section);
	}
	return (toret);
}

bool ConfigOption::hasSection(string sec)
{
	for(unsigned i = 0; i < configs.size(); i++)
	{
		if(configs[i].section == sec)
		{
			return (true);
		}
	}
	return (false);
}

vector<string> ConfigOption::getSectionValues(string sec)
{
	for(unsigned i = 0; i < configs.size(); i++)
	{
		if(configs[i].section == sec)
		{
			return (configs[i].val);
		}
	}
	throw ConfigException("Section not found in config file: " + sec);
}

string ConfigOption::getSectionOptions(string section, string ref)
{
	for(unsigned int i = 0; i < configs.size(); i++)
	{
		if(configs[i].section == section)
		{
			for(unsigned int j = 0; j < configs[i].refs.size(); j++)
			{
				if(configs[i].refs[j] == ref)
				{
					return (configs[i].val[j]);
				}
			}
		}
	}
#ifdef DEBUG
	writeWarning("No reference found for " + section + ", ");
#endif
	return ("null");
}

string ConfigOption::getSectionOptions(string section, string ref, string def)
{
	for(unsigned int i = 0; i < configs.size(); i++)
	{
		if(configs[i].section == section)
		{
			for(unsigned int j = 0; j < configs[i].refs.size(); j++)
			{
				if(configs[i].refs[j] == ref)
				{
					return (configs[i].val[j]);
				}
			}
		}
	}
	return def;
}

int ConfigOption::importConfig(vector<string> &comargs)
{
	// Check that the previous arguments have already been imported.
	if(bMain)
	{
		if(comargs.size() != 3)
		{
			throw ConfigException(
					"ERROR_CONF_003: Number of command line arguments not correct before importSpatialParameters.");
		}
	}
	ifstream is_file;
	try
	{
		is_file.open(configfile);
	}
	catch(...)
	{
		throw ConfigException(
				"ERROR_CONF_004a: Could not open the config file. Check file exists and is readable.");
	}
	if(!is_file.fail())
	{
		string line;
		while(getline(is_file, line))
		{
			istringstream is_line(line);
			string key;
			is_line >> skipws;
			if(line[0] == '[')
			{
				continue;
			}
			if(getline(is_line, key, '='))
			{
				// Could implement proper data parsing based on the key object.
				is_line >> skipws;
				string value;
				if(getline(is_line, value))
				{
					value.erase(std::remove(value.begin(), value.end(), ' '), value.end());
					if(!is_line)
					{
						stringstream os;
						os << value << endl;
						writeWarning(os.str());
						throw ConfigException("ERROR_CONF_001: Read error in config file.");
					}
					char *tmp = new char[value.length() + 1];
					strcpy(tmp, value.c_str());
					comargs.push_back(tmp);
				}
			}
		}
	}
	else
	{
		throw ConfigException(
				"ERROR_CONF_004d: Could not open the config file. Check file exists and is readable.");
	}
	if(is_file.eof())
	{
		is_file.close();
	}
	else
	{
		throw ConfigException("ERROR_CONF_002: End of file not reached. Check input file formating.");
	}
	if(bMain)
	{
		// remove the file name from the command line arguments to maintain the vector format.
		comargs.erase(comargs.begin() + 2);
	}
	return comargs.size();
}

ostream &operator<<(ostream &os, const ConfigOption &c)
{
	os << c.configfile << "\n" << c.bConfig << "\n" << c.bMain << "\n" << c.bFullParse << "\n" << c.configs.size()
	   << "\n";
	for(unsigned int i = 0; i < c.configs.size(); i++)
	{
		os << c.configs[i];
	}
	return os;
}

istream &operator>>(istream &is, ConfigOption &c)
{
	unsigned int configsize;
	is.ignore();
	getline(is, c.configfile);
	is >> c.bConfig >> c.bMain >> c.bFullParse >> configsize;
//		os << "file: " << c.configfile << endl;
//		os << "bconf: " << c.bConfig << endl;
//		os << "bmain: " << c.bMain << endl;
//		os << "fullp: " << c.bFullParse << endl;
	SectionOption tmpoption;
//		 os << "configsize: " << configsize << endl;
//		cout << os.str() << endl;
	// check that the config size isn't completely stupid!
	if(configsize > 10000)
	{
		throw runtime_error("Config size extremely large, check file: " + to_string(configsize));
	}
	if(configsize > 0)
	{
		for(unsigned int i = 0; i < configsize; i++)
		{
			is >> tmpoption;
			c.configs.push_back(tmpoption);
		}
	}
//		os << "end config" << endl;
	return is;
}
