//This file is part of NECSim project which is released under MIT license.
//See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.

/**
 * @author Samuel Thompson
 * @file SpeciesList.h
 * @brief  Contains the SpeciesList class for usage in coalescence simulations.
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 */
#include <iostream>
#include "SpeciesList.h"


SpeciesList::SpeciesList() : list_size(0), maxsize(0), next_active(0), nwrap(0)
{
	list.setSize(0);
}

void SpeciesList::fillList()
{
	if(maxsize!=0)
	{
		for(unsigned int i = 0; i < maxsize; i++)
		{
			list[i] = 0;
		}
	}
}

void SpeciesList::initialise(unsigned long maxsizein)
{
	maxsize = maxsizein;
	nwrap = 0;
	list_size = 0;
	list.setSize(maxsize);
}

void SpeciesList::setMaxsize(unsigned long maxsizein)
{
	maxsize = maxsizein;
}

void SpeciesList::setSpecies(unsigned long index, unsigned long new_val)
{
	if(list[index] == 0)
	{
		throw runtime_error("ERROR_MOVE_027: List position to be replaced is zero. Check list assignment.");
	}
	list[index] = new_val;
}

void SpeciesList::setSpeciesEmpty(unsigned long index, unsigned long new_val)
{
	if(list[index] != 0)
	{
		throw runtime_error("ERROR_MOVE_027b: List position to be replaced is not zero. Check list assignment.");
	}
	list[index] = new_val;
}

void SpeciesList::setNext(unsigned long n)
{
	next_active = n;
}

void SpeciesList::setNwrap(unsigned long nr)
{
	nwrap = nr;
}

unsigned long SpeciesList::addSpecies(unsigned long new_spec)
{
#ifdef DEBUG
	if(list_size + 1 > maxsize)
	{
		stringstream ss;
		ss << "list size: " << list_size << endl;
		ss << "max size: " << maxsize << endl;
		writeLog(10, ss.str());
		throw out_of_range("Could not add species - list size greater than max size.");
	}
#endif
	// First loop from the list size value
	for(unsigned long i = list_size; i<list.size(); i++)
	{
		if(list[i] == 0)
		{
			list_size++;
			list[i] = new_spec;
			return i;
		}
	}
	// Now loop over the rest of the lineages
	for(unsigned long i = 0; i < list_size; i ++)
	{
		if(list[i] == 0)
		{
			list_size++;
			list[i] = new_spec;
			return i;
		}
	}
	throw out_of_range("Could not add species - no empty space");
}

void SpeciesList::addSpeciesSilent(unsigned long new_spec)
{
	for(unsigned long i =0;i<maxsize;i++)
	{
		if(list[i] == 0)
		{
			list_size++;
			list[i] = new_spec;
			return;
		}

	}
	throw out_of_range("Could not add species - no empty space");
}

void SpeciesList::deleteSpecies(unsigned long index)
{
	list[index] = 0;
	list_size --;
}

void SpeciesList::decreaseNwrap()
{
	if(nwrap == 0)
	{
		throw runtime_error("Nwrap should never be decreased less than 0");
	}
	else if(nwrap == 1)
	{
		if( next_active != 0)
		{
			throw runtime_error("Nwrap is being set at 0 when an wrapped lineage is still present");
		}
	}
	nwrap --;
}

void SpeciesList::increaseListSize()
{
	list_size ++;
}

void SpeciesList::increaseNwrap()
{
	nwrap ++;
}

void SpeciesList::changePercentCover(unsigned long newmaxsize)
{
	Row<unsigned long> templist(list);
	maxsize = newmaxsize;
	list.setSize(newmaxsize);
	for(unsigned int i=0;i<max(newmaxsize, templist.size());i++)
	{
		if(i<templist.size())
		{
			list[i] = templist[i];
		}
		else
		{
			list[i] = 0;
		}
	}
#ifdef DEBUG
	if(list.size() > maxsize)
		{
			throw out_of_range("List size not equal to maxsize");
		}
#endif
}

unsigned long SpeciesList::getRandLineage(NRrand &rand_no)
{
	double rand_index;
	if(maxsize <= list_size)
	{
		// Then the list size is larger than the actual size. This means we must return a lineage.
		try
		{
			do
			{
				rand_index = rand_no.d01();
				rand_index *= list.size();
				//os << "ref: " << rand_index << ", " << list[round(rand_index)] << endl;
			} while(list[floor(rand_index)] == 0);
			//os << "RETURNING!" << endl;
			return(list[floor(rand_index)]);
		}
		catch(out_of_range &oor)
		{
			throw runtime_error("ERROR_MOVE_001b: Listpos outside maxsize.");
		}
	}
	else
	{
		rand_index =  rand_no.d01();
//		os << "rand_index: " << rand_index << endl;
		rand_index  *= maxsize;
		// Dynamically resize the list if required. Otherwise, to save memory, the list will not be resized;
		if(rand_index>=list.size())
		{
			changePercentCover(maxsize);
		}

		unsigned long i = static_cast<unsigned long>(floor(rand_index));

#ifdef DEBUG
		if(rand_index>maxsize)
			{
				stringstream ss;
				ss << "Random index is greater than the max size. Fatal error, please report this bug." << endl;
				throw runtime_error(ss.str());
			}
#endif // DEBUG
		return list[i];
	}
}

unsigned long SpeciesList::getSpecies(unsigned long index)
{
	return list[index];
}

unsigned long SpeciesList::getNext()
{
	return next_active;
}

unsigned long SpeciesList::getNwrap()
{
	return nwrap;
}

unsigned long SpeciesList::getListSize()
{
	return list_size;
}

unsigned long SpeciesList::getMaxSize()
{
	return maxsize;
}

unsigned long SpeciesList::getListLength()
{
	return list.size();
}

void SpeciesList::wipeList()
{
	fillList();
	next_active=0;
	nwrap =0;
	list_size=0;
}

ostream &operator<<(ostream &os, const SpeciesList &r)
{
	//os << m.numRows<<" , "<<m.numCols<<" , "<<endl;
	os << r.list << ",";
	os << r.list_size << ",";
	os << r.maxsize << ",";
	os << r.next_active << ",";
	os << r.nwrap << ",";
	return os;
}

istream &operator>>(istream &is, SpeciesList &r)
{
	char delim;
	//double temp1,temp2;
	//is << m.numRows<<" , "<<m.numCols<<" , "<<endl;
	is >> r.list;
	is >> delim;
	is >> r.list_size;
	is >> delim;
	is >> r.maxsize;
	is >> delim;
	is >> r.next_active;
	is >> delim;
	is >> r.nwrap;
	is >> delim;
	return is;
}
