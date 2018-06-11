// This file is part of NECSim project which is released under MIT license.
// See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details
/**
 * @author Samuel Thompson
 * @file Map.h
 * @brief Contains Map for importing .tif files and obtaining a variety of information from them.
 *
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 */

#ifndef MAP_H
#define MAP_H
#ifdef with_gdal

#include <string>
#include <cstring>
#include <sstream>
#include "Logger.h"
#include <gdal_priv.h>
#include <cpl_conv.h> // for CPLMalloc()
#include <sstream>
#include "Matrix.h"
#include "Logger.h"
#include "CustomExceptions.h"
#include "CPLCustomHandlerNecsim.h"

using namespace std;
#ifdef DEBUG
#include "CustomExceptions.h"
#endif // DEBUG

/**
 * @brief Read a a tif file to a matrix and obtain spatial metadata.
 * @tparam T The type of the Matrix to create.
 */
template<class T>
class Map : public virtual Matrix<T>
{
protected:
	GDALDataset *poDataset;
	GDALRasterBand *poBand;
	unsigned long blockXSize, blockYSize;
	double noDataValue;
	string filename;
	GDALDataType dt{};
	// Contains the error object to check for any gdal issues
	CPLErr cplErr{CE_None};
	double upper_left_x{}, upper_left_y{}, x_res{}, y_res{};
	using Matrix<T>::matrix;
	using Matrix<T>::numCols;
	using Matrix<T>::numRows;
public:
	using Matrix<T>::setSize;
	using Matrix<T>::getCols;
	using Matrix<T>::getRows;
	using Matrix<T>::operator*;
	using Matrix<T>::operator+;
	using Matrix<T>::operator-;
	using Matrix<T>::operator-=;
	using Matrix<T>::operator+=;
	using Matrix<T>::operator[];

	Map() : Matrix<T>(0, 0)
	{
		GDALAllRegister();
		poDataset = nullptr;
		poBand = nullptr;
		filename = "";
		blockXSize = 0;
		blockYSize = 0;
		noDataValue = 0.0;
		CPLSetErrorHandler(cplNecsimCustomErrorHandler);
	}

	~Map()
	{
		close();
	}

	/**
	 * @brief Opens the provided filename to the poDataset object.
	 * @param filename file to open in read-only mode.
	 */
	void open(const string &filename_in)
	{
		if(!poDataset)
		{
			filename = filename_in;
			poDataset = (GDALDataset *) GDALOpen(filename.c_str(), GA_ReadOnly);
		}
		else
		{
			throw FatalException("File already open at " + filename);
		}
		if(!poDataset)
		{
			string s = "File " + filename + " not found.";
			throw runtime_error(s);
		}
	}

	/**
	 * @brief Overloaded open for using the preset file name.
	 */
	void open()
	{
		open(filename);
	}

	/**
	 * @brief Checks if the connection to the map file has already been opened.
	 *
	 * All this does is check if poDataset is a null pointer.
	 * @return true if poDataset is a null pointer.
	 */
	bool isOpen()
	{
		return poDataset != nullptr;
	}

	/**
	 * @brief Destroys the connection to the dataset.
	 */
	void close()
	{
		if(poDataset)
		{
			GDALClose(poDataset);
//			if(poDataset)
//			{
//				throw FatalException("poDataset not nullptr after closing, please report this bug.");
//			}
			poDataset = nullptr;
			poBand = nullptr;
		}
	}

	/**
	 * @brief Sets the raster band to the first raster.
	 */
	void getRasterBand()
	{
		poBand = poDataset->GetRasterBand(1);
	}

	/**
	 * @brief Obtains the x and y dimensions from the tif file for reading in blocks.
	 */
	void getBlockSizes()
	{
		blockXSize = static_cast<unsigned long>(poDataset->GetRasterXSize());
		blockYSize = static_cast<unsigned long>(poDataset->GetRasterYSize());
	}

	/**
	 * @brief Sets the no data, data type and data type name values from the tif file.
	 */
	void getMetaData()
	{
		try
		{
			int pbSuccess;
			noDataValue = poBand->GetNoDataValue(&pbSuccess);
			if(!pbSuccess)
			{
				noDataValue = 0.0;
			}
		}
		catch(out_of_range &out_of_range1)
		{
			noDataValue = 0.0;
		}
		stringstream ss;
		ss << "No data value is: " << noDataValue << endl;
		writeInfo(ss.str());
		// Check sizes match
		dt = poBand->GetRasterDataType();
		double geoTransform[6];
		cplErr = poDataset->GetGeoTransform(geoTransform);
		if(cplErr >= CE_Warning)
		{
			CPLError(cplErr, 6, "No transform present in dataset for %s.", filename.c_str());
			CPLErrorReset();
		}
		upper_left_x = geoTransform[0];
		upper_left_y = geoTransform[3];
		x_res = geoTransform[1];
		y_res = -geoTransform[5];
//		checkTifImportFailure();
#ifdef DEBUG
		printMetaData();
#endif // DEBUG
	}

#ifdef DEBUG
	void printMetaData()
	{
		stringstream ss;
		const char *dt_name = GDALGetDataTypeName(dt);
		ss << "Filename: " << filename << endl;
		writeLog(10, ss.str());
		ss.str("");
		ss << "data type: " << dt << "(" << dt_name << ")" << endl;
		writeLog(10, ss.str());
		ss.str("");
		ss << "Geo-transform (ulx, uly, x res, y res): " << upper_left_x << ", " << upper_left_y << ", ";
		ss << x_res << ", " << y_res << ", " <<endl;
		writeLog(10, ss.str());
		ss.str("");
		ss << "No data value: " << noDataValue << endl;
		writeLog(10, ss.str());

	}
#endif //DEBUG

	/**
	 * @brief Gets the upper left x (longitude) coordinate
	 * @return upper left x of the map
	 */
	double getUpperLeftX()
	{
		return upper_left_x;
	}

	/**
	 * @brief Gets the upper left y (latitude) coordinate
	 * @return upper left y of the map
	 */
	double getUpperLeftY()
	{
		return upper_left_y;
	}

	/**
	 * @brief Imports the matrix from a csv file.
	 *
	 * @throws runtime_error: if type detection for the filename fails.
	 * @param filename the file to import.
	 */
	void import(const string &filename) override
	{
		if(!importTif(filename))
		{
			Matrix<T>::import(filename);
		}
	}

	/**
	 * @brief Imports the matrix from a tif file using the gdal library functions.
	 * @note Opens a connection to the file object, which should be closed.
	 * @param filename the path to the file to import.
	 */
	bool importTif(const string &filename)
	{

		if(filename.find(".tif") != string::npos)
		{
			stringstream ss;
			ss << "Importing " << filename << " " << flush;
			writeInfo(ss.str());
			open(filename);
			getRasterBand();
			getBlockSizes();
			getMetaData();
			// If the sizes are 0 then use the raster sizes
			if(numCols == 0 || numRows == 0)
			{
				setSize(blockYSize, blockXSize);
			}
			// Check sizes
			if((numCols != blockXSize || numRows != blockYSize) || numCols == 0 ||
			   numRows == 0)
			{
				stringstream stringstream1;
				stringstream1 << "Raster data size does not match inputted dimensions for " << filename
							  << ". Using raster sizes."
							  << endl;
				stringstream1 << "Old dimensions: " << numCols << ", " << numRows << endl;
				stringstream1 << "New dimensions: " << blockXSize << ", " << blockYSize << endl;
				writeWarning(stringstream1.str());
				setSize(blockYSize, blockXSize);
			}
			// Check the data types are support
			const char *dt_name = GDALGetDataTypeName(dt);
			if(dt == 0 || dt > 7)
			{
				throw FatalException("Data type of " + string(dt_name) + " is not supported.");
			}
#ifdef DEBUG
			if(sizeof(T) * 8 != gdal_data_sizes[dt])
			{
				stringstream ss2;
				ss2 << "Object data size: " << sizeof(T) * 8 << endl;
				ss2 << "Tif data type: " << dt_name << ": " << gdal_data_sizes[dt] << " bytes" << endl;
				ss2 << "Tif data type does not match object data size in " << filename << endl;
				writeWarning(ss2.str());
			}
#endif
			// Just use the overloaded method for importing between types
			internalImport();
			writeInfo("done!\n");
			return true;
		}
		return false;
	}

	/**
	 * @brief Opens the offset map and fetches the metadata.
	 * @param offset_map the offset map to open (should be the larger map).
	 * @return true if the offset map is opened within this function
	 */
	bool openOffsetMap(Map &offset_map)
	{
		bool opened_here = false;
		if(!offset_map.isOpen())
		{
			opened_here = true;
			offset_map.open();
		}
		offset_map.getRasterBand();
		offset_map.getMetaData();
		return opened_here;
	}

	void closeOffsetMap(Map &offset_map, const bool &opened_here)
	{
		if(opened_here)
		{
			offset_map.close();
		}
	}

	/**
	 * @brief Calculates the offset between the two maps.
	 *
	 * The offset_map should be larger and contain this map, otherwise returned values will be negative
	 *
	 * @note Opens a connection to the tif file (if it has not already been opened), which is then closed. If the
	 * 		 connection is already open, then it will not be closed and it is assumed logic elsewhere achieves this.
	 *
	 * @note Offsets are returned as rounded integers at the resolution of the smaller map.
	 *
	 * @param offset_map the offset map to read from
	 * @param offset_x the x offset variable to fill
	 * @param offset_y the y offset variable to fill
	 */
	void calculateOffset(Map &offset_map, long &offset_x, long &offset_y)
	{
		auto opened_here = openOffsetMap(offset_map);
		offset_x = static_cast<long>(round((upper_left_x - offset_map.upper_left_x) / x_res));
		offset_y = static_cast<long>(round((offset_map.upper_left_y - upper_left_y )/ y_res));
		closeOffsetMap(offset_map, opened_here);
	}

	/**
	 * @brief Calculates the relative scale of this map compared to the offset map.
	 *
	 * The offset map should be larger and contain this map.
	 *
	 * @note Only the x resolution is checked, it is assumed the x and y resolutions of both maps is the same (i.e. each
	 * 		 cell on the map is a square.
	 *
	 * @param offset_map the offset map object to read from
	 * @return the relative scale of the offset map
	 */
	unsigned long roundedScale(Map &offset_map)
	{
		auto opened_here = openOffsetMap(offset_map);
		closeOffsetMap(offset_map, opened_here);
		return static_cast<unsigned long>(floor(offset_map.x_res / x_res));
	}

	/**
	 * @brief Default importer when we rely on the default gdal method of converting between values.
	 * Note that importing doubles to ints results in the values being rounded down.
	 * @return true if a tif file exists and can be imported, false otherwise.
	 */
	void internalImport()
	{
		writeWarning("No type detected for Map type. Attempting default importing (potentially undefined behaviour).");
		defaultImport();
	}

	/**
	 * @brief Default import routine for any type. Provided as a separate function so implementation can be called from
	 * any template class type.
	 */
	void defaultImport()
	{
		unsigned int number_printed = 0;
		for(uint32_t j = 0; j < numRows; j++)
		{
			printNumberComplete(j, number_printed);
			cplErr = poBand->RasterIO(GF_Read, 0, j, static_cast<int>(blockXSize), 1, &matrix[j][0],
									  static_cast<int>(blockXSize), 1, dt, 0, 0);
			checkTifImportFailure();
			// Now convert the no data values to 0
			for(uint32_t i = 0; i < numCols; i++)
			{
				if(matrix[j][i] == noDataValue)
				{
					matrix[j][i] = 0;
				}
			}
		}
	}

	/**
	 * @brief Imports from the supplied filename into the GeoTiff object, converting doubles to booleans.
	 * The threshold for conversion is x>0.5 -> true, false otherwise.
	 */
	void importFromDoubleAndMakeBool()
	{
		unsigned int number_printed = 0;
		// create an empty row of type float
		double *t1;
		t1 = (double *) CPLMalloc(sizeof(double) * numCols);
		// import the data a row at a time, using our template row.
		for(uint32_t j = 0; j < numRows; j++)
		{
			printNumberComplete(j, number_printed);
			cplErr = poBand->RasterIO(GF_Read, 0, j, static_cast<int>(blockXSize), 1, &t1[0],
									  static_cast<int>(blockXSize), 1, GDT_Float64, 0, 0);
			checkTifImportFailure();
			// now copy the data to our Map, converting float to int. Round or floor...? hmm, floor?
			for(unsigned long i = 0; i < numCols; i++)
			{
				if(t1[i] == noDataValue)
				{
					matrix[j][i] = false;
				}
				else
				{
					matrix[j][i] = t1[i] >= 0.5;
				}
			}
		}
		CPLFree(t1);
	}

	/**
	 * @brief Imports from the supplied filename into the GeoTiff object, converting doubles to booleans.
	 * The threshold for conversion is x>0.5 -> true, false otherwise.
	 *
	 * @param dt_buff: the buffer type for the data
	 * @tparam T2 the template type for data reading.
	 */
	template<typename T2>
	void importUsingBuffer(GDALDataType dt_buff)
	{
		unsigned int number_printed = 0;
		// create an empty row of type float
		T2 *t1;
		t1 = (T2 *) CPLMalloc(sizeof(T2) * numCols);
		// import the data a row at a time, using our template row.
		for(uint32_t j = 0; j < numRows; j++)
		{
			printNumberComplete(j, number_printed);
			cplErr = poBand->RasterIO(GF_Read, 0, j, static_cast<int>(blockXSize), 1, &t1[0],
									  static_cast<int>(blockXSize), 1, dt_buff, 0, 0);
			checkTifImportFailure();
			// now copy the data to our Map, converting float to int. Round or floor...? hmm, floor?
			for(unsigned long i = 0; i < numCols; i++)
			{
				if(t1[i] == noDataValue)
				{
					matrix[j][i] = static_cast<T>(0);
				}
				else
				{
					matrix[j][i] = static_cast<T>(t1[i]);
				}
			}
		}
		CPLFree(t1);
	}

	/**
	 * @brief Print the percentage complete during import
	 * @param j the reference for the counter
	 * @param number_printed the number of previously printed lines
	 */
	void printNumberComplete(const uint32_t &j, unsigned int &number_printed)
	{
		double dComplete = ((double) j / (double) numRows) * 20;
		if(number_printed < dComplete)
		{
			stringstream os;
			os << "\rImporting " << filename << " ";
			number_printed = 0;
			while(number_printed < dComplete)
			{
				os << ".";
				number_printed++;
			}
			os << flush;
			writeInfo(os.str());
		}
	}

	/**
	 * @brief Checks the error code of the CPLErr object and formats the error
	 */
	void checkTifImportFailure()
	{
		if(cplErr >= CE_Warning)
		{
			CPLError(cplErr, 3, "CPL error thrown during import of %s\n", filename.c_str());
			CPLErrorReset();
		}
	}

	friend ostream &operator>>(ostream &os, const Map &m)
	{
		return Matrix<T>::writeOut(os, m);
	}

	friend istream &operator<<(istream &is, Map &m)
	{
		return Matrix<T>::readIn(is, m);
	}

};

/**
 * @brief Overloaded imported for handling conversion of types to boolean. This function should only be once
 * elsewhere, so inlining is fine, allowing this file to remain header only.
 */
template<>
inline void Map<bool>::internalImport()
{
	if(dt <= 7)
	{
		// Then the tif file type is an int/byte
		// we can just import as it is
		importUsingBuffer<uint8_t>(GDT_Byte);
	}
	else
	{
		// Conversion from double to boolean
		importFromDoubleAndMakeBool();
	}
}

/**
 * @brief Overloaded functions for importing from tifs and matching between gdal and C types.
 * @param filename the path to the filename
 * @param poBand the GDALRasterBand pointer to import from
 * @param nBlockXSize the number of elements per row
 * @param dt the datatype (not required, exists for function overloading)
 * @param r the error reference object
 * @param ndv the no data value
 */
template<>
inline void Map<int8_t>::internalImport()
{
	importUsingBuffer<int16_t>(GDT_Int16);
}

template<>
inline void Map<uint8_t>::internalImport()
{
	dt = GDT_Byte;
	defaultImport();
}

template<>
inline void Map<int16_t>::internalImport()
{
	dt = GDT_Int16;
	defaultImport();
}

template<>
inline void Map<uint16_t>::internalImport()
{
	dt = GDT_UInt16;
	defaultImport();
}

template<>
inline void Map<int32_t>::internalImport()
{
	dt = GDT_Int32;
	defaultImport();
}

template<>
inline void Map<uint32_t>::internalImport()
{
	dt = GDT_UInt32;
	defaultImport();
}

template<>
inline void Map<float>::internalImport()
{
	dt = GDT_Float32;
	defaultImport();
}

template<>
inline void Map<double>::internalImport()
{
	dt = GDT_Float64;
	defaultImport();
}

#endif // with_gdal
#endif //MAP_H
