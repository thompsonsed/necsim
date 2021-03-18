//This file is part of necsim project which is released under MIT license.
//See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.

/**
 * @author Samuel Thompson
 * @file Matrix.h
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 * @brief Contains a template for a matrix with all the basic matrix operations overloaded.
 *
 * @details Provides an efficient, general purpose 2D matrix object with an efficient indexing system designed for
 * modern CPUs (where memory access times are often much longer than compute times for mathematical operations).
 * Most operations are low-level, but some higher level functions remain, such as importCsv().
 *
 * Contact: thompsonsed@gmail.com
 */


#ifndef MATRIX
#define MATRIX
#define null 0

#include <cstdio>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <vector>

#ifdef use_csv
#include<cmath>
#include <stdexcept>
#include "fast-cpp-csv-parser/csv.h"
#endif

#include <cstdint>
#include "Logging.h"

using std::vector;
namespace necsim
{
    // Array of data sizes for importing tif files.
    const int gdal_data_sizes[] = {0, 8, 16, 16, 32, 32, 32, 64};

    /**
     * @brief A class containing the Matrix object, set up as an array of Row objects.
     * Includes basic operations, as well as the importCsv() function for more advanced reading from file.
     * @tparam T the type of the values in the matrix
     */
    template<class T>
    class Matrix
    {

    protected:

        // number of rows and columns
        unsigned long num_cols{};
        unsigned long num_rows{};
        // a matrix is an array of rows
        vector<T> matrix;
    public:

        Matrix() : num_cols(0), num_rows(0), matrix()
        {

        }

        /**
         * @brief The standard constructor
         * @param rows optionally provide the number of rows.
         * @param cols optionally provide the number of columns.
         */
        explicit Matrix(unsigned long rows, unsigned long cols) : num_cols(cols), num_rows(rows),
                                                                  matrix(rows * cols, T())
        {
        }

        Matrix(Matrix &&m) noexcept
        {
            *this = std::move(m);
        }

        /**
         * @brief The copy constructor.
         * @param m a Matrix object to copy from.
         */
        Matrix(const Matrix &m) noexcept: num_cols(0), num_rows(0), matrix()
        {
            *this = m;
        }

        /**
        * @brief The destructor.
        */
        virtual ~Matrix() = default;

        /**
         * @brief Sets the matrix size.
         * Similar concept to that for Rows.
         * @param rows the number of rows.
         * @param cols the number of columns.
         */
        void setSize(unsigned long rows, unsigned long cols)
        {
            if(!matrix.empty())
            {
                matrix.clear();
            }
            matrix.resize(rows * cols);
            num_cols = cols;
            num_rows = rows;
        }

        /**
         * @brief Getter for the number of columns.
         * @return the number of columns.
         */
        unsigned long getCols() const
        {
            return num_cols;
        }

        /**
         * @brief Getter for the number of rows.
         * @return the number of rows.
         */
        unsigned long getRows() const
        {
            return num_rows;
        }

        /**
         * @brief Fills the matrix with the given value
         * @param val the value to fill
         */
        void fill(T val)
        {
            std::fill(matrix.begin(), matrix.end(), val);
        }

        /**
         * @brief Gets the index of a particular row and column in the matrix.
         * @param row the row number to index
         * @param col the column number to index
         * @return the index of row and column within the matrix
         */
        unsigned long index(const unsigned long &row, const unsigned long &col) const
        {
#ifdef DEBUG
            if(num_cols == 0 || num_rows == 0)
            {
                throw std::out_of_range ("Matrix has 0 rows and columns for indexing from.");
            }
            if(col + num_cols * row > matrix.size())
            {
                std::stringstream ss;
                ss << "Index of " << col + num_cols * row << ", (" << row << ", " << col << ")";
                ss << " is out of range of matrix vector with size " << matrix.size() << std::endl;
                throw std::out_of_range (ss.str());
            }
#endif // DEBUG
            return col + num_cols * row;
        }

        /**
         * @brief Gets the value at a particular index.
         * @param row the row number to get the value at
         * @param col the column number to get the value at
         * @return the value at the specified row and column
         */
        T &get(const unsigned long &row, const unsigned long &col)
        {
#ifdef DEBUG
            if(row < 0 || row >= num_rows || col < 0 || col >= num_cols)
            {
                std::stringstream ss;
                ss << "Index of " << row << ", " << col << " is out of range of matrix with size " << num_rows;
                ss << ", " << num_cols << std::endl;
                throw std::out_of_range (ss.str());
            }
#endif
            return matrix[index(row, col)];
        }

        /**
         * @brief Gets the value at a particular index.
         * @param row the row number to get the value at
         * @param col the column number to get the value at
         * @return the value at the specified row and column
         */
        const T &get(const unsigned long &row, const unsigned long &col) const
        {
#ifdef DEBUG
            if(row < 0 || row >= num_rows || col < 0 || col >= num_cols)
            {
                std::stringstream ss;
                ss << "Index of " << row << ", " << col << " is out of range of matrix with size " << num_rows;
                ss << ", " << num_cols << std::endl;
                throw std::out_of_range (ss.str());
            }
#endif
            return matrix[index(row, col)];
        }

        /**
         * @brief Gets the value at a particular index.
         * @param row the row number to get the value at
         * @param col the column number to get the value at
         * @return the value at the specified row and column
         */
        T getCopy(const unsigned long &row, const unsigned long &col) const
        {
#ifdef DEBUG
            if(row < 0 || row >= num_rows || col < 0 || col >= num_cols)
            {
                std::stringstream ss;
                ss << "Index of " << row << ", " << col << " is out of range of matrix with size " << num_rows;
                ss << ", " << num_cols << std::endl;
                throw std::out_of_range (ss.str());
            }
#endif
            return matrix[index(row, col)];
        }

        /**
         * @brief Returns iterators for range-based for loops.
         * @return iterator to the start of the vector
         */
        typename vector<T>::iterator begin()
        {
            return matrix.begin();
        }

        /**
         * @brief Returns end iterators for range-based for loops.
         * @return iterator to the end of the vector
         */
        typename vector<T>::iterator end()
        {
            return matrix.end();
        }

        /**
         * @brief Returns iterators for range-based for loops.
         * @return iterator to the start of the vector
         */
        typename vector<T>::const_iterator begin() const
        {
            return matrix.begin();
        }

        /**
         * @brief Returns end iterators for range-based for loops.
         * @return iterator to the end of the vector
         */
        typename vector<T>::const_iterator end() const
        {
            return matrix.end();
        }

        /**
         * @brief Gets the arithmetic mean of the Matrix
         * @return the mean value in the matrix
         */
        double getMean() const
        {
            if(matrix.empty())
            {
                return 0.0;
            }
            T total = sum();
            return double(total) / matrix.size();

        }

        T sum() const
        {
            if(matrix.empty())
            {
                return T();
            }
            else
            {
                T total = 0;
                for(const auto &item : matrix)
                {
                    total += item;
                }
                return total;
            }
        }

        /**
         * @brief Overloading the = operator.
         * @param m the matrix to copy from.
         */
        Matrix &operator=(const Matrix &m) noexcept
        {
            if(m.matrix.empty())
            {
                matrix.clear();
                num_cols = 0;
                num_rows = 0;
            }
            else
            {
                matrix = m.matrix;
                num_cols = m.num_cols;
                num_rows = m.num_rows;
            }
            return *this;
        }

        Matrix &operator=(Matrix &&m) noexcept
        {
            if(m.matrix.empty())
            {
                matrix.clear();
                num_cols = 0;
                num_rows = 0;
            }
            else
            {
                matrix = std::move(m.matrix);
                num_cols = m.num_cols;
                num_rows = m.num_rows;
            }
            return *this;
        }

        /**
         * @brief Overloading the + operator.
         * @note If matrices are of different sizes, the operation is performed on the 0 to minimum values of each
         *       dimension.
         * @param m the matrix to add to this matrix.
         * @return the matrix object which is the sum of the two matrices.
         */
        Matrix operator+(const Matrix &m) const
        {
            //Since addition creates a new matrix, we don't want to return a reference, but an actual matrix object.
            unsigned long new_num_cols = findMinCols(this, m);
            unsigned long new_num_rows = findMinRows(this, m);
            Matrix result(new_num_rows, new_num_cols);
            for(unsigned long r = 0; r < new_num_rows; r++)
            {
                for(unsigned long c = 0; c < new_num_cols; c++)
                {
                    result.get(r, c) = get(r, c) + m.get(r, c);
                }
            }
            return result;
        }

        /**
         * @brief Overloading the - operator.
         * @note If matrices are of different sizes, the operation is performed on the 0 to minimum values of each
         *       dimension.
         * @param m the matrix to subtract from this matrix.
         * @return the matrix object which is the subtraction of the two matrices.
         * */
        Matrix operator-(const Matrix &m) const
        {
            unsigned long new_num_cols = findMinCols(this, m);
            unsigned long new_num_rows = findMinRows(this, m);
            Matrix result(new_num_rows, new_num_cols);
            for(unsigned long r = 0; r < new_num_rows; r++)
            {
                for(unsigned long c = 0; c < new_num_cols; c++)
                {
                    result.get(r, c) = get(r, c) - m.get(r, c);
                }
            }
            return result;
        }

        /**
         * @brief Overloading the += operator so that the new object is written to the current object.
         * @note If matrices are of different sizes, the operation is performed on the 0 to minimum values of each
         *       dimension.
         * @param m the Matrix object to add to this matrix.
         */
        Matrix &operator+=(const Matrix &m)
        {
            unsigned long new_num_cols = findMinCols(this, m);
            unsigned long new_num_rows = findMinRows(this, m);
            for(unsigned long r = 0; r < new_num_rows; r++)
            {
                for(unsigned long c = 0; c < new_num_cols; c++)
                {
                    get(r, c) += m.get(r, c);
                }
            }
            return *this;
        }

        /**
         * @brief Overloading the -= operator so that the new object is written to the current object.
         * @note If matrices are of different sizes, the operation is performed on the 0 to minimum values of each
         *       dimension.
         * @param m the Matrix object to subtract from this matrix.
         */
        Matrix &operator-=(const Matrix &m)
        {
            unsigned long new_num_cols = findMinCols(this, m);
            unsigned long new_num_rows = findMinRows(this, m);
            for(unsigned long r = 0; r < new_num_rows; r++)
            {
                for(unsigned long c = 0; c < new_num_cols; c++)
                {
                    matrix.get(r, c) -= m.get(r, c);
                }
            }
            return *this;
        }

        /**
         * @brief Overloading the * operator for scaling.
         * @note If matrices are of different sizes, the operation is performed on the 0 to minimum values of each
         *       dimension.
         * @param s the constant to scale the matrix by.
         * @return the scaled matrix.
         */
        Matrix operator*(const double s) const
        {
            Matrix result(num_rows, num_cols);
            for(unsigned long r = 0; r < num_rows; r++)
            {
                for(unsigned long c = 0; c < num_cols; c++)
                {
                    result.get(r, c) = get(r, c) * s;
                }
            }
            return result;
        }

        /**
         * @brief Overloading the * operator for matrix multiplication.
         * @note If matrices are of different sizes, the operation is performed on the 0 to minimum values of each
         *       dimension.
         * Multiplies each value in the matrix with its corresponding value in the other matrix.
         * @param m the matrix to multiply with
         * @return the product of each ith,jth value of the matrix.
         */
        Matrix operator*(Matrix &m) const
        {
            unsigned long new_num_cols = findMinCols(this, m);
            unsigned long new_num_rows = findMinRows(this, m);

            Matrix result(num_rows, m.num_cols);
            for(unsigned long r = 0; r < new_num_rows; r++)
            {
                for(unsigned long c = 0; c < new_num_cols; c++)
                {
                    result.get(r, c) = get(r, c) * m.get(r, c);
                }
            }
            return result;
        }

        /**
         * @brief Overloading the *= operator so that the new object is written to the current object.
         * @note If matrices are of different sizes, the operation is performed on the 0 to minimum values of each
         *       dimension.
         * @param m the Matrix object to add to this matrix.
         */
        Matrix &operator*=(const double s)
        {
            for(unsigned long r = 0; r < num_rows; r++)
            {
                for(unsigned long c = 0; c < num_cols; c++)
                {
                    get(r, c) *= s;
                }
            }
            return *this;
        }

        /**
         * @brief Overloading the *= operator so that the new object is written to the current object.
         * @note If matrices are of different sizes, the operation is performed on the 0 to minimum values of each
         *       dimension.
         * @param m the Matrix object to add to this matrix.
         */
        Matrix &operator*=(const Matrix &m)
        {
            unsigned long new_num_cols = findMinCols(this, m);
            unsigned long new_num_rows = findMinRows(this, m);
            for(unsigned long r = 0; r < new_num_rows; r++)
            {
                for(unsigned long c = 0; c < new_num_cols; c++)
                {
                    get(r, c) *= m.get(r, c);
                }
            }
            return *this;
        }

        /**
         * @brief Overloading the / operator for scaling.
         * @note If matrices are of different sizes, the operation is performed on the 0 to minimum values of each
         *       dimension.
         * @param s the constant to scale the matrix by.
         * @return the scaled matrix.
         */
        Matrix operator/(const double s) const
        {
            Matrix result(num_rows, num_cols);
            for(unsigned long r = 0; r < num_rows; r++)
            {
                for(unsigned long c = 0; c < num_cols; c++)
                {
                    result.get(r, c) = get(r, c) / s;
                }
            }
            return result;
        }

        /**
         * @brief Overloading the /= operator so that the new object is written to the current object.
         * @note If matrices are of different sizes, the operation is performed on the 0 to minimum values of each
         *       dimension.
         * @param m the Matrix object to add to this matrix.
         */
        Matrix &operator/=(const double s)
        {
            for(unsigned long r = 0; r < num_rows; r++)
            {
                for(unsigned long c = 0; c < num_cols; c++)
                {
                    get(r, c) /= s;
                }
            }
            return *this;
        }

        /**
         * @brief Overloading the /= operator so that the new object is written to the current object.
         * @note If matrices are of different sizes, the operation is performed on the 0 to minimum values of each
         *       dimension.
         * @param m the Matrix object to add to this matrix.
         */
        Matrix &operator/=(const Matrix &m)
        {
            unsigned long new_num_cols = findMinCols(this, m);
            unsigned long new_num_rows = findMinRows(this, m);
            for(unsigned long r = 0; r < new_num_rows; r++)
            {
                for(unsigned long c = 0; c < new_num_cols; c++)
                {
                    get(r, c) /= m.get(r, c);
                }
            }
            return *this;
        }

        /**
         * @brief Writes the object to the output stream.
         * @note This is done slightly inefficiently to preserve the output taking the correct form.
         * @param os the output stream to write to
         * @param m the object to write out
         * @return the output stream
         */
        friend std::ostream &writeOut(std::ostream &os, const Matrix &m)
        {
            for(unsigned long r = 0; r < m.num_rows; r++)
            {
                for(unsigned long c = 0; c < m.num_cols; c++)
                {
                    os << m.getCopy(r, c) << ",";
                }
                os << "\n";
            }
            return os;
        }

        /**
         * @brief Reads in from the input stream.
         * @param is the input stream to read from
         * @param m the object to read into
         * @return
         */
        friend std::istream &readIn(std::istream &is, Matrix &m)
        {
            char delim;
            for(unsigned long r = 0; r < m.num_rows; r++)
            {
                for(unsigned long c = 0; c < m.num_cols; c++)
                {
                    is >> m.get(r, c);
                    is >> delim;
                }
            }
            return is;
        }

        /**
         * @brief Overloading the << operator for outputting to an output stream.
         * This can be used for writing to console or storing to file.
         * @param os the output stream.
         * @param m the matrix to output.
         * @return the output stream.
         */
        friend std::ostream &operator<<(std::ostream &os, const Matrix &m)
        {
            return writeOut(os, m);
        }

        /**
         * @brief Overloading the >> operator for inputting from an input stream.
         * This can be used for writing to console or storing to file.
         * @param is the input stream.
         * @param m the matrix to input to.
         * @return the input stream.
         */
        friend std::istream &operator>>(std::istream &is, Matrix &m)
        {
            return readIn(is, m);
        }

        /**
         * @brief Sets the value at the specified indices, including handling type conversion from char to the template
         * class.
         * @param row the row index.
         * @param col the column index.
         * @param value the value to set
         */
        void setValue(const unsigned long &row, const unsigned long &col, const char* value)
        {
            matrix[index(row, col)] = static_cast<T>(*value);
        }

        /**
         * @brief Sets the value at the specified indices, including handling type conversion from char to the template
         * class.
         * @param row the row index.
         * @param col the column index.
         * @param value the value to set
         */
        void setValue(const unsigned long &row, const unsigned long &col, const T &value)
        {
            matrix[index(row, col)] = value;
        }

        /**
         * @brief Imports the matrix from a csv file.
         *
         * @throws runtime_error: if type detection for the filename fails.
         * @param filename the file to import.
         */
        virtual void import(const string &filename)
        {
            if(!importCsv(filename))
            {
                string s = "Type detection failed for " + filename + ". Check file_name is correct.";
                throw std::runtime_error(s);
            }
        }

        /**
         * @brief Imports the matrix from a csv file using the fast-csv-parser method.
         * @param filename the path to the file to import.
         */
#ifdef use_csv
        bool importCsv(const string &file_name)
        {
        if(file_name.find(".csv") != string::npos)
            {
                std::stringstream os;
                os  << "Importing " << file_name << " " << std::flush;
                writeInfo(os.str());
                // LineReader option
                io::LineReader in(file_name);
                // Keep track of whether we've printed to terminal or not.
                bool bPrint = false;
                // Initialies empty variable so that the setValue operator overloading works properly.
                unsigned int number_printed = 0;
                for(unsigned long i =0; i<num_rows; i++)
                {
                    char* line = in.next_line();
                    if(line == nullptr)
                    {
                        if(!bPrint)
                        {
                            writeError("Input dimensions incorrect - read past end of file.");
                            bPrint = true;
                        }
                        break;
                    }
                    else
                    {
                        char *dToken;
                        dToken = strtok(line,",");
                        for(unsigned long j = 0; j<num_cols; j++)
                        {
                            if(dToken == nullptr)
                            {
                                if(!bPrint)
                                {
                                writeError("Input dimensions incorrect - read past end of file.");
                                    bPrint = true;
                                }
                                break;
                            }
                            else
                            {
                                // This function is overloaded to correctly determine the type of the template
                                setValue(i,j,dToken);
                                dToken = strtok(NULL,",");
                            }
                        }
                        // output the percentage complete
                        double dComplete = ((double)i/(double)num_rows)*20;
                        if( number_printed < dComplete)
                        {
                            std::stringstream os;
                            os  << "\rImporting " << file_name << " ";
                            number_printed = 0;
                            while(number_printed < dComplete)
                            {
                                os << ".";
                                number_printed ++;
                            }
                            os << std::flush;
                            writeInfo(os.str());
                        }

                    }
                }
                writeInfo("done.\n");
                return true;
            }
            return false;
        }
#endif
#ifndef use_csv

        /**
         * @brief Imports the matrix from a csv file using the standard, slower method.
         * @deprecated this function should not be used any more as it is much slower.
         * @param filename the path to the file to import.
         * @return true if the csv can be imported.
         */
        bool importCsv(const string &filename)
        {
            if(filename.find(".csv") != string::npos)
            {
                std::stringstream os;
                os << "Importing" << filename << " " << std::flush;
                std::ifstream inputstream;
                inputstream.open(filename.c_str());
                unsigned long number_printed = 0;
                for(uint32_t j = 0; j < num_rows; j++)
                {
                    string line;
                    getline(inputstream, line);
                    std::istringstream iss(line);
                    for(uint32_t i = 0; i < num_cols; i++)
                    {
                        char delim;
                        T val;
                        iss >> val >> delim;
                        this->setValue(j, i, val);
                    }
                    double dComplete = ((double) j / (double) num_rows) * 5;
                    if(number_printed < dComplete)
                    {
                        os << "\rImporting " << filename << " " << std::flush;
                        while(number_printed < dComplete)
                        {
                            os << ".";
                            number_printed++;
                        }
                        os << std::flush;
                        writeInfo(os.str());

                    }
                }
                std::stringstream os2;
                os2 << "\rImporting" << filename << "..." << "done." << "                          " << std::endl;
                inputstream.close();
                writeInfo(os2.str());
                return true;
            }
            return false;
        }

#endif // use_csv
    };

    /**
     * @brief Find the minimum columns of the two objects.
     * @tparam T The type of the Matrix class
     * @param matrix1 the first matrix
     * @param matrix2 the second matrix
     * @return the minimum number of columns between the two matrices
     */
    template<typename T> unsigned long findMinCols(const Matrix<T> &matrix1, const Matrix<T> &matrix2)
    {
        if(matrix1.getCols() < matrix2.getCols())
        {
            return matrix1.getCols();
        }
        return matrix2.getCols();
    }

    /**
     * @brief Find the minimum rows of the two objects.
     * @tparam T The type of the Matrix class
     * @param matrix1 the first matrix
     * @param matrix2 the second matrix
     * @return the minimum number of rows between the two matrices
     */
    template<typename T> unsigned long findMinRows(const Matrix<T> &matrix1, const Matrix<T> &matrix2)
    {
        if(matrix1.getRows() < matrix2.getRows())
        {
            return matrix1.getRows();
        }
        return matrix2.getRows();
    }
}
#endif // MATRIX
