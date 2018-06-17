//This file is part of NECSim project which is released under MIT license.
//See file **LICENSE.txt** or visit https://opensource.org/licenses/MIT) for full license details.

/**
 * @author James Rosindell
 * @file NRrand.h
 * @brief Contains a generic random number generator.
 * Provided by James Rosindell (j.rosindell@imperial.ac.uk) with moderate modifications by
 * Samuel Thompson (thomsonsed@gmail.com).
 * 
 * The definitions for the constants defined here should not be altered.
 * @copyright <a href="https://opensource.org/licenses/MIT"> MIT Licence.</a>
 */
#ifndef FATTAIL_H
#define FATTAIL_H

//
//#define IM1 2147483563
//#define IM2 2147483399
//#define AM (1.0/IM1)
//#define IMM1 (IM1-1)
//#define IA1 40014
//#define IA2 40692
//#define IQ1 53668
//#define IQ2 5277
//#define IR1 12211
//#define IR2 3791
//#define NTAB 32
//#define NDIV (1+IMM1/NTAB)
//#define EPS 1.2e-8
//#define RNMX (1.0-EPS)

#include <cstdio>
#include <string>
#include <iomanip>

#define _USE_MATH_DEFINES

#include <cmath>
#include <algorithm>
#include <vector>
#include <iostream>
#include <fstream>
#include <climits>
#include "Logger.h"

using namespace std;
/**
 * @brief These variables contain the special numbers for random number generation.
 */
const long IM1 = 2147483563;
const long IM2 = 2147483399;
const double AM = (1.0 / IM1);
const long IMM1 = (IM1 - 1);
const long IA1 = 40014;
const long IA2 = 40692;
const long IQ1 = 53668;
const long IQ2 = 5277;
const long IR1 = 12211;
const long IR2 = 3791;
const long NTAB = 32;
const double NDIV = (1 + IMM1 / NTAB);
const double EPS = 1.2e-8;
const double RNMX(1.0 - EPS);

/**
 * @brief Contains the functions for random number generation.
 */
class NRrand {

private:
	long idum{};
	int j{};
	long k{};
	long idum2{};
	long iy{};
	long iv[NTAB]{};
	double temp{};
	bool seeded;

	double lastresult{};
	bool normflag;
	// for the L value of the dispersal kernel (the width - does not affect the shape).
	double tau;
	// for the sigma value of the dispersal kernel (the variance of a normal distribution).
	double sigma;

	typedef double (NRrand::*fptr)(); // once setup will contain the dispersal function to use for this simulation.
	fptr dispersalFunction;
	// the probability that dispersal comes from the uniform distribution. This is only relevant for uniform dispersals.
	double m_prob{};
	// the cutoff for the uniform dispersal function i.e. the maximum value to be drawn from the uniform distribution.
	double cutoff{};
public:

	/**
	 * @brief Standard constructor.
	 */
	NRrand() {
		seeded = false;
		normflag = true;
		dispersalFunction = nullptr;
		sigma = 0;
		tau = 0;
	}

	/**
	 * @brief Sets the seed to the given input.
	 * Is only seeded if the seed hasn't already been provided.
	 * @param seed the input seed.
	 */
	void setSeed(long seed) {
		if (!seeded) {
			idum2 = 123456789;
			iy = 0;
			idum = seed;
			if (idum < 1) idum = 1;
			//Be sure to prevent idum = 0.
			idum2 = (idum);
			for (j = NTAB + 7; j >= 0; j--) {
				//Load the shuffle table (after 8 warm-ups).
				k = (idum) / IQ1;
				idum = IA1 * (idum - k * IQ1) - k * IR1;
				if (idum < 0) idum += IM1;
				if (j < NTAB) iv[j] = idum;
			}
			iy = iv[0];
			seeded = true;
		} else {
			throw runtime_error("Trying to set the seed again: this can only be set once.");
		}
	}

	/**
	 * @brief Clears the seed, if it has already been set.
	 * Keeps other simulation parameters, such as sigma and tau.
	 */
	void wipeSeed() {
		seeded = false;
	}

	/**
	 * @brief The random number generator.
	 * Uses Schrage's method and a shuffle table to generate the output.
	 * @return the random number (a double between 0 and 1).
	 */
	double d01() {
		k = (idum) / IQ1;
		//Start here when not initializing.
		idum = IA1 * (idum - k * IQ1) - k * IR1;
		//Compute idum=(IA1*idum) % IM1 without overflows by Schrage's method. 
		if (idum < 0) idum += IM1;
		k = idum2 / IQ2;
		idum2 = IA2 * (idum2 - k * IQ2) - k * IR2;
		//Compute idum2=(IA2*idum) % IM2 likewise.
		if (idum2 < 0) idum2 += IM2;
		j = iy / NDIV;
		//Will be in the range 0..NTAB-1.
		iy = iv[j] - idum2;
		//Here idum is shuffled, idum and idum2 are combined to generate output. 
		iv[j] = idum;
		if (iy < 1) iy += IMM1;
		if ((temp = AM * iy) > RNMX) {
			//os << "random call = " << "RNMAX" << "\n";
			return RNMX; //Because users don't expect endpoint values.
		}
		return temp;

	}

	/**
	 * @brief Generates a random number uniformly from 0 to the maximum value provided.
	 * @param max the maximum number.
	 * @return an integer of the produced random number.
	 */
	unsigned long i0(unsigned long max) {
		return (unsigned long) (d01() * (max + 1));
	}

	/**
	 * @brief Generates a normally distributed number
	 * Uses the standard normal distribution from a Box-Muller transform.
	 * @return the random number from a normal distribution.
	 */
	double norm() {
		if (normflag) {
			double r2 = 2;
			double xx = 0;
			double yy = 0;
			while (r2 > 1) {
				xx = 2.0 * d01() - 1.0;
				yy = 2.0 * d01() - 1.0;
				r2 = (xx * xx) + (yy * yy);
			}
			double fac = sqrt(-2.0 * log(r2) / r2);
			lastresult = xx * fac;
			double result = yy * fac;
			normflag = false;
			return sigma * result;
		} else {
			normflag = true;
			return sigma * lastresult;
		}
	}

	/**
	 * @brief Returns a 2 dimensional call from a normal distribution, giving a distance in cartesian space
	 * This way is slightly inefficient for normal distributions, but it's kept this way to make the 
	 * Map::runDispersal() function applicable with any dispersal type.
	 * @return dispersal distance of a normal distribution
	 */
	double norm2D() {
		double distx, disty;
		distx = norm();
		disty = norm();
		return pow(pow(distx, 2) + pow(disty, 2), 0.5);
	}

	/**
	 * @brief Sets the dispersal parameters, avoiding requirement to provide these numbers each function call.
	 * This is only relevant for fat-tailed dispersal calls.
	 * @param sigmain the fatness of the fat-tailed dispersal kernel.
	 * @param tauin the width of the fat-tailed dispersal kernel.
	 */
	void setDispersalParams(const double sigmain, const double tauin) {
		sigma = sigmain;
		tau = tauin; // used to invert the sign here, doesn't any more.
	}


	/**
	 * @brief Call from the fat-tailed dispersal kernel with the provided sigma.
	 * @deprecated This is the original version used in J Rosindell's codebase, and has been altered for
	 * a version which approximates the gaussian distribution at extreme limits.
	 * @param z the desired sigma.
	 * @return a random number drawn from the fat-tailed dispersal kernel.
	 */
	double fattail(double z) {
		double result;
		result = pow((pow(d01(), (1.0 / (1.0 - z))) - 1.0), 0.5);
		return result;
	}

	// this new version corrects the 1.0 to 2.0 and doesn't require the values to be passed every time.
	/**
	 * @brief Call from fat-tailed dispersal kernel.
	 * This function requires setDispersalParams() has already been called.
	 * @deprecated deprecated, kept for testing purposes only
	 * @return a random number drawn from the fat-tailed dispersal kernel.
	 */
	double fattail() {
		double result;
		// old function version (kept for reference)
//		result = (tau * pow((pow(d01(),(2.0/(2.0-sigma)))-1.0),0.5));
		result = (sigma * pow((tau * (pow(d01(), -2.0 / tau)) - 1.0), 0.5));
		return result;
	}

	/**
	 * @brief Old version of the function call reparameterised for different nu and sigma.
	 * @deprecated Kept only for testing purposes.
	 * @return a random number drawn from the fat-tailed dispersal kernel.
	 */
	double fattail_old() {
		double result;
		result = (sigma * pow((pow(d01(), (2.0 / (2.0 + tau))) - 1.0), 0.5));
		return result;
	}

	/**
	 * @brief Generates a direction in radians.
	 * @return the direction in radians
	 */
	double direction() {
		return (d01() * 2 * M_PI);
	}

	/**
	 * @brief For a given event probability, returns the probability that the event has occured.
	 * @param event_probability the event probability.
	 * @return whether or not the event has occured.
	 */
	bool event(double event_probability) {
		if (event_probability < 0.000001) {
			if (d01() <= 0.000001) {
				return (event(event_probability * 1000000.0));
			}
			return false;
		}
		if (event_probability > 0.999999) {
			return (!(event(1.0 - event_probability)));
		}
		return (d01() <= event_probability);


	}

	/**
	 * @brief Normal distribution, with percentage chance to choose a uniform distribution instead. 
	 * @note This function will not produce the same output as norm() for the same parameters, even with a
	 * zero chance of picking from the uniform distribution (due to random number draws).
	 * @return normally (or uniformly) distributed number
	 */
	double normUniform() {
		// Check if the dispersal event comes from the uniform distribution
		if (d01() < m_prob) {
			// Then it does come from the uniform distribution
			return (d01() * cutoff);
		}
		return norm2D();
	}


	/**
	 * @brief Two uniform distributions, the first between 0 and 0.1*cutoff, and the second between 0.9*cutoff and
	 * cutoff. Selects from both distributions equally.
	 * @note The mean for this function should be identical to a uniform distribution between 0 and cutoff.
	 * @return uniformly distributed number
	 */
	double uniformUniform() {
		if (d01() < 0.5) {
			// Then value comes from the first uniform distribution
			return (d01() * cutoff * 0.1);
		}
		// Then the value comes from the second uniform distribution
		return 0.9 * cutoff + (d01() * cutoff * 0.1);
	}

	/**
	 * @brief Sets the dispersal method by creating the link between dispersalFunction() and the correct
	 * dispersal character
	 * @param dispersal_method string containing the dispersal type. Can be one of [normal, fat-tail, norm-uniform]
	 * @param m_probin the probability of drawing from the uniform distribution. Only relevant for uniform dispersals.
	 * @param cutoffin the maximum value to be drawn from the uniform dispersal. Only relevant for uniform dispersals.
	 */
	void setDispersalMethod(const string &dispersal_method, const double &m_probin, const double &cutoffin) {
		if (dispersal_method == "normal") {
			dispersalFunction = &NRrand::norm2D;
			if (sigma < 0) {
				throw invalid_argument("Cannot have negative sigma with normal dispersal");
			}
		} else if (dispersal_method == "fat-tail" || dispersal_method == "fat-tailed") {
			dispersalFunction = &NRrand::fattail;
			if (tau < 0 || sigma < 0) {
				throw invalid_argument("Cannot have negative sigma or tau with fat-tailed dispersal");
			}
		} else if (dispersal_method == "norm-uniform") {
			dispersalFunction = &NRrand::normUniform;
			if (sigma < 0) {
				throw invalid_argument("Cannot have negative sigma with normal dispersal");
			}
		} else if (dispersal_method == "uniform-uniform") {
			// This is just here for testing purposes
			dispersalFunction = &NRrand::uniformUniform;
		}
			// Also provided the old version of the fat-tailed dispersal kernel
		else if (dispersal_method == "fat-tail-old") {
			dispersalFunction = &NRrand::fattail_old;
			if (tau > -2 || sigma < 0) {
				throw invalid_argument(
						"Cannot have sigma < 0 or tau > -2 with fat-tailed dispersal (old implementation).");
			}
		} else {
			throw runtime_error("Dispersal method not detected. Check implementation exists");
		}
		m_prob = m_probin;
		cutoff = cutoffin;
	}


	/**
	 * @brief Runs the dispersal with the allocated dispersal function.
	 *
	 * @note This function will never return a value larger than the size of LONG_MAX to avoid issues of converting
	 * doubles to integers. For dispersal distance within coalescence simulations, this is seemed a reasonable
	 * assumption, but may cause issues if code is re-used in later projects.
	 *
	 * @return distance the dispersal distance
	 */
	double dispersal() {
		return min(double(LONG_MAX), (this->*dispersalFunction)());
	}

	// to reconstruct distribution, use x = fattail/squrt(1+direction) , y = fattail/squrt(1+(direction^-1))

	/**
	 * @brief Outputs the NRrand object to the output stream.
	 * Used for saving the object to file.
	 * @param os the output stream.
	 * @param r the NRrand object to output.
	 * @return the output stream.
	 */
	friend ostream &operator<<(ostream &os, const NRrand &r) {
		//os << m.numRows<<" , "<<m.numCols<<" , "<<endl; 
		os << setprecision(64);
		os << r.idum << ",";
		os << r.j << ",";
		os << r.k << ",";
		os << r.idum2 << ",";
		os << r.iy << ",";
		for (long i : r.iv) {
			os << i << ",";
		}
		os << r.temp << ",";
		os << r.seeded << ",";
		os << r.lastresult << ",";
		os << r.normflag << "," << r.tau << "," << r.sigma << "," << r.m_prob << "," << r.cutoff;
		return os;
	}

	/**
	 * @brief Inputs the NRrand object from the input stream.
	 * Used for reading the NRrand object from a file.
	 * @param is the input stream.
	 * @param r the NRrand object to input to.
	 * @return the input stream.
	 */
	friend istream &operator>>(istream &is, NRrand &r) {
//		os << "starting NR read" << endl;
		char delim;
		//double temp1,temp2;
		//is << m.numRows<<" , "<<m.numCols<<" , "<<endl; 
		is >> r.idum;
//		os << r.idum << endl;
//		string tmp;
//		is >> delim >> tmp;
//		os << tmp << endl;
//		os << delim;
		is >> delim;
		is >> r.j;
		is >> delim;
		is >> r.k;
		is >> delim;
		is >> r.idum2;
		is >> delim;
		is >> r.iy;
		is >> delim;
		for (long &i : r.iv) {
			is >> i;
			is >> delim;
		}
		is >> r.temp;
		is >> delim;
		is >> r.seeded;
		is >> delim;
		is >> r.lastresult;
		is >> delim;
		is >> r.normflag;
		is >> delim >> r.tau >> delim >> r.sigma >> delim >> r.m_prob >> delim >> r.cutoff;
		return is;
	}
};

#endif