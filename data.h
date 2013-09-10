/*
 *  data.h
 *
 *  Holds full data of a simulation run
 *
 */ 

#ifndef _DATA_H
#define _DATA_H

#include <ctime>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>

// boost numeric libraries
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
// boost property tree for saving and loading runs
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

namespace ublas = boost::numeric::ublas;


struct Body
{
	std::string name;
	std::string color;
	int is_probe;				// 0 = planet, 1 = probe
	ublas::matrix<double> init;	// [5,3] Attr, X, V, A, Aext
	
	/* Both:	Attr[0] = mass [em] 
	 * Planet:	Attr[1] = radius [au]
	 * Probe:	Attr[1] = fuel [d]
	 */
};


struct RunData
{
	std::vector< Body > B;					// [N] Body	
	std::vector< ublas::matrix<double> > V;	// [i,N,3] Velocity
	std::vector< ublas::vector<double> > W;	// [i,3] Wtot, Wkin, Wpot
	std::vector< ublas::matrix<double> > X; // [i,N,3] Position
	int Ndays, N, Nplanets;
	double Tstep;

	std::streampos file_size(const std::string &filename = "run.dat");	
	void load(const std::string &filename = "run.dat");
	void save(const std::string &filename = "run.dat");
};

#endif

