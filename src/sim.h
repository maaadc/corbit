/*
 *  sim.h
 *
 *  Simulation of an N-body system under the influence of gravitational force
 *
 */ 

#ifndef _SIM_H
#define _SIM_H

#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>

// boost numeric libraries
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
// boost property tree for saving and loading runs
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
// boost timer 
#include <boost/timer.hpp>

#include "data.h"

namespace ublas = boost::numeric::ublas;


class Simulation
{

public:
	/* variables */
	RunData data;	// stores full data of the run

	/* methods */
	Simulation();
	void add_probes(int Nprobes, double r0, double v0);
	void read_planets(std::string file_name);
	void run(int Ndays, double Tstep);

private:
    /* physics constants */
    static const double unit_length = 149.60e9;	// [m/au] mean distance earth-sun
    static const double unit_time = 86400.;     // [s/d] earth day
    static const double unit_mass = 5.9736e24;  // [kg/em] earth mass 
    static const double G = 8.8897235e-10;  	// [au^3 em^-1 d^-2] gravity const
    
	static const double ad_threshold = 2.e-8;	// threshold acceleration of a probe

    /* variables */
    double dt, dt2;      			// full and half time step
	int ad_index;					// adaptive index
    ublas::vector<double> M;		// [N] mass
	ublas::matrix<double> X;		// [N,3] position
	ublas::matrix<double> V;		// [N,3] velocity
	ublas::matrix<double> A;		// [N,3] acceleration
	ublas::matrix<double> Aext;		// [N,3] external acceleration
	ublas::vector<double> W;		// [3] Wtot, Wkin, Wpot
	ublas::matrix<double> D;		// [N,N] distance
   
   	/* temporary variables, used by few functions */
	double temp;
	ublas::vector<double> x_ij;
    
    /* methods */
    void acceleration();
	void adaptive();
	void collision();
    void energy();
    void integrate(int Nsteps);
	void prepare();
};


#endif

