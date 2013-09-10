/*
 *  corbit.cc
 *
 */ 

#include "sim.h"


int main()
{
    // set up system
    Simulation sim;
    
    // add planets as of 2012-01-01 00:00 UTC
    // units: length [au], time [d], mass [em]
    // data from http://ssd.jpl.nasa.gov/horizons.cgi
	sim.data.load("planets.dat");
	
	// add probes starting from earth
//	sim.add_probes(25, 1.e9, 1.e4);
	
    // run simulation
    sim.run(365, 10.);

	// save data
	sim.data.save();

    return 0;
}

