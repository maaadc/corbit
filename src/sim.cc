/*
 *  sim.cc
 *
 */ 


#include "sim.h"


Simulation::Simulation()
/* Constructor */
{        
    // initialize variables
    data.N = 0;
	data.Nplanets = 0;
	ad_index = 1;
   	
	x_ij = ublas::vector<double> (3);
	 
    // set runtime variables
	M = ublas::vector<double> (1);
	X = ublas::matrix<double> (1,3);
	V = ublas::matrix<double> (1,3);
	A = ublas::matrix<double> (1,3);
	Aext = ublas::zero_matrix<double> (1,3);
	W = ublas::vector<double> (1,3);
}


void Simulation::acceleration()
/* calculate acceleration of each object based on Newton's law */
{
    int i,j,d;
	int j_max;
	
	/* Modes (highest accuracy and lowest speed at bottom):
	 * ad_index = ...
	 * 1				Each object only interacts with sun
	 * data.Nplanets	Each object only interacts with sun and planets (use this)
	 * data.N			Each object interacts with each other
	 */
	    
    // External acceleration:
	// Planets have zero, probes can have finite through engines 
	A = Aext;

    // Gravitational acceleration: For each object i ...
    for (i = 0; i < data.N; i++)
	{
		if (i >= data.Nplanets && ad_index == 1)
			j_max = 1;
		else
			j_max = i;

        // ... add up the accelerations towards objects j
        for (j = 0; j < j_max; j++)
        // calculate only the lower half up to the diagional
		// and neglect probe-probe interaction
		{
			// vector x_ij from object i to j
			for (d = 0; d < 3; d++)
				x_ij(d) = X(j,d) - X(i,d);
			// add G and inverse distance |x_ij|^-3
			x_ij *= G * pow( ublas::norm_2(x_ij) ,-3);
			// Newton's law: a_ij = - G m_j x_ij / |x_ij|^3
			for (d = 0; d < 3; d++)
			{
				// add to acceleration of objects i and j
				// by using x_ij = - x_ji
				A(i,d) += M(j) * x_ij(d);
				A(j,d) -= M(i) * x_ij(d);
			}
		}
	}
}


void Simulation::adaptive()
/* decide whether to use a faster, less accurate algorithm or the slower one */
{
	int i;
	double norm = 0;
	double maxacc = 0;
	const double threshold = 1.e-3;	
	ublas::matrix<double> Atemp;

	// get maximum gravitational acceleration of probes caused by planets:
	// get acc w/o planet-probe interaction
	ad_index = 1;
	acceleration();
	Atemp = A;

	// get acc w/ planet-probe interaction
	ad_index = data.Nplanets;
	acceleration();
		
	// calculate acc caused by planets over acc with sun and get maximum
	for (i = data.Nplanets; i < data.N; i++) 
	{
		norm = ublas::norm_2( ublas::row(A-Atemp, i) ); 
		norm /= ublas::norm_2( ublas::row(A-Aext, i) );
		if (norm > maxacc)
			maxacc = norm;
	}
	
	// compare maximum acceleration with threshold value
	std::cout << " AD: " << maxacc*100 <<  " % => ";
	if (maxacc > threshold)
	{
		ad_index = data.Nplanets;
		std::cout << "slow";
	}
	else
	{
		ad_index = 1;
		std::cout << "fast";
	}


}


void Simulation::add_probes(int Nprobes, double r0, double v0)
/* add some probes to the simulation */
{
	int i;
	double pi = 3.14159;

	std::cout << "> Probe start: r0 = " << r0 << " m,  v0 = " << v0 << " m/s." << std::endl;

	// r0 [m] => [au]  and v0 [m/s] => [au/d]
	r0 /= unit_length;
	v0 /= unit_length / unit_time;

	for (i = 0; i < Nprobes; i++)
	{
		Body bod;

		bod.name = "pr";
		bod.name += (char) (i+65);	// prA, prB, ...
		bod.color = "#aaaaaa";
		bod.is_probe = 1;
		
		// initialize init
		bod.init = ublas::zero_matrix<double> (5,3);

		// set mass to 750 kg
		bod.init(0,0) = 750. / unit_mass;

		// set fuel to 100 d
		bod.init(0,1) = 100.;

		// start probes in a circle around earth (index 3)
		// with radius r0 [au] and velocity v0 [au/d]
		bod.init(1,0) = data.B[3].init(1,0) + r0 * cos(2.*pi*i/Nprobes);
		bod.init(1,1) = data.B[3].init(1,1) + r0 * sin(2.*pi*i/Nprobes);
		bod.init(1,2) = data.B[3].init(1,2);

		bod.init(2,0) = data.B[3].init(2,0) + v0 * cos(2.*pi*i/Nprobes);
		bod.init(2,1) = data.B[3].init(2,1) + v0 * sin(2.*pi*i/Nprobes);
		bod.init(2,2) = data.B[3].init(2,2);
		
		// save
		data.B.push_back( bod );
			
		// update N
		data.N++;
	}

	std::cout << "> " << data.N-data.Nplanets << " probes added." << std::endl;
}


void Simulation::collision()
/* Checks if a probe collided into a planet */
{
	int i,j,d;
	double dist, coll_radius;
	double radius_factor = 2.;

	// fill matrix D(i, j) with distance between objects i and j, j < i
	D = ublas::zero_matrix<double> (data.N, data.N);
	for (i = 0; i < data.N; i++)
		for (j = 0; j < i; j++)
			D(i,j) = ublas::norm_2( ublas::row(X, j) - ublas::row(X, i) );
	
	// check for collisions of probes
	for (i = 0; i < data.Nplanets; i++)
		for (j = data.Nplanets; j < data.N; j++) 
		{
			// distance between planet i and probe j
			for (d = 0; d < 3; d++)
				x_ij(d) = X(j,d) - X(i,d);
			dist = ublas::norm_2(x_ij);

			// set collision radius as a multiple of planets radius
			coll_radius = radius_factor * data.B[i].init(0,1);
			
			// compare them and print a warning
			if (dist < coll_radius)
			{
				std::cout << "> Probe " << i << " collided with planet " << j << ": d = ";
				std::cout << dist << " au." << std::endl;
			}
		}

}


void Simulation::energy()
/* Return kinetic, potential and total energy */
{
    int i,j,d; 
	
	// set W to zero
	W = ublas::zero_vector<double> (3);
    
    // W_kin = 1/2 Sum_i m_i |v_i|^2
    for (i = 0; i < data.N; i++)
        for (d = 0; d < 3; d++)
			W(1) += M(i) * V(i,d) * V(i,d);
    W(1) *= .5;
    
    // W_pot = - G Sum_i,j m_i m_j / r_ij for i != j
    for (i = 0; i < data.N; i++)
        for (j = 0; j < i; j++)
        {
            // Vector x_ij from object i to j
            for (d = 0; d < 3; d++)
				x_ij(d) = X(j,d) - X(i,d);
			// add to potential energy
			W(2) += M(i) * M(j) / ublas::norm_2(x_ij);
        }
    W(2) *= -1. * G;
    
    // W_tot = W_kin + W_pot
	W(0) = W(1) + W(2);
}


void Simulation::integrate(int Nsteps)
/* Computes Nsteps integration steps */
{
    int i;
   
    for (i = 0; i < Nsteps; i++)
    {
        // Velocity Verlet algorithm: Force must not depend on velocity.
        
        // (1) x(t+dt) = x(t) + v(t) dt + 1/2 a(t) dt^2
		X += V * dt + A * dt * dt2;
		// (2) v(t+dt/2) = v(t) + 1/2 a(t) dt
		V += A * dt2;
        // (3) derive a(t+dt) from x(t+dt)
        acceleration();
        // (4) v(t+dt) = v(t+dt/2) + 1/2 a(t+dt) dt
		V += A * dt2;
    }
}


void Simulation::prepare() 
/* prepares a run by setting all needed variables */
{
	int i;

	// set time step in simulation units	
	dt = data.Tstep / unit_time;
	dt2 = 0.5 * dt;

	// initialize runtime variables to zero
	M = ublas::zero_vector<double> (data.N);
	X = ublas::zero_matrix<double> (data.N,3);
	V = ublas::zero_matrix<double> (data.N,3);
	A = ublas::zero_matrix<double> (data.N,3);
	Aext = ublas::zero_matrix<double> (data.N,3);
	W = ublas::zero_vector<double> (3);
	
	// fill with bodies 
	for (i = 0; i < data.N; i++)
	{
		M(i) = data.B[i].init(0,0);
		ublas::row(X, i) = ublas::row(data.B[i].init, 1);
		ublas::row(V, i) = ublas::row(data.B[i].init, 2);
		ublas::row(A, i) = ublas::row(data.B[i].init, 3);
		ublas::row(Aext, i) = ublas::row(data.B[i].init, 4);
	}

}


void Simulation::run(int Ndays, double Tstep)
{
    int i;

	data.Ndays = Ndays; // earth days to simulate 
	data.Tstep = Tstep;	// time step [s]

	// prepare run
	prepare();
	
   	// console output 
	std::cout << "> Starting simulation of " << data.N << " objects for " << std::setprecision(2);
	std::cout << std::fixed << (double) data.Ndays/365 << " year(s)." << std::endl;
	 
    // set acceleration for 1st time step
    acceleration();
	
	// start timer 
	boost::timer timer;

    for (i = 0; i < data.Ndays; i++)
    {
		// save data
		data.V.push_back( V );
		data.W.push_back( W );
		data.X.push_back( X );
		
		// adaptive decision
		adaptive();

        // integrate a whole day
        integrate( (int) unit_time/data.Tstep );
		
		// calculate energy 
		energy();

		// collision check
		collision();
			
		// console output
		std::cout << "\r> Running: " << (double) 100*(i+1)/data.Ndays << " % at ";
		std::cout << (double) i/timer.elapsed()  " d/s." << std::flush;
	}
    
    // close output
    std::cout << std::endl;
}

