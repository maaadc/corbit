corbit
======
The aim of this hobby simulation is to find a way for a Planetary Grand Tour, see http://en.wikipedia.org/wiki/Planetary_Grand_Tour

Overview
--------
Simulates an N-body system with gravitational interaction using the Velocity Verlet algorithm. The sun, planets and probes (e. g. Satellites) can be declared via a text file. The interaction between planets and probes is neglected, if a probe is far away from a planet.<br />
The output can be visualized in 3D along with some other plots.


Dependencies
------------
Compiles and runs fine on Debian Wheezy with the following dependencies:<br />
Simulation (corbit): g++ 4.7.2, boost 1.49<br />
Plotting (plot.py): python 2.7.3, numpy 1.6.2, matplotlib 1.1.1

