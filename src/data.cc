/*
 *  data.cc
 *
 */ 

#include "data.h"
#include "sim.h"


std::streampos RunData::file_size(const std::string &filename)
{
	std::ifstream file( filename.c_str(), std::ios::binary | std::ios::ate);
	return file.tellg()/1024; 
}


void RunData::load(const std::string &filename)
{
	int section;
	int i,j,d;

	ublas::vector<double> vec;
	ublas::matrix<double> mat;

	// clean (maybe) used vectors
	B.clear();
	V.clear();
	W.clear();
	X.clear();

	// open input file
	std::ifstream file(filename.c_str());
	std::string line;
	
	while (std::getline(file, line))
	{
		// skip comments
		if (line[0] == '#')
			continue;
		
		// determine section by *
		if (line[0] == '*')
		{
			section = line[1];
			continue;
		}
		
		// prepare reading
		std::istringstream ss(line);

		switch (section)
		{
			// parameters
			case 'P':
				ss >> Ndays >> N >> Nplanets >> Tstep;
				break;
		
			// body
			case 'B':
			{
				Body bod;
				ss >> bod.name >> bod.color >> bod.is_probe;
				bod.init = ublas::zero_matrix<double> (5,3);
				for (j = 0; j < 5; j++)
					for (d = 0; d < 3; d++)
						ss >> bod.init(j,d);
				B.push_back(bod);
				break;
			}
			
			// velocity
			case 'V':
				mat = ublas::zero_matrix<double> (N,3);
				for (j = 0; j < N; j++)
					for (d = 0; d < 3; d++)
						ss >> mat(j,d);
				V.push_back(mat);
				break;
				
			// energy
			case 'W':
				vec = ublas::zero_vector<double> (3);	
				ss >> vec(0) >> vec(1) >> vec(2);
				W.push_back(vec);
				break;

			// positons
			case 'X':
				mat = ublas::zero_matrix<double> (N,3);
				for (j = 0; j < N; j++)
					for (d = 0; d < 3; d++)
						ss >> mat(j,d);
				X.push_back(mat);
				break;
		}

	}

	// close file
	file.close();

	// console output
	std::cout << "> File '" << filename << "' (" << file_size(filename)  << "k) loaded." << std::endl;
	std::cout << "> " << Nplanets << " planets and " << N-Nplanets << " probes found." << std::endl;
}


void RunData::save(const std::string &filename)
{
	int i,j,d;
	time_t now = time(0);

    // open output file and write header with current date and time
    std::ofstream file(filename.c_str());
    file << "# Simulation: " << ctime(&now);
	file << std::fixed << std::setprecision(10); // ALERT: set precision back to 10
	file << "# Units: P [1], X [au], W [au^2 em d^-2]" << std::endl;
    
	// parameters
	file << "*P\n";
	file << Ndays << "\t" << N << "\t" << Nplanets << "\t" << Tstep << "\n";

	// body
	file << "*B\n";
	for (i = 0; i < N; i++)
	{
		file << B[i].name << "\t" << B[i].color << "\t" << B[i].is_probe << "\t";
		for (j = 0; j < 5; j++)
			for (d = 0; d < 3; d++)
				file << B[i].init(j,d) << "\t";
		file << "\n";
	}
	
	// velocity
	file << "*V\n";
	for (i = 0; i < Ndays; i++)
	{
		for (j = 0; j < N; j++)
			file << V[i](j,0) << "\t" << V[i](j,1) << "\t" << V[i](j,2) << "\t";

		file << "\n";
	}

	// energy	
	file << "*W\n";
	for (i = 0; i < Ndays; i++)
		file << W[i](0) << "\t" << W[i](1) << "\t" << W[i](2) << "\n";	

	// positions	
	file << "*X\n";
    for (i = 0; i < Ndays; i++)
	{
        for (j = 0; j < N; j++)
			file << X[i](j,0) << "\t" << X[i](j,1) << "\t" << X[i](j,2) << "\t";

		file << "\n";
	}
 
	// close file
	file.close();

	// console output
	std::cout << "> File '" << filename << "' (" << file_size(filename)  << "k) written." << std::endl;
}

