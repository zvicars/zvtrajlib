#pragma once
#include "basic_structs.h"
struct input_file
{
	std::string inputfile;
	std::string trajectoryfile; //xtc compressed trajectory
	std::string outputfile;
	std::string outputtype = "stl_ascii"; //file type for output mesh, currently stl_binary and stl_ascii are valid
	std::string indexfile; //.ndx file from gromacs
	std::string mode; //average or sequence
	std::vector<std::string> groupname; //group entries in index files
	std::vector<double> sigma; //coarse-graining parameter for each group
	double tmin;
	double tmax;
	std::vector<double> density; //number density/nm
	double boxsize[3]; //absolute size of box
	double gridspacing;
	double gridpoints[3];
	double isovalue = 0.5;
	int ngroups; //number of idx file groups
	int skip; //gap between processed frames
};
void read_input_file(input_file& in)
{
	std::ifstream ifile(in.inputfile);
	std::string line;
	int total_args = 0;
	while(getline(ifile, line))
	{
		if(line.find("output file:") != std::string::npos)
		{
			in.outputfile = trim(line.substr(line.find(':') + 1));
			total_args++;
			continue;
		}
		if(line.find("trajectory file:") != std::string::npos)
		{
			in.trajectoryfile = trim(line.substr(line.find(':') + 1));
			total_args++;
			continue;
		}
		if(line.find("index file:") != std::string::npos)
		{
			in.indexfile = trim(line.substr(line.find(':') + 1));
			total_args++;
			continue;
		}
		if(line.find("group name:") != std::string::npos)
		{
			std::stringstream ss(line.substr(line.find(':') + 1));
			while(ss >> line)
			{
				in.groupname.push_back(line);
			}
			total_args++;
			continue;
		}
		if(line.find("sigma:") != std::string::npos)
		{
			std::stringstream ss(line.substr(line.find(':') + 1));
			while(ss >> line)
			{
				in.sigma.push_back(std::stod(line));
			}
			total_args++;
			continue;
		}
		if(line.find("density:") != std::string::npos)
		{
			std::stringstream ss(line.substr(line.find(':') + 1));
			while(ss >> line)
			{
				in.density.push_back(std::stod(line));
			}
			total_args++;
			continue;
		}
		if(line.find("mode:") != std::string::npos)
		{
			in.mode = trim(line.substr(line.find(':') + 1));
			total_args++;
			continue;		
		}
		if(line.find("tmin:") != std::string::npos)
		{
			in.tmin = std::stod(line.substr(line.find(':') + 1));
			total_args++;
			continue;
		}
		if(line.find("tmax:") != std::string::npos)
		{
			in.tmax = std::stod(line.substr(line.find(':') + 1));
			total_args++;
			continue;
		}
		if(line.find("skip:") != std::string::npos)
		{
			in.skip = std::stoi(line.substr(line.find(':') + 1));
			total_args++;
			continue;
		}
		if(line.find("box size:") != std::string::npos)
		{
			std::stringstream ss(line.substr(line.find(':') + 1));
			ss >> in.boxsize[0] >> in.boxsize[1] >> in.boxsize[2];
			total_args++;
			continue;
		}
		if(line.find("grid spacing:") != std::string::npos)
		{
			in.gridspacing = std::stod(line.substr(line.find(':') + 1));
			total_args++;
			continue;
		}
		if(line.find("isovalue:") != std::string::npos)//optional default 0.5
		{
			in.isovalue = std::stod(line.substr(line.find(':') + 1));
			continue;
		}
		if(line.find("output type:") != std::string::npos) //optional, default stl ascii which can be opened in vmd
		{
			in.outputtype = trim(line.substr(line.find(':') + 1));
			continue;
		}
	}
	if(total_args != 12 )
	{
		std::cout << "Not all arguments specified in input file... closing" << std::endl;
		throw;
	}
	if(in.density.size() != in.groupname.size() || 
		in.groupname.size() != in.sigma.size())
	{
		std::cout << "Size mismatch between vector parameters.";
		throw;
	}
	in.gridpoints[0] = round(in.boxsize[0]/in.gridspacing);
	in.gridpoints[1] = round(in.boxsize[1]/in.gridspacing);
	in.gridpoints[2] = round(in.boxsize[2]/in.gridspacing);
	in.ngroups = in.groupname.size();
	ifile.close();
	return;
}