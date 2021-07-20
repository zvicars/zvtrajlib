#pragma once
#include "basic_structs.h"
#include "xdr_wrapper.h"

std::vector<int> get_index_list(std::string index_file_name, std::string atom_name)
{
	std::vector<int> output;
	std::ifstream ifile(index_file_name);
	std::string line;
	bool read_flag = 0;
	while(getline(ifile, line))
	{
		if(line.find("[") != std::string::npos)
		{
			if(read_flag == 0)
			{
				line = line.substr(line.find("[") + 1, line.find("]") - line.find("[") - 1);
				line = trim(line);
				bool result = line == atom_name;
				if(line == atom_name)
				{
					read_flag = 1;
					
				}
				continue;
			}
			else break;
		}
		std::stringstream ss(line);
		int temp_var;
		while(ss >> temp_var && read_flag == 1)
		{
			output.push_back(temp_var);
		}
		
	}
	return output;
}