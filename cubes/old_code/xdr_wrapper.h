#pragma once
#include "xdrfile_xtc.h"
#include "basic_structs.h"
//XDRFILE* xdrfile_open(const char* path, const char* mode);
//read_trr(XDRFILE *xd,int natoms,int *step,float *t,float *lambda, matrix box,rvec *x,rvec *v,rvec *f);
const int DIM_  = 3;
using Rvec   = float[DIM_];

class file_wrapper
{
	XDRFILE* file_handle = 0;
	int natoms = 0;
	int step = 0;
	float t = 0;
	float prec;
	matrix box;
	Rvec* x = 0;
public:
	file_wrapper(std::string filename, std::string mode = "r")
	{
		file_handle = xdrfile_open(filename.c_str(), mode.c_str());
		read_xtc_natoms( (char *)filename.c_str(), &natoms);
		x = new Rvec[natoms];
	}
	~file_wrapper()
	{
		xdrfile_close(file_handle);
		delete x;

	}
	bool next_frame()
	{
		std::cout << "Time = " << t << std::endl;
		return exdrOK == read_xtc(file_handle,natoms,&step,&t,box,x,&prec);
	}
	float get_x(int a_idx, int dim_idx)
	{
		return x[a_idx][dim_idx];
	}
	float get_t()
	{
		return t; 
	}
	float get_natoms()
	{
		return natoms;
	}
	void get_box(int (&x_out)[3])
	{
		x_out[0] = box[0][0];
		x_out[1] = box[1][1];
		x_out[2] = box[2][2];
		return;
	}
	
};

