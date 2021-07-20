#include "basic_structs.h"
#include "xdr_wrapper.h"
#include "xdr_interface.h"
#include "parse_input.h"
#include "./marching_cubes/TriangulateGlobal.h"
extern "C" int Triangulate(float *vol_data, float **vertex_data,  int **triangle_data,
		int *n, float *s, int *nvtx, int *ntri, float thresh, int cmp);
int ri(int* i, int x, int y)
{
	return i[3*x+y];
}
float rt(float* r, int x, int y)
{
	return r[6*x+y];
}
void output_stl(int numTriangles, float* t, int* idx, std::string filename, double box_vec[3])
{
	std::ofstream ofile(filename + ".stl");
	ofile << "solid " << filename << "\n";
	for(int i = 0; i < numTriangles; i++)
	{
		ofile << "facet normal " << rt(t, ri(idx, i, 0),3) << "  " << rt(t, ri(idx, i, 0),4) << "  " << rt(t, ri(idx, i, 0),5) << "\n"; 
		ofile << "    outer loop\n";
		ofile << "vertex " << rt(t, ri(idx, i, 0),0)+0.5*box_vec[0] << " " << rt(t, ri(idx, i, 0),1)+0.5*box_vec[1] << " " << rt(t, ri(idx, i, 0),2)+0.5*box_vec[2] << "\n";
		ofile << "vertex " << rt(t, ri(idx, i, 1),0)+0.5*box_vec[0] << " " << rt(t, ri(idx, i, 1),1)+0.5*box_vec[1] << " " << rt(t, ri(idx, i, 1),2)+0.5*box_vec[2] << "\n";
		ofile << "vertex " << rt(t, ri(idx, i, 2),0)+0.5*box_vec[0] << " " << rt(t, ri(idx, i, 2),1)+0.5*box_vec[1] << " " << rt(t, ri(idx, i, 2),2)+0.5*box_vec[2] << "\n";
		ofile << "    endloop\n";
		ofile << "endfacet\n";
	}
	ofile << "endsolid " << filename << std::endl;
	ofile.close();
}


void write_float_binary(std::ofstream& out, float f)
{
   out.write( reinterpret_cast<const char*>( &f ), sizeof( float ) );
}
void write_uint32_binary(std::ofstream& out, uint32_t f)
{
   out.write( reinterpret_cast<const char*>( &f ), sizeof( uint32_t ) );
}
void write_char_array_binary(std::ofstream& out, char* f, int size)
{
	out.write(f, size);
}
void output_stl_binary(int numTriangles, float* t, int* idx, std::string ofilename, double box_vec[3])
{
	std::ofstream ofile(ofilename + ".stl", std::ios::binary);
	char header[80] = {0};
	char term[2] = {0,0};
	uint32_t num_faces = numTriangles;
	write_char_array_binary(ofile, header, 80);
	write_uint32_binary(ofile, num_faces);
	for(unsigned int i = 0; i < numTriangles; i++)
	{
		write_float_binary(ofile, rt(t, ri(idx, i, 0),3));
		write_float_binary(ofile, rt(t, ri(idx, i, 0),4));
		write_float_binary(ofile, rt(t, ri(idx, i, 0),5));
		
		write_float_binary(ofile, rt(t, ri(idx, i, 0),0) + 0.5*box_vec[0]);
		write_float_binary(ofile, rt(t, ri(idx, i, 0),1) + 0.5*box_vec[1]);
		write_float_binary(ofile, rt(t, ri(idx, i, 0),2) + 0.5*box_vec[2]);
		
		write_float_binary(ofile, rt(t, ri(idx, i, 1),0) + 0.5*box_vec[0]);
		write_float_binary(ofile, rt(t, ri(idx, i, 1),1) + 0.5*box_vec[1]);
		write_float_binary(ofile, rt(t, ri(idx, i, 1),2) + 0.5*box_vec[2]);	
		
		write_float_binary(ofile, rt(t, ri(idx, i, 2),0) + 0.5*box_vec[0]);
		write_float_binary(ofile, rt(t, ri(idx, i, 2),1) + 0.5*box_vec[1]);
		write_float_binary(ofile, rt(t, ri(idx, i, 2),2) + 0.5*box_vec[2]);		
		
		write_char_array_binary(ofile, term, 2);
	}
	ofile.close();
}

void marching_cubes_driver(input_file in, voxel_grid& v, std::string out_mod)
{
	double gs = v.get_gs();
	int nx = v.getsize(0);
	int ny = v.getsize(1);
	int nz = v.getsize(2);
	double box_vec[3] = {gs*(nx), gs*(ny), gs*(nz)};
	std::vector<float> volume_data(nx*ny*nz, 0);
	int iterator = 0;
	for(int k = 0; k < nz; k++)
	for(int j = 0; j < ny; j++)
	for(int i = 0; i < nx; i++)
	{
		volume_data[iterator] = v.get_grid_val(i, j, k);;
		iterator++;
	}
	float* vertex_data;
	int* triangle_data;
	int nvtx;
	int ntri;
	int n[3] = {nx, ny, nz};
	float s[3] = {gs, gs, gs};
	float thresh = in.isovalue; //isovalue
	Triangulate(&volume_data[0], &vertex_data, &triangle_data,
		n, s, &nvtx, &ntri, thresh, 0);	
	if(in.outputfile == "stl_binary")
	output_stl_binary(ntri, vertex_data, triangle_data, in.outputfile + out_mod, box_vec);	
	else //defaults to stl_ascii
	output_stl(ntri, vertex_data, triangle_data, in.outputfile + out_mod, box_vec);
	return;
}

void create_voxel_grid(input_file& in)
{

	read_input_file(in);
	file_wrapper traj_file(in.trajectoryfile, "r");
	int natoms = traj_file.get_natoms();
	std::vector<std::vector<int> > index_lists;
	for(int i = 0; i < in.ngroups; i++)
	{
		std::vector<int> index_list = get_index_list(in.indexfile, in.groupname[i]);
		index_lists.push_back(index_list);
	}
	if(in.mode == "average")
	{
		voxel_grid v(in.gridpoints[0], in.gridpoints[1], in.gridpoints[2], 
					in.gridspacing, in.sigma, in.density);
		double t = 0;
		int skip_iterator = 0;
		while(traj_file.next_frame())
		{
			t = traj_file.get_t();
			if(t < in.tmin || t > in.tmax) continue;
			skip_iterator++;
			if(skip_iterator%in.skip != 1) continue;
			else v.it_counter();
			for(int k = 0; k < in.ngroups; k++)
			for(int i = 0; i < (int)index_lists[k].size(); i++)
			{
				v.add_gaussian(
							traj_file.get_x(index_lists[k][i], 0),
							traj_file.get_x(index_lists[k][i], 1),
							traj_file.get_x(index_lists[k][i], 2), k);
			}
		}
		marching_cubes_driver(in, v, "_avg");
	}
	else if(in.mode == "sequence")
	{
		double t = 0;
		int iterator = 0;
		int skip_iterator = 0;
		while(traj_file.next_frame())
		{
		voxel_grid v(in.gridpoints[0], in.gridpoints[1], in.gridpoints[2], 
					in.gridspacing, in.sigma, in.density);
			t = traj_file.get_t();
			if(t < in.tmin || t > in.tmax) continue;
			skip_iterator++;
			if(skip_iterator%in.skip != 1) continue;
			v.it_counter();
			for(int k = 0; k < in.ngroups; k++)
			for(int i = 0; i < (int)index_lists[k].size(); i++)
			{
				v.add_gaussian(
							traj_file.get_x(index_lists[k][i], 0),
							traj_file.get_x(index_lists[k][i], 1),
							traj_file.get_x(index_lists[k][i], 2), k);
			}
			marching_cubes_driver(in, v, "_"+std::to_string(iterator));
			iterator++;
		}
	}
	
	return;
}


int main(int argc, char **argv)
{
	if(argc == 1)
	{
		std::cout << "Not enough input arguments... please specify an input file" << std::endl;
		return 0;
	}
	input_file in;
	std::string infile_name = argv[1];
	std::ifstream test_file(infile_name);
	if(!test_file.is_open())
	{
		std::cout << "Could not open file... closing. \n";
		return 0;		
	}
	
	in.inputfile = infile_name;
	create_voxel_grid(in);
	return 0;
}
