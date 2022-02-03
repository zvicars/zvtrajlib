#include "actions.hpp"
#include "../tools/pbcfunctions.hpp"
#include "../tools/cellgrid.hpp"
#include "../tools/stlmath.hpp"
#include "../tools/StringTools.hpp"
#include "Eigen/Eigen"
#include <set>
struct UnitCell
{
	double a, b, c, alpha, beta, gamma;
	std::string resname;
	std::vector<Atom> atoms;
	UnitCell(){return;}
	UnitCell(double ai, double bi, double ci, double alphai, double betai, double gammai, std::string rname, const std::vector<Atom>& atoms_in){
		a = ai; b=bi; c=ci; alpha=alphai; beta=betai; gamma=gammai; resname = rname;
		int counter = 1;
		for(auto atom : atoms_in){
			atom.index = counter;
			add_atom(atom);
			counter++;
		}
		return;		
	}
	void add_atom(Atom atom)
	{
		atom.resnr = 1;
		atom.resname = resname;
		atoms.push_back(atom);
		return;
	}
};
std::vector<Atom> get_shifted_relative(double x, double y, double z, const UnitCell& uc, const Eigen::Matrix3d& m_it )
{
	std::vector<Atom> atoms;
	for(std::size_t i = 0; i < uc.atoms.size(); i++)
	{
		Atom a1 = uc.atoms[i];
		Eigen::Vector3d v, xv;
		v << a1.x[0], a1.x[1], a1.x[2];
		xv = m_it*v;
		a1.x[0] = xv[0] + x;
		a1.x[1] = xv[1] + y;
		a1.x[2] = xv[2] + z;
		atoms.push_back(a1);
	}
	return atoms;
}
UnitCell load_crystal(std::string filename)
{
	double ca, cb, cc, calpha, cbeta, cgamma;
	std::string rname = "CRY";
	bool ba=0, bb=0, bc=0, balpha=0, bbeta=0, bgamma=0;
	std::vector<Atom> atoms;
	
	std::ifstream ifile;
	std::string line;
	ifile.open(filename);
	FANCY_ASSERT(ifile.is_open(), "Failed to open file in boxtools::actions::load_crystal");
	while(std::getline(ifile, line))
	{
		line = StringTools::trimWhitespace(line);
		if(line.at(0) == '_') //unit cell parameter is being specified
		{
			std::stringstream ss(line);
			std::string param_name;
			double param_value;
			if(!(ss >> param_name)){
				throw std::string("Failed to load param name");
			}
			if(!(ss >> param_value)){
				throw std::string("Failed to load parameter value for " + param_name);
			}	
			if(param_name == "_cell_length_a"){
				ca = 0.1*param_value; ba = 1; continue;
			}
			if(param_name == "_cell_length_b"){
				cb = 0.1*param_value; bb = 1; continue;
			}
			if(param_name == "_cell_length_c"){
				cc = 0.1*param_value; bc = 1; continue;
			}
			if(param_name == "_cell_angle_alpha"){
				calpha = param_value; balpha = 1; continue;
			}
			if(param_name == "_cell_angle_beta"){
				cbeta = param_value; bbeta = 1; continue;
			}
			if(param_name == "_cell_angle_gamma"){
				cgamma = param_value; bgamma = 1; continue;
			}			
		}
		else if(line.substr(0,7) == "resname"){
			rname = StringTools::trimWhitespace(line.substr(8));
		}
		else if(line == "BEGIN LIST")
		{
			while(1)
			{
				std::getline(ifile, line);
				line = StringTools::trimWhitespace(line);
				if(line == "END LIST") break;
				std::string name, atype;
				double x_rel, y_rel, z_rel;
				std::stringstream ss(line);
				if(!(ss >> name)) throw std::string("Failed to load atom name");
				if(!(ss >> atype)) throw std::string("Failed to read atom type");
				if(!(ss >> x_rel)) throw std::string("Failed to read x value for " + name + " and " + atype);
				if(!(ss >> y_rel)) throw std::string("Failed to read y value" + name + " and " + atype);
				if(!(ss >> z_rel)) throw std::string("Failed to read z value" + name + " and " + atype);
				Atom a1;
				a1.name = name;
				a1.type = atype;
				a1.x[0] = x_rel;
				a1.x[1] = y_rel;
				a1.x[2] = z_rel;
				a1.v.fill(0.0);
				atoms.push_back(a1);
			}
		}
	}
	FANCY_ASSERT(ba == 1 && bb == 1 && bc == 1 && balpha == 1 && bbeta == 1 && bgamma == 1, "Not all parameters specified in .zcif file.");
	return UnitCell(ca, cb, cc, calpha, cbeta, cgamma, rname, atoms);
}
std::vector<Atom> replicate_cells(int nx, int ny, int nz, const UnitCell& uc, std::vector<double>& box_dims)
{
	std::vector<Atom> atoms;
	Eigen::Matrix3d m;
	Eigen::Matrix3d m_it;
	double alpha = uc.alpha* 3.141592/180.0;
	double beta = uc.beta* 3.141592/180.0;
	double gamma = uc.gamma* 3.141592/180.0;
	double omega = uc.a*uc.b*uc.c*sqrt(1 - cos(alpha)*cos(alpha) - cos(beta)*cos(beta) - cos(gamma)*cos(gamma) + 2*cos(alpha)*cos(beta)*cos(gamma));
	m << uc.a, uc.b*cos(gamma), uc.c*cos(beta), 0, uc.b*sin(gamma),  uc.c*(cos(alpha)-cos(beta)*cos(gamma))/sin(gamma), 0, 0, omega/(uc.a*uc.b*sin(gamma));
	//m_it = m.inverse().transposeInPlace();
	Eigen::Vector3d av, bv, cv;
	av << uc.a, 0, 0;
	bv << uc.b*cos(gamma), uc.b*sin(gamma), 0;
	cv << uc.c*cos(beta), uc.c*(cos(alpha)-cos(gamma)*cos(beta))/sin(gamma), uc.c*sqrt(1 - cos(beta)*cos(beta) - pow((cos(alpha) - cos(gamma)*cos(beta)),2)	/ (sin(gamma)*sin(gamma)) );
	av = nx*av;
	bv = ny*bv;
	cv = nz*cv;
	int index_counter = 0;
	int res_counter = 0;
	for(int i = 0; i < nx; i++)
	for(int j = 0; j < ny; j++)
	for(int k = 0; k < nz; k++)
	{
		Eigen::Vector3d vc, xc;
		vc << (double)i, (double)j, (double)k;
		xc = m*vc;
		std::vector<Atom> atom_step = get_shifted_relative(xc[0], xc[1], xc[2], uc, m);
		for(auto& atom : atom_step){
			atom.resnr += res_counter;
			atom.index += index_counter;
		}
		atoms.insert(atoms.end(), atom_step.begin(), atom_step.end());
		index_counter += atom_step.size();
		res_counter += 1;
	}
	
	box_dims.push_back(av[0]);
	box_dims.push_back(bv[1]);
	box_dims.push_back(cv[2]);
	box_dims.push_back(av[1]);
	box_dims.push_back(av[2]);
	box_dims.push_back(bv[0]);
	box_dims.push_back(bv[2]);
	box_dims.push_back(cv[0]);
	box_dims.push_back(cv[1]);
	
	return atoms;
}
//expand box such that all of the atoms fit within an orthorhombic box
void makeOrthorhombic(Box& box){
	Vec3<double> min_x; min_x.fill(std::numeric_limits<double>::max());
	Vec3<double> max_x; max_x.fill(std::numeric_limits<double>::min());
	for(const auto& atom : box.atoms){
		for(int i = 0; i < 3; i++){
			if(atom.x[i] < min_x[i]) min_x[i] = atom.x[i];
			if(atom.x[i] > max_x[i]) max_x[i] = atom.x[i];
		}
	}
	//remove off-diagonals and recale box
	for(int i = 0; i < 3; i++){
		for(int j = 0; j < 3; j++){
			if(j == i){
				box.boxvec[i][i] = max_x[i] - min_x[i];
			}
			else box.boxvec[i][j] = 0.0;
		}
	}
	//place all atoms in box
	for(auto& atom : box.atoms){
		for(int i = 0; i < 3; i++){
			atom.x[i] -= min_x[i];
		}
	}
	return;
}

void boxtools::actions::supercell(GroManipData& data, const std::vector<std::string>& args){
  //takes a filename and an output box name
  FANCY_ASSERT(args.size() == 5, "Invalid call to boxtools::actions::supercell(), requires zcif filename, nx, ny, nz, and output_box_name");
  Box* b_out = data.findBox(args[4]);
  if(b_out == 0){
    b_out = new Box;
  }

	auto ucell = load_crystal(args[0]);
	std::vector<double> box_dims;
	std::vector<Atom> atoms = replicate_cells(std::stoi(args[1]), std::stoi(args[2]), std::stoi(args[3]), ucell, box_dims);
	std::cout << "Atom count: " << atoms.size() << std::endl;
	b_out->atoms = atoms;
	b_out->hasNamedAtoms = 1;
	makeOrthorhombic(*b_out);
  data.addBox(args[4], b_out);
  return;	
}
//special procedure for decorating a feldspar surface with appropriate OH groups and assigning atom types
//Arg Order: box_in [atom names to consider] smearing_length density_threshold bond_distance box_out
void boxtools::actions::decoratefeldspar(GroManipData& data, const std::vector<std::string>& args){
  //takes a an input box name and a filename
  FANCY_ASSERT(args.size() == 6, "Invalid call to boxtools::actions::decoratefeldspar(), requires box_in [atom names to consider] nn_length density_threshold bond_distance box_out");
  Box* b1 = data.findBox(args[0]);
  FANCY_ASSERT(b1 != 0, "Failed to find input box in boxtools::actions::writegro()"); 
	std::string output_box_name = args[5];
  Box* b_out = data.findBox(output_box_name);
  if(b_out == 0){
    b_out = new Box;
  }
	Box box_out = *b1;

	//other input args
	std::set<std::string> name_set;
	std::string name_data = args[1].substr(args[1].find('[')+1, args[1].find(']') - args[1].find('[') - 1);
	std::stringstream ss(name_data);
	while(std::getline(ss, name_data, ',')){
		name_set.insert(name_data);
	}
	double cutoff = std::stod(args[2]);
	int threshold = std::stoi(args[3]);
	double bond_distance = std::stod(args[4]);
	//create cellgrid
	std::array<double,3> box_size = {box_out.boxvec[0][0], box_out.boxvec[1][1], box_out.boxvec[2][2]};
	CellGrid c1(cutoff, box_size);
	//cell grid uses actual true indices of atoms in atoms vector
	int index = 0;
	std::array<double,3> com; com.fill(0.0);
	for(auto& atom : box_out.atoms){	
		c1.addIndexToGrid(index, atom.x);
		com = com + atom.x;
		index++;
	}
	com = com * (1.0/(double)box_out.atoms.size());
	index = 0;
	std::vector<int> hydroxyl_idx;
	for(const auto& atom : box_out.atoms){
		if(name_set.find(atom.name) != name_set.end()){
			int nn_count = 0;
			auto indices = c1.getNearbyIndices(index, atom.x);
			for(auto i : indices){
				double dist = getDistance(box_out.atoms[i].x, atom.x, box_size);
				if(dist <= cutoff) nn_count++;
			}
			if(nn_count <= threshold) hydroxyl_idx.push_back(index);
			std::cout << nn_count << std::endl;
		}
		index++;
	}
	//now I have a list of oxygen atoms with hydroxyls, need to add hydrogens to residues in approximately the correct location
	//need to draw a vector from the com to the oxygen atom, then go 1 bond length past that
	//also need to change atom names to reflect this and make sure not to have name clashes
	int activeResNum = -1;
	int resCnt = 0;
	std::vector<int> insert_idx;
	std::vector<Atom> insert_atoms;
	for(auto i : hydroxyl_idx){
		auto& atom = box_out.atoms[i];
		int resNum = atom.resnr;
		if(resNum == activeResNum) resCnt++;
		else{
			//we have a new residue, atoms need unique names for bonds within a residue
			activeResNum = resNum;
			resCnt = 0;
		}
		atom.name = "OH" + std::to_string(resCnt);
		atom.type = "OH";
		Atom hatom;
		Vec3<double> dx = atom.x - com;
		Vec3<double> h_pos = (1.0+(bond_distance/norm2(dx)))*dx + com;
		hatom.name = "HO" + std::to_string(resCnt);
		hatom.type = "HO";
		hatom.resname = atom.resname;
		hatom.resnr = activeResNum;
		hatom.x = h_pos;
		hatom.v.fill(0.0);
		insert_atoms.push_back(hatom);
		insert_idx.push_back(i+1);
	}
	for(int i = insert_idx.size()-1; i >= 0; i--){
		box_out.atoms.insert(box_out.atoms.begin() + insert_idx[i], insert_atoms[i]);
	}
	shrinkWrap(box_out);
  *b_out = box_out;
  data.addBox(output_box_name, b_out); 
  return;	
}