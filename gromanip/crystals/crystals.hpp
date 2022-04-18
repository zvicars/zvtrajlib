#pragma once
#include "../actions.hpp"
#include "../../tools/pbcfunctions.hpp"
#include "../../tools/cellgrid.hpp"
#include "../../tools/stlmath.hpp"
#include "../../tools/StringTools.hpp"
#include "Eigen/Eigen"
#include <cstdlib>
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