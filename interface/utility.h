#pragma once
#include "datatypes.h"
//contains files that operate on datatypes to extract useful information from the box, i.e. some simple order parameters
inline int OP_Nv(double pp[6], const Box& box, std::string idx = ""){
    int count = 0;
    for(int i = 0; i < box.atoms.size(); i++)
    {
        if(!(box.atoms[i].x[0] >= pp[0] && box.atoms[i].x[0] <= pp[3])) continue;
        if(!(box.atoms[i].x[1] >= pp[1] && box.atoms[i].x[1] <= pp[4])) continue;
        if(!(box.atoms[i].x[2] >= pp[2] && box.atoms[i].x[2] <= pp[5])) continue;        
        count++;
    }
    return count;
}