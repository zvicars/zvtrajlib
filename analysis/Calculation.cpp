#include "Calculation.hpp"

Calculation::Calculation(InputPack& input)
{
    input.params().readString("name", KeyType::Required, name_);
    input.params().readString("type", KeyType::Required, type_);

    base_ = name_;
    input.params().readString("base", KeyType::Optional, base_);
    
    output_freq_ = -1; //will never output
    input.params().readNumber("frequency", KeyType::Optional, output_freq_);
    
    equilibration_ = 0;
    //equilibration in ps
    input.params().readNumber("equilibration", KeyType::Optional, equilibration_);
    
    input.addCalculation(name_, this);
    return;
}