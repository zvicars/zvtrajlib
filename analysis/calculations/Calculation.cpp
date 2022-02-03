#include "Calculation.hpp"
#include <limits>
Calculation::Calculation(InputPack& input)
{
    input.params().readString("name", KeyType::Required, name_);
    input.params().readString("type", KeyType::Required, type_);
    base_ = name_;
    input.params().readString("base", KeyType::Optional, base_);
    output_freq_ = -1; //will never output
    input.params().readNumber("output_frequency", KeyType::Optional, output_freq_);
    if(output_freq_ > 0) calc_freq_ = output_freq_;
    else calc_freq_ = 1;
    input.params().readNumber("calculation_frequency", KeyType::Optional, calc_freq_);
    equilibration_ = 0;
    FANCY_ASSERT( output_freq_%calc_freq_ == 0 || output_freq_ < 0, "Can only output data on a calculation step! Output frequency should be a multiple of calculation frequency.");
    //equilibration in ps
    input.params().readNumber("equilibration", KeyType::Optional, equilibration_);
    end_ = std::numeric_limits<double>::max();
    input.params().readNumber("end", KeyType::Optional, end_);
    input.addCalculation(name_, this);
    box = input.getBox();
    return;
}

bool Calculation::doCalculate(){
    if(current_time_ < equilibration_ || current_time_ > end_) return 0;
    if(current_frame_%calc_freq_ != 0) return 0;
    return 1;
}
bool Calculation::doOutput(){
    if(current_time_ < equilibration_ || current_time_ > end_) return 0;
    if(current_frame_%output_freq_ != 0) return 0;
    return 1;
}

void Calculation::update(){
    current_time_ = box->time;
    current_frame_ = box->frame_counter; 
    return;
}