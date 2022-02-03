#include "Calc_Histogram.hpp"
Calculation_Histogram::Calculation_Histogram(InputPack& input):Calculation{input}
{
    input.params().readFlag("histogram", KeyType::Optional, doHistogram);
    input.params().readFlag("timeseries", KeyType::Optional, doTimeseries);
    bool found_min = input.params().readNumber("min_bin", KeyType::Optional, min_bin_);
    if(found_min) forceMin = 1;
    bool found_max = input.params().readNumber("max_bin", KeyType::Optional, max_bin_);
    if(found_max) forceMax = 1;
    bool found_bs = input.params().readNumber("bin_size", KeyType::Optional, bin_size_);
    if(found_bs) forceBS = 1;
    if(found_min || found_max || found_bs) doHistogram = 1;
    return;
}