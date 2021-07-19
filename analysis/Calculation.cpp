#include "Calculation.hpp"

Calculation::Calculation(InputPack& input)
{
    input.params().readString("name", KeyType::Required, name_);
    input.params().readString("type", KeyType::Required, type_);
    input.addCalculation(name_, this);
    return;
}