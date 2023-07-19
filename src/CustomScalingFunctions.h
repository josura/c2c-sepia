#pragma once
#include <functional>

//add scaling function as lambda function for dissipation
std::function<double(double)> getDissipationScalingFunction();

//add scaling function as lambda function for conservation
std::function<double(double)> getConservationScalingFunction();