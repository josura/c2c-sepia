#pragma once
#include <functional>
#include "utilities.h"

//add scaling function as lambda function for dissipation
std::function<double(double)> getDissipationScalingFunction();

//add scaling function as lambda function for conservation
std::function<double(double)> getConservationScalingFunction();

//add scaling function as lambda function for propagation
std::function<double(double)> getPropagationScalingFunction();