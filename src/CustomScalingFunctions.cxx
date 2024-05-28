#include "CustomScalingFunctions.h"


//add scaling function as lambda function for dissipation, value should be between 0 and 1 to not get negative values
std::function<double(double)> getDissipationScalingFunction(){
    return [](double x) {
        //rewrite accordingly
        return 0.5;
        
    };
}

//add scaling function as lambda function for conservation, should also be between 0 and 1 to not get negative values
std::function<double(double)> getConservationScalingFunction(){
    return [](double x) {
        //rewrite accordingly
        return 0.5;      
    };
}

std::function<double(double)> getPropagationScalingFunction(){
    return [](double x) {
        //rewrite accordingly
        return 0.5;
    };
}