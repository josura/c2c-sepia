#include "CustomScalingFunctions.h"

//add scaling function as lambda function for dissipation
std::function<double(double)> getDissipationScalingFunction(){
    return [](double x) {
        //rewrite accordingly
        return 1.0;      
    };
}

//add scaling function as lambda function for conservation
std::function<double(double)> getConservationScalingFunction(){
    return [](double x) {
        //rewrite accordingly
        return 1.0;      
    };
}

std::function<double(double)> getPropagationScalingFunction(){
    return [](double x) {
        //rewrite accordingly
        return 1.0;      
    };
}