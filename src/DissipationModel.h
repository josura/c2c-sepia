#pragma once
#include <armadillo>

class DissipationModel{
    DissipationModel();
    ~DissipationModel();
    public:
        arma::Col<double> dissipatePow2Self(arma::Col<double> input);
        arma::Col<double> dissipateSelfScaled(arma::Col<double> input, double scale);
        arma::Col<double> dissipateSelfScaled(arma::Col<double> input, arma::Col<double> scales);
        arma::Col<double> dissipateSelfPeriodic(arma::Col<double> input, double period, double amplitude, double phase);
        arma::Col<double> dissipateSelfPeriodic(arma::Col<double> input, arma::Col<double> periods, arma::Col<double> amplitudes, arma::Col<double> phases);
        arma::Col<double> dissipateSelfPeriodicShift(arma::Col<double> input, double shiftY, double period, double amplitude, double phase);
        arma::Col<double> dissipateSelfPeriodicShift(arma::Col<double> input, arma::Col<double> shiftY, arma::Col<double> periods, arma::Col<double> amplitudes, arma::Col<double> phases);
        arma::Col<double> dissipateSelfRandom(arma::Col<double> input, double min, double max);
};