#ifndef SIRS_H
#define SIRS_H

#include <vector>
#include <iostream>
#include <fstream>
#include <armadillo>
#include <iomanip>
#include <cmath>
#include <string>

using namespace arma;
using namespace std;

class SIRS {

    vec p; // population (susceptible, infected, recovered)
    vec p0; // initial values for S, I, R

    char taskCase; // a, b, c, d, e
    string filename;

    int N; // whole population

    double a; // the rate of transmission
    double b; // the rate of recovery
    double c; // the rate of immunity loss
    double t; // total time
    double dt; // time step

    double e; // birth rate
    double d; // death rate
    double di; // death rate of infected people

    // Number of death from all groups: dS, dI, dR, dIi
    int dS; // death from susceptible group
    int dI; // death from infected group
    int dR; // death from recovery group
    int dIi; // death from infected group due to the disease

    int bornCount; // all born babies are initially in susceptible group (S)

    string method; // Monte Carlo or Runge Kutta

    double w; // frequency of oscillation
    double f; // vaccination rate

    public:
        SIRS(char taskCase, string filename, vec p0, double a, double b, double c, double t, double dt, double e=0., double d=0., double di=0., double w = 0., double f = 0. );
        void rungeKutta(char populationCase);
        void monteCarlo(char populationCase, int num_sim);

    private:
        void rungeKuttaSingleStep();
        vec func(vec dp);
        void printToFile(ofstream &out, double t, double S, double I, double R, double deathRegular = 0., double deathDisease=0., double bornInS = 0 );
        void monteCarloSingleStep();
        double randNumber();
        void SeasonalVariatio(int i);
        void harmonicFunction(int i);
};

#endif