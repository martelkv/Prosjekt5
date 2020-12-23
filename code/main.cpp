#include "sirs.h"


int main(int argc, char const *argv[]) {
    vec P0{ 300.0, 100.0, 0.0 }; // initial values for S, I, R

    // constants
    double a = 4.; // the rate of transmission
    double b = 1; // the rate of recovery
    double c = 0.5; // the rate of immunity loss
    double t = 30; // total time
    double dt = 0.005; // time step

    int num_sim = 600; // Monte Carlo - number of simulations

    // constants for tack c.)
    double e = 0.007; // birth rate
    double d = 0.002; // death rate
    double di = 0.03; // death rate of infected people

    double w = 0.5; // frekvency of oscillation
    double f = 100; // vaccination rate

    char population[] = {'A','B','C','D'};

    for(char i : population) {
        cout << "Population: " << i << endl;

        // task a.)
        SIRS rk_a('a', "rungeKutta", P0, a, b, c, t, dt);
        rk_a.rungeKutta(i);

        // task b.)
        SIRS mc_b('b', "monteCarlo", P0, a, b, c, t, dt);
        mc_b.monteCarlo( i, num_sim);

        // task c.)
        SIRS rk_c('c', "rungeKutta_vital_dynamics", P0, a, b, c, t, dt, e, d, di);
        rk_c.rungeKutta(i);

        SIRS mc_c('c', "monteCarlo_vital_dynamics", P0, a, b, c, t, dt, e, d, di);
        mc_c.monteCarlo( i, num_sim);

        // task d.)
        SIRS rk_d('d', "rungeKutta_seasonal_variation", P0, a, b, c, t, dt, e, d, di, w);
        rk_d.rungeKutta(i);

        SIRS mc_d('d', "monteCarlo_seasonal_variation", P0, a, b, c, t, dt, e, d, di, w);
        mc_d.monteCarlo( i, num_sim);

        // task e.)
        SIRS rk_e('e', "rungeKutta_vaccination", P0, a, b, c, t, dt, e, d, di, w, f);
        rk_e.rungeKutta(i);

        SIRS mc_e('e', "monteCarlo_vaccination", P0, a, b, c, t, dt, e, d, di, w, f);
        mc_e.monteCarlo( i, num_sim);

        // task f.)
        SIRS rk_f('f', "rungeKutta_vaccination_dynamic_f", P0, a, b, c, t, dt, e, d, di, w, f);
        rk_f.rungeKutta(i);

        SIRS mc_f('f', "monteCarlo_vaccination_dynamic_f", P0, a, b, c, t, dt, e, d, di, w, f);
        mc_f.monteCarlo( i, num_sim);

        ++b; // parameter b increases for each population by 1
        cout << endl;
    }


    return 0;
}



