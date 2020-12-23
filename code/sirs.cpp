#include "sirs.h"

SIRS::SIRS(char taskCase, string filename, vec p0, double a, double b, double c, double t, double dt, double e, double d, double di, double w, double f) {
    p = p0; // initialize vector p with 3 elements
    this->p0 = p0;
    this->taskCase = taskCase;
    this->filename = filename;

    N = p(0) + p(1) + p(2);

    this->a = a;
    this->b = b;
    this->c = c;
    this->t = t;
    this->dt = dt;
    this->e = e;
    this->d = d;
    this->di = di;
    this->w = w;
    this->f = f;
}

vec SIRS::func(vec dp) {
    vec t = vec(3, fill::zeros);

    switch(taskCase) {
        case 'a': case 'd': 
        {
            // S' = c*(N-S-I) -a*S*I/N
            // I' = a*S*I/N - b*I
            // R = N - S - I => S' = c*R - a*S*I/N
            t(0) = c*dp(2) - a*dp(1)*dp(0)/N;
            t(1) = a*dp(0)*dp(1)/N - b*dp(1);
            break;
         }
         case 'c': 
         {
            // S' = c*R - a*S*I/N -d*S + e*N
            // I' = a*S*I/N - b*I - d*I - di*I
            // R' = b*I - c*R - d*R
            t(0) = c*dp(2) - dp(0) * (a * dp(1) / N  + d)  + e * N;
            t(1) = dp(1)*(a*dp(0)/N - b-d-di);
            t(2) = b*dp(1) - dp(2)*(c+d);
            break;
         }
         case 'e': case 'f':
         {
            // S' = c*R - a*S*I/N - f
            // I' = a*S*I/N - b*I
            // R' = b*I - c*R + f
            t(0) = c * dp(2) - a * dp(0) * dp(1) / N - f;
            t(1) = a * dp(0) * dp(1) / N - b * dp(1);
            t(2) = b * dp(1) - c * dp(2) + f;

            break;
         }
         default: {
            cout << "Case does not exist!!!" << endl;
            break;
         }
    }

    return t;
}

void SIRS::rungeKutta(char populationCase){
    cout << "Runge Kutta - task " << taskCase << flush;
    method = "RK";

    ofstream out;
    out.open(filename+"_" + populationCase+".txt");

    int n = int(t/dt);
    printToFile(out, 0, p(0), p(1), p(2));

    int dots = n / 10; // number of dots - simulate progress bar

    for(int i = 0; i < n; i++){
        // calculating parameter a(t) - Seasonal Variation
        if(taskCase == 'd') {
            SeasonalVariatio(i);
        }
        if(taskCase == 'f') {
            harmonicFunction(i);
        }

        rungeKuttaSingleStep();
        printToFile(out, (i+1)*dt, p(0), p(1), p(2));

         if(i % dots == 0) { // progress bar
            cout << "." << flush; // force showing '.' on the screen
         }
    }

    cout << endl << flush;

    out.close();
}

void SIRS::rungeKuttaSingleStep() {
    vec K1(3), K2(3), K3(3), K4(3);

    K1 = dt*func(p);
    K2 = dt*func(p + K1/2);
    K3 = dt*func(p + K2/2);
    K4 = dt*func(p+K3);
    p += (K1 + 2*K2 + 2*K3 + K4)/6.0;
    if(taskCase == 'a') {
        p(2) = N - p(1) - p(0);
    } else {
        N =  p(0) + p(1) + p(2); // p(2) is calculatinf in the function f
    }
}

void SIRS::printToFile(ofstream &out, double t, double S, double I, double R, double deathRegular, double deathDisease, double bornInS) {
    out << setw(15) << setprecision(8) << t;
    out << setw(15) << setprecision(8) << S;
    out << setw(15) << setprecision(8) << I;
    out << setw(15) << setprecision(8) << R;
    if((taskCase == 'c') && method == "MC") {
        out << setw(15) << setprecision(8) << deathRegular;
        out << setw(15) << setprecision(8) << deathDisease;
        out << setw(15) << setprecision(8) << bornInS;
    }
    out << endl;
}

void SIRS::monteCarlo(char populationCase, int num_sim){
    cout << "Monte Carlo - task " << taskCase << flush;
    method = "MC";

    int dots = num_sim / 10; // number of dots - simulate progress bar
    ofstream out;
    out.open(filename+"_" + populationCase+".txt");

    dt = min({4/a/N , 1/b/N, 1/c/N});
    int n = int(t/dt);

    printToFile(out, 0, p(0), p(1), p(2));

    vec S = vec(n, fill::zeros);
    vec I = vec(n, fill::zeros);
    vec R = vec(n, fill::zeros);

    // caseTask c
    vec deathRegular = vec(n, fill::zeros);
    vec deathDisease = vec(n, fill::zeros);
    vec bornInS = vec(n, fill::zeros);

    for (int i = 0; i < num_sim; i ++) {
        for (int j = 1; j < n; j++){
            // calculating parameter a(t) - Seasonal Variation
            if(taskCase == 'd') {
                SeasonalVariatio(j);
            }
            if(taskCase == 'f') {
                harmonicFunction(j);
            }

            monteCarloSingleStep();

            S(j) += p(0);
            I(j) += p(1);
            R(j) += p(2);

            deathRegular(j) += (dS + dI + dR);
            deathDisease(j) += dIi;
            bornInS(j) += bornCount;
        }

        p = p0; // set to initial values

        if(i % dots == 0) { // progress bar
            cout << "." << flush; // force showing '.' on the screen
        }
    }

    cout << endl << flush;

    S /= num_sim;
    I /= num_sim;
    R /= num_sim;
    deathRegular /= num_sim;
    deathDisease /= num_sim;
    bornInS /= num_sim;

    for(int i =1; i < n; i++) {
        printToFile(out, i*dt, S(i), I(i), R(i), deathRegular(i), deathDisease(i), bornInS(i));
    }
}

void SIRS::monteCarloSingleStep() {

    //Calculate tansition probabilities
    double S_I = a * p(0) * p(1) / N * dt;
    double I_R = b * p(1) * dt;
    double R_S = c * p(2) * dt;

    // count transition probabilities and transfer from/to relevant group
    int n_SI = 0; // from S->I
    int n_IR = 0; // from I->R
    int n_RS = 0; // from R->S

    dS = dI = dR = dIi = bornCount = 0;

    if(randNumber() < S_I) {
        ++n_SI;
    }
    if(randNumber() < I_R) {
        ++n_IR;
    }
    if(randNumber() < R_S) {
        ++n_RS;
    }

    if (taskCase == 'c') {
        if(randNumber() < d*p(0)*dt) {
            ++dS; // death from susceptible group
        }
        if(randNumber() < d*p(1)*dt) {
            ++dI; // death from infected group
        }
        if(randNumber() < d*p(2)*dt) {
            ++dR; // death from recovery group
        }
        if(randNumber() < di*p(1)*dt) {
            ++dIi; // death from infected group due to the disease
        }
        if(randNumber() < e*N*dt) {
            ++bornCount; // born are initially in susceptible group
        }
    }

    if(taskCase == 'e' || taskCase == 'f') {
        if(randNumber() < f*dt && p(0) > 0) {
            p(0) -= 1;
            p(2) += 1;
        }
    }


    p(0) += (n_RS - n_SI + bornCount - dS);
    p(1) += (n_SI - n_IR - dI - dIi);
    p(2) += (n_IR - n_RS - dR);

    N = p(0) + p(1) + p(2);
}

// random number generator between 0 and 1
double SIRS::randNumber() {
    double num = rand() % 10001;
    num /= 10000.0;

    return num;
}

// calculate parameter a(t) - Seasonal Variation
void SIRS::SeasonalVariatio(int i) {
    double A = 1.5; // max deviattion from a0
    double a0 = 4; // avarage transmission rate

    a = A * cos(w * dt * i) + a0;
}

// calculate parameter f(t) - Harmonic function
void SIRS::harmonicFunction(int i) {
    double F = 40; // max deviattion from f0
    double f0 = 100; // avarage transmission rate

    f = F * cos(w * dt * i+ 3) + f0;
}
