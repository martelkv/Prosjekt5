import numpy as np
import matplotlib.pyplot as plt

a = False
b = False
c_ODE = False
c_MonteCarlo = False
d_ODE = False
d_MonteCarlo = False
e_ODE = False
e_MonteCarlo = False
std = False
std_rungakutta= True
mean = False


if a:
    tid, S, I, R = np.loadtxt("../outputs/rungeKutta_A.txt", usecols=(0,1,2,3), unpack =True)
    plt.plot(tid,S, label="S", color='#F5B14C')
    plt.plot(tid,I, label ="I", color="forestgreen")
    plt.plot(tid,R, label="R", color="salmon")
    plt.title("Runga Kutta, Population A")
    plt.ylabel("Number of people")
    plt.xlabel("Time [days]")
    plt.legend()
    plt.show()

    tid, S, I, R = np.loadtxt("../outputs/rungeKutta_B.txt", usecols=(0,1,2,3), unpack =True)
    plt.plot(tid,S, label="S", color='#F5B14C')
    plt.plot(tid,I, label ="I", color="forestgreen")
    plt.plot(tid,R, label="R", color="salmon")
    plt.title("Runga Kutta, Population B")
    plt.ylabel("Number of people")
    plt.xlabel("Time [days]")
    plt.legend()
    plt.show()

    tid, S, I, R = np.loadtxt("../outputs/rungeKutta_C.txt", usecols=(0,1,2,3), unpack =True)
    plt.plot(tid,S, label="S", color='#F5B14C')
    plt.plot(tid,I, label ="I", color="forestgreen")
    plt.plot(tid,R, label="R", color="salmon")
    plt.title("Runga Kutta, Population C")
    plt.ylabel("Number of people")
    plt.xlabel("Time [days]")
    plt.legend()
    plt.show()

    tid, S, I, R = np.loadtxt("../outputs/rungeKutta_D.txt", usecols=(0,1,2,3), unpack =True)
    plt.plot(tid,S, label="S", color='#F5B14C')
    plt.plot(tid,I, label ="I", color="forestgreen")
    plt.plot(tid,R, label="R", color="salmon")
    plt.title("Runga Kutta, Population D")
    plt.ylabel("Number of people")
    plt.xlabel("Time [days]")
    plt.legend()
    plt.show()

if b:
    tid, S, I, R = np.loadtxt("../outputs/monteCarlo_A.txt", usecols=(0,1,2,3), unpack =True)
    plt.plot(tid,S, label="S, MC", color='tab:orange')
    plt.plot(tid,I, label ="I, MC", color="tab:green")
    plt.plot(tid,R, label="R, MC", color="tab:red")
    tid, S, I, R = np.loadtxt("../outputs/rungeKutta_A.txt", usecols=(0,1,2,3), unpack =True)
    plt.plot(tid,S, label="S, RK4", color='#F5B14C', linestyle= "--")
    plt.plot(tid,I, label ="I, RK4", color="forestgreen",  linestyle= "--")
    plt.plot(tid,R, label="R, RK4", color="salmon",  linestyle= "--")
    plt.title("Monte Carlo, Population A")
    plt.ylabel("Number of people")
    plt.xlabel("Time [days]")
    plt.legend()
    plt.show()

    tid, S, I, R = np.loadtxt("../outputs/monteCarlo_B.txt", usecols=(0,1,2,3), unpack =True)
    plt.plot(tid,S, label="S, MC", color='tab:orange')
    plt.plot(tid,I, label ="I, MC", color="tab:green")
    plt.plot(tid,R, label="R, MC", color="tab:red")
    tid, S, I, R = np.loadtxt("../outputs/rungeKutta_B.txt", usecols=(0,1,2,3), unpack =True)
    plt.plot(tid,S, label="S, RK4", color='#F5B14C', linestyle= "--")
    plt.plot(tid,I, label ="I, RK4", color="forestgreen",  linestyle= "--")
    plt.plot(tid,R, label="R, RK4", color="salmon",  linestyle= "--")
    plt.title("Monte Carlo, Population B")
    plt.ylabel("Number of people")
    plt.xlabel("Time [days]")
    plt.legend()
    plt.show()

    tid, S, I, R = np.loadtxt("../outputs/monteCarlo_C.txt", usecols=(0,1,2,3), unpack =True)
    plt.plot(tid,S, label="S, MC", color='tab:orange')
    plt.plot(tid,I, label ="I, MC", color="tab:green")
    plt.plot(tid,R, label="R, MC", color="tab:red")
    tid, S, I, R = np.loadtxt("../outputs/rungeKutta_C.txt", usecols=(0,1,2,3), unpack =True)
    plt.plot(tid,S, label="S, RK4", color='#F5B14C', linestyle= "--")
    plt.plot(tid,I, label ="I, RK4", color="forestgreen",  linestyle= "--")
    plt.plot(tid,R, label="R, RK4", color="salmon",  linestyle= "--")
    plt.title("Monte Carlo, Population C")
    plt.ylabel("Number of people")
    plt.xlabel("Time [days]")
    plt.legend()
    plt.show()

    tid, S, I, R = np.loadtxt("../outputs/monteCarlo_D.txt", usecols=(0,1,2,3), unpack =True)
    plt.plot(tid,S, label="S, MC", color='tab:orange')
    plt.plot(tid,I, label ="I, MC", color="tab:green")
    plt.plot(tid,R, label="R, MC", color="tab:red")
    tid, S, I, R = np.loadtxt("../outputs/rungeKutta_D.txt", usecols=(0,1,2,3), unpack =True)
    plt.plot(tid,S, label="S, RK4", color='#F5B14C', linestyle= "--")
    plt.plot(tid,I, label ="I, RK4", color="forestgreen",  linestyle= "--")
    plt.plot(tid,R, label="R, RK4", color="salmon",  linestyle= "--")
    plt.title("Monte Carlo, Population D")
    plt.ylabel("Number of people")
    plt.xlabel("Time [days]")
    plt.legend()
    plt.show()

if c_ODE:
    tid, S, I, R = np.loadtxt("../outputs/rungeKutta_vital_dynamics_A.txt", usecols=(0,1,2,3), unpack =True)
    plt.plot(tid,S, label="S", color='#F5B14C')
    plt.plot(tid,I, label ="I", color="forestgreen")
    plt.plot(tid,R, label="R", color="salmon")
    plt.title("Vital dynamics with Runge-Kutta, Population A")
    plt.ylabel("Number of people")
    plt.xlabel("Time [days]")
    plt.legend()
    plt.show()

    tid, S, I, R = np.loadtxt("../outputs/rungeKutta_vital_dynamics_B.txt", usecols=(0,1,2,3), unpack =True)
    plt.plot(tid,S, label="S", color='#F5B14C')
    plt.plot(tid,I, label ="I", color="forestgreen")
    plt.plot(tid,R, label="R", color="salmon")
    plt.title("Vital dynamics with Runge-Kutta, Population B")
    plt.ylabel("Number of people")
    plt.xlabel("Time [days]")
    plt.legend()
    plt.show()

    tid, S, I, R = np.loadtxt("../outputs/rungeKutta_vital_dynamics_C.txt", usecols=(0,1,2,3), unpack =True)
    plt.plot(tid,S, label="S", color='#F5B14C')
    plt.plot(tid,I, label ="I", color="forestgreen")
    plt.plot(tid,R, label="R", color="salmon")
    plt.title("Vital dynamics with Runge-Kutta, Population C")
    plt.ylabel("Number of people")
    plt.xlabel("Time [days]")
    plt.legend()
    plt.show()

    tid, S, I, R = np.loadtxt("../outputs/rungeKutta_vital_dynamics_D.txt", usecols=(0,1,2,3), unpack =True)
    plt.plot(tid,S, label="S", color='#F5B14C')
    plt.plot(tid,I, label ="I", color="forestgreen")
    plt.plot(tid,R, label="R", color="salmon")
    plt.title("Vital dynamics with Runge-Kutta, Population D")
    plt.ylabel("Number of people")
    plt.xlabel("Time [days]")
    plt.legend()
    plt.show()

if c_MonteCarlo:
    tid, S, I, R = np.loadtxt("../outputs/monteCarlo_vital_dynamics_A.txt", usecols=(0,1,2,3), unpack =True)
    plt.plot(tid,S, label="S, MC", color='tab:orange')
    plt.plot(tid,I, label ="I, MC", color="tab:green")
    plt.plot(tid,R, label="R, MC", color="tab:red")
    tid, S, I, R = np.loadtxt("../outputs/rungeKutta_vital_dynamics_A.txt", usecols=(0,1,2,3), unpack =True)
    plt.plot(tid,S, label="S, RK4", color='#F5B14C', linestyle= "--")
    plt.plot(tid,I, label ="I, RK4", color="forestgreen",  linestyle= "--")
    plt.plot(tid,R, label="R, RK4", color="salmon",  linestyle= "--")
    plt.title("Vital dynamics with Runge-Kutta, Population A")
    plt.ylabel("Number of people")
    plt.xlabel("Time [days]")
    plt.legend()
    plt.show()

    tid, S, I, R = np.loadtxt("../outputs/monteCarlo_vital_dynamics_B.txt", usecols=(0,1,2,3), unpack =True)
    plt.plot(tid,S, label="S, MC", color='tab:orange')
    plt.plot(tid,I, label ="I, MC", color="tab:green")
    plt.plot(tid,R, label="R, MC", color="tab:red")
    tid, S, I, R = np.loadtxt("../outputs/rungeKutta_vital_dynamics_B.txt", usecols=(0,1,2,3), unpack =True)
    plt.plot(tid,S, label="S, RK4", color='#F5B14C', linestyle= "--")
    plt.plot(tid,I, label ="I, RK4", color="forestgreen",  linestyle= "--")
    plt.plot(tid,R, label="R, RK4", color="salmon",  linestyle= "--")
    plt.title("Vital dynamics, Population B")
    plt.ylabel("Number of people")
    plt.xlabel("Time [days]")
    plt.legend()
    plt.show()

    tid, S, I, R = np.loadtxt("../outputs/monteCarlo_vital_dynamics_C.txt", usecols=(0,1,2,3), unpack =True)
    plt.plot(tid,S, label="S, MC", color='tab:orange')
    plt.plot(tid,I, label ="I, MC", color="tab:green")
    plt.plot(tid,R, label="R, MC", color="tab:red")
    tid, S, I, R = np.loadtxt("../outputs/rungeKutta_vital_dynamics_C.txt", usecols=(0,1,2,3), unpack =True)
    plt.plot(tid,S, label="S, RK4", color='#F5B14C', linestyle= "--")
    plt.plot(tid,I, label ="I, RK4", color="forestgreen",  linestyle= "--")
    plt.plot(tid,R, label="R, RK4", color="salmon",  linestyle= "--")
    plt.title("Vital dynamics, Population C")
    plt.ylabel("Number of people")
    plt.xlabel("Time [days]")
    plt.legend()
    plt.show()

    tid, S, I, R = np.loadtxt("../outputs/monteCarlo_vital_dynamics_D.txt", usecols=(0,1,2,3), unpack =True)
    plt.plot(tid,S, label="S, MC", color='tab:orange')
    plt.plot(tid,I, label ="I, MC", color="tab:green")
    plt.plot(tid,R, label="R, MC", color="tab:red")
    tid, S, I, R = np.loadtxt("../outputs/rungeKutta_vital_dynamics_D.txt", usecols=(0,1,2,3), unpack =True)
    plt.plot(tid,S, label="S, RK4", color='#F5B14C', linestyle= "--")
    plt.plot(tid,I, label ="I, RK4", color="forestgreen",  linestyle= "--")
    plt.plot(tid,R, label="R, RK4", color="salmon",  linestyle= "--")
    plt.title("Vital dynamics, Population D")
    plt.ylabel("Number of people")
    plt.xlabel("Time [days]")
    plt.legend()
    plt.show()

if d_ODE:
    tid, S, I, R = np.loadtxt("../outputs/rungeKutta_seasonal_variation_A.txt", usecols=(0,1,2,3), unpack =True)
    plt.plot(tid,S, label="S", color='#F5B14C')
    plt.plot(tid,I, label ="I", color="forestgreen")
    plt.plot(tid,R, label="R", color="salmon")
    plt.title("Runga Kutta, Sesonal Variation, Population A")
    plt.ylabel("Number of people")
    plt.xlabel("Time [days]")
    plt.legend()
    plt.show()

    tid, S, I, R = np.loadtxt("../outputs/rungeKutta_seasonal_variation_B.txt", usecols=(0,1,2,3), unpack =True)
    plt.plot(tid,S, label="S", color='#F5B14C')
    plt.plot(tid,I, label ="I", color="forestgreen")
    plt.plot(tid,R, label="R", color="salmon")
    plt.title("Runga Kutta, Sesonal Variation, Population B")
    plt.ylabel("Number of people")
    plt.xlabel("Time [days]")
    plt.legend()
    plt.show()

    tid, S, I, R = np.loadtxt("../outputs/rungeKutta_seasonal_variation_C.txt", usecols=(0,1,2,3), unpack =True)
    plt.plot(tid,S, label="S", color='#F5B14C')
    plt.plot(tid,I, label ="I", color="forestgreen")
    plt.plot(tid,R, label="R", color="salmon")
    plt.title("Runga Kutta, Sesonal Variation, Population C")
    plt.ylabel("Number of people")
    plt.xlabel("Time [days]")
    plt.legend()
    plt.show()

    tid, S, I, R = np.loadtxt("../outputs/rungeKutta_seasonal_variation_D.txt", usecols=(0,1,2,3), unpack =True)
    plt.plot(tid,S, label="S", color='#F5B14C')
    plt.plot(tid,I, label ="I", color="forestgreen")
    plt.plot(tid,R, label="R", color="salmon")
    plt.title("Runga Kutta, Sesonal Variation, Population D")
    plt.ylabel("Number of people")
    plt.xlabel("Time [days]")
    plt.legend()
    plt.show()

if d_MonteCarlo:
    tid, S, I, R = np.loadtxt("../outputs/monteCarlo_seasonal_variation_A.txt", usecols=(0,1,2,3), unpack =True)
    plt.plot(tid,S, label="S, MC", color='tab:orange')
    plt.plot(tid,I, label ="I, MC", color="tab:green")
    plt.plot(tid,R, label="R, MC", color="tab:red")
    tid, S, I, R = np.loadtxt("../outputs/rungeKutta_seasonal_variation_A.txt", usecols=(0,1,2,3), unpack =True)
    plt.plot(tid,S, label="S, RK4", color='#F5B14C', linestyle= "--")
    plt.plot(tid,I, label ="I, RK4", color="forestgreen",  linestyle= "--")
    plt.plot(tid,R, label="R, RK4", color="salmon",  linestyle= "--")
    plt.title("Seasonal Variation, Population A")
    plt.ylabel("Number of people")
    plt.xlabel("Time [days]")
    plt.legend()
    plt.show()

    tid, S, I, R = np.loadtxt("../outputs/monteCarlo_seasonal_variation_B.txt", usecols=(0,1,2,3), unpack =True)
    plt.plot(tid,S, label="S, MC", color='tab:orange')
    plt.plot(tid,I, label ="I, MC", color="tab:green")
    plt.plot(tid,R, label="R, MC", color="tab:red")
    tid, S, I, R = np.loadtxt("../outputs/rungeKutta_seasonal_variation_B.txt", usecols=(0,1,2,3), unpack =True)
    plt.plot(tid,S, label="S, RK4", color='#F5B14C', linestyle= "--")
    plt.plot(tid,I, label ="I, RK4", color="forestgreen",  linestyle= "--")
    plt.plot(tid,R, label="R, RK4", color="salmon",  linestyle= "--")
    plt.title("Seasonal Variation, Population B")
    plt.ylabel("Number of people")
    plt.xlabel("Time [days]")
    plt.legend()
    plt.show()

    tid, S, I, R = np.loadtxt("../outputs/monteCarlo_seasonal_variation_C.txt", usecols=(0,1,2,3), unpack =True)
    plt.plot(tid,S, label="S, MC", color='tab:orange')
    plt.plot(tid,I, label ="I, MC", color="tab:green")
    plt.plot(tid,R, label="R, MC", color="tab:red")
    tid, S, I, R = np.loadtxt("../outputs/rungeKutta_seasonal_variation_C.txt", usecols=(0,1,2,3), unpack =True)
    plt.plot(tid,S, label="S, RK4", color='#F5B14C', linestyle= "--")
    plt.plot(tid,I, label ="I, RK4", color="forestgreen",  linestyle= "--")
    plt.plot(tid,R, label="R, RK4", color="salmon",  linestyle= "--")
    plt.title("Seasonal Variation, Population C")
    plt.ylabel("Number of people")
    plt.xlabel("Time [days]")
    plt.legend()
    plt.show()

    tid, S, I, R = np.loadtxt("../outputs/monteCarlo_seasonal_variation_D.txt", usecols=(0,1,2,3), unpack =True)
    plt.plot(tid,S, label="S, MC", color='tab:orange')
    plt.plot(tid,I, label ="I, MC", color="tab:green")
    plt.plot(tid,R, label="R, MC", color="tab:red")
    tid, S, I, R = np.loadtxt("../outputs/rungeKutta_seasonal_variation_D.txt", usecols=(0,1,2,3), unpack =True)
    plt.plot(tid,S, label="S, RK4", color='#F5B14C', linestyle= "--")
    plt.plot(tid,I, label ="I, RK4", color="forestgreen",  linestyle= "--")
    plt.plot(tid,R, label="R, RK4", color="salmon",  linestyle= "--")
    plt.title("Seasonal Variation, Population D")
    plt.ylabel("Number of people")
    plt.xlabel("Time [days]")
    plt.legend()
    plt.show()

if e_ODE:
    tid, S, I, R = np.loadtxt("../outputs/rungeKutta_vaccination_A.txt", usecols=(0,1,2,3), unpack =True)
    plt.plot(tid,S, label="S", color='#F5B14C')
    plt.plot(tid,I, label ="I", color="forestgreen")
    plt.plot(tid,R, label="R", color="salmon")
    plt.title("Runga Kutta, Vaccination, Population A")
    plt.ylabel("Number of people")
    plt.xlabel("Time [days]")
    plt.legend()
    plt.show()

    tid, S, I, R = np.loadtxt("../outputs/rungeKutta_vaccination_B.txt", usecols=(0,1,2,3), unpack =True)
    plt.plot(tid,S, label="S", color='#F5B14C')
    plt.plot(tid,I, label ="I", color="forestgreen")
    plt.plot(tid,R, label="R", color="salmon")
    plt.title("Runga Kutta, Vaccination, Population B")
    plt.ylabel("Number of people")
    plt.xlabel("Time [days]")
    plt.legend()
    plt.show()

    tid, S, I, R = np.loadtxt("../outputs/rungeKutta_vaccination_C.txt", usecols=(0,1,2,3), unpack =True)
    plt.plot(tid,S, label="S", color='#F5B14C')
    plt.plot(tid,I, label ="I", color="forestgreen")
    plt.plot(tid,R, label="R", color="salmon")
    plt.title("Runga Kutta, Vaccination, Population C")
    plt.ylabel("Number of people")
    plt.xlabel("Time [days]")
    plt.legend()
    plt.show()

    tid, S, I, R = np.loadtxt("../outputs/rungeKutta_vaccination_D.txt", usecols=(0,1,2,3), unpack =True)
    plt.plot(tid,S, label="S", color='#F5B14C')
    plt.plot(tid,I, label ="I", color="forestgreen")
    plt.plot(tid,R, label="R", color="salmon")
    plt.title("Runga Kutta, Vaccination, Population D")
    plt.ylabel("Number of people")
    plt.xlabel("Time [days]")
    plt.legend()
    plt.show()

if e_MonteCarlo:
    tid, S, I, R = np.loadtxt("../outputs/monteCarlo_vaccination_A.txt", usecols=(0,1,2,3), unpack =True)
    plt.plot(tid,S, label="S, MC", color='tab:orange')
    plt.plot(tid,I, label ="I, MC", color="tab:green")
    plt.plot(tid,R, label="R, MC", color="tab:red")
    tid, S, I, R = np.loadtxt("../outputs/rungeKutta_vaccination_A.txt", usecols=(0,1,2,3), unpack =True)
    plt.plot(tid,S, label="S, RK4", color='#F5B14C', linestyle= "--")
    plt.plot(tid,I, label ="I, RK4", color="forestgreen",  linestyle= "--")
    plt.plot(tid,R, label="R, RK4", color="salmon",  linestyle= "--")
    plt.title("Vaccination, Population A")
    plt.ylabel("Number of people")
    plt.xlabel("Time [days]")
    plt.legend()
    plt.show()

    tid, S, I, R = np.loadtxt("../outputs/monteCarlo_vaccination_B.txt", usecols=(0,1,2,3), unpack =True)
    plt.plot(tid,S, label="S, MC", color='tab:orange')
    plt.plot(tid,I, label ="I, MC", color="tab:green")
    plt.plot(tid,R, label="R, MC", color="tab:red")
    tid, S, I, R = np.loadtxt("../outputs/rungeKutta_vaccination_B.txt", usecols=(0,1,2,3), unpack =True)
    plt.plot(tid,S, label="S, RK4", color='#F5B14C', linestyle= "--")
    plt.plot(tid,I, label ="I, RK4", color="forestgreen",  linestyle= "--")
    plt.plot(tid,R, label="R, RK4", color="salmon",  linestyle= "--")
    plt.title("Vaccination, Population B")
    plt.ylabel("Number of people")
    plt.xlabel("Time [days]")
    plt.legend()
    plt.show()

    tid, S, I, R = np.loadtxt("../outputs/monteCarlo_vaccination_C.txt", usecols=(0,1,2,3), unpack =True)
    plt.plot(tid,S, label="S, MC", color='tab:orange')
    plt.plot(tid,I, label ="I, MC", color="tab:green")
    plt.plot(tid,R, label="R, MC", color="tab:red")
    tid, S, I, R = np.loadtxt("../outputs/rungeKutta_vaccination_C.txt", usecols=(0,1,2,3), unpack =True)
    plt.plot(tid,S, label="S, RK4", color='#F5B14C', linestyle= "--")
    plt.plot(tid,I, label ="I, RK4", color="forestgreen",  linestyle= "--")
    plt.plot(tid,R, label="R, RK4", color="salmon",  linestyle= "--")
    plt.title("Vaccination, Population C")
    plt.ylabel("Number of people")
    plt.xlabel("Time [days]")
    plt.legend()
    plt.show()

    tid, S, I, R = np.loadtxt("../outputs/monteCarlo_vaccination_D.txt", usecols=(0,1,2,3), unpack =True)
    plt.plot(tid,S, label="S, MC", color='tab:orange')
    plt.plot(tid,I, label ="I, MC", color="tab:green")
    plt.plot(tid,R, label="R, MC", color="tab:red")
    tid, S, I, R = np.loadtxt("../outputs/rungeKutta_vaccination_D.txt", usecols=(0,1,2,3), unpack =True)
    plt.plot(tid,S, label="S, RK4", color='#F5B14C', linestyle= "--")
    plt.plot(tid,I, label ="I, RK4", color="forestgreen",  linestyle= "--")
    plt.plot(tid,R, label="R, RK4", color="salmon",  linestyle= "--")
    plt.title("Vaccination, Population D")
    plt.ylabel("Number of people")
    plt.xlabel("Time [days]")
    plt.legend()
    plt.show()

if std:
    tid, S, I, R = np.loadtxt("../outputs/monteCarlo_A.txt", usecols=(0,1,2,3), unpack =True)
    std_S = np.std(S[1200:])/(len(S[1200:]))
    std_I = np.std(I[1200:])/(len(S[1200:]))
    std_R = np.std(R[1200:])/(len(S[1200:]))
    print("Monte Carlo______Population A _______")
    print("std S: ",std_S)
    print("std I: ",std_I)
    print("std R: ",std_R)
    tid, S, I, R = np.loadtxt("../outputs/monteCarlo_B.txt", usecols=(0,1,2,3), unpack =True)
    std_S = np.std(S[1600:])/(len(S[1600:]))
    std_I = np.std(I[1600:])/(len(S[1600:]))
    std_R = np.std(R[1600:])/(len(S[1600:]))
    print("Monte Carlo______Population B _______")
    print("std S: ",std_S)
    print("std I: ",std_I)
    print("std R: ",std_R)
    tid, S, I, R = np.loadtxt("../outputs/monteCarlo_C.txt", usecols=(0,1,2,3), unpack =True)
    std_S = np.std(S[1600:])/(len(S[1600:]))
    std_I = np.std(I[1600:])/(len(S[1600:]))
    std_R = np.std(R[1600:])/(len(S[1600:]))
    print("Monte Carlo______Population C _______")
    print("std S: ",std_S)
    print("std I: ",std_I)
    print("std R: ",std_R)
    tid, S, I, R = np.loadtxt("../outputs/monteCarlo_D.txt", usecols=(0,1,2,3), unpack =True)
    std_S = np.std(S[1600:])/(len(S[1600:]))
    std_I = np.std(I[1600:])/(len(S[1600:]))
    std_R = np.std(R[1600:])/(len(S[1600:]))
    print("Monte Carlo______Population B _______")
    print("std S: ",std_S)
    print("std I: ",std_I)
    print("std R: ",std_R)

if mean:
    tid, S, I, R = np.loadtxt("../outputs/monteCarlo_A.txt", usecols=(0,1,2,3), unpack =True)
    mean_S = np.mean(S[1200:])#/(len(S[1200:]))
    mean_I = np.mean(I[1200:])#/(len(S[1200:]))
    mean_R = np.mean(R[1200:])#/(len(S[1200:]))
    print("Mean MC______Population A _______")
    print("mean S: ",mean_S)
    print("mean I: ",mean_I)
    print("mean R: ",mean_R)
    tid, S, I, R = np.loadtxt("../outputs/monteCarlo_B.txt", usecols=(0,1,2,3), unpack =True)
    mean_S = np.mean(S[1600:])#/(len(S[1200:]))
    mean_I = np.mean(I[1600:])#/(len(S[1200:]))
    mean_R = np.mean(R[1600:])#/(len(S[1200:]))
    print("Mean MC______Population B _______")
    print("mean S: ",mean_S)
    print("mean I: ",mean_I)
    print("mean R: ",mean_R)
    tid, S, I, R = np.loadtxt("../outputs/monteCarlo_C.txt", usecols=(0,1,2,3), unpack =True)
    mean_S = np.mean(S[1600:])#/(len(S[1200:]))
    mean_I = np.mean(I[1600:])#/(len(S[1200:]))
    mean_R = np.mean(R[1200:])#/(len(S[1200:]))
    print("Mean MC______Population C _______")
    print("mean S: ",mean_S)
    print("mean I: ",mean_I)
    print("mean R: ",mean_R)
    tid, S, I, R = np.loadtxt("../outputs/monteCarlo_D.txt", usecols=(0,1,2,3), unpack =True)
    mean_S = np.mean(S[1600:])#/(len(S[1200:]))
    mean_I = np.mean(I[1600:])#/(len(S[1200:]))
    mean_R = np.mean(R[1600:])#/(len(S[1200:]))
    print("Mean MC______Population D _______")
    print("mean S: ",mean_S)
    print("mean I: ",mean_I)
    print("mean R: ",mean_R)


if std_rungakutta:
    tid, S, I, R = np.loadtxt("../outputs/rungeKutta_seasonal_variation_A.txt", usecols=(0,1,2,3), unpack =True)
    std_S = np.std(S[1200:])/(len(S[1200:]))
    std_I = np.std(I[1200:])/(len(S[1200:]))
    std_R = np.std(R[1200:])/(len(S[1200:]))
    print("RK______Population A _______")
    print("std S: ",std_S)
    print("std I: ",std_I)
    print("std R: ",std_R)
    tid, S, I, R = np.loadtxt("../outputs/rungeKutta_seasonal_variation_B.txt", usecols=(0,1,2,3), unpack =True)
    std_S = np.std(S[1600:])/(len(S[1600:]))
    std_I = np.std(I[1600:])/(len(S[1600:]))
    std_R = np.std(R[1600:])/(len(S[1600:]))
    print("RK______Population B _______")
    print("std S: ",std_S)
    print("std I: ",std_I)
    print("std R: ",std_R)
    tid, S, I, R = np.loadtxt("../outputs/rungeKutta_seasonal_variation_C.txt", usecols=(0,1,2,3), unpack =True)
    std_S = np.std(S[1600:])/(len(S[1600:]))
    std_I = np.std(I[1600:])/(len(S[1600:]))
    std_R = np.std(R[1600:])/(len(S[1600:]))
    print("RK______Population C _______")
    print("std S: ",std_S)
    print("std I: ",std_I)
    print("std R: ",std_R)
    tid, S, I, R = np.loadtxt("../outputs/rungeKutta_seasonal_variation_D.txt", usecols=(0,1,2,3), unpack =True)
    std_S = np.std(S[1600:])/(len(S[1600:]))
    std_I = np.std(I[1600:])/(len(S[1600:]))
    std_R = np.std(R[1600:])/(len(S[1600:]))
    print("RK______Population D _______")
    print("std S: ",std_S)
    print("std I: ",std_I)
    print("std R: ",std_R)
