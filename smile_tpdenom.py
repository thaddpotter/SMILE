#%%
import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import odeint, quad, dblquad
from scipy import integrate
from scipy.interpolate import interp1d
import numba as nb
from numba import njit
from numba import cfunc
from numba.types import intc, CPointer, float64
from scipy import LowLevelCallable

#This here is the big SMILE function.  This is the model itself
def Light_curve_final(m1, m2, rstar, ecc, phi, omega, alpha, gamma, g, Periods, ti):
    times= np.linspace(0*86400,14*86400,2000)
    t0 = 3.2 #Pulse Mid-Point
    mu = (m1* (1.988e30) * m2* (1.988e30)) / (m1* (1.988e30) + m2* (1.988e30))  # ratio
    G = 6.6743 * 10 ** (-11)  # Gravitational Constant
    c = 2.99792458e8  # speed of light
    rs = 2 * G * (m1* (1.988e30) + m2* (1.988e30)) / (c ** 2)  # Total mass Schwarzschild Radius
    P = Periods * 86400  # orbital period in seconds
    a = ((P ** 2 * G * (m1* (1.988e30) + m2* (1.988e30))) / (4 * 3.14159265358 ** 2)) ** (1 / 3)  # system semi-major axis
    aph = a * (1 + ecc)  # aphelion

    # initial conditions
    x1_0 = (m2* (1.988e30) * aph) / (m1* (1.988e30) + m2* (1.988e30))  # starting x position of mass 1
    x2_0 = -(m1* (1.988e30) * aph) / (m1* (1.988e30) + m2* (1.988e30))  # starting x position of mass 2
    y1_0 = 0  # starting y position of mass 1
    y2_0 = 0  # starting y position of mass 2
    vel1 = ((G * mu / a) * (m2 / m1) * ((1 - ecc) / (1 + ecc))) ** 0.5  # initial y velocity of mass 1
    vel2 = (-1) * ((G * mu / a) * (m2 / m1) * ((1 - ecc) / (1 + ecc))) ** 0.5 * (m1 / m2)  # initial y velocity of mass 2

    # system of equations
    #@njit(fastmath=True)
    def orbits(S, t, t1):
        X1, Y1, X2, Y2, VX1, VY1, VX2, VY2 = S

        bottom = (((((((X1 - X2) ** 2) + ((Y1 - Y2) ** 2)) ** (1 / 2)) - rs) ** 2) * (
                        (((X1 - X2) ** 2) + ((Y1 - Y2) ** 2)) ** (1 / 2)))
        return [VX1,
                VY1,
                VX2,
                VY2,
                -((G * m2* (1.988e30)) * (X1 - X2)) / bottom,
                -((G * m2* (1.988e30)) * (Y1 - Y2)) / bottom,
                ((G * m1* (1.988e30)) * (X1 - X2)) / bottom,
                ((G * m1* (1.988e30)) * (Y1 - Y2)) / bottom]

    # initial conditions
    X1_0 = x1_0
    Y1_0 = y1_0
    X2_0 = x2_0
    Y2_0 = y2_0
    VX1_0 = 0
    VY1_0 = vel1
    VX2_0 = 0
    VY2_0 = vel2

    # Shift time coordinate by ti
    t_shifted = np.linspace(-3*P, 3*P, 1000)-ti*86400
    S_0 = (X1_0, Y1_0, X2_0, Y2_0, VX1_0, VY1_0, VX2_0, VY2_0)

    ans = odeint(orbits, y0=S_0, t=t_shifted, args=(ti,))

    xp1 = ans.T[0]
    yp1 = ans.T[1]
    vxp1 = ans.T[4]
    vyp1 = ans.T[5]

    xp2 = ans.T[2]
    yp2 = ans.T[3]
    vxp2 = ans.T[6]
    vyp2 = ans.T[7]

    f_xp10 = interp1d(t_shifted, xp1)
    xp10 = f_xp10((t0)*86400)
    f_xp20 = interp1d(t_shifted, xp2)
    xp20 = f_xp20((t0)*86400)
    f_yp10 = interp1d(t_shifted, yp1)
    yp10 = f_yp10((t0)*86400)
    f_yp20 = interp1d(t_shifted, yp2)
    yp20 = f_yp20((t0)*86400)

    f_xp1 = interp1d(t_shifted, xp1)
    xp1 = f_xp1(times)
    f_yp1 = interp1d(t_shifted, yp1)
    yp1 = f_yp1(times)
    f_vxp1 = interp1d(t_shifted, vxp1)
    vxp1 = f_vxp1(times)
    f_vyp1 = interp1d(t_shifted, vyp1)
    vyp1 = f_vyp1(times)
    f_xp2 = interp1d(t_shifted, xp2)
    xp2 = f_xp2(times)
    f_yp2 = interp1d(t_shifted, yp2)
    yp2 = f_yp2(times)
    f_vxp2 = interp1d(t_shifted, vxp2)
    vxp2 = f_vxp2(times)
    f_vyp2 = interp1d(t_shifted, vyp2)
    vyp2 = f_vyp2(times)



    M = m1 * (1.988 * 10 ** 30)  # Mass of the Black Hole (kg)
    ap = np.sqrt((xp2 - xp1) ** 2 + (yp2 - yp1) ** 2)
    ap0 = np.sqrt((xp20 - xp10) ** 2 + (yp20 - yp10) ** 2)
    vtran = np.sqrt(vxp1 ** 2 + vyp1 ** 2) + np.sqrt(vxp2 ** 2 + vyp2 ** 2)
    Rs = (2 * G * M) / c ** 2  # Schwarzschild Radius
    Rein = (2 * Rs * ap) ** 0.5  # Einstein Radius
    Rein0 = (2 * Rs * ap0) ** 0.5 
    te = (Rein / vtran) / 86400  # Einstein Time
    u0 = (ap0 / Rein0) * np.sin(phi * (np.pi / 180)) # Closest Angular Approach
    pstar = ((rstar * (6.957e8))/Rein)

    u=((u0**2)+((((times/86400)-t0)/te))**2)**0.5  # Angular Separation
    
    beta = 0.15 * (1+g) * ((15+gamma)/(3-gamma))
    
    sev = beta * ((m1/m2)*(rstar* (6.957e8)/(np.sqrt((xp2-xp1)**2+(yp2-yp1)**2)))**3*np.cos(2*((2*np.pi/P)*(times+ti*86400)+np.pi+omega))*np.sin((np.pi/2-phi* (np.pi / 180)))**2)
    delf = ((3-(alpha))*((((1/c)*(-1*(vxp2*np.sin(omega)*np.sin((np.pi/2-phi* (np.pi / 180)))+vyp2*np.cos(omega)*np.sin((np.pi/2-phi* (np.pi / 180))))))*np.sin((np.pi/2-phi* (np.pi / 180))))))

    variables=[pstar,u,times]

    
    # take integral of top for all time values

    #This is where the very involved integration begins
    def jit_integrand_function(integrand_function):
        jitted_function = nb.jit(integrand_function, nopython=True)

        #error_model="numpy" -> Don't check for division by zero
        @cfunc(float64(intc, CPointer(float64)),error_model="numpy",fastmath=True)
        def wrapped(n, xx):
            ar = nb.carray(xx, n)
            return jitted_function(ar[0], ar[1], ar[2], ar[3], ar[4])
        return LowLevelCallable(wrapped.ctypes)

    @jit_integrand_function
    def pre_top(r, th, u, pstar, gamma):
        return (1 - gamma * ((1 - (3 / 2) * ((pstar ** 2 - r ** 2) / pstar ** 2) ** 0.5))) * ((r * (u ** 2 + r ** 2 - 2 * u * r * np.cos(th) + 2)) / ( (((u ** 2 + r ** 2 - 2 * u * r * np.cos(th))) ** 0.5) * (((u ** 2 + r ** 2 - 2 * u * r * np.cos(th) + 4)) ** 0.5)))
    
    def function(variables,gamma):
        [pstar,u,times]=variables
        A = np.zeros(len(times))
        for i, time in enumerate(times):
            bottom = np.pi * pstar[i]**2
            top_res, top_err = integrate.dblquad(pre_top, 0, (2 * np.pi), 0, pstar[i], args=(u[i], pstar[i], gamma))
            A[i] = (top_res / bottom)
        return A
    return function(variables,gamma)*(1+(delf+sev))





#These are our input parameters
m1 = 2.4
m2 = 0.5
rstar = 0.95
ecc = 0.1
phi = 0.001
omega = 1
alpha = -1.3
gamma = 0.4
g = 0.4
Periods = 14.1
ti = 12.5


#Lastly, we plot the model over the data
final=(Light_curve_final(m1, m2, rstar, ecc, phi, omega, alpha, gamma, g, Periods, ti));
timed=np.linspace(0*86400,14*86400,2000)
plt.xlabel("Time (Days)",fontweight='bold',fontsize=14)
plt.ylabel("Relative Flux",fontweight='bold',fontsize=14)
plt.plot(timed/86400, final, "#FF6C0C")
#plt.xlim(3.1,3.3)
plt.show()
# %%