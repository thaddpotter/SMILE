import ctypes as ct
import numpy as np
import matplotlib
matplotlib.use('TkAgg')

import matplotlib.pyplot as plt
import os


#Get function
smilepath = os.getcwd()
lib = ct.CDLL(os.path.join(smilepath,'smile.so'))


#Input paramters
m1 = ct.c_float(2.4)
m2 = ct.c_float(0.5)
rstar = ct.c_float(0.95)
ecc = ct.c_float(0.1)
phi = ct.c_float(0.001)
omega = ct.c_float(1)
alpha = ct.c_float(-1.3)
gamma = ct.c_float(0.4)
g = ct.c_float(0.4)
Periods = ct.c_float(14.1)
ti = ct.c_float(12.5)

sday = ct.c_double(86400)

#Array Buffers
t_array = (ct.c_double * 2048)()
curve = (ct.c_double * 2048)()
times = np.zeros(2048)


#Initialize t_array, and fill curve array to make sure it gets written to
for i in range(2048):
    t_array[i] = i * (14 * 86400/(2048-1))
    curve[i] = i * (14 * 86400/(2048-1))
    times[i] = i * (14 /(2048-1))

print(curve)



a = lib.gen_light_curve(t_array,curve,m1,m2,rstar,ecc,phi,omega,alpha,gamma,g,Periods,ti)

np.ctypeslib.as_array(t_array)
#np.ctypeslib.as_array(curve)

plt.xlabel("Time (Days)",fontweight='bold',fontsize=14)
plt.ylabel("Relative Flux",fontweight='bold',fontsize=14)
plt.plot(times, curve, "#FF6C0C")
#plt.xlim(3.1,3.3)
plt.show()