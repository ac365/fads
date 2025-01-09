import numpy as np
import matplotlib.pyplot as plt

from ifs.earth.ussa1976 import USSA1976

ctr      = 0
temp     = np.empty([84999,1])
pressure = np.empty([84999,1])
density  = np.empty([84999,1])
atmos    = USSA1976()

for h in range(1,85000):
    T,p,rho,a = atmos.getAtmosParams(h)
    temp[ctr] = T
    pressure[ctr] = p
    density[ctr]  = rho

    ctr += 1

plt.figure()
plt.plot(temp-273.15,range(1,85000))
plt.title('Altitude vs. Temperature')
plt.ylabel('Altitude (m)')
plt.xlabel('Temperature (C)')
plt.grid(visible=True, which='major')

plt.figure()
plt.plot(pressure/1e3,range(1,85000))
plt.title('Altitude vs. Pressure')
plt.ylabel('Altitude (m)')
plt.xlabel('Pressure (kPa)')
plt.grid(visible=True, which='major')

plt.figure()
plt.plot(density,range(1,85000))
plt.title('Altitude vs. Density (kPa)')
plt.ylabel('Altitude (m)')
plt.xlabel('Density (kg/m^3)')
plt.grid(visible=True, which='major')
