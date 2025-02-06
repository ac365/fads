import math
import numpy as np
import matplotlib.pyplot as plt

from ifs.eom.aero.fleemanAero import Aero
from ifs.eom.atmos.ussa1976   import USSA1976

#body params, from Fleeman pg. 591
body = {}
body["bodyLen"]        = 3.65506 #143.9in
body["bodyDiam"]       = 0.2032  #8in
body["noseFineness"]   = 2.4
body["bodyWidth"]      = 0.2032  #8in, circle cross-section
body["bodyHeight"]     = 0.2032  #8in, circle cross-section
body["Ae"]             = 0.007238695 #11.22sq in 
body["inertMass"]      = 166.4684    #367lb
body["propellantMass"] = 0

#init aero and atmosphere
aero      = Aero(body)
atmos     = USSA1976()
T,p,rho,a = atmos.getAtmosParams(6096)
sRef      = math.pi/4*body["bodyWidth"]*body["bodyHeight"]

#zero lift drag unit test
machBkpts  = np.arange(0.1,5,0.1)
Cd0        = np.empty(np.size(machBkpts))

ctr = 0
for mach in machBkpts:
    qBar       = 0.5*rho*(mach*a)*(mach*a)
    
    aeroForces = aero.computeForces(phi=0, alpha=0, mach=mach, qBar=qBar, 
                       powerOn=False)
    aeroCoeff  = aeroForces/(qBar*sRef)

    Cn       = -aeroCoeff[1]/math.cos(0)
    Cd0[ctr] = -(aeroCoeff[0] + Cn*math.sin(0))/math.cos(0)
        
    ctr += 1

plt.figure()
plt.plot(machBkpts,Cd0)
plt.title('Cd0 vs. Mach')
plt.ylabel('Cd0')
plt.xlabel('Mach')
plt.grid(visible=True, which='major')


#normal force unit test
alphaBkpts = np.arange(0,90,0.25)*math.pi/180
Cn         = np.empty(np.size(alphaBkpts))

ctr = 0
for alpha in alphaBkpts:
    mach = 0.8
    qBar = 0.5*rho*(mach*a)*(mach*a)
    
    aeroForces = aero.computeForces(phi=0, alpha=alpha, mach=mach, 
                                    qBar=qBar, powerOn=False)
    aeroCoeff  = aeroForces/(qBar*sRef)

    Cn[ctr] = -aeroCoeff[1]/math.cos(0)

    ctr += 1

plt.figure()
plt.plot(alphaBkpts*180/math.pi,Cn)
plt.title('Cn vs. Angle of Attack')
plt.ylabel('Cn')
plt.xlabel('Angle of Attack [deg]')
plt.grid(visible=True, which='major')
