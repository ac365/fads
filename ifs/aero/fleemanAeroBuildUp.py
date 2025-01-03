import numpy as np
import math

class Aero:
    def __init__(self, noseFineness, bodyLenFt, bodyDiamFt, sRef,
                 bodyWidth, bodyHeight):
        self.sRef         = sRef
        self.bodyLenFt    = bodyLenFt
        self.bodyDiamFt   = bodyDiamFt
        self.noseFineness = noseFineness
        self.bodyWidth    = bodyWidth
        self.bodyHeight   = bodyHeight

    def zeroLiftDrag(self, mach:np.double, qBar:np.double,
                     powerOn:bool, Ae:np.double):
        
        supersonic = mach >= 1.0
        if supersonic:
            baseDrag = 0.25/mach*(powerOn*(1-Ae/self.sRef))
        else:
            baseDrag = (0.12 + 0.13*mach*mach)*(powerOn*(1-Ae/self.sRef))
        waveDrag   = ((1.586 + 1.834/(mach*mach))*
            math.pow(math.atan(0.5/self.noseFineness),1.69))
        skinDrag = (0.053*(self.bodyLenFt/self.bodyDiamFt)*
            math.pow(mach/(qBar*self.bodyLenFt)))

        Cd0 = waveDrag + baseDrag + skinDrag
        return Cd0 
    
    def computeForces(self, rollAng:np.double, alpha:np.double,
                     mach:np.double, qBar:np.double, powerOn:bool,
                     Ae:np.double):
        Cd0 = self._zeroLiftDrag(mach,qBar,powerOn,Ae)

        Cn  = self.bodyWidth/self.bodyHeight*math.cos(rollAng)**2
        Cn += self.bodyHeight/self.bodyWidth*math.sin(rollAng)**2
        Cn *= (abs(math.sin(2*alpha)*math.cos(alpha/2)) + 
               1.3*self.bodyLenFt/self.bodyDiamFt*math.sin(alpha)**2)
        
        Cl = Cn*math.cos(alpha) - Cd0*math.sin(alpha)
        Cd = Cn*math.sin(alpha) + Cd0*math.cos(alpha)

        return Cl,Cd
