import numpy as np
import math

class Aero:
    def __init__(self, noseFineness, bodyLen, bodyDiam, sRef,
                 bodyWidth, bodyHeight):
        self.sRef         = sRef
        self.bodyLen      = bodyLen
        self.bodyDiam     = bodyDiam
        self.noseFineness = noseFineness
        self.bodyWidth    = bodyWidth
        self.bodyHeight   = bodyHeight
    
    def computeForces(self, phi:np.double, alpha:np.double,
                      mach:np.double, qBar:np.double, powerOn:bool,
                      Ae:np.double, rho:np.double, velMag:np.double):
        #axial
        if mach >= 1.0:
            baseDrag = 0.25/mach
        else:
            baseDrag = (0.12 + 0.13*mach*mach)
        if powerOn:
            baseDrag *= (1-Ae/self.sRef)

        if mach >= 0.8:
            waveDrag = ((1.586 + 1.834/(mach*mach))*
                        math.pow(math.atan(0.5/self.noseFineness),1.69))
        else:
            waveDrag = 0.0
        
        skinDrag = (0.053*(self.bodyLen/self.bodyDiam)*
                    math.pow(mach/(qBar*self.bodyLen)))

        Cd0 = waveDrag + baseDrag + skinDrag
        
        #normal
        Cn  = self.bodyWidth/self.bodyHeight*math.cos(phi)**2
        Cn += self.bodyHeight/self.bodyWidth*math.sin(phi)**2
        Cn *= (abs(math.sin(2*alpha)*math.cos(alpha/2)) + 
               1.3*self.bodyLen/self.bodyDiam*math.sin(alpha)**2)
        
        #body
        Cx = -Cd0*math.cos(alpha)
        Cy = -Cn*math.cos(phi)
        Cz = -Cn*math.sin(phi)

        aeroForces  = np.array([Cx,Cy,Cz])
        aeroForces *= qBar*self.sRef

        return aeroForces
