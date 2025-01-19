import numpy as np
import math

class Aero:
    def __init__(self, body:dict):
        self._bodyLen      = body["bodyLen"]
        self._bodyDiam     = body["bodyDiam"]
        self._noseFineness = body["noseFineness"]
        self._bodyWidth    = body["bodyWidth"]
        self._bodyHeight   = body["bodyHeight"]
        self._Ae           = body["Ae"]
        self._sRef         = math.pi/4*self._bodyWidth*self._bodyHeight
    
    def computeForces(self, phi:np.double, alpha:np.double,
                      mach:np.double, qBar:np.double, powerOn:bool):
        #zero lift drag
        if mach >= 1.0:
            baseDrag = 0.25/mach
        else:
            baseDrag = (0.12 + 0.13*mach*mach)
        if powerOn:
            baseDrag *= (1-self._Ae/self._sRef)

        if mach >= 0.8:
            waveDrag = ((1.586 + 1.834/(mach*mach))*
                        math.pow(math.atan(0.5/self._noseFineness),1.69))
        else:
            waveDrag = 0.0
        
        skinDrag = (0.053*(self._bodyLen/self._bodyDiam)*
                    math.pow(mach/(qBar*self._bodyLen),0.2))

        Cd0 = waveDrag + baseDrag + skinDrag
        
        #normal
        Cn  = self._bodyWidth/self._bodyHeight*math.cos(phi)**2
        Cn += self._bodyHeight/self._bodyWidth*math.sin(phi)**2
        Cn *= (abs(math.sin(2*alpha)*math.cos(alpha/2)) + 
               1.3*self._bodyLen/self._bodyDiam*math.sin(alpha)**2)
        
        #drag
        Cd = Cn*math.sin(alpha) + Cd0*math.cos(alpha)
        
        #body
        Cx = -Cd*math.cos(alpha)
        Cy = -Cn*math.cos(phi)
        Cz = -Cn*math.sin(phi) + Cd*math.sin(alpha)

        aeroForces  = np.array([Cx,Cy,Cz])
        aeroForces *= qBar*self._sRef

        return aeroForces
