import math
import numpy as np

from ifs.eom.aero.fleemanAero import Aero
from ifs.eom.atmos.ussa1976   import USSA1976

class EOM:
    def __init__(self, initEuler:np.array, initPos:np.array,
                 initVel:np.array, body:dict):
        #initial conditions
        self._quat = self._initQuaternion(initEuler)
        self._pos  = initPos
        self._vel  = initVel

        #mass
        self._mass = body["inertMass"] + body["propellantMass"]

        #objects
        self._atmos = USSA1976()
        self._aero  = Aero(body)

    def _initQuaternion(self, euler:np.array):
        psi   = euler[0]
        theta = euler[1]
        phi   = euler[2]

        q0 = (np.cos(psi/2)*np.cos(theta/2)*np.cos(phi/2) + 
              np.sin(psi/2)*np.sin(theta/2)*np.sin(phi/2))
        q1 = (np.cos(psi/2)*np.cos(theta/2)*np.sin(phi/2) - 
              np.sin(psi/2)*np.sin(theta/2)*np.cos(phi/2))
        q2 = (np.cos(psi/2)*np.sin(theta/2)*np.cos(phi/2) +
              np.sin(psi/2)*np.cos(theta/2)*np.sin(phi/2))
        q3 = (np.sin(psi/2)*np.cos(theta/2)*np.cos(phi/2) -
              np.cos(psi/2)*np.sin(theta/2)*np.sin(phi/2))
        
        return np.array([q0, q1, q2, q3])

    def updateQuaternion(self, omega:np.array, dt:np.double):
        p = omega[0]
        q = omega[1]
        r = omega[2]

        q0 = self._quat[0]
        q1 = self._quat[1]
        q2 = self._quat[2]
        q3 = self._quat[3]

        k = 0.5
        l = 1 - (q0*q0 + q1*q1 + q2*q2 + q3*q3)
        m = np.matrix([0, -p, -q, -r], [p, 0, r, -q],
                      [q, -r, 0, p], [r, q, -p, 0])
        
        self._quat += (0.5*np.matmul(m,self._quat) + k*l*self._quat)*dt

    def calculateBodyToEarth(self):
        q0 = self._quat[0]
        q1 = self._quat[1]
        q2 = self._quat[2]
        q3 = self._quat[3]

        b2e = np.zeros([3,3])
        b2e[0][0] = q0*q0 + q1*q1 - q2*q2 - q3*q3
        b2e[0][1] = 2*(q1*q2 + q0*q3)
        b2e[0][2] = 2*(q1*q3 - q0*q2)
        b2e[1][0] = 2*(q1*q2 - q0*q3)
        b2e[1][1] = q0*q0 - q1*q1 + q2*q2 - q3*q3
        b2e[1][2] = 2*(q2*q3 + q0*q1)
        b2e[2][0] = 2*(q1*q3 + q0*q2)
        b2e[2][1] = 2*(q2*q3 - q0*q1)
        b2e[2][2] = q0*q0 - q1*q1 - q2*q2 + q3*q3

        return b2e
    
    def calculateBodyForces(self):
        #quaternion
        q0 = self._quat[0]
        q1 = self._quat[1]
        q2 = self._quat[2]
        q3 = self._quat[3]

        #thrust
        thrust = np.zeros(3) #TODO: call to propulsion object here

        #aero forces
        rho, a = self._atmos.getAtmosParams(-self._pos[3])[2:]
        velMag = np.linalg.norm(self._vel)
        
        qBar  = 0.5*rho*velMag*velMag
        mach  = velMag/a
        aoa   = math.atan(self._vel[2]/self.vel[0])
        phi   = math.atan2(2*(q2*q3 + q0*q1),
                          (q0*q0 - q1*q1 - q2*q2 + q3*q3))
        
        #TODO: powerOn should be dependent on propulsion...
        aeroForces = self._aero.computeForces(phi, aoa, mach, qBar,
                                              powerOn=False) 
        
        return aeroForces + thrust

    def calculateBodyRates(self, accCmd:np.array): 
        #NOTE: this function is actually a pseudo-autopilot        
        
        #quaternion
        q0 = self._quat[0]
        q1 = self._quat[1]
        q2 = self._quat[2]
        q3 = self._quat[3]

        #aero params
        rho, a = self._atmos.getAtmosParams(-self._pos[2])[2:]
        velMag = np.linalg.norm(self._vel)
        
        qBar  = 0.5*rho*velMag*velMag
        mach  = velMag/a
        phi   = math.atan2(2*(q2*q3 + q0*q1),
                           (q0*q0 - q1*q1 - q2*q2 + q3*q3))

        #command should be perpendicular to velocity and g-limited
        velHat  = self._vel/np.linalg.norm(self._vel)
        accCmd -= np.dot(accCmd,velHat)*velHat
        cmdHat  = accCmd/np.linalg.norm(accCmd)
        if np.linalg.norm(accCmd) > 45:
            accCmd = cmdHat*45
            #TODO: g-limit shouldn't be hard-coded

        #convert acceleration command to aoa command
        aoaMin = 0
        aoaMax = math.radians(25)
        tol    = math.radians(0.001)

        hiVel   = np.array([math.cos(aoaMax),0,math.sin(aoaMax)])
        hiLift  = self._aero.computeForces(phi,aoaMax,mach,qBar,False)
        hiLift -= np.dot(hiLift,hiVel)*hiVel
        hiLift  = np.linalg.norm(hiLift)
        loVel   = np.array([math.cos(aoaMin),0,math.sin(aoaMin)])
        loLift  = self._aero.computeForces(phi,aoaMin,mach,qBar,False)
        loLift -= np.dot(loLift,loVel)*loVel
        loLift  = np.linalg.norm(loLift)
            #TODO: powerOn should depend on thrust object

        ctr = 1
        while ctr < 50 and (aoaMax - aoaMin) > tol:
            guess  = (aoaMin + aoaMax)/2
            gVel   = np.array([math.cos(guess),0,math.sin(guess)])
            gLift  = self._aero.computeForces(phi,guess,mach,qBar,False)
                #TODO: powerOn should depend on thrust object
            gLift -= np.dot(gLift,gVel)*gVel
            gLift  = np.linalg.norm(gLift)
            
            lErr = loLift/self._mass - np.linalg.norm(accCmd)
            gErr = gLift/self._mass - np.linalg.norm(accCmd) 

            if np.sign(gErr) == np.sign(lErr):
                aoaMin = guess
                loLift = gLift
            else:
                aoaMax = guess

            ctr += 1
        aoaDes = guess
        aoa    = math.atan(self._vel[2]/self._vel[0])
        aoaCmd = aoaDes - aoa

        #convert acceleration command to phi command
        phiDes = math.atan2(cmdHat[2],cmdHat[1])
        phiCmd = (phi-phiDes+math.pi)%(2*math.pi) - math.pi
        if phiCmd < -math.pi:
            phiCmd += 2*math.pi