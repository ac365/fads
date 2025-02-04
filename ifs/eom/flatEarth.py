import math
import numpy as np

from ifs.eom.aero.fleemanAero import Aero
from ifs.eom.atmos.ussa1976   import USSA1976
from core.components.eomBase  import eomBase
from utils.utils              import eulerToQuaternion
class EOM(eomBase):
    def __init__(self, initEuler:np.array, initPos:np.array,
                 initVel:np.array, initOmega:np.array, body:dict):
        #initial conditions
        self._quat   = eulerToQuaternion(initEuler)
        self._e2b    = self._calculateEarthToBody()
        self._posNed = initPos
        self._velNed = np.linalg.inv(self._e2b)*np.array([initVel,0,0])
        self._omega  = initOmega

        #mass
        self._mass = body["inertMass"] + body["propellantMass"]

        #objects
        self._atmos = USSA1976()
        self._aero  = Aero(body)

    def _updateQuaternion(self, dt:np.double):
        p = self._omega[0]
        q = self._omega[1]
        r = self._omega[2]

        q0 = self._quat[0]
        q1 = self._quat[1]
        q2 = self._quat[2]
        q3 = self._quat[3]

        k = 0.5
        l = 1 - (q0*q0 + q1*q1 + q2*q2 + q3*q3)
        m = np.matrix([0, -p, -q, -r], [p, 0, r, -q],
                      [q, -r, 0, p], [r, q, -p, 0])
        
        self._quat += (0.5*np.matmul(m,self._quat) + k*l*self._quat)*dt

    def _calculateEarthToBody(self):
        q0 = self._quat[0]
        q1 = self._quat[1]
        q2 = self._quat[2]
        q3 = self._quat[3]

        e2b = np.zeros([3,3])
        e2b[0][0] = q0*q0 + q1*q1 - q2*q2 - q3*q3
        e2b[0][1] = 2*(q1*q2 + q0*q3)
        e2b[0][2] = 2*(q1*q3 - q0*q2)
        e2b[1][0] = 2*(q1*q2 - q0*q3)
        e2b[1][1] = q0*q0 - q1*q1 + q2*q2 - q3*q3
        e2b[1][2] = 2*(q2*q3 + q0*q1)
        e2b[2][0] = 2*(q1*q3 + q0*q2)
        e2b[2][1] = 2*(q2*q3 - q0*q1)
        e2b[2][2] = q0*q0 - q1*q1 - q2*q2 + q3*q3

        self._e2b = e2b
    
    def _calculateBodyForces(self):
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
        aoa   = math.atan(-self._vel[2]/self.vel[0])
        phi   = math.atan2(2*(q2*q3 + q0*q1),
                          (q0*q0 - q1*q1 - q2*q2 + q3*q3))
        
        #TODO: powerOn should be dependent on propulsion...
        aeroForces = self._aero.computeForces(phi, aoa, mach, qBar,
                                              powerOn=False) 
        
        return aeroForces + thrust

    def _eulerEquations(self, omegaDot:np.array, dt:np.double): 
        #NOTE: this is 5dof... angular acceleration from psuedo-autopilot
        self._omega += omegaDot*dt
            #TODO: should dt be an input or member var...?

    def _newtonsEquations(self, bdyForces:np.array):
        gBdy = np.matmul(self._e2b,np.array([0,0,9.81]))
        fx = bdyForces[0]
        fy = bdyForces[1]
        fz = bdyForces[2]

        velBdy = self._e2b*self._velNed
        u = velBdy[0]
        v = velBdy[1]
        w = velBdy[2]

        p = self._omega[0]
        q = self._omega[1]
        r = self._omega[2]

        ax = r*v - q*w + fx/self._mass + gBdy[0]
        ay = p*w - r*u + fy/self._mass + gBdy[1]
        az = q*u - p*v + fz/self._mass + gBdy[2]
        accBdy  = np.array([ax,ay,az])

        return accBdy 

    def step(self, omegaDot:np.array, dt:np.double):
        bdyForces = self._calculateBodyForces()
        self._eulerEquations(omegaDot,dt)
        accBdy = self._newtonsEquations(bdyForces)
        
        #update states
        self._updateQuaternion()
        self._calculateEarthToBody()
        accNed = np.linalg.inv(self._e2b)*accBdy
        self._velNed += accNed*dt
        self._posNed += self._velNed*dt

        #update event
        self.event.quat      = self._quat
        self.event.angRate   = self._omega
        self.event.angAcc    = omegaDot
        self.event.posNed    = self._posNed
        self.event.velNed    = self._velNed
        self.event.accNed    = accNed
        self.event.e2b       = self._e2b
        self.event.mass      = self._mass
        self.event.bdyForces = bdyForces

        return self.event