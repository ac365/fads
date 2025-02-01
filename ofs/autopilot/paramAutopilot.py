import math
import numpy as np

from ifs.eom.aero.fleemanAero import Aero
from ifs.eom.atmos.ussa1976   import USSA1976
from core.events.guideEvent   import guideEvent
from core.events.navEvent     import navEvent

class Autopilot():
    def __init__(self, body:dict):
        self._aero  = Aero(body)#TODO:should be trim aero model not truth
        self._atmos = USSA1976()

    def accCmdToIncidenceAngleCmd(self, ge:guideEvent, ne:navEvent):
        #quaternion
        q0 = ne.quat[0]
        q1 = ne.quat[1]
        q2 = ne.quat[2]
        q3 = ne.quat[3]

        #aero params
        rho, a = self._atmos.getAtmosParams(-ne.posNed[2])[2:]
        velMag = np.linalg.norm(ne.velNed)
        
        qBar  = 0.5*rho*velMag*velMag
        mach  = velMag/a
        phi   = math.atan2(2*(q2*q3 + q0*q1),
                           (q0*q0 - q1*q1 - q2*q2 + q3*q3))

        #command should be perpendicular to velocity and g-limited
        velBdy  = ne.e2b*ne.velNed
        velHat  = velBdy/np.linalg.norm(velBdy)
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
            
            lErr = loLift/ne.mass - np.linalg.norm(accCmd)
            gErr = gLift/ne.mass - np.linalg.norm(accCmd) 

            if np.sign(gErr) == np.sign(lErr):
                aoaMin = guess
                loLift = gLift
            else:
                aoaMax = guess

            ctr += 1
        aoaDes = guess
        aoa    = math.atan(velBdy[2]/velBdy[0])
        aoaCmd = aoaDes - aoa

        #convert acceleration command to phi command
        phiDes = math.atan2(cmdHat[2],cmdHat[1])
        phiCmd = (phi-phiDes+math.pi)%(2*math.pi) - math.pi
        if phiCmd < -math.pi:
            phiCmd += 2*math.pi

        return aoaCmd, phiCmd