import math
import numpy as np

class EOM:
    def __init__(self, initEuler, initPos, initVel):
        #TODO: init properly
        self._quat = self._initQuaternion(initEuler)
        self._pos  = initPos
        self._vel  = initVel

    def _initQuaternion(self, euler):
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

    def updateQuaternion(self, omega, dt):
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
        b2e[1][0] = 2*(q1*q2 + q0*q3)
        b2e[1][1] = q0*q0 - q1*q1 + q2*q2 - q3*q3
        b2e[1][2] = 2*(q2*q3 + q0*q1)
        b2e[2][0] = 2*(q1*q3 + q0*q2)
        b2e[2][1] = 2*(q2*q3 - q0*q1)
        b2e[2][2] = q0*q0 - q1*q1 - q2*q2 + q3*q3

        return b2e