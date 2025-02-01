import math
import numpy as np
 
def rodriguesRotation(vector,axis,angle):
    if np.linalg.norm(axis) != 1:
        axis /= np.linalg.norm(axis)

    rotVec  = vector*math.cos(angle)
    rotVec += np.cross(axis,vector)*math.sin(angle)
    rotVec += axis*np.dot(axis,vector)*(1-math.cos(angle))

    return rotVec

def eulerToQuaternion(euler:np.array):
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
