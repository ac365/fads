import math
import numpy as np
 
def rodriguesRotation(vector,axis,angle):
    if np.linalg.norm(axis) != 1:
        axis /= np.linalg.norm(axis)

    rotVec  = vector*math.cos(angle)
    rotVec += np.cross(axis,vector)*math.sin(angle)
    rotVec += axis*np.dot(axis,vector)*(1-math.cos(angle))

    return rotVec