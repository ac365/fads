import numpy as np
from dataclasses import dataclass

@dataclass
class eomEvent:
    quat:np.array      = np.zeros(4,)
    angRate:np.array   = np.zeros(3,)
    angAcc:np.array    = np.zeros(3,)
    posNed:np.array    = np.zeros(3,)
    velNed:np.array    = np.zeros(3,)
    accNed:np.array    = np.zeros(3,)
    e2b:np.ndarray     = np.zeros([3,3])
    mass:np.double     = 0.
    bdyForces:np.array = np.zeros(3,)

    name:str = "eom"