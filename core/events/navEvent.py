import numpy as np
from dataclasses import dataclass

@dataclass
class navEvent:
    quat:np.array
    angRate:np.array
    angAcc:np.array
    posNed:np.array
    velNed:np.array
    accNed:np.array
    e2b:np.matrix
    mass:np.double
    bdyForces:np.array

    name:str = "nav"





