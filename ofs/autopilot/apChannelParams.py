import numpy as np
from dataclasses import dataclass

@dataclass
class autopilotLateral:
    tau:np.double = 0.1
    xi:np.double  = 0.71
    w:np.double   = 16.05
    wz:np.double  = 20.5

@dataclass
class autopilotRoll:
    xi:np.double = 0.71
    w:np.double  = 7