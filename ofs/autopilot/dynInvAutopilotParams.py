import numpy as np
from dataclasses import dataclass

@dataclass
class apLateral:
    tau:np.double = 0.1
    xi:np.double  = 0.71
    w:np.double   = 16.05
    wz:np.double  = 20.5

@dataclass
class apRoll:
    xi:np.double = 0.71
    w:np.double  = 7