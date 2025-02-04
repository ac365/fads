import numpy as np
from dataclasses import dataclass

@dataclass
class apEvent:
    angAcc:np.array = np.zeros(3,)
    name:str = "ap"