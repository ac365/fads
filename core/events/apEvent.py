import numpy as np
from dataclasses import dataclass

@dataclass
class apEvent:
    angAcc:np.array
    name:str = "ap"