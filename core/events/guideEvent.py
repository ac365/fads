import numpy as np
from dataclasses import dataclass

@dataclass
class guideEvent:
    accCmd:np.array = np.zeros(3,)
    name:str = "guide"