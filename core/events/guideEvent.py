import numpy as np
from dataclasses import dataclass

@dataclass
class guideEvent:
    accCmd:np.array
    
    name:str = "guide"