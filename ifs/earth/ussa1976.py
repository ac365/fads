import math

class USSA1976:
    def __init__(self):
        self._R     = 287.053
        self._g     = 9.80665
        self._gamma = 1.401

    def getAtmosParams(self,height):
        if height <= 11e3:
            lapseRate = -0.0065
            refHeight = 0.
            refTemp   = 15 + 273.15
            refPressure = 101325
        elif height <= 20e3:
            lapseRate = 0.
            refHeight = 11e3
            refTemp   = -56.5 + 273.15
            refPressure = 22632.1
        elif height <= 32e3:
            lapseRate = 0.001
            refHeight = 20e3
            refTemp   = -56.5 + 273.15
            refPressure = 5474.88
        elif height <= 47e3:
            lapseRate = 0.0028
            refHeight = 32e3
            refTemp   = -44.5 + 273.15
            refPressure = 868.018 
        elif height <= 51e3:
            lapseRate = 0
            refHeight = 47e3
            refTemp   = -2.5 + 273.15
            refPressure = 110.906 
        elif height <= 71e3:
            lapseRate = -0.0028
            refHeight = 51e3
            refTemp   = -2.5 + 273.15
            refPressure = 66.9388
        elif height <= 85e3:
            lapseRate   = -0.002
            refHeight   = 71e3
            refTemp     = -58.5 + 273.15
            refPressure = 3.95641
    
        #atmospheric params
        temp = refTemp + (height - refHeight)*lapseRate
        if lapseRate != 0:
            pressure  = 1 + (height - refHeight)*lapseRate/refTemp
            pressure  = math.pow(pressure,-self._g/(self._R*lapseRate))
            pressure *= refPressure
        else:
            exp      = -self._g/(self._R*refTemp)*(height - refHeight)
            pressure = refPressure*math.pow(math.e,exp)
        density = pressure/(self._R*temp)

        #speed of sound
        a = math.sqrt(self._gamma*self._R*temp)

        return temp,pressure,density,a