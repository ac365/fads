class Autopilot():
    def __init__(self, ap:dict):
        self.tauRoll   = ap["tauRoll"]
        self.tauPitch  = ap["tauPitch"]
        self.gLimit    = ap["gLimit"]
        self.maxAngAcc = ap["maxAngAcc"]

    #TODO: autopilot logic should be done here instead of in EOM