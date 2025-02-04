from core.components.navBase import navBase
from core.events.eomEvent    import eomEvent

class TruthNav(navBase):
    def __init__(self):
        pass

    def step(self, ee:eomEvent):
        self.event.quat      = ee.quat
        self.event.angRate   = ee.angRate
        self.event.angAcc    = ee.angAcc
        self.event.posNed    = ee.posNed
        self.event.velNed    = ee.velNed
        self.event.accNed    = ee.accNed
        self.event.e2b       = ee.e2b
        self.event.mass      = ee.mass
        self.event.bdyForces = ee.bdyForces

        return self.event