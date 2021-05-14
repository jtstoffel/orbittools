from dataclasses import dataclass
from orbittools.constants import planets


class Body:
    def __init__(self, name):
        self.name = name
        try:
            self.data = planets[name]
            print('Loaded data for ' + name)
        except KeyError:
            print('Unable to load planet data for: ' + name)


@dataclass
class Satellite:
    mass: float
    orbit: Orbit
    body: Body
    reference_frame: str


@dataclass
class EarthSystem:
    pass
