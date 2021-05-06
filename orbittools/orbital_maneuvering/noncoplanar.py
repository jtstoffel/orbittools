import numpy as np
from orbittools.orbital_maneuvering.coplanar import vis_viva
import orbittools.constants as oc

radius_earth = oc.planets['earth']['radius']
mu_earth = oc.planets['earth']['gravitational_constant']


class CircularOrbit3D:
    def __init__(self, radius, inclination, right_ascension, mu=mu_earth):
        self.radius = radius
        self.vel = vis_viva(radius, radius)
        self.inclination = np.deg2rad(inclination)
        self.right_ascension = np.deg2rad(right_ascension)

    def inclination_change(self, delta_inc, flight_path_angle=0):
        return 2 * self.vel * np.cos(np.deg2rad(flight_path_angle)) * np.sin(np.deg2rad(delta_inc)/2)

    def ascending_node_change(self, delta_asc_node):
        return 2 * self.vel * np.sin(np.arccos((np.cos(self.inclination))**2 + ((np.sin(self.inclination))**2 * np.cos(np.deg2rad(delta_asc_node)))) / 2)

    def asc_node_inc_change(self, delta_inc, delta_asc_node):
        return 2 * self.vel * np.sin((np.arccos((np.cos(self.inclination) * np.cos(self.inclination + np.deg2rad(delta_inc))) + (np.sin(self.inclination) * np.sin(self.inclination + np.deg2rad(delta_inc)) * np.cos(np.deg2rad(delta_asc_node))))) / 2)
