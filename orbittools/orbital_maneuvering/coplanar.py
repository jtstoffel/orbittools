import numpy as np


def vis_viva(radius, semimajor_axis, mu):
    return np.sqrt((2*mu/radius) - (mu/semimajor_axis))


class Orbit2D:
    def __init__(self, semimajor_axis, eccentricity, mu):
        self.semimajor_axis = semimajor_axis
        self.eccentricity = eccentricity
        self.mu = mu

        self.determine_orbit_type()

        if self.elliptic:
            self.radius_p = semimajor_axis * (1 - eccentricity)
            self.radius_a = semimajor_axis * (1 + eccentricity)
            self.vel_p = vis_viva(radius_p, semimajor_axis, mu)
            self.vel_a = vis_viva(radius_a, semimajor_axis, mu)

        if self.circular:
            self.radius = semimajor_axis
            self.vel = vis_viva(self.radius, semimajor_axis, mu)

    def determine_orbit_type(self):
        if self.eccentricity < 0:
            raise ValueError
        if self.eccentricity == 0:
            self.circular = True
        if self.eccentricity == 1:
            self.parabolic = True
        if self.eccentricity > 1:
            self.hyperbolic = True
        else:
            self.elliptic = True


class CoplanarCoaxialTransfer:
    def __init__(self, orbit_i, orbit_f):
        self.orbit_i = orbit_i
        self.orbit_f = orbit_f

    def hohmann_transfer(self):
        pass

    def bielliptic_transfer(self, radius_trans):
        pass

    def one_tangent_burn(self, ):
        pass


class CoplanarPhasing:
    pass
