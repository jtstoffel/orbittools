import numpy as np
import orbittools.constants as oc

radius_earth = oc.planets['earth']['radius']
mu_earth = oc.planets['earth']['gravitational_parameter']


def vis_viva(radius, semimajor_axis, mu=mu_earth):
    return np.sqrt((2*mu/radius) - (mu/semimajor_axis))


class Orbit2D:
    def __init__(self, semimajor_axis, eccentricity, mu=mu_earth):
        self.semimajor_axis = semimajor_axis
        self.eccentricity = eccentricity
        self.mu = mu
        self.circular = False
        self.elliptic = False
        self.hyperbolic = False
        self.parabolic = False
        self.determine_orbit_type()

        if self.elliptic or self.circular:
            self.radius_p = semimajor_axis * (1 - eccentricity)
            self.radius_a = semimajor_axis * (1 + eccentricity)
            self.vel_p = vis_viva(self.radius_p, semimajor_axis, mu)
            self.vel_a = vis_viva(self.radius_a, semimajor_axis, mu)

        if self.hyperbolic or self.parabolic:
            # Add quantities of interest for these orbit types
            pass

    def determine_orbit_type(self):
        if self.eccentricity < 0:
            raise ValueError
        elif self.eccentricity == 0:
            self.circular = True
        elif self.eccentricity == 1:
            self.parabolic = True
        elif self.eccentricity > 1:
            self.hyperbolic = True
        elif 0 < self.eccentricity < 1:
            self.elliptic = True


class CoplanarCoaxialTransfer:
    def __init__(self, orbit_i, orbit_f):
        self.orbit_i = orbit_i
        self.orbit_f = orbit_f
        assert orbit_i.mu == orbit_f.mu
        self.mu = orbit_i.mu
        self.radius_i = orbit_i.radius_p
        self.vel_i = orbit_i.vel_p
        self.radius_f = orbit_f.radius_a
        self.vel_f = self.orbit_f.vel_a

    def hohmann_transfer(self):
        semimajor_axis_trans = (self.radius_i + self.radius_f) / 2.0
        vel_i_trans = vis_viva(self.radius_i, semimajor_axis_trans, self.mu)
        vel_f_trans = vis_viva(self.radius_f, semimajor_axis_trans, self.mu)
        deltaV1 = vel_i_trans - self.vel_i
        deltaV2 = self.vel_f - vel_f_trans
        deltaV_total = np.abs(deltaV1) + np.abs(deltaV2)
        duration = np.pi * np.sqrt(semimajor_axis_trans**3 / self.mu)
        return {'deltaV1': deltaV1,
                'deltaV2': deltaV2,
                'deltaV_total': deltaV_total,
                'semimajor_axis_trans': semimajor_axis_trans,
                'duration': duration}

    def bielliptic_transfer(self, radius_trans):
        semimajor_axis_trans1 = (self.radius_i + radius_trans) / 2.0
        semimajor_axis_trans2 = (self.radius_f + radius_trans) / 2.0
        vel_i_trans1 = vis_viva(self.radius_i, semimajor_axis_trans1, self.mu)
        vel_f_trans1 = vis_viva(radius_trans, semimajor_axis_trans1, self.mu)
        vel_i_trans2 = vis_viva(radius_trans, semimajor_axis_trans2, self.mu)
        vel_f_trans2 = vis_viva(self.radius_f, semimajor_axis_trans2, self.mu)
        deltaV1 = vel_i_trans1 - self.vel_i
        deltaV2 = vel_i_trans2 - vel_f_trans1
        deltaV3 = self.vel_f - vel_f_trans2
        deltaV_total = np.abs(deltaV1) + np.abs(deltaV2) + np.abs(deltaV3)
        duration = (np.pi * np.sqrt(semimajor_axis_trans1**3 / self.mu)) + \
            (np.pi * np.sqrt(semimajor_axis_trans2**3 / self.mu))
        return {'deltaV1': deltaV1,
                'deltaV2': deltaV2,
                'deltaV3': deltaV3,
                'deltaV_total': deltaV_total,
                'semimajor_axis_trans1': semimajor_axis_trans1,
                'semimajor_axis_trans2': semimajor_axis_trans2,
                'duration': duration}

    def one_tangent_burn(self, true_anomaly_trans, at_periapse=True):
        radius_ratio = self.radius_i / self.radius_f
        if at_periapse:
            eccentricity_trans = (radius_ratio - 1) / \
                (np.cos(true_anomaly_trans) - radius_ratio)
            semimajor_axis_trans = self.radius_i / (1 - eccentricity_trans)
        else:
            eccentricity_trans = (radius_ratio - 1) / \
                (np.cos(true_anomaly_trans) + radius_ratio)
            semimajor_axis_trans = self.radius_i / (1 + eccentricity_trans)
        vel_trans_i = vis_viva(self.radius_i, semimajor_axis_trans, self.mu)
        vel_trans_f = vis_viva(self.radius_f, semimajor_axis_trans, self.mu)

        burn_angle = np.arctan((eccentricity_trans * np.sin(true_anomaly_trans)) /
                               (1 + eccentricity_trans * np.cos(true_anomaly_trans)))
        deltaV1 = vel_trans_i - self.vel_i
        deltaV2 = np.sqrt(vel_trans_f**2 + self.vel_f**2 -
                          (2 * self.vel_f * vel_trans_f * np.cos(burn_angle)))
        deltaV_total = np.abs(deltaV1) + np.abs(deltaV2)
        E = np.arccos((eccentricity_trans + np.cos(true_anomaly_trans)) /
                      (1 + (eccentricity_trans * np.cos(true_anomaly_trans))))
        duration = np.sqrt(semimajor_axis_trans**3 / self.mu) * \
            ((E - eccentricity_trans * np.sin(E)))
        return {
            'deltaV1': deltaV1,
            'deltaV2': deltaV2,
            'deltaV_total': deltaV_total,
            'semimajor_axis_trans': semimajor_axis_trans,
            'eccentricity_trans': eccentricity_trans,
            'burn_angle': np.rad2deg(burn_angle),
            'duration': duration
        }
