import numpy as np
import orbittools.constants as oc

mu_earth = oc.planets['earth']['gravitational_constant']


def circular_coplanar_phasing(radius_target, radius_interceptor, phase_offset, k=1, mu=mu_earth):
    ang_vel_target = np.sqrt(mu / radius_target**3)
    ang_vel_interceptor = np.sqrt(mu / radius_interceptor**3)
