import numpy as np
import orbittools.constants as oc
from orbittools.propagators.kepler_eqn import *
from orbittools.math_tools import *
from dataclasses import dataclass
from orbittools.system.coordinates import PQW_to_ECI, ECI_to_classical_elements, classical_elements_to_ECI


mu_earth = oc.planets['earth']['gravitational_parameter']


def kepler_orbit_propagate(radius, velocity, duration, mu=mu_earth):
    elements = ECI_to_classical_elements(radius, velocity)
    elements_next = elements
    e = elements['eccentricity']
    a = elements['semimajor_axis']
    p = elements['semilatus_rectum']
    if e != 0:
        v = np.deg2rad(elements['true_anomaly'])
        anomaly = true_to_anomaly(e, v)
    else:
        if elements['circular_inclined']:
            anomaly = np.deg2rad(elements['arg_of_latitude'])
        elif elements['circular_equatorial']:
            anomaly = np.deg2rad(elements['true_longitude'])

    # Propogating anomaly
    if e < 1:
        M0 = anomaly - e * np.sin(anomaly)
        n = np.sqrt(mu / (a**3))
        M = M0 + n * duration
        anomaly_next = kepler_eqn_elliptic(M, e)
    if e > 1:
        M0 = e * np.sinh(anomaly) - anomaly
        n = np.sqrt(mu / (-a**3))
        M = M0 + n * duration
        anomaly_next = kepler_eqn_hyperbolic(M, e)
    if e == 1:
        anomaly_next = kepler_eqn_parabolic(duration, p)

    if e != 0:
        v_next = anomaly_to_true(anomaly_next, e)
        elements_next['true_anomaly'] = np.rad2deg(v_next)
    else:
        if elements['circular_inclined']:
            elements_next['arg_of_latitude'] = np.rad2deg(anomaly_next)
        elif elements['circular_equatorial']:
            elements_next['true_longitude'] = np.rad2deg(anomaly_next)

    radius_next, velocity_next = classical_elements_to_ECI(
        elements_next)
    return radius_next, velocity_next


def find_time_of_flight(radius, radius_next, semilatus_rectum, mu=mu_earth):
    p = semilatus_rectum
    r0 = np.sqrt(np.dot(radius, radius))
    r = np.sqrt(np.dot(radius_next, radius_next))
    cosdv = np.dot(radius, radius_next) / (r0 * r)
    dv = np.arccos(cosdv)
    k = r0 * r * (1 - cosdv)
    l = r0 + r
    m = r0 * r * (1 + cosdv)
    a = m * k * p / (((2 * m - l**2) * p**2) + (2 * k * l * p) - k**2)
    f = 1 - ((r/p) * (1 - cosdv))
    g = r0 * r * np.sin(dv) / np.sqrt(mu * p)
    if a > 0:
        fdot = np.sqrt(mu/p) * np.tan(dv/2) * \
            (((1 - cosdv) / p) - (1/r0) - (1/r))
        cosdE = 1 - (r0 * (1 - f) / a)
        sindE = -r0 * r * fdot / np.sqrt(mu * a)
        dE = np.arctan2(sindE, cosdE)
        tof = g + np.sqrt(a**3 / mu) * (dE - sindE)
    if a < 0:
        coshdH = 1 + (f - 1) * (r0 / a)
        dH = np.arccosh(coshdH)
        tof = g + np.sqrt((-a)**3 / mu) * (np.sinh(dH) - dH)
    else:
        c = np.sqrt(r0**2 + r**2 - (2 * r0 * r * cosdv))
        s = (r0 + r + c)/2
        tof = (2/3) * np.sqrt(s**3/(2*mu)) * (1 - (((s - c) / s)**(3/2)))
    return tof
