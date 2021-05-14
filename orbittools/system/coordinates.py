import numpy as np
import orbittools.constants as oc
from orbittools.math_tools import rotx, roty, rotz
from dataclasses import dataclass

mu_earth = oc.planets['earth']['gravitational_parameter']


@dataclass
class Orbit:
    semiparameter: float
    semimajor_axis: float
    eccentricity: float
    inclination: float


def SEZ_to_ECEF():
    pass


def PQW_to_ECI(inclination, right_ascension, arg_of_perigee, mu=mu_earth):
    return rotz(right_ascension) @ rotx(inclination) @ rotz(arg_of_perigee)


def PQW_to_RSW(true_anomaly):
    return rotz(-true_anomaly)


def RSW_to_PQW(true_anomaly):
    return rotz(true_anomaly)


def ECI_to_EQW(inclination, right_ascension):
    pass


def SEZ_to_BODY(roll, pitch, yaw):
    return rotz(-yaw) @ rotx(-roll) @ roty(-pitch)


def ECI_to_classical_elements(radius, velocity, mu=mu_earth):
    h = np.cross(radius, velocity)
    h_mag = np.sqrt(np.dot(h, h))
    n = np.cross(np.array([0, 0, 1]), h)
    n_mag = np.sqrt(np.dot(n, n))
    v_mag = np.sqrt(np.dot(velocity, velocity))
    r_mag = np.sqrt(np.dot(radius, radius))
    e = (((v_mag**2 - (mu/r_mag)) * radius) -
         (np.dot(radius, velocity) * velocity)) / mu
    e_mag = np.sqrt(np.dot(e, e))
    E = (v_mag**2 / 2) - (mu/r_mag)
    i = np.arccos(h[2] / h_mag)

    ce = False
    ci = False
    nce = False
    if e_mag == 0:
        if i == 0:
            ce = True
        else:
            ci = True
    elif i == 0:
        nce = True

    if e_mag != 1:
        a = -mu/(2*E)
        p = a * (1-e_mag**2)
    else:
        p = h_mag**2 / mu
        a = np.inf

    omega = np.arccos(n[0] / n_mag)
    if n[1] < 0:
        omega = np.deg2rad(360) - omega

    w = np.arccos(np.dot(n, e)/(n_mag * e_mag))
    if e[2] < 0:
        w = np.deg2rad(360) - w

    v = np.arccos(np.dot(e, radius)/(e_mag * r_mag))
    if np.dot(radius, velocity) < 0:
        v = np.deg2rad(360) - v

    w_true = np.arccos(e[0] / e_mag)
    if e[1] < 0:
        w_true = np.deg2rad(360) - w_true

    u = np.arccos(np.dot(n, radius) / (n_mag * r_mag))
    if radius[2] < 0:
        u = np.deg2rad(360) - u

    lon = np.arccos(radius[0] / r_mag)
    if radius[1] < 0:
        lon = np.deg2rad(360) - lon

    return {
        'semilatus_rectum': p,
        'semimajor_axis': a,
        'eccentricity': e_mag,
        'inclination': np.rad2deg(i),
        'right_ascension': np.rad2deg(omega),
        'true_anomaly': np.rad2deg(v),
        'arg_of_perigee': np.rad2deg(w),
        'true_lon_of_perigee': np.rad2deg(w_true),
        'arg_of_latitude': np.rad2deg(u),
        'true_longitude': np.rad2deg(lon),
        'noncircular_equatorial': nce,
        'circular_inclined': ci,
        'circular_equatorial': ce}


def classical_elements_to_ECI(elements, mu=mu_earth):
    e = elements['eccentricity']
    i = np.deg2rad(elements['inclination'])
    p = elements['semilatus_rectum']

    if elements['circular_equatorial']:
        w = 0
        omega = 0
        v = np.deg2rad(elements['true_longitude'])
    elif elements['circular_inclined']:
        w = 0
        omega = np.deg2rad(elements['right_ascension'])
        v = np.deg2rad(elements['arg_of_latitude'])
    elif elements['noncircular_equatorial']:
        w = np.deg2rad(elements['true_lon_of_perigee'])
        omega = 0
        v = np.deg2rad(elements['true_anomaly'])
    else:
        w = np.deg2rad(elements['arg_of_perigee'])
        omega = np.deg2rad(elements['right_ascension'])
        v = np.deg2rad(elements['true_anomaly'])

    radius_orbital = np.array(
        [p * np.cos(v) / (1 + e*np.cos(v)), p * np.sin(v) / (1 + e*np.cos(v)), 0.0])

    velocity_orbital = np.array(
        [-np.sqrt(mu/p) * np.sin(v), np.sqrt(mu/p) * (e + np.cos(v)), 0])

    rot = PQW_to_ECI(i, omega, w)
    radius = rot @ radius_orbital
    velocity = rot @ velocity_orbital
    return radius, velocity
