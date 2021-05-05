import numpy as np
import orbittools.constants as oc
from orbittools.kepler.kepler_eqn import *

mu_earth = oc.planets['earth']['gravitational_constant']


def radius_velocity_to_classical_elements(radius, velocity, mu=mu_earth):
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


'''
r = np.array([6524.834, 6862.875, 6448.296])
v = np.array([4.901327, 5.533756, -1.976341])
print(radius_velocity_to_elements(r, v))
'''


def rotx(theta):
    c = np.cos(theta)
    s = np.sin(theta)
    return np.array([[1, 0, 0],
                     [0, c, -s],
                     [0, s, c]])


def roty(theta):
    c = np.cos(theta)
    s = np.sin(theta)
    return np.array([[c, 0, s],
                     [0, 1, 0],
                     [-s, 0, c]])


def rotz(theta):
    c = np.cos(theta)
    s = np.sin(theta)
    return np.array([[c, -s, 0],
                     [s, c, 0],
                     [0, 0, 1]])


def orbit_to_3D_transform(inclination, right_ascension, arg_of_perigee, mu=mu_earth):
    return rotz(right_ascension) @ rotx(inclination) @ rotz(arg_of_perigee)


def classical_elements_to_radius_velocity(elements, mu=mu_earth):
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

    rot = orbit_to_3D_transform(i, omega, w)
    radius = rot @ radius_orbital
    velocity = rot @ velocity_orbital
    return radius, velocity


'''
elements = {
    'semilatus_rectum': 11067.79,
    'semimajor_axis': None,
    'eccentricity': 0.83285,
    'inclination': 87.87,
    'right_ascension': 227.89,
    'true_anomaly': 92.335,
    'arg_of_perigee': 53.38,
    'true_lon_of_perigee': None,
    'arg_of_latitude': None,
    'true_longitude': None,
    'noncircular_equatorial': False,
    'circular_inclined': False,
    'circular_equatorial': False}
print(classical_elements_to_radius_velocity(elements))
'''


def kepler_orbit_propagate(radius, velocity, duration, mu=mu_earth):
    elements = radius_velocity_to_classical_elements(radius, velocity)
    elements_next = elements
    print(elements)
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

    radius_next, velocity_next = classical_elements_to_radius_velocity(
        elements_next)
    return radius_next, velocity_next


'''
r = np.array([1131.34, -2282.343, 6672.423])
v = np.array([-5.64305, 4.30333, 2.42879])
dt = 40 * 60
print(kepler_orbit_propagate(r, v, dt))
'''


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
