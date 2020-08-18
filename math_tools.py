import numpy as np
import constants


def classical_to_equinoctial(x):
    semimajor_axis = a = x[0]
    eccentricity = e = x[1]
    inclination = i = x[2]
    argument_of_perigee = w = x[3]
    right_ascension_of_ascending_node = raas = x[4]
    true_anomaly = v = x[5]
    p = a * (1 - e**2)
    f = e * np.cos(w + raas)
    g = e * np.sin(w + raas)
    h = np.tan(i / 2.0) * np.sin(raas)
    k = np.tan(i / 2.0) * np.cos(raas)
    L = raas + w + v
    return np.array([p, f, g, h, k, L])


def equinoctial_to_ECI(x):
    p = x[0]
    f = x[1]
    g = x[2]
    h = x[3]
    k = x[4]
    L = x[5]
    mu = constants.planets['earth']['gravitational_constant']
    alpha = h**2 - k**2
    beta = 1 + h**2 + k**2
    w = 1 + (f * np.cos(L)) + (g * np.sin(L))
    r = p/w
    x = (r / beta) * (np.cos(L) + (alpha * np.cos(L)) + (2 * h * k * np.sin(L)))
    y = (r / beta) * (np.sin(L) - (alpha * np.sin(L)) + (2 * h * k * np.cos(L)))
    z = (2 * r / beta) * ((h * np.sin(L)) - (k * np.cos(L)))
    vx = (-1 / beta) * np.sqrt(mu / p) * (np.sin(L) + (alpha * np.sin(L)
                                                       ) - (2 * h * k * np.cos(L)) + g - (2 * f * h * k) + (alpha * g))
    vy = (-1 / beta) * np.sqrt(mu / p) * (-np.cos(L) + (alpha * np.cos(L)
                                                        ) + (2 * h * k * np.cos(L)) - f + (2 * g * h * k) + (alpha * f))
    vz = (2 / beta) * np.sqrt(mu / p) * \
        ((h * np.cos(L)) + (k * np.sin(L)) + (f * h) + (g * k))
    return np.array([x, y, z , vx, vy, vz])


