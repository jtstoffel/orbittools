import numpy as np
from numpy import pi
import orbittools.constants as oc

mu_earth = oc.planets['earth']['gravitational_constant']


def kepler_eqn_elliptic(mean_anomaly, eccentricity):
    if -pi < mean_anomaly < 0 or mean_anomaly > pi:
        eccentric_anomaly_prev = mean_anomaly - eccentricity
    else:
        eccentric_anomaly_prev = mean_anomaly + eccentricity
    while True:
        eccentric_anomaly = eccentric_anomaly_prev + \
            ((mean_anomaly - eccentric_anomaly_prev +
             (eccentricity * np.sin(eccentric_anomaly_prev))) / (1 - (eccentricity * np.cos(eccentric_anomaly_prev))))
        if np.abs(eccentric_anomaly - eccentric_anomaly_prev) < 10e-8:
            break
        eccentric_anomaly_prev = eccentric_anomaly
    return eccentric_anomaly


def kepler_eqn_parabolic(duration, semilatus_rectum, mu=mu_earth):
    mean_angular_vel = 2 * np.sqrt(mu / (semilatus_rectum**3))
    s = 0.5 * np.arctan(1/(3 * mean_angular_vel * duration / 2))
    w = np.arctan((np.tan(s))**(1/3))
    parabolic_anomaly = 2 / np.tan(2 * w)
    return parabolic_anomaly


def kepler_eqn_hyperbolic(mean_anomaly, eccentricity):
    if eccentricity < 1.6:
        if -pi < mean_anomaly < 0 or mean_anomaly > pi:
            hyperbolic_anomaly_prev = mean_anomaly - eccentricity
        else:
            hyperbolic_anomaly_prev = mean_anomaly + eccentricity
    else:
        if eccentricity < 3.6 and np.abs(mean_anomaly) > pi:
            hyperbolic_anomaly_prev = mean_anomaly - \
                np.sign(mean_anomaly) * eccentricity
        else:
            hyperbolic_anaomaly_prev = mean_anomaly / (eccentricity - 1)
    while True:
        hyperbolic_anomaly = hyperbolic_anomaly_prev + \
            ((mean_anomaly + hyperbolic_anomaly_prev -
             (eccentricity * np.sinh(hyperbolic_anomaly_prev))) / (-1 + (eccentricity * np.cosh(hyperbolic_anomaly_prev))))
        if np.abs(hyperbolic_anomaly - hyperbolic_anomaly_prev) < 10e-8:
            break
        hyperbolic_anomaly_prev = hyperbolic_anomaly
    return hyperbolic_anomaly


def true_to_anomaly(true_anomaly, eccentricity):
    if eccentricity < 1:
        eccentric_anomaly = np.arctan2((np.sin(true_anomaly) * np.sqrt(1 - eccentricity**2)) / (1 + eccentricity * np.cos(
            true_anomaly)), (eccentricity + np.cos(true_anomaly)) / (1 + eccentricity * np.cos(true_anomaly)))
        return eccentric_anomaly
    if eccentricity > 1:
        hyperbolic_anomaly = np.arctanh(((np.sin(true_anomaly) * np.sqrt(eccentricity**2 - 1)) / (1 + eccentricity * np.cos(
            true_anomaly))) / ((eccentricity + np.cos(true_anomaly)) / (1 + eccentricity * np.cos(true_anomaly))))
        return hyperbolic_anomaly
    else:
        parabolic_anomaly = np.tan(true_anomaly / 2)
        return parabolic_anomaly


def anomaly_to_true(anomaly, eccentricity):
    if eccentricity < 1:
        true_anomaly = np.arctan2((np.sin(
            anomaly) * np.sqrt(1 - eccentricity**2)) / (1 - eccentricity * np.cos(anomaly)), (np.cos(anomaly) - eccentricity) / (1 - eccentricity * np.cos(anomaly)))
        return true_anomaly
    if eccentricity > 1:
        hyperbolic_anomaly = np.arctan2((-np.sinh(anomaly) * np.sqrt(eccentricity**2 - 1)) / (
            1 - eccentricity * np.cosh(anomaly)), (np.cosh(anomaly) - eccentricity) / (1 - eccentricity * np.cosh(anomaly)))
        return true_anomaly
    else:
        print("parablic eccentricity (e = 1), use parabolic_anomaly_to_true() instead")
        return None


def parabolic_anomaly_to_true(anomaly, semilatus_rectum, radius):
    return np.arctan2(semilatus_rectum * anomaly / radius, (semilatus_rectum - radius) / radius)
