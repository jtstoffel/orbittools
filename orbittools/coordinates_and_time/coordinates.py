import numpy as np
import orbittools.constants as oc


def ecef_to_lat_long(r_ecef):
    r_eq = np.sqrt(r_ecef[0]**2 + r_ecef[1]**2)
    sina = r_ecef[1] / r_eq
    cosa = r_ecef[0] / r_eq
    lon = np.arctan2(sina, cosa)
