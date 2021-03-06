import numpy as np
from orbittools.system.coordinates import *


def test_ECI_to_elements():
    r = np.array([6524.834, 6862.875, 6448.296])
    v = np.array([4.901327, 5.533756, -1.976341])
    answer = {'semilatus_rectum': 11067.798350991809, 'semimajor_axis': 36127.337763974785, 'eccentricity': 0.8328533990836885, 'inclination': 87.86912617702644, 'right_ascension': 227.8982603572737, 'true_anomaly': 92.3351567104033,
              'arg_of_perigee': 53.38493067019384, 'true_lon_of_perigee': 247.8064482284738, 'arg_of_latitude': 145.72008738059714, 'true_longitude': 55.282707981472676, 'noncircular_equatorial': False, 'circular_inclined': False, 'circular_equatorial': False}
    assert ECI_to_classical_elements(r, v) == answer


def test_classical_elements_to_ECI():
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
    r_answer = np.array([6525.36812099, 6861.5318349, 6449.11861416])
    v_answer = np.array([4.90227864,  5.53313957, -1.9757101])

    r, v = classical_elements_to_ECI(elements)
    assert np.allclose(r, r_answer)
    assert np.allclose(v, v_answer)
