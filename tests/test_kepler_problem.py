import numpy as np
from orbittools.propagators.kepler_problem import *


def test_kepler_orbit_propagate():
    r = np.array([1131.34, -2282.343, 6672.423])
    v = np.array([-5.64305, 4.30333, 2.42879])
    dt = 40 * 60
    r_next_answer = np.array([-4191.18587486,  4348.12233743, -4005.749023])
    v_next_answer = np.array([3.72381066, -1.951891, -6.08036812])
    r_next, v_next = kepler_orbit_propagate(r, v, dt)
    assert np.allclose(r_next, r_next_answer)
    assert np.allclose(v_next, v_next_answer)
