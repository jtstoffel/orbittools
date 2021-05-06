from orbittools.orbital_maneuvering.noncoplanar import *

orbit = CircularOrbit3D(11480.649364225166, 55, 0)


def test_inclination_change():
    di = 15
    assert orbit.inclination_change(di) == 1.5382018364126486


def test_asc_node_change():
    dasc = 45
    assert orbit.ascending_node_change(dasc) == 3.694195175425932


def test_asc_node_inc_change():
    di = -15
    dasc = 45
    assert orbit.asc_node_inc_change(di, dasc) == 3.6159245484963183
