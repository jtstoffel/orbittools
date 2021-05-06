from orbittools.orbital_maneuvering.coplanar import *


def test_hohmann():
    orbit1 = Orbit2D(191.344 + radius_earth, 0)
    orbit2 = Orbit2D(35781.348 + radius_earth, 0)
    transfer = CoplanarCoaxialTransfer(orbit1, orbit2)
    answer = {'deltaV1': 2.4570367587766535, 'deltaV2': 1.4781863179340886, 'deltaV_total': 3.9352230767107423,
              'semimajor_axis_trans': 24364.485999999997, 'duration': 18924.171087666222}
    assert transfer.hohmann_transfer() == answer


def test_bielliptic():
    orbit1 = Orbit2D(191.34411 + radius_earth, 0)
    orbit2 = Orbit2D(376310 + radius_earth, 0)
    trans_radius = 503873 + radius_earth
    transfer = CoplanarCoaxialTransfer(orbit1, orbit2)
    answer = {'deltaV1': 3.156232460620802, 'deltaV2': 0.6773579575497282, 'deltaV3': -0.07046593618058172, 'deltaV_total': 3.9040563543511118,
              'semimajor_axis_trans1': 258410.31205500002, 'semimajor_axis_trans2': 446469.64, 'duration': 2138112.0337858316}
    assert transfer.bielliptic_transfer(trans_radius) == answer


def test_one_tangent_burn():
    orbit1 = Orbit2D(191.344 + radius_earth, 0)
    orbit2 = Orbit2D(35781.348 + radius_earth, 0)
    v = np.deg2rad(160)
    transfer = CoplanarCoaxialTransfer(orbit1, orbit2)
    answer = {'deltaV1': 2.5753941000853002, 'deltaV2': 2.1239359698163502, 'deltaV_total': 4.6993300699016505, 'semimajor_axis_trans':
              28634.249800512778, 'eccentricity_trans': 0.7705725120871735, 'burn_angle': 43.688796032297425, 'duration': 12446.730369829895}
    assert transfer.one_tangent_burn(v) == answer


def test_determine_orbit():
    orbit = Orbit2D(10000, 0)
    assert orbit.circular
    orbit = Orbit2D(10000, 0.3)
    assert orbit.elliptic
    orbit = Orbit2D(10000, 1.5)
    assert orbit.hyperbolic
    orbit = Orbit2D(10000, 1)
    assert orbit.parabolic
