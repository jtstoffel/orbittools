from orbittools.system.time import *


def test_julian_date():
    assert julian_date(1996, 10, 26, 14, 20, 0) == 2450383.097222222


'''
def test_get_gregorian_date():
    assert get_gregorian_date(2449877.3458762)

'''
