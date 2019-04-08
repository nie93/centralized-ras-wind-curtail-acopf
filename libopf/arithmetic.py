from numpy import pi


__author__ = "Zhijie Nie"
__copyright__ = "Copyright (c) 2019, Zhijie Nie, Shyam Gopal, Washington State University"
__credits__ = []
__version__ = "0.1"


def deg2rad(d):
    return d / 180 * pi

def rad2deg(r):
    return r / pi * 180    