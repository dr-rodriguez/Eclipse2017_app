# Function to calculate the distance to the moon

import numpy as np
from astropy.constants import R_earth, R_sun
from astropy import units as u
from astropy.coordinates import get_sun, get_moon, EarthLocation, AltAz


def moon_distance(loc, t2, t3):
    # Get the distance to the moon via the umbra diameter estimate

    if isinstance(loc, type(u'')) or isinstance(loc, type('')):
        loc = EarthLocation.of_address(loc)

    delta_t = t3-t2

    # Location of Sun & Moon
    sun = get_sun(t2)
    moon = get_moon(t2, loc)
    sun_altaz = sun.transform_to(AltAz(obstime=t2, location=loc))

    # Umbra size, from time and location
    s_umbra = 2 * np.pi * R_earth * np.cos(loc.latitude.to(u.rad)) / u.day * delta_t
    s_umbra /= np.tan(sun_altaz.alt.to(u.rad))  # correction for altitude

    # Constants in km
    re = R_earth.to(u.km).value
    ds = sun.distance.to(u.km).value
    rs = R_sun.to(u.km).value
    rm = 1738.1
    su = s_umbra.to(u.km).value

    dm = (2 * rm - su) * ds / (2 * (rs - rm)) + re
    print('Distance estimate: {:.7} km (Actual: {:.7}). Ratio: {:.2}. Difference: {:.7}'.format(dm, moon.distance,
                                                                                                moon.distance / (
                                                                                                dm * u.km),
                                                                                                dm * u.km - moon.distance))
    return dm, moon.distance.to(u.km).value


def estimate_umbra(loc, t2, t3):
    # Estimate the size of the umbra from the contact times
    delta_t = t3 - t2

    # Location of Sun
    sun = get_sun(t2)
    sun_altaz = sun.transform_to(AltAz(obstime=t2, location=loc))

    # Umbra size, from time and location
    s_umbra = 2 * np.pi * R_earth * np.cos(loc.latitude.to(u.rad)) / u.day * delta_t
    s_umbra /= np.tan(sun_altaz.alt.to(u.rad))  # correction for altitude

    return s_umbra.to(u.km).value


def get_umbra(loc, t2):
    # Get the real size of the umbra

    sun = get_sun(t2)
    moon = get_moon(t2, loc)
    sun_altaz = sun.transform_to(AltAz(obstime=t2, location=loc))

    # Constants in km
    re = R_earth.to(u.km).value
    ds = sun.distance.to(u.km).value
    rs = R_sun.to(u.km).value
    rm = 1738.1

    alpha = np.arcsin((rs - rm) / (ds - moon.distance.to(u.km).value))
    su = 2 * (rm / np.sin(alpha) - (moon.distance.to(u.km).value - re)) * np.tan(alpha)
    su / np.tan(sun_altaz.alt.to(u.rad))

    return su