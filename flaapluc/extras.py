#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Time-stamp: "2017-07-29 22:22:51 jlenain"

"""
Extras functions for FLaapLUC

@author Jean-Philippe Lenain <mailto:jlenain@in2p3.fr>
"""

import numpy as np


def met2mjd(met):
    """
    Converts Mission Elapsed Time (MET, in seconds) in Modified Julian Day.
    Cf. http://fermi.gsfc.nasa.gov/ssc/data/analysis/documentation/Cicerone/Cicerone_Data/Time_in_ScienceTools.html
    to see how the time is handled in the Fermi Science Tools.

    Input: time in MET (s)
    Output: time in MJD (fraction of a day)
    """
    MJDREFI = 51910.0
    MJDREFF = 7.428703703703703e-4
    return (MJDREFI+MJDREFF+met/24./60./60.)


def mjd2met(mjd):
    """
    Converts Modified Julian Day in Mission Elapsed Time (MET, in seconds).
    Cf. http://fermi.gsfc.nasa.gov/ssc/data/analysis/documentation/Cicerone/Cicerone_Data/Time_in_ScienceTools.html
    to see how the time is handled in the Fermi Science Tools.

    Input:  time in MJD (fraction of a day)
    Output: time in MET (s)
    """
    MJDREFI = 51910.0
    MJDREFF = 7.428703703703703e-4
    return (24. * 60. * 60 * (mjd - MJDREFI - MJDREFF))


def unixtime2mjd(unixtime):
    """
    Converts a UNIX time stamp in Modified Julian Day

    Input:  time in UNIX seconds
    Output: time in MJD (fraction of a day)

    """

    # unixtime gives seconds passed since "The Epoch": 1.1.1970 00:00
    # MJD at that time was 40587.0

    result = 40587.0 + unixtime / (24.*60.*60.)
    return result


def jd2gd(x):
    """
    Compute gregorian date out of julian date

    input: julian date x (float)
    return value: string of gregorian date

    based on/copied from script jd2dg.py from Enno Middelberg
    http://www.atnf.csiro.au/people/Enno.Middelberg/python/jd2gd.py

    task to convert a list of julian dates to gregorian dates
    description at http://mathforum.org/library/drmath/view/51907.html
    Original algorithm in Jean Meeus, "Astronomical Formulae for Calculators"
    """

    jd = float(x)

    jd = jd+0.5
    Z = int(jd)
    F = jd-Z
    alpha = int((Z-1867216.25)/36524.25)
    A = Z + 1 + alpha - int(alpha/4)

    B = A + 1524
    C = int((B-122.1)/365.25)
    D = int(365.25*C)
    E = int((B-D)/30.6001)

    dd = B - D - int(30.6001*E) + F

    if E < 13.5:
        mm = E-1

    if E > 13.5:
        mm = E-13

    if mm > 2.5:
        yyyy = C-4716

    if mm < 2.5:
        yyyy = C-4715

    daylist = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    daylist2 = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

    h = int((dd-int(dd))*24)
    min = int((((dd-int(dd))*24)-h)*60)
    sec = 86400*(dd-int(dd))-h*3600-min*60

    # Now calculate the fractional year. Do we have a leap year?
    if (yyyy % 4 != 0):
        days = daylist2
    elif (yyyy % 400 == 0):
        days = daylist2
    elif (yyyy % 100 == 0):
        days = daylist
    else:
        days = daylist2

    string = "%04d-%02d-%02d %02d:%02d:%04.1f" % (yyyy, mm, dd, h, min, sec)

    return string


def mjd2gd(time):
    """
    Converts Modified Julian Day in Gregorian Date.

    Under the hood, it calls jd2gd().
    """

    return jd2gd(time+2400000.5)


def rad2deg(angle):
    """
    Convert an angle from radians to degrees.

    @param angle in radians
    """
    return angle*180./pi


def deg2rad(angle):
    """
    Convert an angle from degrees to radians.

    @param angle in degrees
    """
    return angle*pi/180.


def angsep((ra1,dec1), (ra2,dec2), deg=True):
    """
    Calculates the angular separation between two points on the sky.

    @param (ra1,dec1) coordinates of 1st source
    @param (ra2,dec2) coordinates of 2nd source
    @param deg flag whether inputs/outputs are in degrees or radians
    """
    if deg:
        ra1 = deg2rad(ra1)
        dec1 = deg2rad(dec1)
        ra2 = deg2rad(ra2)
        dec2 = deg2rad(dec2)

    SEP = np.arccos(np.cos(dec1)*np.cos(dec2)*np.cos(ra1-ra2)+np.sin(dec1)*np.sin(dec2))
    if deg:
        SEP = rad2deg(SEP)
    return SEP
