import datetime
import pytz
import numpy as np
import time
from datetime import datetime, date, tzinfo
import math
from math import pi, sin, asin, cos, acos, sqrt, atan2

from math import floor

SIDEREAL_DAY = 86164.0905  # Sidereal day in seconds
SOLAR_DAY = SIDEREAL_DAY - 235.9095
LUNAR_DAY = SIDEREAL_DAY - 2089.2292
SIDERAL_ANGLE_RATE = (2 * pi) / SIDEREAL_DAY
SOLAR_ANGLE_RATE = (2 * pi) / SOLAR_DAY
LUNAR_ANGLE_RATE = (2 * pi) / LUNAR_DAY


def JulianDate(year, month, day, utc=0):
    """
    Returns the Julian date, number of days since 1 January 4713 BC 12:00.
    utc is UTC in decimal hours. If utc=0, returns the date at 12:00 UTC.
    """
    if month > 2:
        y = year
        m = month
    else:
        y = year - 1
        m = month + 12
    d = day
    h = utc / 24
    if year <= 1582 and month <= 10 and day <= 4:
        # Julian calendar
        b = 0
    elif year == 1582 and month == 10 and day > 4 and day < 15:
        # Gregorian calendar reform: 10 days (5 to 14 October 1582) were skipped.
        # In 1582 after 4 October follows the 15 October.
        d = 15
        b = -10
    else:
        # Gregorian Calendar
        a = int(y / 100)
        b = 2 - a + int(a / 4)
    jd = int(365.25 * (y + 4716)) + int(30.6001 * (m + 1)) + d + h + b - 1524.5
    return jd


def SiderialTime(year, month, day, utc=0, longitude=0):
    # Returns the siderial time in decimal hours. Longitude (long) is in decimal degrees.
    # If long=0, return value is Greenwich Mean Siderial Time (GMST).

    jd = JulianDate(year, month, day)
    t = (jd - 2451545.0) / 36525
    # Greenwich siderial time at 0h UTC (hours)
    st = (
        24110.54841 + 8640184.812866 * t + 0.093104 * t ** 2 - 0.0000062 * t ** 3
    ) / 3600
    # Greenwich siderial time at given UTC
    st = st + 1.00273790935 * utc
    # Local siderial time at given UTC (longitude in degrees)
    st = st + longitude / 15
    st = st % 24
    return st


# // The trig-value matrix on page 103 of the Explanatory Supplement is NOT THE SAME as
# // the one at http://chsfpc5.chem.ncsu.edu/~franzen/CH795Z/math/lab_frame/lab_frame.html
# // Further, I had to use the transpose to get the right answer (matching other 1950->2000
# // results), presumably because the Supplement does a post- not a pre- multiply or
# // the other way around.
def MakeMatrixSupplement(a, b, c):
    # zeta(supp) is A, theta(supp) is B; and the third greek letter is C
    m = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
    cA = cos(a)
    sA = sin(a)
    cB = cos(b)
    sB = sin(b)
    cC = cos(c)
    sC = sin(c)
    m[0][0] = cA * cB * cC - sA * sC
    m[1][0] = -cA * cB * sC - sA * cC
    m[2][0] = -cA * sB
    m[0][1] = sA * cB * cC + cA * sC
    m[1][1] = -sA * cB * sC + cA * cC
    m[2][1] = -sA * sB
    m[0][2] = sB * cC
    m[1][2] = -sB * sC
    m[2][2] = cB
    return m


#  Calculate the three key angles (A B and C in the above)
#  Supplement page 104
#  To match the above, theta(supp) is B; zeta(supp) is A, and the third greek letter is C.
def Supplement(fixed, date):  # Here years, maybe should be JDays
    T = (fixed - 2000.0) / 100.0
    t = (date - fixed) / 100.0
    # should be Julian days / 36535
    asec = (2306.218 + 1.397 * T) * t + 1.095 * t ** 2
    bsec = (2004.311 - 0.853 * T) * t - 0.427 * t ** 2
    csec = (2306.218 + 1.397 * T) * t + 0.302 * t ** 2
    SecondsPerRadian = (180.0 / pi) * 3600.0
    m = MakeMatrixSupplement(
        asec / SecondsPerRadian, bsec / SecondsPerRadian, csec / SecondsPerRadian
    )
    return m


def Transform(ra, dec, matrix):  # returns  ra  & dec
    r0 = [cos(ra) * cos(dec), sin(ra) * cos(dec), sin(dec)]
    s0 = [
        r0[0] * matrix[0][0] + r0[1] * matrix[0][1] + r0[2] * matrix[0][2],
        r0[0] * matrix[1][0] + r0[1] * matrix[1][1] + r0[2] * matrix[1][2],
        r0[0] * matrix[2][0] + r0[1] * matrix[2][1] + r0[2] * matrix[2][2],
    ]
    r = sqrt(s0[0] ** 2 + s0[1] ** 2 + s0[2] ** 2)
    r_dec = asin(s0[2] / r)  # New dec in range -90.0 -- +90.0
    # or use sin^2 + cos^2 = 1.0
    cosaa = (s0[0] / r) / cos(r_dec)
    sinaa = (s0[1] / r) / cos(r_dec)
    r_ra = atan2(sinaa, cosaa)
    if r_ra < 0.0:
        r_ra = r_ra + pi + pi
    return r_ra, r_dec


# Fractional part
def frac(x):
    x -= int(x)
    if x < 0:
        return x + 1
    else:
        return x


# Map a time in hours to the range  0  to 24
def Map24(hour):
    if hour < 0.0:
        n = int(hour / 24.0) - 1
        return hour - n * 24.0
    elif hour >= 24.0:
        n = (int)(hour / 24.0)
        return hour - n * 24.0
    else:
        return hour


# Compute Greenwich Mean Sidereal Time (gmst)
# TU is number of Julian centuries since 2000 January 1.5
# Expression for gmst from the Astronomical Almanac Supplement
def CalcLST(year, month, day, ut, SiteLongitude):
    TU = (CalcJD(year, month, day, 0.0) - 2451545.0) / 36525.0
    TU2 = TU * TU
    TU3 = TU2 * TU
    T0 = (
        (24110.54841 / 3600.0)
        + 8640184.812866 / 3600.0 * TU
        + 0.093104 / 3600.0 * TU2
        - 6.2e-6 / 3600.0 * TU3
    )
    T0 = Map24(T0)
    gmst = Map24(T0 + ut * 1.002737909)
    lmst = 24.0 * frac((gmst - SiteLongitude / 15.0) / 24.0)
    return lmst


# Calculate the local sidereal time with millisecond resolution
def LSTNow(SiteLongitude):
    now = datetime.utcnow()
    year = now.year
    month = now.month
    day = now.day
    minutes = now.minute
    hours = now.hour
    seconds = now.second
    usecs = now.microsecond
    milliseconds = usecs / 1000
    # Calculate floating point ut in hours
    ut = (milliseconds / 3600000) + (seconds / 3600) + (minutes / 60) + hours
    lst = CalcLST(year, month, day, ut, SiteLongitude)
    return lst


def CalcJD(ny, nm, nd, ut):
    day = nd + ut / 24.0
    if nm == 1 or nm == 2:
        ny = ny - 1
        nm = nm + 12

    if ny + nm / 12.0 + day / 365.25 >= (1582.0 + 10.0 / 12.0 + 15.0 / 365.25):
        A = int(ny / 100.0)
        B = 2.0 - A + int(A / 4.0)
    else:
        B = 0.0

    if ny < 0.0:
        C = int(365.25 * ny - 0.75)
    else:
        C = int(365.25 * ny)
    D = int(30.6001 * (nm + 1))
    jd = B + C + D + day + 1720994.5
    return jd


# Calculate the Julian date with millisecond resolution
def JDNow():
    now = datetime.utcnow()
    year = now.year
    month = now.month
    day = now.day
    minutes = now.minute
    hours = now.hour
    seconds = now.second
    usecs = now.microsecond
    milliseconds = usecs / 1000
    # Calculate floating point ut in hours
    ut = (milliseconds / 3600000) + (seconds / 3600) + (minutes / 60) + hours
    jd = CalcJD(year, month, day, ut)
    return jd


# Precession from J2000 to EOD or back
def Precession(ra_rad, dec_rad, dirflag):
    if dirflag > 0:
        return PrecessToEOD(2000.0, ra_rad, dec_rad)
    else:
        return PrecessToEpoch(2000.0, ra_rad, dec_rad)


# Precess in place from epoch to EOD
# Coordinates ra in RAD and dec in RAD
# Call this in the form PrecessToEOD(epoch,ra_rad,dec_rad)
def PrecessToEOD(epoch, ra_rad, dec_rad):
    jdfixed = (epoch - 2000.0) * 365.25 + 2451545.0
    # JD for epoch of date
    jdnow = JDNow()
    # Julian centuries for the fixed epoch from a base epoch 2000.0
    T = (jdfixed - 2451545.0) / 36525.0
    # Julian centuries for the epoch of date from the fixed epoch
    t = (jdnow - jdfixed) / 36525.0
    # Evaluate the constants in arc seconds
    zeta = (
        (2306.2181 + 1.39656 * T - 0.000139 * T * T) * t
        + (0.30188 - 0.000344 * T) * t * t
        + (0.017998) * t * t * t
    )
    z = (
        (2306.2181 + 1.39656 * T - 0.000139 * T * T) * t
        + (1.09468 + 0.000066 * T) * t * t
        + (0.018203) * t * t * t
    )
    theta = (
        (2004.3109 - 0.85330 * T - 0.000217 * T * T) * t
        + (-0.42665 - 0.000217 * T) * t * t
        + (-0.041833) * t * t * t
    )
    # Convert to radians
    zeta_rad = zeta * pi / (180.0 * 3600.0)
    z_rad = z * pi / (180.0 * 3600.0)
    theta_rad = theta * pi / (180.0 * 3600.0)
    # Calculate the precession
    a = sin(ra_rad + zeta_rad) * cos(dec_rad)
    b = cos(ra_rad + zeta_rad) * cos(theta_rad) * cos(dec_rad) - sin(theta_rad) * sin(
        dec_rad
    )
    c = cos(ra_rad + zeta_rad) * sin(theta_rad) * cos(dec_rad) + cos(theta_rad) * sin(
        dec_rad
    )
    if c > 0.9:
        dec_out = math.acos(sqrt(a * a + b * b))
    elif c < -0.9:
        dec_out = -math.acos(sqrt(a * a + b * b))
    else:
        dec_out = asin(c)
    ra_out = atan2(a, b) + z_rad
    if dec_out > pi / 2:
        dec_out = pi - dec_out
        ra_out = ra_out + pi

    if dec_out < -pi / 2:
        dec_out = -pi - dec_out
        ra_out = ra_out + pi
    return ra_out, dec_out


# Precess from EOD to epoch
# Mean coordinates ra in rad and dec in rad returned in place
# Will not remove nutation and aberration
# Call this in the form PrecessToEpoch(epoch,ra_rad,dec_rad)
def PrecessToEpoch(epoch, ra_rad, dec_rad):
    # JD for the fixed epoch
    jdfixed = (epoch - 2000.0) * 365.25 + 2451545.0
    # JD for epoch of date
    jdnow = JDNow()
    # Julian centuries for the fixed epoch from a base epoch 2000.0
    T = (jdnow - 2451545.0) / 36525.0
    # Julian centuries for the epoch of date from the fixed epoch
    t = (jdfixed - jdnow) / 36525.0
    # Evaluate the constants in arc seconds

    zeta = (
        (2306.2181 + 1.39656 * T - 0.000139 * T * T) * t
        + (0.30188 - 0.000344 * T) * t * t
        + (0.017998) * t * t * t
    )
    z = (
        (2306.2181 + 1.39656 * T - 0.000139 * T * T) * t
        + (1.09468 + 0.000066 * T) * t * t
        + (0.018203) * t * t * t
    )
    theta = (
        (2004.3109 - 0.85330 * T - 0.000217 * T * T) * t
        + (-0.42665 - 0.000217 * T) * t * t
        + (-0.041833) * t * t * t
    )
    # Convert to radians
    zeta_rad = zeta * pi / (180.0 * 3600.0)
    z_rad = z * pi / (180.0 * 3600.0)
    theta_rad = theta * pi / (180.0 * 3600.0)

    # Calculate the precession
    a = sin(ra_rad + zeta_rad) * cos(dec_rad)
    b = cos(ra_rad + zeta_rad) * cos(theta_rad) * cos(dec_rad) - sin(theta_rad) * sin(
        dec_rad
    )
    c = cos(ra_rad + zeta_rad) * sin(theta_rad) * cos(dec_rad) + cos(theta_rad) * sin(
        dec_rad
    )
    if c > 0.9:
        dec_out = math.acos(sqrt(a * a + b * b))
    elif c < -0.9:
        dec_out = -math.acos(sqrt(a * a + b * b))
    else:
        dec_out = asin(c)
    ra_out = atan2(a, b) + z_rad
    if dec_out > pi / 2:
        dec_out = pi - dec_out
        ra_out = ra_out + pi
    if dec_out < -pi / 2:
        dec_out = -pi - dec_out
        ra_out = ra_out + pi
    return ra_out, dec_out


def gms(grad_f):
    grad = int(grad_f)
    if grad < 0:
        grad_f = grad_f * -1
    minuten = int((grad_f - int(grad_f)) * 60)
    sekunden = (((grad_f - int(grad_f)) * 60) - int((grad_f - int(grad_f)) * 60)) * 60
    return "{grad}°{minuten}´{sekunden:.3f}´´".format(
        grad=grad, minuten=minuten, sekunden=sekunden
    )


# Now we have the RA, DEC and HA for the object,
# and the Latitude (LAT) of the observing site,
# the following formulas will give us the ALT and
# AZ of the object at the current LST.


def HA_DEC_to_AZ_ALT(ha_rad, dec_rad, latitude_rad):
    # sin(ALT) = sin(DEC)*sin(LAT)+cos(DEC)*cos(LAT)*cos(HA)
    sin_ALT_rad = sin(dec_rad) * sin(latitude_rad) + cos(dec_rad) * cos(
        latitude_rad
    ) * cos(ha_rad)
    ALT_rad = asin(sin_ALT_rad)
    #            sin(DEC) - sin(ALT)*sin(LAT)
    # cos(A)   =   ---------------------------------
    #                cos(ALT)*cos(LAT)
    cos_a_rad = (sin(dec_rad) - sin_ALT_rad * sin(latitude_rad)) / (
        cos(ALT_rad) * cos(latitude_rad)
    )
    AZ_rad = math.acos(cos_a_rad)
    return AZ_rad, ALT_rad
