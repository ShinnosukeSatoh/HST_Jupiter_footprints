""" ftpS3.py

Created on Mar 15, 2023
@author: Shin Satoh

Description:


Version
1.0.0 (Mar 15, 2023)
2.0.0 (Apr 4, 2023) Io's and Ganymede's footprint indicated.
"""

import spiceypy as spice
import numpy as np


class ftpS3():
    def __init__(self):
        return None

    def FP(self, utc, satmodel, MOON: str):
        # Position of the Galilean moons in the "IAU_JUPITER" frame
        # Calculating position of the moon by spiceypy.
        # Then converting to the System III longitude.

        # spice.furnsh('kernel/cassMetaK.txt')
        et_hst = spice.str2et(utc)

        # HST's position seen from the HST in IAU_JUPITER coordinate.
        _, lightTimes = spice.spkpos(
            targ='HST', et=et_hst, ref='IAU_JUPITER', abcorr='LT+S', obs='JUPITER'
        )

        # Moon's position seen from Jupiter in IAU_JUPITER coordinate.
        pos, _ = spice.spkpos(
            targ=MOON, et=et_hst-lightTimes, ref='IAU_JUPITER', abcorr='none', obs='JUPITER'
        )

        posx, posy, posz = pos[0], pos[1], pos[2]
        posr = np.sqrt(posx**2 + posy**2 + posz**2)
        # postheta = np.arccos(posz/posr)
        posphi = np.arctan2(posy, posx)
        if posphi < 0:
            Sys3 = np.degrees(-posphi)
        else:
            Sys3 = np.degrees(2*np.pi - posphi)
        print('R', posr/71492, ', PHI', np.degrees(posphi), ', SYS3', Sys3)

        # Search the System III index
        s3list = np.arange(0, 361, 5)
        argsorted = np.argsort(np.abs(s3list-Sys3), axis=0)
        s3_idx0, s3_idx1 = argsorted[0], argsorted[1]

        MOONs3 = [s3list[s3_idx0], s3list[s3_idx1]]
        # print('MOONS3', MOONs3)     # [350, 355]のときに計算が狂う

        if MOON == 'IO':
            s3wlon = [satmodel.iowlon[s3_idx0], satmodel.iowlon[s3_idx1]]
            s3lat = [satmodel.iolat[s3_idx0], satmodel.iolat[s3_idx1]]

        if MOON == 'EUROPA':
            s3wlon = [satmodel.euwlon[s3_idx0], satmodel.euwlon[s3_idx1]]
            s3lat = [satmodel.eulat[s3_idx0], satmodel.eulat[s3_idx1]]

        if MOON == 'GANYMEDE':
            s3wlon = [satmodel.gawlon[s3_idx0], satmodel.gawlon[s3_idx1]]
            s3lat = [satmodel.galat[s3_idx0], satmodel.galat[s3_idx1]]

        # print('s3wlon', s3wlon)
        # print('s3lat', s3lat)

        if (s3wlon[0]-s3wlon[1]) > 300:
            ds3wlon10 = (s3wlon[1]+360)-s3wlon[0]
            y0 = s3wlon[0]
        elif (s3wlon[1]-s3wlon[0]) > 300:
            ds3wlon10 = s3wlon[1]-(s3wlon[0]+360)
            y0 = s3wlon[0]+360
        else:
            ds3wlon10 = s3wlon[1]-s3wlon[0]
            y0 = s3wlon[0]

        s3wlon_lin = (ds3wlon10/(MOONs3[1]-MOONs3[0]))*(Sys3-MOONs3[0]) + y0
        s3lat_lin = (
            (s3lat[1]-s3lat[0])/(MOONs3[1]-MOONs3[0]))*(Sys3-MOONs3[0]) + s3lat[0]

        if s3wlon_lin > 360:
            s3wlon_lin += -360

        print('EFP LAT', s3lat_lin, ', SYS3', s3wlon_lin)

        return Sys3, s3wlon_lin, s3lat_lin

    def EFP(self, utc, satmodel):
        # Position of the Galilean moons in the "IAU_JUPITER" frame
        # Calculating position of the moon by spiceypy.
        # Then converting to the System III longitude.
        spice.furnsh('kernel/cassMetaK.txt')
        et_hst = spice.str2et(utc)

        # Europa's position seen from the HST in IAU_JUPITER coordinate.
        _, lightTimes = spice.spkpos(
            targ='HST', et=et_hst, ref='IAU_JUPITER', abcorr='LT+S', obs='JUPITER'
        )

        # Europa's position seen from Jupiter in IAU_JUPITER coordinate.
        pos, _ = spice.spkpos(
            targ='EUROPA', et=et_hst-lightTimes, ref='IAU_JUPITER', abcorr='none', obs='JUPITER'
        )

        posx, posy, posz = pos[0], pos[1], pos[2]
        posr = np.sqrt(posx**2 + posy**2 + posz**2)
        # postheta = np.arccos(posz/posr)
        posphi = np.arctan2(posy, posx)
        if posphi < 0:
            Sys3 = np.degrees(-posphi)
        else:
            Sys3 = np.degrees(2*np.pi - posphi)
        print('R', posr/71492, ', PHI', np.degrees(posphi), ', SYS3', Sys3)

        # Search the System III index
        s3list = np.arange(0, 361, 5)
        argsorted = np.argsort(np.abs(s3list-Sys3), axis=0)
        s3_idx0, s3_idx1 = argsorted[0], argsorted[1]

        EUs3 = [s3list[s3_idx0], s3list[s3_idx1]]
        print('EUS3', EUs3)     # [350, 355]のときに計算が狂う

        s3wlon = [satmodel.euwlon[s3_idx0], satmodel.euwlon[s3_idx1]]
        s3lat = [satmodel.eulat[s3_idx0], satmodel.eulat[s3_idx1]]
        print('s3wlon', s3wlon)
        print('s3lat', s3lat)

        if (s3wlon[0]-s3wlon[1]) > 300:
            ds3wlon10 = (s3wlon[1]+360)-s3wlon[0]
            y0 = s3wlon[0]
        elif (s3wlon[1]-s3wlon[0]) > 300:
            ds3wlon10 = s3wlon[1]-(s3wlon[0]+360)
            y0 = s3wlon[0]+360
        else:
            ds3wlon10 = s3wlon[1]-s3wlon[0]
            y0 = s3wlon[0]

        s3wlon_lin = (ds3wlon10/(EUs3[1]-EUs3[0]))*(Sys3-EUs3[0]) + y0
        s3lat_lin = (
            (s3lat[1]-s3lat[0])/(EUs3[1]-EUs3[0]))*(Sys3-EUs3[0]) + s3lat[0]

        if s3wlon_lin > 360:
            s3wlon_lin += -360

        print('FP LAT', s3lat_lin, ', SYS3', s3wlon_lin)

        return s3wlon_lin, s3lat_lin
