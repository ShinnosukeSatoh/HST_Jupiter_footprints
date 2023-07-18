""" localtime.py

Created on Jun 13, 2023
@author: Shin Satoh

Description:


Version
1.0.0 (Jun 13, 2023)
"""

import spiceypy as spice
import numpy as np
import math
import datetime


class LT():
    def __init__(self):
        return None

    def S3(self, pos_arr):
        """_summary_

        Args:
            pos_arr (1d-ndarray): object position in IAU_JUPITER [m]

        Returns:
            Sys3: System III longitude of the object [deg]
            posphi: longitude in S3RH [deg]
        """

        posx, posy, posz = pos_arr[0], pos_arr[1], pos_arr[2]
        # posr = np.sqrt(posx**2 + posy**2 + posz**2)
        # postheta = np.arccos(posz/posr)
        posphi = np.arctan2(posy, posx)
        if posphi < 0:
            Sys3 = np.degrees(-posphi)
        else:
            Sys3 = np.degrees(2*np.pi - posphi)

        return Sys3, np.degrees(posphi)

    def lt(self, utc):
        """_summary_

        Args:
            utc: observation date

        Returns:
           elong: east longitude of Europa (0 is defined at the anti-solar position) [deg]
           td: local time [sec]
        """
        et_hst = spice.str2et(utc)

        # The sun's position seen from the Jupiter in IAU_JUPITER coordinate.
        posSUN, _ = spice.spkpos(
            targ='SUN', et=et_hst, ref='IAU_JUPITER', abcorr='NONE', obs='JUPITER'
        )

        S3wlon_SUN, S3RH_SUN = self.S3(posSUN)

        # Europa's position seen from the Jupiter in IAU_JUPITER coordinate.
        posEUR, _ = spice.spkpos(
            targ='EUROPA', et=et_hst, ref='IAU_JUPITER', abcorr='NONE', obs='JUPITER'
        )

        S3wlon_EUR, S3RH_EUR = self.S3(posEUR)

        dot = posSUN[0]*posEUR[0]+posSUN[1]*posEUR[1]
        R_SUN = math.sqrt(posSUN[0]**2+posSUN[1]**2+posSUN[2]**2)
        R_EUR = math.sqrt(posEUR[0]**2+posEUR[1]**2+posEUR[2]**2)
        arg = math.degrees(math.acos(dot/(R_SUN*R_EUR)))

        # print(posEUR)
        # print(arg)
        # print('SUN:', S3wlon_SUN, S3RH_SUN)
        # print('EUR:', S3wlon_EUR, S3RH_EUR)

        # Dawn or dusk
        yaxis = np.array([-posSUN[1], posSUN[0]])
        dot = yaxis[0]*posEUR[0]+yaxis[1]*posEUR[1]
        if dot >= 0:
            # DUSK側
            # print('DUSK')
            elong = arg + 180
        else:
            # DAWN側
            # print('DAWN')
            elong = 180 - arg

        sec = (3600*24/360)*elong    # [sec]
        # print(elong)
        # print(td)
        return elong, sec


spice.furnsh('kernel/cassMetaK.txt')
for i in range(8):
    elong, sec = LT().lt('2015-01-01T0'+str(i)+':00:00')
    td = datetime.timedelta(seconds=sec)

# a = np.load('ref/refmon.npz')
# print(a.files)
# print(a['hilat'])
