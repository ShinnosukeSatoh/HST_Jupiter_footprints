""" spiceypy_test.py

Created on Mar 15, 2023
@author: Shin Satoh

Description:


Version
1.0.0 (Mar 15, 2023) Determination of Jupiter's coordinate

"""

#
# Solution convtm
#
import matplotlib.pyplot as plt
import numpy as np
import spiceypy as spice
import datetime


spice.furnsh('kernel/cassMetaK.txt')

# 木星半径
RJ = 71492   # [km]

# 時間刻み数
step = 2
# we are going to get positions between these two dates
utc = ['2014-03-09T01:03:30', '2022-05-22']

# get et values one and two, we could vectorize str2et
etOne = spice.str2et(utc[0])
etTwo = spice.str2et(utc[1])
print("ET One: {}, ET Two: {}".format(etOne, etTwo))

print(spice.et2utc(etOne, "ISOC", 6))

# get times
times = [x*(etTwo-etOne)/step + etOne for x in range(step)]

# check first few times:
# print(times[0:3])

#
#
# 木星からみたEuropaの座標
# Run spkpos as a vectorized function
positions, lightTimes = spice.spkpos(
    targ='EUROPA', et=etOne, ref='IAU_JUPITER', abcorr='NONE', obs='JUPITER')

# Positions is a 3xN vector of XYZ positions
print("Positions [km]: ")
print(positions)
print("Positions [RJ]: ")
print(np.sqrt(positions[0]**2 + positions[1]**2 + positions[2]**2)/RJ)

"""
# Run spkpos as a vectorized function
positions, lightTimes = spice.spkpos(
    targ='Cassini', et=times, ref='J2000', abcorr='NONE', obs='SATURN BARYCENTER')

# Positions is a 3xN vector of XYZ positions
# print("Positions: ")
# print(positions[0])

# Light times is a N vector of time
# print("Light Times: ")
# print(lightTimes[0])

# positions is shaped (4000, 3), let's transpose to (3, 4000) for easier indexing
positions = positions.T
fig = plt.figure(figsize=(9, 9))
ax = fig.add_subplot(111, projection='3d')
ax.plot(positions[0], positions[1], positions[2])
plt.title('SpiceyPy Cassini Position Example from Jun 20, 2004 to Dec 1, 2005')
plt.show()
"""
