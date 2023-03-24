
# Some examples of the use of HSTProjImage.
# Filenames will have to be set to your system setup

import matplotlib.pyplot as plt
# (note I have a directory in my pythonpath called hst)
import hstprojimage as hst

EXTRACTDIR = 'data/sample/'

h = hst.HSTProjImage(
    EXTRACTDIR+'jup_22-142-10-30-41_0030_v06_stis_f25srf2_flatproj.fits')
h.readHSTFile()
plt.figure(1)
plt.clf()
f, axs = plt.subplots(1, 2, num=1)
ax = axs[0]
hax = h.tvPolar(ax, vmin=10, vmax=2000, refmainoval=True,
                draw_labels=False, reflon=h.alm.cml)

plt.show()
h = hst.HSTProjImage(
    EXTRACTDIR+'jup_22-140-15-38-53_0030_v03_stis_f25srf2_flatproj.fits')
h.readHSTFile()
ax = axs[1]
h.tvPolar(ax, vmin=10, vmax=1000, norm='lin',
          draw_labels=True, refmainoval=False)
plt.show()


h = hst.HSTProjImage(
    EXTRACTDIR+'jup_22-142-10-30-41_0030_v06_stis_f25srf2_flatproj.fits')
h.readHSTFile()
plt.figure(2)
plt.clf()
f, axs = plt.subplots(1, 1, num=2)
ax = axs
hax = h.tvProj(ax, vmin=10, vmax=2000, refmainoval=True,  draw_labels=True)
