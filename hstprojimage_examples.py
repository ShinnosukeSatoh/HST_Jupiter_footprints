
# Some examples of the use of HSTProjImage.
# Filenames will have to be set to your system setup

import matplotlib.pyplot as plt
# (note I have a directory in my pythonpath called hst)
import hstprojimage as hst
from mpl_toolkits.axes_grid1 import make_axes_locatable

EXTRACTDIR = 'data/sample/'

fontname = 'Nimbus Sans'
plt.rcParams.update({'font.sans-serif': fontname,
                     'font.family': 'sans-serif',
                     'mathtext.fontset': 'custom',
                     'mathtext.rm': fontname,
                     'mathtext.it': fontname+':italic',
                     # 'mathtext.bf': 'Nimbus Sans:italic:bold',
                     'mathtext.bf': fontname+':bold'
                     })

"""
h = hst.HSTProjImage(
    EXTRACTDIR+'jup_22-142-10-30-41_0030_v06_stis_f25srf2_flatproj.fits')
h.readHSTFile()
plt.figure(1)
plt.clf()
f, axs = plt.subplots(1, 2, num=1)
ax = axs[0]
hax = h.tvPolar(ax, vmin=10, vmax=2000, refmainoval=True,
                draw_labels=False, reflon=h.alm.cml)

h = hst.HSTProjImage(
    EXTRACTDIR+'jup_22-142-10-30-41_0030_v06_stis_f25srf2_flatproj.fits')
h.readHSTFile()
ax = axs[1]
h.tvPolar(ax, vmin=10, vmax=1000, norm='lin',
          draw_labels=True, refmainoval=False)
plt.show()
print('two done')

h = hst.HSTProjImage(
    EXTRACTDIR+'jup_22-142-10-30-41_0030_v06_stis_f25srf2_flatproj.fits')
h.readHSTFile()
plt.figure(2)
plt.clf()
f, axs = plt.subplots(1, 1, num=2)
ax = axs
hax = h.tvProj(ax, vmin=10, vmax=2000, refmainoval=True,  draw_labels=True)
plt.show()
print('three done')
"""

# fitsname = 'jup_14-001-03-03-03_0030_v01_stis_f25srf2_flatproj.fits'
# fitsname = 'jup_22-140-15-38-53_0030_v03_stis_f25srf2_flatproj.fits'
fitsname = 'jup_22-142-10-30-41_0030_v06_stis_f25srf2_flatproj.fits'

plt.figure(1, figsize=(5, 6), dpi=150)
plt.clf()
fig, axs = plt.subplots(2, 1, num=1, gridspec_kw={
                        'height_ratios': [20, 1.6]})
# 上下を結合
plt.subplots_adjust(top=1, bottom=0.07, hspace=0.4)

# Colorbar axes
divider = make_axes_locatable(axs[0])
cax = divider.append_axes('bottom', size='3%', pad=0)

# HST data
h = hst.HSTProjImage(EXTRACTDIR+fitsname)
h.readHSTFile()
ax = axs[0]
h.tvPolar(ax, vmin=10, vmax=2000,
          draw_labels=True, refmainoval=False, reflon=h.alm.cml)
# h.tvProj(ax, vmin=10, vmax=2000,
#          draw_labels=True, refmainoval=False, reflon=h.alm.cml)

# Colorbar axes plot
axpos = ax.get_position()
pp = fig.colorbar(h.tvim, cax=cax, orientation='horizontal')
pp.set_label('Intensity [kR]', fontsize=11)

# Observation information
ax = axs[1]
fontsize = 12
title = 'STIS FUV-MAMA F25SRF2'
rootid = 'Root ID: '+h.obsid
filename = fitsname
obsdate = 'HST date: '+h.datetime
integ = 'Integration: '+h.exptime+' s'
OBSDAY, OBSTIME = h.datetime.split('T')
doy = OBSDAY + ' (DoY '+h.DOY + ')'
visit = OBSTIME + ' (Visit '+h.VISIT+')'

if h.s3lat_lin >= 0:
    NS = 'N'
else:
    NS = 'S'
eurftp = 'Europa footprint: '+str(round(h.s3wlon_lin, 2)) + \
    '°W, ' + str(round(abs(h.s3lat_lin), 2)) + '°'+NS

cml = 'CML: '+h.CML+'°'
info_txtL = doy + '\n' + cml + '\n' + eurftp
info_txtR = visit + '\n' + integ
ax.set_title(title, weight='bold', fontsize=fontsize, loc='left')
ax.set_title('Jupiter North', weight='bold', fontsize=fontsize, loc='right')

ax.text(0, 1, info_txtL,
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax.transAxes,
        fontsize=fontsize)
ax.text(1, 1, info_txtR,
        horizontalalignment='right',
        verticalalignment='top',
        transform=ax.transAxes,
        fontsize=fontsize)
ax.axis('off')                                # 全部消す

plt.savefig('img/HST_plot_'+h.obsid+'.jpg')
plt.show()
