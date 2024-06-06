""" hstprojimage_plots.py

Created on Apr 3, 2023
@author: Shin Satoh

Description:


Version
1.0.0 (Apr 3, 2023)
2.0.0 (Apr 4, 2023) Io's and Ganymede's footprint indicated.
"""

import matplotlib.pyplot as plt
import hst_unproj as hst
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os
import spiceypy as spice

fontname = 'Nimbus Sans'
plt.rcParams.update({'font.sans-serif': fontname,
                    'font.family': 'sans-serif',
                     'mathtext.fontset': 'custom',
                     'mathtext.rm': fontname,
                     'mathtext.it': fontname+':italic',
                     # 'mathtext.bf': 'Nimbus Sans:italic:bold',
                     'mathtext.bf': fontname+':bold'
                     })

moon = 'EUROPA'
year = '2014'
yearlydata_dir = 'data/red/'+year
spice.furnsh('kernel/cassMetaK.txt')

doy_visit_list = sorted(os.listdir(yearlydata_dir))
print(doy_visit_list)
for doyvisit in doy_visit_list[3:4]:
    savedir = 'img/red/'+moon+'/'+year+'/'+doyvisit

    fits30s = sorted(os.listdir(yearlydata_dir+'/'+doyvisit))

    for fitsname in fits30s:

        fig, ax = plt.subplots()

        # Colorbar axes
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('bottom', size='3%', pad=0)

        # HST data
        h = hst.HSTProjImage(yearlydata_dir+'/'+doyvisit+'/'+fitsname)
        h.readHSTFile()
        h.MOON = moon
        h.tvPolar(ax, vmin=10, vmax=2000,
                  # h.tvPolar(ax, vmin=50, vmax=60000,
                  draw_labels=True, refmainoval=False,
                  satovals=['all'],
                  # reflon=None,
                  grid=True,
                  reflon=h.alm.cml,
                  ext=None,
                  )

        # Colorbar axes plot
        axpos = ax.get_position()
        pp = fig.colorbar(h.tvim, cax=cax, orientation='horizontal')
        pp.set_label('Intensity [kR]', fontsize=11)

        ax.axis('off')

        savename, _ = fitsname.split('_stis_')
        savename = savename+'_'+moon.lower()[0]+'fp'
        print(savename)
        # plt.savefig(savedir+'/'+savename+'.jpg')
        # plt.pause(10)
        plt.show()
        # plt.close()

        del h

        break
