""" hstprojimage_plots.py

Created on Apr 3, 2023
@author: Shin Satoh

Description:


Version
1.0.0 (Apr 3, 2023)

"""

import matplotlib.pyplot as plt
import hstprojimage as hst
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os

fontname = 'Nimbus Sans'
plt.rcParams.update({'font.sans-serif': fontname,
                    'font.family': 'sans-serif',
                     'mathtext.fontset': 'custom',
                     'mathtext.rm': fontname,
                     'mathtext.it': fontname+':italic',
                     # 'mathtext.bf': 'Nimbus Sans:italic:bold',
                     'mathtext.bf': fontname+':bold'
                     })

year = '2022'
yearlydata_dir = 'data/red/'+year

doy_visit_list = sorted(os.listdir(yearlydata_dir))
for doyvisit in doy_visit_list[12:]:
    savedir = 'img/red/'+year+'/'+doyvisit
    try:
        os.makedirs(savedir)
    except FileExistsError:
        savedir += '_r'
        os.makedirs(savedir)

    fits30s = sorted(os.listdir(yearlydata_dir+'/'+doyvisit))

    for fitsname in fits30s:
        plt.figure(1, figsize=(5, 6), dpi=200)
        plt.clf()
        fig, axs = plt.subplots(2, 1, num=1,
                                gridspec_kw={'height_ratios': [20, 1.6]})

        plt.subplots_adjust(top=1, bottom=0.07, hspace=0.4)

        # Colorbar axes
        divider = make_axes_locatable(axs[0])
        cax = divider.append_axes('bottom', size='3%', pad=0)

        # HST data
        h = hst.HSTProjImage(yearlydata_dir+'/'+doyvisit+'/'+fitsname)
        h.readHSTFile()
        ax = axs[0]
        h.tvPolar(ax, vmin=10, vmax=2000,
                  draw_labels=True, refmainoval=False, reflon=h.alm.cml)

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
            NORTHSOUTH = 'North'
        else:
            NS = 'S'
            NORTHSOUTH = 'South'
        eurftp = 'Europa footprint: '+str(round(h.s3wlon_lin, 2)) + \
            '°W, ' + str(round(abs(h.s3lat_lin), 2)) + '°'+NS

        cml = 'CML: '+h.CML+'°'
        info_txtL = doy + '\n' + cml + '\n' + eurftp
        info_txtR = visit + '\n' + integ
        ax.set_title(title, weight='bold', fontsize=fontsize, loc='left')
        ax.set_title('Jupiter '+NORTHSOUTH, weight='bold',
                     fontsize=fontsize, loc='right')

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
        ax.axis('off')

        savename, _ = fitsname.split('_stis_')
        print(savename)
        plt.savefig(savedir+'/'+savename+'.jpg')
        # plt.pause(10)
        plt.close()
