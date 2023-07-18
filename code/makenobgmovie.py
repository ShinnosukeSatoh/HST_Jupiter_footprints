import numpy as np
import matplotlib.pyplot as plt
import pdb
import glob
from datetime import timedelta
import datetime
from datetime import datetime as dt
from astropy.time import Time
from importlib import reload
import hst.hstimage as hst
reload(hst)
import matplotlib.gridspec as gridspec
# from multiprocessing import Pool 
from pathos.multiprocessing import Pool
import subprocess
import warnings
import spacecraft.junopjajtimes as pjaj
from pathlib import Path
from numba import njit
plt.ioff()

def createMovieFrame(i):
    """docstring for createMovieFrame"""
    
    print('Processing image {:0>3d}'.format(i),':', files[i])
    # return
    fim = hst.NewHSTImage(files[i])
    fim.readHSTFile(cr=2.5)
    pim = hst.HSTProjImage(piles[i])
    pim.readHSTFile(cr=2.5)

# New figure
    f = plt.figure(1, figsize=(4.2,  7.1))
    f.clf()
    f = plt.figure(1, figsize=(4.2,  7.1))
    # dpi = 120
    dpi = 200
    f.set_dpi(dpi)


# Plot unprojected image
    if fim.alm.hemisph == 'north':
        shx=-200.*np.sin(np.deg2rad(cml0))
        # exparams = [(314, 900), (822,699+shx)]
        exparams = [(314, 900), (175,699+shx)]
    else:
        shx=-80.*np.cos(np.deg2rad(cml0))
        # exparams = [(314, 900), (572,699+shx)]
        exparams = [(314, 900), (180,699+shx)]
        
    fim.extractImage(exparams[0],exparams[1])

    # ax = f.add_axes([0,0.1,1,0.6])
    ax = plt.subplot(gs[0], aspect='equal')
    fim.tvImage(vmin=10,vmax=maxint,delon=20, grid=True)
    alm = fim.alm

    ax2 = plt.subplot(gs[1], aspect='equal')
    #   shift slightly upward
    oldpos = ax2.get_position().bounds
    newpos = tuple(np.array(oldpos))+np.array([0,0.042,0,0])
    ax2.set_position(newpos)
    mo = True#{'north':True, 'south':False}[pim.alm.hemisph]
    yc = {'north':85, 'south':90}[pim.alm.hemisph]
    if reflonstr == '_cmls':
        reflon = pim.alm.cmls
    else:
        reflon = None
    pimax, xyimg, gl = pim.tvPolar(ax2, vmin=10,vmax=maxint, refmainoval=mo, yclip=yc, badcol='k', satovals=['all'], draw_labels=True, reflon=reflon)
    gl.rotate_labels=False
    gl.xlabel_style = {'c':'k', 'size':'xx-small'}
    gl.ylabel_style = {'c':'k', 'size':'xx-small'}
    latmod = {'north':1, 'south':-1}[pim.alm.hemisph]
    pimax.plot([pim.alm.cmls,pim.alm.cmls], [48*latmod,50*latmod], 'r', transform=pim.geodetic, lw=1.5) 

    
 # Colourbar
    ticks = {300:[1,10,100,300],
             800:[1,10,100,800],
             1000:[1,10,100,1000],
             2000:[1,10,100,1000],
             3000:[1,10,100,1000,3000],
             5000:[1,10,100,1000,3000],
             10000:[1,10,100,1000,10000]}[maxint]
    cax = f.add_axes([0.05, 0.105, 0.9, 0.02])
    cbar = f.colorbar(fim.tvim, orientation='horizontal',format='%i', cax=cax,
                      extend='both', ticks=ticks)
    # cbar.locator = MaxNLocator( nbins = 3)
    cbar.ax.tick_params(labelsize=12, length=5)
    # cbar.ax.tick_params(which='minor', length=2.5)
    # minorticks = fim.tvim.norm(np.arange(1, 10, 2))
    # cbar.ax.xaxis.set_ticks(minorticks, minor=True)
    cbar.ax.set_xticks(cbar.ax.xaxis.get_ticklocs()[2:])
    cbar.set_label('Intensity (kR)', size=12, labelpad=2)
    cbar.solids.set_edgecolor("face")
    # cbar.ax.minorticks_on()
    # cbar = fig.colorbar(fim.tvim, )
    cbar.update_ticks()
    cml = alm.cml % 360.
    lcml = (360. - cml) % 360.
    dece = alm.dece
    ud = alm.udate
    lt = round(alm.lightime)


# Operations on the image time
    t = dt.strptime(ud, '%Y-%m-%d %H:%M:%S')
    for pj, ajdt in enumerate(ajt.datetime):
        if t < ajdt:
            break
    datestr = t.date().isoformat()
    doy = t.utctimetuple().tm_yday
    doystr = str(doy)

    delt = timedelta(seconds = lt)
    ts = t - delt
    tmids = tmid - delt
    tstarts = tstart - delt
    tends = tend - delt
    tsdate = ts.date()
    timestr = t.time().isoformat()
    stimstr = ts.time().isoformat()
    pjstr = fim.alm.tardescr[-2:]
    hemistr = pim.alm.hemisph.capitalize()





    ta = ax.transAxes
    tsize = 10.
    ax.text(0.0, 1.42, 'Jupiter '+hemistr, transform = ta, size=13)
    # datestr = datestr +' (DoY ' +doystr+ ') PJ {:>02d}'.format(pj)
    datestr = datestr +' (DoY ' +doystr+ ')'
    ax.text(0.0, 1.24, datestr, transform = ta, size=tsize)
    ax.text(0.0, 1.07, timestr + ' (Earth) ', transform = ta, size=tsize)
    ax.text(0.34, 1.07, stimstr + ' (Jupiter)', transform = ta, size=tsize)
    ax.text(1., 1.07,
        r'$\mathrm{{CML}} (\lambda_{{III}}) = {:.1f}^\circ{{}}$'.format(alm.cml), ha = 'right',
        transform = ta, size=tsize)

    ax.text(1., 1.42, 'Visit ' + alm.visit[:2], ha = 'right', transform = ta, size=tsize)
    ax.text(1., 1.24,  alm.filter, ha = 'right', transform = ta, size=tsize)
    # phistr = r'$\lambda_{{III}} = {:.0f}^\circ{{}}$'.format(alm.cml)
    # ax.text(1, 1.08, phistr, ha = 'right', transform = ta)
    sigtext = 'HST/STIS.  NASA, ESA, J. D. Nichols'

    ax2.text(0.5, -0.21, sigtext, ha = 'center', size=11,
            transform = ax2.transAxes)

    plt.savefig('/Users/jdn4/data/HST/jupiter/'+datestr[:4]+'/movies/ims_'+'{:0>3d}'.format(i)+'.jpg', transparent=False, dpi=dpi, pil_kwargs={'quality':90})

    # plt.close()
    del(f)
    del(fim)
    del(pim)
    # plt.draw()
    # plt.show()
    # import pdb; pdb.set_trace()
        
        
"""docstring for showim"""

# Backend needs to be 'agg' for the multiprocessing to work

if __name__ == '__main__':

    # imglist = np.recfromtxt('/Users/jdn4/data/HST/jupiter/junoimglist.dat',
    # names=['filename', 'rootname', 'hemisph', 'filter', 'cml1', 'cml2'], encoding=None)
    # inx = imglist.hemisph == 'north'
    # files = imglist.filename[inx]

    visits = [
            #   '001_v01', 2014
            #   '002_v02',
            #   '003_v03',
            #   '004_v04',
            #   '005_v05',
            #   '006_v06',
            #   '007_v07',
            #   '010_v09',
            #   '011_v10',
            #   '011_v11',
            #   '013_v08',
            #   '013_v13',
            #   '013_v14',
            #   '016_v12',

              '140_v03', #2022
              '141_v04',
              '141_v05',
              '142_v06',
              '142_v07',
              '184_v08',
              '185_v09',
              '185_v10',
              '186_v11',
              '228_v13',
              '229_v14',
              '271_v18',
              '273_v16',
              '274_v17',
              '309_v20',
              '310_v19',
              '311_v21',
              '348_v22',
              '349_v23',
              '349_v24',
    ]
        
    # visits = ['309_v20', '310_v19']
    # yrsearchstr = '202[12]'
    yrsearchstr = '2022'
   

    ajs =  np.array(pjaj.getAJTimes())
    ajt = Time.strptime(ajs, '%Y-%m-%d %H:%M:%S')
        
    rseq = 60268E3
    rspl = 54364E3

    maxint = 2000

    hisatstr = '_hisat'


    parallel = True

    for v in visits:

        files = sorted(np.array(glob.glob('/Users/jdn4/data/HST/jupiter/'+yrsearchstr+'/extract/' + v + '/*_0030_*2.fits')))
        piles = sorted(np.array(glob.glob('/Users/jdn4/data/HST/jupiter/'+yrsearchstr+'/extract/' + v + '/*_0030_*_flatproj.fits')))

        yrstr = files[0][29:33]
            
        nf = len(files)
        # nf = 8
        iarr = list(range(nf))

        fim = hst.NewHSTImage(files[nf//2])
        fim.readHSTFile()
        cml0 = fim.alm.cml+10
        ud = fim.alm.udate
        tmid = dt.strptime(ud, '%Y-%m-%d %H:%M:%S')

        fim = hst.NewHSTImage(files[0])
        fim.readHSTFile()
        ud = fim.alm.udate
        tstart = dt.strptime(ud, '%Y-%m-%d %H:%M:%S')

        fim = hst.NewHSTImage(files[-1])
        fim.readHSTFile()
        ud = fim.alm.udate
        tend = dt.strptime(ud, '%Y-%m-%d %H:%M:%S')

        
        gs = gridspec.GridSpec(2, 1, height_ratios=[1,3])
        gs.update(top=0.87, hspace = 0.00001, bottom=0.05, left=0.05, right=0.95)

        # reflon_list = ['_siii', '_cmls']
        reflon_list = ['_siii']

        for reflonstr in reflon_list:

            with warnings.catch_warnings():
                warnings.simplefilter("ignore")            

                if parallel is True:
                    pool = Pool()
                    pool.map(createMovieFrame, iarr)
                    pool.close()
                    pool.join()
                else:
                    for i in iarr:
                        createMovieFrame(i)
                        # import pdb; pdb.set_trace()
                        
            
            fim = hst.NewHSTImage(files[0])
            fim.readHSTFile()
            ud = fim.alm.udate
            t = dt.strptime(ud, '%Y-%m-%d %H:%M:%S')
            datestr = t.date().isoformat()
            doy = t.utctimetuple().tm_yday
            doystr = str(doy)
            timestr = t.time().isoformat()
            shtyr = datestr[2:4]
            yrstr = datestr[:4]
            sfx = 'jup_'+shtyr+'_'+doystr+'_'+timestr[0:2]+timestr[3:5]+'_v'+fim.alm.visit[:2].lower()
            path = '/Users/jdn4/data/HST/jupiter/'+yrstr+'/movies/'
        
            print('Creating projection movie in file...')
            outfile = path+sfx+'_nobg'+hisatstr+reflonstr+'.mp4'
            # outfile = path+sfx+'_p180_hisat.mp4'
            print(outfile)
            mp4str = 'ffmpeg -y -loglevel -8 -framerate 18 -i "'+path+'ims_%03d.jpg" -pix_fmt yuv420p -b:v 7M -vf "fps=18" '+outfile
            try:
                retcode = subprocess.call(mp4str, shell=True)
            except OSError as e:
                print("Execution failed:", e)


            deleteims = True
            if deleteims is True:
                delstr = '/bin/rm '+path+'ims_*.jpg'
                try:
                    subprocess.call(delstr, shell=True)
                except OSError as e:
                    print("Execution failed:", e)

            openstr = 'open '+outfile
            subprocess.call(openstr, shell=True)
            
    
