#!/usr/bin/env python
# encoding: utf-8
"""
hstprojimage.py

Version 1.0 Created by jdn on 2023-03-24.
Version 1.1 Modified by jdn on 2023-03-24 to run with Cartopy 0.21
Copyright (c) 2013 University of Leicester. All rights reserved.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize, LinearSegmentedColormap
from numba import njit
from pathlib import Path
import scipy.constants as c
from datetime import datetime as dt
from scipy.io.idl import AttrDict
import astropy.io.fits as pf
import cartopy.crs as ccrs
from cartopy.img_transform import warp_array
from cartopy.mpl.ticker import LongitudeFormatter
import warnings
import copy
import sys
import spiceypy as spice


class PlanetLonFormatter(LongitudeFormatter):
    """Object that modifies the Cartopy LongitudeFormatter to formats the
       longitude labels in the customary manner for planetary observations,
       i.e. from 0-360 deg, rather than +/- 180 E/W """

    def __init__(self, *args, **kwargs):
        super(PlanetLonFormatter, self).__init__(*args, **kwargs)

    def _format_value(self, value, original_value):
        number_format = self._degrees_number_format
        value = f"{(value+360)%360:{number_format}}{self._degree_symbol}"
        return value


class HSTProjImage(object):

    # Define reference file directory
    # This will be different for each setup
    # for Macbook Pro
    BASEDIR = Path('/Users/shin/Documents/Research/Jupiter/Codes/HST/')
    # for Macbook Air
    BASEDIR = Path('/Users/satoshin/Documents/Research/Jupiter/Codes/HST/')
    if BASEDIR.exists() is False:
        print('Please set the base directory to an existing directory')
    REFDIR = BASEDIR.joinpath('ref')
    if REFDIR.exists() is False:
        print('Please set the reference file directory to an existing directory')

    def __init__(self, filename):
        """
        Object for a projected HST auroral image. Primary functions are to read
        in and plot the image.

        Parameters
        ----------
        filename
            A string pointing to the projection FITS file
        """

        # Check filename is a string
        if isinstance(filename, str) is False:
            print('filename must be a string. Returning...')
            return

        # Check that the file is a projection
        check = filename.find('proj')
        if check == -1:
            print('Warning, file is not a projection!')

        self.filename = filename

        # This is because Cartopy raises a bunch of irritating warnings.
        # Comment out if preferred
        warnings.filterwarnings("ignore")

    # Convert the FITS header to an AttrDict for ease of element access
    def h2Alm(self, h):

        self.alm = AttrDict()
        for i in range(len(h)):
            self.alm.__setitem__(list(h.keys())[i], list(h.values())[i])

    def readHSTFile(self, ext=1, hext=0, cr=2.5, domu=True):
        """
        Reads in FITS file and header.
        Converts the almanac in the FITS header to AttrDict
        Sets kR values based on defined colour ratio
        Inflates the image from 1D to 2D array if file is *_flatproj

        Parameters
        ----------
        ext: optional
            The FITS extension holding the array to be read. Usually 1
        hext: optional
           The FITS extension holding the header to be read. Usually 0
        cr: optional
            The value of the colour ratio used when converting counts to kR
            Defaults to 2.5
            The values used are those in Table 1 of Gustin et al. (2012).
        domu: bool
            When True, calculate the cosine of the observation angle, mu.
            Usually needed, so defaults to True.
        """

        with pf.open(self.filename) as hdul:
            image = hdul[ext].data.astype(np.float32)
            h = hdul[hext].header

        self.h2Alm(h)
        try:
            self.alm.deltarpkm = self.alm.delrpkm
            self.alm.lightime = self.alm.dist_org*c.au/c.c
        except Exception:
            pass

        # Convert to new Gustin kR values at this stage to make sure it's done
        # The values in Table 1 of Gustin et al. (2012) are the reciprocals of the
        # equivalent values previously used. These are given in the arrays below

        if cr != 2.5:

            colr = [1.04, 1.10, 1.50, 2.00, 2.50, 3.00, 3.50, 4.00, 4.50, 5.00,
                    6.00, 7.00, 8.00, 9.00, 10.00, 12.00, 14.00, 16.00, 18.00,
                    20.00, 25.00]
            f1 = [469., 488, 596, 701, 789, 852, 908, 964, 1016, 1043, 1097,
                  1151, 1205, 1245, 1259, 1288, 1317, 1346, 1375, 1403, 1476]
            f2 = [815., 835, 950, 1049, 1127, 1176, 1218, 1261, 1300, 1318,
                  1355, 1391, 1427, 1453, 1462, 1480, 1498, 1516, 1534, 1552,
                  1597]
            cl = [1027., 1072, 1335, 1602, 1833, 2002, 2157, 2313, 2458, 2538,
                  2698, 2859, 3019, 3137, 3183, 3274, 3366, 3458, 3550, 3641,
                  3871]
            sr = [3948., 3994, 4215, 4391, 4523, 4605, 4675, 4746, 4810, 4841,
                  4902, 4963, 5024, 5069, 5086, 5120, 5155, 5189, 5224, 5258,
                  5344]
            # oldfactor = {'F115LP': 2.103E-3,'F125LP':1.473E-3}[self.alm.filter]
            oldfactor = 1./self.alm.cts2kr
            newfactors = {'F115LP': f1, 'F125LP': f2,
                          'CLEAR': cl, 'F25SRF2': sr}[self.alm.filter]
            newfactor = newfactors[colr.index(cr)]
            if np.round(newfactor) != np.round(oldfactor):
                image *= (newfactor/oldfactor)
                self.alm.cts2kr = 1./newfactor

        # Inflate the image if 1D
        @njit
        def _inflate_im(flatim, iinx, jinx, nx, ny):
            n = len(flatim)
            im = np.zeros((ny, nx), dtype=np.float32)
            for i in range(n):
                im[jinx[i], iinx[i]] = flatim[i]
            return im

        if image.ndim == 1:

            dat = np.load(HSTProjImage.REFDIR.joinpath(self.alm.ijinxfle))
            iinx = dat['iinx']
            jinx = dat['jinx']
            ny, nx = dat['shape']
            image = _inflate_im(image, iinx, jinx, nx, ny)

        self.image = image
        self.ny, self.nx = np.shape(self.image)

        # Save image time as a DateTime object
        self.t = dt.strptime(self.alm.udate, '%Y-%m-%d %H:%M:%S')

        # Define the coordinate boundaries based on the image shape and hemisphere
        # x and y are logical coordinates based on the full 1440 x 720 array (0.25 deg resolution)
        if self.alm.hemisph == 'north':
            if self.ny == 180:
                y1, y2, x1, x2 = 540, 720, 0, 1440  # 45 -> 89.75 deg lat
                lon1, lon2, lat1, lat2 = 360, 0, 45, 90
            elif self.ny == 200:
                y1, y2, x1, x2 = 520, 720, 0, 1440  # 40 -> 89.75 deg lat
                lon1, lon2, lat1, lat2 = 360, 0, 40, 90
            elif self.ny == 720:
                y1, y2, x1, x2 = 0, 720, 0, 1440    # 0 -> 89.75 deg lat
                lon1, lon2, lat1, lat2 = 360, 0, -90, 90
            else:
                print('Unrecognised cylindrical projection array shape:',
                      self.nx, self.ny)
                sys.exit()

        elif self.alm.hemisph == 'south':
            if self.ny == 181:
                y1, y2, x1, x2 = 0, 181, 0, 1440    # 0 -> -45 deg lat
                # longitudes flipped so aurora plotted through planet
                lon1, lon2, lat1, lat2 = 360, 0, -90, -45
            elif self.ny == 720:
                y1, y2, x1, x2 = 0, 720, 0, 1440    # 0 -> 89.75 deg lat
                lon1, lon2, lat1, lat2 = 360, 0, -90, 90
            else:
                print('Unrecognised cylindrical projection array shape:',
                      self.nx, self.ny)
                sys.exit()

        else:
            if self.ny == 720:
                y1, y2, x1, x2 = 0, 720, 0, 1440    # 0 -> 89.75 deg lat
                lon1, lon2, lat1, lat2 = 360, 0, -90, 90
            else:
                print('Unrecognised cylindrical projection array shape:',
                      self.nx, self.ny)
                sys.exit()

        self.y1, self.y2, self.x1, self.x2 = y1, y2, x1, x2
        self.cylproj_extent = lon1, lon2, lat1, lat2

        # Instantiate the Cartopy cylindrical projection
        # Planetary radius set to 1 R_p. Oblateness was taken care of in reduction
        self.globe = ccrs.Globe(
            semimajor_axis=1, semiminor_axis=1, ellipse=None)
        self.cylproj = ccrs.PlateCarree(central_longitude=0, globe=self.globe)
        self.geodetic = ccrs.Geodetic(self.globe)

        # Compute mu = the cosine of the observation angle

        @njit
        def _mu(y1, y2, x1, x2, selat, selon):
            mu = np.zeros((y2-y1, x2-x1))
            d2r = np.pi/180.
            sinselt = np.sin(selat*d2r)
            cosselt = np.cos(selat*d2r)

            for j in range(y1, y2):
                lat = (j/4 - 90.) * d2r
                sinlat = np.sin(lat)
                coslat = np.cos(lat)
                for i in range(x1, x2):
                    lon = (360. - i/4) * d2r
                    cosdellon = np.cos(lon - selon*d2r)
                    mu[j - y1, i] = sinlat*sinselt + coslat*cosselt*cosdellon
            return mu
        if domu is True:
            self.mu = _mu(y1, y2, x1, x2, self.alm.dece, self.alm.cml)

    def _getGistCool(self):
        """Returns a colour map based on gist_heat,
        but with the R and B channels swapped"""

        gh = copy.deepcopy(plt.cm.gist_heat)
        ghd = gh._segmentdata
        ghd['blue'], ghd['red'] = ghd['red'], ghd['blue']
        coolmap = LinearSegmentedColormap('gist_cool', ghd, 256)

        return coolmap

    # Simple limb brightening correction by multiplication by the cosine of the observation angle mu
    def limbBrighteningCorrection(self, image=None):
        if image is None:
            image = self.image
        image = image[:] * clip(self.mu, 0., 1.)
        return image

    # Set to zero if the observation angle is larger than a given value
    def clipLimb(self, image=None, yclip=89., val=0.):
        if image is None:
            image = self.image
        image = np.where(self.mu > np.cos(np.deg2rad(yclip)), image, val)
        return image

###############################################################################

    def tvPolar(self, axis=None, vmin=None, vmax=None, yclip=88, yclipval=None,
                norm='log', grid=True, reflon=None, zero=True, badcol='0.',
                cmap=None, satovals=['all'], refmainoval=False, noplot=False,
                draw_labels=False, grid_spacing=(20, 10), ylabellon=220, **kwargs):
        """
        Converts the cylindrical projection of the image to
        Polar Stereographic projections using Cartopy and plots.
        The southern aurora is plotted as if looking through the planet.
        Parameters control some aspects of the plot

        Parameters
        ----------
        axis: optional
            A class:`matplotlib.axes.Axes` instance, used to define the position
            of the  :class:`cartopy.mpl.geoaxes.GeoAxes` object that is created
            to plot the image. axis is destroyed after its position is used.
        vmin, vmax: float, optional
           Defines the data range the colormap covers. Defaults to the full
           range of the data
        yclip: float, optional
            The maximum value of the observation angle in degrees to be shown.
            Used to clip stretched pixels near the limb. The default is 88 deg.
        yclipval: float, optional
            The values that pixels beyond the clip limit are set to. Defaults to
            vmin
        norm: str, optional
            The normalisation scheme used for the colormap.
                - 'lin': Linear
                - 'log': Logarithmic'
        grid: bool, optional
            When True, draw gridlines using :class`cartopy.mpl.gridliner.Gridliner`
            Defaults to True.
        reflon: float, optional
            Defines the longitude oriented toward the bottom of the image.
            Defaults to 180 for the north, and 0 for the south.
        zero: bool, optional
            When True, draw the prime meridian in red dotted lines
        badcol: optional
            An object that can be interpreted as a colour by matplotlib.colors
        cmap: 'str' or :class:`matplotlib.colors.Colormap`
            The :class:`matplotlib.colormap` instance or registered colormap
            name used to map scalar data to colors.
            Defaults to "gist_cool", which is like gist_heat but with
            the blue and red channels swapped
        satovals: list, optional
            List of strings that sets which (if any) JRM33 satellite contours
            are plotted. Defaults to Io, Europa, Ganymede
                - ['io'] Just Io
                - ['io', 'eu', 'ga'] Io, Europa, Ganymede
                - ['all'] Io, Europa, Ganymede
        refmainoval: bool, optional
            When True, plot the reference main oval from Nichols et al. (2017)
            in red
        noplot: bool, optional
            When True, return before plotting the projected image
        draw_labels: bool, optional
            Passed to :class`cartopy.mpl.gridliner.Gridliner`. When True, draw
            grid line labels. Defaults to False
        grid_spacing: tuple, optional
            A 2-tuple (lon_spacing, lat_spacing), which defines the spacing
            in degrees of drawn gridlines.
        ylabellon: float, optional
            The longitude at which northern inline latitude labels are shown
        """
        # Position of the Galilean moons in the "IAU_JUPITER" frame
        # Calculating position of the moon by spiceypy.
        # Then converting to the System III longitude.
        spice.furnsh("cassMetaK.txt")
        utc = '2022-05-22T10:30:41'
        et = spice.str2et(utc)

        # Europa's position seen from the HST in IAU_JUPITER coordinate.
        pos, lightTimes = spice.spkpos(
            targ='EUROPA', et=et, ref='IAU_JUPITER', abcorr='LT+S', obs='HST'
        )

        # Jupiter's position seen from the HST in IAU_JUPITER coordinate.
        pos1, lightTimes = spice.spkpos(
            targ='JUPITER', et=et, ref='IAU_JUPITER', abcorr='LT+S', obs='HST'
        )

        # Europa's position seen from Jupiter in IAU_JUPITER coordinate.
        pos = pos - pos1

        posx, posy = pos[0], pos[1]
        posz = pos[2]
        posr = np.sqrt(posx**2 + posy**2 + posz**2)
        # postheta = np.arccos(posz/posr)
        posphi = np.arctan2(posy, posx)
        if posphi < 0:
            Sys3 = np.degrees(-posphi)
        else:
            Sys3 = np.degrees(2*np.pi - posphi)
        print('R', posr/71492, 'PHI', posphi, 'SYS3', Sys3)

        s3list = np.arange(0, 360, 5, dtype=int)
        s3_idx = np.argmin(np.abs(np.arange(0, 360, 5)-Sys3))
        print('S3_IDX', s3_idx)
        # s3_candidate = s3list[s3_idx]   # Closest Sys3 long.

        # Output projection
        hem = self.alm.hemisph
        if hem == 'north':
            if reflon is None:
                reflon = 180
            outproj = ccrs.NorthPolarStereo(
                central_longitude=reflon, globe=self.globe)
        else:
            if reflon is None:
                reflon = 0
            outproj = ccrs.SouthPolarStereo(
                central_longitude=reflon + 180, globe=self.globe)

        # Copy the image array for a version to operate on for display
        dimage = self.image.copy()

        # Colour scale limits
        if vmin is None:
            vmin = np.nanmin(dimage)
        if vmax is None:
            vmax = np.nanmax(dimage)

        dimage = np.clip(dimage, vmin, vmax)

        # Load the colour table and set bad point colour
        if cmap is None:
            cmap = self._getGistCool()
            cmap.set_bad(badcol, 1.)

        # Set values near the limb to zero if required
        if yclip is not None:
            if yclipval is None:
                yclipval = vmin
            dimage = self.clipLimb(image=dimage, yclip=yclip, val=yclipval)

        # This creates a new axis - this needs to be done here as its type is
        # determined by the output projection. Copy the position of the input axis
        if axis is None:
            axis = plt.gca()
        ax = plt.gcf().add_axes(axis.get_position(), projection=outproj)
        axis.remove()
        del axis
        ax.set_aspect('equal', adjustable='box')

        # Set the plot extents in the output coordinate system
        def _shift_centre(aplat, aplon, out_extent):
            # aplat and aplon are the coords of approx. centroid of auroral region
            apx, apy = outproj.transform_point(aplon, aplat, self.cylproj)
            midx = out_extent[0] + (out_extent[1] - out_extent[0])/2.
            midy = out_extent[2] + (out_extent[3] - out_extent[2])/2.
            shiftx = apx - midx
            shifty = apy - midy
            out_extent += np.array([shiftx, shiftx, shifty, shifty])
            return out_extent

        if hem == 'north':
            out_extent = np.array([-0.51, 0.49, -0.77, 0.23])
            if reflon != 180:
                out_extent = _shift_centre(70, 180, out_extent)
        else:
            out_extent = np.array([-0.60, 0.4, -0.6, 0.4])
            if reflon != 0:
                out_extent = _shift_centre(-82, 42, out_extent)

        # Perform the image transformation
        dimage, dum = warp_array(dimage, outproj, source_proj=self.cylproj,
                                 target_res=(500, 500), source_extent=self.cylproj_extent,
                                 target_extent=out_extent)

        # Return if no plotting is required
        if noplot is True:
            return ax, dimage

        # Show the image and plot gridlines aand tick labels if required
        self.normfunc = {'log': LogNorm(vmin=vmin, vmax=vmax),
                         'lin': Normalize(vmin=vmin, vmax=vmax)}[norm]
        self.tvim = ax.imshow(dimage, origin='lower',
                              extent=out_extent, norm=self.normfunc, cmap=cmap, zorder=0.8)

        if grid is True:
            gl = ax.gridlines(crs=self.cylproj, xlocs=np.arange(-180, 180, grid_spacing[0]),
                              ylocs=np.arange(-90, 90, grid_spacing[1]),
                              ylim=(-80, 80), linestyle=':', draw_labels=draw_labels, y_inline=True, x_inline=False,
                              xformatter=PlanetLonFormatter(direction_label=False))

        else:
            gl = []

        # Load and plot the reference main oval
        if refmainoval is True:
            if hem == 'north':
                mofile = HSTProjImage.REFDIR.joinpath('refmon.npz')
                latmod = 1.
            else:
                mofile = HSTProjImage.REFDIR.joinpath('refmoshires.npz')
                latmod = -1.
            npzfile = np.load(mofile)
            molat = npzfile['hilat'] * latmod
            molon = npzfile['hilon']
            ax.plot(molon, molat, 'r', transform=self.geodetic)

        if zero is True:
            ax.plot([0, 0], [-80, 80], 'r:', transform=self.geodetic)

        # Load and plot the satellite contours
        if satovals == ['all']:
            satovals = ['io', 'eu', 'ga']
        if len(satovals) > 0:
            num = {'north': 2, 'south': 3}[hem]
            satoval = np.recfromtxt(HSTProjImage.REFDIR.joinpath(
                '2021je007055-sup-000'+str(1+num)+'-table si-s0'+str(num)+'.txt'), skip_header=3,
                names=['wlon', 'amlat', 'amwlon', 'iolat', 'iowlon', 'eulat', 'euwlon', 'galat', 'gawlon'])
            euftp_lon = (satoval.euwlon[s3_idx-1]+satoval.euwlon[s3_idx])/2
            euftp_lat = (satoval.eulat[s3_idx-1]+satoval.eulat[s3_idx])/2
            print('EUROPA FTP', euftp_lon, euftp_lat)
            if 'io' in satovals:
                ax.plot(satoval.iowlon, satoval.iolat, 'w',
                        transform=self.geodetic, lw=0.4)
            if 'eu' in satovals:
                ax.plot(satoval.euwlon, satoval.eulat, 'w',
                        transform=self.geodetic, lw=0.4, zorder=0.8)
                ax.plot(euftp_lon, euftp_lat, markersize=18, marker='o', markerfacecolor='none', markeredgecolor='red',
                        transform=self.geodetic, zorder=5)
            if 'ga' in satovals:
                ax.plot(satoval.gawlon, satoval.galat, 'w',
                        transform=self.geodetic, lw=0.4)

        ax.set_extent(out_extent, crs=outproj)
        if hem == 'north':
            ax.invert_xaxis()

        # Finally draw labels if required
        if draw_labels is True:
            gl.xlabel_style = {'c': 'w', 'size': 6}
            gl.ylabel_style = {'c': 'w', 'size': 6}
            gl.rotate_labels = False
            gl.xpadding = -4
            gl.geo_labels = False
            gl.right_labels = True  # Gets overridden by Gridliner for some reason
            gl.left_label = True    # Ditto
            plt.gcf().canvas.draw()
            for i, lbl in enumerate(gl._labels):
                tx, ty = lbl.artist.get_position()
                r = np.sqrt(tx**2 + ty**2)
                if lbl.loc == 'right' or lbl.loc == 'left':  # Manually set
                    lbl.artist.set_visible(True)
                if hem == 'north':
                    if lbl.xy == 'y' and lbl.loc == 'inline':
                        lbl.artist.set_position((ylabellon, ty))
                    gl.right_labels = True
        self.gl = gl

        return ax, dimage, gl

    def tvProj(self, axis=None, vmin=None, vmax=None, norm='log', grid=True, reflon=180,
               badcol='0.', satovals=['all'], cmap=None, refmainoval=False,
               draw_labels=False, grid_spacing=(20, 10), fulldisc=True, **kwargs):
        """
        Plots the cylindrical projection image as a PlateCarree projection

        Parameters
        ----------
        axis: optional
            A class:`matplotlib.axes.Axes` instance, used to define the position
            of the  :class:`cartopy.mpl.geoaxes.GeoAxes` object that is created
            to plot the image. This is destroyed after its position is used.
        vmin, vmax: float, optional
           Defines the data range the colormap covers. Defaults to the full
           range of the data
        norm: str, optional
            The normalisation scheme used for the colormap.
                - 'lin': Linear
                - 'log': Logarithmic'
        grid: bool, optional
            When True, draw gridlines using :class`cartopy.mpl.gridliner.Gridliner`
            Defaults to True.
        reflon: float, optional
            Defines the longitude oriented toward the bottom of the image.
            Defaults to 180 for the north, and 0 for the south.
        badcol: optional
            An object that can be interpreted as a colour by matplotlib.colors
        cmap: 'str' or :class:`matplotlib.colors.Colormap`
            The :class:`matplotlib.colormap` instance or registered colormap
            name used to map scalar data to colors.
            Defaults to "gist_cool", which is like gist_heat but with
            the blue and red channels swapped
        refmainoval: bool, optional
            When True, plot the reference main oval from Nichols et al. (2017)
            in red
        draw_labels: bool, optional
            Passed to :class`cartopy.mpl.gridliner.Gridliner`. When True, draw
            grid line labels. Defaults to False
        grid_spacing: tuple, optional
            A 2-tuple (lon_spacing, lat_spacing), which defines the spacing
            in degrees of drawn gridlines.
        fulldisc: bool, optional
            When True, plot the full disc, i.e. 0->360 lon, -90->90 lat,
            otherwise plot only the region covered by the image
        """

        # Output projection
        hem = self.alm.hemisph
        outproj = ccrs.PlateCarree(central_longitude=reflon, globe=self.globe)

        # Copy the image array for a version to operate on for display
        dimage = self.image.copy()

        # Colour scale limits
        if vmin is None:
            vmin = np.nanmin(dimage)
        if vmax is None:
            vmax = np.nanmax(dimage)

        dimage = np.clip(dimage, vmin, vmax)

        # Load the colour table and set bad point colour
        if cmap is None:
            cmap = self._getGistCool()
            cmap.set_bad(badcol, 1.)

        # This creates a new axis - this needs to be done here as its type is
        # determined by the output projection. Copy the position of the input axis
        if axis is None:
            axis = plt.gca()
        ax = plt.gcf().add_axes(axis.get_position(), projection=outproj)
        axis.remove()
        del axis
        ax.set_aspect('auto', adjustable='box')

        out_extent = [-180, 180, self.cylproj_extent[2],
                      self.cylproj_extent[3]]

        # Perform the image transformation
        dimage, dum = warp_array(dimage, outproj, source_proj=self.cylproj,
                                 target_res=(
                                     self.nx, self.ny), source_extent=self.cylproj_extent,
                                 target_extent=out_extent)

        # Emplace the image array in one covering the full disc if required
        if fulldisc is True:
            temp = np.zeros((720, 1440))
            temp[self.y1:self.y2, self.x1:self.x2] = dimage
            dimage = temp
            out_extent = [-180, 180, -90, 90]

        # Show the image and plot gridlines and labels if required
        self.normfunc = {'log': LogNorm(vmin=vmin, vmax=vmax),
                         'lin': Normalize(vmin=vmin, vmax=vmax)}[norm]
        self.tvim = ax.imshow(dimage, origin='lower',
                              extent=out_extent, norm=self.normfunc, cmap=cmap)
        ax.set_aspect('equal', adjustable='box')
        if grid is True:
            gl = ax.gridlines(crs=self.cylproj, xlocs=np.arange(-180, 180, grid_spacing[0]),
                              ylocs=np.arange(-90, 90, grid_spacing[1]),
                              ylim=(-90, 90), linestyle=':', draw_labels=draw_labels,
                              xformatter=PlanetLonFormatter(direction_label=False))
            gl.rotate_labels = False
            gl.xlabel_style = {'c': 'k', 'size': 'xx-small'}
            gl.ylabel_style = {'c': 'k', 'size': 'xx-small'}
        else:
            gl = []
        self.gl = gl

        # Load and plot the reference main oval
        if refmainoval is True:
            if hem == 'north':
                mofile = HSTProjImage.REFDIR.joinpath('refmon.npz')
                latmod = 1.
            else:
                mofile = HSTProjImage.REFDIR.joinpath('refmoshires.npz')
                latmod = -1.
            npzfile = np.load(mofile)
            molat = npzfile['hilat'] * latmod
            molon = npzfile['hilon']
            ax.plot(molon, molat, 'r', transform=self.geodetic)

        # Load and plot the satellite contours
        if satovals == ['all']:
            satovals = ['io', 'eu', 'ga']
        if len(satovals) > 0:
            num = {'north': 2, 'south': 3}[hem]
            satoval = np.recfromtxt(HSTProjImage.REFDIR.joinpath(
                '2021je007055-sup-000'+str(1+num)+'-table si-s0'+str(num)+'.txt'), skip_header=3,
                names=['wlon', 'amlat', 'amwlon', 'iolat', 'iowlon', 'eulat', 'euwlon', 'galat', 'gawlon'])
            if 'io' in satovals:
                ax.plot(satoval.iowlon, satoval.iolat, 'w',
                        transform=self.geodetic, lw=0.4)
            if 'eu' in satovals:
                ax.plot(satoval.euwlon, satoval.eulat, 'w',
                        transform=self.geodetic, lw=0.4)
            if 'ga' in satovals:
                ax.plot(satoval.gawlon, satoval.galat, 'w',
                        transform=self.geodetic, lw=0.4)

        ax.set_extent(out_extent, crs=outproj)
        ax.invert_xaxis()

        return ax, dimage, gl


"""
# plt.ion()
plt.figure(1)
plt.clf()
EXTRACTDIR = 'data/sample/'
h = HSTProjImage(
    EXTRACTDIR+'jup_22-142-10-30-41_0030_v06_stis_f25srf2_flatproj.fits')
# h = HSTProjImage(EXTRACTDIR+'140_v03/jup_22-140-15-38-53_0030_v03_stis_f25srf2_flatproj.fits')
h.readHSTFile()
ax = plt.subplot(111)
h.tvPolar(ax, vmin=10, vmax=2000,  draw_labels=True,
          yclip=88, refmainoval=True, reflon=180)
plt.show()
"""
