""" F25SR2F.py

F25SR2Fバンドパスフィルターのthroughputを計算するプログラム

"""

import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
import matplotlib.patches as patches
from matplotlib.ticker import AutoMinorLocator
from astropy import units as u
import spiceypy as spice
spice.furnsh('kernel/cassMetaK.txt')


fontname = 'Nimbus Sans'
plt.rcParams.update({'font.sans-serif': fontname,
                    'font.family': 'sans-serif',
                     'mathtext.fontset': 'custom',
                     'mathtext.rm': fontname,
                     'mathtext.it': fontname+':italic',
                     # 'mathtext.bf': 'Nimbus Sans:italic:bold',
                     'mathtext.bf': fontname+':bold',
                     'text.kerning_factor': 0,
                     # 'font.size': 19,
                     })
params = {
    # 'lines.markersize': 1,
    # 'lines.linewidth': 1,
    'axes.linewidth': 2,
    'xtick.major.size': 5,
    'xtick.minor.size': 3.5,
    'xtick.major.width': 2.0,
    'xtick.minor.width': 1.25,
    'ytick.major.size': 5,
    'ytick.minor.size': 3,
    'ytick.major.width': 2.0,
    'ytick.minor.width': 1.25,
}
plt.rcParams.update(params)

# Color universal design
cud4 = ['#FF3300', '#FFF100', '#03AF7A', '#005AFF',
        '#4DC4FF', '#FF8082', '#F6AA00', '#990099', '#804000']
cud4bs = ['#FFCABF', '#FFFF80', '#D8F255', '#BFE4FF',
          '#FFCA80', '#77D9A8', '#C9ACE6', '#84919E']


# Constance
A_HST = 45239           # Area of HST's main mirror [cm2]
Sigma_px = 0.0246       # Plate scale of STIS/FUV-MAMA [arcsec/px]
Sigma_px *= 1/3600      # [deg/px]
Sigma_px *= np.pi/180   # [rad/px]
Sigma_px = Sigma_px**2  # [str/px]


# Set legend shadow
def legend_shadow(fig, ax, legend, dx, dy):

    frame = legend.get_window_extent()

    xmin, ymin = fig.transFigure.inverted().transform((frame.xmin, frame.ymin))
    xmax, ymax = fig.transFigure.inverted().transform((frame.xmax, frame.ymax))

    # plot patch shadow
    rect = patches.Rectangle((xmin+dx, ymin+dy), xmax-xmin, ymax-ymin,
                             transform=fig.transFigure,
                             edgecolor='k', facecolor='k',
                             clip_on=False)
    ax.add_patch(rect)

    return None


# Data load (Brightness)
def load(csvname0):
    df = pd.read_csv(csvname0, sep='\t')
    utc = df.loc[:, 'date']                                 # UTC date
    # EFP latitude [deg]
    efplat = df.loc[:, 'lat [deg]'].values
    # EFP System III longitude [deg]
    efpwlong = df.loc[:, 'wlong [deg]'].values
    # Europa's System III longitude [deg]
    moons3 = df.loc[:, 'Moon S3 [deg]'].values
    # EFP brightness [kR]
    final_phot_ave = df.loc[:, 'spot brightness [kR]']
    # local background [kR]
    annulus_median = df.loc[:, 'local background [kR]']

    b0_arr = np.zeros(len(final_phot_ave))
    b1_arr = np.zeros(len(final_phot_ave))
    efplat0_arr = np.zeros(len(final_phot_ave))
    efpwlong0_arr = np.zeros(len(final_phot_ave))
    moons30_arr = np.zeros(len(final_phot_ave))
    for i in range(len(final_phot_ave)):
        if final_phot_ave[i] != '0':
            b0_arr[i] = final_phot_ave[i]
            b1_arr[i] = annulus_median[i]
            efplat0_arr[i] = efplat[i]
            efpwlong0_arr[i] = efpwlong[i]
            moons30_arr[i] = moons3[i]
        else:
            continue

    return utc, b0_arr, b1_arr, efplat0_arr, efpwlong0_arr, moons30_arr


# Data load (STIS throughput)
a = np.loadtxt('data/F25SRF2/HST_STIS_FUV.F25SRF2.dat')
wl_data = a[:, 0]/10    # [nm]
th_data = a[:, 1]*100   # [%]


# Effective range of wavelength
wl_e0 = 128       # [nm]
wl_e1 = 178+0.5   # [nm]


def WL_PIVOT(wl_e0: float, wl_e1: float,):
    """
    Args:
        wl_e0 (float): min wavelength [nm]
        wl_e1 (float): max wavelength [nm]
    """
    wl_erange = np.where(((wl_data >= wl_e0) & (wl_data <= wl_e1)))

    # Calculate an effective throughput
    wl_e = wl_data[wl_erange]
    th_e = th_data[wl_erange]
    # print(wl_e)

    # Pivot wavelength
    S2_u = 0
    S2_l = 0
    for i in range(wl_e.size-1):
        dwl = wl_e[i+1]-wl_e[i]
        S2_u += (th_e[i+1]+th_e[i])*dwl/2
        S2_l += ((th_e[i+1]/wl_e[i+1]**2)+(th_e[i]/wl_e[i]**2))*dwl/2

    wl_pivot = math.sqrt(S2_u/S2_l)

    return wl_pivot


def TH_LINEAR(wl_input: float, ):
    """
    Args:
        wl_input (float): input wavelength [nm]
    """
    # Calculate throughput
    argsorted = np.argsort(np.abs(wl_data-wl_input))
    idx0, idx1 = argsorted[0], argsorted[1]

    # Linearly interpolated
    dx = wl_data[idx1]-wl_data[idx0]
    dy = th_data[idx1]-th_data[idx0]
    th_linear = (dy/dx)*(wl_input-wl_data[idx0]) + th_data[idx0]

    return th_linear


def TH_EQUIV(wl_e0: float, wl_e1: float):
    """
    Args:
        wl_e0 (float): min wavelength [nm]
        wl_e1 (float): max wavelength [nm]
    """

    S = 0
    wl_erange = np.where(((wl_data >= wl_e0) & (wl_data <= wl_e1)))

    # Calculate an effective throughput
    wl_e = wl_data[wl_erange]
    th_e = th_data[wl_erange]
    S = 0
    for i in range(th_e.size-1):
        dwl = wl_e[i+1]-wl_e[i]
        S += 0.5*(th_e[i+1]+th_e[i])*dwl

    th_equiv = S/(wl_e[-2]-wl_e[0])
    print('EQUIV: ', th_equiv)

    return th_equiv


wl_pivot = WL_PIVOT(wl_e0, wl_e1)
th_pivot = TH_LINEAR(wl_pivot)
th_equiv = TH_EQUIV(wl_e0, wl_e1)
print('Pivot [nm]: ', wl_pivot)
print('Pivot throughpuut: ', th_pivot)
print('Equiv throughpuut: ', th_equiv)


# Plot (1)
fontsize = 19
fig, ax = plt.subplots(figsize=(5, 3.5), dpi=150)
fig.tight_layout()
ax.set_title('F25SRF2 Throughput', fontsize=fontsize, weight='bold')
ax.tick_params(axis='both', labelsize=fontsize)
ax.set_xlabel('Wavelength [nm]', fontsize=fontsize)
ax.set_ylabel('Throughput [%]', fontsize=fontsize)
# ax.set_yscale('log')
ax.set_xlim(100, 300)
ax.set_ylim(1E-6, 3)
ax.plot(wl_data, th_data, color='k')
ax.scatter(wl_pivot, th_pivot, color=cud4[3], marker='s', s=20)
ax.axhline(y=th_pivot, color=cud4[3], label='Pivot')
ax.axhline(y=th_equiv, color=cud4[0], label='Equiv')
ax.axvline(x=wl_e0, color='gray', linestyle='dashed')
ax.axvline(x=wl_e1, color='gray', linestyle='dashed')
legend = ax.legend(loc='upper center',
                   ncol=1,
                   markerscale=3.5,
                   bbox_to_anchor=(0.85, 1.0),
                   fancybox=False,
                   facecolor='white',
                   framealpha=1,
                   edgecolor='k',
                   fontsize=fontsize*0.75,
                   labelspacing=0.34,
                   handlelength=0.5,
                   scatterpoints=1, )
legend_shadow(fig, ax, legend, 0.0001, -0.078)
fig.tight_layout()
plt.show()


# Brightness data
north_doy = ['14/006_v06', '14/013_v13',
             '14/016_v12', '22/271_v18', '22/274_v17']
doy = north_doy
csvname0 = 'img/red3_half2/EUROPA/20'+doy[0]+'/brightness.csv'


def EP_CALC(HSTUTC, B_kR, CR):
    # HST's position seen from the HST in IAU_JUPITER coordinate.
    et_hst = spice.str2et(HSTUTC)
    pos, _ = spice.spkpos(
        targ='HST', et=et_hst, ref='IAU_JUPITER', abcorr='LT+S', obs='JUPITER'
    )
    r_hst = np.sqrt(pos[:, 0]**2+pos[:, 1]**2+pos[:, 2]**2)  # [km]
    # print('Distance [km] and [AU]: ', r_hst, r_hst/AU)

    # Count rate
    Ct_convert = 1/4523     # Gustin+2012 CR=2.5
    Ct = Ct_convert*B_kR
    # print(Ct)

    # Emitted power
    if CR == 1.5:
        EP_70180 = (4215/3948)*(9.54E-10)*(r_hst**2)*Ct    # [W]
    if CR == 2.0:
        EP_70180 = (4391/3948)*(9.94E-10)*(r_hst**2)*Ct    # [W]
    if CR == 2.5:
        EP_70180 = (4523/3948)*(1.02E-9)*(r_hst**2)*Ct    # [W]
    elif CR == 5.0:
        EP_70180 = (4841/3948)*(1.10E-9)*(r_hst**2)*Ct    # [W]
    elif CR == 10:
        EP_70180 = (5086/3948)*(1.15E-9)*(r_hst**2)*Ct    # [W]

    return EP_70180


def EP_CALC2(B_kR):
    # Assumed: 10 kR = 1 mW m-2 electron energy flux (Gustin+2012)
    # Assumed: radius of EFP = 300 km on Jupiter (Moirano+2021)
    Rf = 300*1E+3   # radius of EFP [m]
    EP_E = B_kR * 0.1 * np.pi * (Rf**2) * 1E-3      # [W]

    return EP_E


# Plot (3)
fontsize = 18
fig, ax = plt.subplots(figsize=(5.8, 3.4), dpi=326)
ax.set_zorder(1)
ax.patch.set_visible(False)  # prevents ax1 from hiding ax2
fig.tight_layout()
# ax.set_title('FUV Power (70-180 nm)', fontsize=fontsize, weight='bold')
ax.tick_params(axis='both', labelsize=fontsize)
ax.set_xlim(0, 360)
ax.set_ylim(0, 7)
ax.set_xticks(np.arange(0, 361, 45))
ax.set_xticklabels(np.arange(0, 361, 45), fontsize=fontsize)
ax.xaxis.set_minor_locator(AutoMinorLocator(3))  # minor ticks
ax.set_yticks(np.arange(0, 7, 2))
# ax.set_yticklabels(np.arange(0, 23, 5), fontsize=fontsize)
ax.yaxis.set_minor_locator(AutoMinorLocator(4))  # minor ticks
ax.set_xlabel(
    r'Moon System III longitude $\lambda_{\rm III}$ [deg]', fontsize=fontsize)
ax.set_ylabel('Total emission [GW]', fontsize=fontsize)
ax.text(0.01, 1.020, r'Total UV emission in 70-180 nm',
        color='k',
        horizontalalignment='left',
        verticalalignment='bottom',
        transform=ax.transAxes,
        fontsize=fontsize*0.42)
ax.text(0.035, 0.93, 'STIS/SrF2',
        weight='bold',
        color='k',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax.transAxes,
        fontsize=fontsize*0.8)
ax.text(0.70, 0.93, '2014',
        weight='bold',
        color=cud4[0],
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax.transAxes,
        fontsize=fontsize*0.8)
ax.text(0.85, 0.93, '2022',
        weight='bold',
        color=cud4[3],
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax.transAxes,
        fontsize=fontsize*0.8)

ax2 = ax.twinx()
ax2.tick_params(axis='both', labelsize=fontsize)
ax2.set_ylim(0, 450)
ax2.set_yticks(np.arange(0, 450, 100))
# ax2.set_yticklabels(np.arange(0, 101, 25), fontsize=fontsize)
ax2.yaxis.set_minor_locator(AutoMinorLocator(4))  # minor ticks
ax2.set_ylabel(r'Poynting flux $S$ [GW]', fontsize=fontsize)

# OBSERVED POWER
cud4_N = [cud4[0], cud4[2], cud4[3], cud4[5], cud4[7]]
cud4_N2 = [cud4[0], cud4[0], cud4[0], cud4[3], cud4[3]]
doyname = ['2014-01-06', '2014-01-13', '2014-01-16'] + \
    ['2022-09-28', '2022-10-01']
colorratio = [2.0, 2.0, 2.0, 2.0, 2.0]
for i in range(len(doy)):
    # print(doyname[i])
    csvname0 = 'img/red3_half2/EUROPA/20'+doy[i]+'/brightness.csv'
    utc, b0_arr, b1_arr, efplat0_arr, efpwlong0_arr, moons30_arr = load(
        csvname0)
    idx = np.where(b0_arr > 0)
    utc = utc[idx[0]]
    b0_arr = b0_arr[idx]       # [kR]
    b1_arr = b1_arr[idx]       # [kR]
    efplat0_arr = efplat0_arr[idx]
    efpwlong0_arr = efpwlong0_arr[idx]
    moons30_arr = moons30_arr[idx]

    # リード角補正
    lead_data = np.loadtxt('data/red3_leadangle/EUROPA/20'+doy[i]+'_eq.txt')
    eq_leadangle = lead_data[1, :]   # lead angle [deg]
    if doy[i] == '22/271_v18':
        # なぜかリード角のデータと1行目が合わない
        utc = utc[1:]
        efplat0_arr = efplat0_arr[1:]
        efpwlong0_arr = efpwlong0_arr[1:]
        moons30_arr = moons30_arr[1:]
        b0_arr = b0_arr[1:]
        b1_arr = b1_arr[1:]
    moons30_ave = np.average(moons30_arr)
    moons30_std = np.std(moons30_arr)
    b0_ave = np.average(b0_arr)
    b1_ave = np.average(b1_arr)
    b0_std = np.std(b0_arr)
    moon_leadback = moons30_arr-eq_leadangle
    moon_leadback_ave = np.average(moon_leadback)
    moon_leadback_std = np.std(moon_leadback)

    EP_70180 = EP_CALC(HSTUTC=utc, B_kR=b0_arr,
                       CR=colorratio[i])     # Emitted power [W]
    EP_70180 = EP_CALC2(B_kR=b0_arr)
    EP_std = np.std(EP_70180)           # Standard deviation [W]
    print('Ave. emission [GW]', np.average(EP_70180)*1E-9)

    """
    ax.scatter(moon_leadback, EP_70180*1E-6, marker='s', s=1,
               color=cud4_N[i], label=doyname[i], zorder=10-i, )
    """
    ax.errorbar(moon_leadback_ave, np.average(EP_70180)*1E-9,
                xerr=moon_leadback_std, yerr=EP_std*1E-9,
                marker='s', markersize=3.5, mfc=cud4_N2[i],
                mec=cud4_N2[i], linestyle='none', ecolor=cud4_N2[i],
                elinewidth=1.2, zorder=1.5)


# ESTIMATED POWER
pwr_14 = np.loadtxt('data/Poyntingflux/PY_2014_R4.txt')
pwr_22 = np.loadtxt('data/Poyntingflux/PY_2022_R4.txt')

pwr_14[1, :] *= 1E-9    # [GW]
pwr_22[1, :] *= 1E-9    # [GW]
print('Average S 14 & 22 [GW]', np.average(
    pwr_14[1, :]), np.average(pwr_22[1, :]),)

ax2.plot(pwr_14[0, :], pwr_14[1, :],
         linewidth=1.0, color=cud4[0], zorder=0.8)
ax2.plot(pwr_22[0, :], pwr_22[1, :],
         linewidth=1.0, color=cud4[3], zorder=0.8)
ax.yaxis.set_major_formatter(
    ptick.ScalarFormatter(useMathText=True))    # 指数表記
pwr_14_arr = np.zeros((12, pwr_14[1, :].size))
pwr_22_arr = np.zeros((12, pwr_22[1, :].size))
for i in range(12):
    pwr_14_1 = np.loadtxt('data/Poyntingflux/PY_2014_R4_edge_'+str(i)+'.txt')
    pwr_22_1 = np.loadtxt('data/Poyntingflux/PY_2022_R4_edge_'+str(i)+'.txt')
    pwr_14_arr[i, :] = pwr_14_1[1, :]*1E-9   # [GW]
    pwr_22_arr[i, :] = pwr_22_1[1, :]*1E-9   # [GW]
ax2.fill_between(x=pwr_14[0, :],
                 y1=pwr_14_arr.max(axis=0),
                 y2=pwr_14_arr.min(axis=0),
                 alpha=0.6,
                 color=cud4bs[0],
                 edgecolor=None,
                 zorder=0.6,)
ax2.fill_between(x=pwr_22[0, :],
                 y1=pwr_22_arr.max(axis=0),
                 y2=pwr_22_arr.min(axis=0),
                 alpha=0.6,
                 color=cud4bs[3],
                 edgecolor=None,
                 zorder=0.6,)

fig.tight_layout()

plt.savefig('Total_power_CR5.png')
plt.show()
