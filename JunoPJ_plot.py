"""
JunoPJ_plot.py

Equatorial lead angle obtained by Juno/UVS observations.
The data is available in Hue+2023.

"""


import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import matplotlib.patches as patches
import astropy.io.fits as fits
from matplotlib.colors import LinearSegmentedColormap  # colormapをカスタマイズする
import pandas as pd
import ftpS3

# Color universal design
cud4 = ['#FF3300', '#FFF100', '#03AF7A', '#005AFF',
        '#4DC4FF', '#FF8082', '#F6AA00', '#990099', '#804000']
cud4bs = ['#FFCABF', '#FFFF80', '#D8F255', '#BFE4FF',
          '#FFCA80', '#77D9A8', '#C9ACE6', '#84919E']

# matplotlib フォント設定
fontname = 'Nimbus Sans'
plt.rcParams.update({'font.sans-serif': fontname,
                     'font.family': 'sans-serif',
                     'mathtext.fontset': 'custom',
                     'mathtext.rm': fontname,
                     'mathtext.it': fontname+':italic',
                     # 'mathtext.bf': 'Nimbus Sans:italic:bold',
                     'mathtext.bf': fontname+':bold'
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


hem = 'south'

refnum = 0
if hem == 'south':
    refnum = 1
satoval = np.recfromtxt('ref/2021je007055-sup-000'+str(3+refnum)+'-table si-s0'+str(2+refnum)+'.txt', skip_header=3,
                        names=['wlon', 'amlat', 'amwlon', 'iolat', 'iowlon', 'eulat', 'euwlon', 'galat', 'gawlon'])

pdcsv = pd.read_csv('data/Hue23_EFP_'+hem+'.txt', skiprows=2,
                    names=['PJ', 'UTC_TIME', 'Moon_SIII_LON', 'FP_LON', 'FP_LAT', 'FP_LON_ERR', 'FP_LAT_ERR', 'EMISSION_ANGLE'], sep='\t')
# print(pdcsv['PJ'])
# print(pdcsv['PJ'].str.contains('PJ07'))

PJ_list = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10',
           '11', '12', '13', '14', '15', '16', '17', '18', '19', '20',
           '21', '22', '23', '24', '25', '26', '27', '28', '29', '30',
           '31', '32', '33', '34', '35', '36', '37', '38', '39', '40',
           '41', '42', '43']

for j in range(len(PJ_list)):
    # Specify PJ number
    PJnum = 'PJ'+PJ_list[j]
    PJdata = pdcsv[pdcsv['PJ'].str.contains(PJnum)]
    # print(PJdata)

    h23_moonS3 = np.array(PJdata['Moon_SIII_LON'])
    h23_efpWlon = np.array(PJdata['FP_LON'])
    h23_efpWlon_ERR = np.array(PJdata['FP_LON_ERR'])

    h23_efpWlonEQ = np.zeros(h23_efpWlon.shape)
    h23_efpWlonEQ_P = np.zeros(h23_efpWlon.shape)
    h23_efpWlonEQ_M = np.zeros(h23_efpWlon.shape)
    for i in range(h23_efpWlonEQ.size):
        h23_efpWlonEQ[i] = ftpS3.ftpS3().S3EQ(
            h23_efpWlon[i], satoval, 'EUROPA')
        h23_efpWlonEQ_P[i] = ftpS3.ftpS3().S3EQ(
            h23_efpWlon[i]+h23_efpWlon_ERR[i], satoval, 'EUROPA')
        h23_efpWlonEQ_M[i] = ftpS3.ftpS3().S3EQ(
            h23_efpWlon[i]-h23_efpWlon_ERR[i], satoval, 'EUROPA')

    h23_leadA = h23_moonS3-h23_efpWlonEQ
    ERR_M = h23_moonS3 - h23_efpWlonEQ_M
    ERR_P = h23_moonS3 - h23_efpWlonEQ_P

    h23_leadA_ERR = np.zeros((2, h23_efpWlon.size))
    h23_leadA_ERR[0, :] = ERR_M - h23_leadA
    h23_leadA_ERR[1, :] = h23_leadA - ERR_P

    fontsize = 18
    fig, ax = plt.subplots(figsize=(5, 3.5), dpi=226)
    ax.tick_params(axis='x', labelsize=fontsize)
    ax.tick_params(axis='y', labelsize=fontsize)
    ax.set_title(PJnum+' ('+hem+')', weight='bold', fontsize=fontsize)
    ax.set_xlabel('Moon System III long. [deg]', fontsize=fontsize)
    ax.set_ylabel('Eq. lead angle [deg]', fontsize=fontsize)
    ax.set_xlim(0, 360)
    ax.set_ylim(-2.05, 14)
    ax.set_xticks(np.arange(0, 361, 45))
    ax.set_xticklabels(np.arange(0, 361, 45))
    ax.xaxis.set_minor_locator(AutoMinorLocator(3))  # minor ticks
    ax.set_yticks(np.arange(0, 14, 5))
    ax.set_yticklabels(np.arange(0, 14, 5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))  # minor ticks

    ax.scatter(h23_moonS3, h23_leadA, marker=',', s=0.95, c=cud4[0],
               linewidths=0.7, zorder=1)
    ax.errorbar(h23_moonS3, h23_leadA, yerr=h23_leadA_ERR,
                linestyle='none', ecolor=cud4[0], elinewidth=0.55, marker='none', zorder=1.5)

    fig.tight_layout()
    plt.savefig('img/Juno_PJ/'+hem+'/Juno_'+PJnum+'_'+hem+'_EUR.jpg')
    plt.close()
