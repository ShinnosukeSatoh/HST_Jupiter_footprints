""" F25SR2F.py

F25SR2Fバンドパスフィルターのthroughputを計算するプログラム

"""


import numpy as np
import math
import matplotlib.pyplot as plt

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


# List of wavelength [A]
wl = np.arange(1300, 1825+1, 25)

# List of throughput [%]
th = np.array([
    2.10, 2.59, 2.54, 2.42,
    2.27, 2.07, 1.85, 1.63,
    1.42, 1.25, 1.10, 0.94,
    0.81, 0.71, 0.62, 0.52,
    0.42, 0.35, 0.29, 0.24,
    0.20, 0.15
])

# Data load
a = np.loadtxt('data/F25SRF2/HST_STIS_FUV.F25SRF2.dat')
# a[:, 0] wavelength
# a[:, 1] throughput [ratio]


# Plot
fontsize = 19
fig, ax = plt.subplots()
ax.tick_params(axis='both', labelsize=fontsize)
ax.set_xlabel('Wavelength [nm]', fontsize=fontsize)
ax.set_ylabel('Throughput [%]', fontsize=fontsize)
ax.set_yscale('log')
ax.set_xlim(100, 340)
ax.set_ylim(1E-6, 10)
ax.plot(wl/10, th, color='k')
ax.plot(a[:, 0]/10, a[:, 1]*100, color='r')

fig.tight_layout()
plt.show()


# Calculate throughput
wl_input = 142.1          # [nm]
wl_data = a[:, 0]/10    # [nm]
th_data = a[:, 1]*100   # [%]
wl_del = np.abs(wl_data-wl_input)
argsorted = np.argsort(wl_del)
idx0, idx1 = argsorted[0], argsorted[1]
print(wl_data[idx0], wl_data[idx1])
print(th_data[idx0], th_data[idx1])

# Linearly interpolated
dx = wl_data[idx1]-wl_data[idx0]
dy = th_data[idx1]-th_data[idx0]
y = (dy/dx)*(wl_input-wl_data[idx0]) + th_data[idx0]

print('Linear')
print(y)


# Average throughput
S = 0
for i in range(wl_data.size-1):
    S += (th_data[i+1]+th_data[i])*(wl_data[i+1]-wl_data[i])/2
th_ave = S/(wl_data[-2]-wl_data[0])
print(S, th_ave, (wl_data[-2]-wl_data[0]))

# Plot
fontsize = 19
fig, ax = plt.subplots()
ax.set_title('F25SRF2 Throughput', fontsize=fontsize, weight='bold')
ax.tick_params(axis='both', labelsize=fontsize)
ax.set_xlabel('Wavelength [nm]', fontsize=fontsize)
ax.set_ylabel('Throughput [%]', fontsize=fontsize)
# ax.set_yscale('log')
ax.set_xlim(100, 340)
ax.set_ylim(1E-6, 3)
ax.plot(wl/10, th, color='k')
ax.plot(a[:, 0]/10, a[:, 1]*100, color='r')
ax.axhline(y=th_ave, color='b')
ax.axvline(x=wl_data[0], color='gray')
ax.axvline(x=wl_data[-2], color='gray')

fig.tight_layout()
plt.show()
