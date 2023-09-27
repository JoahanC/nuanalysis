from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np


tot_data = fits.getdata("./../bifrost_data/60502017008/event_cl/nu60502017008_2.96-7.96_keV_stacked.evt", 1)
events = []
times = []
for datum in tot_data:
    #print(datum)
    #if float(datum[0]) > float(tstart) and float(datum[0]) < float(tstop):
    #    datum_xpix = datum[13]
    #    datum_ypix = datum[14]
    #    if np.sqrt((float(xpix) - float(datum_xpix))**2 + (float(ypix) - float(datum_ypix))**2) < self.rlimit / 2.5:
    #        events.append(datum)
    times.append(float(datum[0]))
times = np.array(times)
tstart = np.min(times)
tstop = np.max(times)
lc_bins = np.arange(float(tstart), float(tstop), 10)
lc, bines = np.histogram(times, lc_bins)
t = []
for idx, l in enumerate(lc):
    t.append(idx)
t = np.array(t) * 10

fig, axes = plt.subplots(ncols=1, nrows=3)
fig.patch.set_facecolor('black')
ax1 = axes[0]
ax2 = axes[1]
ax3 = axes[2]
ax1.plot(t, lc, ms= 1)
ax2.plot(t, lc, ms= 1)
ax3.plot(t, lc, ms= 1)
ax1.tick_params(color='white', labelcolor='white')
ax2.tick_params(color='white', labelcolor='white')
ax3.tick_params(color='white', labelcolor='white')
fig.text(0.5, 0.04, 'Elapsed Time', ha='center', color="white", fontsize=15)
fig.text(0.07, 0.5, 'Counts', va='center', rotation='vertical', color="white", fontsize=15)
plt.setp(ax1.get_xticklabels(), visible=False)
plt.setp(ax2.get_xticklabels(), visible=False)
ax1.tick_params(axis='both', which='major', labelsize=15)
ax1.tick_params(axis='both', which='minor', labelsize=15)
ax2.tick_params(axis='both', which='major', labelsize=15)
ax2.tick_params(axis='both', which='minor', labelsize=15)
ax3.tick_params(axis='both', which='major', labelsize=15)
ax3.tick_params(axis='both', which='minor', labelsize=15)
ax1.set_xlim(46000, 50000)
ax2.set_xlim(46000, 50000)
ax3.set_xlim(46000, 50000)

telapsed = tstop - tstart
indices = np.arange(0, telapsed, 500)
intervals = []
last = 0
for i in range(len(indices) - 1):
    last = indices[i + 1]
    intervals.append([indices[i], indices[i + 1]])
intervals.append([last, telapsed])
#print
#print(indices)
intervals2 = [[0, 500]]
indices2 = np.arange(250, telapsed, 500)
last = 0
for i in range(len(indices2) - 1):
    last = indices2[i + 1]
    intervals2.append([indices2[i], indices2[i + 1]])
intervals2.append([last, telapsed])
for idx, interval in enumerate(intervals):
    if idx % 2 == 0:
        ax1.axvspan(interval[0], interval[1], alpha = 0.5, color='red')
        ax3.axvline(interval[0], color='r')
        ax3.axvline(interval[1], color='r')
    else:
        ax1.axvspan(interval[0], interval[1], alpha = 0.3, color='red')
for idx, interval in enumerate(intervals2):
    if idx % 2 == 0:
        ax2.axvspan(interval[0], interval[1], alpha = 0.5, color='green')
        ax3.axvline(interval[0], color='g')
        ax3.axvline(interval[1], color='g')
    else:
        ax2.axvspan(interval[0], interval[1], alpha = 0.3, color='green')
ax1.set_title('dt = 500 seconds', color='white', loc='right', fontsize=15)
plt.show()
plt.close()

"""lc_bins = np.arange(float(tstart), float(tstop), 10)
lc, bines = np.histogram(times, lc_bins)
t = []
for idx, l in enumerate(lc):
    t.append(idx)
t = np.array(t) * 10
plt.plot(t, lc, ms= 1)
plt.show()
plt.close()"""