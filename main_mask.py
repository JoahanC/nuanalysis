"""
A file for testing out methods for masking out main sources and identifying them.
"""
from astropy.stats import sigma_clipped_stats
from photutils.detection import find_peaks
from astropy.io.fits import getdata 
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from astropy.visualization import simple_norm
from astropy.visualization.mpl_normalize import ImageNormalize
from photutils.aperture import CircularAperture

# Generic case
data = getdata("../bifrost_data/80902312004/science.fits")
mean, median, std = sigma_clipped_stats(data, sigma=3.0)
threshold = median + (5.0 * std)
tbl = find_peaks(data, threshold, box_size=11)
tbl['peak_value'].info.format = '%.8g'  # for consistent table output
print(tbl[:10])  # print only the first 10 peaks

fig, ax = plt.subplots()

ax.imshow(data, origin='lower', norm=matplotlib.colors.LogNorm())
ax.set_xlim(250, 750)
ax.set_ylim(250, 750)

plt.show()
plt.close()

positions = np.transpose((tbl['x_peak'], tbl['y_peak']))
apertures = CircularAperture(positions, r=5.0)
norm = simple_norm(data, 'log', percent=99.9)
plt.imshow(data, cmap='Greys_r', origin='lower', norm=norm,
           interpolation='nearest')
apertures.plot(color='#0547f9', lw=1.5)
plt.xlim(0, data.shape[1] - 1)
plt.ylim(0, data.shape[0] - 1)
plt.show()

