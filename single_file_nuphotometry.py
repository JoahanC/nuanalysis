import astropy
import numpy as np
from wrappers import *
from astropy.coordinates import SkyCoord
from astropy.wcs.utils import skycoord_to_pixel
from astropy.io.fits import getdata
from photutils.aperture import RectangularAperture, CircularAperture, CircularAnnulus
from tqdm import tqdm
from nuphotometry import *


#make_det1_image("./../bifrost_data/60701033002/event_cl/nu60701033002A01_cl.evt", elow=1.6, ehigh=79, outpath="./../bifrost_data/60701033002/event_cl")
stacked_data = astropy.io.fits.getdata("./../bifrost_data/60701033002/event_cl/nu60701033002A01_cl_1.6to79keV_det1.fits")
stacked_a_data = astropy.io.fits.getdata("./../bifrost_data/60701033002/event_cl/nu60701033002A01_cl.evt")
wcs = astropy.wcs.WCS(astropy.io.fits.getheader("./../bifrost_data/60701033002/science.fits"))
header = astropy.io.fits.getheader("./test/nu60701033002A01_cl_1.6to79keV_det1.fits")
ra = header['RA_OBJ']
dec = header['DEC_OBJ']
obj_pos = SkyCoord(ra, dec, unit='deg')
pix_coordinates = [skycoord_to_pixel(obj_pos, wcs)]
print(pix_coordinates)



x_vals = []
y_vals = []
maping = getdata("./../bifrost_data/60701033002/event_cl/nu60701033002A_det1.fits")
for datum in stacked_a_data:
    x = datum[13]
    y = datum[14]
    if x > pix_coordinates[0][0] - 0.5 and x < pix_coordinates[0][0]+ 0.5:
        if y > pix_coordinates[0][1] - 0.5 and y < pix_coordinates[0][1] + 0.5:
            x_vals.append(datum[9])
            y_vals.append(datum[10])

det_x = np.mean(x_vals)
det_y = np.mean(y_vals)
src_coordinates = [[det_x, det_y]]
      
total_aper = RectangularAperture([180.0, 180.0], 330.0, 330.0)
total_sums, total_error = total_aper.do_photometry(stacked_data)
count_ratios = []
min_deviance = 1
optimal_radius = 0
for i in tqdm(np.arange(30, 290, 1)):
    source_aper = CircularAperture(src_coordinates, i)
    src_sums, src_error = source_aper.do_photometry(stacked_data)
    count_ratio = src_sums / total_sums
    difference = np.abs(0.90 - count_ratio)
    if difference < min_deviance:
        optimal_radius = i
        min_deviance = difference
    count_ratios.append((i, count_ratio))


background_area = calculate_background_area(stacked_data, det_x, det_y, optimal_radius)
arcsecond_area = background_area * u.arcsecond * u.arcsecond
print(f"Background area: {arcsecond_area}")

print(optimal_radius, min_deviance)
source_aper = CircularAperture(src_coordinates, optimal_radius)
bkg_ann = CircularAnnulus(src_coordinates, optimal_radius, optimal_radius + 30)
bkg_count, bkg_error = source_aper.do_photometry(stacked_data)
bkg_count = bkg_count * u.ct
area_count_rate = bkg_count / arcsecond_area
print(f"Area Count Rate: {area_count_rate}")


