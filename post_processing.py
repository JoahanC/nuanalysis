import numpy as np 
import matplotlib.pyplot as plt 
from astropy.utils.data import get_pkg_data_filename
from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
from astropy.coordinates import SkyCoord
from regions import PixCoord
from regions import CircleSkyRegion, CirclePixelRegion
import radial_profile

from helpers import *
from nuanalysis import NuAnalysis

seqid = "60965003002"
#seqid = "70901004002"
path = f"../data/{seqid}/"
evdir = f"{path}event_cl/"
out_path = f"{path}products/"
test = NuAnalysis(10, 3, path=path, evdir=evdir, seqid=seqid, out_path=out_path, clean=True)
test.test2_gen()
#test.process_detections()

'''for bound in test._phi_bounds[2:3]:
    data = test.read_unique_detections(bound)
    outfile = "test2.fits"
    hdu = fits.open(outfile, uint=True)[0]
    wcs = WCS(hdu.header)
    infile = "./../data/70901004002/event_cl/nu70901004002A01_cl.evt , ./../data/70901004002/event_cl/nu70901004002B01_cl.evt"
    outfile = "./../data/70901004002/qa/test2.fits"'''
    
"""outfile = "test2.fits"
    hdu = fits.open(outfile, uint=True)[0]
    wcs = WCS(hdu.header)
    coords = (hdu.header["RA_OBJ"], hdu.header["DEC_OBJ"])
    coordinates = radial_profile.find_source(outfile, show_image=False, filt_range=3)
    rind, rad_profile, radial_err, psf_profile = radial_profile.make_radial_profile(outfile, show_image=False,
                                                                 coordinates = coordinates)
    rlimit = radial_profile.optimize_radius_snr(rind, rad_profile, radial_err, psf_profile, show=False)
    
    plt.subplot(projection=wcs)
    ax = plt.gca()
    ax.imshow(hdu.data, origin='lower')
    center_sky = SkyCoord(coords[0], coords[1], unit='deg', frame='fk5')
    region_sky = CircleSkyRegion(center=center_sky, radius=50*u.arcsecond)
    pixel_region = region_sky.to_pixel(wcs)
    pixel_region.plot(ax=ax, color="green")
    
    for idx in data["INDEX"]:
        tmp_coords = test.ra_dec_todeg(data["RA"][int(idx) - 1], data["DEC"][int(idx) - 1])
        center_sky = SkyCoord(tmp_coords[0], tmp_coords[1], unit='deg', frame='fk5')
        region_sky = CircleSkyRegion(center=center_sky, radius=5*u.arcsecond)
        pixel_region = region_sky.to_pixel(wcs)
        pixel_region.plot(ax=ax, color="white")
    plt.xlabel('')
    plt.ylabel('')
    plt.show()"""



#infile = "./../data/30501002002/event_cl/nu30501002002A01_cl.evt , ./../data/30501002002/event_cl/nu30501002002B01_cl.evt"
#outfile = "./../data/30501002002/qa/test2.fits"
#coordinates = radial_profile.find_source(outfile, show_image=False, filt_range=3)
#rind, rad_profile, radial_err, psf_profile = radial_profile.make_radial_profile(outfile, show_image=False,
#                                                                 coordinates = coordinates)
#rlimit = radial_profile.optimize_radius_snr(rind, rad_profile, radial_err, psf_profile, show=False)
#print(rlimit)

"""hdu = fits.open(outfile, uint=True)[0]
wcs = WCS(hdu.header)
coords = (hdu.header["RA_OBJ"], hdu.header["DEC_OBJ"])

print(wcs.proj_plane_pixel_scales())
center_sky = SkyCoord(176.869, -61.953728, unit='deg', frame='fk5')
region_sky = CircleSkyRegion(center=center_sky, radius=rlimit*u.arcsecond)
pixel_region = region_sky.to_pixel(wcs)



plt.subplot(projection=wcs)
ax = plt.gca()
ax.imshow(hdu.data, origin='lower')
pixel_region.plot(ax=ax, color="white")
plt.grid(color='white', ls='solid')
plt.xlabel('Galactic Longitude')
plt.ylabel('Galactic Latitude')
plt.show()"""
