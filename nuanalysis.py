"""
Main class for performing X-ray FOV source detections on NuSTAR observations. 
"""

import os
import subprocess
import glob
import matplotlib
import radial_profile
import numpy as np
import matplotlib.pyplot as plt

from shapely.geometry import Point
from shapely.plotting import plot_polygon
from astropy.io import fits
from scipy.stats import poisson_means_test
from tqdm import tqdm

from nustar import *
from helpers import *
from wrappers import *
from photometry import *

from astropy import units as u
from astropy.wcs import WCS
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.wcs.utils import skycoord_to_pixel, pixel_to_skycoord
from regions import CircleSkyRegion
from photutils.aperture import RectangularAperture, CircularAperture, CircularAnnulus



class NuAnalysis(Observation):
    """
    This class defines an object to be used for performing analysis on a NuSTAR observation.
    
    Arguments
    ---------
    
    dtime : float
        The characteristic timescale to perform analysis and binning
    snr_threshold : float
        The snr limit that is used in ximage calls during slidecell detection routines
    low_phi_file : str
        The file defining the lower bounds for PI channel bands
    high_phi_file : str
        The file defining the higher bounds for PI channel bands
    path : str
        The high level directory pointing to the seqid directory 
    seqid : str
        The sequence id used to label this observation within NuSTAR's catalog
    bifrost : bool
        A flag indicating if analysis is being performed on the bifrost server
    object_name : str
        The object_name of the observation
    """

    def __init__(self, dtime, snr_threshold, low_phi_file, high_phi_file, 
                 path=False, seqid=False, clean=False, bifrost=False, object_name=None):
        
        # Initialize NuSTAR object
        self.ns = NuSTAR()

        # Define class parameters and core file locations
        self._object = object_name
        self._seqid = seqid
        self._mainpath = os.path.abspath(path)
        self._evtpath = os.path.join(self._mainpath, "event_cl")
        self._outpath = os.path.join(self._mainpath, "outputs")
        self._impath = os.path.join(self._mainpath, "images")
        self._detpath = os.path.join(self._mainpath, "detections")
        self._logpath = os.path.join(self._mainpath, "logs")
        self._low_phi_file = os.path.abspath(low_phi_file)
        self._high_phi_file = os.path.abspath(high_phi_file)
        self._maincontents = os.listdir(self._mainpath)

        # Ensure core folders are present
        if self._evtpath not in self._maincontents:
            generate_directory(self._evtpath)
        if self._outpath not in self._maincontents:
            generate_directory(self._outpath)
        if self._impath not in self._maincontents:
            generate_directory(self._impath)
        if self._detpath not in self._maincontents:
            generate_directory(self._detpath)
        if self._logpath not in self._maincontents:
            generate_directory(self._logpath)
        
        # Record all present files
        self._maincontents = os.listdir(self._mainpath)
        self._evtcontents = os.listdir(self._evtpath)
        self._outcontents = os.listdir(self._outpath)
        self._imcontents = os.listdir(self._impath)
        self._detectcontents = os.listdir(self._detpath)
        self._sessionid = generate_random_id()
        os.chdir(self._mainpath)

        # Read in analysis parameters
        self._snr = snr_threshold
        self._dtime = dtime
        self._clean = clean
        self._phi_channels, self._kev_levels = self._read_in_phi_channels(self._low_phi_file, 
                                                                         self._high_phi_file)
        
        # Check for reprocessing
        if not clean:

            print("Running nuproducts on observation.")
            # Special case for bifrost
            if bifrost:
                rawfilepath = f"./../../../../../nustar/fltops/{self._object}/{self._seqid}"
                os.system(f"nupipeline {rawfilepath} nu{self._seqid} ./event_cl saamode=STRICT clobber=yes")
                with open(self._evtpath + "processing_flag.txt", "w") as file:
                    file.write("PROCESSING COMPLETE")

                # Update known files
                self._maincontents = os.listdir(self._mainpath)
                self._evtcontents = os.listdir(self._evtpath)
                self._clean = True
            
            else:
                os.system(f"nupipeline ./ nu{self._seqid} ./event_cl saamode=STRICT clobber=yes")

                # Update known files
                self._maincontents = os.listdir(self._mainpath)
                self._evtcontents = os.listdir(self._evtpath)
                self._clean = True

        # Perform file sanity check
        self.completeness_flag = True
        if f"nu{self._seqid}A01_cl.evt" not in self._evtcontents:
            print("Clean events file for FPMA missing!")
            self.completeness_flag = False
        if f"nu{self._seqid}B01_cl.evt" not in self._evtcontents:
            print("Clean events file for FPMB missing!")
            self.completeness_flag = False
        
        if self.completeness_flag:
        
            self._fpma_eventpath = os.path.join(self._evtpath, f"nu{self._seqid}A01_cl.evt")
            self._fpmb_eventpath = os.path.join(self._evtpath, f"nu{self._seqid}B01_cl.evt")

            # Initialize the Observation superclass
            super().__init__(self._mainpath, self._seqid, self._evtpath, self._outpath)
            

            # Generate subsets of original event files
            self._process_evt_files()
            
            # Generate optimal radius for snr 
            self._im_skycoord = pixel_to_skycoord(self._source_pix_coordinates[0][0], self._source_pix_coordinates[0][1], self.wcs)
            rind, rad_profile, radial_err, psf_profile = radial_profile.make_radial_profile(self._im_paths[f"{self._kev_levels[0][0]}-{self._kev_levels[0][1]}"][2],
                                                                                            show_image=False,
                                                                                            coordinates = self._source_pix_coordinates)
            self.rlimit = radial_profile.optimize_radius_snr(rind, rad_profile, radial_err, psf_profile, show=False)
            
            # Buffer radius to reduce main source polluting detections
            if self.rlimit == 0:
                self.rlimit = 120
            #if self.rlimit <= 90:
            #    self.rlimit = 90
            
            # Initialize detection parameters
            self._time_bins = self.generate_timebins()
            self._detections = None
            
            # Ensure detector images are generated
            if "det1_coords.txt" not in self._imcontents:
                print("Detector coordinates not found! Estimating coordinates.")
                self.source_det1_coords()
            #self.generate_background_images()
            
            # Read in detector coordinates
            self._source_det1_coords = {}
            
            with open(os.path.join(self._impath, "det1_coords.txt"), 'r') as file:
                data_line = file.readline().split()
                self._source_det1_coords['A'] = [float(data_line[0]), float(data_line[1])]
                self._source_det1_coords['B'] = [float(data_line[2]), float(data_line[3])]

            self._display_terminal()

    # Mutable properties begin below

    @property
    def snr(self):
        """
        Returns the set SNR threshold for astrometry analysis.
        """
        return self._snr


    @property
    def dtime(self):
        """
        Returns the time cutting parametry used for astrometry analysis.
        """
        return self._dtime


    @property
    def phi_channels(self):
        """
        Returns the phi_channels set for analysis.
        """
        return self._phi_channels


    @property 
    def time_bins(self):
        """
        Returns the time bins set for analysis.
        """
        return self._time_bins
    

    @property
    def detections(self):
        """
        Returns the detections found for this object.
        """
        return self._detections
    
    # Initialization methods below

    def _read_in_phi_channels(self, pilow_file, pihi_file):
        """
        Makes use of two files `nustar_pilow.txt` and `nustar_pihi.txt` within the `nustar` directory to set PI channel thresholds 
        for data analysis using NUSTARDAS

        Arguments
        ---------
        pilow_file : `str`
            The name of the file containing the pilow bounds.
        pihi_file : `str`
            The name of the file containing the pihi bounds.
        
        Raises
        ------
        IndexError
            when the PI channel bounds are not of the same length
        """

        phi_channels = []
        kev_levels = []
        low_phi = open(pilow_file, 'r')
        high_phi = open(pihi_file, 'r')
        low_phis = low_phi.readlines()
        high_phis = high_phi.readlines()
        
        if len(low_phis) != len(high_phis):
            raise IndexError("PHI bounds must have the same length!")
        
        for idx in range(len(low_phis)):
            pi_low = int(low_phis[idx].replace('\n', ''))
            pi_high = int(high_phis[idx].replace('\n', ''))
            kev_low = round(self.ns.channel_to_energy(float(pi_low)), 3)
            kev_high = round(self.ns.channel_to_energy(float(pi_high)), 3)
            phi_channels.append((pi_low, pi_high))
            kev_levels.append((kev_low, kev_high))
        
        return phi_channels, kev_levels


    def generate_timebins(self):
        """
        This method executes a simple splitting algorithm to dither up the NuSTAR observation by a set `dt` value 
        defined at initialization. 

        Arguments
        ---------
        None
        
        Returns
        -------
        list
            Returns a list of two dictionaries which correspond to the dithered time segments. The first dictionary begins 
            at `TSTART` while the second dictionary begins at `TSTART` + dt / 2. The structure of each dictionary is 

            key : [interval_start, interval_end]
            value : list of all events inclusively contained in each interval
        
        Raises
        ------
        None
        """

        # generate intervals for PASS 1
        intervals_p1 = []
        cur_int = self._event_times[0]
        max_int = np.max(self._event_times)
        while cur_int < max_int:
            intervals_p1.append(cur_int)
            cur_int += self._dtime
        intervals_p1.append(max_int)
        self._start_intervals_pass1 = len(intervals_p1)

        # generate intervals for PASS 2
        intervals_p2 = []
        intervals_p2.append(self._event_times[0])
        cur_int = self._event_times[0] + self._dtime / 2
        max_int = np.max(self._event_times)
        while cur_int < max_int:
            intervals_p2.append(cur_int)
            cur_int += self._dtime 
        intervals_p2.append(max_int)
        self._start_intervals_pass2 = len(intervals_p2)

        # split data for PASS 1
        idx = 0
        first = True
        data_intervals_1 = {}
        last = self._event_times[-1]
        for datapoint in self._event_times:
            if first:
                data = []
                data_intervals_1[idx] = [intervals_p1[idx], intervals_p1[idx + 1], data]
                data_intervals_1[idx][2].append(datapoint)
                first = False
                continue
            if datapoint == last:
                if idx not in data_intervals_1:
                    data = []
                    data_intervals_1[idx] = [intervals_p1[idx], intervals_p1[idx + 1], data]
                    data_intervals_1[idx][2].append(datapoint)
                    break
                else:
                    data_intervals_1[idx][2].append(datapoint)
                    break
            if datapoint > intervals_p1[idx] and datapoint <= intervals_p1[idx + 1]:
                data_intervals_1[idx][2].append(datapoint)
            else:
                if datapoint > intervals_p1[idx + 1] and datapoint < intervals_p1[idx + 2]:
                    idx += 1
                    data = []
                    data_intervals_1[idx] = [intervals_p1[idx], intervals_p1[idx + 1], data]
                    data_intervals_1[idx][2].append(datapoint)
                else:
                    idx += 1

        # split data for PASS 2
        idx = 0
        first = True
        data_intervals_2 = {}
        last = self._event_times[-1]
        for datapoint in self._event_times:
            if first:
                data = []
                data_intervals_2[idx] = [intervals_p2[idx], intervals_p2[idx + 1], data]
                data_intervals_2[idx][2].append(datapoint)
                first = False
                continue
            if datapoint == last:
                if idx not in data_intervals_2:
                    data = []
                    data_intervals_2[idx] = [intervals_p2[idx], intervals_p2[idx + 1], data]
                    data_intervals_2[idx][2].append(datapoint)
                    break
                else:
                    data_intervals_2[idx][2].append(datapoint)
                    break
            if datapoint > intervals_p2[idx] and datapoint <= intervals_p2[idx + 1]:
                data_intervals_2[idx][2].append(datapoint)
            else:
                if datapoint > intervals_p2[idx + 1] and datapoint < intervals_p2[idx + 2]:
                    idx += 1
                    data = []
                    data_intervals_2[idx] = [intervals_p2[idx], intervals_p2[idx + 1], data]
                    data_intervals_2[idx][2].append(datapoint)
                else:
                    idx += 1
        
        return [data_intervals_1, data_intervals_2]


    def _process_evt_files(self):
        """
        Performs energy binning on clean event files and generates stacked fits image 
        files for photometry.
        """
        from astropy.io.fits import getdata

        # Filtered event files
        self._evt_files = {}

        # keV level focused FPMA evt files
        for level in self._kev_levels:
            infiles = os.path.relpath(self._fpma_eventpath)
            evpath = os.path.relpath(self._evtpath)
            outfile = os.path.join(evpath, f"nu{self._seqid}_{level[0]}-{level[1]}_keV_A01.evt")
            if f"nu{self._seqid}_{level[0]}-{level[1]}_keV_A01.evt" not in self._evtcontents:
                print(f"Generating nu{self._seqid}_{level[0]}-{level[1]}_keV_A01.evt")
                script_id = generate_random_id()
                evt_xselect_filter(infiles, outfile, '.', level[0], level[1], script_id)
                os.system(f"xselect @{script_id}xselect.xco > {os.path.join(os.path.relpath(self._logpath), f'xselect_{level[0]}-{level[1]}_keV_A01_evt')}.log")
                os.system(f"rm {script_id}xselect.xco")
            self._evt_files[f"{level[0]}-{level[1]}"] = [os.path.abspath(outfile)]

        # keV level focused FPMB evt files
        for level in self._kev_levels:
            infiles = os.path.relpath(self._fpmb_eventpath)
            evpath = os.path.relpath(self._evtpath)
            outfile = os.path.join(evpath, f"nu{self._seqid}_{level[0]}-{level[1]}_keV_B01.evt")
            if f"nu{self._seqid}_{level[0]}-{level[1]}_keV_B01.evt" not in self._evtcontents:
                print(f"Generating nu{self._seqid}_{level[0]}-{level[1]}_keV_B01.evt")
                script_id = generate_random_id()
                evt_xselect_filter(infiles, outfile, '.', level[0], level[1], script_id)
                os.system(f"xselect @{script_id}xselect.xco > {os.path.join(os.path.relpath(self._logpath), f'xselect_{level[0]}-{level[1]}_keV_B01_evt')}.log")
                os.system(f"rm {script_id}xselect.xco")
            self._evt_files[f"{level[0]}-{level[1]}"].append(os.path.abspath(outfile))

        # keV level focused stacked evt files
        for level in self._kev_levels:
            infiles = f"{os.path.relpath(self._fpma_eventpath)} , {os.path.relpath(self._fpmb_eventpath)}"
            evpath = os.path.relpath(self._evtpath)
            outfile = os.path.join(evpath, f"nu{self._seqid}_{level[0]}-{level[1]}_keV_stacked.evt")
            if f"nu{self._seqid}_{level[0]}-{level[1]}_keV_stacked.evt" not in self._evtcontents:
                print(f"Generating nu{self._seqid}_{level[0]}-{level[1]}_keV_stacked.evt")
                script_id = generate_random_id()
                evt_xselect_filter(infiles, outfile, '.', level[0], level[1], script_id)
                os.system(f"xselect @{script_id}xselect.xco > {os.path.join(os.path.relpath(self._logpath), f'xselect_{level[0]}-{level[1]}_keV_stacked_evt')}.log")
                os.system(f"rm {script_id}xselect.xco")
            self._evt_files[f"{level[0]}-{level[1]}"].append(os.path.abspath(outfile))

        # Filtered image files
        self._im_paths = {}

        # keV level focused stacked img files
        for level in self._kev_levels:
            infiles = os.path.relpath(self._fpma_eventpath)
            impath = os.path.relpath(self._impath)
            outfile = os.path.join(impath, f"nu{self._seqid}_{level[0]}-{level[1]}_keV_A01.fits")
            if f"nu{self._seqid}_{level[0]}-{level[1]}_keV_A01.fits" not in self._imcontents:
                print(f"Generating nu{self._seqid}_{level[0]}-{level[1]}_keV_A01.fits")
                script_id = generate_random_id()
                evt_to_fits_image(infiles, outfile, '.', level[0], level[1], script_id)
                os.system(f"xselect @{script_id}xselect.xco > {os.path.join(os.path.relpath(self._logpath), f'xselect_{level[0]}-{level[1]}_keV_A01_fits')}.log")
                os.system(f"rm {script_id}xselect.xco")
            self._im_paths[f"{level[0]}-{level[1]}"] = [os.path.abspath(outfile)]

        # keV level focused stacked img files
        for level in self._kev_levels:
            infiles = os.path.relpath(self._fpmb_eventpath)
            impath = os.path.relpath(self._impath)
            outfile = os.path.join(impath, f"nu{self._seqid}_{level[0]}-{level[1]}_keV_B01.fits")
            if f"nu{self._seqid}_{level[0]}-{level[1]}_keV_B01.fits" not in self._imcontents:
                print(f"Generating nu{self._seqid}_{level[0]}-{level[1]}_keV_B01.fits")
                script_id = generate_random_id()
                evt_to_fits_image(infiles, outfile, '.', level[0], level[1], script_id)
                os.system(f"xselect @{script_id}xselect.xco > {os.path.join(os.path.relpath(self._logpath), f'xselect_{level[0]}-{level[1]}_keV_B01_fits')}.log")
                os.system(f"rm {script_id}xselect.xco")
            self._im_paths[f"{level[0]}-{level[1]}"].append(os.path.abspath(outfile))

        # keV level focused stacked img files
        for level in self._kev_levels:
            infiles = f"{os.path.relpath(self._fpma_eventpath)} , {os.path.relpath(self._fpmb_eventpath)}"
            impath = os.path.relpath(self._impath)
            outfile = os.path.join(impath, f"nu{self._seqid}_{level[0]}-{level[1]}_keV_stacked.fits")
            if f"nu{self._seqid}_{level[0]}-{level[1]}_keV_stacked.fits" not in self._imcontents:
                print(f"Generating nu{self._seqid}_{level[0]}-{level[1]}_keV_stacked.fits")
                script_id = generate_random_id()
                evt_to_fits_image(infiles, outfile, '.', level[0], level[1], script_id)
                os.system(f"xselect @{script_id}xselect.xco > {os.path.join(os.path.relpath(self._logpath), f'xselect_{level[0]}-{level[1]}_keV_stacked_fits')}.log")
                os.system(f"rm {script_id}xselect.xco")
            self._im_paths[f"{level[0]}-{level[1]}"].append(os.path.abspath(outfile))

        # Store relevant parameters
        print(outfile)
        outfile = f"images/nu{self.seqid}_2.96-78.96_keV_stacked.fits"
        hdu = fits.open(outfile, uint=True)[0]
        self.stacked_data = getdata(outfile)
        self.wcs = WCS(hdu.header)
        self._source_pix_coordinates = [skycoord_to_pixel(self._source_position, self.wcs)]
        self.n_cuts = self.exposure['A01'] / self._dtime
        return
    
    
    def generate_detector_images(self):
        """
        Generates fits images which are in DET1 coordinates for specified energy ranges.
        """
        filtered_image_files = {}
        filtered_image_files['A'] = []
        filtered_image_files['B'] = []
        for channels in self._phi_channels:
            
            elow = round(self.ns.channel_to_energy(float(channels[0])), 3)
            ehigh = round(self.ns.channel_to_energy(float(channels[1])), 3)
            photometry_files = os.listdir(self._impath)
            
            if f"nu{self._seqid}A01_cl_{elow}to{ehigh}keV_det1.fits" not in photometry_files:
                a_file = make_det1_image(os.path.join(self._evtpath, f"nu{self._seqid}A01_cl.evt"), 
                                elow=elow,
                                ehigh=ehigh,
                                outpath=self._impath)
                filtered_image_files['A'].append(a_file)
            
            if f"nu{self._seqid}B01_cl_{elow}to{ehigh}keV_det1.fits" not in photometry_files:
                b_file = make_det1_image(os.path.join(self._evtpath, f"nu{self._seqid}B01_cl.evt"), 
                                elow=elow,
                                ehigh=ehigh,
                                outpath=self._impath)
                filtered_image_files['B'].append(b_file)
    
    
    def recalculate_center(self):
        import astropy
        from astropy.stats import sigma_clipped_stats
        from photutils.detection import find_peaks
        data = astropy.io.fits.getdata(self._im_paths["2.96-78.96"][2])
        mean, median, std = sigma_clipped_stats(data, sigma=3)
        if std == 0:
            std = 0.1
        tbl = find_peaks(data, std*10, box_size=20)
        tbl['peak_value'].info.format = '%.8g'  # for consistent table output
        tbl.sort('peak_value')
        tbl.reverse()
        print(tbl[0])
        #print(self._source_pix_coordinates)
        self._source_pix_coordinates = [[np.array(tbl['x_peak'][0]), np.array(tbl['y_peak'][0])]]
        #print(self._source_pix_coordinates)
    
    def event_extraction(self):
        """
        Performs the event extraction procedures for binning the event file for this observation into 
        time scales of self._dtime
        """
        
        # Display terminal
        print('#' * 90)
        print(f"Performing Event Extraction for SEQID: {self._seqid}")
        print(f"Time scale: {self._dtime} seconds.")
        print('#' * 90)

        if f"nu{self._seqid}A01_cl.evt" not in self._evtcontents:
            print("Core event files missing! Skipping!")
            return 0
        if f"nu{self._seqid}B01_cl.evt" not in self._evtcontents:
            print("Core event files missing! Skipping!")
            return 0
        
        os.chdir(self._evtpath)

        # Iterate through all PHI channel bounds
        for channel in self._phi_channels:

            # Generate detection folder for this PI channel
            channel_path = os.path.relpath(os.path.join(self._detpath, f"{channel[0]}-{channel[1]}_{self._dtime}-{self._snr}"))
            generate_directory(channel_path, overwrite=False)

            # Begin first pass through of data
            print(f"Binning data: PHI Channels {channel[0]}-{channel[1]}; Cycle 1")
            for interval in tqdm(self._time_bins[0]):
                if len(self._time_bins[0][interval][2]) == 0:
                    continue
                else:
                    outfile = os.path.join(channel_path, f"nu_{self._time_bins[0][interval][0]}-{self._time_bins[0][interval][1]}.evt")
                    logfile = os.path.join(self._logpath, f"binning_{channel[0]}-{channel[1]}_{self._time_bins[0][interval][0]}-{self._time_bins[0][interval][1]}.log")
                    with open("xselect.xco", 'w') as script:
                        script.write(f'{self._sessionid}\n')
                        script.write(f"read events\n")
                        script.write(".\n")
                        script.write(f"nu{self._seqid}A01_cl.evt , nu{self._seqid}B01_cl.evt\n")
                        script.write('yes\n')
                        script.write("filter PHA_CUTOFF\n")
                        script.write(f"{channel[0]}\n")
                        script.write(f"{channel[1]}\n")
                        script.write("filter time scc\n")
                        script.write(f"{self._time_bins[0][interval][0]} , {self._time_bins[0][interval][1]}\n")
                        script.write("x\n")
                        script.write("extract events\n")
                        script.write("\n")
                        script.write(f"save events {outfile}\n")
                        script.write('no\n')
                        script.write("exit no\n")
                    os.system(f"xselect @xselect.xco > {logfile}")
                    os.system("rm xselect.xco")
            
            # Begin second pass through of data
            print(f"Binning data: PHI Channels {channel[0]}-{channel[1]}; Cycle 2")
            for interval in tqdm(self._time_bins[1]):
                if len(self._time_bins[1][interval][2]) == 0:
                    continue
                else:
                    outfile = os.path.join(channel_path, f"nu_{self._time_bins[1][interval][0]}-{self._time_bins[1][interval][1]}.evt")
                    logfile = os.path.join(self._logpath, f"binning_{channel[0]}-{channel[1]}_{self._time_bins[1][interval][0]}-{self._time_bins[1][interval][1]}.log")
                    with open("xselect.xco", 'w') as script:
                        script.write(f'{self._sessionid}\n')
                        script.write(f"read events\n")
                        script.write(".\n")
                        script.write(f"nu{self._seqid}A01_cl.evt , nu{self._seqid}B01_cl.evt\n")
                        script.write('yes\n')
                        script.write("filter PHA_CUTOFF\n")
                        script.write(f"{channel[0]}\n")
                        script.write(f"{channel[1]}\n")
                        script.write("filter time scc\n")
                        script.write(f"{self._time_bins[1][interval][0]} , {self._time_bins[1][interval][1]}\n")
                        script.write("x\n")
                        script.write("extract events\n")
                        script.write("\n")
                        script.write(f"save events {outfile}\n")
                        script.write('no\n')
                        script.write("exit no\n")
                    os.system(f"xselect @xselect.xco > {logfile}")
                    os.system("rm xselect.xco")
        
        # Write completion flag
        with open(f"{self._dtime}_binning_flag.txt", "w") as file:
            file.write("PROCESSING COMPLETE")
        os.chdir(self._mainpath)
    
    
    def verify_event_extraction(self):
        """
        This method performs a check to verify if all of the relevant .evt files 
        have been created within the event extraction phase
        """

        # Display terminal
        print('#' * 90)
        print(f"Verifying Event Extraction for SEQID: {self._seqid}")
        print(f"Time scale: {self._dtime} seconds.")
        print('#' * 90)

        os.chdir(self._detpath)
        complete_vals = []
        
        # Iterate through all PHI channel bounds
        for channel in self._phi_channels:

            channel_path = os.path.relpath(os.path.join(self._detpath, f"{channel[0]}-{channel[1]}_{self._dtime}-{self._snr}"))
            if not os.path.isdir(channel_path):
                complete_vals.append(0)
            if os.path.isdir(channel_path):
                evt_files = glob.glob(os.path.join(channel_path, "*.evt"))
                trim_files = []
                for file in evt_files:
                    trim_files.append(file.replace(channel_path, '').replace('/', ''))

                # PASS 1 verification
                count = 0
                for idx in self._time_bins[0]:
                    test_file = f"nu_{self._time_bins[0][idx][0]}-{self._time_bins[0][idx][1]}.evt"
                    if test_file in trim_files:
                        count += 1
                    
                # PASS 2 verification
                for idx in self._time_bins[1]:
                    test_file = f"nu_{self._time_bins[1][idx][0]}-{self._time_bins[1][idx][1]}.evt"
                    if test_file in trim_files:
                        count += 1

                percent = round(count / (len(self._time_bins[0]) + len(self._time_bins[1])), 2)
                print(f"Percent of evt files properly processed: {percent}")
                complete_vals.append(round(percent, 4))

        return complete_vals
    
    
    def source_det1_coords(self):
        """
        This method is DEPRECATED
        This method estimates the position of the main source within the detector coordinate 
        frame.
        """
        from astropy.io.fits import getdata
        
        # FPMA
        fpma_event_data = getdata(self._fpma_eventpath)
        x_vals = []
        y_vals = []

        # Loop through all events
        for datum in fpma_event_data:
            x = datum[13]
            y = datum[14]

            # Record the detector coordinate if event is near the pixel coordinates of main source
            if x >= self._source_pix_coordinates[0][0] - 1 and x <= self._source_pix_coordinates[0][0] + 1:
                if y >= self._source_pix_coordinates[0][1] - 1 and y <= self._source_pix_coordinates[0][1] + 1:
                    x_vals.append(datum[9])
                    y_vals.append(datum[10])
        
        # Take a mean to find the pixel to within a 2x2 grid of the actual position
        det_x = np.mean(x_vals)
        det_y = np.mean(y_vals)
        a_src_coordinates = [[det_x, det_y]]
        
        # FPMB
        fpmb_event_data = getdata(self._fpmb_eventpath)
        x_vals = []
        y_vals = []

        # Loop through all events
        for datum in fpmb_event_data:
            x = datum[13]
            y = datum[14]

            # Record the detector coordinate if event is near the pixel coordinates of main source
            if x >= self._source_pix_coordinates[0][0] - 1 and x <= self._source_pix_coordinates[0][0] + 1:
                if y >= self._source_pix_coordinates[0][1] - 1 and y <= self._source_pix_coordinates[0][1] + 1:
                    x_vals.append(datum[9])
                    y_vals.append(datum[10])
        
        # Take a mean to find the pixel to within a 2x2 grid of the actual position
        det_x = np.mean(x_vals)
        det_y = np.mean(y_vals)
        b_src_coordinates = [[det_x, det_y]]

        with open(os.path.join(self._impath, "det1_coords.txt"), 'w') as file:
            file.write(f"{a_src_coordinates[0][0]} {a_src_coordinates[0][1]} \
                       {b_src_coordinates[0][0]} {b_src_coordinates[0][1]}")
    
    
    def read_final_detections(self, detectiontype="basic"):
        """
        Reads in all of the detections associated with this object

        Arguments
        ---------
        detectiontype : str
            The type or stage of detections desired to be read in. Supported 
            types include: ximage, poisson, and classify
        
        Returns
        -------
        """

        if detectiontype == "basic":
            detpath = os.path.relpath(os.path.join(self._detpath, f"{self._dtime}_detections.tbl"))
            filename = f"{self._dtime}_detections.tbl"
        if detectiontype == "poisson":
            detpath = os.path.relpath(os.path.join(self._detpath, f"{self._dtime}_poisson.tbl"))
            filename = f"{self._dtime}_poisson.tbl"
        if detectiontype == "classify":
            detpath = os.path.relpath(os.path.join(self._detpath, f"{self._dtime}_classify.tbl"))
            filename = f"{self._dtime}_classify.tbl"
        if detectiontype == "basicdet1":
            detpath = os.path.relpath(os.path.join(self._detpath, f"{self._dtime}_basic_det1.tbl"))
            filename = f"{self._dtime}_basic_det1.tbl"
        if detectiontype == "straycutbasic":
            detpath = os.path.relpath(os.path.join(self._detpath, f"{self._dtime}_straycut_basicdet1.tbl"))
            filename = f"{self._dtime}_straycut_basicdet1.tbl"
        if detectiontype == "straycut":
            detpath = os.path.relpath(os.path.join(self._detpath, f"{self._dtime}_straycut.tbl"))
            filename = f"{self._dtime}_straycut.tbl"
        if not os.path.isfile(detpath):
            print("No final detections found!")
            return None, False, None

        detect_info = Table.read(detpath, format='ipac')
        if len(detect_info['INDEX']) == 0:
            return detect_info, False, filename
        else:
            return detect_info, True, filename
        

    # Display methods below

    def _display_terminal(self):
        """
        Displays basic information of a newly loaded observation for ease-of-reading in the terminal
        """
        print("#" * 90)
        print(f"Succesfully read in {self._object}: {self._seqid}")
        print(f"Main Source located: ({self._source_pix_coordinates[0][0]}, {self._source_pix_coordinates[0][1]})")
        print("#" * 90)


    def display_image(self, savefig=False, display=True, style="dark"):
        """
        Displays the stacked image for this observation with the 
        source labeled with a circular region.

        Arguments
        ---------
        savefig : bool, optional
            A flag determining whether to save the plot
        display : bool, optional
            A flag determining whether to display the plot in an interactive window
        style : str, optional
            A string describing the image display style. Options include ``dark`` 
            which sets a black facecolor, and ``bright`` which sets a white 
            facecolor
        """

        # Initialize figure and plot elements
        fig = plt.figure()
        ax = fig.add_subplot(projection=self.wcs)
        im = ax.imshow(self.stacked_data, origin='lower', norm=matplotlib.colors.LogNorm())
        object_region = CircleSkyRegion(center=self._source_position, radius=(self.rlimit - 50) * u.arcsecond)
        plot_region = object_region.to_pixel(self.wcs)
        plot_region.plot(ax=ax, color="yellow")
        plt.grid(True)
        plt.xlim(250, 750)
        plt.ylim(250, 750)

        # Apply style parameters
        if style == "dark":
            fig.patch.set_facecolor("black")
            ax.set_facecolor("gray")
            ax.tick_params(color="white", labelcolor="white")
            ax.tick_params(axis="both", which="major", labelsize=15)
            im.axes.tick_params(color="white", labelcolor="white")
            plt.xlabel("Right Ascension", color="white", fontsize=18)
            plt.ylabel("Declination", color="white", fontsize=18)
        
        if style == "bright":
            fig.patch.set_facecolor("white")
            ax.set_facecolor("white")
            ax.tick_params(color="black", labelcolor="black")
            ax.tick_params(axis="both", which="major", labelsize=15)
            im.axes.tick_params(color="black", labelcolor="black")
            plt.xlabel("Right Ascension", color="black", fontsize=18)
            plt.ylabel("Declination", color="black", fontsize=18)

        # IO
        if savefig:
            plt.savefig("stacked_full.pdf", dpi=1000)
        if display:
            plt.show()


    def display_detections(self, savefig=False, display=True, style="dark"):
        """
        Displays an image of the stacked full exposure image for this observation with the 
        source labeled with a circular region and all observations for this specific 
        observation labeled.
        
        Arguments
        ---------
        savefig : bool, optional
            A flag determining whether to save the plot
        display : bool, optional
            A flag determining whether to display the plot in an interactive window
        style : str, optional
            A string describing the image display style. Options include ``dark`` 
            which sets a black facecolor, and ``bright`` which sets a white 
            facecolor
        """
        
        # Import current detections
        NoneType = type(None)
        detections, flag, filey = self.read_final_detections()
        
        # If no detections are found, run the default image display method
        if type(detections) == NoneType:
            self.display_image()
            return
        if len(detections["INDEX"]) == 0:
            self.display_image()
            return
        
        # Initialize figure and plot elements
        fig = plt.figure()
        ax = fig.add_subplot(projection=self.wcs)
        stacked_data = fits.getdata(os.path.join(self._impath, f"nu{self._seqid}_{2.96}-{78.96}_keV_stacked.fits"))
        im = ax.imshow(stacked_data, origin='lower', norm=matplotlib.colors.LogNorm())
        center = Point(self._source_pix_coordinates[0][0], self._source_pix_coordinates[0][1]).buffer((self.rlimit) / 2.46, resolution=1000)
        plot_polygon(center, ax=ax, add_points=False, color="yellow")
        
        # Plot all detections
        det_pos_x = []
        det_pos_y = []
        det_pos_prob = []
        for idx, val in enumerate(detections['XPIX']):
            x_pix = detections['XPIX'][idx]
            y_pix = detections['YPIX'][idx]
            
            # Sanity check for main source mask
            if np.sqrt((float(x_pix) - float(self._source_pix_coordinates[0][0]))**2 + (float(y_pix) - float(self._source_pix_coordinates[0][1]))**2) > (self.rlimit) / 2.5:
                det_pos_x.append(float(x_pix))
                det_pos_y.append(float(y_pix)) 
                det_pos_prob.append(float(detections['PROB'][idx]))
        det_scatter = ax.scatter(det_pos_x, det_pos_y, marker="d", c=det_pos_prob, s=120, linewidths=1, edgecolors= "black", cmap="spring", norm=matplotlib.colors.LogNorm()) 
        #cb = plt.colorbar(det_scatter, orientation="vertical")
        plt.title("Probability", x=1.11)
        plt.grid(True)
        plt.xlim(250, 750)
        plt.ylim(250, 750)
        
        # Apply style parameters
        if style == "dark":
            fig.patch.set_facecolor("black")
            ax.set_facecolor("gray")
            ax.tick_params(color="white", labelcolor="white")
            ax.tick_params(axis="both", which="major", labelsize=15)
            im.axes.tick_params(color="white", labelcolor="white")
            #cb.set_label('Probability', color='white', fontsize=15)
            #cb.ax.yaxis.set_tick_params(color='white')
            #cb.outline.set_edgecolor('white')
            #plt.setp(plt.getp(cb.ax.yaxis, 'ticklabels'), color='white', fontsize=15)
            plt.xlabel("Right Ascension", color="white", fontsize=18)
            plt.ylabel("Declination", color="white", fontsize=18)
            
        if style == "bright":
            fig.patch.set_facecolor("white")
            ax.set_facecolor("white")
            ax.tick_params(color="black", labelcolor="black")
            ax.tick_params(axis="both", which="major", labelsize=15)
            im.axes.tick_params(color="black", labelcolor="white")
            #cb.set_label('Probability', color="black", fontsize=15)
            #cb.ax.yaxis.set_tick_params(color="black")
            #cb.outline.set_edgecolor("black")
            #plt.setp(plt.getp(cb.ax.yaxis, 'ticklabels'), color="black", fontsize=15)
            plt.xlabel("Right Ascension", color="black", fontsize=18)
            plt.ylabel("Declination", color="black", fontsize=18)

        # IO
        if savefig:
            savepath = os.path.relpath(os.path.join(self._outpath, "final_detections.pdf"))
            plt.savefig(savepath, dpi=1000)
        if display:
            plt.show()

    # Analysis methods

    def estimate_background(self):
        """
        This method estimates the background for an observation using the DET1 coordinate system. It uses an 
        annulus around the main source mask to estimate the quiessant background rate.
        """
        from astropy.io.fits import getdata
        
        background_rates = {}
        effective_exposure_time = self.calculate_effective_exposure()
        for channel in self._phi_channels:
            # Fine R_90 for FPMA
            elow = round(self.ns.channel_to_energy(float(channel[0])), 3)
            ehigh = round(self.ns.channel_to_energy(float(channel[1])), 3)
            if f"{elow}-{ehigh}" not in background_rates:
                background_rates[f"{elow}-{ehigh}"] = []
            det1_a_data = getdata(self._refpath + f"photometry/nu{self._seqid}A01_cl_{elow}to{ehigh}keV_det1.fits")
            total_aper = RectangularAperture([180.0, 180.0], 330.0, 330.0)
            total_sums, total_error = total_aper.do_photometry(det1_a_data)
            count_ratios = []
            min_deviance = 1
            optimal_radius = 0
            det_x = self._source_det1_coords['A'][0]
            det_y = self._source_det1_coords['A'][1]
            for i in tqdm(np.arange(30, 290, 1)):
                source_aper = CircularAperture(self._source_det1_coords['A'], i)
                src_sums, src_error = source_aper.do_photometry(det1_a_data)
                count_ratio = src_sums / total_sums
                difference = np.abs(0.93 - count_ratio)
                if difference < min_deviance:
                    optimal_radius = i
                    min_deviance = difference
                count_ratios.append((i, count_ratio))
            
            filepath = self._refpath + f"photometry/nu{self._seqid}A01_cl_{elow}to{ehigh}keV_det1_bkg.pdf"
            background_area = calculate_background_area(det1_a_data, det_x, det_y, optimal_radius, filepath=filepath)
            arcsecond_area = background_area * 6.0516 * u.arcsecond * u.arcsecond

            bkg_ann = CircularAnnulus(self._source_det1_coords['A'], optimal_radius, optimal_radius + 30)
            bkg_count, bkg_error = bkg_ann.do_photometry(det1_a_data)
            if bkg_count == 0:
                bkg_ann = CircularAnnulus(self._source_det1_coords['A'], optimal_radius - 50, optimal_radius + 30 - 50)
            bkg_count = bkg_count * u.ct
            area_count_rate = bkg_count / arcsecond_area
            final_count_rate = area_count_rate / effective_exposure_time / u.second
            background_rates[f"{elow}-{ehigh}"].append(final_count_rate)


        for channel in self._phi_channels:
            # Fine R_90 for FPMB
            elow = round(self.ns.channel_to_energy(float(channel[0])), 3)
            ehigh = round(self.ns.channel_to_energy(float(channel[1])), 3)
            det1_a_data = getdata(self._refpath + f"photometry/nu{self._seqid}B01_cl_{elow}to{ehigh}keV_det1.fits")
            total_aper = RectangularAperture([180.0, 180.0], 330.0, 330.0)
            total_sums, total_error = total_aper.do_photometry(det1_a_data)
            count_ratios = []
            min_deviance = 1
            optimal_radius = 0
            det_x = self._source_det1_coords['B'][0]
            det_y = self._source_det1_coords['B'][1]
            for i in tqdm(np.arange(30, 290, 1)):
                source_aper = CircularAperture(self._source_det1_coords['B'], i)
                src_sums, src_error = source_aper.do_photometry(det1_a_data)
                count_ratio = src_sums / total_sums
                difference = np.abs(0.93 - count_ratio)
                if difference < min_deviance:
                    optimal_radius = i
                    min_deviance = difference
                count_ratios.append((i, count_ratio))
            
            filepath = self._refpath + f"photometry/nu{self._seqid}B01_cl_{elow}to{ehigh}keV_det1_bkg.pdf"
            background_area = calculate_background_area(det1_a_data, det_x, det_y, optimal_radius, filepath=filepath)
            arcsecond_area = background_area * 6.0516 * u.arcsecond * u.arcsecond

            bkg_ann = CircularAnnulus(self._source_det1_coords['B'], optimal_radius, optimal_radius + 30)
            bkg_count, bkg_error = bkg_ann.do_photometry(det1_a_data)
            if bkg_count == 0:
                bkg_ann = CircularAnnulus(self._source_det1_coords['B'], optimal_radius - 50, optimal_radius + 30 - 50)
            bkg_count = bkg_count * u.ct
            area_count_rate = bkg_count / arcsecond_area
            final_count_rate = area_count_rate / self.exposure['B01'] / u.second
            background_rates[f"{elow}-{ehigh}"].append(final_count_rate)
        
        rate_string = ""
        for bound in background_rates:
            for val in background_rates[bound]:
                rate_string += f"{val.value} "
        with open(self._refpath + "photometry/background_rates.txt", 'w') as file:
            file.write(rate_string)


    def sliding_cell_detection(self):
        """
        The base sliding cell detection call on the temporal binned data.
        """
        from astropy.io.fits import getdata

        # Display terminal
        print('#' * 90)
        print(f"Performing Sliding Cell Source Detection Search for SEQID: {self._seqid}")
        print(f"Time scale: {self._dtime} seconds.")
        #print(f"Optimized Source Radius: {self.rlimit}")
        print(f"Source Position: {self._source_pix_coordinates[0][0]}, {self._source_pix_coordinates[0][1]}")
        print('#' * 90)

        #if f"{self._dtime}_flag.txt" in self._evtcontents:
        #    print("Old detection pass detected! Skipping!")
        #    return 0
        if f"{self._dtime}_ximage_flag.txt" in self._evtcontents:
            print("New detection pass detected! Skipping!")
            return 0
        
        for channel in self._phi_channels:
            
            # Establish file paths and directories
            channel_path = os.path.relpath(os.path.join(self._detpath, f"{channel[0]}-{channel[1]}_{self._dtime}-{self._snr}"))
            os.chdir(channel_path)
            evt_files = glob.glob("*.evt")
            
            # Remove old detection files
            det_files = glob.glob("*.det")
            fits_files = glob.glob("*.fits")
            for file in det_files:
                os.remove(file)
            for file in fits_files:
                os.remove(file)

            print(f"Processing PHI Channels: {channel[0]}-{channel[1]}")
            for file in tqdm(evt_files):
                if len(getdata(file)) != 0:
                    filestub = file.replace('.evt', '')
                    logfile = os.path.relpath(os.path.join(self._logpath, f"slidecell_{filestub}.log"))
                    with open("ximage.xco", 'w') as script:
                        script.write(f"read/fits/size=800/{file}\n")
                        script.write(f"detect/snr={self._snr}/source_box_size=8/filedet={filestub}.det/fitsdet={filestub}.fits\n")
                        script.write("exit")
                    os.system(f"ximage @ximage.xco > {logfile}")
                    os.system(f"rm -r -f ximage.xco")
        
        # Create completion flag
        os.chdir(self._evtpath)
        with open(f"{self._dtime}_ximage_flag.txt", 'w') as file:
            file.write("PROCESSING COMPLETE")
        os.chdir(self._mainpath)


    def detection_merging(self):
        """
        Performs the first round of detection redundancy culling on detections caught by the sliding cell algorithm.
        """

        for channel in self._phi_channels:
            detpath = os.path.relpath(os.path.join(self._detpath, f"{channel[0]}-{channel[1]}_{self._dtime}-{self._snr}"))
            if not os.path.isdir(detpath):
                continue
            os.chdir(detpath)
            detfiles = glob.glob("*.det")
            
            # Construct xselect script string
            script_string = f"srcmrg/out=mrg.txt/tolerance=13" 
            for file in tqdm(detfiles):
                script_string += f" {file}"
            script_string += "\n"
            
            if len(detfiles) == 0:
                continue
            else:
                logfile = os.path.relpath(os.path.join(self._logpath, f"merge_{channel[0]}-{channel[1]}_{self._dtime}-{self._snr}.log"))
                with open("src_merge.xco", 'w') as script:
                    script.write(script_string)
                    script.write("exit\n")
                os.system(f"ximage @src_merge.xco > {logfile}")
                os.system(f"rm -r -f src_merge.xco")
        
        # Create completion flag
        os.chdir(self._evtpath)
        with open(f"{self._dtime}_merging_flag.txt", 'w') as file:
            file.write("PROCESSING COMPLETE")
        os.chdir(self._mainpath)


    def run_detection_pipeline(self):
        """
        Runs the event extraction, sliding cell detection, and detection merging steps for an observation.
        """

        self.event_extraction()
        self.sliding_cell_detection()
        self.detection_merging()


    def detection_summary(self):
        """
        Provides a brief text summary on the detection results for each PI channel.
        """
        for channel in self._phi_channels:
            detections = self.detection_dir_processing(channel)
            if detections == None:
                n_det = 0
            if detections != None:
                n_det = len(detections["INDEX"]) 
            print(f"PI Channels: {channel[0]}-{channel[1]} -- {n_det} detections found.")


    def process_detections(self):
        """
        Calls ``nuproducts`` on detections.
        """
        for channel in self._phi_channels:
            detections = self.detection_dir_processing(channel)
            if detections == None:
                n_det = 0
            if detections != None:
                n_det = len(detections["INDEX"]) 
                self.nuproducts(detections, channel)
            print(f"PI Channels: {channel[0]}-{channel[1]} -- {n_det} detections found.")


    def read_unique_detections(self, bounds):
        """
        This method reads in the merged detection instances and populates a dictionary with all of the 
        recorded values for future use.

        Arguments
        ---------
        None

        Returns
        -------
        dict
            A dictionary containing the following keys: `INDEX`, `COUNTS`, `XPIX`, `YPIX`, `VIGCOR`, 
            `RA`, `DEC`, `ERR`, `HBOX`, `PROB`, `SNR`. A description of these keys is below:

            INDEX : the number associated with this detection as a sort of 'detection id'
            COUNTS : the counts associated with this detection
            XPIX : the x pixel on which this detection is centered
            YPIX : the y pixel on which this detection is centered
            VIGCOR : 
            RA : the right ascension for this detection
            DEC : the declination for this detection
            ERR : the count error for this detection
            HBOX : 
            PROB : 
            SNR : the snr calculated for this detection
            
        None
            In the event that the unique merger file is not found.

        RAISES
        ------
        None
        """

        detpath = os.path.join(self._detpath, f"{bounds[0]}-{bounds[1]}_{self._dtime}-{self._snr}")
        filepath = os.path.join(detpath, "mrg.txt")
        
        if not os.path.isdir(detpath):
            return None
        if not os.path.isfile(filepath):
            return None           
        
        with open(filepath) as detections:
            detect_info = {}
            detect_info["INDEX"] = []
            detect_info["COUNTS"] = []
            detect_info["XPIX"] = []
            detect_info["YPIX"] = []
            detect_info["VIGCOR"] = []
            detect_info["RA"] = []
            detect_info["DEC"] = []
            detect_info["ERR"] = []
            detect_info["HBOX"] = []
            detect_info["PROB"] = []
            detect_info["SNR"] = []
            
            for line in detections:
                line_info = line.split()
                detect_info["INDEX"].append(line_info[0])
                detect_info["COUNTS"].append(line_info[1])
                detect_info["XPIX"].append(line_info[2])
                detect_info["YPIX"].append(line_info[3])
                detect_info["VIGCOR"].append(line_info[4])
                detect_info["RA"].append(f"{line_info[5]} {line_info[6]} {line_info[7]}")
                detect_info["DEC"].append(f"{line_info[8]} {line_info[9]} {line_info[10]}")
                detect_info["ERR"].append(line_info[11])
                detect_info["HBOX"].append(line_info[12])
                detect_info["PROB"].append(line_info[13])
                detect_info["SNR"].append(line_info[14])
        
        return detect_info


    def read_detection_dir(self, bounds):
        """
        Reads in the individual detection files for the given PI bounds. The detections are then aggregated into an 
        Astropy table and returned.
        """
        detpath = os.path.relpath(os.path.join(self._detpath, f"{bounds[0]}-{bounds[1]}_{self._dtime}-{self._snr}"))
        det_files = glob.glob(os.path.join(detpath, "*.det"))
        
        if len(det_files) == 0:
            return None
        
        detect_info = {}
        
        for file in det_files:
            with open(file) as detections:
                for i in range(14):
                    detections.readline()
                times = tuple(file.replace(detpath, '').replace('/', '').replace(".det", '').split(sep="-"))
                detect_info[times] = {}
                detect_info[times]["INDEX"] = []
                detect_info[times]["COUNTS"] = []
                detect_info[times]["XPIX"] = []
                detect_info[times]["YPIX"] = []
                detect_info[times]["VIGCOR"] = []
                detect_info[times]["RA"] = []
                detect_info[times]["DEC"] = []
                detect_info[times]["ERR"] = []
                detect_info[times]["HBOX"] = []
                detect_info[times]["PROB"] = []
                detect_info[times]["SNR"] = []
                detect_info[times]["FILE"] = []

                for line in detections:
                    line_info = line.split()
                    detect_info[times]["FILE"].append(file)
                    detect_info[times]["INDEX"].append(line_info[0])
                    detect_info[times]["COUNTS"].append(line_info[1])
                    detect_info[times]["XPIX"].append(line_info[2])
                    detect_info[times]["YPIX"].append(line_info[3])
                    detect_info[times]["VIGCOR"].append(line_info[4])
                    detect_info[times]["RA"].append(f"{line_info[5]} {line_info[6]} {line_info[7]}")
                    detect_info[times]["DEC"].append(f"{line_info[8]} {line_info[9]} {line_info[10]}")
                    detect_info[times]["ERR"].append(line_info[11])
                    detect_info[times]["HBOX"].append(line_info[12])
                    detect_info[times]["PROB"].append(line_info[13])
                    detect_info[times]["SNR"].append(line_info[14])

        return detect_info
    
    
    def extract_detections(self):
        """
        Runs the `nuproducts` extraction procedure for all unique detections.
        """

        for channel in self._phi_channels:
            self.nuproducts(self.detection_dir_processing(channel), channel)

    
    def remove_main_source(self, detect_info):
        """
        Removes all detections which where found to be within ``self.rlimit`` arcseconds of the main source.
        """

        trimmed_detect_info = {}        
        for key in detect_info:
            trimmed_detect_info[key] = []
        
        for i in detect_info["INDEX"]:
            ra = detect_info["RA"][int(i) - 1]
            dec = detect_info["DEC"][int(i) - 1]
            detect_position = SkyCoord(f"{ra} {dec}", unit=(u.hourangle, u.deg), frame='fk5')
            if self._source_position.separation(detect_position).arcsec > self.rlimit:
                for key in detect_info:
                    trimmed_detect_info[key].append(detect_info[key][int(i) - 1])
                    
        if len(trimmed_detect_info["INDEX"]) == 0:
            return None
        return trimmed_detect_info
    
    
    def detection_dir_processing(self, bounds):
        """ 
        Trims down set of detections via removing duplicates and detections which 
        are within the main source mask for this target.
        """
        # Read in unique detection information, return None if not present
        unique_detect_info = self.read_unique_detections(bounds)
        if unique_detect_info == None:
            return None

        # Read in all detections for a given PI channel range, return None if not present
        all_detect_info = self.read_detection_dir(bounds)
        if all_detect_info == None:
            return None
        
        # Eliminated detections within main source, if none are left, return None
        trimmed_detect_info = self.remove_main_source(unique_detect_info)
        if trimmed_detect_info == None:
            return None

        # Now remove duplicates and save time slots
        trimmed_all_info = {}
        for key in trimmed_detect_info:
            trimmed_all_info[key] = []
        trimmed_all_info["TIMES"] = []

        for time in all_detect_info:
            for idx in all_detect_info[time]["INDEX"]:

                if all_detect_info[time]["XPIX"][int(idx) - 1] in trimmed_detect_info["XPIX"]:
                    for key in trimmed_detect_info:
                        trimmed_all_info[key].append(all_detect_info[time][key][int(idx) - 1])
                    trimmed_all_info["TIMES"].append(time)
        
        self._detections = trimmed_all_info
        return trimmed_all_info
    

    def write_net_detections(self):
        """ 
        Writes all screened detections to a file.
        """
        
        print("Merging all valid detections")

        # Loop through all channels
        for channel in self._phi_channels:
            trimmed_all_info = self.detection_dir_processing(channel)
            if trimmed_all_info == None:
                print("No valid ximage detections found!")
                return None
            if trimmed_all_info != None:
                n_obj = len(trimmed_all_info["INDEX"])

        reduced_detections, tkeys = self.basic_detection_stats()
        
        if reduced_detections == None:
            return None
        
        # Populate dictionary with aggregated information
        trimmed_all_info = {}
        if reduced_detections != None:
            for key in tkeys:
                trimmed_all_info[key] = []
            key_map = {}
            for idx in range(len(tkeys)):
                key_map[idx] = tkeys[idx]
            for idx in range(len(tkeys)):
                trimmed_all_info[tkeys[idx]] = []
            for det in range(len(reduced_detections)):
                for idx in range(len(tkeys)):
                    trimmed_all_info[tkeys[idx]].append(reduced_detections[det][idx])

            # Applying table indexing corrections.
            n_obj = len(trimmed_all_info["RA"])
            trimmed_all_info["INDEX"] = list(range(1, n_obj + 1))
            
            # Write table to an Astropy table object
            detect_table = Table()
            for key in trimmed_all_info:
                detect_table[key] = trimmed_all_info[key]
            detect_table.write(os.path.join(self._detpath, f"{self._dtime}_detections.tbl"), format='ipac', overwrite=True)
        
        # Create completion flag
        os.chdir(self._evtpath)
        with open(f"{self._dtime}_netdetection_flag.txt", 'w') as file:
            file.write("PROCESSING COMPLETE")
        os.chdir(self._mainpath)



    def nuproducts(self, detect_info, pi_bounds):
        """
        Low level method for running the command lines necessary for calling 
        `nuproducts`
        """
        if detect_info == None:
            return None

        detect_info = self._detections
        generate_directory(self._out_path)
        gti_files = glob.glob(self._refpath + f"products/*.hnd_gti")
        
        for times in set(detect_info["TIMES"]):
            if self._refpath + f"products/{times[0]}_{times[1]}_gti.hnd_gti" not in gti_files:
                generate_gti_files(self._refoutpath, times[0], times[1])
            

        for idx in range(len(detect_info["INDEX"])):
            coords = self.ra_dec_todeg(detect_info["RA"][int(idx)], detect_info["DEC"][int(idx)])
            string = f"nuproducts indir=./event_cl instrument=FPMA steminputs=nu{self._seqid} outdir=./products/{pi_bounds[0]}_{pi_bounds[1]}-{self._dtime}_{self._snr}/{idx}A srcra={coords[0]} srcdec={coords[1]} bkgra={coords[0]} bkgdec={coords[1]} binsize=1 usrgtifile=./products/{detect_info['TIMES'][idx][0]}_{detect_info['TIMES'][idx][1]}_gti.hnd_gti infile=science.evt"
            subprocess.run(string.split(), cwd=self._refpath)
            string = f"nuproducts indir=./event_cl instrument=FPMB steminputs=nu{self._seqid} outdir=./products/{pi_bounds[0]}_{pi_bounds[1]}-{self._dtime}_{self._snr}/{idx}B srcra={coords[0]} srcdec={coords[1]} bkgra={coords[0]} bkgdec={coords[1]} binsize=1 usrgtifile=./products/{detect_info['TIMES'][idx][0]}_{detect_info['TIMES'][idx][1]}_gti.hnd_gti infile=science.evt"
            subprocess.run(string.split(), cwd=self._refpath)


    def ra_dec_todeg(self, ra, dec):
        """
        A simple conversion function which converts from RA DEC convention to degrees:

        Arguments
        ---------
        ra : string

        dec : string

        Returns
        -------
        tuple
            A tuple containing the degree values of `(RA, DEC)`

        Raises
        ------
        ValueError
            if an unacceptable RA or DEC string format is inputted.
        """
        c = SkyCoord(f"{ra} {dec}", unit=(u.hourangle, u.deg))
        return (c.ra.deg, c.dec.deg)


    def ds9_detections(self, mod='B'):
        """
        Displays all of the detections that were processed by nuproducts into an instance of ds9.
        """
        ds9_string = f"ds9 nu{self._seqid}{mod}01_cl.evt "
        detection_dirs = os.listdir(self._refpath + "products/")
        for dir in detection_dirs:
            ds9_string += f"-regions ../products/{dir}/nu{self._seqid}{mod}01_bkg.reg "
            ds9_string += f"-regions ../products/{dir}/nu{self._seqid}{mod}01_src.reg "            
        subprocess.run(ds9_string.split(), cwd=self._refpath + "event_cl/")
        return ds9_string
    

    def list_product_files(self):
        """
        Returns all of the product files 
        """
        all_files = []
        for channel in self._phi_channels:
            files = glob.glob(self._refpath + f"products/{channel[0]}_{channel[1]}-{self._dtime}_{self._snr}/**/*.flc", recursive=True)
            for file in files:
                all_files.append(file)
        return all_files


    def acquire_lightcurve(self, ids, tstart, tstop, xpix, ypix):
        """
        Returns a dictionary containing the background-subtracted lightcurve data.

        Arguments
        ---------
        n : `int`
            The index of the detection.
        
        Returns
        -------

        Raises
        ------
        """

        tot_data = fits.getdata(self._refpath + "science.evt")
        events = []
        times = []
        for datum in tot_data:
            if float(datum[0]) > float(tstart) and float(datum[0]) < float(tstop):
                datum_xpix = datum[13]
                datum_ypix = datum[14]
                if np.sqrt((float(xpix) - float(datum_xpix))**2 + (float(ypix) - float(datum_ypix))**2) < self.rlimit / 2.5:
                    events.append(datum)
                    times.append(float(datum[0]))
        """lc_bins = np.arange(float(tstart), float(tstop), 5)
        lc, bines = np.histogram(times, lc_bins)
        t = []
        for idx, l in enumerate(lc):
            t.append(idx)
        t = np.array(t) * 10
        plt.plot(t, lc, ms= 1)
        plt.savefig(self._refpath + f"products/{ids}lc_5.pdf")
        plt.close()

        lc_bins = np.arange(float(tstart), float(tstop), 10)
        lc, bines = np.histogram(times, lc_bins)
        t = []
        for idx, l in enumerate(lc):
            t.append(idx)
        t = np.array(t) * 10
        plt.plot(t, lc, ms= 1)
        plt.savefig(self._refpath + f"products/{ids}_10lc.pdf")
        plt.close()"""

        lc_bins = np.arange(float(tstart), float(tstop), 20)
        lc, bines = np.histogram(times, lc_bins)
        t = []
        for idx, l in enumerate(lc):
            t.append(idx)
        t = np.array(t) * 10
        plt.plot(t, lc, ms= 1)
        plt.savefig(self._refpath + f"products/{ids}_20lc.pdf")
        plt.close()

        lc_bins = np.arange(float(tstart), float(tstop), 50)
        lc, bines = np.histogram(times, lc_bins)
        t = []
        for idx, l in enumerate(lc):
            t.append(idx)
        t = np.array(t) * 10
        plt.plot(t, lc, ms= 1)
        plt.savefig(self._refpath + f"products/{ids}_50lc.pdf")
        plt.close()

    
    def acquire_event_curves(self):
        """ 
        Low level method for accessing lightcurves for detections
        """
        detections, flag, filey = self.read_final_detections()
        if len(detections["INDEX"]) != 0:
            for idx in range(len(detections['INDEX'])):
                ids = detections['INDEX'][idx]
                xpix = detections['XPIX'][idx]
                ypix = detections['YPIX'][idx]
                tstart = detections['TSTART'][idx]
                tstop = detections['TSTOP'][idx]
                self.acquire_lightcurve(ids, tstart, tstop, xpix, ypix)


    def plot_lightcurves(self):
        """
        Plots the background subtracted lightcurve of all detections.

        Arguments
        ---------
        n : `int`
            The index of the detection being plotted.

        Returns
        -------

        Raises
        ------
        """

        for file in self.list_product_files():
            lc_data = self.acquire_lightcurve(file)
            fig, ax = plt.subplots()
            ax.errorbar(lc_data["TIME"], lc_data["RATE"], lc_data["RATEERR"], fmt='o', lw=0.5, markersize=2)
            ax.set_xlabel("Time (s)")
            ax.set_ylabel("Counts/Time (c/s)")
            #pilow = file.replace()
            #pihi = file.replace()
            #kev_low = chan_to_energy(pilow)
            #kev_high = chan_to_energy(pihi)
            ax.set_title("Background Subtracted Lightcurve", loc="left")
            #ax.set_title(f"Energy Range: {kev_low} - {kev_high}", loc="right")
            plt.savefig(file.replace(".flc", ".pdf"), dpi=800)
            plt.close()


    def basic_detection_stats(self):
        """
        Performs some simple statistics on each detection and corrects some Astropy table formatting
        """
        
        total_detections = {}
        for channel in self._phi_channels:
            trimmed_all_info = self.detection_dir_processing(channel)
            
            if trimmed_all_info == None:
                print("basic_detection_stats: No valid detections inputted!")
                continue
            if trimmed_all_info != None:
                
                # Reindexing all detections
                n_obj = len(trimmed_all_info["INDEX"])
                trimmed_all_info["INDEX"] = list(range(1, n_obj + 1))
                
                # Formatting time ranges into separate columns
                t_starts = [val[0].replace("nu_", '') for val in trimmed_all_info["TIMES"]]
                t_stops = [val[1] for val in trimmed_all_info["TIMES"]]
                trimmed_all_info["TSTART"] = t_starts
                trimmed_all_info["TSTOP"] = t_stops
                
                # Formatting counts and count errors into separate columns
                count_vals = [float(counts.split('+/-')[0]) for counts in trimmed_all_info["COUNTS"]]
                count_err_vals = [float(counts.split('+/-')[1]) for counts in trimmed_all_info["COUNTS"]]
                trimmed_all_info["COUNTS"] = count_vals
                trimmed_all_info["COUNTSERR"] = count_err_vals
                
                # Adding seqid and PI channel information for each detection
                trimmed_all_info["SEQID"] = [self._seqid for i in range(n_obj)]
                trimmed_all_info["BOUND"] = [f"{channel[0]}-{channel[1]}" for i in range(n_obj)]
                
                # Remove deprecated column
                del trimmed_all_info["TIMES"]
                
                # Calculate separation distance from main source
                seps = []
                for i in range(len(trimmed_all_info["INDEX"])):
                    ra = trimmed_all_info['RA'][i]
                    dec = trimmed_all_info['DEC'][i]
                    detect_position = SkyCoord(f"{ra} {dec}", unit=(u.hourangle, u.deg), frame='fk5')
                    seps.append(self._source_position.separation(detect_position).arcsec)
                trimmed_all_info["SEP"] = np.array(seps)
                
                # Adding source position information for future use
                trimmed_all_info["SRCRA"] = [self.source_position.ra.deg for _ in range(n_obj)]
                trimmed_all_info["SRCDEC"] = [self.source_position.dec.deg for _ in range(n_obj)]
                
                # Acquire all keys in final dictionary
                tkeys = []
                for key in trimmed_all_info:
                    tkeys.append(key)
                
                # Aggregating detection details across all PI channels
                for idx in tqdm(range(len(trimmed_all_info["INDEX"]))):
                    values = []
                    for key in trimmed_all_info:
                        values.append(trimmed_all_info[key][idx])
                    total_detections[len(total_detections)] = values

                    
        if len(total_detections.keys()) == 0:
            return None, None
        return total_detections, tkeys


    def slide_cell_verification(self, det_xpix, det_ypix, tstart, tstop, channel):
        
        ## FPMA PASS
        os.chdir(self._evtpath)
        indpath = os.path.relpath(os.path.join(self._detpath, f"{channel[0]}-{channel[1]}_{self._dtime}-{self._snr}/individual"))
        generate_directory(indpath, overwrite=False)
        
        
        # Generate the gti filtered files.  
        logfile = os.path.relpath(os.path.join(self._logpath, f"individual_{channel[0]}-{channel[1]}_{self._dtime}-{self._snr}_{tstart}-{tstop}A"))      
        outfile = os.path.relpath(os.path.join(self._detpath, f"{channel[0]}-{channel[1]}_{self._dtime}-{self._snr}/individual/nu_{tstart}-{tstop}A.evt"))
        with open("xselect.xco", 'w') as script:
            script.write(f'{self._sessionid}\n')
            script.write(f"read events\n")
            script.write(".\n")
            script.write(f"nu{self._seqid}A01_cl.evt\n")
            script.write('yes\n')
            script.write("filter PHA_CUTOFF\n")
            script.write(f"{channel[0]}\n")
            script.write(f"{channel[1]}\n")
            script.write("filter time scc\n")
            script.write(f"{tstart} , {tstop}\n")
            script.write("x\n")
            script.write("extract events\n")
            script.write("\n")
            script.write(f"save events {outfile}\n")
            script.write('no\n')
            script.write("exit no\n")
        os.system(f"xselect @xselect.xco > {logfile}")
        os.system("rm xselect.xco")

        # FPMB PASS
        logfile = os.path.relpath(os.path.join(self._logpath, f"individual_{channel[0]}-{channel[1]}_{self._dtime}-{self._snr}_{tstart}-{tstop}B"))      
        outfile = os.path.relpath(os.path.join(self._detpath, f"{channel[0]}-{channel[1]}_{self._dtime}-{self._snr}/individual/nu_{tstart}-{tstop}B.evt"))
        with open("xselect.xco", 'w') as script:
            script.write(f'{self._sessionid}\n')
            script.write(f"read events\n")
            script.write(".\n")
            script.write(f"nu{self._seqid}B01_cl.evt\n")
            script.write('yes\n')
            script.write("filter PHA_CUTOFF\n")
            script.write(f"{channel[0]}\n")
            script.write(f"{channel[1]}\n")
            script.write("filter time scc\n")
            script.write(f"{tstart} , {tstop}\n")
            script.write("x\n")
            script.write("extract events\n")
            script.write("\n")
            script.write(f"save events {outfile}\n")
            script.write('no\n')
            script.write("exit no\n")
        os.system(f"xselect @xselect.xco > {logfile}")
        os.system("rm xselect.xco")
        os.chdir(self._mainpath)

        # Perform sliding cell detection

        indpath = os.path.relpath(os.path.join(self._detpath, f"{channel[0]}-{channel[1]}_{self._dtime}-{self._snr}/individual"))   
        filestub = os.path.relpath(os.path.join(indpath, f"nu_{tstart}-{tstop}"))
        evtfiles = glob.glob(filestub + '*.evt')
        passing_A = False
        passing_B = False
        counts_A = 0
        counts_B = 0
        for file in evtfiles:
            counts = 0
            data_events = fits.getdata(file)
            for event in data_events:
                xpix = event[13]
                ypix = event[14]
                if np.sqrt((float(xpix) - float(det_xpix))**2 + (float(ypix) - float(det_ypix))**2) < self.rlimit / 2.46:
                    counts += 1
            if counts > 0:
                letter = file.replace(filestub, '').replace(".evt", '').replace('/', '')
                if letter == 'A':
                    passing_A = True
                    counts_A = counts
                if letter == 'B':
                    passing_B = True
                    counts_B = counts
        if passing_A and passing_B:
            return True, counts_A, counts_B
        else:
            return False, counts_A, counts_B    

    def further_processing(self):
        tot_data = fits.getdata(self._refpath + "science.evt", 1)
        times = []
        for datum in tot_data:
            times.append(float(datum[0]))
        
        lc_bins = np.arange(np.min(times), np.max(times), 1)
        lc, bines = np.histogram(times, lc_bins)
        t = []
        for idx, l in enumerate(lc):
            t.append(idx)
        t = np.array(t) * 10
        plt.plot(t, lc, ms= 1)
        
        plt.title(f"Exposure time: {self._exposure['A01']}", loc="left")
        plt.show()
        plt.close()

        detections, flag, filey = self.read_final_detections()
        if not flag:
            for idx in range(len(detections['INDEX'])):
                for interval in self.time_bins[0]:
                    phi_low = detections['BOUND'][idx].split('-')[0]
                    phi_high = detections['BOUND'][idx].split('-')[1]
                    #tstart = detections['TSTART'][idx]
                    #tstop = detections['TSTOP'][idx]
                    tstart = self.time_bins[0][interval][0]
                    tstop = self.time_bins[0][interval][1]
                    tot_data = fits.getdata(self._refpath + f"detections/{phi_low}-{phi_high}_{self._dtime}-{self._snr}/nu_{tstart}-{tstop}.evt", 1)
                    times = []
                    for datum in tot_data:
                        times.append(float(datum[0]))
                    
                    lc_bins = np.arange(np.min(times), np.max(times), 1)
                    lc, bines = np.histogram(times, lc_bins)
                    t = []
                    for idx, l in enumerate(lc):
                        t.append(idx)
                    t = np.array(t) * 10
                    plt.plot(t, lc, ms= 1)
                    
                    plt.title(f"Exposure time: {self._exposure['A01']}", loc="left")
                    plt.show()
                    plt.close()


    def recalculate_poisson(self, override=False):
        """ 
        This method is used to recalculate the poisson statistic implemented by 
        `ximage` as there are known issues with the background estimation. These 
        newly filtered candidates are then written to a new table, labeled:
        ./detections/{dtime}_poisson.tbl
        """
        from astropy.io.fits import getdata
        
        # Check if poisson has been done
        if os.path.isfile(os.path.join(self._evtpath, f"{self._dtime}_poisson_flag.txt")):
            return
        # Check if detections exist
        filepath = os.path.relpath(os.path.join(self._detpath, f"{self._dtime}_detections.tbl"))
        if not os.path.isfile(filepath):
            # Return completion flag, there are no detections to do poisson stats on
            os.chdir(self._evtpath)
            with open(f"{self._dtime}_poisson_flag.txt", 'w') as file:
                file.write("PROCESSING COMPLETE")
            os.chdir(self._mainpath)
            return 
        
        # Read in detections
        data_table = Table.read(filepath, format='ipac')
        n_obj = len(data_table['INDEX'])
        passing = 0
        valid_indices = []
        net_counts = []
        pvals = []
        det_xpix = []
        det_ypix = []
        
        self.adata = getdata(self._fpma_eventpath)

        # Look through each index
        for idx in tqdm(range(len(data_table['INDEX']))):
            
            if len(data_table['INDEX']) > 800:
                print("WARNING: Over 800 flagged locations found!")
                print("Straylight region/main source mask failure may be culprit!")
                print("Set override=True to continue processing!")
                if override:
                    pass
                else:
                    return
            
            # Get detector coordinates
            detect_position = SkyCoord(f"{data_table['RA'][idx]} {data_table['DEC'][idx]}", unit=(u.hourangle, u.deg), frame='fk5')
            det_x = detect_position.ra.deg
            det_y = detect_position.dec.deg
            outfile = self.sky_to_detector(detect_position.ra.deg, detect_position.dec.deg)
            data = fits.getdata(outfile)
            mean_detx = []
            mean_dety = []
            for datum in data:
                mean_detx.append(float(datum[1]))
                mean_dety.append(float(datum[2]))
            det_x = np.mean(mean_detx)
            det_y = np.mean(mean_dety)

            # Acquire flagged pixel information
            blow = int(data_table['BOUND'][idx].split('-')[0])
            bhigh = int(data_table['BOUND'][idx].split('-')[1])
            elow = round(self.ns.channel_to_energy(float(blow)), 3)
            ehigh = round(self.ns.channel_to_energy(float(bhigh)), 3)
            tlow = data_table['TSTART'][idx]
            thigh = data_table['TSTOP'][idx]
            datafile = os.path.relpath(os.path.join(self._detpath, f"{blow}-{bhigh}_{self._dtime}-{self._snr}/nu_{tlow}-{thigh}.evt"))
            
            evt_data = getdata(datafile)
            effective_exposure_time = self.calculate_effective_time(evt_data)
            outpath = os.path.relpath(os.path.join(self._impath, f"{blow}-{bhigh}_{self._dtime}-{self._snr}"))
            generate_directory(outpath)
            
            # Perform photometry
            src_img = make_det1_image(datafile, elow=elow, ehigh=ehigh, outpath=outpath)
            det1_data = getdata(os.path.abspath(src_img))
            src_area = calculate_source_area(det1_data, det_x, det_y, 10) * 6.0516 * u.arcsecond * u.arcsecond
            bkg_area = calculate_background_area(det1_data, det_x, det_y, 10) * 6.0516 * u.arcsecond * u.arcsecond

            actual_counts = source_counts(det1_data, det_x, det_y, 10)
            bk_counts = bkg_counts(det1_data, det_x, det_y, 10)
            src_predict = bk_counts * src_area / bkg_area
            
            if np.isnan(actual_counts.value):
                continue    
            if effective_exposure_time <= 0:
                continue
            
            # Conduct poisson test
            res = poisson_means_test(int(actual_counts.value), effective_exposure_time, int(src_predict.value), effective_exposure_time)
            if float(data_table["XPIX"][idx]) < 520 and float(data_table["XPIX"][idx]) > 480:
                if float(data_table["YPIX"][idx]) < 389 and float(data_table["YPIX"][idx]) > 360:
                    print(res.pvalue)
            if res.pvalue < 0.01:
                passing += 1
                valid_indices.append(idx)
                det_xpix.append(det_x)
                det_ypix.append(det_y)
                pvals.append(res.pvalue)
                net_counts.append(actual_counts)
        if n_obj != 0:
            print(f"Percent of objects kept: {passing/n_obj}")
        
        # Rewrite remaining detections
        filtered_detections = {}
        for col in data_table.colnames:
            filtered_detections[col] = []
        for col in data_table.colnames:
            for idx in valid_indices:
                filtered_detections[col].append(data_table[col][idx])
        
        filtered_detections["PVAL"] = pvals
        filtered_detections["NETCOUNTS"] = net_counts
        filtered_detections["DET1X"] = det_xpix
        filtered_detections["DET1Y"] = det_ypix
        
        detect_table = Table()
        for key in filtered_detections:
            detect_table[key] = filtered_detections[key]
        finalpath = os.path.relpath(os.path.join(self._detpath, f"{self._dtime}_poisson.tbl"))
        detect_table.write(finalpath, format='ipac', overwrite=True)

        # Create completion flag
        os.chdir(self._evtpath)
        with open(f"{self._dtime}_poisson_flag.txt", 'w') as file:
            file.write("PROCESSING COMPLETE")
        os.chdir(self._mainpath)


    def calculate_effective_time(self, data, evt=True):
        """ 
        This method returns an estimate for the exposure time of an 
        observation, accounting for geometry-related dead gaps in 
        lightcurves
        """
        times = []
        for datum in data:
            if evt:
                times.append(float(datum[0]))
            else:
                times.append(float(datum))
        tstart = np.min(times)
        tstop = np.max(times)
        lc_bins = np.arange(float(tstart), float(tstop), 20)
        lc, bines = np.histogram(times, lc_bins)
        count = 0
        for idx, l in enumerate(lc):
            if l != 0:
                count += 1
        return count * 20
    
    
    def detection_lightcurve(self, xpix, ypix, low_pi, high_pi, tstart, tstop):
        """ 
        Low level method for plotting lightcurve data.
        """
        from astropy.io.fits import getdata
        dt = 30

        evt_data = getdata(os.path.join(self._mainpath, "science.evt"))
        times = []
        for datum in evt_data:
            x = float(datum[13])
            y = float(datum[14])
            if np.sqrt((x - xpix) ** 2 + (y - ypix)**2) <= 15:
                times.append(float(datum[0]))
        totstart = np.min(times)
        totstop = np.max(times)
        lc_bins = np.arange(float(totstart), float(totstop), dt)
        lc, bines = np.histogram(times, lc_bins)
        t = []
        count = 0
        for idx, l in enumerate(lc):
            t.append(idx)
            if l != 0:
                count += 1
        t = np.array(t) * dt

        intervals = []
        for i in range(len(lc_bins) - 1):
            intervals.append((lc_bins[i], lc_bins[i + 1]))
        count_rates = {}
        for interval in intervals:
            count_rates[interval] = []
        for time in times:
            for interval in intervals:
                if time > interval[0] and time <= interval[1]:
                    count_rates[interval].append(time)
        
        rates = []
        for interval in count_rates:
            if len(count_rates[interval]) == 0:
                rates.append(0)
                continue
            eff_time = self.calculate_effective_time(count_rates[interval], evt=False)
            if eff_time == 0:
                rates.append(0)
            else:
                count_rate = len(count_rates[interval]) / eff_time
                rates.append(count_rate)
            
        fig, axs = plt.subplots(nrows=2, ncols=1)
        axs[0].plot(t, lc, ms= 1)
        axs[0].axvspan(tstart - totstart, tstop - totstart, alpha=0.5, color="red")
        axs[0].set_xlabel("Elapsed Time (sec)")
        axs[0].set_ylabel("Total Counts")
        axs[0].set_title(f"{round(xpix)}, {round(ypix)}")
        axs[1].plot(t, rates, ms= 1)
        axs[1].axvspan(tstart - totstart, tstop - totstart, alpha=0.5, color="red")
        axs[1].set_xlabel("Elapsed Time (sec)")
        axs[1].set_ylabel("Adjusted Count Rate")
        plt.show()
        plt.close()


    
    def plot_detection_lightcurves(self, override=False):
        """ 
        User-interface method for mass plotting lightcurves
        """
        filepath = self._refpath + f"detections/{self._dtime}_3.tbl"
        data_table = Table.read(filepath, format='ipac')
        if len(data_table['INDEX']) > 5:
            if override:
                pass 
            else:
                print(f"WARNING: {len(data_table['INDEX'])} candidates found!")
                print("Set override=True to plot them")
                return
        for idx in range(len(data_table['INDEX'])):
            xpix = float(data_table['XPIX'][idx])
            ypix = float(data_table['YPIX'][idx])
            lbound = data_table['BOUND'][idx].split('-')[0]
            hbound = data_table['BOUND'][idx].split('-')[1]
            tstart = float(data_table['TSTART'][idx])
            tstop = float(data_table['TSTOP'][idx])
            self.detection_lightcurve(xpix, ypix, lbound, hbound, tstart, tstop)


    def classify_sources(self):
        """
        This method applies a simple algorithm to classify persistent vs 
        transient sources.
        """
        
        print("#" * 90)
        print(f"Measuring the persistence value for poisson verified sources.")
        print("#" * 90)
        detections, flag, filey = self.read_final_detections()
        
        NoneType = type(None)
        if type(detections) == NoneType:
            return
        if len(detections["INDEX"]) == 0:
            return
        
        tpersist = []
        cpersist = []
        for index, value in enumerate(detections['INDEX']):
            
            tstarts = []
            tstops = []
            channels = []
            xpix = float(detections["XPIX"][index])
            ypix = float(detections["YPIX"][index])
            tstarts.append(detections["TSTART"][index])
            tstops.append(detections["TSTOP"][index])
            channels.append(detections["BOUND"][index])
            for index2, value2 in enumerate(detections['INDEX']):
                if np.sqrt((float(detections["XPIX"][index2]) - xpix) ** 2 + (float(detections["YPIX"][index2]) - ypix) ** 2) <= 8:
                    tstarts.append(detections["TSTART"][index2])
                    tstops.append(detections["TSTOP"][index2])
                    channels.append(detections["BOUND"][index2])
            present_percent = len(set(tstarts)) / (2 * self.n_cuts)
            channels_percent = len(set(channels)) / 4
            tpersist.append(present_percent)
            cpersist.append(channels_percent)
        detections['TPERSIST'] = tpersist
        detections['CPERSIST'] = cpersist
        detect_table = Table(detections)
        detect_table.write(os.path.join(self._detpath, f"{self._dtime}_classify.tbl"), format='ipac', overwrite=True)


    def determine_bright(self):
        """
        Checks if the source is bright
        """
        from astropy.io.fits import getdata
        fpma_data = getdata(self._fpma_eventpath)
        fpmb_data = getdata(self._fpmb_eventpath)
        
        times_fpma = []
        fpma_x_center = self._source_pix_coordinates[0][0]
        fpma_y_center = self._source_pix_coordinates[0][1]
        for datum in fpma_data:
            x = float(datum[13])
            y = float(datum[14])
            if np.sqrt((x - fpma_x_center) ** 2 + (y - fpma_y_center)**2) <= 30:
                times_fpma.append(float(datum[0]))
        if len(times_fpma) == 0:
            return 0
        totstart = np.min(times_fpma)
        totstop = np.max(times_fpma)
        lc_bins = np.arange(float(totstart), float(totstop), 1)
        lc, bines = np.histogram(times_fpma, lc_bins)
        t = []
        count = 0
        for idx, l in enumerate(lc):
            t.append(idx)
            if l != 0:
                count += 1
        t = np.array(t)
        if len(lc) == 0:
            return 0
        return max(lc)


    def sky_to_detector(self, ra, dec, mod='A'):
        """
        Calculate the position of an object in detector coordinates from sky 
        coordinates using ``nuskytodet``

        Parameters
        ----------
        ra : float
            The RA of the object

        dec : float
            The dec of the object

        Returns
        """
        cur_path = os.getcwd()
        os.chdir(self._evdir)

        mastaspectfile = f"nu{self.seqid}_mast.fits"
        attfile = f"./../auxil/nu{self.seqid}_att.fits"
        outfile = f"nu{self.seqid}_detpos.fits"

        command_string = f"nuskytodet instrument=FPMA mastaspectfile={mastaspectfile} "
        command_string += f"attfile={attfile} pntra={ra} pntdec={dec} skydetfile={outfile} "
        command_string += f"clobber=yes >/dev/null 2>&1"
        os.system(command_string)
        os.chdir(cur_path)
        return os.path.abspath(os.path.join("event_cl", outfile))


    def convert_detection_coords(self):
        """ 
        Runs the sky to DET1 coordinate conversion routine on an observation.
        """
        detections, flag, filename = self.read_final_detections()
        
        NoneType = type(None)
        if type(detections) == NoneType:
            return
        if len(detections["INDEX"]) == 0:
            with open(os.path.join(self._detpath, filename.replace(".tbl", "_detc.txt")), 'w') as file:
                    file.write(f" ")
            return
        
        ras = detections["RA"]
        decs = detections["DEC"]
        ra_deg = []
        dec_deg = []
        for idx in range(len(detections["INDEX"])):
            detect_position = SkyCoord(f"{ras[idx]} {decs[idx]}", unit=(u.hourangle, u.deg), frame='fk5')
            ra_deg.append(detect_position.ra.deg)
            dec_deg.append(detect_position.dec.deg)
        for idx in range(len(detections["INDEX"][0:5])):
            outfile = self.sky_to_detector(ra_deg[idx], dec_deg[idx])
            print(outfile)
            data = fits.getdata(outfile)
            mean_detx = []
            mean_dety = []
            for datum in data:
                mean_detx.append(float(datum[1]))
                mean_dety.append(float(datum[2]))
            detx = np.mean(mean_detx)
            dety = np.mean(mean_dety)
            if not os.path.isfile(os.path.join(self._detpath, filename.replace(".tbl", "_detc.txt"))):
                with open(os.path.join(self._detpath, filename.replace(".tbl", "_detc.txt")), 'w') as file:
                    file.write(f"{idx} {detx} {dety}\n")
            else:
                with open(os.path.join(self._detpath, filename.replace(".tbl", "_detc.txt")), 'a') as file:
                    file.write(f"{idx} {detx} {dety}\n")

        # Create completion flag
        os.chdir(self._evtpath)
        with open(f"{self._dtime}_det1convert_flag.txt", 'w') as file:
            file.write("PROCESSING COMPLETE")
        os.chdir(self._mainpath)
            

    def display_detections_det(self, savefig=False, display=True):
        """
        Displays an image of the stacked full exposure image for this observation with the 
        source labeled with a circular region and all observations for this specific 
        observation labeled.
        """
        
        # Import current detections
        NoneType = type(None)
        detections, flag, filey = self.read_final_detections(detectiontype="basic")
        
        # If no detections are found, run the default image display method
        if type(detections) == NoneType:
            self.display_image()
            return
        if len(detections["INDEX"]) == 0:
            self.display_image()
        else:
            fig = plt.figure()
            ax = fig.add_subplot(projection=self.wcs)
            ax.set_facecolor('gray')
            #ax = plt.subplot(projection=self.wcs)
            fig.patch.set_facecolor('black')
            im = ax.imshow(self.stacked_data, origin='lower', norm=matplotlib.colors.LogNorm())
            circle = Point(self._source_pix_coordinates[0][0], self._source_pix_coordinates[0][1]).buffer(self.rlimit / 2.46, resolution=1000)
            plot_polygon(circle, ax=ax, add_points=False, color="yellow")
            #object_region = CircleSkyRegion(center=self._source_position, radius=self.rlimit*u.arcsecond)
            #plot_region = object_region.to_pixel(self.wcs)
            #plot_region.plot(ax=ax, color="green")

            # Loop through detections and plot
            
            x_pix = detections['DET1X']
            x_pi = [float(pix) for pix in x_pix]
            y_pix = detections['DET1Y']
            y_pi = [float(pix) for pix in y_pix]
            probs = np.array(detections['PROB'])
            probs_pi = [float(prob) for prob in probs]

            ploty = ax.scatter(x_pi, y_pi, marker="d", c=probs_pi, s=80, linewidths=1, edgecolors= "black", cmap='spring', norm=matplotlib.colors.LogNorm())               
            ax.tick_params(color='white', labelcolor='white')
            ax.tick_params(axis='both', which='major', labelsize=15)
            im.axes.tick_params(color='white', labelcolor='white')
            #ax.get_coords_overlay(self.wcs)
            #cb = plt.colorbar(ploty, orientation='horizontal', shrink=0.5)
            #cb.set_label('Probability', color='white', fontsize=15)
            #cb.ax.xaxis.set_tick_params(color='white')
            #cb.outline.set_edgecolor('white')
            #plt.setp(plt.getp(cb.ax.axes, 'xticklabels'), color='white', fontsize=15)
            plt.grid(True)
            plt.xlabel('Right Ascension', color='white', fontsize=18)
            plt.ylabel('Declination', color='white', fontsize=18)
            plt.xlim(250, 750)
            plt.ylim(250, 750)
            

            if savefig:
                savepath = os.path.relpath(os.path.join(self._outpath, "final_detections.pdf"))
                plt.savefig(savepath, dpi=1000)
            if display:
                plt.show()


    def collect_det1coords(self, dtime):
        """ 
        Returns all DET1 coordinates for a given detection set.
        """
        with open(os.path.join(self._detpath, f"{dtime}_detections_detc.txt")) as file:
            data = file.readlines()
        coords = {}
        for datum in data:
            parts = datum.split()
            detx = float(parts[1])
            dety = float(parts[2])
            coords[parts[0]] = (detx, dety)
        return coords
    
    
    def straycat_removal(self):
        """
        Checks if this Observation was analyzed in any of the StrayCats surveys 
        and flags detections which where located within any identified regions.
        """
        from astropy.io.fits import getdata
        from regions import Regions
        
        # Import current detections
        NoneType = type(None)
        detections, flag, filename = self.read_final_detections(detectiontype="basicdet1")
                
        # If no detections are found, stop program
        if type(detections) == NoneType:
            return
        if len(detections["INDEX"]) == 0:
            return
                
        with open("/Users/jcj/Documents/research/nustar/nuanalysis/straydir.txt") as file:
            data = file.readlines()

        # Make a clone of the table for record keeping and iteration            
        detect_table = Table()
        for key in detections.columns:
            detect_table[key] = detections[key]
        detect_table.write(os.path.join(self._detpath, f"{self._dtime}_straycut_basicdet1.tbl"), format='ipac', overwrite=True)
        
        # Iterate through strayfiles
        for line in tqdm(data):
            strayinfo = line.split()
            if strayinfo[1] == self._seqid:
                print(f"Match found! {strayinfo[1]}")
                
                # Read in region file
                reg_file = os.path.join("/Volumes/data_ssd_2/straycat/", strayinfo[0])
                region_list = Regions.read(reg_file, format="ds9")
                with open(reg_file, 'r') as filey:
                    reg_info = filey.readlines()
                
                # Ensure det1 images are made
                # self.generate_detector_images()
                
                # Read in mask constructables
                construct = []
                flag = False
                for line in reg_info:
                    if flag:
                        if line[0] == '-':
                            construct.append('-')
                        else:
                            construct.append('+')
                    if line == "image\n":
                        flag = True

                # Trivial single case
                if len(region_list) == 1:
                    pos = region_list[0].center.xy
                    radius = region_list[0].radius
                    src_x = pos[0]
                    src_y = pos[1]
                    mask = Point(src_x, src_y).buffer(radius, resolution=1000)
                    
                # Multi region file case
                if len(region_list) > 1:
                    pos = region_list[0].center.xy
                    radius = region_list[0].radius
                    src_x = pos[0]
                    src_y = pos[1]
                    mask = Point(src_x, src_y).buffer(radius, resolution=1000)
                    
                    # Construct the full region
                    for idx, region in enumerate(region_list):
                        if idx == 0:
                            continue 
                        print(type(region))
                        if construct[idx] == '+':
                            pos = region.center.xy
                            radius = region.radius
                            src_x = pos[0]
                            src_y = pos[1]
                            mask_add = Point(src_x, src_y).buffer(radius, resolution=1000)
                            mask = mask.union(mask_add)
                            
                        if construct[idx] == '-':
                            pos = region.center.xy
                            radius = region.radius
                            src_x = pos[0]
                            src_y = pos[1]
                            mask_subtract = Point(src_x, src_y).buffer(radius, resolution=1000)
                            mask = mask.difference(mask_subtract)

                # Optional plotting mechanism
                #fig, ax = plt.subplots()
                #plot_polygon(mask, ax=ax, add_points=False)
                #plt.xlim(0, 300)
                #plt.ylim(0, 300)
                #plt.show()
                
                # Read in all detections
                detpath = os.path.relpath(os.path.join(self._detpath, f"{self._dtime}_straycut_basicdet1.tbl"))
                detections = Table.read(detpath, format='ipac')
                
                # Terminate early if iteration weeds out all detections
                if len(detections["INDEX"]) == 0:
                    print("All detections eliminated!")
                    return
                
                # Initialize new table with empty columns
                detect_table = Table()
                for key in detections.columns:
                    detect_table[key] = []
                    
                # Test
                filtered_detections = {}
                for col in detections.colnames:
                    filtered_detections[col] = []
                
                # Iterate through detections and record those that pass inspection
                count = 0
                for idx in tqdm(range(len(detections["INDEX"]))):
                    detx = detections["DET1X"][idx]
                    dety = detections["DET1Y"][idx]
                    det_point = Point(detx, dety).buffer(1, resolution=1000)
                    
                    if mask.contains(det_point):
                        count += 1
                        continue
                    for col in detections.colnames:
                        filtered_detections[col].append(detections[col][idx])
                    
                detect_table = Table()
                for key in filtered_detections:
                    detect_table[key] = filtered_detections[key]
                
                print(f"Ratio of detections removed: {count/len(detections['INDEX'])}")
                detect_table.write(os.path.join(self._detpath, f"{self._dtime}_straycut_basicdet1.tbl"), format='ipac', overwrite=True)
                
        # Create completion flag
        os.chdir(self._evtpath)
        with open(f"{self._dtime}_straycutbasic_flag.txt", 'w') as file:
            file.write("PROCESSING COMPLETE")
        os.chdir(self._mainpath)
        return
    
    
    def convert_to_det1(self):
        """ 
        Evaluates the DET1 coordinates for recalculated poisson detections.
        """
        # Display terminal
        print('#' * 90)
        print(f"Converting Poisson Detections to Det1 for SEQID: {self._seqid}")
        print(f"Time scale: {self._dtime} seconds.")
        print(f"Source Position: {self._source_pix_coordinates[0][0]}, {self._source_pix_coordinates[0][1]}")
        print('#' * 90)
        
        detections, flag, filename = self.read_final_detections(detectiontype="poisson")
        
        NoneType = type(None)
        if type(detections) == NoneType:
            return
        if len(detections["INDEX"]) == 0:
            # Create completion flag
            os.chdir(self._evtpath)
            with open(f"{self._dtime}_det1poisson_flag.txt", 'w') as file:
                file.write("PROCESSING COMPLETE")
            os.chdir(self._mainpath)
            return
        
        detx_list = []
        dety_list = []
        for idx in tqdm(range(len(detections["INDEX"]))):
            ra = detections["RA"][idx]
            dec = detections["DEC"][idx]
            detect_position = SkyCoord(f"{ra} {dec}", unit=(u.hourangle, u.deg), frame='fk5')
            ra = detect_position.ra.deg
            dec = detect_position.dec.deg
            outfile = self.sky_to_detector(ra, dec)
            data = fits.getdata(outfile)
            mean_detx = []
            mean_dety = []
            for datum in data:
                mean_detx.append(float(datum[1]))
                mean_dety.append(float(datum[2]))
            detx = np.mean(mean_detx)
            dety = np.mean(mean_dety)
            detx_list.append(detx)
            dety_list.append(dety)
            os.remove(outfile)
        detections["DET1X"] = detx_list
        detections["DET1Y"] = dety_list
        detections.write(os.path.join(self._detpath, f"{self._dtime}_poisson_det1.tbl"), format='ipac', overwrite=True)
        
        # Create completion flag
        os.chdir(self._evtpath)
        with open(f"{self._dtime}_det1poisson_flag.txt", 'w') as file:
            file.write("PROCESSING COMPLETE")
        os.chdir(self._mainpath)
        return
    
    
    def convert_basic_to_det1(self):
        """ 
        Calculates and assigns DET1 coordinates to detection flagged pixels.
        """
        
        # Display terminal
        print('#' * 90)
        print(f"Converting Basic Detections to Det1 for SEQID: {self._seqid}")
        print(f"Time scale: {self._dtime} seconds.")
        print(f"Source Position: {self._source_pix_coordinates[0][0]}, {self._source_pix_coordinates[0][1]}")
        print('#' * 90)
        
        detections, flag, filename = self.read_final_detections(detectiontype="basic")
        
        NoneType = type(None)
        if type(detections) == NoneType:
            # Create completion flag
            os.chdir(self._evtpath)
            with open(f"{self._dtime}_det1basic_flag.txt", 'w') as file:
                file.write("PROCESSING COMPLETE")
            os.chdir(self._mainpath)
            return
        if len(detections["INDEX"]) == 0:
            # Create completion flag
            os.chdir(self._evtpath)
            with open(f"{self._dtime}_det1basic_flag.txt", 'w') as file:
                file.write("PROCESSING COMPLETE")
            os.chdir(self._mainpath)
            return
        
        detx_list = []
        dety_list = []
        for idx in tqdm(range(len(detections["INDEX"]))):
            #if len(detections['INDEX']) > 800:
            #    return
            ra = detections["RA"][idx]
            dec = detections["DEC"][idx]
            detect_position = SkyCoord(f"{ra} {dec}", unit=(u.hourangle, u.deg), frame='fk5')
            ra = detect_position.ra.deg
            dec = detect_position.dec.deg
            outfile = self.sky_to_detector(ra, dec)
            data = fits.getdata(outfile)
            mean_detx = []
            mean_dety = []
            for datum in data:
                mean_detx.append(float(datum[1]))
                mean_dety.append(float(datum[2]))
            detx = np.mean(mean_detx)
            dety = np.mean(mean_dety)
            print(detx)
            detx_list.append(detx)
            dety_list.append(dety)
            os.remove(outfile)
        detections["DET1X"] = detx_list
        detections["DET1Y"] = dety_list
        detections.write(os.path.join(self._detpath, f"{self._dtime}_basic_det1.tbl"), format='ipac', overwrite=True)
        
        # Create completion flag
        os.chdir(self._evtpath)
        with open(f"{self._dtime}_det1basic_flag.txt", 'w') as file:
            file.write("PROCESSING COMPLETE")
        os.chdir(self._mainpath)
        return


    def main_source_removal(self):
        """
        This method loops through a set of detections and eliminates 
        those which are located within the main source mask for this 
        observation.
        """
        
        # Import current detections
        NoneType = type(None)
        detections, flag, filename = self.read_final_detections(detectiontype="basicdet1")
        
        # If no detections are found, stop program
        if type(detections) == NoneType:
            return
        if len(detections["INDEX"]) == 0:
            return
        
        # Sanity check for pixels not in the main source mask
        valid_indices = []
        for idx, val in tqdm(enumerate(detections['XPIX'])):
            pix_x = detections['XPIX'][idx]
            pix_y = detections['YPIX'][idx]
            if np.sqrt((float(pix_x) - float(self._source_pix_coordinates[0][0]))**2 + (float(pix_y) - float(self._source_pix_coordinates[0][1]))**2) > (180) / 2.5:
                valid_indices.append(idx)    
        
        # Recreate dictionary for remaining detections
        filtered_detections = {}
        for col in detections.colnames:
            filtered_detections[col] = []
        for col in detections.colnames:
            for idx in valid_indices:
                filtered_detections[col].append(detections[col][idx])
        
        # Write newly filtered dictionary to astropy table
        detect_table = Table()
        for key in filtered_detections:
            detect_table[key] = filtered_detections[key]
        finalpath = os.path.relpath(os.path.join(self._detpath, f"{self._dtime}_basicdet_main3.tbl"))
        detect_table.write(finalpath, format='ipac', overwrite=True)
        
        # Create completion flag
        os.chdir(self._evtpath)
        with open(f"{self._dtime}_mainsourceremoval_flag.txt", 'w') as file:
            file.write("PROCESSING COMPLETE")
        os.chdir(self._mainpath)
        return
    
    
    def classify_persistent_transients(self):
        """
        This method applies a simple algorithm to classify persistent vs 
        transient sources.
        """
        
        print("#" * 90)
        print(f"Measuring the persistence value for sources.")
        print("#" * 90)
        detections, flag, filey = self.read_final_detections(detectiontype="basicdet1")
        
        # Check if det1
        test_flag = os.path.join(self._evdir , f"{self._dtime}_det1basic_flag.txt")
        if not os.path.isfile(test_flag):
            return
        
        # Check if atleast one detection is present
        NoneType = type(None)
        if type(detections) == NoneType:
            # Create completion flag
            #os.chdir(self._evtpath)
            #with open(f"{self._dtime}_classify_det1_flag.txt", 'w') as file:
            #    file.write("PROCESSING COMPLETE")
            #os.chdir(self._mainpath)
            os.chdir(self._evtpath)
            with open(f"{self._dtime}_classify_det1_flag.txt", 'w') as file:
                file.write("PROCESSING COMPLETE")
            os.chdir(self._mainpath)
            print("No det1 detections found")
            return
        if len(detections["INDEX"]) == 0:
            # Create completion flag
            os.chdir(self._evtpath)
            with open(f"{self._dtime}_classify_det1_flag.txt", 'w') as file:
                file.write("PROCESSING COMPLETE")
            os.chdir(self._mainpath)
            print("No det1 detections were confirmed")
            return
        
        tpersist = []
        cpersist = []
        
        for index, value in tqdm(enumerate(detections['INDEX'])):
            tstarts = []
            tstops = []
            channels = []
            xpix = float(detections["XPIX"][index])
            ypix = float(detections["YPIX"][index])
            tstarts.append(detections["TSTART"][index])
            tstops.append(detections["TSTOP"][index])
            channels.append(detections["BOUND"][index])
            for index2, value2 in enumerate(detections['INDEX']):
                if np.sqrt((float(detections["XPIX"][index2]) - xpix) ** 2 + (float(detections["YPIX"][index2]) - ypix) ** 2) <= 20:
                    tstarts.append(detections["TSTART"][index2])
                    tstops.append(detections["TSTOP"][index2])
                    channels.append(detections["BOUND"][index2])
            present_percent = min(len(set(tstarts)), len(set(tstops))) / (self._start_intervals_pass1 + self._start_intervals_pass2)
            channels_percent = len(set(channels)) / 4
            tpersist.append(present_percent)
            cpersist.append(channels_percent)            
        detections['TPERSIST'] = tpersist
        detections['CPERSIST'] = cpersist
        detect_table = Table(detections)
        detect_table.write(os.path.join(self._detpath, f"{self._dtime}_basicdet1_classify.tbl"), format='ipac', overwrite=True)
        os.chdir(self._evtpath)
        with open(f"{self._dtime}_classify_det1_flag.txt", 'w') as file:
            file.write("PROCESSING COMPLETE")
        os.chdir(self._mainpath)