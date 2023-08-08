import os
import subprocess
import glob
import random
import numpy as np
import matplotlib.pyplot as plt
from nustar import *
from helpers import *
from astropy.io import fits
from astropy import units as u
from astropy.wcs import WCS
from astropy.table import QTable, Table, Column
from tqdm import tqdm
import radial_profile
from astropy.coordinates import SkyCoord
from regions import CircleSkyRegion


class NuAnalysis(Observation):
    """
    This class defines an object to be used for performing analysis on a NuSTAR observation.
    """

    def __init__(self, dtime, snr_threshold, path=False, seqid=False, evdir=False, out_path=False, clean=False, bifrost=False, object_name=None, sessionid=None):
        # Define analysis parameters
        self._snr = snr_threshold
        self._dtime = dtime
        self._clean = clean
        self._phi_bounds = self.read_in_phi_bounds("nustar_pilow.txt", "nustar_pihi.txt")
        self._refpath = path
        self._refoutpath = out_path
        self._sessionid = sessionid
        
        # Modifies for server performance and fidelity
        if not bifrost:
            self._contents = os.listdir(path)
        if bifrost:
            if not os.path.isdir(path.replace(seqid, '')):
                os.mkdir(path.replace(seqid, ''))
            if not os.path.isdir(path):
                os.mkdir(path)
            if not os.path.isdir(evdir):
                os.mkdir(evdir)
            self._contents = os.listdir(path)
        generate_directory(self._refoutpath, overwrite=False)
        
        if not clean:
            if bifrost:
                self._object = object_name
                if self._object == None:
                    raise AssertionError("OBJECT NAME NOT PROVIDED")
                generate_directory(evdir, overwrite=True)
                generate_directory(evdir + "event_cl/", overwrite=True)
                datapath = f"../../../../../nustar/fltops/{self._object}/{seqid}"
                subprocess.run(["nupipeline", datapath, f"nu{seqid}", "./event_cl/", "saamode=STRICT", "tentacle=yes", "clobber=yes"], cwd=path)
                with open(path + "event_cl/processing_flag.txt", "w") as file:
                    file.write("PROCESSING COMPLETE")
                self._contents = os.listdir(path)
                self._clean = True
            else:
                #generate_directory(evdir, overwrite=True)
                subprocess.run(["nupipeline", path, f"nu{seqid}", evdir, "saamode=STRICT", "tentacle=yes", "clobber=yes"])
                self._clean = True
        
        super().__init__(path, seqid, evdir, out_path)
        
        if not self._clean:
            self.run_cleaning_script()

        #if "science.fits" not in self._contents:
        #    print("Generating Science FITS image")
        #    print(f"FILES: {self.science_files['A'][0]}, {self.science_files['B'][0]}")
        #    infiles = f"{self.science_files['A'][0].replace(self._refpath, '')} , {self.science_files['B'][0].replace(self._refpath, '')}"
        #    outfile = self._refpath + "science.fits"
        #    name, number = make_xselect_commands(infiles, outfile, self._refpath, 1.6, 79, sessionid + 1000, evt_extract=True)
        #    subprocess.run(["xselect", f"@{number}xsel.xco"])
        #    subprocess.run(["rm", f"{number}xsel.xco"])

        #hdu = fits.open(self._refpath + "science.fits", uint=True)[0]
        #self.wcs = WCS(hdu.header)
        #self.data = hdu.data
        #self.coords = (hdu.header["RA_OBJ"], hdu.header["DEC_OBJ"])
        #coordinates = radial_profile.find_source(self._refpath + "science.fits", show_image=False, filt_range=3)
        
        #rind, rad_profile, radial_err, psf_profile = radial_profile.make_radial_profile(self._refpath + "science.fits", show_image=False,
        #                                                         coordinates = coordinates)
        #self.rlimit = radial_profile.optimize_radius_snr(rind, rad_profile, radial_err, psf_profile, show=False)
        #print(self.rlimit)
        self._time_bins = self.generate_timebins()
        #self._detections = None

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
    def phi_bounds(self):
        """
        Returns the phi_bounds set for analysis.
        """
        return self._phi_bounds


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
        
    # Methods begin below

    def display_image(self):
        ax = plt.subplot(projection=self.wcs)
        im = ax.imshow(self.data, origin='lower')
        center_sky = SkyCoord(self.coords[0], self.coords[1], unit='deg', frame='fk5')
        region_sky = CircleSkyRegion(center=center_sky, radius=50*u.arcsecond)
        pixel_region = region_sky.to_pixel(self.wcs)
        pixel_region.plot(ax=ax, color="green")


        ax.get_coords_overlay()

        plt.colorbar(im)
        plt.grid(True)
        plt.xlabel('')
        plt.ylabel('')
        plt.show()

    
    def display_detections(self):
        ax = plt.subplot(projection=self.wcs)
        im = ax.imshow(self.data, origin='lower')
        center_sky = SkyCoord(self.coords[0], self.coords[1], unit='deg', frame='fk5')
        region_sky = CircleSkyRegion(center=center_sky, radius=50*u.arcsecond)
        pixel_region = region_sky.to_pixel(self.wcs)
        pixel_region.plot(ax=ax, color="yellow")

        for bound in self.phi_bounds:
            detections = self.detection_dir_processing(bound)
            if detections == None:
                continue
            if detections != None:
                for idx in range(len(detections["RA"])):
                    ra, dec = self.ra_dec_todeg(detections["RA"][idx], detections["DEC"][idx])
                    detect_coord = SkyCoord(ra, dec, unit='deg', frame='fk5')
                    detect_circle = CircleSkyRegion(center=detect_coord, radius=5*u.arcsecond)
                    detect_region = detect_circle.to_pixel(self.wcs)
                    detect_region.plot(ax=ax, color="white")                
        
        ax.get_coords_overlay(self.wcs)
        plt.colorbar(im)
        plt.grid(True)
        plt.xlabel('')
        plt.ylabel('')
        plt.show()


    def read_in_phi_bounds(self, pilow_file, pihi_file):
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

        phi_bounds = []
        low_phi = open(f"./ref_files/{pilow_file}", 'r')
        high_phi = open(f"./ref_files/{pihi_file}", 'r')
        low_phis = low_phi.readlines()
        high_phis = high_phi.readlines()
        
        if len(low_phis) != len(high_phis):
            raise IndexError("PHI bounds must have the same length!")
        
        for idx in range(len(low_phis)):
            phi_bounds.append((low_phis[idx].replace('\n', ''), high_phis[idx].replace('\n', '')))
        
        return phi_bounds


    def run_cleaning_script(self):
        """
        This script will run `nupipeline` on data that either lacks clean event files, or is determined 
        to require reprocessing.
        """
        
        generate_directory(self._evdir, overwrite=True)
        #subprocess.run(["nupipeline", self._refpath, f"nu{self._seqid}", self._evdir])


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

        # generate intervals for PASS 2
        intervals_p2 = []
        intervals_p2.append(self._event_times[0])
        cur_int = self._event_times[0] + self._dtime / 2
        max_int = np.max(self._event_times)
        while cur_int < max_int:
            intervals_p2.append(cur_int)
            cur_int += self._dtime 
        intervals_p2.append(max_int)

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


    def event_extraction(self):
        """
        The event extraction and temporal binning process for finding XRTs.
        """
        if "detections" not in self._contents:
            os.mkdir(self._refpath + "detections/")

        for bound in self._phi_bounds:
            detect_path = self._refpath + f"detections/{bound[0]}-{bound[1]}_{self._dtime}-{self._snr}/"
            generate_directory(detect_path, overwrite=False)
            print(f"Binning data: PHI Bound {bound[0]}-{bound[1]}; Cycle 1")
            for interval in tqdm(self._time_bins[0]):
                if len(self._time_bins[0][interval][2]) == 0:
                    continue
                else:
                    with open(self._refpath + "/event_cl/xselect.xco", 'w') as script:
                        script.write(f'{self._sessionid}\n')
                        script.write(f"read events\n")
                        script.write(".\n")
                        script.write(f"nu{self._seqid}A01_cl.evt , nu{self._seqid}B01_cl.evt\n")
                        script.write('yes\n')
                        script.write("filter PHA_CUTOFF\n")
                        script.write(f"{bound[0]}\n")
                        script.write(f"{bound[1]}\n")
                        script.write("filter time scc\n")
                        script.write(f"{self._time_bins[0][interval][0]} , {self._time_bins[0][interval][1]}\n")
                        script.write("x\n")
                        script.write("extract events\n")
                        script.write("\n")
                        script.write(f"save events ./../detections/{bound[0]}-{bound[1]}_{self._dtime}-{self._snr}/nu_{self._time_bins[0][interval][0]}-{self._time_bins[0][interval][1]}.evt\n")
                        script.write('no\n')
                        script.write("exit no\n")
                    subprocess.run(["xselect", "@xselect.xco"], cwd=self._evdir, capture_output=True)
                    subprocess.run(["rm", "xselect.xco"], cwd=self._evdir, capture_output=True)
                    with open(self._evdir + f"{self._dtime}_binning_flag.txt", "w") as file:
                        file.write("PROCESSING COMPLETE")
            print(f"Binning data: PHI Bound {bound[0]}-{bound[1]}; Cycle 2")
            for interval in tqdm(self._time_bins[1]):
                if len(self._time_bins[1][interval][2]) == 0:
                    continue
                else:
                    with open(self._refpath + "/event_cl/xselect.xco", 'w') as script:
                        script.write(f'{self._sessionid}\n')
                        script.write(f"read events\n")
                        script.write(".\n")
                        script.write(f"nu{self._seqid}A01_cl.evt , nu{self._seqid}B01_cl.evt\n")
                        script.write('yes\n')
                        script.write("filter PHA_CUTOFF\n")
                        script.write(f"{bound[0]}\n")
                        script.write(f"{bound[1]}\n")
                        script.write("filter time scc\n")
                        script.write(f"{self._time_bins[1][interval][0]} , {self._time_bins[1][interval][1]}\n")
                        script.write("x\n")
                        script.write("extract events\n")
                        script.write("\n")
                        script.write(f"save events ./../detections/{bound[0]}-{bound[1]}_{self._dtime}-{self._snr}/nu_{self._time_bins[1][interval][0]}-{self._time_bins[1][interval][1]}.evt\n")
                        script.write('no\n')
                        script.write("exit no\n")
                    subprocess.run(["xselect", "@xselect.xco"], cwd=self._evdir, capture_output=True)
                    subprocess.run(["rm", "xselect.xco"], cwd=self._evdir, capture_output=True)
                    with open(self._evdir + f"{self._dtime}_binning_flag.txt", "w") as file:
                        file.write("PROCESSING COMPLETE")


    def sliding_cell_detection(self):
        """
        The sliding cell detection call on the temporal binned data.
        """

        print("Performing sliding cell source detection on stacked images.")
        for bound in self._phi_bounds:
            print(f"Processing PHI Channels: {bound[0]}-{bound[1]}")
            running_directory = self._refpath + f"detections/{bound[0]}-{bound[1]}_{self._dtime}-{self._snr}/"
            stacked_images = glob.glob(running_directory + "*.evt")
            stacked_files = [file.replace(running_directory, '') for file in stacked_images]
            for file in tqdm(stacked_files):
                if len(fits.getdata(running_directory + file)) != 0:
                    with open(self._refpath + f"detections/{bound[0]}-{bound[1]}_{self._dtime}-{self._snr}/ximage.xco", 'w') as script:
                        script.write(f"read/fits/size=800/{file}\n")
                        script.write(f"detect/snr={self._snr}/source_box_size=32/filedet={file.replace('.evt', '')}.det/fitsdet={file.replace('.evt', '')}.fits\n")
                        script.write("exit")
                    subprocess.run(["ximage", "@ximage.xco"], cwd=running_directory, capture_output=True)
                    subprocess.run(["rm", "ximage.xco"], cwd=running_directory, capture_output=True)
                    with open(self._evdir + f"{self._dtime}_flag.txt", 'w') as file:
                        file.write("DONE")

    
    def detection_merging(self):
        """
        Performs the first round of detection redundancy culling on detections caught by the sliding cell algorithm
        """
        for bound in self.phi_bounds:
            det_files = glob.glob(self._refpath + f"detections/{bound[0]}-{bound[1]}_{self._dtime}-{self._snr}/*.det")
            file_strings = []
            for file in det_files:
                file_strings.append(file.replace(self._refpath + f"detections/{bound[0]}-{bound[1]}_{self._dtime}-{self._snr}/", ''))
            script_string = "srcmrg/out=mrg.txt/tolerance=5e1"
            for file in file_strings:
                script_string += f" {file}"
            script_string += "\n"
            if len(file_strings) == 0:
                continue
                #print("No detections found!")
            else:
                with open(self._refpath + f"detections/{bound[0]}-{bound[1]}_{self._dtime}-{self._snr}/src_merge.xco", 'w') as script:
                    script.write(script_string)
                    script.write("exit\n")
                subprocess.run(["ximage", "@src_merge.xco"], cwd=self._refpath + f"detections/{bound[0]}-{bound[1]}_{self._dtime}-{self._snr}", capture_output=True)
                subprocess.run(["rm", "src_merge.xco"], cwd=self._refpath + f"detections/{bound[0]}-{bound[1]}_{self._dtime}-{self._snr}", capture_output=True)


    def run_detection_pipeline(self):
        self.event_extraction()
        self.sliding_cell_detection()
        self.detection_merging()


    def test2_gen(self):
        """
        This method generates a series of detection searches after first stacking and restricting the PI channels 
        of the clean event files using `xselect`, afterwhich it runs a sliding cell detection algorithm using `ximage` 
        The detections are calculated for each time bin defined in the `generate_timebins` method and are then aggregated.
        The directories underwhich the detection information is stored has the following syntax:

        dir = "./seqid/detections/phi_low-phi_high-dtime-snr/"
        """

        if "detections" not in self._contents:
            os.mkdir(self._refpath + "detections/")
        
        # Begin by selecting the data which lies in the appropriate energy range.
        
        for bound in self._phi_bounds:
            detect_path = self._refpath + f"detections/{bound[0]}-{bound[1]}_{self._dtime}-{self._snr}/"
            generate_directory(detect_path, overwrite=False)
            print(f"Binning data: PHI Bound {bound[0]}-{bound[1]}; Cycle 1")
            for interval in tqdm(self._time_bins[0]):
                if len(self._time_bins[0][interval][2]) == 0:
                    continue
                else:
                    with open(self._refpath + "/event_cl/xselect.xco", 'w') as script:
                        n1 = random.random()
                        n2 = random.random()
                        n3 = random.random()
                        n4 = random.choice([100, 1000, 10000])
                        number = n1*n2*n3*n4
                        script.write(f'{number}\n')
                        script.write('\n')
                        script.write(f"read events\n")
                        script.write(".\n")
                        script.write(f"nu{self._seqid}A01_cl.evt , nu{self._seqid}B01_cl.evt\n")
                        script.write('yes\n')
                        script.write("filter PHA_CUTOFF\n")
                        script.write(f"{bound[0]}\n")
                        script.write(f"{bound[1]}\n")
                        script.write("filter time scc\n")
                        script.write(f"{self._time_bins[0][interval][0]} , {self._time_bins[0][interval][1]}\n")
                        script.write("x\n")
                        script.write("extract events\n")
                        script.write("\n")
                        script.write(f"save events ./../detections/{bound[0]}-{bound[1]}_{self._dtime}-{self._snr}/nu_{self._time_bins[0][interval][0]}-{self._time_bins[0][interval][1]}.evt\n")
                        script.write('no\n')
                        script.write("exit no\n")
                    subprocess.run(["xselect", "@xselect.xco"], cwd=self._evdir, capture_output=True)
                    subprocess.run(["rm", "xselect.xco"], cwd=self._evdir, capture_output=True)
            print(f"Binning data: PHI Bound {bound[0]}-{bound[1]}; Cycle 2")
            for interval in tqdm(self._time_bins[1]):
                if len(self._time_bins[1][interval][2]) == 0:
                    continue
                else:
                    with open(self._refpath + "/event_cl/xselect.xco", 'w') as script:
                        n1 = random.random()
                        n2 = random.random()
                        n3 = random.random()
                        n4 = random.choice([100, 1000, 10000])
                        number = n1*n2*n3*n4
                        script.write(f'{number}\n')
                        script.write('\n')
                        script.write(f"read events\n")
                        script.write(".\n")
                        script.write(f"nu{self._seqid}A01_cl.evt , nu{self._seqid}B01_cl.evt\n")
                        script.write('yes\n')
                        script.write("filter PHA_CUTOFF\n")
                        script.write(f"{bound[0]}\n")
                        script.write(f"{bound[1]}\n")
                        script.write("filter time scc\n")
                        script.write(f"{self._time_bins[1][interval][0]} , {self._time_bins[1][interval][1]}\n")
                        script.write("x\n")
                        script.write("extract events\n")
                        script.write("\n")
                        script.write(f"save events ./../detections/{bound[0]}-{bound[1]}_{self._dtime}-{self._snr}/nu_{self._time_bins[1][interval][0]}-{self._time_bins[1][interval][1]}.evt\n")
                        script.write('no\n')
                        script.write("exit no\n")
                    subprocess.run(["xselect", "@xselect.xco"], cwd=self._evdir, capture_output=True)
                    subprocess.run(["rm", "xselect.xco"], cwd=self._evdir, capture_output=True)
        
        print("Performing sliding cell source detection on stacked images.")
        for bound in self._phi_bounds:
            print(f"Processing PHI Channels: {bound[0]}-{bound[1]}")
            running_directory = self._refpath + f"detections/{bound[0]}-{bound[1]}_{self._dtime}-{self._snr}/"
            stacked_images = glob.glob(running_directory + "*.evt")
            stacked_files = [file.replace(running_directory, '') for file in stacked_images]
            for file in tqdm(stacked_files):
                if len(fits.getdata(running_directory + file)) != 0:
                    with open(self._refpath + f"detections/{bound[0]}-{bound[1]}_{self._dtime}-{self._snr}/ximage.xco", 'w') as script:
                        script.write(f"read/fits/size=800/{file}\n")
                        script.write(f"detect/snr={self._snr}/filedet={file.replace('.evt', '')}.det/fitsdet={file.replace('.evt', '')}.fits\n")
                        script.write("exit")
                    subprocess.run(["ximage", "@ximage.xco"], cwd=running_directory, capture_output=True)
                    subprocess.run(["rm", "ximage.xco"], cwd=running_directory, capture_output=True)
                
            
        print("Merging together results to find unique detections.")
        for bound in self.phi_bounds:
            det_files = glob.glob(self._refpath + f"detections/{bound[0]}-{bound[1]}_{self._dtime}-{self._snr}/*.det")
            file_strings = []
            for file in det_files:
                file_strings.append(file.replace(self._refpath + f"detections/{bound[0]}-{bound[1]}_{self._dtime}-{self._snr}/", ''))
            script_string = "srcmrg/out=mrg.txt/tolerance=5e1"
            for file in file_strings:
                script_string += f" {file}"
            script_string += "\n"
            if len(file_strings) == 0:
                print("No detections found!")
            else:
                with open(self._refpath + f"detections/{bound[0]}-{bound[1]}_{self._dtime}-{self._snr}/src_merge.xco", 'w') as script:
                    script.write(script_string)
                    script.write("exit\n")
                subprocess.run(["ximage", "@src_merge.xco"], cwd=self._refpath + f"detections/{bound[0]}-{bound[1]}_{self._dtime}-{self._snr}")
                subprocess.run(["rm", "src_merge.xco"], cwd=self._refpath + f"detections/{bound[0]}-{bound[1]}_{self._dtime}-{self._snr}")
        

    def detection_summary(self):
        for bound in self.phi_bounds:
            detections = self.detection_dir_processing(bound)
            if detections == None:
                n_det = 0
            if detections != None:
                n_det = len(detections["INDEX"]) 
            print(f"PI Channels: {bound[0]}-{bound[1]} -- {n_det} detections found.")



    def process_detections(self):
        for bound in self.phi_bounds:
            detections = self.detection_dir_processing(bound)
            if detections == None:
                n_det = 0
            if detections != None:
                n_det = len(detections["INDEX"]) 
                self.nuproducts(detections, bound)
            print(f"PI Channels: {bound[0]}-{bound[1]} -- {n_det} detections found.")


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
            a dictionary containing the following keys: `INDEX`, `COUNTS`, `XPIX`, `YPIX`, `VIGCOR`, 
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

        RAISES
        ------
        None
        """

        if not os.path.isdir(self._refpath + f"detections/{bounds[0]}-{bounds[1]}_{self._dtime}-{self._snr}/"):
            return None
        if "mrg.txt" not in os.listdir(self._refpath + f"detections/{bounds[0]}-{bounds[1]}_{self._dtime}-{self._snr}/"):
            return None           

        with open(self._refpath + f"detections/{bounds[0]}-{bounds[1]}_{self._dtime}-{self._snr}/mrg.txt") as detections:
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
        det_files = glob.glob(self._refpath + f"detections/{bounds[0]}-{bounds[1]}_{self._dtime}-{self._snr}/*.det")
        detect_info = {}
        #print(f"READING ALL DETECTIONS: {bounds}")
        #print(len(det_files))
        for file in det_files:
            with open(file) as detections:
                for i in range(14):
                    detections.readline()
                times = tuple(file.replace(self._refpath + f"detections/{bounds[0]}-{bounds[1]}_{self._dtime}-{self._snr}/", '').replace(".det", '').split(sep="-"))
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

        for bounds in self.phi_bounds:
            print(bounds)
            self.nuproducts(self.detection_dir_processing(bounds), bounds)

    
    def remove_main_source(self, detect_info):
        trimmed_detect_info = {}        
        for key in detect_info:
            trimmed_detect_info[key] = []
        
        for i in detect_info["INDEX"]:
            ra = detect_info["RA"][int(i) - 1]
            dec = detect_info["DEC"][int(i) - 1]
            c = SkyCoord(f"{ra} {dec}", unit=(u.hourangle, u.deg))
            if self._source_position.separation(c).arcsec > 50:
                for key in detect_info:
                    trimmed_detect_info[key].append(detect_info[key][int(i) - 1])
        return trimmed_detect_info
    
    
    def detection_dir_processing(self, bounds):
        
        # Read in unique detection information
        unique_detect_info = self.read_unique_detections(bounds)
        if unique_detect_info == None:
            return None

        # Read in all detections for a given PI channel range
        all_detect_info = self.read_detection_dir(bounds)
        n_all_det = 0
        for time in all_detect_info:
            n_all_det += len(all_detect_info[time]["INDEX"])

        # Begin by eliminating the main source
        trimmed_detect_info = self.remove_main_source(unique_detect_info)
        
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
        
        #for time in all_detect_info:
        #    print(len(all_detect_info[time]["INDEX"]))

        
        self._detections = trimmed_all_info
        return trimmed_all_info
    

    def write_net_detections(self):
        for bound in self._phi_bounds:
            trimmed_all_info = self.detection_dir_processing(bound)
            if trimmed_all_info != None:

                # Applying table corrections.
                n_obj = len(trimmed_all_info["INDEX"])
                trimmed_all_info["INDEX"] = list(range(1, n_obj + 1))
                t_starts = [val[0].replace("nu_", '') for val in trimmed_all_info["TIMES"]]
                t_stops = [val[1] for val in trimmed_all_info["TIMES"]]
                trimmed_all_info["TSTART"] = t_starts
                trimmed_all_info["TSTOP"] = t_stops
                count_vals = [counts.split('+/-')[0] for counts in trimmed_all_info["COUNTS"]]
                count_err_vals = [counts.split('+/-')[1] for counts in trimmed_all_info["COUNTS"]]
                trimmed_all_info["COUNTS"] = count_vals
                trimmed_all_info["COUNTSERR"] = count_err_vals
                trimmed_all_info["SEQID"] = [self._seqid for i in range(n_obj)]

                del trimmed_all_info["TIMES"]
                # Construct and write table
                detect_table = Table()
                for key in trimmed_all_info:
                    detect_table[key] = trimmed_all_info[key]
                detect_table.write(self._refpath + f"detections/{self._dtime}.tbl", format='ipac', overwrite=True)



    def nuproducts(self, detect_info, pi_bounds):
        print('here')
        if detect_info == None:
            print("Oh no! No detections!")
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
        all_files = []
        for bound in self._phi_bounds:
            files = glob.glob(self._refpath + f"products/{bound[0]}_{bound[1]}-{self._dtime}_{self._snr}/**/*.flc", recursive=True)
            for file in files:
                all_files.append(file)
        return all_files


    def acquire_lightcurve(self, fits_file):
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

        lc_data = {}
        lc_data["TIME"] = []
        lc_data["RATE"] = []
        lc_data["RATEERR"] = []
        with fits.open(fits_file) as hdul:
            fits_data = hdul[1].data
        for datum in fits_data:
            lc_data["TIME"].append(datum[0])
            lc_data["RATE"].append(datum[2])
            lc_data["RATEERR"].append(datum[3])
        return lc_data


    def acquire_pha(self, pilow, pihi, n):
        """
        """
        pha_file = self._refpath + f"products/{pilow}_{pihi}/{n}/nu{self._seqid}B01_sr.pha"
        with fits.open(pha_file) as hdul:
            pha_data = hdul[1].data
        return pha_data


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


    def plot_pha(self, pilow, pihi, n):
        """
        """

        pha_data = self.acquire_pha(pilow, pihi, n)
        fig, ax = plt.subplots()
        y_dat = []
        x_dat = []
        first_key = list(self._exposure.keys())[0]
        for datum in pha_data:
            x_dat.append(chan_to_energy(datum[0]))
            y_dat.append(datum[1]/self._exposure[first_key])
        ax.plot(x_dat, y_dat)
        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Counts/Time (c/s)")
        plt.show()
    

    def generate_region(self, path, xpix, ypix):
        """
        This method generates a region file to be used for data analysis with a fits image.
        """
        pass

    

