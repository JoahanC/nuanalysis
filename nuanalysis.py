import os
import subprocess
import glob
import numpy as np
import matplotlib.pyplot as plt
from nustar import *
from helpers import *
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS


class NuAnalysis(Observation):
    """
    This class defines an object to be used for performing analysis on a NuSTAR observation.
    """

    def __init__(self, dtime, snr_threshold, path=False, seqid=False, evdir=False, out_path=False, clean=False):
        self._snr = snr_threshold
        self._dtime = dtime
        self._clean = clean
        self._phi_bounds = self.read_in_phi_bounds("nustar_pilow.txt", "nustar_pihi.txt")
        self._refpath = path
        self._contents = os.listdir(path)
        generate_directory(out_path, overwrite=True)
        if "event_cl" not in self._contents or not clean:
            generate_directory(evdir)
            subprocess.run(["nupipeline", path, f"nu{seqid}", evdir])
            self._clean = True
        
        super().__init__(path, seqid, evdir, out_path)

        if not self._clean:
            self.run_cleaning_script()

        self._time_bins = self.generate_timebins()
        self._detections = None

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
        subprocess.run(["nupipeline", self._refpath, f"nu{self._seqid}", self._evdir])


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
                data_intervals_2[idx][2].append(datapoint)
                break
            if datapoint > intervals_p2[idx] and datapoint <= intervals_p2[idx + 1]:
                data_intervals_2[idx][2].append(datapoint)
            else:
                idx += 1
                data = []
                data_intervals_2[idx] = [intervals_p2[idx], intervals_p2[idx + 1], data]
                data_intervals_2[idx][2].append(datapoint)

        return [data_intervals_1, data_intervals_2]


    def generate_detections(self):
        """
        This method generates a series of detection searches after first stacking and restricting the PI channels 
        of the clean event files using `xselect`, afterwhich it runs a sliding cell detection algorithm using `ximage` 
        The detections are calculated for each time bin defined in the `generate_timebins` method and are then aggregated.
        The directories underwhich the detection information is stored has the following syntax:

        dir = "./obsid/detections/phi_low-phi_high-snr-t_low-t_high"

        Arguments
        ---------
        None

        Returns
        -------
        None

        Raises
        ------
        None
        """

        if "detections" not in self._contents:
            os.mkdir(self._refpath + "/detections/")
        detect_path = self._refpath + f"detections/{self._dtime}-{self._snr}/"
        if f"{self._dtime}-{self._snr}" in os.listdir(self._refpath + "detections/"):
            self._detections = self.read_detections()
        else:
            generate_directory(detect_path, overwrite=True)
            # Begin by selecting the data which lies in the appropriate energy range.
            for mod in self.modules:
                # PASS 1
                for interval in self._time_bins[0]:
                    if len(self._time_bins[0][interval][2]) == 0:
                        continue
                    with open(self._refpath + "/event_cl/xselect.xco", 'w') as script:
                        script.write('\n')
                        script.write('\n')
                        script.write(f"read events nu{self._seqid}{mod}01_cl.evt\n")
                        script.write("./\n")
                        script.write('yes\n')
                        script.write("filter time scc\n")
                        script.write(f"{self._time_bins[0][interval][0]} , {self._time_bins[0][interval][1]}\n")
                        script.write("x\n")
                        script.write("extract events\n")
                        script.write(f"save events ./../detections/{self._dtime}-{self._snr}/nu{mod}_{self._time_bins[0][interval][0]}-{self._time_bins[0][interval][1]}.evt\n")
                        script.write('\n')
                        script.write("exit no")
                        script.write('\n')
                    subprocess.run(["xselect", "@xselect.xco"], cwd=self._evdir)
                    subprocess.run(["rm", "xselect.xco"], cwd=self._evdir)

                # PASS 2
                for interval in self._time_bins[1]:
                    if len(self._time_bins[1][interval][2]) == 0:
                        continue
                    with open(self._refpath + "/event_cl/xselect.xco", 'w') as script:
                        script.write('\n')
                        script.write('\n')
                        script.write(f"read events nu{self._seqid}{mod}01_cl.evt\n")
                        script.write("./\n")
                        script.write('yes\n')
                        script.write("filter time scc\n")
                        script.write(f"{self._time_bins[1][interval][0]} , {self._time_bins[1][interval][1]}\n")
                        script.write("x\n")
                        script.write("extract events\n")
                        script.write(f"save events ./../detections/{self._dtime}-{self._snr}/nu{mod}_{self._time_bins[1][interval][0]}-{self._time_bins[1][interval][1]}.evt\n")
                        script.write('\n')
                        script.write("exit no")
                        script.write('\n')
                    subprocess.run(["xselect", "@xselect.xco"], cwd=self._evdir)
                    subprocess.run(["rm", "xselect.xco"], cwd=self._evdir)

            # Now we merge the FPMA and FPMB data together into one file structure.
            # PASS 1
            for interval in self._time_bins[0]:
                if len(self._time_bins[0][interval][2]) == 0:
                    continue
                with open(self._refpath + f"/detections/{self._dtime}-{self._snr}/ximage.xco", 'w') as script:
                    script.write(f"read/fits/size=800/nuA_{self._time_bins[0][interval][0]}-{self._time_bins[0][interval][1]}.evt\n")
                    script.write("save_image\n")
                    script.write(f"read/fits/size=800/nuB_{self._time_bins[0][interval][0]}-{self._time_bins[0][interval][1]}.evt\n")
                    script.write("sum_image\n")
                    script.write("save_image\n")
                    script.write(f"write/fits nu_{self._time_bins[0][interval][0]}-{self._time_bins[0][interval][1]}.evt\n")
                    script.write("exit\n")
                subprocess.run(["ximage", "@ximage.xco"],cwd=self._refpath+f"./detections/{self._dtime}-{self._snr}/")
                subprocess.run(["rm", "ximage.xco"], cwd=self._refpath+f"./detections/{self._dtime}-{self._snr}/")
                # Now we perform detections.
                with open(self._refpath + f"/detections/{self._dtime}-{self._snr}/src_detect.xco", "w") as script:
                    script.write(f"read/fits/size=800/nuA_{self._time_bins[0][interval][0]}-{self._time_bins[0][interval][1]}.evt\n")
                    script.write("save_image\n")
                    script.write(f"read/fits/size=800/nuB_{self._time_bins[0][interval][0]}-{self._time_bins[0][interval][1]}.evt\n")
                    script.write("sum_image\n")
                    script.write("save_image\n")
                    script.write(f"detect/snr={self._snr}/filedet={self._time_bins[0][interval][0]}-{self._time_bins[0][interval][1]}.det/fitsdet={self.time_bins[0][interval][0]}-{self.time_bins[0][interval][1]}.fits\n")
                    script.write("exit\n")
                subprocess.run(["ximage", "@src_detect.xco"], cwd=self._refpath+f"./detections/{self._dtime}-{self._snr}/")#, 
                subprocess.run(["rm", "src_detect.xco"], cwd=self._refpath+f"./detections/{self._dtime}-{self._snr}/")
            
            # PASS 2
            for interval in self._time_bins[1]:
                if len(self._time_bins[1][interval][2]) == 0:
                    continue
                with open(self._refpath + f"/detections/{self._dtime}-{self._snr}/ximage.xco", 'w') as script:
                    script.write(f"read/fits/size=800/nuA_{self._time_bins[1][interval][0]}-{self._time_bins[1][interval][1]}.evt\n")
                    script.write("save_image\n")
                    script.write(f"read/fits/size=800/nuB_{self._time_bins[1][interval][0]}-{self._time_bins[1][interval][1]}.evt\n")
                    script.write("sum_image\n")
                    script.write("save_image\n")
                    script.write(f"write/fits nu_{self._time_bins[1][interval][0]}-{self._time_bins[1][interval][1]}.evt\n")
                    script.write("exit\n")
                subprocess.run(["ximage", "@ximage.xco"],cwd=self._refpath+f"./detections/{self._dtime}-{self._snr}/")
                subprocess.run(["rm", "ximage.xco"], cwd=self._refpath+f"./detections/{self._dtime}-{self._snr}/")
                # Now we perform detections.
                with open(self._refpath + f"/detections/{self._dtime}-{self._snr}/src_detect.xco", "w") as script:
                    script.write(f"read/fits/size=800/nuA_{self._time_bins[1][interval][0]}-{self._time_bins[1][interval][1]}.evt\n")
                    script.write("save_image\n")
                    script.write(f"read/fits/size=800/nuB_{self._time_bins[1][interval][0]}-{self._time_bins[1][interval][1]}.evt\n")
                    script.write("sum_image\n")
                    script.write("save_image\n")
                    script.write(f"detect/snr={self._snr}/filedet={self._time_bins[1][interval][0]}-{self._time_bins[1][interval][1]}.det/fitsdet={self.time_bins[1][interval][0]}-{self.time_bins[1][interval][1]}.fits\n")
                    script.write("exit\n")
                subprocess.run(["ximage", "@src_detect.xco"], cwd=self._refpath+f"./detections/{self._dtime}-{self._snr}/")#, 
                                            #capture_output=True, text=True)
                subprocess.run(["rm", "src_detect.xco"], cwd=self._refpath+f"./detections/{self._dtime}-{self._snr}/")
        det_files = glob.glob(f"{self._refpath}/detections/{self._dtime}-{self._snr}/*.det")
        file_strings = []
        for file in det_files:
            file_strings.append(file.replace(f"{self._refpath}/detections/{self._dtime}-{self._snr}/", ''))
        script_string = "srcmrg/out=mrg.txt/tolerance=5e1"
        for file in file_strings:
            script_string += f" {file}"
        script_string += "\n"
        if len(file_strings) == 0:
            print("No detections found!")
        else:
            with open(self._refpath + f"/detections/{self._dtime}-{self._snr}/src_merge.xco", 'w') as script:
                script.write(script_string)
                script.write("exit\n")
            subprocess.run(["ximage", "@src_merge.xco"], cwd=self._refpath+f"detections/{self._dtime}-{self._snr}/")
            subprocess.run(["rm", "src_merge.xco"], cwd=self._refpath+f"detections/{self._dtime}-{self._snr}/")
            self._detections = self.read_detections()


    def read_detections(self):
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

        with open(self._refpath + f"/detections/{self._dtime}-{self._snr}/mrg.txt") as detections:
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


    #def plot_temporal_differences(self, save_fig=False, display_fig=True):
    #    """
    #    Plots the elapsed livetime for each event within the observation data.

    #    Arguments
    #    ---------
    #    save_fig : bool
    #        a flag which determines whether or not to save an image under the file name:
    #        `elapsed_livetime.pdf`
    #    display_fig : bool
    #        a flag which determines whether of not to display an image when calling this method

    #    Returns
    #    -------
    #    None

    #    Raises
    #    ------
    #    None
    #    """

    #    plt.rcParams["font.family"] = "monospace"
    #    fig, ax = plt.subplots()
    #    ax.hist(self.header_info["DTIMES"], bins=10000)
    #    ax.set_xlim(xmin=0, xmax = 0.1)
    #    ax.set_title("Time since last detection")
    #    ax.set_xlabel("Time (s)")
    #    ax.set_ylabel("Number")
    #    if save_fig:
    #        plt.savefig(self.path + f"event_cl/elapsed_livetime.pdf", dpi=1000)
    #    if display_fig:
    #        plt.show()
    

    def extract_detections(self):
        """
        This method runs `nuproducts` on all valid detections found within the /detections/ directory

        Arguments
        ---------
        None

        Returns
        -------
        None

        Raises
        ------
        None
        """

        self._detections = self.read_detections()
        for idx in self.detect_info["INDEX"]:
            coords = self.ra_dec_todeg(self.detections["RA"][int(idx) - 1], self.detections["DEC"][int(idx) - 1])
            string = f"nuproducts indir=./event_cl instrument=FPMA steminputs=nu{self.obsid} outdir=./products/{idx} srcra={coords[0]} srcdec={coords[1]} bkgra={coords[0]} bkgdec={coords[1]} binsize=5"
            string = f"nuproducts indir=./event_cl instrument=FPMB steminputs=nu{self.obsid} outdir=./products/{idx} srcra={coords[0]} srcdec={coords[1]} bkgra={coords[0]} bkgdec={coords[1]} binsize=5"
            subprocess.run(string.split(), cwd=self._refpath)


    def remove_main_source(self):
        pass


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


    def ds9_detection(self, mod='B', n=1):
        """
        Spawns an instance of ds9 which showcases a specific detection.

        Arguments
        ---------
        n : int
            The index associated with this element

        Returns:
        --------

        Raises
        ------
        """
        pass


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
    

    def acquire_lightcurve(self, n=1):
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

        fits_file = self._refpath + f"products/{n}/nu{self._seqid}B01.flc"
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


    def acquire_pha(self, n=1):
        """
        """
        pha_file = self._refpath + f"products/{n}/nu{self._seqid}B01_sr.pha"
        with fits.open(pha_file) as hdul:
            pha_data = hdul[1].data
        return pha_data


    def plot_lightcurve(self, n=1):
        """
        Plots the background subtracted lightcurve of a detection.

        Arguments
        ---------
        n : `int`
            The index of the detection being plotted.

        Returns
        -------

        Raises
        ------
        """

        lc_data = self.acquire_lightcurve(n)
        fig, ax = plt.subplots()
        ax.errorbar(lc_data["TIME"], lc_data["RATE"], lc_data["RATEERR"], fmt='o', lw=0.5, markersize=2)
        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Counts/Time (c/s)")
        plt.show()


    def plot_pha(self, n=1):
        """
        """

        pha_data = self.acquire_pha(n)
        fig, ax = plt.subplots()
        y_dat = []
        x_dat = []
        for datum in pha_data:
            x_dat.append(datum[0] * 0.04 + 1.6)
            y_dat.append(datum[1]/44034)
        ax.plot(x_dat, y_dat)
        ax.set_xlim(3, 79)
        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Counts/Time (c/s)")
        plt.show()
    

    def generate_background_rate(self):
        """
        This method calls an algorithm for identifying the background rates for a given fits image 
        by energy range. Four background regions are defined for a given FITS image and these 
        background regions individually are analyzed to determine a distribution of count rates. 
        The median of the four regions is then taken and turned into a background rate distribution.
        """
        pass

    
    def generate_region(self, path, xpix, ypix):
        """
        This method generates a region file to be used for data analysis with a fits image.
        """
        pass

    def generate_directory(self, overwrite=False):
        """
        Erases and creates a new directory
        """
        pass

