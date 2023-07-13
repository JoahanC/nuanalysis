import os
import subprocess
import glob
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS


class NuObs:
    """
    This class defines a basic shell object for performing analysis on NuSTAR observations using an installation of HEASoft, 
    DS9, and the NuSTAR calibration database.
    """

    def __init__(self, obsid, dt=1000, snr=3):
        self.obsid = obsid
        self.path = f"./{obsid}/"
        self.contents = os.listdir(self.path)
        if "event_cl" not in self.contents:
            print("Clean science events not found. Running nupipeline cleaning process.")
            self.run_cleaning_script()
        self.header_info = self.generate_headerinfo()
        self.read_in_phi()
        self.dt = dt
        self.snr = snr
        self.time_bins = self.generate_timebins()
        #self.ra = self.ra_dec_todeg()
        #self.dec = self.ra_dec_todeg()
        self.detect_info = {}
        
        
    def read_in_phi(self):
        """
        Makes use of two files `nustar_pilow.txt` and `nustar_pihi.txt` within the `nustar` directory to set PI channel thresholds 
        for data analysis using NUSTARDAS

        Arguments
        ---------
        None

        Returns
        -------
        None

        Raises
        ------
        IndexError
            when the PI channel bounds are not of the same length
        """

        self.phi_bounds = []
        low_phi = open("nustar_pilow.txt", 'r')
        high_phi = open("nustar_pihi.txt", 'r')

        low_phis = low_phi.readlines()
        high_phis = high_phi.readlines()
        
        if len(low_phis) != len(high_phis):
            raise IndexError("PHI bounds must have the same length!")
        
        for idx in range(len(low_phis)):
            self.phi_bounds.append((low_phis[idx].replace('\n', ''), high_phis[idx].replace('\n', '')))


    def run_cleaning_script(self):
        """
        In the event that a NuSTAR observation is recently taken, or has not had a proper data reduction process enacted, this 
        method will run `nupipeline` on the unprocessed data to generate an `event_cl` directory with reduced data following 
        Stage 1 and Stage 2 calibration.

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
        
        output_path = f"./{self.obsid}/event_cl/"
        os.mkdir(output_path)
        subprocess.run(["nupipeline", self.path, f"nu{self.obsid}", output_path])

    
    def generate_headerinfo(self):
        """
        This function reads in the header information of the A01 clean event file for the NuSTAR observation and catalogs
        values which may be of relevance during data analysis.

        Arguments
        ---------
        None

        Returns
        -------
        dict
            a dictionary containing the following keys: `OBJECT`, `DATE-OBS`, `RA_OBJ`, `DEC_OBJ`, `TSTART`, `TSTOP`, `TELAPSE`.
            These keys represent the following values:
            
            OBJECT : accepted name for this object
            DATE-OBS : the date during which the observation of this object was taken
            RA_OBJ : the right-ascension at which the object was detected
            DEC_OBJ : the declination at which the object was detected
            TSTART : the time at which the observation began
            TSTOP : the time at which the observation ended
            TELAPSE : the total time which has elapsed during the observation
            TIMES : the times at which an event occured
            DTIMES : the elapsed time between two events

        Raises
        ------
        None
        """

        header_info = {}
        keys = ["OBJECT", "DATE-OBS", "RA_OBJ", "DEC_OBJ", "TSTART", "TSTOP", "TELAPSE"]
        times = []
        time_diffs = []

        with fits.open(self.path + f"event_cl/nu{self.obsid}A01_cl.evt") as hdul:
            for array in hdul[1].data:
                times.append(array[0])
                time_diffs.append(array[1])
            for key in keys:
                header_info[key] = hdul[1].header[key]
            self.wcs = WCS(hdul[1].header)
            header_info["TIMES"] = times
            header_info["DTIMES"] = time_diffs
        
        return header_info


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

        obs_times = self.header_info["TIMES"]


        # generate intervals for PASS 1
        intervals_p1 = []
        cur_int = obs_times[0]
        max_int = np.max(obs_times)
        while cur_int < max_int:
            intervals_p1.append(cur_int)
            cur_int += self.dt
        intervals_p1.append(max_int)

        # generate intervals for PASS 2
        intervals_p2 = []
        intervals_p2.append(obs_times[0])
        cur_int = obs_times[0] + self.dt / 2
        max_int = np.max(obs_times)
        while cur_int < max_int:
            intervals_p2.append(cur_int)
            cur_int += self.dt 
        intervals_p2.append(max_int)


        # split data for PASS 1
        idx = 0
        first = True
        data_intervals_1 = {}
        last = obs_times[-1]
        for datapoint in obs_times:
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
        last = obs_times[-1]
        for datapoint in obs_times:
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

        if "detections" not in self.contents:
            os.mkdir(self.path + "detections/")
        detect_path = self.path + f"detections/{self.dt}-{self.snr}/"
        if f"{self.dt}-{self.snr}" in os.listdir(self.path + "detections/"):
            self.detect_info = self.read_detections()
        else:
            os.mkdir(detect_path)
            # Begin by selecting the data which lies in the appropriate energy range.
            for letter in ["A", "B"]:
                # PASS 1
                for interval in self.time_bins[0]:
                    if len(self.time_bins[0][interval][2]) == 0:
                        continue
                    with open(self.path + "event_cl/xselect.xco", 'w') as script:
                        script.write('\n')
                        script.write('\n')
                        script.write(f"read events nu{self.obsid}{letter}01_cl.evt\n")
                        script.write("./\n")
                        script.write('yes\n')
                        script.write("filter time scc\n")
                        script.write(f"{self.time_bins[0][interval][0]} , {self.time_bins[0][interval][1]}\n")
                        script.write("x\n")
                        script.write("extract events\n")
                        script.write(f"save events ./../detections/{self.dt}-{self.snr}/nu{letter}_{self.time_bins[0][interval][0]}-{self.time_bins[0][interval][1]}.evt\n")
                        script.write('\n')
                        script.write("exit no")
                        script.write('\n')
                    subprocess.run(["xselect", "@xselect.xco"], cwd=self.path+"event_cl/")
                    subprocess.run(["rm", "xselect.xco"], cwd=self.path+"event_cl/")

                # PASS 2
                for interval in self.time_bins[1]:
                    if len(self.time_bins[1][interval][2]) == 0:
                        continue
                    with open(self.path + "event_cl/xselect.xco", 'w') as script:
                        script.write('\n')
                        script.write('\n')
                        script.write(f"read events nu{self.obsid}{letter}01_cl.evt\n")
                        script.write("./\n")
                        script.write('yes\n')
                        script.write("filter time scc\n")
                        script.write(f"{self.time_bins[1][interval][0]} , {self.time_bins[1][interval][1]}\n")
                        script.write("x\n")
                        script.write("extract events\n")
                        script.write(f"save events ./../detections/{self.dt}-{self.snr}/nu{letter}_{self.time_bins[1][interval][0]}-{self.time_bins[1][interval][1]}.evt\n")
                        script.write('\n')
                        script.write("exit no")
                        script.write('\n')
                    subprocess.run(["xselect", "@xselect.xco"], cwd=self.path+"event_cl/")
                    subprocess.run(["rm", "xselect.xco"], cwd=self.path+"event_cl/")

            # Now we merge the FPMA and FPMB data together into one file structure.
            # PASS 1
            for interval in self.time_bins[0]:
                if len(self.time_bins[0][interval][2]) == 0:
                    continue
                with open(self.path + f"detections/{self.dt}-{self.snr}/ximage.xco", 'w') as script:
                    script.write(f"read/fits/size=800/nuA_{self.time_bins[0][interval][0]}-{self.time_bins[0][interval][1]}.evt\n")
                    script.write("save_image\n")
                    script.write(f"read/fits/size=800/nuB_{self.time_bins[0][interval][0]}-{self.time_bins[0][interval][1]}.evt\n")
                    script.write("sum_image\n")
                    script.write("save_image\n")
                    script.write(f"write/fits nu_{self.time_bins[0][interval][0]}-{self.time_bins[0][interval][1]}.evt\n")
                    script.write("exit\n")
                subprocess.run(["ximage", "@ximage.xco"],cwd=self.path+f"./detections/{self.dt}-{self.snr}/")
                subprocess.run(["rm", "ximage.xco"], cwd=self.path+f"./detections/{self.dt}-{self.snr}/")
                # Now we perform detections.
                with open(self.path + f"detections/{self.dt}-{self.snr}/src_detect.xco", "w") as script:
                    script.write(f"read/fits/size=800/nuA_{self.time_bins[0][interval][0]}-{self.time_bins[0][interval][1]}.evt\n")
                    script.write("save_image\n")
                    script.write(f"read/fits/size=800/nuB_{self.time_bins[0][interval][0]}-{self.time_bins[0][interval][1]}.evt\n")
                    script.write("sum_image\n")
                    script.write("save_image\n")
                    script.write(f"detect/snr={self.snr}/filedet={self.time_bins[0][interval][0]}-{self.time_bins[0][interval][1]}.det/fitsdet={self.time_bins[0][interval][0]}-{self.time_bins[0][interval][1]}.fits\n")
                    script.write("exit\n")
                subprocess.run(["ximage", "@src_detect.xco"], cwd=self.path+f"./detections/{self.dt}-{self.snr}/")#, 
                                            #capture_output=True, text=True)
                subprocess.run(["rm", "src_detect.xco"], cwd=self.path+f"./detections/{self.dt}-{self.snr}/")
            
            # PASS 2
            for interval in self.time_bins[1]:
                if len(self.time_bins[1][interval][2]) == 0:
                    continue
                with open(self.path + f"detections/{self.dt}-{self.snr}/ximage.xco", 'w') as script:
                    script.write(f"read/fits/size=800/nuA_{self.time_bins[1][interval][0]}-{self.time_bins[1][interval][1]}.evt\n")
                    script.write("save_image\n")
                    script.write(f"read/fits/size=800/nuB_{self.time_bins[1][interval][0]}-{self.time_bins[1][interval][1]}.evt\n")
                    script.write("sum_image\n")
                    script.write("save_image\n")
                    script.write(f"write/fits nu_{self.time_bins[1][interval][0]}-{self.time_bins[1][interval][1]}.evt\n")
                    script.write("exit\n")
                subprocess.run(["ximage", "@ximage.xco"],cwd=self.path+f"./detections/{self.dt}-{self.snr}/")
                subprocess.run(["rm", "ximage.xco"], cwd=self.path+f"./detections/{self.dt}-{self.snr}/")
                # Now we perform detections.
                with open(self.path + f"detections/{self.dt}-{self.snr}/src_detect.xco", "w") as script:
                    script.write(f"read/fits/size=800/nuA_{self.time_bins[1][interval][0]}-{self.time_bins[1][interval][1]}.evt\n")
                    script.write("save_image\n")
                    script.write(f"read/fits/size=800/nuB_{self.time_bins[1][interval][0]}-{self.time_bins[1][interval][1]}.evt\n")
                    script.write("sum_image\n")
                    script.write("save_image\n")
                    script.write(f"detect/snr={self.snr}/filedet={self.time_bins[1][interval][0]}-{self.time_bins[1][interval][1]}.det/fitsdet={self.time_bins[1][interval][0]}-{self.time_bins[1][interval][1]}.fits\n")
                    script.write("exit\n")
                subprocess.run(["ximage", "@src_detect.xco"], cwd=self.path+f"./detections/{self.dt}-{self.snr}/")#, 
                                            #capture_output=True, text=True)
                subprocess.run(["rm", "src_detect.xco"], cwd=self.path+f"./detections/{self.dt}-{self.snr}/")
        det_files = glob.glob(f"{self.path}/detections/{self.dt}-{self.snr}/*.det")
        file_strings = []
        for file in det_files:
            file_strings.append(file.replace(f"{self.path}/detections/{self.dt}-{self.snr}/", ''))
        script_string = "srcmrg/out=mrg.txt/tolerance=5e1"
        for file in file_strings:
            script_string += f" {file}"
        script_string += "\n"
        if len(file_strings) == 0:
            print("No detections found!")
        else:
            with open(self.path + f"detections/{self.dt}-{self.snr}/src_merge.xco", 'w') as script:
                script.write(script_string)
                script.write("exit\n")
            subprocess.run(["ximage", "@src_merge.xco"], cwd=self.path+f"detections/{self.dt}-{self.snr}/")
            subprocess.run(["rm", "src_merge.xco"], cwd=self.path+f"detections/{self.dt}-{self.snr}/")
            self.detect_info = self.read_detections()


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

        with open(self.path + f"detections/{self.dt}-{self.snr}/mrg.txt") as detections:
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


    def plot_temporal_differences(self, save_fig=False, display_fig=True):
        """
        Plots the elapsed livetime for each event within the observation data.

        Arguments
        ---------
        save_fig : bool
            a flag which determines whether or not to save an image under the file name:
            `elapsed_livetime.pdf`
        display_fig : bool
            a flag which determines whether of not to display an image when calling this method

        Returns
        -------
        None

        Raises
        ------
        None
        """

        plt.rcParams["font.family"] = "monospace"
        fig, ax = plt.subplots()
        ax.hist(self.header_info["DTIMES"], bins=10000)
        ax.set_xlim(xmin=0, xmax = 0.1)
        ax.set_title("Time since last detection")
        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Number")
        if save_fig:
            plt.savefig(self.path + f"event_cl/elapsed_livetime.pdf", dpi=1000)
        if display_fig:
            plt.show()
    

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

        self.detect_info = self.read_detections()
        for idx in self.detect_info["INDEX"]:
            coords = self.ra_dec_todeg(self.detect_info["RA"][int(idx) - 1], self.detect_info["DEC"][int(idx) - 1])
            string = f"nuproducts indir=./event_cl instrument=FPMA steminputs=nu{self.obsid} outdir=./products/{idx} srcra={coords[0]} srcdec={coords[1]} bkgra={coords[0]} bkgdec={coords[1]} binsize=5"
            string = f"nuproducts indir=./event_cl instrument=FPMB steminputs=nu{self.obsid} outdir=./products/{idx} srcra={coords[0]} srcdec={coords[1]} bkgra={coords[0]} bkgdec={coords[1]} binsize=5"
            subprocess.run(string.split(), cwd=self.path)


    def remove_main_source(self):
        print(self.header_info(""))


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


    def ds9_detections(self, mod='B'):
        """
        Displays all of the detections that were processed by nuproducts into an instance of ds9.
        """
        ds9_string = f"ds9 nu{self.obsid}{mod}01_cl.evt "
        detection_dirs = os.listdir(self.path + "products/")
        for dir in detection_dirs:
            ds9_string += f"-regions ../products/{dir}/nu{self.obsid}{mod}01_bkg.reg "
            ds9_string += f"-regions ../products/{dir}/nu{self.obsid}{mod}01_src.reg "            
        subprocess.run(ds9_string.split(), cwd=self.path + "event_cl/")
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

        fits_file = self.path + f"products/{n}/nu{self.obsid}B01.flc"
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
        pha_file = self.path + f"products/{n}/nu{self.obsid}B01_sr.pha"
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
        with open(self.path/)

    
    def generate_region(self, path, xpix, ypix):
        """
        This method generates a region file to be used for data analysis with a fits image.
        """

    def generate_directory(self, overwrite=False)
        



    



#recover_events()
#current_obs = [glob.glob("*/")[0]]
#for obs in current_obs:
#    print(obs[:-1])
#    test = NuObs(obs[:-1], dt=10000)
#    test.generate_detections()

test = NuObs("30501002002", 10000)
test.plot_pha(7)
