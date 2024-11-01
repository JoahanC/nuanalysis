from astropy.time import Time, TimeDelta
from astropy.wcs import WCS
import astropy.units as u
import os
import subprocess
import glob
import numpy as np


class NuSTAR():
    '''
    Class for holding constant attributes about NuSTAR, time conversion from
    MET to 'TIME' objects and vice-versa, and energy conversions from PI to keV 
    and vice-versa
    '''
    
    def __init__(self):
        self._ref_epoch = Time('2010-01-01T00:00:00', format='fits', scale='utc')
        self._raw_pixel = 604.8 * u.micron
        self._pixel_um = self._raw_pixel / 5.
        self._pixel = 2.54 * u.arcsec
        self._launch = Time('2012-06-13T00:00:00')
        self._tick = 16./14745600. 

    # Mutable properties begin below

    @property
    def launch(self):
        '''
        Returns the launch date
        '''
        
        return self._launch
        

    @property
    def ref_epoch(self):
        '''
        Returns MET time reference epoch (Jan 1, 2010, UTC)
        '''
        
        return self._ref_epoch


    @property
    def pixel(self):
        '''
        Returns the size of a sky pixel (2.54 arcsec)
        '''
        
        return self._pixel


    @property
    def pixel_um(self):
        '''
        Returns the physical size of a NuSTAR sub-pixel (raw_pixel / 5)
        '''
        
        return self._pixel_um


    @property
    def tick(self):
        '''
        Returns the clock resolution, nominally 16 samples at 14.7456 Mhz
        '''
        
        return self._tick

    # Methods begin below

    def time_to_met(self, time):
        ''' 
        Convert a Time object to a unitless NuSTAR MET second. Note that the time scale
        for MET seconds is *always* 'TT' seconds.
        
        Parameters
        ----------
        time : Astropy Time object
            The time or array of times that you want to convert to MET.
        
        Returns
        -------
        met : float
            Seconds since Jan 1, 2010 UTC in TT seconds.
        
        Examples
        --------
        >>> ns = NuSTAR()
        >>> time = Time('2020-01-01T12:34:42', format = 'fits', scale = 'utc')
        >>> met = ns.time_to_met(time)
        >>> print(met)
        315578085.0
        '''
        
        met = (time.tt - self.ref_epoch).sec
        return met
    

    def met_to_time(self, met):
        '''
        Assumes unitless MET seconds input.
        
        Parameters
        ----------
        met: float
            Unitless NuSTAR MET seconds.
        
        Examples
        --------
        >>> ns = NuSTAR()
        >>> met = 315578085.0
        >>> time = ns.met_to_time(met)
        >>> print(time.fits)
        2020-01-01T12:34:42.000
        '''

        # Hacky way to check and see if the entries are floats
        from numpy import asarray
        foo = asarray(met)
        if foo.ndim == 0:
            assert isinstance(met, float), "met_to_time: met must be a float"
        else:
            assert isinstance(met[0], float), "met_to_time: met must be a float"
        
        this_time = TimeDelta(met, format='sec', scale ='tt') + self.ref_epoch
        return this_time


    def rate_conversion(self, rate, incident=False):
        '''
        Converts observed counts per second into an incident rate
    
        Parameters
        ----------
        rate: float
            Unitless NuSTAR MET seconds.
        incident: boolean
            If False, converts between measured and incident rates
            If True, converts between incident and measured rates.
    
        Returns
        -------
        float
            Measured/incident rate if the incident parameters is False/True
    
        Examples
        --------
        >>> ns = NuSTAR()
        >>> met = 315578085.0
        >>> time = ns.met_to_time(met)
        >>> print(time.fits)
        2020-01-01T12:34:42.000
        '''

        if incident is True:
            # Convert from incident to meaured rate
            result = rate / (1.0 + rate * 2.5e-3)
        else:
            # Convert from measured to incident rate
            result = rate / (1.0 - rate * 2.5e-3)
        return result
    
    
    def energy_to_channel(self, keV):
        """
        Convert keV to NuSTAR PI channels: PI = (keV - 1.6) / 0.04
        
        Parameters
        ----------
        
        keV: float
            Unitless float giving the energy you want to convert to channels.
        
        Example
        ---------
        >>> from nustar_gen.utils import energy_to_chan
        >>> chan = energy_to_chan(10.)
        >>> np.isclose(chan, 210)
        True
        """
        
        return int((keV - 1.6) / 0.04)


    def channel_to_energy(self, chan):
        """
        Convert NuSTAR PI channels to keV: keV = PI * 0.04 + 1.6
        
        
        Parameters
        ----------
        
        channel: int
            Integer PI channel to be converted into keV    
        
        Example
        -------
        >>> from nustar_gen.utils import chan_to_energy
        >>> chan = chan_to_energy(210)
        >>> np.isclose(chan, 10.)
        True
        """

        try:
            en = [float(x) *  0.04 + 1.6 for x in chan]
        except:
            en = chan * 0.04 + 1.6
        return en



class Observation():
    '''
    Class for storing meta-data about a given NuSTAR observation.
    
    Parameters
    ----------
    path : str
        The top-level working directory. All paths are assumed to be relative to this
        location. This should be one level above the sequence ID location.
    
    seqid : str
        The sequence id for the observation
 
    evdir : str
        Full path to the directory containing the event files
        
    out_path : str
        Full path to the desired output location
    '''
        
    def __init__(self, path='./', seqid=False, evdir=False,
                out_path=False):
                                
        self._evdir_lock = False
        self.modules = ['A', 'B']
        self.set_path(path)

        # If you specify the evdir location, make sure nothing can change this
        if evdir is not False:
            evdir = evdir#os.path.abspath(evdir)
            self._set_evdir(evdir, lock=True)
        if seqid is False:
            self._seqid=False
        else:
            self._set_seqid(seqid)
        
        if out_path is False:
            self.set_outpath(self.evdir)
        else:
            out_path = os.path.abspath(out_path)
            self.set_outpath(out_path)

    # Mutable properties begin below

    @property
    def path(self):
        '''
        Returns the top-level path
        '''
        return self._path


    @property
    def observation_date(self):
        '''
        Returns the date of the observation parsed from the file headers
        '''
        return self._observation_date

    
    @property
    def evdir(self):
        '''
        Returns the event file directory
        '''
        return self._evdir


    @property
    def seqid(self):
        '''
        Returns the current sequence ID
        '''
        return self._seqid
        

    @property
    def exposure(self):
        '''
        Returns an dict (with 'A' and 'B' as keys) with lists of exposures for
        all event files.
        '''
        return self._exposure        
    

    @property
    def source_position(self):
        '''
        Returns the current source RA/Dec from the FITS headers as Astropy
        SkyCoord object
        '''
        return self._source_position


    @property
    def science_files(self):
        '''
        Returns a list of science (01) event files
        '''
        
        if self._seqid is False:
            raise ValueError(f"Set sequence ID and path first!")
        
        science_files = {}
        for mod in self.evtfiles:
            science_files[mod] = []
            for file in self.evtfiles[mod]:
                if ( (f'{mod}01' in file) | (f'{mod}06' in file) ) &  \
                    ( (file.endswith('gz') ) | \
                    (file.endswith('evt') ) ):
                    science_files[mod].extend([file])
            
            # If you have 06 data, make sure you either have CHU split data
            # or not:
            # See if you have any CHU mode 6 data:
            for file in science_files[mod]:
                if f'{mod}06_chu' in file:
                    # You do, so find the "non CHU" mode06 file:
                    for ind, clfile in enumerate(science_files[mod]):
                        if f'{mod}06_cl' in clfile:
                            science_files[mod].pop(ind)
                            break
        return science_files


    @property
    def out_path(self):
        '''
        Returns the output path.
        '''
        return self._out_path
    
    # Methods begin below

    def _set_evdir(self, value, lock=False):
        '''
        Set the event directory.
        '''
        
        if self._evdir_lock is False:
            self._evdir = value
        else:
            # Try the default:
            self._evdir=os.path.join(self._path, )
            
        assert os.path.isdir(self._evdir), f"Event file path does not exist! {self._evdir}"
       
        # Set the lock flag
        if lock is True:
            self._evdir_lock=True
        return


    def set_path(self, path):
        '''
        Sets the path. Makes sure that path ends with a '/'
        '''
        path = os.path.abspath(path)
        self._path = path
        return


    def _set_seqid(self, value):
        '''Set the sequence ID. Raise error if obs.path+obs.seqid doesn't exist. Finds
        clean event files and parses the input fits header to populate the other
        Observation() attributes.'''
        
        self._seqid=value
        
        if self._evdir_lock is False:
            default_datapath = os.path.join(self._path, self._seqid)
            self._set_evdir(os.path.join(default_datapath, 'event_cl'))
        
        # Check to make sure this exiss:
        assert os.path.isdir(self.evdir), f'Event file path does not exist! {self.evdir}'
        self._find_cleaned_files()
        self._parse_header()


    def set_outpath(self, value):
        '''
        Set the output path.
        '''
        
        self._out_path=value
        assert os.path.isdir(self._out_path), f"Output path does not exist! {self._out_path}"
        
        return


    def _find_cleaned_files(self):
        '''
        Uses self.evdir to find all of the cleaned event files.
        '''
        
        self.evtfiles = {}
        for mod in self.modules:
            self.evtfiles[mod] = sorted(glob.glob(os.path.join(self._evdir, f'nu*{mod}*cl.evt')))            
        
        return
        

    def _parse_header(self):
        from astropy.io.fits import getheader, getdata
        from astropy.coordinates import SkyCoord
        
        self._exposure = {}
        self._event_times = []
        for mod in self.modules:
            for evtfile in self.science_files[mod]:
                hdr = getheader(evtfile)
                hdr2 = getheader(evtfile, 1)
                self.wcs = WCS(hdr2)
                self._observation_date = Time(hdr['DATE-OBS'], format='fits')
                self._source_position = \
                        SkyCoord(hdr['RA_OBJ'], hdr['DEC_OBJ'], unit='deg')
                
                for ti in range(7):
                    keystr = f'{mod}'+f'{ti+1}'.zfill(2)
                    if keystr in evtfile:
                        if keystr not in self._exposure:
                            self._exposure[keystr] = hdr['EXPOSURE']
                        else:
                            self._exposure[keystr] += hdr['EXPOSURE']
        evt_array = getdata(self.science_files[mod][0])
        self._data = evt_array
        for evt in evt_array:
            self._event_times.append(evt[0])
        
        return
        

    def exposure_report(self):
        '''
        Make the output report on the exposure for various observation types
        '''
        
        if self._seqid is False:
            raise ValueError(f"Set sequence ID and path first!")
        
        for mod in self.modules:
            for ti in range(7):
                keystr = f'{mod}'+f'{ti+1}'.zfill(2)
                if keystr in self._exposure:
                    print(f'Exposure for FPM{mod}, mode '+f'{ti+1}'.zfill(2)+f' is: {1e-3*self._exposure[keystr]:10.4} ks')
            print()

        return


    def download_bgd_report(self):
        '''
        Downloads the background report from the SOC:
        '''
        
        base_html = 'http://www.srl.caltech.edu/NuSTAR_Public/NuSTAROperationSite/'
        base_html += 'SAA_Filtering/nulyses_reports/'
        for mod in self.modules:
            pdf_html = base_html + f'{self.seqid}/nu{self.seqid}_SAA_Report_{mod}.pdf'
            subprocess.run((f'wget {pdf_html}').split())
            subprocess.run((f'mv {os.path.basename(pdf_html)} {self.out_path}').split())
        return
    

    