#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Time-stamp: "2017-10-24 12:33:35 jlenain"

"""
FLaapLUC (Fermi/LAT automatic aperture photometry Light C<->Urve)

Automatic generation of aperture photometric light curves of
high energy sources, for a given source.

No likelihood fit is performed, the results solely rely on the
4FGL spectral fits, if available.

More information are available at:
http://fermi.gsfc.nasa.gov/ssc/data/analysis/scitools/aperture_photometry.html

@author Jean-Philippe Lenain <mailto:jlenain@in2p3.fr>
"""

import datetime
import glob
import logging
import matplotlib
matplotlib.use('Agg')

import os
import sys
import time
import numpy as np
from ConfigParser import ConfigParser
from matplotlib import pyplot as plt
from matplotlib.ticker import FuncFormatter
from matplotlib.dates import date2num, HourLocator, DayLocator, DateFormatter
from matplotlib import patheffects

import ephem
from astropy.io import ascii
from astropy.io import fits
from astropy import coordinates as c
from astropy import units as u
from astropy import time as t

import gt_apps as fermi
from flaapluc import extras


# Flags
BATCH = True
# Flag to know whether Gamma is assumed to be ASSUMEDGAMMA
# or taken from the 4FGL.
FLAGASSUMEDGAMMA = False

# Global variables
TOFFSET = 54000.  # offset in MJD for plot creation
# assumed photon index for a source not belonging to the 4FGL
ASSUMEDGAMMA = -2.5

# Visibility parameters
ASTRONOMICAL_TWILIGHT_ANGLE = -16. * u.deg
CIVIL_TWILIGHT_ANGLE = -12. * u.deg
MOON_TWILIGHT_ANGLE = -0.83 * u.deg

def getConfigList(option, sep=','):
    return [stuff for stuff in option.split(sep)]


def processSrc(mysrc=None, useThresh=False, daily=False, mail=True, longTerm=False, test=False, yearmonth=None,
               mergelongterm=False, withhistory=False, update=False, configfile='default.cfg', force_daily=False,
               stopmonth=None, stopday=None, forcealert=False, too=False, log=logging.INFO):
    """
    Process a given source.
    """

    logging.basicConfig(format='[%(levelname)s] %(message)s', level=log)

    if mysrc is None:
        logging.error('Missing input source !')
        sys.exit(1)

    if too=="True":
        logging.info('[{src:s}] \033[94m*** Dealing with a ToO, forcing daily-binned light curve, and alert\033[0m'.format(src=mysrc))
        forcealert = True
        force_daily = True

    # If we asked for a daily light curve, first make sure that the long time-binned data already exists, otherwise this script will crash, since the daily-binned PNG needs the long time-binned data to be created. No mail alert is sent at this step.
    # We automatically recreate here any missing long time-binned data.
    if daily and not longTerm and not force_daily:
        logging.info('[{src:s}] Daily light curve asked for, I will first process the long time-binned one'.format(src=mysrc))
        longtermactive, visible = processSrc(mysrc=mysrc,
                                             useThresh=useThresh,
                                             daily=False,
                                             mail=False,
                                             longTerm=longTerm,
                                             yearmonth=yearmonth,
                                             mergelongterm=mergelongterm,
                                             withhistory=withhistory,
                                             update=update,
                                             configfile=configfile,
                                             stopmonth=stopmonth,
                                             stopday=stopday,
                                             forcealert=forcealert,
                                             too=too,
                                             log=log)
        if longtermactive and visible:
            logging.info('[{src:s}] Source {src:s} is active and visible in long time-binned data, processing daily-binned light curve...'.format(src=mysrc))
        elif longtermactive and not visible:
            logging.info('[{src:s}] \033[91mSource {src:s} is active but not visible. Daily-binned light curve aborted...\033[0m'.format(src=mysrc))
            return False
        elif not longtermactive and visible:
            logging.info('[{src:s}] \033[91mSource {src:s} is visible but not active. Daily-binned light curve aborted...\033[0m'.format(src=mysrc))
            return False
        elif not longtermactive and not visible:
            logging.info('[{src:s}] \033[91mSource {src:s} is neither active nor visible. Daily-binned light curve aborted...\033[0m'.format(src=mysrc))
            return False
        else:
            logging.info('[{src:s}] \033[91mDaily-binned light curve aborted, for unknown reason...\033[0m'.format(src=mysrc))
            return False
    elif force_daily:
        logging.info('[{src:s}] Forcing daily light curve, I will first process the long time-binned one'.format(src=mysrc))
        longtermactive, visible = processSrc(mysrc=mysrc,
                                             useThresh=useThresh,
                                             daily=False,
                                             mail=False,
                                             longTerm=longTerm,
                                             yearmonth=yearmonth,
                                             mergelongterm=mergelongterm,
                                             withhistory=withhistory,
                                             update=update,
                                             configfile=configfile,
                                             stopmonth=stopmonth,
                                             stopday=stopday,
                                             forcealert=forcealert,
                                             log=log)
    else:
        logging.info('[{src:s}] Processing long time-binned light curve...'.format(src=mysrc))

    auto = automaticLightCurve(customThreshold=useThresh, daily=daily, longTerm=longTerm, yearmonth=yearmonth,
                               mergelongterm=mergelongterm, withhistory=withhistory, configfile=configfile, stopmonth=stopmonth,
                               stopday=stopday, forcealert=forcealert, log=log)
    auto.readSourceList(mysrc)

    if longTerm is True and mergelongterm is True:

        # Remove all the old merged file for this source, before reprocessing the merged data
        if not auto.daily:
            for file in glob.glob(auto.workDir + '/' + auto.src + '*'):
                os.remove(file)
        if auto.daily:
            for file in glob.glob(auto.workDir + '/' + auto.src + '*daily*'):
                os.remove(file)

        # TO BE CHANGED !!!
        startyearmonth = '200808'
        # Hardcoded !!! Beurk, not good, ugly, bad !!!

        if auto.stopmonth is not None:
            thisyearmonth = auto.stopmonth
        else:
            thisyearmonth = datetime.date.today().strftime('%Y%m')
        thisyear = thisyearmonth[:-2]
        thismonth = thisyearmonth[-2:]
        startyear = startyearmonth[:-2]
        startmonth = startyearmonth[-2:]

        # First make sure that all the month-by-month long-term data have been processed
        #
        # Loop on month from 2008/08 to this month
        for year in range(int(startyear), int(thisyear) + 1):
            for month in range(1, 12 + 1):
                # To retrieve the correct results directories, 'month' should be made of 2 digits
                month = '%02d' % month
                tmpyearmonth = str(year) + str(month)
                if (year == int(startyear) and int(month) < int(startmonth)) or (
                        year == int(thisyear) and int(month) > int(thismonth)):
                    continue

                # If year=thisyear and month=thismonth, we should remove all data for this source and reprocess everything again with fresh, brand new data !
                # BUT only if update=True
                if year == int(thisyear) and int(month) == int(thismonth) and update is True:
                    tmpworkdir = auto.baseOutDir + "/longTerm/" + str(year) + str(month)
                    if not auto.daily:
                        for file in glob.glob(tmpworkdir + '/' + auto.src + '*'):
                            os.remove(file)
                    if auto.daily:
                        for file in glob.glob(tmpworkdir + '/' + auto.src + '*daily*'):
                            os.remove(file)

                processSrc(mysrc=auto.src, useThresh=useThresh, daily=auto.daily, mail=False, longTerm=True, test=False,
                           yearmonth=tmpyearmonth, mergelongterm=False, update=update, configfile=configfile,
                           stopmonth=stopmonth, forcealert=forcealert, too=auto.too)

        # Then merge the GTI files together, and run createXML, photoLC, exposure, createDAT, createLCfig, createEnergyTimeFig. No mail is sent here.
        auto.mergeGTIfiles()
        if auto.fglName is not None:
            auto.createXML()
            mygamma = None
        else:
            mygamma = ASSUMEDGAMMA
            logging.info('[{src:s}] \033[93mNo 4FGL counterpart given in the list of sources, assuming photon index of {gamma:.2f} for the light curve generation.\033[0m'.format(src=auto.src, gamma=mygamma))
        auto.photoLC()
        auto.exposure(gamma=mygamma)
        auto.createDAT()
        auto.createLCfig()
        auto.createEnergyTimeFig()
        # Exit here
        return False
    # End mergelongterm


    # When mergelongterm is False, we do the following:
    auto.selectSrc()
    auto.makeTime()

    # If we are in --long-term mode, but not in --merge-long-term mode, we can stop here, since the --merge-long-term mode then starts at the mergeGTIfiles level
    if longTerm:
        return False

    global FLAGASSUMEDGAMMA
    if auto.fglName is not None:
        auto.createXML()
        mygamma = None
        FLAGASSUMEDGAMMA = False
    else:
        mygamma = ASSUMEDGAMMA
        logging.info('[{src:s}] \033[93mNo 4FGL counterpart given in the list of sources, assuming photon index of {gamma:.2f} for the light curve generation.\033[0m'.format(src=auto.src, gamma=mygamma))
        FLAGASSUMEDGAMMA = True
    try:
        auto.photoLC()
        auto.exposure(gamma=mygamma)
        auto.createDAT()
        auto.createLCfig()
        auto.createEnergyTimeFig()
        alertSent = auto.sendAlert(nomailall=test, sendmail=mail)
    except KeyError:
        return False, False

    return auto.active, auto.visible


class automaticLightCurve:
    """
    FLaapLUC

    Automatic aperture photometry light curve generation.
    Main class, for a given of source.
    """

    def __init__(self, file=None, customThreshold=False, daily=False,
                 longTerm=False, yearmonth=None, mergelongterm=False,
                 withhistory=False, stopmonth=None, stopday=None,
                 configfile='default.cfg', forcealert=False,
                 log=logging.INFO):

        self.config = self.getConfig(configfile=configfile)
        self.allskyDir = self.config.get('InputDirs', 'AllskyDir')
        self.archiveDir = self.config.get('InputDirs', 'ArchiveDir')
        self.templatesDir = self.config.get('InputDirs', 'TemplatesDir')
        self.catalogFile = self.config.get('InputFiles', 'CatalogFile')
        if file is None:
            self.file = self.config.get('InputFiles', 'SourceList')
        else:
            self.file = file
        self.baseOutDir = self.config.get('OutputDirs', 'OutputResultsDir')
        self.allskyFile = self.allskyDir + "/" + self.config.get('InputFiles', 'WholeAllskyFile')
        self.lastAllskyFile = self.allskyDir + "/" + self.config.get('InputFiles', 'LastAllskyFile')
        self.spacecraftFile = self.allskyDir + "/" + self.config.get('InputFiles', 'SpacecraftFile')
        self.forcealert = forcealert

        try:
            self.longtimebin = float(self.config.get('AlertTrigger', 'LongTimeBin'))
        except:
            # Take 7 days by default
            self.longtimebin = 7.
            logging.warning('\033[93mCan not read LongTimeBin in config file, taking {lt:.1f} days as default.\033[0m'.format(lt=self.longtimebin))

        try:
            self.sigma = float(self.config.get('AlertTrigger', 'Sigma'))
        except:
            # Take 2 sigma by default
            self.sigma = 3.
            logging.warning('\033[93mCan not read Sigma in config file, taking {sig:.1f} as default.\033[0m'.format(sig=self.sigma))

        try:
            self.sigmaLT = float(self.config.get('AlertTrigger', 'SigmaLT'))
        except:
            # Take 2 sigma by default
            self.sigmaLT = 1.5
            logging.warning('\033[93mCan not read SigmaLT in config file, taking {siglt:.1f} as default.\033[0m'.format(siglt=self.sigmaLT))

        try:
            self.lowerExposure = float(self.config.get('AlertTrigger', 'LowerExposure'))
        except:
            self.lowerExposure = 0.

        # Read maxz and maxZA as lists, not as single floats
        self.maxz = [float(i) for i in getConfigList(self.config.get('AlertTrigger', 'MaxZ'))]
        self.maxZA = [float(i) for i in getConfigList(self.config.get('AlertTrigger', 'MaxZA'))]
        try:
            self.checkVisibility = self.config.get('Site', 'CheckVisibility')
            self.siteLon = float(self.config.get('Site', 'SiteLongitude'))
            self.siteLat = float(self.config.get('Site', 'SiteLatitude'))
            self.siteAlt = float(self.config.get('Site', 'SiteAltitude'))
        except:
            # Don't check the source visibility, by default
            self.checkVisibility = False

        # Is the source currently observable from the site under reasonable Moonlight ?
        self.observableUnderMoonlight = False
        self.visibleUnderMoonlight = False
        self.visibleUnderMoonlightOnly = False
        self.visibleInDarkness = False
        self.visibilityFig = None

        self.daily = daily
        self.withhistory = withhistory

        # Mail sender and recipients
        self.usualRecipients = getConfigList(self.config.get('MailConfig', 'UsualRecipients'))
        self.testRecipients = getConfigList(self.config.get('MailConfig', 'TestRecipients'))
        self.mailSender = self.config.get('MailConfig', 'MailSender')

        today = datetime.date.today().strftime('%Y%m%d')
        self.stopday = stopday
        if self.stopday is not None:
            today = self.stopday.replace('-', '')
            self.lastAllskyFile = self.allskyFile

        # Setting file names and directories
        if longTerm:
            self.allsky = self.allskyFile
            if not mergelongterm:
                self.workDir = self.baseOutDir + "/longTerm/" + yearmonth
            else:
                self.workDir = self.baseOutDir + "/longTerm/merged"
        else:
            self.allsky = self.lastAllskyFile
            self.workDir = self.baseOutDir + "/" + today

        self.spacecraft = self.spacecraftFile
        if not os.path.isdir(self.workDir):
            try:
                os.makedirs(self.workDir)
            except OSError:
                pass

        self.fermiDir = os.getenv('FERMI_DIR')

        # Setting default parameters
        self.roi = 1.  # degrees (http://fermi.gsfc.nasa.gov/ssc/data/analysis/scitools/aperture_photometry.html: "For aperture photometry we select a very small aperture (rad=1 degree), because we are not fitting the background.")
        try:
            self.emin = float(self.config.get('Erange', 'Emin'))
        except:
            # Take 100 MeV by default
            self.emin = 1.e2  # E min
            logging.warning('\033[93mCan not read Emin in config file, taking {emin:.1g} MeV as default.\033[0m'.format(emin=self.emin))
        try:
            self.emax = float(self.config.get('Erange', 'Emax'))
        except:
            # Take 500 GeV by default
            self.emax = 5.e5  # E max
            logging.warning('\033[93mCan not read Emax in config file, taking {emax:.1g} MeV as default.\033[0m'.format(emax=self.emax))
        self.zmax = 90.  # degrees
        self.rockangle = 52.  # maximal allowed rocking angle

        if self.daily:
            self.tbin = 24. * 60. * 60.  # seconds, daily bins
        else:
            self.tbin = self.longtimebin * 24. * 60. * 60.  # seconds, longtimebin by defaults

        self.threshold = 1.e-6  # ph cm^-2 s^-1
        self.customThreshold = customThreshold

        self.stopmonth = stopmonth

        # Open allsky file to get the start and stop dates
        try:
            hdu = fits.open(self.allsky)
        except IOError as e:
            logging.error("""I/O error ({errno}): can not open file {f}: {err}
I will create the allsky file on the fly for you, for the last month of available data, using enrico.
First, retrieving the last photon files...
            """.format(errno=e.errno, f=self.allsky, err=e.strerror))
            cmd = 'enrico_download --download_data'
            r = os.system(cmd)
            assert (r == 0), "Could not properly download the last data."
            logging.error('Second, retrieving the last spacecraft file...')
            cmd = 'enrico_download --download_spacecraft'
            r = os.system(cmd)
            assert (r == 0), "Could not properly download the last spacecraft file."
            logging.error('Third, creating the allsky file with enrico...')
            cmd = 'enrico_download --preprocess_data --steps=gtselect --event_classes=source --selections=all --emins=100'
            r = os.system(cmd)
            assert (r == 0), "Could not properly generate the allsky file."
            import enrico.data
            self.allsky = enrico.data.PREPROCESSED_DIR + '/source/all/emin_000100/gtselect.fits'
            hdu = fits.open(self.allsky)
        header = hdu[0].header

        self.filtermoon = self.hasmoon()
        
        if not longTerm:
            self.tstart = header['TSTART']
            self.tstop = header['TSTOP']
            if self.stopday is not None:
                logging.debug('Input stop day is {}'.format(self.stopday))
                logging.debug('i.e. MJD {}'.format(t.Time('%s 00:00:00' % self.stopday, format='iso', scale='utc').mjd))
                stopdaymet = extras.mjd2met(t.Time('%s 00:00:00' % self.stopday, format='iso', scale='utc').mjd)
                logging.debug('i.e. {} MET'.format(stopdaymet))
                logging.debug('while tstop is {}'.format(self.tstop))
                if stopdaymet > self.tstop:
                    logging.error('FLaapLUC can not predict the future !')
                    sys.exit(1)
                self.tstop = stopdaymet
                self.tstart = self.tstop - 30 * 24 * 3600  # stop - 30 days
        else:
            missionStart = header['TSTART']  # in MET
            missionStop = header['TSTOP']  # in MET

            if not mergelongterm:
                # Need to convert 'yearmonth' in MET
                # self.tstart is the first day of yearmonth at 00:00:00, or missionStart
                # self.tstop  is the first day of next month at 00:00:00, or missionStop
                year = yearmonth[:-2]
                month = yearmonth[-2:]

                # Get date of first day of yearmonth at 00:00:00, in UNIX time (timetuple transform a datetime object in time object ???)
                #                                                        year         month   day   hour   minute  second  microsecond
                yearmonthStart = time.mktime(datetime.datetime(int(year), int(month), 1, 0, 0, 0, 0).timetuple())
                if int(month) < 12:
                    yearmonthStop = time.mktime(datetime.datetime(int(year), int(month) + 1, 1, 0, 0, 0, 0).timetuple())
                else:
                    yearmonthStop = time.mktime(datetime.datetime(int(year) + 1, 1, 1, 0, 0, 0, 0).timetuple())

                # Convert these from UNIX time to MET
                tmptstart = extras.mjd2met(extras.unixtime2mjd(yearmonthStart))
                tmptstop = extras.mjd2met(extras.unixtime2mjd(yearmonthStop))

                logging.debug('INIT yearmonthStart={}'.format(yearmonthStart))
                logging.debug('INIT yearmonthStop={}'.format(yearmonthStop))

                # Make sure that start of yearmonth is after the launch of Fermi, and that stop of yearmonth is before the very last data we have from NASA servers !
                if tmptstart > missionStart:
                    self.tstart = tmptstart
                else:
                    self.tstart = missionStart

                if tmptstop < missionStop:
                    self.tstop = tmptstop
                else:
                    self.tstop = missionStop

            if mergelongterm is True:
                self.tstart = missionStart
                self.tstop = missionStop

    def getConfig(self, configfile='./default.cfg'):
        """Get configuration from a configuration file."""
        self.config = ConfigParser()
        self.config.readfp(open(configfile))
        return self.config

    def readSourceList(self, mysrc=None):
        """
        Read the list of sources.

        @todo Use a mySQL database instead of an ASCII file for the list of sources ?
        """

        try:
            srcList = ascii.read(self.file)
        except IOError:
            logging.error('Can not open {f:s}'.format(f=self.file))
            sys.exit(1)

        src = srcList['Name']
        ra = srcList['RA']
        dec = srcList['Dec']
        z = srcList['z']
        fglName = srcList['4FGLname']
        # Read the threshold for the source from the source list, if we asked to process with custom thresholds when instanciating the class
        myThreshold = srcList['Threshold']
        too = srcList['ToO']

        # If we ask for a particular source, return the parameters for that source
        if mysrc != None:
            # Find our input src in the list of sources
            for i in range(len(src)):
                if src[i] == mysrc:
                    # Redefine the threshold if we provided a custom threshold
                    if self.customThreshold and myThreshold[i] != 0.:
                        try:
                            self.threshold = float(myThreshold[i])
                            logging.debug('[{src:s}] Read threshold={th} from source list'.format(src=mysrc,th=self.threshold))
                        except ValueError:
                            logging.warning('The threshold of the source {src:s} is ill-defined (not a float). Please, check the list of sources !'.format(src=mysrc))
                            sys.exit(2)
                    self.src = src[i]
                    self.ra = ra[i]
                    self.dec = dec[i]
                    self.z = z[i]
                    self.fglName = fglName[i]
                    if self.fglName == 'None':
                        self.fglName = None
                    self.too = too[i]
                    return

            logging.warning('\033[92mCan\'t find your source {src:s} in the list of sources !\033[0m'.format(src=mysrc))
            self.src = None
            self.ra = None
            self.dec = None
            self.z = None
            self.fglName = None
            self.too = None
            return

        # Otherwise, return the whole list of parameters for all the sources
        else:
            return src, ra, dec, z, fglName, too

    def selectSrc(self):
        """
        Filter a given source, running gtselect
        """
        # Do we have to deal with a FITS file or an ASCII list of FITS file ?
        allskyext = os.path.splitext(self.allsky)[1]
        if allskyext in [".fit", ".fits"]:
            fermi.filter['infile'] = self.allsky
        else:
            fermi.filter['infile'] = '@%s' % self.allsky
        if self.daily:
            outfile = self.workDir + '/' + str(self.src) + '_daily.fits'
        else:
            outfile = self.workDir + '/' + str(self.src) + '.fits'
        fermi.filter['outfile'] = outfile

        # If outfile already exists, we don't do anything
        if os.path.isfile(outfile):
            return True

        fermi.filter['ra'] = self.ra
        fermi.filter['dec'] = self.dec
        fermi.filter['rad'] = self.roi
        fermi.filter['emin'] = self.emin
        fermi.filter['emax'] = self.emax
        fermi.filter['tmin'] = self.tstart
        fermi.filter['tmax'] = self.tstop
        fermi.filter['zmax'] = self.zmax
        fermi.filter['evclass'] = 128
        logging.info('Running gtselect')
        fermi.filter.run()

    def hasmoon(self):
        """
        Return True if the FT2 spacecraft file has the Moon coordinates (using the moonpos user contributed script), False otherwise
        """
        try:
            d = fits.open(self.spacecraft)[1].data
        except:
            liste = np.loadtxt(self.spacecraft,unpack=True,dtype=str)
            d = fits.open(liste[0])[1].data
        try:
            m = d.field('RA_MOON')
            r = True
        except KeyError:
            r = False
        return r

    def makeTime(self):
        """
        Filter the GTI for a given source
        """
        fermi.maketime['scfile'] = self.spacecraft

        if self.daily:
            fermi.maketime['evfile'] = self.workDir + '/' + str(self.src) + '_daily.fits'
            outfile = self.workDir + '/' + str(self.src) + '_daily_gti.fits'
        else:
            fermi.maketime['evfile'] = self.workDir + '/' + str(self.src) + '.fits'
            outfile = self.workDir + '/' + str(self.src) + '_gti.fits'
        fermi.maketime['outfile'] = outfile

        # If outfile already exists, we don't do anything
        if os.path.isfile(outfile):
            return True

        # cf. http://fermi.gsfc.nasa.gov/ssc/data/analysis/scitools/aperture_photometry.html
        gtifilter = "LAT_CONFIG==1 && DATA_QUAL>0 && (angsep(" + str(self.ra) + "," + str(
            self.dec) + ",RA_SUN,DEC_SUN)>5.)"
        if self.filtermoon:
            gtifilter += " && (angsep(" + str(self.ra) + "," + str(
                self.dec) + ",RA_MOON,DEC_MOON)>5.)"
        fermi.maketime['filter'] = gtifilter
        fermi.maketime['roicut'] = 'no'
        fermi.maketime['tstart'] = self.tstart
        fermi.maketime['tstop'] = self.tstop
        logging.info('Running gtmktime')
        fermi.maketime.run()

    def mergeGTIfiles(self):
        """
        Merge multiple GTI files when mergelongterm is True.
        Use gtselect.
        Assume the current workDir is longTerm/merged.
        """

        # Create list of GTI files
        if not self.daily:
            listname = self.workDir + '/' + self.src + '_gti.list'
        else:
            listname = self.workDir + '/' + self.src + '_daily_gti.list'
        filelist = open(listname, 'w')
        list = []
        if not self.daily:
            for file in glob.glob(self.workDir + '/../20????/' + self.src + '_gti.fits'):
                list.append(file)
        else:
            for file in glob.glob(self.workDir + '/../20????/' + self.src + '_daily_gti.fits'):
                list.append(file)
        # Sort the list of GTI files
        list = sorted(list)
        for item in list:
            filelist.write(item + '\n')
        filelist.close()

        fermi.filter['infile'] = '@' + listname
        if not self.daily:
            outfile = self.workDir + '/' + str(self.src) + '_gti.fits'
        else:
            outfile = self.workDir + '/' + str(self.src) + '_daily_gti.fits'
        fermi.filter['outfile'] = outfile

        # If outfile already exists, we re-create it
        if os.path.isfile(outfile):
            os.remove(outfile)

        fermi.filter['ra'] = self.ra
        fermi.filter['dec'] = self.dec
        fermi.filter['rad'] = self.roi
        fermi.filter['emin'] = self.emin
        fermi.filter['emax'] = self.emax
        fermi.filter['tmin'] = self.tstart
        fermi.filter['tmax'] = self.tstop
        fermi.filter['zmax'] = self.zmax
        fermi.filter['evclass'] = 128
        logging.info('Running gtselect')
        fermi.filter.run()

    def createXML(self):
        """
        Create an XML model file based on the 4FGL catalogue
        """

        if self.daily:
            evfile = self.workDir + '/' + str(self.src) + '_daily_gti.fits'
            modelfile = self.workDir + '/' + str(self.src) + '_daily.xml'
        else:
            evfile = self.workDir + '/' + str(self.src) + '_gti.fits'
            modelfile = self.workDir + '/' + str(self.src) + '.xml'

        # If modelfile already exists, we don't do anything
        if os.path.isfile(modelfile):
            return True

        import make4FGLxml

        mymodel = make4FGLxml.srcList(self.catalogFile, evfile, modelfile)
        logging.info('Running makeModel')
        mymodel.makeModel(GDfile=self.fermiDir + '/refdata/fermi/galdiffuse/gll_iem_v07.fits', GDname='GalDiffuse',
                          ISOfile=self.fermiDir + '/refdata/fermi/galdiffuse/iso_P8R3_SOURCE_V2.txt',
                          ISOname='IsotropicDiffuse', extDir=self.templatesDir, makeRegion=False)

    def photoLC(self):
        """
        Compute the photometric light curve for a given source
        """

        if self.daily:
            fermi.evtbin['evfile'] = self.workDir + '/' + str(self.src) + '_daily_gti.fits'
            outfile = self.workDir + '/' + str(self.src) + '_daily_lc.fits'
        else:
            fermi.evtbin['evfile'] = self.workDir + '/' + str(self.src) + '_gti.fits'
            outfile = self.workDir + '/' + str(self.src) + '_lc.fits'

        # If outfile already exists, we don't do anything
        if os.path.isfile(outfile):
            return True

        fermi.evtbin['outfile'] = outfile
        fermi.evtbin['scfile'] = self.spacecraft
        fermi.evtbin['algorithm'] = 'LC'
        fermi.evtbin['tbinalg'] = 'LIN'
        fermi.evtbin['tstart'] = self.tstart
        fermi.evtbin['tstop'] = self.tstop
        fermi.evtbin['dtime'] = self.tbin
        logging.info('Running gtbin')
        try:
            fermi.evtbin.run()
        except:
            fermi.evtbin['scfile'] = '@'+self.spacecraft   
            fermi.evtbin.run()

    def exposure(self, gamma=None):
        """
        Compute exposure on source src, to add a flux column for the photometric light curve.

        Warning: the input file is modified in place, with an additional exposure column added to the file !
        """

        if self.daily:
            infile = self.workDir + '/' + str(self.src) + '_daily_lc.fits'
            srcmdl = self.workDir + '/' + str(self.src) + '_daily.xml'
        else:
            infile = self.workDir + '/' + str(self.src) + '_lc.fits'
            srcmdl = self.workDir + '/' + str(self.src) + '.xml'

        # If infile already contains an EXPOSURE column, we don't do anything
        hdu = fits.open(infile)
        if hdu[1].header.get('TTYPE5') == 'EXPOSURE':
            return True

        scfile = self.spacecraft
        irfs = 'P8R3_SOURCE_V2'
        rad = str(self.roi)

        options = 'infile=' + infile + ' scfile=' + scfile + ' irfs=' + irfs + ' rad=' + rad
        if self.fglName is not None:
            target = self.fglName.replace('4FGLJ', '4FGL J')
            # Special case for IC component of Crab, which is extended, and which name is replaced by "Crab IC" by make4FGLxml:
            if self.fglName == '4FGLJ0534.5+2201i':
                target = 'Crab IC'
            logging.debug('exposure: target=%s', target)
            options += ' srcmdl=' + srcmdl + ' target="' + target + '"'
        else:
            options += ' srcmdl="none" specin=' + str(gamma)
        cmd = 'time -p gtexposure ' + options
        logging.info('Running gtexposure')
        logging.debug('Running gtexposure with command: {cmd}'.format(cmd=cmd))
        os.system(cmd)

    def createDAT(self):
        """
        Create a data file with the light curve of a given source.
        """

        # Read LC file
        if self.daily:
            infile = self.workDir + '/' + str(self.src) + '_daily_lc.fits'
            outfile = self.workDir + '/' + str(self.src) + '_daily_lc.dat'
        else:
            infile = self.workDir + '/' + str(self.src) + '_lc.fits'
            outfile = self.workDir + '/' + str(self.src) + '_lc.dat'

        # If outfile already exists, we don't do anything
        if os.path.isfile(outfile):
            return True

        try:
            hdu = fits.open(infile)
        except:
            logging.critical('Exception: can not open file %s', infile)
            raise
        data = hdu[1].data

        file = open(outfile, 'w')
        file.write("MET\tMJD\tFlux\tFluxError\n")
        file.write("#[MET]\t[MJD]\t[ph cm^-2 s^-1]\t[ph cm^-2 s^-1]\n")
        time = data.field('TIME')  # MET
        counts = data.field('COUNTS')
        countsErr = data.field('ERROR')  # error on counts
        exposure = data.field('EXPOSURE')  # cm^2 s^1
        flux = counts / exposure  # approximate flux in ph cm^-2 s^-1
        fluxErr = countsErr / exposure  # approximate flux error in ph cm^-2 s^-1

        timeMjd = extras.met2mjd(time)
        # We can do this because time is NOT a list, but a numpy.array

        for i in range(len(time)):
            # Exposure can be 0 if longTerm=True and TSTOP in photon file > TSTOP in spacecraft file, or if Fermi operated in pointed mode for a while.
            if exposure[i] != 0.:
                file.write(str(time[i]) + "\t" + str(timeMjd[i]) + "\t" + str(flux[i]) + "\t" + str(fluxErr[i]) + "\n")
        file.close()

    def getBAT(self):
        import urllib2

        # daily fits example url:
        # http://swift.gsfc.nasa.gov/docs/swift/results/transients/CygX-3.lc.fits

        # Some sources need replacement names to match the BAT names
        urls = {
            '4U1907+09': 'H1907+097',
            '1FGLJ1018.6-5856': '0FGLJ1018.2-5858',
            'H1743-322': 'IGRJ17464-3213',
            'V4641Sgr': 'SAXJ1819.3-2525',
            '1E1841-04.5': 'Kes73',
        }

        # Remove '+', add file ending
        if urls.has_key(self.src):
            file = urls[self.src].replace('+', 'p') + ".lc.fits"
        else:
            file = self.src.replace('+', 'p') + ".lc.fits"
        urlprefix = "http://swift.gsfc.nasa.gov/docs/swift/results/transients/"

        # lc files can be in a weak/ subdir for weak sources, we try both
        try:
            baturl = urlprefix + file
            webfile = urllib2.urlopen(baturl)
        except (urllib2.HTTPError, urllib2.URLError) as e:
            try:
                baturl = urlprefix + 'weak/' + file
                webfile = urllib2.urlopen(baturl)
            except (urllib2.HTTPError, urllib2.URLError) as e:
                return False, None
            except socket.error:
                return False, None
        except socket.error:
            return False, None

        # save lc to local file
        localfile = open(file, 'w')
        localfile.write(webfile.read())
        webfile.close()
        localfile.close()
        # read local file with fits into batlc
        batfits = fits.open(file)
        batlc = np.array(batfits[1].data)
        batfits.close()
        # delete local file
        os.unlink(file)

        return True, batlc

    def createLCfig(self):
        """
        Create a PNG figure with the light curve of a given source. Any existing PNG file is overwritten !
        """

        # Read the .dat LC file
        if self.daily:
            infile = self.workDir + '/' + str(self.src) + '_daily_lc.dat'
            outfig = self.workDir + '/' + str(self.src) + '_daily_lc.png'
            infileLongTimeBin = self.workDir + '/' + str(self.src) + '_lc.dat'
            duration = 1.  # duration of a time bin, in days
        else:
            infile = self.workDir + '/' + str(self.src) + '_lc.dat'
            outfig = self.workDir + '/' + str(self.src) + '_lc.png'
            duration = self.longtimebin  # duration of a time bin, in days

        data = ascii.read(infile)
        # the times are already read as MJD, cf createDAT function.
        timelc = data['MJD']
        flux = data['Flux']
        fluxErr = data['FluxError']

        if self.daily:
            dataLongTimeBin = ascii.read(infileLongTimeBin)
            # the times are already read as MJD, cf createDAT function.
            timeLongTimeBin = dataLongTimeBin['MJD']
            fluxLongTimeBin = dataLongTimeBin['Flux']
            fluxErrLongTimeBin = dataLongTimeBin['FluxError']
            durationLongTimeBin = self.longtimebin  # duration of a time bin, in days

        # Download Swift/BAT data if available
        # xray is boolean flag indicating that X-ray BAT data is available
        xray, batlc = self.getBAT()

        # Redefine the trigger threshold if withhistory=True
        if self.withhistory:
            self.dynamicalTrigger()
        else:
            self.LTfluxAverage = self.threshold
            self.LTfluxRMS = 0.

        fig = plt.figure()

        if xray:
            ax = fig.add_subplot(211)
            axbat = fig.add_subplot(212, sharex=ax)
        else:
            ax = fig.add_subplot(111)

        if self.fglName is not None:
            title = str(self.src) + ', ' + str(self.fglName).replace('_2FGLJ', '2FGL J').replace('3FGLJ', '3FGL J').replace('4FGLJ', '4FGL J')
        else:
            title = str(self.src) + ', no known 4FGL counterpart'
        if self.z == 'None':
            title = title + ' (z unknown)'
        else:
            title = title + ' (z=' + str(self.z) + ')'

        ax.set_title(title)

        # Force the y-axis ticks to use 1e-6 as a base exponent
        ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: ('%.2f') % (x * 1e6)))
        ax.set_ylabel('F (%.0f MeV-%.0f GeV) (%s 10$^{-6}$ ph cm$^{-2}$ s$^{-1}$)' % (
        self.emin, self.emax / 1000., r'$\times$'))  # , size='x-small')

        ## Make the x-axis ticks shifted by some value
        ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: '%.0f' % (x - TOFFSET)))
        ax.set_xlabel('MJD-' + str(TOFFSET))
        # ax.set_xlabel('MJD')

        # Plot the Fermi/LAT light curve
        if self.daily:
            # Also plot the long time-binned light curve
            ax.errorbar(x=timelc, xerr=duration / 2., y=flux, yerr=fluxErr, fmt='ro')
            ax.errorbar(x=timeLongTimeBin, xerr=durationLongTimeBin / 2., y=fluxLongTimeBin, yerr=fluxErrLongTimeBin,
                        fmt='bo')
            # The last plot called is on top of the others in matplotlib (are you sure ???). Here, we want the long time-binned LC on top, for visibility.
        else:
            ax.errorbar(x=timelc, xerr=duration / 2., y=flux, yerr=fluxErr, fmt='bo')

        # Plot a line at the threshold value
        ax.axhline(y=self.threshold, linewidth=3, linestyle='--', color='r')
        if self.withhistory:
            ax.axhline(y=self.LTfluxAverage, linewidth=1, linestyle='-', color='b')
            ax.axhline(y=self.LTfluxAverage + self.LTfluxRMS, linewidth=1, linestyle='--', color='b')
            ax.axhline(y=self.LTfluxAverage - self.LTfluxRMS, linewidth=1, linestyle='--', color='b')

        # Plot a line at flux=0, for visibility/readibility
        ax.axhline(y=0., color='k')

        # Add a label for the creation date of this figure
        # x,y in relative 0-1 coords in figure
        plt.figtext(0.98, 0.95,
                    'plot creation date: %s (UTC)' % (time.strftime("%a, %d %b %Y %H:%M:%S", time.gmtime())),
                    horizontalalignment="right",
                    rotation='vertical',
                    size='xx-small'
                    )

        # Plot Swift/BAT lightcurve
        if xray:
            axbat.errorbar(batlc['TIME'] + 0.5, batlc['RATE'], batlc['ERROR'], fmt="none", capsize=0, elinewidth=1,
                           ecolor='b', color='b')
            axbat.set_xlabel('MJD-' + str(TOFFSET))
            # axbat.set_xlabel('MJD')
            axbat.set_ylabel('F (15-50 keV) (count cm^-2 s^-1)', size='x-small')
            try:
                axbat.set_xlim(xmin=timelc[0] - duration / 2. - 1., xmax=timelc[-1:] + duration / 2. + 1.)
                axbat.set_ylim(ymin=0.)
            except:
                pass

        # Need to zoom in or not, at the very end, after any call to other matplotlib functions
        NEEDTOZOOMIN = False
        for i in range(len(flux)):
            if fluxErr[i] > 5. * flux[i]:
                NEEDTOZOOMIN = True
        if NEEDTOZOOMIN:
            maxy = 1.5 * max(flux)
            if maxy > self.threshold:
                ax.set_ylim(ymin=-1.e-7, ymax=maxy)
            else:
                ax.set_ylim(ymin=-1.e-7, ymax=self.threshold)

        # Don't show the figure in batch mode
        if not BATCH:
            plt.show()
        # Save the figure
        fig.savefig(outfig)
        plt.close(fig)

    def createEnergyTimeFig(self, eThresh=1.e2):
        """
        Create a PNG figure with the energy vs time of a given source, above eThresh MeV. Any existing PNG file is overwritten !
        """

        # Read the GTI FITS file
        infile = self.workDir + '/' + str(self.src) + '_gti.fits'
        outfig = self.workDir + '/' + str(self.src) + '_energyTime.png'

        hdu = fits.open(infile)
        data = hdu[1].data
        mask = data.field('ENERGY') > eThresh
        datac = data[mask]
        if not datac.size:
            logging.warning('[{src:s}] \033[92mEmpty energy vs time plot above {eth:0.2f} GeV\033[0m'.format(src=self.src, eth=eThresh/1.e3))
            return

        t = extras.met2mjd(datac['TIME'])
        e = datac['ENERGY']

        fig = plt.figure()
        ax = fig.add_subplot(111)

        if self.fglName is not None:
            title = str(self.src) + ', ' + str(self.fglName).replace('_2FGLJ', '2FGL J').replace('3FGLJ', '3FGL J').replace('4FGLJ', '4FGL J')
        else:
            title = str(self.src) + ', no known 4FGL counterpart'
        if self.z == 'None':
            title = title + ' (z unknown)'
        else:
            title = title + ' (z=' + str(self.z) + ')'

        ax.set_title(title)

        ylabel = 'Energy (MeV)'
        if eThresh > self.emin:
            ylabel += ' -- only data above %.1f GeV are shown' % (eThresh / 1.e3)
        # ax.set_ylabel(ylabel, size='x-small')
        ax.set_ylabel(ylabel)

        ## Make the x-axis ticks shifted by some value
        ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: '%.0f' % (x - TOFFSET)))
        ax.set_xlabel('MJD-' + str(TOFFSET))
        try:
            ax.set_xlim(xmin=t[0] - 1., xmax=t[-1:] + 1.)
        except:
            pass

        # Plot the energy vs time distribution
        try:
            # cf. http://stackoverflow.com/questions/20105364/how-can-i-make-a-scatter-plot-colored-by-density-in-matplotlib
            from scipy.stats import gaussian_kde
            xy = np.vstack([t, e])
            z = gaussian_kde(xy)(xy)
            # Re-normalize the density
            z = z / max(z)
            idx = z.argsort()
            t, e, z = t[idx], e[idx], z[idx]
            pcm = ax.scatter(t, e, c=z, s=100, edgecolor='')
            cbar = plt.colorbar(pcm, ax=ax)
            cbar.set_label('Kernel-density estimates (arb. unit)', rotation=90)
        except ImportError:
            ax.plot(t, e, 'bo')
        ax.set_yscale('log')

        # Add a label for the creation date of this figure
        # x,y in relative 0-1 coords in figure
        plt.figtext(0.98, 0.95,
                    'plot creation date: %s (UTC)' % (time.strftime("%a, %d %b %Y %H:%M:%S", time.gmtime())),
                    horizontalalignment="right",
                    rotation='vertical',
                    size='xx-small'
                    )

        # Don't show the figure in batch mode
        if not BATCH:
            plt.show()
        # Save the figure
        fig.savefig(outfig)
        plt.close(fig)

    def createVisibilityPlot(self, begindate):
        '''
        Create a PNG figure with the visibility plot for the source from the site.
        '''

        logging.info('[{src:s}] Creating visibility plot'.format(src=self.src))
        logging.debug('[createVisibiliytPlot] begin date is {}'.format(begindate))

        fig = plt.figure(figsize=(10, 5))
        ax = plt.subplot(1, 1, 1)
        ax.xaxis_date()
        ax.xaxis.set_minor_formatter(DateFormatter("%Hh"))
        ax.xaxis.set_minor_locator(HourLocator(interval=2))
        ax.xaxis.set_major_formatter(DateFormatter("%Y-%m-%d     "))
        ax.xaxis.set_major_locator(DayLocator())
        ax.tick_params(pad=20)
        if self.fglName is not None:
            title = str(self.src) + ', ' + str(self.fglName).replace('_2FGLJ', '2FGL J').replace('3FGLJ', '3FGL J').replace('4FGLJ', '4FGL J')
        else:
            title = str(self.src) + ', no known 4FGL counterpart'
        if self.z == 'None':
            title += ' (z unknown)'
        else:
            title += ' (z=' + str(self.z) + ')'
        ax.set_title(title)

        self.visibilityFig = self.workDir + '/' + str(self.src) + '_visibilityPlot.png'

        site = c.EarthLocation(lon=self.siteLon * u.deg, lat=self.siteLat * u.deg,
                               height=self.siteAlt * u.m)
        # Time range for plot
        hrange = 10.
        centerhour = 24 - int(round(site.longitude.hour))
        begindate = begindate.datetime()  # ephem.Date to datetime.datetime
        begindate = t.Time(datetime.datetime(begindate.year, begindate.month, begindate.day, centerhour),
                           scale='utc')
        dt = np.linspace(-int(hrange), int(hrange),
                         100 * int(hrange)) * u.hour
        times = begindate + dt
        fdates = np.array([date2num(x.datetime) for x in times])

        # Sun, Moon and source altitudes
        horizframe = c.AltAz(obstime=times, location=site)
        sun = c.get_sun(times).transform_to(horizframe)

        srcCoords = c.SkyCoord(ra=self.ra * u.degree, dec=self.dec * u.degree, frame='icrs')
        altaz = srcCoords.transform_to(horizframe)

        obs = ephem.Observer()
        moonHorizon = ephem.degrees(str(MOON_TWILIGHT_ANGLE.value))  # Moon twilight angle to -0.83°
        obs.lon = ephem.degrees(str(self.siteLon))  # ephem needs this as string
        obs.lat = ephem.degrees(str(self.siteLat))  # ephem needs this as string
        obs.elev = self.siteAlt
        obs.compute_pressure()

        moon = ephem.Moon()
        moon_alts = np.zeros_like(times)

        for ii, tt in enumerate(times):
            obs.date = ephem.Date(tt.datetime)
            moon.compute(obs)
            moon_alts[ii] = c.Angle(moon.alt, unit=u.rad).degree
            # logging.debug('[createVisibilityPlot] It is {date}, {src} alt. {srcalt:.1f}, Sun alt. {sun:.1f}, Moon alt. {moon:.1f} deg.'.format(date=obs.date, src=self.src, srcalt=altaz.alt[ii], sun=sun.alt[ii], moon=moon_alts[ii]))

        moon_alts = moon_alts * u.deg

        # plt.fill_between(fdates, 0, 90, sun.alt < -0 * u.deg,
        #                  color='0.5', zorder=0)

        # nautical twilight
        plt.fill_between(fdates, 0, 90, sun.alt < CIVIL_TWILIGHT_ANGLE,
                         color='0.3', zorder=0, label='Twilight')

        # astronomical twilight
        darkness = (sun.alt < ASTRONOMICAL_TWILIGHT_ANGLE)
        plt.fill_between(fdates, 0, 90, darkness,
                         color='k', zorder=0, label='Darkness')

        # plot the moon window. Note that we don't use -18 degrees
        moonup = (moon_alts > MOON_TWILIGHT_ANGLE) & (sun.alt < ASTRONOMICAL_TWILIGHT_ANGLE)
        plt.fill_between(fdates, 0, 90, moonup,
                         color='sienna', zorder=0, label='Moon')

        srcObsDark = (darkness) & (~moonup) & (altaz.alt.degree > self.srcMinAlt)
        logging.debug('[{src}] Observable under dark conditions between {start} and {stop}'.format(src=self.src, start=times[srcObsDark][0], stop=times[srcObsDark][-1]))
        if self.visibleUnderMoonlight:
            srcObsMoon = (darkness) & (moonup) & (altaz.alt.degree > self.srcMinAlt)
            logging.debug('[{src}] Observable under moonlight conditions between {start} and {stop}'.format(src=self.src, start=times[srcObsMoon][0], stop=times[srcObsMoon][-1]))

        ax.axhline(y=self.srcMinAlt,
                   linestyle='--',
                   color='red',
                   linewidth=2.0,
                   label='ZA={za} deg'.format(za=(90.-self.srcMinAlt)))

        fig.autofmt_xdate(rotation=0, ha='center')
        plt.xlabel("Time ({})".format(times.scale.upper()), fontsize='small')
        plt.ylabel("Altitude (deg)", fontsize='small')
        plt.ylim(0, 90)

        plt.grid(color="0.8", which="minor")
        plt.grid(color="0.8", which="major")
        plt.grid(color="0.8", which="major", lw=4, axis="x")

        leg = plt.legend(loc="upper right", frameon=True)

        plt.plot(fdates, altaz.alt, lw=3, label=self.src,
             path_effects=[patheffects.Stroke(linewidth=4, foreground='white'),
                           patheffects.Normal()])
        for tt in leg.texts:
            tt.set_path_effects([patheffects.Stroke(linewidth=3, foreground="w"),
                                 patheffects.Normal()])

        # Add a label for the creation date of this figure
        # x,y in relative 0-1 coords in figure
        plt.figtext(0.98, 0.95,
                    'plot creation date: %s (UTC)' % (time.strftime("%a, %d %b %Y %H:%M:%S", time.gmtime())),
                    horizontalalignment="right",
                    rotation='vertical',
                    size='xx-small')

        # Add information on start and stop times of observability under dark and/or moonlight conditions
        if self.visibleInDarkness:
            obstext = '''{src} observable under dark conditions between {start} and {stop}, for ZA<{maxza} deg
'''.format(src=self.src,
           start=t.Time(times[srcObsDark][0], format='iso', scale='utc', out_subfmt='date_hm'),
           stop=t.Time(times[srcObsDark][-1], format='iso', scale='utc', out_subfmt='date_hm'),
           maxza=(90.-self.srcMinAlt))
            if self.visibleUnderMoonlight:
                obstext += '''{src} observable under reasonable moonlight conditions between {start} and {stop}
'''.format(src=self.src,
           start=t.Time(times[srcObsMoon][0], format='iso', scale='utc', out_subfmt='date_hm'),
           stop=t.Time(times[srcObsMoon][-1], format='iso', scale='utc', out_subfmt='date_hm'))
            else:
                obstext += "but not observable under reasonable moonlight conditions."
        else:
            obstext = "{src} is not observable (nor under dark or moonlight conditions) during the night".format(src=self.src)
        plt.figtext(0.05, 0.01,
                    obstext,
                    size='x-small')

        # Don't show the figure in batch mode
        if not BATCH:
            plt.show()
        # Save the figure
        fig.savefig(self.visibilityFig)
        plt.close(fig)


    def zaAtCulmination(self):
        """
        Returns the zenith angle of a source at culmination, for the provided site.
        """
        za = np.abs(self.dec - self.siteLat)
        return za

    def is_visible(self):
        '''
        Check whether the current source is visible at the site provided.
        '''

        # Define site for pyephem
        site = ephem.Observer()
        astroHorizon = ephem.degrees(str(ASTRONOMICAL_TWILIGHT_ANGLE.value))  # astronomical twilight
        moonHorizon = ephem.degrees(str(MOON_TWILIGHT_ANGLE.value))  # Moon twilight angle to -0.83°
        site.horizon = astroHorizon
        site.lon = ephem.degrees(str(self.siteLon))  # ephem needs this as string
        site.lat = ephem.degrees(str(self.siteLat))  # ephem needs this as string
        site.elev = self.siteAlt  # meters
        site.compute_pressure()

        srcCoords = c.SkyCoord(ra=self.ra * u.degree, dec=self.dec * u.degree, frame='icrs')

        # If input z is None, make it believe it is 0, otherwise msk crashes:
        if self.z == 'None':
            z = 0.
            logging.warning('[{src:s}] \033[93mUnknown z, taking z=0.\033[0m'.format(src=self.src))
        else:
            z = float(self.z)

        # We also want the max allowed ZA for the given z of the source
        maxz = np.array(self.maxz)
        maxZA = np.array(self.maxZA)
        if z > np.max(maxz):
            thismaxZA = np.min(maxZA)
            logging.warning('[{src}] \033[93mz={z} is greater than max probed z={maxz} !\033[0m'.format(src=self.src, z=z, maxz=np.max(maxz)))
        else:
            msk = np.where(z < maxz)
            # Get the first item in the mask, to get the corresponding ZA:
            thismaxZA = maxZA[msk[0][0]]
            thismaxz = maxz[msk[0][0]]

        # Convert ZA to Alt
        self.srcMinAlt = np.abs(90. - thismaxZA)

        ephemSrc = ephem.FixedBody()
        ephemSrc._ra = ephem.hours(str(srcCoords.ra.to_string(unit=u.hourangle,
                                                              sep=':')))  # Careful: ephem should be given hours here, but only for RA !
        ephemSrc._dec = ephem.degrees(str(srcCoords.dec.to_string(unit=u.degree, sep=':')))

        visibleFlag = False

        zaAtCulmin = self.zaAtCulmination()
        if zaAtCulmin > 90.:
            # the source is basically NEVER visible at the site
            logging.info('[{src:s}] \033[91mNEVER above horizon at the site, consider discarding this source from your source list...\033[0m'.format(src=self.src))
            return False

        if thismaxZA < zaAtCulmin:
            # the source is never above maxZA set by 2D mask on Dec/z
            logging.info('[{src:s}]\033[91m Never above allowed max ZA, consider relaxing the Dec/z cuts or discarding this source from your source list...\033[0m'.format(src=self.src))
            logging.debug('[{src:s}] thismaxZA={maxZA:f}, zaAtCulmin={zaculmin:f}, thismaxz={maxz:f}'.format(src=self.src, maxZA=thismaxZA, zaculmin=zaAtCulmin, maxz=thismaxz))
            return False

        # All times are handled here in UTC (pyEphem only uses UTC)
        if self.stopday is not None:
            site.date = ephem.Date(self.stopday.replace('-','/')+' 12:00:00')
        else:
            site.date = datetime.datetime.utcnow()
        sun = ephem.Sun()
        sun.compute(site)
        nextSunset = site.next_setting(sun)
        nextSunrise = site.next_rising(sun)
        # so far, so good. All of this is OK if we execute the program during day time.

        # However, if the program is run during dark time, we should look at the ephemerids of next night (not current night):
        if nextSunrise < nextSunset:
            logging.info('Looking at visibility for tomorrow')
            # we just put the current time at next sunrise + 10 min., to be sure to fall on tomorrow's morning day time
            site.date = nextSunrise.datetime() + datetime.timedelta(minutes=10)
            sun.compute(site)
            nextSunset = site.next_setting(sun)
            nextSunrise = site.next_rising(sun)

        logging.debug('[is_visible] Date is set to {date}'.format(date=site.date))

        # The Moon just needs to be below the horizon, not below astronomical twilight angle
        site.horizon = moonHorizon
        moon = ephem.Moon()
        moon.compute(site)
        nextMoonset = site.next_setting(moon)
        nextMoonrise = site.next_rising(moon)
        site.horizon = astroHorizon

        logging.debug('[is_visible] Next sun rise at {0}'.format(nextSunrise))
        logging.debug('[is_visible] Next sun set at {0}'.format(nextSunset))
        logging.debug('[is_visible] Next moon rise at {0}'.format(nextMoonrise))
        logging.debug('[is_visible] Next moon set at {0}'.format(nextMoonset))

        ephemSrc.compute(site)
        moonPhase = moon.phase  # phase of the Moon between 0 and 100
        targetMoonSeparation = c.Angle(ephem.separation(ephemSrc, moon), unit=u.rad).degree  # ephem computes angles in radians !!
        logging.debug('[is_visible] Moon phase is {phase:.1f}% and the Moon/{src} separation is {sep:.1f}°'.format(phase=moonPhase, src=self.src, sep=targetMoonSeparation))

        # Add criteria on Moon phase and Moon/target separation
        if moonPhase < 40. and targetMoonSeparation >= 45. and targetMoonSeparation <= 145.:
            self.observableUnderMoonlight = True
        logging.debug('[is_visible] Is {src} observable under reasonable Moonlight conditions ? {ans}'.format(src=self.src, ans='Yes' if self.observableUnderMoonlight else 'No'))

        srcTransitTime = site.next_transit(ephemSrc)
        site.date = srcTransitTime
        ephemSrc.compute(site)
        srcAltAtTransit = c.Angle(ephemSrc.alt, unit=u.rad).degree

        # If srcAltAtTransit is below self.srcMinAlt, the source is just not optimally visible and we stop here
        logging.debug('[is_visible] Acceptable minimium altitude for {src}: {alt:.1f}°'.format(src=self.src, alt=self.srcMinAlt))
        if srcAltAtTransit < self.srcMinAlt:
            return False

        # Compute start and end of darkness time
        # (from here, as if we're computing these during the day)
        if nextMoonset > nextSunset and nextMoonset < nextSunrise and not self.observableUnderMoonlight:
            beginDarkness = nextMoonset
        else:
            beginDarkness = nextSunset

        if nextMoonrise < nextSunrise and nextMoonrise > nextSunset and not self.observableUnderMoonlight:
            endDarkness = nextMoonrise
        else:
            endDarkness = nextSunrise

        site.date = beginDarkness
        ephemSrc.compute(site)
        srcAltAtStartDarkTime = c.Angle(ephemSrc.alt, unit=u.rad).degree

        site.date = endDarkness
        ephemSrc.compute(site)
        srcAltAtEndDarkTime = c.Angle(ephemSrc.alt, unit=u.rad).degree

        darknessDuration = (endDarkness - beginDarkness) * 24. * 60.  # day to minutes
        logging.debug('[is_visible] Darkness begins at {0}'.format(beginDarkness))
        logging.debug('[is_visible] Darkness ends at {0}'.format(endDarkness))
        logging.debug('[is_visible] Darkness duration: {0:.1f} minutes'.format(darknessDuration))
        logging.debug('[is_visible] {src} altitude at beginning of darkness: {alt:.1f}°'.format(src=self.src, alt=srcAltAtStartDarkTime))
        logging.debug('[is_visible] {src} transitting on {transit} at an altitude of {alt:.1f}°'.format(src=self.src, transit=srcTransitTime, alt=srcAltAtTransit))
        logging.debug('[is_visible] {src} altitude at end of darkness: {alt:.1f}°'.format(src=self.src, alt=srcAltAtEndDarkTime))

        # check if source is visible, above minAlt, during this night
        for step in range(0, np.int(darknessDuration)):
            site.date = beginDarkness.datetime() + datetime.timedelta(minutes=step)
            ephemSrc.compute(site)
            srcAlt = c.Angle(ephemSrc.alt, unit=u.rad).degree
            moon.compute(site)
            moonAlt = c.Angle(moon.alt, unit=u.rad).degree
            moonUp = moonAlt > MOON_TWILIGHT_ANGLE.value
            # logging.debug('[is_visible] Time loop: It is now {now} and {src} is at alt. of {alt:.2f}°. Moon is {isup} at alt. of {moonAlt:.1f}°'.format(now=site.date,
            #                                                                                                                       src=self.src,
            #                                                                                                                       alt=srcAlt,
            #                                                                                                                       isup='up' if moonUp else 'down',
            #                                                                                                                       moonAlt=moonAlt))
            if srcAlt > self.srcMinAlt and not visibleFlag:
                visibleFlag = True
                logging.info('[{src:s}] {src:s} starts to be optimally visible, that is above an altitude of {alt:.2f}°, at {date}'.format(src=self.src, alt=self.srcMinAlt, date=site.date))
            if visibleFlag and not moonUp:
                self.visibleInDarkness = True
            if moonUp and visibleFlag and self.observableUnderMoonlight and not self.visibleUnderMoonlight:
                self.visibleUnderMoonlight = True

        if visibleFlag and self.visibleUnderMoonlight and not self.visibleInDarkness:
            self.visibleUnderMoonlightOnly = True

        logging.debug('[is_visible] Is {src} visible at all during the night ? {vis}'.format(src=self.src, vis='Yes' if visibleFlag else 'No'))
        logging.debug('[is_visible] Is {src} visible under reasonable Moonlight ? {vis}'.format(src=self.src, vis='Yes' if self.visibleUnderMoonlight else 'No'))
        logging.debug('[is_visible] Is {src} visible only under reasonable Moonlight and never during dark time ? {vis}'.format(src=self.src, vis='Yes' if self.visibleUnderMoonlightOnly else 'No'))

        if self.daily:
            self.createVisibilityPlot(beginDarkness)

        return visibleFlag


    def killTrigger(self):
        """
        Defines cuts on (RA,Dec,z) before assessing whether a mail alert should be sent for a source which flux is above the trigger threshold.
        We cut on a combination (z, ZenithAngle), using a bit mask.

        @rtype bool
        @todo      Introduce an additional cut on Gal latitude ?

        The 'return' value is a bit counter-intuitive. It answers the question 'Should we kill an imminent mail alert ?', i.e. if a source has the last flux point above the flux threshold, does it also fulfill the requirements on both z (not too far away) and zenith angle (not too low in the sky) ? So if an alert should definitely be sent, this function returns 'False' !
        """

        # Numpy array
        # combination of acceptable
        #                           z         ZA@culmination
        grid = np.array(zip(self.maxz, self.maxZA))

        zaAtCulmin = self.zaAtCulmination()

        # If input z is None, make it believe it is 0, otherwise msk crashes:
        if self.z == 'None':
            z = 0.
        else:
            z = float(self.z)

        # Mask on both (z, ZA at culmin)
        #          z column               ZA column
        msk = (z <= grid[:, 0]) & (zaAtCulmin <= grid[:, 1])

        # Assess whether the source is currently visible at the provided site
        if self.checkVisibility == 'True':
            self.visible = self.is_visible()
        else:
            # The source is assumed to be visible in any case, i.e. we don't care about its visibility status at the provided site to send a potential alert
            self.visible = True

        # if the mask has at least one 'True' element, we should send an alert
        if True in msk and self.visible:
            # An alert should be triggered !
            return False
        else:
            # No alert triggered
            return True

    def dynamicalTrigger(self):
        '''
        If long-term data are available for a source, dynamically computes a flux trigger threshold based on the flux history of the source. Otherwise, fall back with default fixed trigger threshold.
        '''

        # Read the longterm .dat LC file
        infile = self.baseOutDir + '/longTerm/merged/' + str(self.src) + '_lc.dat'
        try:
            data = ascii.read(infile)
        except IOError:
            logging.error('[{src:s}] \033[95m* Long term data file unavailable for source {src:s}\033[0m'.format(src=self.src))
            # Falling back to default fixed trigger threshold
            self.withhistory = False
            self.LTfluxAverage = self.threshold
            self.LTfluxRMS = 0.
            return None

        flux = data['Flux']
        fluxErr = data['FluxError']
        try:
            from uncertainties import unumpy as unp
            logging.info('[{src:s}] The long-term flux average is {f} ph cm^-2 s^-1'.format(src=self.src, f=unp.uarray(flux, fluxErr).mean()))
        except:
            logging.debug('Cannot use uncertainties...')
            pass

        # weighted average of the historical fluxes, weighted by their errors
        try:
            self.LTfluxAverage = np.average(flux, weights=1. / fluxErr)
            self.LTfluxRMS = np.std(flux, dtype=np.float64)
        except ZeroDivisionError:
            logging.error('[{src:s}] \033[95m* Zero division error for computation of long term flux average for source {src:s}\033[0m'.format(src=self.src))
            self.withhistory = False
            self.LTfluxAverage = self.threshold
            self.LTfluxRMS = 0.
            return None

        # Dynamically redefine the flux trigger threshold, using a 2-level criteria depending on whether we are currently looking at short- or long-term data
        if self.daily:
            self.threshold = self.LTfluxAverage + self.sigma * self.LTfluxRMS
        else:
            self.threshold = self.LTfluxAverage + self.sigmaLT * self.LTfluxRMS

    def Triggered(self):
        '''
        Has the source fulfilled the trigger conditions ?

        @return True
        @rtype bool
        '''

        # Read the light curve file
        if self.daily:
            infile = self.workDir + '/' + str(self.src) + '_daily_lc.dat'
            infilefits = self.workDir + '/' + str(self.src) + '_daily_lc.fits'
            self.pngFig = self.workDir + '/' + str(self.src) + '_daily_lc.png'

            # Also take a look in the long time-binned data
            infileLongTimeBin = self.workDir + '/' + str(self.src) + '_lc.dat'
            dataLongTimeBin = ascii.read(infileLongTimeBin)
            timeLongTimeBin = dataLongTimeBin['MET']
            fluxLongTimeBin = dataLongTimeBin['Flux']
            fluxErrLongTimeBin = dataLongTimeBin['FluxError']
            # Catch the last flux point
            self.lastTimeLongTimeBin = timeLongTimeBin[-1:]
            self.lastFluxLongTimeBin = fluxLongTimeBin[-1:]
            self.lastFluxErrLongTimeBin = fluxErrLongTimeBin[-1:]

            # Get the arrival time of the last photon analysed
            photonfileLongTimeBin = self.workDir + '/' + str(self.src) + '_gti.fits'
            photonsLongTimeBin = fits.open(photonfileLongTimeBin)
            photonsLongTimeBinTime = photonsLongTimeBin[1].data.field('TIME')
            self.arrivalTimeLastPhotonLongTimeBin = photonsLongTimeBinTime[-1:]

            photonfile = self.workDir + '/' + str(self.src) + '_daily_gti.fits'
            photons = fits.open(photonfile)
            photonsTime = photons[1].data.field('TIME')
            self.arrivalTimeLastPhoton = photonsTime[-1:]
        else:
            infile = self.workDir + '/' + str(self.src) + '_lc.dat'
            infilefits = self.workDir + '/' + str(self.src) + '_lc.fits'
            self.pngFig = self.workDir + '/' + str(self.src) + '_lc.png'

            photonfile = self.workDir + '/' + str(self.src) + '_gti.fits'
            photons = fits.open(photonfile)
            photonsTime = photons[1].data.field('TIME')
            self.arrivalTimeLastPhoton = photonsTime[-1:]
        data = ascii.read(infile)
        time = data['MET']
        flux = data['Flux']
        fluxErr = data['FluxError']
        self.lastExposure = float(fits.open(infilefits)['RATE'].data.field('EXPOSURE')[-1:])

        # Catch the last flux point
        self.lastTime = time[-1:]
        self.lastFlux = flux[-1:]
        self.lastFluxErr = fluxErr[-1:]

        self.energyTimeFig = self.workDir + '/' + str(self.src) + '_energyTime.png'

        logging.debug('[{src:s}] threshold={th:.2g}, lastFlux={lastfl:.2g} +/- {lastflerr:.2g} ph cm^-2 s^-1'.format(src=self.src,
                                                                                                                     th=float(self.threshold),
                                                                                                                     lastfl=float(self.lastFlux),
                                                                                                                     lastflerr=float(self.lastFluxErr)))

        # Do we kill potential trigger due to (ra, dec, z) cut ?
        self.triggerkilled = self.killTrigger()

        # Assess whether flux is above threshold, looking at the last flux point and at the exposure of the last flux measurement
        if self.lastFlux >= self.LTfluxAverage + self.sigma * np.sqrt(self.LTfluxRMS**2 + self.lastFluxErr**2):
            self.active = True
            if self.lastExposure < self.lowerExposure:
                self.active = False
                logging.warning('[{src}] \033[93mPotential trigger killed due to low exposure: last exposure is {lastExp:.2g} cm^2 s versus cutting on the exposure at {lowExp:.2g} cm^2 s.\033[0m'.format(src=self.src, lastExp=self.lastExposure, lowExp=self.lowerExposure))
        else:
            self.active = False

        # Combine killTrigger and flux above threshold criteria, force an alert in case of ToO
        SENDALERT = False
        if ((not self.triggerkilled and self.active) or self.forcealert or self.too) is True:
            SENDALERT = True

        logging.debug('[{src}] dec={dec}, z={z}, maxZA={za}, maxz={maxz}, active={active}, visible={visible}, triggerkilled={kill}, sendalert={send}, too={too}, forcealert={forcealert}'.format(src=self.src,
                                                                                                                                                                                                 dec=self.dec,
                                                                                                                                                                                                 z=self.z,
                                                                                                                                                                                                 za=self.maxZA,
                                                                                                                                                                                                 maxz=self.maxz,
                                                                                                                                                                                                 active=self.active,
                                                                                                                                                                                                 visible=self.visible,
                                                                                                                                                                                                 kill=self.triggerkilled,
                                                                                                                                                                                                 send=SENDALERT,
                                                                                                                                                                                                 too=self.too,
                                                                                                                                                                                                 forcealert=self.forcealert))
        return SENDALERT

    def sendAlert(self, nomailall=False, sendmail=False):
        '''
        Send a mail alert in case a source fulfills the trigger conditions.

        @param nomailall Boolean, should the mail be sent to a restricted list of recipients ?
        @return True
        @rtype bool
        '''

        # Import modules
        try:
            # Import smtplib to send mails
            import smtplib

            # Here are the email package modules we'll need
            from email.MIMEMultipart import MIMEMultipart
            from email.MIMEText import MIMEText
            from email.MIMEBase import MIMEBase
            from email import Encoders

        except:
            logging.error('sendAlert: Can\'t import mail modules.')
            sys.exit(1)

        SENDALERT = self.Triggered()

        # If trigger condition is met, we send a mail
        if SENDALERT and sendmail:
            # Create the container email message.
            msg = MIMEMultipart()
            sender = self.mailSender

            fhlName = self.search2FHLcounterpart()
            if fhlName is not None:
                fhlmessage = "2FHL counterpart is %s" % fhlName
            else:
                fhlmessage = "No 2FHL counterpart found"

            fglName = self.search4FGLcounterpart()
            if fglName is not None:
                fglmessage = "4FGL counterpart is %s" % fglName
            else:
                fglmessage = "No 4FGL counterpart found"

            # To whom the mail should be sent (cf. __init__ function of the class)
            if self.too == 'True':
                tootitleflag = " *** ToO *** Fermi-LAT report"
            else:
                tootitleflag = "Fermi-LAT flare alert"
            if not nomailall:
                recipient = self.usualRecipients
                msg['Subject'] = '[FLaapLUC] %s on %s [2FHL counterpart: %s]' % (tootitleflag, self.src, fhlName)
            else:
                recipient = self.testRecipients
                msg['Subject'] = '[FLaapLUC TEST MAIL] %s on %s [2FHL counterpart: %s]' % (
                tootitleflag, self.src, fhlName)

            msg['From'] = sender
            COMMASPACE = ', '
            msg['To'] = COMMASPACE.join(recipient)
            msg.preamble = 'You will not see this in a MIME-aware mail reader.\n'
            # Guarantees the message ends in a newline
            msg.epilogue = ''

            if not self.too:
                stateinfo = "exceeds"
            else:
                if self.lastFlux > self.threshold:
                    stateinfo = "*exceeds*"
                else:
                    stateinfo="*does not* exceed"

            mailtext = """
     FLaapLUC (Fermi/LAT automatic aperture photometry Light C<->Urve) report

     *** The Fermi/LAT flux (%.0f MeV-%.0f GeV) of %s (%s, %s) %s the trigger threshold of %.2g ph cm^-2 s^-1 ***

     Source coordinates (RA, Dec): %f %f

     """ % (self.emin, self.emax / 1000., self.src, fhlmessage, fglmessage, stateinfo, self.threshold, self.ra, self.dec)


            if self.daily:
                mailtext = mailtext + """

     The last daily-binned flux is:        %.2g +/- %.2g ph cm^-2 s^-1, centred on MET %.0f (MJD %.5f, i.e. %s) (arrival time of last photon analysed: MET %.0f, MJD %.5f, %s)
     and the last %.0f-day binned flux is: %.2g +/- %.2g ph cm^-2 s^-1, centred on MET %.0f (MJD %.5f, i.e. %s) (arrival time of last photon analysed: MET %.0f, MJD %.5f, %s)

""" % (self.lastFlux,
       self.lastFluxErr,
       self.lastTime, extras.met2mjd(self.lastTime), str(extras.mjd2gd(extras.met2mjd(self.lastTime))),
       self.arrivalTimeLastPhoton, extras.met2mjd(self.arrivalTimeLastPhoton),
       str(extras.mjd2gd(extras.met2mjd(self.arrivalTimeLastPhoton))),
       self.longtimebin,
       self.lastFluxLongTimeBin,
       self.lastFluxErrLongTimeBin,
       self.lastTimeLongTimeBin, extras.met2mjd(self.lastTimeLongTimeBin),
       str(extras.mjd2gd(extras.met2mjd(self.lastTimeLongTimeBin))),
       self.arrivalTimeLastPhotonLongTimeBin, extras.met2mjd(self.arrivalTimeLastPhotonLongTimeBin),
       str(extras.mjd2gd(extras.met2mjd(self.arrivalTimeLastPhotonLongTimeBin))))
                mailtext = mailtext + "The most recent lightcurve (%.0f-day binned in red, and %.0f-day binned in blue) is attached." % (
                self.tbin / 24. / 60. / 60., self.longtimebin)
            else:
                mailtext = mailtext + """

     The last %.0f-day binned flux is:      %.2g +/- %.2g ph cm^-2 s^-1, centred on MET %.0f (MJD %.5f, i.e. %s) (arrival time of last photon analysed: %.0f, MJD %.5f, %s)

""" % (self.longtimebin,
       self.lastFlux,
       self.lastFluxErr,
       self.lastTime, extras.met2mjd(self.lastTime), str(extras.mjd2gd(extras.met2mjd(self.lastTime))),
       self.arrivalTimeLastPhoton, extras.met2mjd(self.arrivalTimeLastPhoton),
       str(extras.mjd2gd(extras.met2mjd(self.arrivalTimeLastPhoton))))
                mailtext = mailtext + "The most recent lightcurve (%.0f-day binned) is attached." % (
                self.tbin / 24. / 60. / 60.)

            if self.visibleUnderMoonlightOnly:
                mailtext = mailtext + """

     *NOTE*: %s is visible from your site, but only observable under *moonlight* conditions.
""" % (self.src)
            elif self.visibleUnderMoonlight and not self.visibleUnderMoonlightOnly:
                mailtext = mailtext + """

     *NOTE*: %s is visible from your site, and observable both under *dark* and *moonlight* conditions.
""" % (self.src)
            elif self.visibleInDarkness and not self.visibleUnderMoonlight:
                mailtext = mailtext + """

     *NOTE*: %s is visible from your site, and observable under *dark* conditions.
""" % (self.src)

            if self.too == 'True':
                mailtext=mailtext+"""

     *NOTE*: a ToO is currently ongoing on %s, this alert is forced.
"""%(self.src)

            if FLAGASSUMEDGAMMA is True:
                mailtext = mailtext + """

     *WARNING*: The source %s is not found in the 4FGL catalogue, its photon index is thus assumed to be %.2f for the light curve computation.
""" % (self.src, ASSUMEDGAMMA)

            mailtext = mailtext + """

     *Disclaimer*: Be careful, though, that these light curves are not computed using the usual, clean, standard (un)binned likelihood procedure one should normally use for a good quality, publication-ready result. Those reported here only rely on a "quick & dirty" aperture photometric analysis (cf. e.g. http://fermi.gsfc.nasa.gov/ssc/data/analysis/scitools/aperture_photometry.html), which basically assumes that the data set, within 1 degree around the source, is background-free.

      Cheers,
      FLaapLUC.
"""

            txt = MIMEText(mailtext)
            msg.attach(txt)

            # Attach the figures
            for fig in [self.pngFig, self.energyTimeFig, self.visibilityFig]:
                try:
                    # Open the files in binary mode.  Let the MIMEImage class automatically guess the specific image type.
                    fp = open(fig, 'rb')
                    img = MIMEBase('application', 'octet-stream')
                    img.set_payload(fp.read())
                    Encoders.encode_base64(img)
                    img.add_header('Content-Disposition',
                                   'attachment; filename="%s"' % os.path.basename(fig))
                    fp.close()
                    msg.attach(img)
                except:
                    pass

            # Send the email via our own SMTP server.
            s = smtplib.SMTP()
            s.set_debuglevel(0)
            s.connect()
            s.sendmail(sender, recipient, msg.as_string())
            s.quit()

            logging.info('\033[94m*** Alert sent for {src:s}\033[0m'.format(src=self.src))

            return True
        else:
            return False


    def search3FGLcounterpart(self):
        """
        Search the 3FGL name of a 2FGL or a 4FGL source name
        """
        if self.fglName is None:
            return None
        if "3FGL" in self.fglName:
            return self.fglName.replace('_3FGLJ', '3FGL J').replace('3FGLJ', '3FGL J')

        # Search 3FGL name from a 2FGL name
        if '2FGL' in self.fglName:
            cat3FGLfile = self.catalogFile.replace('gll_psc_v08', 'gll_psc_v16')
            hdulist = fits.open(cat3FGLfile)
            cat = hdulist[1].data
            logging.debug('[{src:s}] 2FGL name is {fgl}'.format(src=self.src, fgl=self.fglName.replace('_2FGLJ', '2FGL J').replace('2FGLJ', '2FGL J')))

            found = False
            for stuff in cat:
                if stuff.field('2FGL_Name') == self.fglName.replace('_2FGLJ', '2FGL J').replace('2FGLJ', '2FGL J'):
                    threefglName = stuff.field('Source_Name')
                    logging.info('[{src:s}] Found the 3FGL counterpart of {two:s}: {three:s}'.format(src=self.src, two=self.fglName, three=threefglName))
                    found = True
                    break

            if not found:
                threefglName = None
                logging.info('[{src:s}] No 3FGL counterpart found for {two:s}'.format(src=self.src, two=self.fglName))

            hdulist.close()
            return threefglName

        # Search 3FGL name from a 4FGL name
        if '4FGL' in self.fglName:
            cat4FGLfile = self.catalogFile.replace('/3FGL/','/4FGL-DR2/').replace('gll_psc_v16', 'gll_psc_v27')
            hdulist = fits.open(cat4FGLfile)
            cat = hdulist[1].data
            logging.debug('[{src:s}] 4FGL name is {fgl}'.format(src=self.src, fgl=self.fglName.replace('_4FGLJ', '4FGL J').replace('4FGLJ', '4FGL J')))

            found = False
            for stuff in cat:
                if stuff.field('Source_Name') == self.fglName.replace('_4FGLJ', '4FGL J').replace('4FGLJ', '4FGL J') and '3FGL' in stuff.field('ASSOC_FGL'):
                    threefglName = stuff.field('ASSOC_FGL')
                    logging.info('[{src:s}] Found the 3FGL counterpart of {two:s}: {three:s}'.format(src=self.src, two=self.fglName, three=threefglName))
                    found = True
                    break

            if not found:
                threefglName = None
                logging.info('[{src:s}] No 3FGL counterpart found for {two:s}'.format(src=self.src, two=self.fglName))

            hdulist.close()
            return threefglName


    def search4FGLcounterpart(self, threefglName=None):
        """
        Search the 4FGL name of a 3FGL source name
        """
        if self.fglName is None:
            return None
        if "4FGL" in self.fglName:
            return self.fglName.replace('4FGLJ', '4FGL J')

        cat4FGLfile = self.catalogFile.replace('gll_psc_v16', 'gll_psc_v27')
        hdulist = fits.open(cat4FGLfile)
        cat = hdulist[1].data
        logging.debug('[{src:s}] 3FGL name is {fgl}'.format(src=self.src, fgl=self.fglName.replace('3FGLJ', '3FGL J')))

        found = False
        for stuff in cat:
            if threefglName is not None and stuff.field('ASSOC_FGL') == threefglName.replace('3FGLJ', '3FGL J'):
                fourfglName = stuff.field('Source_Name')
                logging.info('[{src:s}] Found the 4FGL counterpart of {three:s}: {four:s}'.format(src=self.src, three=threefglName, four=fourfglName))
                found = True
                break
            elif threefglName is None and stuff.field('ASSOC_FGL') == self.fglName.replace('3FGLJ', '3FGL J'):
                fourfglName = stuff.field('Source_Name')
                logging.info('[{src:s}] Found the 4FGL counterpart of {three:s}: {four:s}'.format(src=self.src, three=self.fglName, four=fourfglName))
                found = True
                break

        if not found:
            fourfglName = None
            logging.info('[{src:s}] No 4FGL counterpart found for {three:s}'.format(src=self.src, three=self.fglName))

        hdulist.close()
        return fourfglName


    def search2FHLcounterpart(self):
        """
        Search the 2FHL name of a 2FGL, a 3FGL or a 4FGL source name
        """
        if self.fglName is not None:
            if "2FHL" in self.fglName:
                return self.fglName.replace('_2FHLJ', '2FHL J').replace('2FHLJ', '2FHL J')

            cat2FHLfile = self.catalogFile.replace('/3FGL/',
                                                   '/2FHL/').replace('/4FGL-DR2',
                                                                     '/2FHL/').replace('psc_v08',
                                                                                       'psch_v08').replace('psc_v16',
                                                                                                           'psch_v08').replace('psc_v27',
                                                                                                                               'psch_v08')
            try:
                hdulist = fits.open(cat2FHLfile)
            except IOError:
                logging.info('2FHL catalog file not found')
                return None
            cat = hdulist[1].data

            found = False
            threefglName = self.search3FGLcounterpart()
            for stuff in cat:
                if stuff.field('3FGL_Name') == self.fglName.replace('_3FGLJ', '3FGL J').replace('3FGLJ',
                                                                                                '3FGL J') or stuff.field(
                        '3FGL_Name') == str(threefglName).replace('3FGLJ', '3FGL J'):
                    fhlName = stuff.field('Source_Name')
                    logging.info('[{src:s}] Found the 2FHL counterpart of {fgl:s}: {fhl:s}'.format(src=self.src, fgl=self.fglName, fhl=fhlName))
                    found = True
                    break

            if not found:
                fhlName = None
                logging.info('[{src:s}] No 2FHL counterpart found for {fgl:s}'.format(src=self.src, fgl=self.fglName))

            hdulist.close()
            return fhlName
        else:
            return None
