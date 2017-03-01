#!/bin/env python
# -*- coding: utf-8 -*-
#
# Time-stamp: "2017-03-01 22:29:48 jlenain"

"""
FLaapLUC (Fermi/LAT automatic aperture photometry Light C<->Urve)

Automatic generation of aperture photometric light curves of
high energy sources, for a given source.

No likelihood fit is performed, the results solely rely on the
3FGL spectral fits, if available.

More information are available at:
http://fermi.gsfc.nasa.gov/ssc/data/analysis/scitools/aperture_photometry.html

@author Jean-Philippe Lenain <mailto:jlenain@in2p3.fr>
"""

import sys
import os
import asciidata
import datetime
import time
import glob
from numpy import *
import pyfits
import ephem
from astLib import astCoords
from optparse import OptionParser
from ConfigParser import ConfigParser

# Import some matplotlib modules
try:
    import matplotlib
    matplotlib.use('Agg')

    from matplotlib.pyplot import *
    from matplotlib.ticker import FuncFormatter
except ImportError:
    print "ERROR Can't import matplotlib"
    sys.exit(1)

# Import the Science Tools modules
try:
    from gt_apps import *
except ImportError:
    print "ERROR Can't import the Fermi Science tools"
    sys.exit(1)

# Flags
DEBUG = False
VERBOSE = False
BATCH = True
FORCE_ALERT=False
FORCE_LIKELIHOOD=False
# Flag to know whether Gamma is assumed to be ASSUMEDGAMMA
# or taken from the 3FGL.
FLAGASSUMEDGAMMA = False

# Global variables
TOFFSET = 54000.  # offset in MJD for plot creation
# assumed photon index for a source not belonging to the 3FGL
ASSUMEDGAMMA = -2.5


def met2mjd(met):
    """
    Converts Mission Elapsed Time (MET, in seconds) in Modified Julian Day.
    Cf. http://fermi.gsfc.nasa.gov/ssc/data/analysis/documentation/Cicerone/Cicerone_Data/Time_in_ScienceTools.html
    to see how the time is handled in the Fermi Science Tools.

    Input: time in MET (s)
    Output: time in MJD (fraction of a day)
    """
    MJDREFI = 51910.0
    MJDREFF = 7.428703703703703e-4
    return(MJDREFI+MJDREFF+met/24./60./60.)


def mjd2met(mjd):
    """
    Converts Modified Julian Day in Mission Elapsed Time (MET, in seconds).
    Cf. http://fermi.gsfc.nasa.gov/ssc/data/analysis/documentation/Cicerone/Cicerone_Data/Time_in_ScienceTools.html
    to see how the time is handled in the Fermi Science Tools.

    Input:  time in MJD (fraction of a day)
    Output: time in MET (s)
    """
    MJDREFI = 51910.0
    MJDREFF = 7.428703703703703e-4
    return(24. * 60. * 60 * (mjd - MJDREFI - MJDREFF))


def unixtime2mjd(unixtime):
    """
    Converts a UNIX time stamp in Modified Julian Day

    Highly inspired from the function 'mjd_now()' from
    Marcus Hauser's ADRAS/ATOM pipeline

    Input:  time in UNIX seconds
    Output: time in MJD (fraction of a day)

    @todo Should make sure we use UTC here
    @todo We do not yet account for leap seconds (2nd order compared to UTC handling !!)
    """

    # unixtime gives seconds passed since "The Epoch": 1.1.1970 00:00
    # MJD at that time was 40587.0

    result = 40587.0 + unixtime / (24.*60.*60.)
    return result


def jd2gd(x):
    """
    Compute gregorian date out of julian date

    input: julian date x (float)
    return value: string of gregorian date

    based on/copied from script jd2dg.py from Enno Middelberg
    http://www.atnf.csiro.au/people/Enno.Middelberg/python/jd2gd.py

    task to convert a list of julian dates to gregorian dates
    description at http://mathforum.org/library/drmath/view/51907.html
    Original algorithm in Jean Meeus, "Astronomical Formulae for Calculators"
    """

    jd = float(x)

    jd = jd+0.5
    Z = int(jd)
    F = jd-Z
    alpha = int((Z-1867216.25)/36524.25)
    A = Z + 1 + alpha - int(alpha/4)

    B = A + 1524
    C = int((B-122.1)/365.25)
    D = int(365.25*C)
    E = int((B-D)/30.6001)

    dd = B - D - int(30.6001*E) + F

    if E < 13.5:
        mm = E-1

    if E > 13.5:
        mm = E-13

    if mm > 2.5:
        yyyy = C-4716

    if mm < 2.5:
        yyyy = C-4715

    # months=["January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"]

    daylist = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    daylist2 = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

    h = int((dd-int(dd))*24)
    min = int((((dd-int(dd))*24)-h)*60)
    sec = 86400*(dd-int(dd))-h*3600-min*60

    # Now calculate the fractional year. Do we have a leap year?
    if (yyyy % 4 != 0):
        days = daylist2
    elif (yyyy % 400 == 0):
        days = daylist2
    elif (yyyy % 100 == 0):
        days = daylist
    else:
        days = daylist2

    # print x+" = "+months[mm-1]+" %i, %i, " % (dd, yyyy), 
    # print string.zfill(h,2)+":"+string.zfill(min,2)+":"+string.zfill(sec,2)+" UTC"

    string = "%04d-%02d-%02d %02d:%02d:%04.1f" % (yyyy, mm, dd, h, min, sec)

    return string


###### modified julian date to gregorian date ########################
def mjd2gd(time):
    """
    Converts Modified Julian Day in Gregorian Date.

    Under the hood, it calls jd2gd().
    """

    return jd2gd(time+2400000.5)


def rad2deg(angle):
    """
    Convert an angle from radians to degrees.

    @param angle in radians
    """
    return angle*180./pi


def deg2rad(angle):
    """
    Convert an angle from degrees to radians.

    @param angle in degrees
    """
    return angle*pi/180.


def angsep((ra1,dec1), (ra2,dec2), deg=True):
    """
    Calculates the angular separation between two points on the sky.

    @param (ra1,dec1) coordinates of 1st source
    @param (ra2,dec2) coordinates of 2nd source
    @param deg flag whether inputs/outputs are in degrees or radians
    """
    if deg:
        ra1 = deg2rad(ra1)
        dec1 = deg2rad(dec1)
        ra2 = deg2rad(ra2)
        dec2 = deg2rad(dec2)

    SEP = arccos(cos(dec1)*cos(dec2)*cos(ra1-ra2)+sin(dec1)*sin(dec2))
    if deg:
        SEP = rad2deg(SEP)
    return SEP


def getConfigList(option, sep=','):
    return [stuff for stuff in option.split(sep)]


class autoLC:
    """
    FLaapLUC

    Automatic aperture photometry light curve generation.
    Main class, for a given of source.
    """

    def __init__(self, file=None, customThreshold=False, daily=False,
                 longTerm=False, yearmonth=None, mergelongterm=False,
                 withhistory=False, stopmonth=None, configfile='default.cfg'):
        
        self.config = self.getConfig(configfile=configfile)
        self.allskyDir = self.config.get('InputDirs', 'AllskyDir')
        self.archiveDir = self.config.get('InputDirs', 'ArchiveDir')
        self.templatesDir = self.config.get('InputDirs', 'TemplatesDir')
        self.ATOMSchedulesDir = self.config.get('InputDirs', 'ATOMSchedulesDir')
        self.catalogFile = self.config.get('InputFiles', 'CatalogFile')
        if file is None:
            self.file = self.config.get('InputFiles', 'SourceList')
        else:
            self.file = file
        self.baseOutDir = self.config.get('OutputDirs', 'OutputResultsDir')
        self.allskyFile = self.allskyDir + "/" + self.config.get('InputFiles', 'WholeAllskyFile')
        self.lastAllskyFile = self.allskyDir + "/"+self.config.get('InputFiles', 'LastAllskyFile')
        self.spacecraftFile = self.allskyDir + "/" + self.config.get('InputFiles', 'SpacecraftFile')
        self.webpageDir = self.config.get('OutputDirs', 'OutputWebpageDir')
        self.url = self.config.get('OutputDirs', 'URL')

        try:
            self.longtimebin = float(self.config.get('AlertTrigger', 'LongTimeBin'))
        except:
            # Take 7 days by default
            self.longtimebin = 7.
            print '\033[93mCan not read LongTimeBin in config file, taking %.1f as default.\033[0m' % (self.longtimebin)

        try:
            self.sigma = float(self.config.get('AlertTrigger', 'Sigma'))
        except:
            # Take 2 sigma by default
            self.sigma = 3.
            print '\033[93mCan not read Sigma in config file, taking %.1f as default.\033[0m' % (self.sigma)

        try:
            self.sigmaLT = float(self.config.get('AlertTrigger', 'SigmaLT'))
        except:
            # Take 2 sigma by default
            self.sigmaLT = 1.5
            print '\033[93mCan not read SigmaLT in config file, taking %.1f as default.\033[0m' % (self.sigmaLT)
     
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
        try:
            self.launchLikeAna = self.config.get('AlertTrigger', 'LaunchLikelihoodAnalysis')
        except:
            self.launchLikeAna = False

        try:
            self.likeHighEThresh = self.config.get('AlertTrigger', 'LikelihoodHighEnergyThreshold')
        except:
            self.likeHighEThresh = None

        self.daily = daily
        self.withhistory = withhistory

        # Mail sender and recipients
        self.usualRecipients = getConfigList(self.config.get('MailConfig', 'UsualRecipients'))
        self.testRecipients = getConfigList(self.config.get('MailConfig', 'TestRecipients'))
        self.mailSender = self.config.get('MailConfig', 'MailSender')
        try:
            self.likelihoodRecipients = self.config.get('MailConfig', 'LikelihoodRecipients')
        except:
            self.likelihoodRecipients = None

        today = datetime.date.today().strftime('%Y%m%d')

        # Setting file names and directories
        if longTerm:
            self.allsky = self.allskyFile
            if not mergelongterm:
                self.workDir = self.baseOutDir + "/longTerm/" + yearmonth
            else:
                self.workDir = self.baseOutDir + "/longTerm/merged"
        else:
            # If lock files exist in the archive directory (e.g. if NASA servers are down), we do not do anything and exit
            photonLock = self.archiveDir + "/photon.lock"
            spacecraftLock = self.archiveDir + "/spacecraft.lock"
            if os.path.isfile(photonLock) or os.path.isfile(spacecraftLock):
                self.sendErrorMail(mailall=True)
                sys.exit(10)

            self.allsky = self.lastAllskyFile
            self.workDir = self.baseOutDir + "/" + today

        self.spacecraft = self.spacecraftFile
        if not os.path.isdir(self.workDir):
            try:
	        os.makedirs(self.workDir)
	    except OSError:
		pass
        
        self.fermiDir   = os.getenv('FERMI_DIR')

        # Setting default parameters
        self.roi       = 1.   # degrees (http://fermi.gsfc.nasa.gov/ssc/data/analysis/scitools/aperture_photometry.html: "For aperture photometry we select a very small aperture (rad=1 degree), because we are not fitting the background.")
        try:
            self.emin  = float(self.config.get('Erange','Emin'))
        except:
            # Take 100 MeV by default
            self.emin      = 1.e2 # E min
            print '\033[93mCan not read Emin in config file, taking %.1g as default.\033[0m' % (self.emin)
        try:
            self.emax  = float(self.config.get('Erange','Emax'))
        except:
            # Take 500 GeV by default
            self.emax      = 5.e5 # E max
            print '\033[93mCan not read Emax in config file, taking %.1g as default.\033[0m' % (self.emax)
        self.zmax      = 100. # degrees
        self.rockangle = 52.  # maximal allowed rocking angle

        if self.daily:
            self.tbin =                  24.*60.*60. # seconds, daily bins
        else:
            self.tbin = self.longtimebin*24.*60.*60. # seconds, longtimebin by defaults

        self.threshold = 1.e-6 # ph cm^-2 s^-1
        self.customThreshold=customThreshold

        self.stopmonth = stopmonth
        
        # Open allsky file to get the start and stop dates
        try:
            hdu=pyfits.open(self.allsky)
        except IOError as e:
            print 'I/O error ({0}): can not open file {1}: {2}'.format(e.errno, self.allsky, e.strerror)
            raise
        header = hdu[0].header

        if not longTerm:
            self.tstart = header['TSTART']
            self.tstop  = header['TSTOP']
            
        else:
            missionStart = header['TSTART'] # in MET
            missionStop  = header['TSTOP']  # in MET

            if not mergelongterm:
                # Need to convert 'yearmonth' in MET
                # self.tstart is the first day of yearmonth at 00:00:00, or missionStart
                # self.tstop  is the first day of next month at 00:00:00, or missionStop
                year=yearmonth[:-2]
                month=yearmonth[-2:]

                # Get date of first day of yearmonth at 00:00:00, in UNIX time (timetuple transform a datetime object in time object ???)
                #                                                        year         month   day   hour   minute  second  microsecond
                yearmonthStart     = time.mktime(datetime.datetime(  int(year),   int(month),   1,     0,       0,      0,           0).timetuple())
                if int(month)<12:
                    yearmonthStop  = time.mktime(datetime.datetime(  int(year), int(month)+1,   1,     0,       0,      0,           0).timetuple())
                else:
                    yearmonthStop  = time.mktime(datetime.datetime(int(year)+1,            1,   1,     0,       0,      0,           0).timetuple())
            
                # Convert these from UNIX time to MET
                tmptstart = mjd2met(unixtime2mjd(yearmonthStart))
                tmptstop  = mjd2met(unixtime2mjd(yearmonthStop))

                if DEBUG:
                    print 'DEBUG: INIT yearmonthStart=',yearmonthStart
                    print 'DEBUG: INIT yearmonthStop=',yearmonthStop

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
                    
 

    def getConfig(self,configfile='./default.cfg'):
        """Get configuration from a configuration file."""
        self.config = ConfigParser()
        self.config.readfp(open(configfile))
        return self.config


    def readSourceList(self,mysrc=None):
        """
        Read the list of sources.

        @todo Use a mySQL database instead of an ASCII file for the list of sources ?
        """

        try:
            import asciidata
        except ImportError:
            print "ERROR Can't import asciidata, needed to read the list of sources. Aborting..."
            sys.exit(1)

        try:
            srcList=asciidata.open(self.file)
        except IOError:
            print "ERROR Can't open "+self.file
            sys.exit(1)

        src = srcList[0]
        ra  = srcList[1].tonumpy()
        dec = srcList[2].tonumpy()
        z   = srcList[3].tonumpy()
        fglName=srcList[4]
        # Read the threshold for the source from the source list, if we asked to process with custom thresholds when instanciating the class
        if self.customThreshold:
            myThreshold=srcList[5].tonumpy()

    
        # If we ask for a particular source, return the parameters for that source
        if mysrc != None:
            # Find our input src in the list of sources
            found=False
            for i in range(len(src)):
                if src[i]==mysrc:
                    found=True

                    # Redefine the threshold if we provided a custom threshold
                    if self.customThreshold and myThreshold[i] != 0.:
                        try:
                            float(myThreshold[i])
                            self.threshold=myThreshold[i]
                        except ValueError:
                            print 'WARNING The threshold of the source %s is not a float. Please, check the list of sources !' % mysrc
                            sys.exit(2)
                    self.src     = src[i]
                    self.ra      = ra[i]
                    self.dec     = dec[i]
                    self.z       = z[i]
                    self.fglName = fglName[i]
                    return
            
            # If we end up without any found source, print out a WARNING
            print 'WARNING Can\'t find your source %s in the list of sources !' % str(mysrc)
            self.src     = None
            self.ra      = None
            self.dec     = None
            self.z       = None
            self.fglName = None
            return
        
        # Otherwise, return the whole list of parameters for all the sources
        else:
            return src,ra,dec,z,fglName

    def selectSrc(self):
        """
        Filter a given source, running gtselect
        """
        # Do we have to deal with a FITS file or an ASCII list of FITS file ?
        allskyext = os.path.splitext(self.allsky)[1]
        if allskyext in [".fit", ".fits"]:
            filter['infile'] = self.allsky
        else:
            filter['infile'] = '@%s' % self.allsky
        if self.daily:
            outfile=self.workDir+'/'+str(self.src)+'_daily.fits'
        else:
            outfile=self.workDir+'/'+str(self.src)+'.fits'
        filter['outfile']=outfile
        
        # If outfile already exists, we don't do anything
        if os.path.isfile(outfile):
            return True

        filter['ra']      = self.ra
        filter['dec']     = self.dec
        filter['rad']     = self.roi
        filter['emin']    = self.emin
        filter['emax']    = self.emax
        filter['tmin']    = self.tstart
        filter['tmax']    = self.tstop
        filter['zmax']    = self.zmax
        filter['evclass'] = 128
        filter.run()


    def makeTime(self):
        """
        Filter the GTI for a given source
        """
        maketime['scfile']=self.spacecraft

        if self.daily:
            maketime['evfile']=self.workDir+'/'+str(self.src)+'_daily.fits'
            outfile=self.workDir+'/'+str(self.src)+'_daily_gti.fits'
        else:
            maketime['evfile']=self.workDir+'/'+str(self.src)+'.fits'
            outfile=self.workDir+'/'+str(self.src)+'_gti.fits'
        maketime['outfile']=outfile
        
        # If outfile already exists, we don't do anything
        if os.path.isfile(outfile):
            return True

        # cf. http://fermi.gsfc.nasa.gov/ssc/data/analysis/scitools/aperture_photometry.html
        #maketime['filter']="IN_SAA!=T && LAT_CONFIG==1 && DATA_QUAL==1 && (angsep(RA_ZENITH,DEC_ZENITH,"+str(ra)+","+str(dec)+")+"+str(self.roi)+" <"+str(self.zmax)+") && (angsep("+str(ra)+","+str(dec)+",RA_SCZ,DEC_SCZ)<180.) && (angsep("+str(ra)+","+str(dec)+",RA_SUN,DEC_SUN)>5.)"
        #maketime['filter'] = "IN_SAA!=T && LAT_CONFIG==1 && DATA_QUAL==1 && ABS(ROCK_ANGLE)<"+str(self.rockangle)+" && (angsep(RA_ZENITH,DEC_ZENITH,"+str(self.ra)+","+str(self.dec)+")+"+str(self.roi)+" <"+str(self.zmax)+") && (angsep("+str(self.ra)+","+str(self.dec)+",RA_SUN,DEC_SUN)>5.)"
        maketime['filter'] = "LAT_CONFIG==1 && DATA_QUAL>0"
        maketime['roicut'] = 'no'
        maketime['tstart'] = self.tstart
        maketime['tstop']  = self.tstop
        maketime.run()


    def mergeGTIfiles(self):
        """
        Merge multiple GTI files when mergelongterm is True.
        Use gtselect.
        Assume the current workDir is longTerm/merged.
        """

        # Create list of GTI files
        if not self.daily:
            listname=self.workDir+'/'+self.src+'_gti.list'
        else:
            listname=self.workDir+'/'+self.src+'_daily_gti.list'
        filelist=open(listname,'w')
        list=[]
        if not self.daily:
            for file in glob.glob(self.workDir+'/../20????/'+self.src+'_gti.fits'):
                list.append(file)
        else:
            for file in glob.glob(self.workDir+'/../20????/'+self.src+'_daily_gti.fits'):
                list.append(file)
        # Sort the list of GTI files
        list=sorted(list)
        for item in list:
            filelist.write(item+'\n')
        filelist.close()
        
        filter['infile']='@'+listname
        if not self.daily:
            outfile=self.workDir+'/'+str(self.src)+'_gti.fits'
        else:
            outfile=self.workDir+'/'+str(self.src)+'_daily_gti.fits'
        filter['outfile']=outfile
        
        # If outfile already exists, we re-create it
        if os.path.isfile(outfile):
            os.remove(outfile)

        filter['ra']      = self.ra
        filter['dec']     = self.dec
        filter['rad']     = self.roi
        filter['emin']    = self.emin
        filter['emax']    = self.emax
        filter['tmin']    = self.tstart
        filter['tmax']    = self.tstop
        filter['zmax']    = self.zmax
        filter['evclass'] = 128
        if VERBOSE:
            print 'INFO Running gtmktime'
        filter.run()


    def createXML(self):
        """
        Create an XML model file based on the 3FGL catalogue
        """
        
        if self.daily:
            evfile=self.workDir+'/'+str(self.src)+'_daily_gti.fits'
            modelfile=self.workDir+'/'+str(self.src)+'_daily.xml'
        else:
            evfile=self.workDir+'/'+str(self.src)+'_gti.fits'
            modelfile=self.workDir+'/'+str(self.src)+'.xml'

        # If modelfile already exists, we don't do anything
        if os.path.isfile(modelfile):
            return True

        try:
            import make3FGLxml
        except ImportError:
            print "ERROR Can't import make3FGLxml."
            sys.exit(1)
        
        mymodel=make3FGLxml.srcList(self.catalogFile,evfile,modelfile)
        if VERBOSE:
            print 'INFO Running makeModel'
        mymodel.makeModel(GDfile=self.fermiDir+'/refdata/fermi/galdiffuse/gll_iem_v06.fits',GDname='GalDiffuse',ISOfile=self.fermiDir+'/refdata/fermi/galdiffuse/iso_P8R2_SOURCE_V6_v06.txt',ISOname='IsotropicDiffuse',extDir=self.templatesDir,makeRegion=False)


    def photoLC(self):
        """
        Compute the photometric light curve for a given source
        """

        if self.daily:
            evtbin['evfile']=self.workDir+'/'+str(self.src)+'_daily_gti.fits'
            outfile=self.workDir+'/'+str(self.src)+'_daily_lc.fits'
        else:
            evtbin['evfile']=self.workDir+'/'+str(self.src)+'_gti.fits'
            outfile=self.workDir+'/'+str(self.src)+'_lc.fits'

        # If outfile already exists, we don't do anything
        if os.path.isfile(outfile):
            return True

        evtbin['outfile']   = outfile
        evtbin['scfile']    = self.spacecraft
        evtbin['algorithm'] = 'LC'
        evtbin['tbinalg']   = 'LIN'
        evtbin['tstart']    = self.tstart
        evtbin['tstop']     = self.tstop
        evtbin['dtime']     = self.tbin
        if VERBOSE:
            print 'INFO Running gtbin'
        evtbin.run()


    def exposure(self,gamma=None):
        """
        Compute exposure on source src, to add a flux column for the photometric light curve.

        Warning: the input file is modified in place, with an additional exposure column added to the file !
        """

        if self.daily:
            infile=self.workDir+'/'+str(self.src)+'_daily_lc.fits'
            srcmdl=self.workDir+'/'+str(self.src)+'_daily.xml'
        else:
            infile=self.workDir+'/'+str(self.src)+'_lc.fits'
            srcmdl=self.workDir+'/'+str(self.src)+'.xml'

        # If infile already contains an EXPOSURE column, we don't do anything
        hdu=pyfits.open(infile)
        if hdu[1].header.get('TTYPE5')=='EXPOSURE':
            return True
 
        scfile=self.spacecraft
        irfs='P8R2_SOURCE_V6'
        rad=str(self.roi)

        options='infile='+infile+' scfile='+scfile+' irfs='+irfs+' rad='+rad
        if self.fglName is not None:
            target=self.fglName.replace('3FGLJ','3FGL J')
            if DEBUG:
                print 'DEBUG: exposure: target=%s' % target
            options+=' srcmdl='+srcmdl+' target="'+target+'"'
        else:
            options+=' srcmdl="none" specin='+str(gamma)
        cmd='time -p '+self.fermiDir+'/bin/gtexposure '+options
        if VERBOSE:
            print 'INFO Running %s' % cmd
        os.system(cmd)


    def createDAT(self):
        """
        Create a data file with the light curve of a given source.
        """

        # Read LC file
        if self.daily:
            infile=self.workDir+'/'+str(self.src)+'_daily_lc.fits'
            outfile=self.workDir+'/'+str(self.src)+'_daily_lc.dat'
        else:
            infile=self.workDir+'/'+str(self.src)+'_lc.fits'
            outfile=self.workDir+'/'+str(self.src)+'_lc.dat'

        # If outfile already exists, we don't do anything
        if os.path.isfile(outfile):
            return True

        import pyfits
        try:
            hdu=pyfits.open(infile)
        except:
            print 'Exception: can not open file '+infile
            raise
        data=hdu[1].data


        file=open(outfile,'w')
        file.write("# Time[MET]\tTime[MJD]\tFlux[ph.cm^-2.s^-1]\tFluxError[ph.cm^-2.s^-1]\n")
        time      = data.field('TIME')     # MET
        counts    = data.field('COUNTS')
        countsErr = data.field('ERROR')    # error on counts
        exposure  = data.field('EXPOSURE') # cm^2 s^1
        flux      = counts/exposure        # approximate flux in ph cm^-2 s^-1
        fluxErr   = countsErr/exposure     # approximate flux error in ph cm^-2 s^-1

        timeMjd=met2mjd(time)
        # We can do this because time is NOT a list, but a numpy.array
        
        for i in range(len(time)):
            # Exposure can be 0 if longTerm=True and TSTOP in photon file > TSTOP in spacecraft file, or if Fermi operated in pointed mode for a while.
            if exposure[i] != 0.:
                file.write(str(time[i])+"\t"+str(timeMjd[i])+"\t"+str(flux[i])+"\t"+str(fluxErr[i])+"\n")
        file.close()


    def getBAT(self):
        import urllib2

        # daily fits example url:
        # http://swift.gsfc.nasa.gov/docs/swift/results/transients/CygX-3.lc.fits

        # Some sources need replacement names to match the BAT names
        urls={
                '4U1907+09':'H1907+097',
                '1FGLJ1018.6-5856':'0FGLJ1018.2-5858',
                'H1743-322': 'IGRJ17464-3213',
                'V4641Sgr':'SAXJ1819.3-2525',
                '1E1841-04.5':'Kes73',
            }

        # Remove '+', add file ending
        if urls.has_key(self.src):
            file=urls[self.src].replace('+','p')+".lc.fits"
        else:
            file=self.src.replace('+','p')+".lc.fits"
        urlprefix="http://swift.gsfc.nasa.gov/docs/swift/results/transients/"

        # lc files can be in a weak/ subdir for weak sources, we try both
        try:
            baturl=urlprefix+file
            webfile=urllib2.urlopen(baturl)
        except (urllib2.HTTPError, urllib2.URLError) as e:
            try:
                baturl=urlprefix+'weak/'+file
                webfile=urllib2.urlopen(baturl)
            except (urllib2.HTTPError, urllib2.URLError) as e:
                return False,None

        # save lc to local file
        localfile=open(file,'w')
        localfile.write(webfile.read())
        webfile.close()
        localfile.close()
        # read local file with pyfits into batlc
        batfits=pyfits.open(file)
        batlc=array(batfits[1].data)
        batfits.close()
        # delete local file
        os.unlink(file)

        return True,batlc


        
    def createLCfig(self):
        """
        Create a PNG figure with the light curve of a given source. Any existing PNG file is overwritten !
        """

        # Read the .dat LC file
        if self.daily:
            infile=self.workDir+'/'+str(self.src)+'_daily_lc.dat'
            outfig=self.workDir+'/'+str(self.src)+'_daily_lc.png'
            infileLongTimeBin=self.workDir+'/'+str(self.src)+'_lc.dat'
            duration = 1. # duration of a time bin, in days
        else:
            infile=self.workDir+'/'+str(self.src)+'_lc.dat'
            outfig=self.workDir+'/'+str(self.src)+'_lc.png'
            duration = self.longtimebin # duration of a time bin, in days

        data    = asciidata.open(infile)
        # the times are already read as MJD, cf createDAT function.
        timelc  = data[1].tonumpy()
        flux    = data[2].tonumpy()
        fluxErr = data[3].tonumpy()

        if self.daily:
            dataLongTimeBin    = asciidata.open(infileLongTimeBin)
            # the times are already read as MJD, cf createDAT function.
            timeLongTimeBin    = dataLongTimeBin[1].tonumpy()
            fluxLongTimeBin    = dataLongTimeBin[2].tonumpy()
            fluxErrLongTimeBin = dataLongTimeBin[3].tonumpy()
            durationLongTimeBin= self.longtimebin # duration of a time bin, in days

        # Download Swift/BAT data if available
        # xray is boolean flag indicating that X-ray BAT data is available
        xray,batlc=self.getBAT()

        # Redefine the trigger threshold if withhistory=True
        if self.withhistory:
            (fluxAverage,fluxRMS) = self.dynamicalTrigger()

        fig=figure()

        if xray:
            ax    = fig.add_subplot(211)
            axbat = fig.add_subplot(212,sharex=ax)
        else:
            ax = fig.add_subplot(111)

        if self.fglName is not None:
            title=str(self.src)+', '+str(self.fglName).replace('_2FGLJ','2FGL J').replace('3FGLJ','3FGL J')
        else:
            title=str(self.src)+', no known 3FGL counterpart'
        if str(self.z)=='--': # this is the result of the conversion of None from asciidata to numpy to str
            title=title+' (z unknown)'
        else:
            title=title+' (z='+str(self.z)+')'

        ax.set_title(title)

        
        # Force the y-axis ticks to use 1e-6 as a base exponent
        ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: ('%.2f')%(x*1e6)))
        ax.set_ylabel('F (%.0f MeV-%.0f GeV) (x 10^-6 ph cm^-2 s^-1)'%(self.emin,self.emax/1000.),size='x-small')

        day=24.*60.*60.

        ## Make the x-axis ticks shifted by some value
        ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: '%.0f'%(x-TOFFSET)))
        ax.set_xlabel('MJD-'+str(TOFFSET))
        #ax.set_xlabel('MJD')

        # Plot the Fermi/LAT light curve
        if self.daily:
            # Also plot the long time-binned light curve
            ax.errorbar(x=timelc, xerr=duration/2., y=flux, yerr=fluxErr, fmt='ro')
            ax.errorbar(x=timeLongTimeBin, xerr=durationLongTimeBin/2., y=fluxLongTimeBin, yerr=fluxErrLongTimeBin, fmt='bo')
            # The last plot called is on top of the others in matplotlib (are you sure ???). Here, we want the long time-binned LC on top, for visibility.
        else:
            ax.errorbar(x=timelc, xerr=duration/2., y=flux, yerr=fluxErr, fmt='bo')

        # TODO if threshold is dynamic and it is computed using the last flux
        # measurement uncertainty, threshold will not be a given horizontal
        # line, but depend on each point. In that case, indicate in different
        # color those flux measurements that trigger

        # Plot a line at the threshold value
        ax.axhline(y=self.threshold,linewidth=3,linestyle='--',color='r')
        if self.withhistory:
            ax.axhline(y=fluxAverage,linewidth=1,linestyle='-',color='b')
            ax.axhline(y=fluxAverage+fluxRMS,linewidth=1,linestyle='--',color='b')
            ax.axhline(y=fluxAverage-fluxRMS,linewidth=1,linestyle='--',color='b')

        # Plot a line at flux=0, for visibility/readibility
        ax.axhline(y=0.,color='k')

        # Add a label for the creation date of this figure (highly inspired from Marcus Hauser's ADRAS/ATOM pipeline)
        # x,y in relative 0-1 coords in figure
        figtext(0.98, 0.95,
                'plot creation date: %s (UTC)'%(time.strftime("%a, %d %b %Y %H:%M:%S", time.gmtime())),
                horizontalalignment="right",
                rotation='vertical',
                size='xx-small'
                )

        # Plot Swift/BAT lightcurve
        if xray:
            axbat.errorbar(batlc['TIME']+0.5,batlc['RATE'],batlc['ERROR'],fmt="none",capsize=0,elinewidth=1,ecolor='b',color='b')
            axbat.set_xlabel('MJD-'+str(TOFFSET))
            #axbat.set_xlabel('MJD')
            axbat.set_ylabel('F (15-50 keV) (count cm^-2 s^-1)',size='x-small')
            try:
                axbat.set_xlim(xmin=timelc[0]-duration/2.-1.,xmax=timelc[-1:]+duration/2.+1.)
                axbat.set_ylim(ymin=0.)
            except:
                pass

        # Need to zoom in or not, at the very end, after any call to other matplotlib functions
        NEEDTOZOOMIN=False
        for i in range(len(flux)):
            if fluxErr[i] > 5.*flux[i]:
                NEEDTOZOOMIN=True
        if NEEDTOZOOMIN:
            maxy=1.5*max(flux)
            if maxy>self.threshold:
                ax.set_ylim(ymin=-1.e-7,ymax=maxy)
            else:
                ax.set_ylim(ymin=-1.e-7,ymax=self.threshold)
        
        # Don't show the figure in batch mode
        if not BATCH:
            show()
        # Save the figure
        fig.savefig(outfig)


    def createEnergyTimeFig(self, eThresh=1.e2):
        """
        Create a PNG figure with the energy vs time of a given source, above eThresh MeV. Any existing PNG file is overwritten !
        """

        # Read the GTI FITS file
        infile=self.workDir+'/'+str(self.src)+'_gti.fits'
        outfig=self.workDir+'/'+str(self.src)+'_energyTime.png'

        hdu   = pyfits.open(infile)
        data  = hdu[1].data
        mask  = data.field('ENERGY')>eThresh
        datac = data[mask]
        if not datac.size:
            print '[%s] \033[92mWARNING Empty energy vs time plot above %0.f GeV\033[0m' % (self.src, eThresh/1.e3)
            return

        t=met2mjd(datac['TIME'])
        e=datac['ENERGY']

        fig=figure()
        ax = fig.add_subplot(111)

        if self.fglName is not None:
            title=str(self.src)+', '+str(self.fglName).replace('_2FGLJ','2FGL J').replace('3FGLJ','3FGL J')
        else:
            title=str(self.src)+', no known 3FGL counterpart'
        if str(self.z)=='--': # this is the result of the conversion of None from asciidata to numpy to str
            title=title+' (z unknown)'
        else:
            title=title+' (z='+str(self.z)+')'

        ax.set_title(title)

        ylabel = 'Energy (MeV)'
        if eThresh > self.emin:
            ylable += ' -- only data above %.1f GeV are shown' % (eThresh/1.e3)
        ax.set_ylabel(ylabel, size='x-small')

        ## Make the x-axis ticks shifted by some value
        ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: '%.0f'%(x-TOFFSET)))
        ax.set_xlabel('MJD-'+str(TOFFSET))
        try:
            ax.set_xlim(xmin=t[0]-1.,xmax=t[-1:]+1.)
        except:
            pass

        # Plot the energy vs time distribution
        try:
            # cf. http://stackoverflow.com/questions/20105364/how-can-i-make-a-scatter-plot-colored-by-density-in-matplotlib
            from scipy.stats import gaussian_kde
            xy = vstack([t, e])
            z = gaussian_kde(xy)(xy)
            # Re-normalize the density
            z = z/max(z)
            idx = z.argsort()
            t, e, z = t[idx], e[idx], z[idx]
            pcm = ax.scatter(t, e, c=z, s=100, edgecolor='')
            cbar = colorbar(pcm, ax=ax)
            cbar.set_label('Kernel-density estimates (arb. unit)', rotation=90)
        except ImportError:
            ax.plot(t, e,  'bo')
        ax.set_yscale('log')

        # Add a label for the creation date of this figure (highly inspired from Marcus Hauser's ADRAS/ATOM pipeline)
        # x,y in relative 0-1 coords in figure
        figtext(0.98, 0.95,
                'plot creation date: %s (UTC)'%(time.strftime("%a, %d %b %Y %H:%M:%S", time.gmtime())),
                horizontalalignment="right",
                rotation='vertical',
                size='xx-small'
                )
        
        # Don't show the figure in batch mode
        if not BATCH:
            show()
        # Save the figure
        fig.savefig(outfig)


    def zaAtCulmination(self):
        """
        Returns the zenith angle of a source at culmination, for the provided site.
        """
        return abs(self.dec-self.siteLat)


    def is_visible(self):
        '''
        Check whether the current source is visible at the site provided.
        '''
        
        # Define site for pyephem
        site    = ephem.Observer()
        site.pressure = 0
        astroHorizon = '-18:00' # astronomical twilight
        civilHorizon = '-0:34'
        site.horizon = astroHorizon
        site.lon = astCoords.decimal2dms(self.siteLon, delimiter=':')
        site.lat = astCoords.decimal2dms(self.siteLat, delimiter=':')
        site.elev = self.siteAlt

        # If input z is None, make it believe it is 0, otherwise msk crashes:
        if str(self.z)=='--': # this is the result of the conversion of None from asciidata to numpy to str
            z = 0.
        else:
            z = self.z

        # We also want the max allowed ZA for the given z of the source
        maxz = array(self.maxz)
        maxZA = array(self.maxZA)
        if z > max(maxz):
            thismaxZA = min(maxZA)
        else:
            msk = where(z<maxz)
            # Get the first item in the mask, to get the corresponding ZA:
            thismaxZA = maxZA[msk[0][0]]

        # Convert ZA to Alt
        thisminAlt=abs(90.-thismaxZA)
        
        ephemSrc = ephem.FixedBody()
        ephemSrc._ra=astCoords.decimal2hms(self.ra,delimiter=':')
        ephemSrc._dec=astCoords.decimal2dms(self.dec,delimiter=':')
        
        visibleFlag=False

        zaAtCulmin = self.zaAtCulmination()
        if zaAtCulmin>90.:
            # the source is basically NEVER visible at the site
            print '[%s] \033[91mNEVER above horizon at the site, consider discarding this source from your source list...\033[0m' % self.src
            return False

        if thismaxZA<zaAtCulmin:
            # the source is never above maxZA set by 2D mask on Dec/z
            print '[%s]\033[91m Never above allowed max ZA, consider relaxing the Dec/z cuts or discarding this source from your source list...\033[0m' % self.src
            return False
        
        
        # All times are handled here in UTC (pyEphem only uses UTC)
        now      = datetime.datetime.utcnow()
        # tomorrow = now + datetime.timedelta(days=1)
        
        site.date      = now
        sun            = ephem.Sun()
        nextSunset     = site.next_setting(sun)
        nextSunrise    = site.next_rising(sun)
        # The Moon just needs to be below the horizon, not below astronomical twilight angle
        site.horizon   = civilHorizon
        moon           = ephem.Moon()
        nextMoonset    = site.next_setting(moon)
        nextMoonrise   = site.next_rising(moon)
        site.horizon   = astroHorizon
        # so far, so good. All of this is OK if we execute the program during day time.

        # However, if the program is run during dark time, we should look at the ephemerids of next night (not current night):
        if nextSunrise < nextSunset:
            if VERBOSE:
                print "INFO: looking at visibility for tomorrow"
            # we just put the current time at next sunrise + 10 min., to be sure to fall on tomorrow's morning day time
            site.date = nextSunrise.datetime() + datetime.timedelta(minutes=10)
            nextSunset    = site.next_setting(sun)
            nextSunrise   = site.next_rising(sun)
            site.horizon = civilHorizon
            nextMoonset   = site.next_setting(moon)
            nextMoonrise  = site.next_rising(moon)
            site.horizon = astroHorizon

        ephemSrc.compute(site)
        srcTransitTime = site.next_transit(ephemSrc)
        
        site.date=srcTransitTime
        ephemSrc.compute(site)
        srcAltAtTransit=astCoords.dms2decimal(ephemSrc.alt,delimiter=':')
        
        # If srcAltAtTransit is below thisminAlt, the source is just not correctly visible and we stop here
        if srcAltAtTransit < thisminAlt:
            return False

        # Compute start and end of darkness time
        if nextMoonset > nextSunset and nextMoonset < nextSunrise:
            beginDarkness=nextMoonset
        else:
            beginDarkness=nextSunset

        if nextMoonrise < nextSunrise and nextMoonrise > nextSunset:
            endDarkness=nextMoonrise
        else:
            endDarkness=nextSunrise

        if DEBUG:
            darknessDuration = endDarkness-beginDarkness
            print "DEBUG: darkness begin=%s" % beginDarkness
            print "DEBUG: darkness ends=%s" % endDarkness
            print "DEBUG: darkness duration=%s minutes" % (darknessDuration*24.*60.)

        site.date=beginDarkness
        ephemSrc.compute(site)
        srcAltAtStartDarkTime=astCoords.dms2decimal(ephemSrc.alt,delimiter=':')
        
        site.date=endDarkness
        ephemSrc.compute(site)
        srcAltAtEndDarkTime=astCoords.dms2decimal(ephemSrc.alt,delimiter=':')
        
        # check if source is visible, above minAlt, during this night
        if (srcTransitTime > beginDarkness and srcTransitTime < endDarkness and srcAltAtTransit > thisminAlt) or srcAltAtStartDarkTime > thisminAlt or srcAltAtEndDarkTime > thisminAlt:
            visibleFlag=True

        if VERBOSE:
            print "INFO: is_visible: "+str(visibleFlag)
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
        #                        z         ZA@culmination
        grid = array(zip(self.maxz,self.maxZA))
    
        try:
            # import Kapteyn module for WCS
            from kapteyn import wcs
        except:
            print "ERROR cutsTrigger: Can't import kapteyn.wcs module."
            sys.exit(1)


        
        # ******** FOR THE FUTURE *********
        ## Equatorial coordinates -> Galactic coordinates transformation
        #equat2gal = wcs.Transformation((wcs.equatorial, wcs.fk5,'J2000.0'),wcs.galactic)
        ## WCS module needs a numpy matrix for coord transformation
        #srcradec=matrix(([ra],[dec]))
        #srcGalCoord=equat2gal(srcradec)
        ## Retrieve the Galactic latitude
        #srcGalLat=srcGalCoord[(1,0)]
        ## To be used later, for an additional cut on Gal Lat ?!
        
        zaAtCulmin = self.zaAtCulmination()

        # If input z is None, make it believe it is 0, otherwise msk crashes:
        if str(self.z)=='--': # this is the result of the conversion of None from asciidata to numpy to str
            z = 0.
        else:
            z = self.z

        # Mask on both (z, ZA at culmin)
        #          z column               ZA column
        msk = (z<=grid[:,0])&(zaAtCulmin<=grid[:,1])

        # Assess whether the source is currently visible at the provided site
        if self.checkVisibility == 'True':
            self.visible = self.is_visible()
        else:
            # The source is assumed to be visible in any case, i.e. we don't care about its visibility status at the provided site to send a potential alert
            self.visible = True

        # if the mask has at least one 'True' element, we should send an alert
        if True in msk and self.visible:
            # print 'An alert should be triggered !'
            return False
        else:
            # print 'No alert triggered'
            return True


    def dynamicalTrigger(self):
        '''
        If long-term data are available for a source, dynamically computes a flux trigger threshold based on the flux history of the source. Otherwise, fall back with default fixed trigger threshold.

        @return (fluxAverage,fluxRMS)
        @rtype tuple
        '''
        
        # Read the longterm .dat LC file
        infile = self.baseOutDir+'/longTerm/merged/'+str(self.src)+'_lc.dat'
        try:
            data    = asciidata.open(infile)
        except IOError:
            print '[%s] \033[95m* Long term data file unavailable for source %s\033[0m' % (self.src, self.src)
            # Falling back to default fixed trigger threshold
            self.withhistory=False
            return (False,False)

        flux        = data[2].tonumpy()
        fluxErr     = data[3].tonumpy()
        if VERBOSE:
            from uncertainties import unumpy as unp
            print('INFO: The long-term flux average is ', unp.uarray(flux, fluxErr).mean())

        # weighted average of the historical fluxes, weighted by their errors
        fluxAverage = average(flux, weights=1./fluxErr)
        fluxRMS     = std(flux, dtype=float64)

        # Dynamically redefine the flux trigger threshold, using a 2-level criteria depending on whether we are currently looking at short- or long-term data
        if self.daily:
            self.threshold = fluxAverage + self.sigma*fluxRMS
        else:
            self.threshold = fluxAverage + self.sigmaLT*fluxRMS

        return (fluxAverage,fluxRMS)
        

    def Triggered(self):
        '''
        Has the source fulfilled the trigger conditions ?

        @return True
        @rtype bool
        '''

        # Read the light curve file
        if self.daily:
            infile  = self.workDir+'/'+str(self.src)+'_daily_lc.dat'
            self.pngFig=self.workDir+'/'+str(self.src)+'_daily_lc.png'

            # Also take a look in the long time-binned data
            infileLongTimeBin=self.workDir+'/'+str(self.src)+'_lc.dat'
            dataLongTimeBin=asciidata.open(infileLongTimeBin)
            timeLongTimeBin=dataLongTimeBin[0].tonumpy()
            fluxLongTimeBin=dataLongTimeBin[2].tonumpy()
            fluxErrLongTimeBin=dataLongTimeBin[3].tonumpy()
            # Catch the last flux point
            self.lastTimeLongTimeBin=timeLongTimeBin[-1:]
            self.lastFluxLongTimeBin=fluxLongTimeBin[-1:]
            self.lastFluxErrLongTimeBin=fluxErrLongTimeBin[-1:]

            # Get the arrival time of the last photon analysed
            photonfileLongTimeBin            = self.workDir+'/'+str(self.src)+'_gti.fits'
            photonsLongTimeBin               = pyfits.open(photonfileLongTimeBin)
            photonsLongTimeBinTime           = photonsLongTimeBin[1].data.field('TIME')
            self.arrivalTimeLastPhotonLongTimeBin = photonsLongTimeBinTime[-1:]

            photonfile                  = self.workDir+'/'+str(self.src)+'_daily_gti.fits'
            photons                     = pyfits.open(photonfile)
            photonsTime                 = photons[1].data.field('TIME')
            self.arrivalTimeLastPhoton       = photonsTime[-1:]
        else:
            infile  = self.workDir+'/'+str(self.src)+'_lc.dat'
            self.pngFig=self.workDir+'/'+str(self.src)+'_lc.png'

            photonfile            = self.workDir+'/'+str(self.src)+'_gti.fits'
            photons               = pyfits.open(photonfile)
            photonsTime           = photons[1].data.field('TIME')
            self.arrivalTimeLastPhoton = photonsTime[-1:]
        data    = asciidata.open(infile)
        time    = data[0].tonumpy()
        flux    = data[2].tonumpy()
        fluxErr = data[3].tonumpy()

        # Catch the last flux point
        self.lastTime    = time[-1:]
        self.lastFlux    = flux[-1:]
        self.lastFluxErr = fluxErr[-1:]

        self.energyTimeFig=self.workDir+'/'+str(self.src)+'_energyTime.png'
            
        if DEBUG:
            print 'DEBUG %s, threshold=%g, lastFlux=%g, lastFluxErr=%g' % (self.src,self.threshold,self.lastFlux,self.lastFluxErr)

        # Do we kill potential trigger due to (ra, dec, z) cut ?
        self.triggerkilled = self.killTrigger()

        # Assess whether flux is above threshold, looking at the last flux point
        #if self.daily:
        #    if ((self.lastFlux - self.lastFluxErr) >= self.threshold or (self.lastFluxLongTimeBin - self.lastFluxErrLongTimeBin) >= self.threshold):
        #        self.active=True
        #    else:
        #        self.active=False
        #else:
        #    if (self.lastFlux - self.lastFluxErr) >= self.threshold:
        #        self.active=True
        #    else:
        #        self.active=False

        if (self.lastFlux - self.lastFluxErr) >= self.threshold:
            self.active=True
        else:
            self.active=False

        # Combine killTrigger and flux above threshold criteria
        if (not self.triggerkilled and self.active) or FORCE_ALERT:
            SENDALERT = True
        else:
            SENDALERT = False

        if VERBOSE:
            print "INFO: triggerkilled="+str(self.triggerkilled)
            print "INFO: active="+str(self.active)
            print "INFO: visible="+str(self.visible)
            print "INFO: SENDALERT="+str(SENDALERT)

        if DEBUG:
            print "DEBUG %s, dec=%f, z=%f, maxZA=[%s], maxz=[%s], triggerkilled=%s, sendalert=%s" % (str(self.src),self.dec,self.z,', '.join(map(str,self.maxZA)),', '.join(map(str,self.maxz)),self.triggerkilled,SENDALERT)

        return SENDALERT


    def sendAlert(self,nomailall=False,sendmail=False):
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
            from email.MIMEImage import MIMEImage
            from email.MIMEMultipart import MIMEMultipart
            from email.MIMEText import MIMEText
            from email.MIMEBase import MIMEBase
            from email import Encoders

        except:
            print "ERROR sendAlert: Can't import mail modules."
            sys.exit(1)


        SENDALERT = self.Triggered()


        # If trigger condition is met, we send a mail
        if SENDALERT and sendmail:
            # Create the container email message.
            msg = MIMEMultipart()
            sender = self.mailSender

            fhlName=self.search2FHLcounterpart()
            if fhlName is not None:
                fhlmessage="2FHL counterpart is %s" % fhlName
            else:
                fhlmessage="No 2FHL counterpart found"

            fglName=self.search3FGLcounterpart()
            if fglName is not None:
                fglmessage="3FGL counterpart is %s" % fglName
            else:
                fglmessage="No 3FGL counterpart found"

            # To whom the mail should be sent (cf. __init__ function of the class)
            if not nomailall:
                recipient = self.usualRecipients
                msg['Subject'] = '[FLaapLUC] Fermi/LAT flare alert on %s [2FHL counterpart: %s]' % (self.src,fhlName)
            else:
                recipient = self.testRecipients
                msg['Subject'] = '[FLaapLUC TEST MAIL] Fermi/LAT flare alert on %s [2FHL counterpart: %s]' % (self.src, fhlName)

            msg['From'] = sender
            COMMASPACE = ', '
            msg['To'] =COMMASPACE.join( recipient )
            msg.preamble = 'You will not see this in a MIME-aware mail reader.\n'
            # Guarantees the message ends in a newline
            msg.epilogue = ''
                       
            mailtext="""
     FLaapLUC (Fermi/LAT automatic aperture photometry Light C<->Urve) report

     *** The Fermi/LAT flux (%.0f MeV-%.0f GeV) of %s (%s, %s) exceeds the trigger threshold of %.2g ph cm^-2 s^-1 ***

     """%(self.emin,self.emax/1000.,self.src,fhlmessage,fglmessage,self.threshold)

            if self.daily:
                mailtext=mailtext+"""

     The last daily-binned flux is:        %.2g +/- %.2g ph cm^-2 s^-1, centred on MET %.0f (MJD %.5f, i.e. %s) (arrival time of last photon analysed: MET %.0f, MJD %.5f, %s)
     and the last %.0f-day binned flux is: %.2g +/- %.2g ph cm^-2 s^-1, centred on MET %.0f (MJD %.5f, i.e. %s) (arrival time of last photon analysed: MET %.0f, MJD %.5f, %s)

""" % (self.lastFlux,
       self.lastFluxErr,
       self.lastTime, met2mjd(self.lastTime), str(mjd2gd(met2mjd(self.lastTime))),
       self.arrivalTimeLastPhoton, met2mjd(self.arrivalTimeLastPhoton), str(mjd2gd(met2mjd(self.arrivalTimeLastPhoton))),
       self.longtimebin,
       self.lastFluxLongTimeBin,
       self.lastFluxErrLongTimeBin,
       self.lastTimeLongTimeBin, met2mjd(self.lastTimeLongTimeBin), str(mjd2gd(met2mjd(self.lastTimeLongTimeBin))),
       self.arrivalTimeLastPhotonLongTimeBin, met2mjd(self.arrivalTimeLastPhotonLongTimeBin), str(mjd2gd(met2mjd(self.arrivalTimeLastPhotonLongTimeBin))))
                mailtext=mailtext+"The most recent lightcurve (%.0f-day binned in red, and %.0f-day binned in blue) is attached."%(self.tbin/24./60./60.,self.longtimebin)
            else:
                mailtext=mailtext+"""

     The last %.0f-day binned flux is:      %.2g +/- %.2g ph cm^-2 s^-1, centred on MET %.0f (MJD %.5f, i.e. %s) (arrival time of last photon analysed: %.0f, MJD %.5f, %s)

""" % (self.longtimebin,
       self.lastFlux,
       self.lastFluxErr,
       self.lastTime, met2mjd(self.lastTime), str(mjd2gd(met2mjd(self.lastTime))),
       self.arrivalTimeLastPhoton, met2mjd(self.arrivalTimeLastPhoton), str(mjd2gd(met2mjd(self.arrivalTimeLastPhoton))))
                mailtext=mailtext+"The most recent lightcurve (%.0f-day binned) is attached."%(self.tbin/24./60./60.)

            if self.launchLikeAna == 'True':
                mailtext=mailtext+"""

     *NOTE*: a likelihood analysis has been automatically launched at CCIN2P3 for the time interval corresponding to the last measurement (MET %i - MET %i). Contact Jean-Philippe Lenain (jlenain@in2p3.fr) to know the outcome.

"""%(self.tstop-(self.longtimebin*24.*3600.), self.tstop)

            if FLAGASSUMEDGAMMA is True:
                mailtext=mailtext+"""

     *WARNING*: The source %s is not found in the 3FGL catalogue, its photon index is thus assumed to be %.2f for the light curve computation.
""" % (self.src,ASSUMEDGAMMA)


            mailtext=mailtext+"""

     All available data can be found on the web site
     %s

     *Disclaimer*: Be careful, though, that these light curves are not computed using the usual, clean, standard (un)binned likelihood procedure one should normally use for a good quality, publication-ready result. Those reported here only rely on a "quick & dirty" aperture photometric analysis (cf. e.g. http://fermi.gsfc.nasa.gov/ssc/data/analysis/scitools/aperture_photometry.html), which basically assumes that the data set, within 1 degree around the source, is background-free.

      For more information, please contact J.-P. Lenain <jlenain@in2p3.fr>.

      Cheers,
      FLaapLUC.
""" %(self.url)
 
            txt = MIMEText(mailtext)
            msg.attach(txt)
            
            # Attach the figures
            for fig  in [self.pngFig, self.energyTimeFig]:
                try:
                    # Open the files in binary mode.  Let the MIMEImage class automatically guess the specific image type.
                    fp = open(fig, 'rb')
                    # img = MIMEImage(fp.read(), name=os.path.basename(fig))
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

            print "\033[94m*** Alert sent for %s\033[0m" % self.src

            return True
        else:
            return False



    def sendErrorMail(self,mailall=False):
        '''
        Send an error mail message.

        @param mailall: boolean, should the mail be sent to UsualRecipients or TestRecipients ?
        '''

        # Import modules
        try:
            # Import smtplib to send mails
            import smtplib

            # Here are the email package modules we'll need
            from email.MIMEImage import MIMEImage
            from email.MIMEMultipart import MIMEMultipart
            from email.MIMEText import MIMEText

        except:
            print "ERROR sendAlert: Can't import mail modules."
            sys.exit(1)

        # Create the container email message.
        msg = MIMEMultipart()
        msg['Subject'] = '[FLaapLUC] ERROR Fermi/LAT automatic light curve'
        sender = self.mailSender
        recipient = self.testRecipients

        msg['From'] = sender
        COMMASPACE = ', '
        msg['To'] =COMMASPACE.join( recipient )
        msg.preamble = 'You will not see this in a MIME-aware mail reader.\n'
        # Guarantees the message ends in a newline
        msg.epilogue = ''
            
        mailtext="""
        The Fermi automatic light curve pipeline has been automatically aborted, because new Fermi data could not be properly downloaded. Maybe the NASA servers are down.

        For more information, please contact J.-P. Lenain <jlenain@in2p3.fr>.

        Cheers,
        FLaapLUC.
"""
 
        txt = MIMEText(mailtext)
        msg.attach(txt)
            
        # Send the email via our own SMTP server.
        s = smtplib.SMTP()
        s.set_debuglevel(0)
        s.connect()
        s.sendmail(sender, recipient, msg.as_string())
        s.quit()

        return True

    def search3FGLcounterpart(self):
        """
        Search the 3FGL name of a 2FGL source name
        """
        if self.fglName is not None:
            if "3FGL" in self.fglName:
                return self.fglName.replace('_3FGLJ','3FGL J').replace('3FGLJ','3FGL J')

            cat3FGLfile = self.catalogFile.replace('gll_psc_v08','gll_psc_v16')
            hdulist = pyfits.open(cat3FGLfile)
            cat=hdulist[1].data
            if DEBUG:
                print 'DEBUG: 2FGL name is %s' % self.fglName.replace('_2FGLJ','2FGL J').replace('2FGLJ','2FGL J')

            found=False
            for stuff in cat:
                if stuff.field('2FGL_Name') == self.fglName.replace('_2FGLJ','2FGL J').replace('2FGLJ','2FGL J'):
                    threefglName=stuff.field('Source_Name')
                    if VERBOSE:
                        print 'INFO: Found the 3FGL counterpart of %s: %s' % (self.fglName,threefglName)
                    found=True
                    break

            if not found:
                threefglName=None
                if VERBOSE:
                    print 'INFO: No 3FGL counterpart found for %s' % self.fglName

            hdulist.close()
            return threefglName
        else:
            return None


    def search2FHLcounterpart(self):
        """
        Search the 2FHL name of a 2FGL or a 3FGL source name
        """
        if self.fglName is not None:
            if "2FHL" in self.fglName:
                return self.fglName.replace('_2FHLJ','2FHL J').replace('2FHLJ','2FHL J')

            cat2FHLfile = self.catalogFile.replace('/3FGL/','/2FHL/').replace('psc_v08','psch_v08').replace('psc_v16','psch_v08')
            try:
                hdulist = pyfits.open(cat2FHLfile)
            except IOErrror:
                if VERBOSE:
                    print 'INFO: 2FHL catalog file not found'
                return None
            cat=hdulist[1].data

            found=False
            threefglName=self.search3FGLcounterpart()
            for stuff in cat:
                if stuff.field('3FGL_Name') == self.fglName.replace('_3FGLJ','3FGL J').replace('3FGLJ','3FGL J') or stuff.field('3FGL_Name') == str(threefglName).replace('3FGLJ','3FGL J'):
                    fhlName=stuff.field('Source_Name')
                    if VERBOSE:
                        print 'INFO: Found the 2FHL counterpart of %s: %s' % (self.fglName,fhlName)
                    found=True
                    break

            if not found:
                fhlName=None
                if VERBOSE:
                    print 'INFO: No 2FHL counterpart found for %s' % self.fglName

            hdulist.close()
            return fhlName
        else:
            return None
        

    def launchLikelihoodAnalysis(self, nomailall=False):
        """
        Launch a clean likelihood analysis in Lyon
        """
        
        srcDir=self.src+'_FLaapLUC_'+str(datetime.date.today().strftime('%Y%m%d'))
        anaDir=os.getenv('FERMIUSER')+'/'+srcDir
        # If this analysis has already been launched, we abort here
        if os.path.isdir(anaDir):
            return False
        os.makedirs(anaDir)
        if self.fglName is not None and '2FGL' in self.fglName:
            threefglName=self.search3FGLcounterpart()
        elif self.fglName is not None and '3FGL' in self.fglName:
            threefglName=self.fglName
        else:
            threefglName=None
        if threefglName is not None:
            fglNameFile=anaDir+'/FermiName.txt'
            file=open(fglNameFile,'w')
            file.write(threefglName.replace('3FGLJ','3FGL J'))
            file.close()
        if VERBOSE:
            print 'INFO: launchLikelihoodAnalysis: 3FGL name is %s' % threefglName
        srcSelectFile=anaDir+'/source_selection.txt'
        srcSelect=open(srcSelectFile,'w')
        srcSelect.write("""Search Center (RA,Dec)  =       (%f,%f)
Radius  =       10 degrees
Start Time (MET)        =       %i seconds (MJD%f)
Stop Time (MET) =       %i seconds (MJD%f)
Minimum Energy  =       %i MeV
Maximum Energy  =       %i MeV
""" % (self.ra, self.dec,
       self.tstop-(self.longtimebin*24.*3600.), met2mjd(self.tstop-(self.longtimebin*24.*3600.)),
       self.tstop, met2mjd(self.tstop),
       int(self.emin),
       int(self.emax)))
        srcSelect.close()

        photonFile=anaDir+'/photon.list'
        photonList=open(photonFile,'w')
        photonList.write(str(self.lastAllskyFile))
        photonList.close()

        def writescript(emin,options=""):
            script = """#!/usr/local/bin/bash
export FERMIUSER=/sps/hess/users/lpnhe/jlenain/fermi
SRC=%s
tmpfile=`mktemp --tmpdir=\$FERMIUSER/tmp myLATanalysis_\${SRC}_XXXX.sh`
        
cat > $tmpfile << EOF
#!/usr/local/bin/bash
#$ -S /bin/bash
#$ -N myLATanalysis
#$ -e /sps/hess/lpnhe/jlenain/fermi/log
#$ -o /sps/hess/lpnhe/jlenain/fermi/log
#$ -j yes
#$ -l ct=5:00:00
#$ -l mem_free=2048M
#$ -l fsize=2048M
#$ -l xrootd=0
#$ -l sps=1
#$ -l os=sl6
#$ -p 0
#$ -P P_hess -cwd

SRC=$SRC
export FERMI_DIR=/sps/hess/users/lpnhe/jlenain/local/fermi/ScienceTools-v10r0p5-fssc-20150518-x86_64-unknown-linux-gnu-libc2.12/x86_64-unknown-linux-gnu-libc2.12
export FERMIUSER=/sps/hess/users/lpnhe/jlenain/fermi
source $FERMI_DIR/fermi-init.sh

LOG=${tmpfile%%.sh}.log
EMIN=%i
EMAX=%i
BOUNDARY=\`date +%%s | md5sum | cut -b1-32\`

cat > \$LOG << EOM
From: %s
Subject: [FLaapLUC likelihood] [\$SRC] Fermi likelihood analysis report (\$EMIN MeV - \$EMAX MeV)
To: %s
Content-Type: multipart/mixed; boundary="\$BOUNDARY"

This is a MIME formatted message.  If you see this text it means that your email software does not support MIME formatted messages.

--\$BOUNDARY
Content-Type: text/plain; charset=ISO-8859-1; format=flowed
Content-Transfer-Encoding: 7bit
Content-Disposition: inline

Source:  \$SRC
E_min:   \$EMIN MeV
E_max:   \$EMAX MeV

A likelihood analysis was automatically launched by FLaapLUC after issuing an alert on %s, which log can be found below. You'll most probably be interested mainly in the lines beginning with 'INFO FINAL', giving the fitted flux and photon index for the source.

Note that if the source has a known redshift, the Fermi/LAT spectrum computed from this likelihood analysis is extrapolated to very high energies, assuming a somewhat ad-hoc spectral break of Delta Gamma=1 at 100 GeV to account for intrinsic spectral curvature, absorbed by the EBL following the EBL model of Franceschini et al. (2008), and compared to the H.E.S.S. II differential sensitivity for 20h of observations (Model++ Stereo Std cuts, Zenith Angle of 18 deg.). The resulting plot should be attached to this mail, otherwise please contact Jean-Philippe Lenain (<jlenain@in2p3.fr>).

*********************************************************************************

EOM

## Define the energy range \$ERANGE string
a=\$EMIN
b=\$EMAX
c=1000.0 # Break energy of 1 GeV, between 'MeV' label and 'GeV' label
case \`echo "a=\$a;c=\$c;r=-1;if(a==c)r=0;if(a>c)r=1;r"|bc\` in
    0) ERANGEMIN="1GeV";;
    1) a=\`echo "\$a/1000.0"|bc\`
	ERANGEMIN=\${a}GeV;;
    *) a=\`echo "\$a/1.0"|bc\`
	ERANGEMIN=\${a}MeV;;
esac
case \`echo "b=\$b;c=\$c;r=-1;if(b==c)r=0;if(b>c)r=1;r"|bc\` in
    0) ERANGEMAX="1GeV";;
    1) b=\`echo "\$b/1000.0"|bc\`
	ERANGEMAX=\${b}GeV;;
    *) b=\`echo "\$b/1.0"|bc\`
	ERANGEMAX=\${b}MeV;;
esac
# De-attribute a,b,c
unset a b c
ERANGE=\${ERANGEMIN}_\${ERANGEMAX}

\$FERMIUSER/myLATanalysis.sh -s \${SRC} %s >> \$LOG 2>&1
ATTACH="\${FERMIUSER}/\${SRC}/BINNED/VHEextrapolation/\${SRC}_\${ERANGE}_SED.pdf"
if [ -e \$ATTACH ]; then
    BASENAME=\$(basename \${ATTACH})
    cat >> \$LOG << EOP
--\$BOUNDARY
Content-Type: application/pdf; name="\${BASENAME}"
Content-Transfer-Encoding: base64
Content-Disposition: attachment; filename="\${BASENAME}"

EOP
    cat \$ATTACH | openssl base64 >> \$LOG
fi
ATTACH="\${FERMIUSER}/\${SRC}/BINNED/\${SRC}_\${ERANGE}_sanity_maps.png"
if [ -e \$ATTACH ]; then
    BASENAME=\$(basename \${ATTACH})
    cat >> \$LOG << EOP
--\$BOUNDARY
Content-Type: image/png; name="\${BASENAME}"
Content-Transfer-Encoding: base64
Content-Disposition: attachment; filename="\${BASENAME}"

EOP
    cat \$ATTACH | openssl base64 >> \$LOG
fi
echo "" >> \$LOG
echo "" >> \$LOG
echo "--\${BOUNDARY}--" >> \$LOG
echo "" >> \$LOG
echo "" >> \$LOG
""" % (srcDir, int(emin), int(self.emax), self.mailSender, self.likelihoodRecipients, self.src, options)
            if self.likelihoodRecipients is not None:
                script += "cat \$LOG | sendmail -t"
            script += """
EOF

qsub ${tmpfile}
"""
            return script


        options=""

        if self.fglName is not None:
            options += " -c "
        options += " -a std -q -l -m BINNED -e %i -E %i" % (int(self.emin), int(self.emax))
        if str(self.z)!='--': # this is the result of the conversion of None from asciidata to numpy to str
            # pass the redshift as argument to myLATanalysis to trigger the generation of the EBL absorbed, VHE extrapolation of the Fermi/LAT likelihood spectral results
            options += ' -z %f' % self.z
    
        command=writescript(emin=self.emin,options=options)
        r=os.system(command)

        # If a high energy threshold has been set in the config file, launch another instance of the likelihood analysis using this higher energy threshold:
        if self.likeHighEThresh is not None:
            options=""

            if self.fglName is not None:
                options += " -c "
            options += " -a std -q -l -m BINNED -e %i -E %i" % (int(float(self.likeHighEThresh)), int(self.emax))
            if str(self.z)!='--' and self.z != 0: # this is the result of the conversion of None from asciidata to numpy to str
            # pass the redshift as argument to myLATanalysis to trigger the generation of the EBL absorbed, VHE extrapolation of the Fermi/LAT likelihood spectral results
                options += ' -z %f' % self.z

            command=writescript(emin=float(self.likeHighEThresh),options=options)
            r2=os.system(command)
            return r, r2
            
        return r


def processSrc(mysrc=None,useThresh=False,daily=False,mail=True,longTerm=False,test=False, yearmonth=None, mergelongterm=False,withhistory=False,update=False,configfile='default.cfg',force_daily=False,stopmonth=None):
    """
    Process a given source.
    """

    if mysrc is None:
        print "ERROR Missing input source !"
        sys.exit(1)

    # If we asked for a daily light curve, first make sure that the long time-binned data already exists, otherwise this script will crash, since the daily-binned PNG needs the long time-binned data to be created. No mail alert is sent at this step.
    # We automatically recreate here any missing long time-binned data.
    if daily and not longTerm and not force_daily:
        print "[%s] Daily light curve asked for, I will first process the long time-binned one" % mysrc
        longtermactive, visible=processSrc(mysrc=mysrc,
                                           useThresh=useThresh,
                                           daily=False,
                                           mail=False,
                                           longTerm=longTerm,
                                           yearmonth=yearmonth,
                                           mergelongterm=mergelongterm,
                                           withhistory=withhistory,
                                           update=update,
                                           configfile=configfile,
                                           stopmonth=stopmonth)
        if longtermactive and visible:
            print "[%s] Source %s is active and visible in long time-binned data, processing daily-binned light curve..." % (mysrc, mysrc)
        elif longtermactive and not visible:
            print "[%s] \033[91mSource %s is active but not visible. Daily-binned light curve aborted...\033[0m" % (mysrc, mysrc)
            return False
        elif not longtermactive and visible:
            print "[%s] \033[91mSource %s is visible but not active. Daily-binned light curve aborted...\033[0m" % (mysrc, mysrc)
            return False
        elif not longtermactive and not visible:
            print "[%s] \033[91mSource %s is neither active nor visible. Daily-binned light curve aborted...\033[0m" % (mysrc, mysrc)
            return False
        else:
            print "[%s] \033[91mDaily-binned light curve aborted, for unknown reason...\033[0m" % (mysrc, mysrc)
            return False
    elif force_daily:
        print "[%s] Forcing daily light curve, I will first process the long time-binned one" % mysrc
        longtermactive, visible=processSrc(mysrc=mysrc,
                                           useThresh=useThresh,
                                           daily=False,
                                           mail=False,
                                           longTerm=longTerm,
                                           yearmonth=yearmonth,
                                           mergelongterm=mergelongterm,
                                           withhistory=withhistory,
                                           update=update,
                                           configfile=configfile,
                                           stopmonth=stopmonth)        
    else:
        print "[%s] Processing long time-binned light curve..." % mysrc

    auto=autoLC(customThreshold=useThresh,daily=daily,longTerm=longTerm,yearmonth=yearmonth,mergelongterm=mergelongterm,withhistory=withhistory,configfile=configfile,stopmonth=stopmonth)
    auto.readSourceList(mysrc)

    if DEBUG:
        print '2FHL counterpart is ', auto.search2FHLcounterpart()

    if longTerm is True and mergelongterm is True:

        # Remove all the old merged file for this source, before reprocessing the merged data
        if not auto.daily:
            for file in glob.glob(auto.workDir+'/'+auto.src+'*'):
                os.remove(file)
        if auto.daily:
            for file in glob.glob(auto.workDir+'/'+auto.src+'*daily*'):
                os.remove(file)


        # TO BE CHANGED !!!
        startyearmonth = '200808'
        # Hardcoded !!! Beurk, not good, ugly, bad !!!

        if auto.stopmonth is not None:
            thisyearmonth = auto.stopmonth
        else:
            thisyearmonth  = datetime.date.today().strftime('%Y%m')
        thisyear       = thisyearmonth[:-2]
        thismonth      = thisyearmonth[-2:]
        startyear      = startyearmonth[:-2]
        startmonth     = startyearmonth[-2:]


        # First make sure that all the month-by-month long-term data have been processed
        #
        # Loop on month from 2008/08 to this month
        for year in range(int(startyear),int(thisyear)+1):
            for month in range(1,12+1):
                # To retrieve the correct results directories, 'month' should be made of 2 digits
                month='%02d'%month
                tmpyearmonth=str(year)+str(month)
                if (year==int(startyear) and int(month) < int(startmonth)) or (year==int(thisyear) and int(month) > int(thismonth)):
                    continue

                # If year=thisyear and month=thismonth, we should remove all data for this source and reprocess everything again with fresh, brand new data !
                # BUT only if update=True
                if year==int(thisyear) and int(month)==int(thismonth) and update is True:
                    tmpworkdir=auto.baseOutDir+"/longTerm/"+str(year)+str(month)
                    if not auto.daily:
                        for file in glob.glob(tmpworkdir+'/'+auto.src+'*'):
                            os.remove(file)
                    if auto.daily:
                        for file in glob.glob(tmpworkdir+'/'+auto.src+'*daily*'):
                            os.remove(file)

                processSrc(mysrc=auto.src,useThresh=useThresh,daily=auto.daily,mail=False,longTerm=True,test=False,yearmonth=tmpyearmonth,mergelongterm=False,update=update,configfile=configfile,stopmonth=stopmonth)

                

        # Then merge the GTI files together, and run createXML, photoLC, exposure, createDAT, createLCfig, createEnergyTimeFig. No mail is sent here.
        auto.mergeGTIfiles()
        if auto.fglName is not None:
            auto.createXML()
            mygamma=None
        else:
            mygamma=ASSUMEDGAMMA
            print '[%s] \033[93mNo 3FGL counterpart given in the list of sources, assuming photon index of %.2f for the light curve generation.\033[0m' % (auto.src, mygamma)
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
        mygamma=None
        FLAGASSUMEDGAMMA=False
    else:
        mygamma=ASSUMEDGAMMA
        print '[%s] \033[93mNo 3FGL counterpart given in the list of sources, assuming photon index of %.2f for the light curve generation.\033[0m' % (auto.src, mygamma)
        FLAGASSUMEDGAMMA=True
    auto.photoLC()
    auto.exposure(gamma=mygamma)
    auto.createDAT()
    auto.createLCfig()
    auto.createEnergyTimeFig()
    alertSent=auto.sendAlert(nomailall=test,sendmail=mail)
    if (alertSent and auto.active and auto.launchLikeAna == 'True') or FORCE_LIKELIHOOD:
        auto.launchLikelihoodAnalysis(nomailall=test)
    
    return auto.active, auto.visible


def main(argv=None):
    """
    Main procedure
    """

    # options parser:
    helpmsg="""%prog [options] <source> [<optional YYYYMM>]

This is the version $Id$

If you call %prog using the -l option, you need to provide a year and a month in input, in the format YYYYMM.

Use '-h' to get the help message

"""

    parser = OptionParser(version="$Id$",
                          usage=helpmsg)

    parser.add_option("-d", "--daily", action="store_true", dest="d", default=False,
                      help='use daily bins for the light curves (defaulted to long time-binned)')
    parser.add_option("--force-daily", action="store_true", dest="force_daily", default=False,
                      help='force daily bins for the light curves')
    parser.add_option("-c", "--custom-threshold", action="store_true", dest="c", default=False,
                      help='use custom trigger thresholds from the master list of sources (defaulted to 1.e-6 ph cm^-2 s^-1)')
    parser.add_option("-l", "--long-term", action="store_true", dest="l", default=False,
                      help='generate a long term light curve, for one given month (defaulted to False). With this option, one should provide a source name as usual, but also a month for which the data should be processed, in the format YYYYMM.')
    parser.add_option("-w", "--with-history", action="store_true", dest="history", default=False,
                      help='use the long-term history of a source to dynamically determine a flux trigger threshold, instead of using a fixed flux trigger threshold as done by default. This option makes use of the long-term data on a source, and assumes that these have been previously generated with the --merge-long-term option.')
    parser.add_option("-m", "--merge-long-term", action="store_true", dest="m", default=False,
                      help='merge the month-by-month long-term light curves together. If those do not exist, they will be created on the fly.')
    parser.add_option("-u","--update",action="store_true",dest="u",default=False,
                      help='update with new data for last month/year when used in conjunction of --merge-long-term. Otherwise, has no effect.')
    parser.add_option("--stop-month", default=None, dest="STOPMONTH", metavar="<STOPMONTH>",
                      help="in conjunction with --merge-long-term, defines the stop year/month (in the format YYYYMM) until which the long-term light curve is generated. '%default' by default.")
    parser.add_option("-n", "--no-mail", action="store_true", dest="n", default=False,
                      help='do not send mail alerts')
    parser.add_option("-t", "--test", action="store_true", dest="t", default=False,
                      help='for test purposes. Do not send the alert mail to everybody if a source is above the trigger threshold, but only to test recipients (by default, mail alerts are sent to everybody, cf. the configuration files).')
    parser.add_option("-f", "--config-file", default='default.cfg', dest="CONFIGFILE", metavar="CONFIGFILE",
                      help="provide a configuration file. Using '%default' by default.")
    parser.add_option("-v", "--verbose", action="store_true", dest="v", default=False,
                      help='verbose output.')
    parser.add_option("--debug", action="store_true", dest="debug", default=False,
                      help='debugging output.')
    parser.add_option("--force-alert", action="store_true", dest="force_alert", default=False,
                      help='force an alert to be issued. To be used along with --force-daily. For test purposes only.')
    parser.add_option("--force-like", action="store_true", dest="force_like", default=False,
                      help='force a likelihood follow-up analysis to be launched. To be used along with --force-daily. For test purposes only.')
    (opt, args) = parser.parse_args()

    CONFIGFILE=opt.CONFIGFILE

    global VERBOSE
    if opt.v:
        VERBOSE=True
    else:
        VERBOSE=False

    global DEBUG
    if opt.debug:
        DEBUG=True
    else:
        DEBUG=False

    # If forcing an alert to be issued
    global FORCE_ALERT
    if opt.force_alert:
        FORCE_ALERT=True
    else:
        FORCE_ALERT=False

    # If forcing a likelihood analysis to be launched
    global FORCE_LIKELIHOOD
    if opt.force_like:
        FORCE_LIKELIHOOD=True
        FORCE_ALERT=True
    else:
        FORCE_LIKELIHOOD=False

    # If daily bins
    if opt.d:
        DAILY=True
    else:
        DAILY=False

    if opt.force_daily:
        FORCE_DAILY=True
        DAILY=True
    else:
        FORCE_DAILY=False

    # If custom thresholds
    if opt.c:
        USECUSTOMTHRESHOLD=True
    else:
        USECUSTOMTHRESHOLD=False

    # If no mail is sent
    if opt.n:
        MAIL=False
    else:
        MAIL=True

    # If test mode
    if opt.t:
        TEST=True
    else:
        TEST=False

    if TEST and not MAIL:
        print "ERROR You asked for both the --test and --no-mail options."
        print "      These options are mutually exclusive."
        sys.exit(1)
    
    # If long term
    if opt.l:
        LONGTERM=True
        # Check that we provided the mandatory argument: a source to process, and a month for which the long-term data should be produced !
        if len(args) != 2:
            print "ERROR Main: wrong number of arguments"
            print "            With the --long-term option, you should provide:"
            print "              - a source name"
            print "              - a yearmonth for which the long term data should be produced"
            sys.exit(1)
        yearmonth=args[1]

    else:
        LONGTERM=False
        # Check that we provided the mandatory argument: a source to process !
        if len(args) != 1:
            print "ERROR Main: wrong number of arguments"
            sys.exit(1)
        yearmonth=None

    # If merge long term light curves
    if opt.m:
        MERGELONGTERM=True
        LONGTERM=True
        #DAILY=False
        MAIL=False
        TEST=False
        # If update long term light curves with brand new data
        if opt.u:
            UPDATE=True
        else:
            UPDATE=False
        if opt.STOPMONTH is not None:
            STOPMONTH = str(opt.STOPMONTH)
        else:
            STOPMONTH = None
    else:
        MERGELONGTERM=False
        UPDATE=False
        STOPMONTH = None

    # If dynamical flux trigger threshold based on source history
    if opt.history:
        WITHHISTORY=True
    else:
        WITHHISTORY=False

    src=args[0]

    processSrc(mysrc=src,useThresh=USECUSTOMTHRESHOLD,daily=DAILY,mail=MAIL,longTerm=LONGTERM,test=TEST,yearmonth=yearmonth,mergelongterm=MERGELONGTERM,withhistory=WITHHISTORY,update=UPDATE,configfile=CONFIGFILE,force_daily=FORCE_DAILY, stopmonth=STOPMONTH)

    return True


if __name__ == '__main__':
    """
    Execute main()
    """

    main()
