#!/usr/bin/python

"""
FLaapLUC (Fermi/LAT automatic aperture photometry Light CUrve)

Automatic generation of aperture photometric light curves of Fermi sources, for a given source.

No likelihood fit is performed, the results solely rely on the 2FGL spectral fits.

More information are available at: http://fermi.gsfc.nasa.gov/ssc/data/analysis/scitools/aperture_photometry.html

@author Jean-Philippe Lenain <mailto:jlenain@in2p3.fr>
@date $Date$
@version $Id$
"""

import sys, os, asciidata, datetime, time, glob
from numpy import *
import pyfits
from optparse import OptionParser
from ConfigParser import ConfigParser

# Import some matplotlib modules
try:
    import matplotlib
    matplotlib.use ('Agg')

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
DEBUG=False # Debugging flag
BATCH=True  # True in batch mode

# Global variables
TOFFSET=54600. # MJD

def met2mjd(met):
    """
    Converts Mission Elapsed Time (MET, in seconds) in Modified Julian Day.
    Cf. http://fermi.gsfc.nasa.gov/ssc/data/analysis/documentation/Cicerone/Cicerone_Data/Time_in_ScienceTools.html
    to see how the time is handled in the Fermi Science Tools.
    
    Input: time in MET (s)
    Output: time in MJD (fraction of a day)
    """
    MJDREFI=51910.0
    MJDREFF=7.428703703703703e-4
    return(MJDREFI+MJDREFF+met/24./60./60.)


def mjd2met(mjd):
    """
    Converts Modified Julian Day in Mission Elapsed Time (MET, in seconds).
    Cf. http://fermi.gsfc.nasa.gov/ssc/data/analysis/documentation/Cicerone/Cicerone_Data/Time_in_ScienceTools.html
    to see how the time is handled in the Fermi Science Tools.
    
    Input:  time in MJD (fraction of a day)
    Output: time in MET (s)
    """
    MJDREFI=51910.0
    MJDREFF=7.428703703703703e-4
    return(24.*60.*60* (mjd - MJDREFI - MJDREFF) )


def unixtime2mjd(unixtime):
    """
    Converts a UNIX time stamp in Modified Julian Day

    Highly inspired from the function 'mjd_now()' from Marcus Hauser's ADRAS/ATOM pipeline

    Input:  time in UNIX seconds
    Output: time in MJD (fraction of a day)

    @todo Should make sure we use UTC here
    @todo We do not yet account for leap seconds (2nd order compared to UTC handling !!)
    """

    # unixtime gives seconds passed since "The Epoch": 1.1.1970 00:00
    # MJD at that time was 40587.0

    result = 40587.0 + unixtime / (24.*60.*60.)
    return result


def zaAtCulmination(dec):
    """
    Returns the zenith angle of a source at culmination, for the H.E.S.S. site.
    """
    siteLat = -23.27167 # latitude of H.E.S.S. site in degree
    return abs(dec-siteLat)



def getConfigList(option,sep=','):
    return [ stuff for stuff in option.split(sep) ]


class autoLC:
    """
    FLaapLUC

    Automatic aperture photometry light curve generation.
    Main class, for a given of source.
    """

    def __init__(self,file='listSources.txt',customThreshold=False,daily=False,longTerm=False,yearmonth=None,mergelongterm=False,withhistory=False,configfile='default.cfg'):
        
        self.config           = self.getConfig(configfile=configfile)
        #self.file=file
        self.allskyDir        = self.config.get('InputDirs','AllskyDir')
        self.archiveDir       = self.config.get('InputDirs','ArchiveDir')
        self.templatesDir     = self.config.get('InputDirs','TemplatesDir')
        self.ATOMSchedulesDir = self.config.get('InputDirs','ATOMSchedulesDir')
        self.catalogFile      = self.config.get('InputFiles','CatalogFile')
        self.file             = self.config.get('InputFiles','SourceList')
        self.baseOutDir       = self.config.get('OutputDirs','OutputResultsDir')
        self.allskyFile       = self.allskyDir+"/"+self.config.get('InputFiles','WholeAllskyFile')
        self.lastAllskyFile   = self.allskyDir+"/"+self.config.get('InputFiles','LastAllskyFile')
        self.spacecraftFile   = self.allskyDir+"/"+self.config.get('InputFiles','SpacecraftFile')
        self.webpageDir       = self.config.get('OutputDirs','OutputWebpageDir')
        self.url              = self.config.get('OutputDirs','URL')
        # Read maxz and maxZA as lists, not as single floats
        self.maxz             = [float(i) for i in getConfigList(self.config.get('AlertTrigger','MaxZ'  ))]
        self.maxZA            = [float(i) for i in getConfigList(self.config.get('AlertTrigger','MaxZA'))]

        self.daily            = daily
        self.withhistory      = withhistory

        today=datetime.date.today().strftime('%Y%m%d')

        # Setting file names and directories
        if longTerm is True:
            self.allsky = self.allskyFile
            if mergelongterm is False:
                self.workDir    = self.baseOutDir+"/longTerm/"+yearmonth
            else:
                self.workDir    = self.baseOutDir+"/longTerm/merged"
        else:
            # If lock files exist in the archive directory (e.g. if NASA servers are down), we do not do anything and exit
            photonLock=self.archiveDir+"/photon.lock"
            spacecraftLock=self.archiveDir+"/spacecraft.lock"
            if os.path.isfile(photonLock) or os.path.isfile(spacecraftLock):
                self.sendErrorMail(mailall=True)
                sys.exit(10)

            self.allsky     = self.lastAllskyFile
            self.workDir    = self.baseOutDir+"/"+today

        self.spacecraft = self.spacecraftFile
        if not os.path.isdir(self.workDir):
            os.makedirs(self.workDir)

        self.fermiDir   = os.getenv('FERMI_DIR')

        # Setting default parameters
        self.roi       = 1.   # degrees (http://fermi.gsfc.nasa.gov/ssc/data/analysis/scitools/aperture_photometry.html: "For aperture photometry we select a very small aperture (rad=1 degree), because we are not fitting the background.")
        self.emin      = 1.e2 # E min
        self.emax      = 3.e5 # E max
        self.zmax      = 100. # degrees
        self.rockangle = 52. # maximal allowed rocking angle

        if daily:
            self.tbin =    24.*60.*60. # seconds, daily bins
        else:
            self.tbin = 7.*24.*60.*60. # seconds, weekly bins

        self.threshold = 1.e-6 # ph cm^-2 s^-1
        self.customThreshold=customThreshold
        
        # Open allsky file to get the start and stop dates
        try:
            hdu=pyfits.open(self.allsky)
        except:
            print 'Exception: can not open file '+self.allsky
            raise
        header = hdu[0].header

        if longTerm is False:
            self.tstart = header['TSTART']
            self.tstop  = header['TSTOP']
            
        else:
            missionStart = header['TSTART'] # in MET
            missionStop  = header['TSTOP']  # in MET

            if mergelongterm is False:
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
                    print 'INIT yearmonthStart=',yearmonthStart
                    print 'INIT yearmonthStop=',yearmonthStop

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
                self.tstop  = missionStop

        # Mail sender and recipients
        self.usualRecipients= getConfigList(self.config.get('MailConfig','UsualRecipients'))
        self.testRecipients = getConfigList(self.config.get('MailConfig','TestRecipients'))
        self.mailSender     = self.config.get('MailConfig','MailSender')
    


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
        if self.customThreshold is True:
            myThreshold=srcList[5].tonumpy()

    
        # If we ask for a particular source, return the parameters for that source
        if mysrc != None:
            # Find our input src in the list of sources
            found=False
            for i in range(len(src)):
                if src[i]==mysrc:
                    found=True

                    # Redefine the threshold if we provided a custom threshold
                    if self.customThreshold is True and myThreshold[i] != 0.:
                        try:
                            float(myThreshold[i])
                            self.threshold=myThreshold[i]
                        except ValueError:
                            print "WARNING The threshold of the source "+str(mysrc)+" is not a float. Please, check the list of sources !"
                            sys.exit(2)
                    return src[i],ra[i],dec[i],z[i],fglName[i]
            
            # If we end up without any found source, print out a WARNING
            print "WARNING Can't find your source "+str(mysrc)+" in the list of sources !"
            return None,None,None,None,None
        
        # Otherwise, return the whole list of parameters for all the sources
        else:
            return src,ra,dec,z,fglName


    def selectSrc(self,src,ra,dec):
        """
        Filter a given source, running gtselect
        """
        filter['infile']=self.allsky
        if self.daily:
            outfile=self.workDir+'/'+str(src)+'_daily.fits'
        else:
            outfile=self.workDir+'/'+str(src)+'.fits'
        filter['outfile']=outfile
        
        # If outfile already exists, we don't do anything
        if os.path.isfile(outfile):
            return True

        filter['ra']=ra
        filter['dec']=dec
        filter['rad']=self.roi
        filter['emin']=self.emin
        filter['emax']=self.emax
        filter['tmin']=self.tstart
        filter['tmax']=self.tstop
        filter['zmax']=self.zmax
        filter['evclass']=2
        filter.run()


    def makeTime(self,src,ra,dec):
        """
        Filter the GTI for a given source
        """
        maketime['scfile']=self.spacecraft

        if self.daily:
            maketime['evfile']=self.workDir+'/'+str(src)+'_daily.fits'
            outfile=self.workDir+'/'+str(src)+'_daily_gti.fits'
        else:
            maketime['evfile']=self.workDir+'/'+str(src)+'.fits'
            outfile=self.workDir+'/'+str(src)+'_gti.fits'
        maketime['outfile']=outfile
        
        # If outfile already exists, we don't do anything
        if os.path.isfile(outfile):
            return True

        # cf. http://fermi.gsfc.nasa.gov/ssc/data/analysis/scitools/aperture_photometry.html
        #maketime['filter']="IN_SAA!=T && LAT_CONFIG==1 && DATA_QUAL==1 && (angsep(RA_ZENITH,DEC_ZENITH,"+str(ra)+","+str(dec)+")+"+str(self.roi)+" <"+str(self.zmax)+") && (angsep("+str(ra)+","+str(dec)+",RA_SCZ,DEC_SCZ)<180.) && (angsep("+str(ra)+","+str(dec)+",RA_SUN,DEC_SUN)>5.)"
        maketime['filter']="IN_SAA!=T && LAT_CONFIG==1 && DATA_QUAL==1 && ABS(ROCK_ANGLE)<"+str(self.rockangle)+" && (angsep(RA_ZENITH,DEC_ZENITH,"+str(ra)+","+str(dec)+")+"+str(self.roi)+" <"+str(self.zmax)+") && (angsep("+str(ra)+","+str(dec)+",RA_SUN,DEC_SUN)>5.)"
        maketime['roicut']='no'
        maketime['tstart']=self.tstart
        maketime['tstop']=self.tstop
        maketime.run()


    def mergeGTIfiles(self,src,ra,dec,daily=False):
        """
        Merge multiple GTI files when mergelongterm is True.
        Use gtselect.
        Assume the current workDir is longTerm/merged.
        """

        # Create list of GTI files
        if not daily:
            listname=self.workDir+'/'+src+'_gti.list'
        else:
            listname=self.workDir+'/'+src+'_daily_gti.list'
        filelist=open(listname,'w')
        list=[]
        if not daily:
            for file in glob.glob(self.workDir+'/../20????/'+src+'_gti.fits'):
                list.append(file)
        else:
            for file in glob.glob(self.workDir+'/../20????/'+src+'_daily_gti.fits'):
                list.append(file)
        # Sort the list of GTI files
        list=sorted(list)
        for item in list:
            filelist.write(item+'\n')
        filelist.close()
        
        filter['infile']='@'+listname
        if not daily:
            outfile=self.workDir+'/'+str(src)+'_gti.fits'
        else:
            outfile=self.workDir+'/'+str(src)+'_daily_gti.fits'
        filter['outfile']=outfile
        
        # If outfile already exists, we re-create it
        if os.path.isfile(outfile):
            os.remove(outfile)

        filter['ra']      = ra
        filter['dec']     = dec
        filter['rad']     = self.roi
        filter['emin']    = self.emin
        filter['emax']    = self.emax
        filter['tmin']    = self.tstart
        filter['tmax']    = self.tstop
        filter['zmax']    = self.zmax
        filter['evclass'] = 2
        filter.run()


    def createXML(self,src):
        """
        Create an XML model file based on the 2FGL catalogue
        """
        
        if self.daily:
            evfile=self.workDir+'/'+str(src)+'_daily_gti.fits'
            modelfile=self.workDir+'/'+str(src)+'_daily.xml'
        else:
            evfile=self.workDir+'/'+str(src)+'_gti.fits'
            modelfile=self.workDir+'/'+str(src)+'.xml'

        # If modelfile already exists, we don't do anything
        if os.path.isfile(modelfile):
            return True

        try:
            import make2FGLxml
        except ImportError:
            print "ERROR Can't import make2FGLxml."
            sys.exit(1)
        
        mymodel=make2FGLxml.srcList(self.catalogFile,evfile,modelfile)
        mymodel.makeModel(self.fermiDir+'/refdata/fermi/galdiffuse/gal_2yearp7v6_v0.fits','Gal_2yearp7v6_v0',self.fermiDir+'/refdata/fermi/galdiffuse/iso_p7v6source.txt','iso_p7v6source',self.templatesDir)


    def photoLC(self,src):
        """
        Compute the photometric light curve for a given source
        """

        if self.daily:
            evtbin['evfile']=self.workDir+'/'+str(src)+'_daily_gti.fits'
            outfile=self.workDir+'/'+str(src)+'_daily_lc.fits'
        else:
            evtbin['evfile']=self.workDir+'/'+str(src)+'_gti.fits'
            outfile=self.workDir+'/'+str(src)+'_lc.fits'

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
        evtbin.run()


    def exposure(self,src,fglName,gamma=None):
        """
        Compute exposure on source src, to add a flux column for the photometric light curve.

        Warning: the input file is modified in place, with an additional exposure column added to the file !
        """

        if self.daily:
            infile=self.workDir+'/'+str(src)+'_daily_lc.fits'
            srcmdl=self.workDir+'/'+str(src)+'_daily.xml'
        else:
            infile=self.workDir+'/'+str(src)+'_lc.fits'
            srcmdl=self.workDir+'/'+str(src)+'.xml'

        # If infile already contains an EXPOSURE column, we don't do anything
        hdu=pyfits.open(infile)
        if hdu[1].header.get('TTYPE5')=='EXPOSURE':
            return True
 

        scfile=self.spacecraft
        irfs='P7SOURCE_V6'
        target=fglName
        rad=str(self.roi)
        
        if gamma is None:
            options='infile='+infile+' scfile='+scfile+' irfs='+irfs+' srcmdl='+srcmdl+' target='+target+' rad='+rad
        else:
            options='infile='+infile+' scfile='+scfile+' irfs='+irfs+' srcmdl="none" specin='+str(gamma)+' rad='+rad
        cmd='time -p '+self.fermiDir+'/bin/gtexposure '+options
        os.system(cmd)




    def createDAT(self,src):
        """
        Create a data file with the light curve of a given source.
        """

        # Read LC file
        if self.daily:
            infile=self.workDir+'/'+str(src)+'_daily_lc.fits'
            outfile=self.workDir+'/'+str(src)+'_daily_lc.dat'
        else:
            infile=self.workDir+'/'+str(src)+'_lc.fits'
            outfile=self.workDir+'/'+str(src)+'_lc.dat'

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
            # Exposure can be 0 if longTerm=True and TSTOP in photon file > TSTOP in spacecraft file.
            if exposure[i] != 0.:
                file.write(str(time[i])+"\t"+str(timeMjd[i])+"\t"+str(flux[i])+"\t"+str(fluxErr[i])+"\n")
        file.close()

    def getBAT(self,src):
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
        if urls.has_key(src):
            file=urls[src].replace('+','p')+".lc.fits"
        else:
            file=src.replace('+','p')+".lc.fits"
        urlprefix="http://swift.gsfc.nasa.gov/docs/swift/results/transients/"

        # lc files can be in a weak/ subdir for weak sources, we try both
        try:
            baturl=urlprefix+file
            webfile=urllib2.urlopen(baturl)
        except urllib2.HTTPError:
            try:
                baturl=urlprefix+'weak/'+file
                webfile=urllib2.urlopen(baturl)
            except urllib2.HTTPError:
                return False,None

        # save lc to local file
        localfile=open(file,'w')
        localfile.write(webfile.read())
        webfile.close()
        localfile.close()
        # read local file with pyfits into batlc
        batfits=pyfits.open(file)
        batlc=np.array(batfits[1].data)
        batfits.close()
        # delete local file
        os.unlink(file)

        return True,batlc


        
    def createPNG(self,src,fglName,z):
        """
        Create a PNG figure with the light curve of a given source. Any existing PNG file is overwritten !
        """

        # Read the .dat LC file
        if self.daily:
            infile=self.workDir+'/'+str(src)+'_daily_lc.dat'
            outfig=self.workDir+'/'+str(src)+'_daily_lc.png'
            infileWeekly=self.workDir+'/'+str(src)+'_lc.dat'
            duration = 1. # duration of a time bin, in days
        else:
            infile=self.workDir+'/'+str(src)+'_lc.dat'
            outfig=self.workDir+'/'+str(src)+'_lc.png'
            duration = 7. # duration of a time bin, in days

        data    = asciidata.open(infile)
        # the times are already read as MJD, cf createDAT function.
        timelc  = data[1].tonumpy()
        flux    = data[2].tonumpy()
        fluxErr = data[3].tonumpy()

        if self.daily:
            dataWeekly    = asciidata.open(infileWeekly)
            # the times are already read as MJD, cf createDAT function.
            timeWeekly    = dataWeekly[1].tonumpy()
            fluxWeekly    = dataWeekly[2].tonumpy()
            fluxErrWeekly = dataWeekly[3].tonumpy()
            durationWeekly= 7. # duration of a time bin, in days

        # Download Swift/BAT data if available
        # xray is boolean flag indicating that X-ray BAT data is available
        xray,batlc=self.getBAT(src)

        # Redefine the trigger threshold if withhistory=True
        if self.withhistory is True:
            (fluxAverage,fluxRMS) = self.dynamicalTrigger(src)


        fig=figure()

        if xray:
            ax    = fig.add_subplot(211)
            axbat = fig.add_subplot(212,sharex=ax)
        else:
            ax = fig.add_subplot(111)

        if fglName is not None:
            title=str(src)+', '+str(fglName).replace('_2FGLJ','2FGL J')
        else:
            title=str(src)+', no known 2FGL counterpart'
        if str(z)=='--': # this is the result of the conversion of None to a float
            title=title+' (z unknown)'
        else:
            title=title+' (z='+str(z)+')'

        ax.set_title(title)

        
        # Force the y-axis ticks to use 1e-6 as a base exponent
        ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: ('%.1f')%(x*1e6)))
        ax.set_ylabel('F (%.0f MeV-%.0f GeV) (x 10^-6 ph cm^-2 s^-1)'%(self.emin,self.emax/1000.))

        day=24.*60.*60.

        # Make the x-axis ticks shifted by some value
        ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: '%.0f'%(x-TOFFSET)))
        ax.set_xlabel('MJD-'+str(TOFFSET))
        
        # Plot the Fermi/LAT light curve
        if self.daily:
            # Also plot the weekly-binned light curve
            ax.errorbar(x=timelc, xerr=duration/2., y=flux, yerr=fluxErr/2., fmt='ro')
            ax.errorbar(x=timeWeekly, xerr=durationWeekly/2., y=fluxWeekly, yerr=fluxErrWeekly/2., fmt='bo')
            # The last plot called is on top of the others in matplotlib (are you sure ???). Here, we want the weekly-binned LC on top, for visibility.
        else:
            ax.errorbar(x=timelc, xerr=duration/2., y=flux, yerr=fluxErr/2., fmt='bo')

        # TODO if threshold is dynamic and it is computed using the last flux
        # measurement uncertainty, threshold will not be a given horizontal
        # line, but depend on each point. In that case, indicate in different
        # color those flux measurements that trigger

        # Plot a line at the threshold value
        ax.axhline(y=self.threshold,linewidth=3,linestyle='--',color='r')
        if self.withhistory is True:
            ax.axhline(y=fluxAverage,linewidth=1,linestyle='-',color='b')
            ax.axhline(y=fluxAverage+fluxRMS,linewidth=1,linestyle='--',color='b')
            ax.axhline(y=fluxAverage-fluxRMS,linewidth=1,linestyle='--',color='b')

        # Plot a line at flux=0, for visibility/readibility
        ax.axhline(y=0.,color='k')

        # Plot Swift/BAT lightcurve
        if xray:
            #axbat.errorbar(batlc['TIME']+0.5-TOFFSET,batlc['RATE'],batlc['ERROR'],fmt=None,capsize=0,elinewidth=1,ecolor='b',color='b')
            axbat.errorbar(batlc['TIME'],batlc['RATE'],batlc['ERROR'],fmt=None,capsize=0,elinewidth=1,ecolor='b',color='b')
            axbat.set_xlabel('MJD-'+str(TOFFSET))
            axbat.set_ylabel('F (15-150 keV) (count cm^-2 s^-1)')
            axbat.set_xlim(xmin=timelc[0]-1,xmax=timelc[-1:]+1)
            axbat.set_ylim(ymin=0.)

        # Add a label for the creation date of this figure (highly inspired from Marcus Hauser's ADRAS/ATOM pipeline)
        # x,y in relative 0-1 coords in figure
        figtext(0.98, 0.05,
                'plot creation date: %s (UTC)'%(time.strftime("%a, %d %b %Y %H:%M:%S", time.gmtime())),
                horizontalalignment="left",
                rotation='vertical',
                size='xx-small'
                )
        
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


    def killTrigger(self,ra,dec,z):
        """
        Defines cuts on (Ra,Dec,z) before assessing whether a mail alert shoudl be sent for a source which flux is above the trigger threshold.
        We cut on a combination (z, ZenithAngle), using a bit mask.

        @todo: Introduce an additional cut on Gal latitude ?
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


        # Equatorial coordinates -> Galactic coordinates transformation
        equat2gal = wcs.Transformation((wcs.equatorial, wcs.fk5,'J2000.0'),wcs.galactic)
    
        # WCS module needs a numpy matrix for coord transformation
        srcradec=matrix(([ra],[dec]))
        srcGalCoord=equat2gal(srcradec)
        # Retrieve the Galactic latitude
        srcGalLat=srcGalCoord[(1,0)]
        
        zaAtCulmin=zaAtCulmination(dec)
    

        # Mask on both (z, ZA at culmin)
        #          z column               ZA column
        msk = (z<=grid[:,0])&(zaAtCulmin<=grid[:,1])
    
        # the 'return' below is a bit counter-intuitive. It answers the question 'Should we kill an imminent mail alert ?', i.e. if a source has the last flux point above the flux threshold, does it also fulfill the requirements on both z (not too far away) and zenith angle (not too low in the sky) ? So if an alert should definitely be sent, this function returns 'False' !

        # if the mask has at least one 'True' element, we should send an alert
        if True in msk:
            # print 'An alert should be triggered !'
            return False
        else:
            # print 'No alert triggered'
            return True


    def dynamicalTrigger(self,src):
        '''
        If long-term data are available for a source, dynamically computes a flux trigger threshold based on the flux history of the source. Otherwise, fall back with default fixed trigger threshold.

        @param src Soure name
        '''
        
        # Read the longterm .dat LC file
        infile = self.baseOutDir+'/longTerm/merged/'+str(src)+'_lc.dat'
        try:
            data    = asciidata.open(infile)
        except IOError:
            print 'Long term data file unavailable for source %s' % src
            # Falling back to default fixed trigger threshold
            self.withhistory=False
            return (False,False)

        flux    = data[2].tonumpy()
        fluxErr = data[3].tonumpy()

        # weighted average of the historical fluxes, weighted by their errors
        fluxAverage = average(flux, weights=1./fluxErr)
        fluxRMS     = std(flux, dtype=np.float64)

        self.threshold = fluxAverage + 3.0*fluxRMS
        return (fluxAverage,fluxRMS)
        


    def sendAlert(self,src,ra,dec,z,nomailall=False):
        '''
        Send a mail alert in case a source fulfills the trigger conditions.

        @param src Source name
        @param ra Source right ascension
        @param dec Source declination
        @param z Source redshift
        @param nomailall Boolean, should the mail be sent to a restricted list of recipients ?
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



        # Read the light curve file
        if self.daily:
            infile  = self.workDir+'/'+str(src)+'_daily_lc.dat'
            pngFig=self.workDir+'/'+str(src)+'_daily_lc.png'

            # Also take a look in the weekly data
            infileWeekly=self.workDir+'/'+str(src)+'_lc.dat'
            dataWeekly=asciidata.open(infileWeekly)
            timeWeekly=dataWeekly[0].tonumpy()
            fluxWeekly=dataWeekly[2].tonumpy()
            fluxErrWeekly=dataWeekly[3].tonumpy()
            # Catch the last flux point
            lastTimeWeekly=timeWeekly[-1:]
            lastFluxWeekly=fluxWeekly[-1:]
            lastFluxErrWeekly=fluxErrWeekly[-1:]
        else:
            infile  = self.workDir+'/'+str(src)+'_lc.dat'
            pngFig=self.workDir+'/'+str(src)+'_lc.png'
        data    = asciidata.open(infile)
        time    = data[0].tonumpy()
        flux    = data[2].tonumpy()
        fluxErr = data[3].tonumpy()

        # Catch the last flux point
        lastTime    = time[-1:]
        lastFlux    = flux[-1:]
        lastFluxErr = fluxErr[-1:]

        if DEBUG:
            print 'DEBUG: ',src,self.threshold
            print
            print "self.threshold=",self.threshold
            print "lastFlux=",lastFlux
            print

        KILLTRIGGER = self.killTrigger(ra,dec,z)

        if not KILLTRIGGER:
            # Assess whether the trigger condition is met, looking at the last flux point
            if self.daily:
                if (lastFlux >= self.threshold or lastFluxWeekly >= self.threshold):
                    SENDALERT=True
                else:
                    SENDALERT=False
            else:
                if lastFlux >= self.threshold:
                    SENDALERT=True
                else:
                    SENDALERT=False
        else:
            SENDALERT=False

        if DEBUG:
            print str(src),dec,z,self.maxZA,self.maxz,KILLTRIGGER,SENDALERT

        # If trigger condition is met, we send a mail
        if SENDALERT:
            # Create the container email message.
            msg = MIMEMultipart()
            msg['Subject'] = '[FLaapLUC] Fermi/LAT flare alert on %s' % src
            sender = self.mailSender
            
            # To whom the mail should be sent (cf. __init__ function of the class)
            if nomailall is False:
                recipient = self.usualRecipients
            else:
                recipient = self.testRecipients

            msg['From'] = sender
            COMMASPACE = ', '
            msg['To'] =COMMASPACE.join( recipient )
            msg.preamble = 'You will not see this in a MIME-aware mail reader.\n'
            # Guarantees the message ends in a newline
            msg.epilogue = ''
            
            mailtext="""
     FLaapLUC (Fermi/LAT automatic aperture photometry Light CUrve) report

     *** The Fermi/LAT flux (%.0f MeV-%.0f GeV) of %s exceeds the trigger threshold of %.2g ph cm^-2 s^-1 ***

     """%(self.emin,self.emax/1000.,src,self.threshold)

            if self.daily:
                mailtext=mailtext+"""

     The last daily-binned flux is:      %.2g +/- %.2g ph cm^-2 s^-1, centred on MET %.0f
     and the last weekly-binned flux is: %.2g +/- %.2g ph cm^-2 s^-1, centred on MET %.0f

"""%(lastFlux,lastFluxErr,lastTime,lastFluxWeekly,lastFluxErrWeekly,lastTimeWeekly)
                mailtext=mailtext+"The most recent lightcurve (%.0f-day binned in red, and weekly binned in blue) is attached."%(self.tbin/24./60./60.)
            else:
                mailtext=mailtext+"""

     The last weekly-binned flux is:      %.2g +/- %.2g ph cm^-2 s^-1, centred on MET %.0f

"""%(lastFlux,lastFluxErr,lastTime)
                mailtext=mailtext+"The most recent lightcurve (%.0f-day binned) is attached."%(self.tbin/24./60./60.)
            
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
            
            # Open the files in binary mode.  Let the MIMEImage class automatically guess the specific image type.
            fp = open(pngFig, 'rb')
            img = MIMEImage(fp.read())
            fp.close()
            msg.attach(img)

            # Send the email via our own SMTP server.
            s = smtplib.SMTP()
            s.set_debuglevel(0)
            s.connect()
            s.sendmail(sender, recipient, msg.as_string())
            s.quit()

            print "Alert sent for %s"%src

        return True



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
            
        if mailall is True:
            recipient = self.usualRecipients
        else:
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




def processSrc(mysrc=None,useThresh=False,daily=False,mail=True,longTerm=False,test=False, yearmonth=None, mergelongterm=False,withhistory=False,configfile='default.cfg'):
    """
    Process a given source.
    """

    if DEBUG:
        print 'src=',mysrc

    
    if mysrc is None:
        print "ERROR Missing input source !"
        sys.exit(1)

    auto=autoLC(customThreshold=useThresh,daily=daily,longTerm=longTerm,yearmonth=yearmonth,mergelongterm=mergelongterm,withhistory=withhistory,configfile=configfile)
    src,ra,dec,z,fglName=auto.readSourceList(mysrc)


    if longTerm is True and mergelongterm is True:

        # Remove all the old merged file for this source, before reprocessing the merged data
        if not daily:
            for file in glob.glob(auto.workDir+'/'+src+'*'):
                os.remove(file)
        if daily:
            for file in glob.glob(auto.workDir+'/'+src+'*daily*'):
                os.remove(file)


        # TO BE CHANGED !!!
        startyearmonth = '200808'
        # Hardcoded !!! Beurk, not good, ugly, bad !!!

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
                if year==int(thisyear) and int(month)==int(thismonth):
                    tmpworkdir=auto.baseOutDir+"/longTerm/"+str(year)+str(month)
                    if not daily:
                        for file in glob.glob(tmpworkdir+'/'+src+'*'):
                            os.remove(file)
                    if daily:
                        for file in glob.glob(tmpworkdir+'/'+src+'*daily*'):
                            os.remove(file)

                processSrc(mysrc=src,useThresh=useThresh,daily=daily,mail=False,longTerm=True,test=False,yearmonth=tmpyearmonth,mergelongterm=False,configfile=configfile)

                

        # Then merge the GTI files together, and run createXML, photoLC, exposure, createDAT and createPNG. No mail is sent here.
        auto.mergeGTIfiles(src,ra,dec,daily=daily)
        if fglName is not None:
            auto.createXML(src)
            mygamma=None
        else:
            mygamma=-2.5
            print 'Your source '+src+' has no 2FGL counterpart given in the list of sources. I will assume a photon index of '+str(mygamma)+' for the light curve generation.'
        auto.photoLC(src)
        auto.exposure(src,fglName,gamma=mygamma)
        auto.createDAT(src)
        auto.createPNG(src,fglName,z)
        # Exit here
        return True



    # When mergelongterm is False, we do the following:
    auto.selectSrc(src,ra,dec)
    auto.makeTime(src,ra,dec)
    # If we are in --long-term mode, but not in --merge-long-term mode, we can stop here, since the --merge-long-term mode then starts at the mergeGTIfiles level
    if longTerm is True:
        return True

    if fglName is not None:
        auto.createXML(src)
        mygamma=None
    else:
        mygamma=-2.5
        print 'Your source '+src+' has no 2FGL counterpart given in the list of sources. I will assume a photon index of '+str(mygamma)+' for the light curve generation.'
    auto.photoLC(src)
    auto.exposure(src,fglName,gamma=mygamma)
    auto.createDAT(src)
    auto.createPNG(src,fglName,z)
    if mail is True:
        auto.sendAlert(src,ra,dec,z,nomailall=test)

    
    return True


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
                      help='use daily bins for the light curves (defaulted to weekly)')
    parser.add_option("-c", "--custom-threshold", action="store_true", dest="c", default=False,
                      help='use custom trigger thresholds from the master list of sources (defaulted to 1.e-6 ph cm^-2 s^-1)')
    parser.add_option("-l", "--long-term", action="store_true", dest="l", default=False,
                      help='generate a long term light curve, for one given month (defaulted to False). With this option, one should provide a source name as usual, but also a month for which the data should be processed, in the format YYYYMM.')
    parser.add_option("--with-history", action="store_true", dest="history", default=False,
                      help='use the long-term history of a source to dynamically determine a flux trigger threshold, instead of using a fixed flux trigger threshold done by default.')
    parser.add_option("-m", "--merge-long-term", action="store_true", dest="m", default=False,
                      help='merge the month-by-month long-term light curves together. If those do not exist, they will be created on the fly.')
    parser.add_option("-n", "--no-mail", action="store_true", dest="n", default=False,
                      help='do not send mail alerts')
    parser.add_option("-t", "--test", action="store_true", dest="t", default=False,
                      help='for test purposes. Do not send the alert mail to everybody if a source is above the trigger threshold, but only to J.-P. Lenain (by default, mail alerts are sent to everybody)')
    parser.add_option("--config-file", default='default.cfg', dest="CONFIGFILE", metavar="CONFIGFILE",
                      help="provide a configuration file. Using '%default' by default.")

    (opt, args) = parser.parse_args()

    CONFIGFILE=opt.CONFIGFILE

    # If daily bins
    if opt.d:
        DAILY=True
    else:
        DAILY=False

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

    if TEST is True and MAIL is False:
        print "ERROR You asked for both the --test and --no-mail options."
        print "      These are mutually exclusive options."
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
    else:
        MERGELONGTERM=False
    
    # If dynamical flux trigger threshold based on source history
    if opt.history:
        WITHHISTORY=True
    else:
        WITHHISTORY=False
    
    src=args[0]


    # If we asked for a daily light curve, first make sure that the weekly-binned data already exsits, otherwise this script will crash, since the daily-binned PNG needs the weekly-binned data to be created. No mail alert is sent at this step.
    # We automatically recreate here any missing weekly-binned data.
    if DAILY:
        processSrc(mysrc=src,useThresh=USECUSTOMTHRESHOLD,daily=False,mail=False,longTerm=LONGTERM,yearmonth=yearmonth,mergelongterm=MERGELONGTERM,withhistory=WITHHISTORY,configfile=CONFIGFILE)

    processSrc(mysrc=src,useThresh=USECUSTOMTHRESHOLD,daily=DAILY,mail=MAIL,longTerm=LONGTERM,test=TEST,yearmonth=yearmonth,mergelongterm=MERGELONGTERM,withhistory=WITHHISTORY,configfile=CONFIGFILE)

    return True


if __name__ == '__main__':
    """
    Execute main()
    """

    main()
