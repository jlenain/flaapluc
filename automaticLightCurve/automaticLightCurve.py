#!/usr/bin/env python

"""
Automatic generation of aperture photometric light curves of Fermi sources, for a given source.

No likelihood fit is performed, the results solely rely on the 2FGL spectral fits.

More information are available at: http://fermi.gsfc.nasa.gov/ssc/data/analysis/scitools/aperture_photometry.html

@author Jean-Philippe Lenain <mailto:jplenain@lsw.uni-heidelberg.de>
@date $Date$
@version $Id$
"""

import sys, os, asciidata, datetime
from numpy import *
#from multiprocessing import Process, Queue
from optparse import OptionParser

# Import some matplotlib modules
try:
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


class autoLC:
    """
    Automatic aperture photometry light curve generation, for a list of sources
    """

    def __init__(self,file="/home/fermi/local/automaticLightCurve/listSources.txt",customThreshold=False,daily=False):
        self.file=file

        # Setting file names and directories
        #self.allsky     = "/data/fermi/allsky/allsky_30MeV_300GeV_diffuse_filtered.fits"
        #self.allsky     = "/data/fermi/allsky/allsky_lastMonth_30MeV_300GeV_diffuse_filtered.fits"
        self.allsky     = "/data/fermi/allsky/allsky_last70days_30MeV_300GeV_diffuse_filtered.fits"
        self.spacecraft = "/data/fermi/allsky/allsky_SC00.fits"
        self.workDir    = "/home/fermi/data/automaticLightCurveOutput/"+datetime.date.today().strftime('%Y%m%d')
        if not os.path.isdir(self.workDir):
            os.makedirs(self.workDir)

        self.fermiDir   = os.getenv('FERMI_DIR')

        # Setting default parameters
        self.roi  = 1.   # degrees (http://fermi.gsfc.nasa.gov/ssc/data/analysis/scitools/aperture_photometry.html: "For aperture photometry we select a very small aperture (rad=1 degree), because we are not fitting the background.")
        self.emin = 1.e2 # E min
        self.emax = 1.e5 # E max
        self.zmax = 100. # degrees
        self.daily=daily

        if not daily:
            self.tbin = 7.*24.*60.*60. # seconds, weekly bins
        else:
            self.tbin = 24.*60.*60. # seconds, daily bins

        self.threshold = 1.e-6 # ph cm^-2 s^-1
        self.customThreshold=customThreshold

        # Open allsky file to get the start and stop dates
        import pyfits
        try:
            hdu=pyfits.open(self.allsky)
        except:
            print 'Exception: can not open file '+self.allsky
            raise
        header      = hdu[0].header
        self.tstart = header['TSTART']
        self.tstop  = header['TSTOP']


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

                    # Redefine the class threshold if we provided a custom threshold
                    if self.customThreshold is True and myThreshold[i] != 0.:
                        self.threshold=myThreshold[i]
                        
                    return src[i],ra[i],dec[i],z[i],fglName[i]
            
            # If we end up without any found source, print an error and exits
            print "WARNING Can't find your source "+str(mysrc)+" in the list of sources !"
            return None,None,None,None,None
            #sys.exit(1)
  
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
        
        # If outfile exsits, we remove it before updating it
        if os.path.isfile(outfile):
            os.remove(outfile)

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
        
        # If outfile exsits, we remove it before updating it
        if os.path.isfile(outfile):
            os.remove(outfile)

        # cf. http://fermi.gsfc.nasa.gov/ssc/data/analysis/scitools/aperture_photometry.html
        maketime['filter']="IN_SAA!=T && LAT_CONFIG==1 && DATA_QUAL==1 && (angsep(RA_ZENITH,DEC_ZENITH,"+str(ra)+","+str(dec)+")+"+str(self.roi)+" <"+str(self.zmax)+") && (angsep("+str(ra)+","+str(dec)+",RA_SCZ,DEC_SCZ)<180.)"
        maketime['roicut']='no'
        maketime['tstart']=self.tstart
        maketime['tstop']=self.tstop
        maketime.run()


    def createXML(self,src):
        """
        Create an XML model file based on the 2FGL catalogue
        """
        
        try:
            import make2FGLxml
        except ImportError:
            print "ERROR Can't import make2FGLxml."
            sys.exit(1)
        
        if self.daily:
            evfile=self.workDir+'/'+str(src)+'_daily_gti.fits'
            modelfile=self.workDir+'/'+str(src)+'_daily.xml'
        else:
            evfile=self.workDir+'/'+str(src)+'_gti.fits'
            modelfile=self.workDir+'/'+str(src)+'.xml'

        # If modelfile exsits, we remove it
        if os.path.isfile(modelfile):
            if DEBUG:
                return
            os.remove(modelfile)
        
        mymodel=make2FGLxml.srcList('./gll_psc_v07.fit',evfile,modelfile)
        mymodel.makeModel(self.fermiDir+'/refdata/fermi/galdiffuse/gal_2yearp7v6_v0.fits','Gal_2yearp7v6_v0',self.fermiDir+'/refdata/fermi/galdiffuse/iso_p7v6source.txt','iso_p7v6source','/home/jplenain/fermi/2FGL/Templates')


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

        # If outfile exsits, we remove it before updating it
        if os.path.isfile(outfile):
            if DEBUG:
                return
            os.remove(outfile)

        evtbin['outfile']=outfile
        evtbin['scfile']=self.spacecraft
        evtbin['algorithm']='LC'
        evtbin['tbinalg']='LIN'
        evtbin['tstart']=self.tstart
        evtbin['tstop']=self.tstop
        evtbin['dtime']=self.tbin
        evtbin.run()

    def exposure(self,src,fglName,gamma=None):
        """
        Compute exposure on source src, to add a flux column for the photometric light curve.

        Warning: the input file is modified in place !
        """

        if self.daily:
            infile=self.workDir+'/'+str(src)+'_daily_lc.fits'
            srcmdl=self.workDir+'/'+str(src)+'_daily.xml'
        else:
            infile=self.workDir+'/'+str(src)+'_lc.fits'
            srcmdl=self.workDir+'/'+str(src)+'.xml'

        scfile=self.spacecraft
        irfs='P7SOURCE_V6'
        target=fglName
        rad=str(self.roi)
        
        if gamma is None:
            options='infile='+infile+' scfile='+scfile+' irfs='+irfs+' srcmdl='+srcmdl+' target='+target+' rad='+rad
        else:
            options='infile='+infile+' scfile='+scfile+' irfs='+irfs+' srcmdl="none" specin='+str(gamma)+' rad='+rad
        os.system(self.fermiDir+'/bin/gtexposure '+options)


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

        import pyfits
        try:
            hdu=pyfits.open(infile)
        except:
            print 'Exception: can not open file '+infile
            raise
        data=hdu[1].data
        duration=data.field('TIMEDEL')[0]/3600./24. # sec -> days

        file=open(outfile,'w')
        file.write("#Time[MET]\tFlux[ph.cm^-2.s^-1]\tFluxError[ph.cm^-2.s^-1]\n")
        time      = data.field('TIME')     # MET
        counts    = data.field('COUNTS')
        countsErr = data.field('ERROR')    # error on counts
        exposure  = data.field('EXPOSURE') # cm^2 s^1
        flux      = counts/exposure        # approximate flux in ph cm^-2 s^-1
        fluxErr   = countsErr/exposure     # approximate flux error in ph cm^-2 s^-1
        
        for i in range(len(time)):
            if exposure[i] != 0.:
                file.write(str(time[i])+"\t"+str(flux[i])+"\t"+str(fluxErr[i])+"\n")
        file.close()


    def met2mjd(self,time):
        """
        Converts Mission Elapsed Time (MET, in seconds) in Modified Julian Day.
        Cf. http://fermi.gsfc.nasa.gov/ssc/data/analysis/documentation/Cicerone/Cicerone_Data/Time_in_ScienceTools.html
        to see how the time is handled in the Fermi Science Tools.
        
        Input: time in MET (s)
        Output: time in MJD (fraction of a day)
        """
        MJDREFI=51910.0
        MJDREFF=7.428703703703703e-4
        return(MJDREFI+MJDREFF+time/24./60./60.)

        
    def createPNG(self,src,fglName):
        """
        Create a PNG figure with the light curve of a given source.
        """

        # Read the .dat LC file
        if self.daily:
            infile=self.workDir+'/'+str(src)+'_daily_lc.dat'
            outfig=self.workDir+'/'+str(src)+'_daily_lc.png'
        else:
            infile=self.workDir+'/'+str(src)+'_lc.dat'
            outfig=self.workDir+'/'+str(src)+'_lc.png'

        data    = asciidata.open(infile)
        time    = data[0].tonumpy()
        flux    = data[1].tonumpy()
        fluxErr = data[2].tonumpy()

        fig=figure()
        ax = fig.add_subplot(111)
        ax.set_title(str(src)+', '+str(fglName).replace('_2FGLJ','2FGL J'))

        # Force the y-axis ticks to use 1e-6 as a base exponent
        ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: ('%.1f')%(x*1e6)))
        ax.set_ylabel('$F_{%.0f MeV-%.0f GeV}$ (x $10^{-6}$ ph cm$^{-2}$ s$^{-1}$)'%(self.emin,self.emax/1000.))

        # Make the x-axis ticks formatted to 0 decimal places
        day=24.*60.*60.
        toffset=54600
        time = self.met2mjd(time)  # Conversion MET -> MJD
        # We can do this because t is NOT a list, but a numpy.array

        ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: '%.0f'%(x-54600.)))
        ax.set_xlabel('MJD-'+str(toffset))
        
        # Plot a line at the threshold value
        axhline(y=self.threshold,linewidth=3,linestyle='--',color='r')

        # Plot a line at flux=0, for visibility/readibility
        axhline(y=0.,color='k')

        # Plot the light curve
        errorbar(x=time, y=flux, yerr=fluxErr/2., fmt='ro')
        

        # Don't show the figure in batch mode
        if not BATCH:
            show()
        ## Save the figure
        fig.savefig(outfig)


    def sendAlert(self,src):

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
        else:
            infile  = self.workDir+'/'+str(src)+'_lc.dat'
            pngFig=self.workDir+'/'+str(src)+'_lc.png'
        data    = asciidata.open(infile)
        time    = data[0].tonumpy()
        flux    = data[1].tonumpy()
        fluxErr = data[2].tonumpy()

        # Catch the last flux point
        lastFlux=flux[-1:]

        if DEBUG:
            print
            print "self.threshold=",self.threshold
            print "lastFlux=",lastFlux
            print
    
        # If the last flux is above the threshold, we send a mail
        if lastFlux >= self.threshold:
            # Create the container email message.
            msg = MIMEMultipart()
            msg['Subject'] = 'Fermi/LAT flare alert on %s' % src
            sender = 'Fermi automatic light curve robot <fermi@hess-lsw.lsw.uni-heidelberg.de>'
            
            #recipient = ['Gabriele Cologna <g.cologna@lsw.uni-heidelberg.de>',
            #             'Sarah Kaufmann <s.kaufmann@lsw.uni-heidelberg.de>',
            #             'Jean-Philippe Lenain <jp.lenain@lsw.uni-heidelberg.de>',
            #             'Mahmoud Mohamed <m.mohamed@lsw.uni-heidelberg.de>',
            #             'Stephanie Schwemmer <s.schwemmer@lsw.uni-heidelberg.de>',
            #             'Stefan Wagner <s.wagner@lsw.uni-heidelberg.de>']

            recipient = ['Jean-Philippe Lenain <jp.lenain@lsw.uni-heidelberg.de>']
            
            msg['From'] = sender
            COMMASPACE = ', '
            msg['To'] =COMMASPACE.join( recipient )
            msg.preamble = 'You will not see this in a MIME-aware mail reader.\n'
            # Guarantees the message ends in a newline
            msg.epilogue = ''
            
            mailtext="""
     *** The Fermi/LAT flux (%.0f MeV-%.0f GeV) of %s exceeds the trigger threshold of %.2g ph cm^-2 s^-1 ***

     The most recent lightcurve (%.0f-day binned) is attached.

     All available data can be found on 'hess-lsw' at
     %s
     
    """ %(self.emin,self.emax/1000.,src,self.threshold,self.tbin/24./60./60.,self.workDir)
            # (PS: The 'fermi' account password is the first name of the eponym famous physicist, in lower case ;-) ).
 
            txt = MIMEText(mailtext)
            msg.attach(txt)
            
            # Open the files in binary mode.  Let the MIMEImage class automatically
            # guess the specific image type.
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




def processSrc(mysrc=None,q=None,useThresh=False,daily=False):
    """
    Process a given source.
    """

    if DEBUG:
        print 'src=',mysrc
    
    if(mysrc != None):
        auto=autoLC(customThreshold=useThresh,daily=daily)
    else:
        print "ERROR Missing input source !"
        sys.exit(1)

    src,ra,dec,z,fglName=auto.readSourceList(mysrc)

    if q==None:
        auto.selectSrc(src,ra,dec)
        auto.makeTime(src,ra,dec)
        if fglName is not None:
            auto.createXML(src)
            mygamma=None
        else:
            mygamma=-2.5
            print 'Your source '+src+' has no 2FGL name in the list of sources. I will assume a photon index of '+str(mygamma)+' for the light curve generation.'
        auto.photoLC(src)
        auto.exposure(src,fglName,gamma=mygamma)
        auto.createDAT(src)
        auto.createPNG(src,fglName)
        auto.sendAlert(src)
    else:
        print "The MULTITHREAD flag is deprecated. Aborting..."
        return False
        #q.put([
        #        auto.selectSrc(src,ra,dec),
        #        auto.makeTime(src,ra,dec),
        #        auto.createXML(src),
        #        auto.photoLC(src),
        #        auto.exposure(src,fglName),
        #        auto.createDAT(src),
        #        auto.createPNG(src,fglName),
        #        auto.sendAlert(src)
        #        ])
    
    return True


def main(argv=None):
    """
    Main procedure
    """

    # options parser:
    helpmsg="""%prog [options] <source>

This is the version $Id$

Use '-h' to get the help message

"""

    parser = OptionParser(version="$Id$",
                          usage=helpmsg)

    parser.add_option("-d", "--daily", action="store_true", dest="d", default=False,
                      help='use daily bins for the light curves (defaulted to weekly)')
    parser.add_option("-c", "--custom-threshold", action="store_true", dest="c", default=False,
                      help='use custom trigger thresholds from the master list of sources (defaulted to 1.e-6 ph cm^-2 s^-1)')

    (opt, args) = parser.parse_args()


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

    # Check that we provided the mandatory argument: a source to process !
    if len(args) != 1:
        print "ERROR Main: wrong number of arguments"
        sys.exit(1)
    
    src=args[0]

    processSrc(mysrc=src,useThresh=USECUSTOMTHRESHOLD,daily=DAILY)

    return True


if __name__ == '__main__':
    """
    Execute main()
    """

    main()
