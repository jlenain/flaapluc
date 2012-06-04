#!/usr/bin/env python

"""
Automatic generation of aperture photometric light curves of Fermi sources, for a given source.

No likelihood fit is performed, the results solely rely on the 2FGL spectral fits.

More information are available at: http://fermi.gsfc.nasa.gov/ssc/data/analysis/scitools/aperture_photometry.html

@author Jean-Philippe Lenain <mailto:jplenain@lsw.uni-heidelberg.de>
@date $Date$
@version $Id$
"""

import sys, os, asciidata
from numpy import *

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

    def __init__(self,file="/home/jplenain/fermi/automaticLightCurve/listSources.txt"):
        self.file=file

        # Setting file names and directories
        self.allsky     = "/data/fermi/allsky/allsky_30MeV_300GeV_diffuse_filtered.fits"
        self.spacecraft = "/data/fermi/allsky/allsky_SC00.fits"
        self.workDir    = "/home/fermi/tmp"
        self.fermiDir   = os.getenv('FERMI_DIR')

        # Setting default parameters
        self.roi  = 1.   # degrees (http://fermi.gsfc.nasa.gov/ssc/data/analysis/scitools/aperture_photometry.html: "For aperture photometry we select a very small aperture (rad=1 degree), because we are not fitting the background.")
        self.emin = 1.e2 # E min
        self.emax = 1.e5 # E max
        self.zmax = 100.

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


    def readSourceList(self,mysrc):
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
    
        # Find our input src in the list of sources
        found=False
        for i in range(len(src)):
            if src[i]==mysrc:
                found=True
                return src[i],ra[i],dec[i],z[i],fglName[i]
            
        # If we end up without any found source, print an error and exits
        print "ERROR Can't find your source "+str(mysrc)+" in the list of sources !"
        sys.exit(1)


    def selectSrc(self,src,ra,dec):
        """
        Filter a given source, running gtselect
        """
        filter['infile']=self.allsky
        outfile=self.workDir+'/'+str(src)+'.fits'
        filter['outfile']=outfile
        
        # If outfile exsits, we remove it before updating it
        if os.path.isfile(outfile):
            if DEBUG:
                return
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

        outfile=self.workDir+'/'+str(src)+'_gti.fits'
        maketime['outfile']=outfile
        
        # If outfile exsits, we remove it before updating it
        if os.path.isfile(outfile):
            if DEBUG:
                return
            os.remove(outfile)

        # cf. http://fermi.gsfc.nasa.gov/ssc/data/analysis/scitools/aperture_photometry.html
        maketime['filter']="IN_SAA!=T && LAT_CONFIG==1 && DATA_QUAL==1 && (angsep(RA_ZENITH,DEC_ZENITH,"+str(ra)+","+str(dec)+")+"+str(self.roi)+" <"+str(self.zmax)+") && (angsep("+str(ra)+","+str(dec)+",RA_SCZ,DEC_SCZ)<180.)"
        maketime['roicut']='no'
        maketime['tstart']=self.tstart
        maketime['tstop']=self.tstop
        maketime['evfile']=self.workDir+'/'+str(src)+'.fits'
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
        evtbin['dtime']=604800 # sec. = 1 week
        evtbin.run()

    def exposure(self,src,fglName):
        """
        Compute exposure on source src, to add a flux column for the photometric light curve.

        Warning: the input file is modified in place !
        """

        infile=self.workDir+'/'+str(src)+'_lc.fits'
        scfile=self.spacecraft
        irfs='P7SOURCE_V6'
        srcmdl=self.workDir+'/'+str(src)+'.xml'
        target=fglName
        rad=str(self.roi)
        
        options='infile='+infile+' scfile='+scfile+' irfs='+irfs+' srcmdl='+srcmdl+' target='+target+' rad='+rad
        os.system(self.fermiDir+'/bin/gtexposure '+options)


    def createDAT(self,src):
        """
        Create a data file with the light curve of a given source.
        """

        # Read LC file
        infile=self.workDir+'/'+str(src)+'_lc.fits'

        import pyfits
        try:
            hdu=pyfits.open(infile)
        except:
            print 'Exception: can not open file '+infile
            raise
        data=hdu[1].data
        duration=data.field('TIMEDEL')[0]/3600./24. # sec -> days

        outfile=self.workDir+'/'+str(src)+'_lc.dat'
        file=open(outfile,'w')
        file.write("#Time[MET]\tFlux[ph.cm^-2.s^-1]\tFluxError[ph.cm^-2.s^-1]\n")
        time      = data.field('TIME')     # MET
        counts    = data.field('COUNTS')
        countsErr = data.field('ERROR')    # error on counts
        exposure  = data.field('EXPOSURE') # cm^2 s^1
        flux      = counts/exposure        # approximate flux in ph cm^-2 s^-1
        fluxErr   = countsErr/exposure     # approximate flux error in ph cm^-2 s^-1
        
        for i in range(len(time)):
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

        try:
            from matplotlib.pyplot import *
            from matplotlib.ticker import FuncFormatter
        except ImportError:
            print "ERROR Can't import matplotlib"
            sys.exit(1)

        # Read the .dat LC file
        infile  = self.workDir+'/'+str(src)+'_lc.dat'
        data    = asciidata.open(infile)
        time    = data[0].tonumpy()
        flux    = data[1].tonumpy()
        fluxErr = data[2].tonumpy()

        outfig=self.workDir+'/'+str(src)+'_lc.png'
        fig=figure()
        ax = fig.add_subplot(111)
        ax.set_title(str(src)+', '+str(fglName).replace('_2FGLJ','2FGL J'))

        # Force the y-axis ticks to use 1e-8 as a base exponent
        ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: ('%.1f')%(x*1e8)))
        ax.set_ylabel('$F_{100 MeV-100 GeV}$ (x $10^{-8}$ ph cm$^{-2}$ s$^{-1}$)')

        # Make the x-axis ticks formatted to 0 decimal places
        day=24.*60.*60.
        toffset=54600
        time = self.met2mjd(time)  # Conversion MET -> MJD
        # We can do this because t is NOT a list, but a numpy.array

        ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: '%.0f'%(x-54600.)))
        ax.set_xlabel('MJD-'+str(toffset))
        
        errorbar(x=time, y=flux, yerr=fluxErr, fmt='ro')
        
        threshold=1.e-6 # ph cm^-2 s^-1
        axhline(y=threshold,color='k')

        # Don't show the figure in batch mode
        if BATCH is False:
            show()
        ## Save the figure
        fig.savefig(outfig)



def main(argv=None):
    """
    Main procedure
    """

    argc    = len(sys.argv)
    argList = sys.argv
    
    if(argc==2):
    #    file=argList[0]
    #    print "Overriding default list of source: using "+file
    #    auto=autoLC(file)
        src=argList[1]
        auto=autoLC()
    else:
        print "ERROR Missing input source !"
        sys.exit(1)

    src,ra,dec,z,fglName=auto.readSourceList(src)

    auto.selectSrc(src,ra,dec)
    auto.makeTime(src,ra,dec)
    auto.createXML(src)
    auto.photoLC(src)
    auto.exposure(src,fglName)
    auto.createDAT(src)
    auto.createPNG(src,fglName)


if __name__ == '__main__':
    """
    Execute main()
    """

    main()
