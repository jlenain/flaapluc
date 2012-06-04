#!/usr/bin/env python

"""
Automatic generation of aperture photometric light curves of Fermi sources.

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
DEBUG=True



class autoLC:
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


    def readSourceList(self):
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
    
        if DEBUG:
            for i in range(len(src)):
                print "DEBUG src=",src[i]
                print "DEBUG ra =",ra[i]
                print "DEBUG dec=",dec[i]
                print "DEBUG z  =",z[i]
                print "DEBUG fglName=",fglName[i]

        return src,ra,dec,z,fglName



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

        try:
            import matplotlib.pyplot as plt
        except ImportError:
            print "ERROR Can't import matplotlib"
            sys.exit(1)

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
        for value in data:
            time      = data.field('TIME')
            counts    = data.field('COUNTS')
            countsErr = data.field('ERROR') # error on counts
            expo      = data.field('EXPOSURE') # cm^2 s^1
            flux      = counts/exposure # approximate flux in ph cm^-2 s^-1
            fluxErr   = countsErr/exposure # approximate flux error in ph cm^-2 s^-1
            file.write("%8d\t%5.5e\t%5.5e\n")%(time,flux,fluxErr)
        file.close()

        
    def createPNG(self,src):
        """
        Create a PNG figure with the light curve of a given source.
        """




def main(argv=None):
    """
    Main procedure
    """

    argc    = len(sys.argv)
    argList = sys.argv
    
    if(argc==2):
        file=argList[0]
        print "Overriding default list of source: using "+file
        auto=autoLC(file)
    else:
        auto=autoLC()

    src,ra,dec,z,fglName=auto.readSourceList()
    auto.selectSrc(src[0],ra[0],dec[0])
    auto.makeTime(src[0],ra[0],dec[0])
    auto.createXML(src[0])
    auto.photoLC(src[0])
    auto.exposure(src[0],fglName[0])


if __name__ == '__main__':
    """
    Execute main()
    """

    main()
