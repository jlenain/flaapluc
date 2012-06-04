#!/usr/bin/env python
#
# @author $Author$ 
# @version $Id$
# @date $Date$

"""
Automatic generation of photometric light curves of Fermi sources.
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
DEBUG=False



class autoLC:
    def __init__(self,file="/home/jplenain/fermi/automaticLightCurve/listSources.txt"):
        self.file=file

        # Setting file names and directories
        self.allsky     = "/data/fermi/allsky/allsky_30MeV_300GeV_diffuse_filtered.fits"
        self.spacecraft = "/data/fermi/allsky/allsky_SC00.fits"
        self.workDir    = "/home/fermi/tmp"

        # Setting default parameters
        self.roi  = 10.  # degrees
        self.emin = 1.e2 # E min
        self.emax = 1.e5 # E max

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
    
        if DEBUG:
            for i in range(len(src)):
                print "DEBUG src=",src[i]
                print "DEBUG ra =",ra[i]
                print "DEBUG dec=",dec[i]
                print "DEBUG z  =",z[i]

        return src,ra,dec,z


    def selectSrc(self,src,ra,dec):
        """
        Filter a given source, running gtselect
        """
        filter['infile']=self.allsky
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
        filter['zmax']=100.
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
            os.remove(outfile)

        maketime['filter']="IN_SAA!=T && LAT_CONFIG==1 && DATA_QUAL==1 && ABS(ROCK_ANGLE)<52 && ANGSEP("+ra+","+dec+",RA_SUN,DEC_SUN)+"+self.roi+">5."
        maketime['roicut']='yes'
        maketime['tstart']=self.tstart
        maketime['tstop']=self.tstop
        maketime['evfile']=self.workDir+'/'+str(src)+'.fits'
        maketime.run()


    def photomLC(src):
        """
        Compute the photometric light curve for a given source
        """
        print "Not yet implemented"




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

    src,ra,dec,z=auto.readSourceList()
    auto.selectSrc(src[0],ra[0],dec[0])
    auto.makeTime(src[0],ra[0],dec[0])
    


if __name__ == '__main__':
    """
    Execute main()
    """

    main()
