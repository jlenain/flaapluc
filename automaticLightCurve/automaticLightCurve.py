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

# Flags
DEBUG=False



class autoLC:
    def __init__(self,file="/home/jplenain/fermi/automaticLightCurve/listSources.txt"):
        self.file=file

        # Import the Science Tools modules
        try:
            from gt_apps import *
        except ImportError:
            print "ERROR Can't import the Fermi Science tools"
            sys.exit(1)

        # Setting file names and directories
        self.allsky     = "/data/fermi/allsky/allsky_30MeV_300GeV_diffuse_filtered.fits"
        self.allskyGti  = "/data/fermi/allsky/allsky_30MeV_300GeV_diffuse_filtered_gti.fits"
        self.spacecraft = "/data/fermi/allsky/allsky_SC00.fits"
        self.workDir    = "/home/fermi/tmp"

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

        return zip(src,ra,dec,z)


    def selectSrc(src,ra,dec):
        """
        Filter a given source, running gtselect
        """
        filter['infile']=self.allskyGti
        filter['outfile']=self.workDir+'/'+str(src)+'_gti.fits'
        filter['ra']=ra
        filter['dec']=dec
        filter['rad']=10
        filter['emin']=100.
        filter['emax']=1.e5

        # Open allsky file to get the start and stop dates
        import pyfits
        try:
            hdu=pyfits.open(self.allskyGti)
        except:
            print 'Exception: can not open file '+allskyGti
            raise
        header=hdu[0].header
        
        filter['tmin']=header['TSTART']
        filter['tmax']=header['TSTOP']
        filter['zmax']=100.
        filter['evclass']=2
        filter.run()


    def makeTime(src):
        """
        Filter the GTI for a given source
        """
        print "Not yet implemented"
        


    def computePhotometricLC(src):
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

    srcLst=auto.readSourceList()
    src1=srcLst[0]
    auto.selectSrc(zip(*src1))
    


if __name__ == '__main__':
    """
    Execute main()
    """

    main()
