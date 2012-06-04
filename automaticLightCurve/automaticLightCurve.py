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
DEBUG=True



class autoLC:
    def __init__(self,file="/home/jplenain/fermi/automaticLightCurve/listSources.txt"):
        self.file=file

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


    def computePhotometricLC():
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

    auto.readSourceList()


if __name__ == '__main__':
    """
    Execute main()
    """

    main()
