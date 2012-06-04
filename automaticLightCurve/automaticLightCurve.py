#!/usr/bin/env python
#
# @author $Author$ 
# Jean-Philippe Lenain, jplenain@lsw.uni-heidelberg.de
# @version $Id$

"""
Automatic generation of photometric light curves of Fermi sources.
"""

import sys, os, asciidata

def readSourceList(file="/home/jplenain/fermi/automaticLightSource/listSources.txt"):
    try:
        import asciidata
    except ImportError:
        print "ERROR Can't import asciidata, needed to read the list of sources. Aborting..."
        sys.exit(1)

    try:
        srcList=asciidata.open(file)
    except IOError:
        print "ERROR Can't open "+file
        sys.exit(1)

    src = srcList[0]
    ra  = srcList[1]
    dec = srcList[2]
    z   = srcList[3]
    
    if DEBUG:
        print "DEBUG src="+src[1]
        print "DEBUG ra ="+ra[1]
        print "DEBUG dec="+dec[1]
        print "DEBUG z  ="+z[1]




def main(argv=None):
    """
    Main procedure
    """

    argc    = len(sys.argv)
    argList = sys.argv
    
    if(argc<2):
        file=argList[0]
        print "Overriding default list of source: using "+file
        readSourceList(file)
    else:
        readSourceList


if __name__ == '__main__':
    """
    Execute main()
    """

    main()
