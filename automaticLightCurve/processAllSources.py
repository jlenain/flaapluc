#!/usr/bin/env python

"""
Process all sources for automatic aperture photometry of interesting 2FGL sources, with parametric batch jobs.

@author Jean-Philippe Lenain <mailto:jplenain@lsw.uni-heidelberg.de>
@date $Date$
@version $Id$
"""

import sys, os

# Import custom module
try:
    from automaticLightCurve import *
except ImportError:
    print "ERROR Can't import automaticLightCurve"
    sys.exit(1)

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
    # Number of sources
    nbSrc=len(src)

    print "I will process ",nbSrc," sources."
    print

    for i in range(nbSrc):
        cmd='echo "./automaticLightCurve "'+str(src[i]+' | batch') # | batch # to launch it
        os.system(cmd)




if __name__ == '__main__':
    """
    Execute main()
    """

    main()
