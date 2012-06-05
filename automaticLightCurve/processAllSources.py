#!/usr/bin/env python

"""
Process all sources for automatic aperture photometry of interesting 2FGL sources, with parametric batch jobs.

@author Jean-Philippe Lenain <mailto:jplenain@lsw.uni-heidelberg.de>
@date $Date$
@version $Id$
"""

import sys, os

# Flags
MULTITHREAD=True


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

    
    ## Do it the dirty way, invoking os.system
    #for i in range(nbSrc):
    #    cmd='echo "./automaticLightCurve.py "'+str(src[i]+' | batch')
    #    os.system(cmd)
    

    if MULTITHREAD:

        # Use the multiprocessing Python module
        from multiprocessing import Process, Queue
        q=Queue()
        proc=[]
        for i in range(nbSrc):
            print 'Starting proc ',i,' for source ',src[i]
            proc.append(Process(target=processSrc, args=(src[i],q)))
            proc[i].start()

    else:

        # Or directly process everything sequentially
        for i in range(nbSrc):
            print 'Starting process ',i,' for source ',src[i]
            processSrc(src[i],q)
    
    return True



if __name__ == '__main__':
    """
    Execute main()
    """

    main()
