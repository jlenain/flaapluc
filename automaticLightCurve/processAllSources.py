#!/usr/bin/env python

"""
Process all sources for automatic aperture photometry of interesting 2FGL sources, with parametric batch jobs.

@author Jean-Philippe Lenain <mailto:jplenain@lsw.uni-heidelberg.de>
@date $Date$
@version $Id$
"""

import sys, os

#os.environ['PYTHONPATH']='/usr/local/fermi/src/ScienceTools-v9r27p1-fssc-20120410/external/x86_64-unknown-linux-gnu-libc2.12.2/lib/python2.7/site-packages:'+os.getenv('PYTHONPATH')
#os.environ['PATH']='/usr/kerberos/sbin:/usr/kerberos/bin:/usr/local/root/5.28.00/bin:/usr/local/bin:/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin:/home/jplenain/hess/bin:/usr/local/hess/bin:/opt/dell/srvadmin/bin:/home/jplenain/local/matlab/bin:/home/jplenain/bin'

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

    
    #for i in range(nbSrc):
    #    cmd='echo "./automaticLightCurve.py "'+str(src[i]+' | batch')
    #    os.system(cmd)
    
    from multiprocessing import Process, Queue
    q=Queue()
    proc=[]
    for i in range(nbSrc):
        p[i]=Process(target=processSrc, args=(src[i]))
        p[i].start()




if __name__ == '__main__':
    """
    Execute main()
    """

    main()
