#!/usr/bin/env python

"""
Process all sources for automatic aperture photometry of interesting 2FGL sources, with parametric batch jobs.

@author Jean-Philippe Lenain <mailto:jplenain@lsw.uni-heidelberg.de>
@date $Date$
@version $Id$
"""

import sys, os, datetime
from optparse import OptionParser

# Flags
MULTITHREAD=False
PARALLEL=True
USETHRESHOLD=True


# Import custom module
try:
    from automaticLightCurve import *
except ImportError:
    print "ERROR Can't import automaticLightCurve"
    sys.exit(1)


def readATOMschedule(infile='/home/fermi/'+datetime.date.today().strftime('%y%m%d')+'.sched'):
    """
    Read the ATOM schedule file for observations of tonight, automatically put on hess-lsw by "copy_schedule_to_attel.sh"

    By default, it takes as argument the schedule file of today.
    """

    if not os.path.isfile(infile):
        print "WARNING readATOMschedule: ",infile," does not exist. I will try to see if there is any ATOM schedule file for the last 10 days."
    
        found=False
        for i in range(1,10):
            infile='/home/fermi/'+(datetime.date.today()-datetime.timedelta(i)).strftime('%y%m%d')+'.sched'
            if os.path.isfile(infile):
                found=True
                print 'I found one ! I will use the file: '+infile
                print
            if found:
                break

        if not found:
            print 'ERROR readATOMschedule: I could not find any ATOM schedule file for the last 10 days. Aborting...'
            sys.exit(1)
        

    f=open(infile,'r')
    ATOMsrcs=[]
    for line in f:
        ATOMsrcs.append(
            (
                (
                    (
                        line.replace('\t','')
                        ).replace(' ','')
                    ).replace('\n','')
                )
            )
    
    # Return a list of sources observed with ATOM last night
    return ATOMsrcs


def main(argv=None):
    """
    Main procedure
    """

    # options parser:
    helpmsg="""%prog [options] <optional list of sources>

This is the version $Id$

Use '-h' to get the help message

If called with '-a', the list of sources will be taken from the last ATOM schedule file.

"""

    parser = OptionParser(version="$Id$",
                          usage=helpmsg)

    parser.add_option("-a", "--atom", action="store_true", dest="a", default=False,
                      help='use ATOM list of sources')
    parser.add_option("-d", "--daily", action="store_true", dest="d", default=False,
                      help='use daily bins for the light curves (defaulted to weekly)')

    (opt, args) = parser.parse_args()

        
    # If daily bins
    if opt.d:
        DAILY=True
    else:
        DAILY=False

    if(len(args)!=0):
        file=args[0]
        print "Overriding default list of source: using "+file
        auto=autoLC(file,customThreshold=USETHRESHOLD,daily=DAILY)
    else:
        auto=autoLC(customThreshold=USETHRESHOLD,daily=DAILY)


        
    # If use ATOM schedule
    if opt.a:
        src=[]
        ra=[]
        dec=[]
        z=[]
        fglName=[]
        ATOMsrcs=readATOMschedule()
        for i in range(len(ATOMsrcs)):
            ATOMsrc=ATOMsrcs[i]
            print 'Searching for the ATOM source '+ATOMsrc+' in the master list of sources...'
            tmpsrc,tmpra,tmpdec,tmpz,tmpfglName=auto.readSourceList(ATOMsrc)
            # Check if the ATOM source is in the master list of sources.
            if tmpsrc is None:
                print "Your source "+ATOMsrc+" can't be found in the master list of sources, skip it."
                print
                continue
            src.append(tmpsrc)
            ra.append(tmpra)
            dec.append(tmpdec)
            z.append(tmpz)
            fglName.append(tmpfglName)
        
    else:
        src,ra,dec,z,fglName=auto.readSourceList()

    # Total number of sources to process
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

        # Or use the shell command "parallel"
        if PARALLEL:
            options=[]
            for i in range(nbSrc):
                options.append('\"./automaticLightCurve.py '+str(src[i])+'\"')
            cmd="parallel --jobs 8 ::: "+" ".join(options)
            # use --dry-run just to test the parallel command
            os.system(cmd)

        else:
            # Or directly process everything sequentially
            for i in range(nbSrc):
                print 'Starting process ',i,' for source ',src[i]
                processSrc(src[i])
    
    return True



if __name__ == '__main__':
    """
    Execute main()
    """

    main()