#!/usr/bin/python

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

    print 'First, we look for the ATOM schedule for today:'+infile
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


    # Now, we check that all the ATOM sources are known to the "master" list of sources
    FermiSrcInATOMSchedule=[]
    auto=autoLC()
    for i in range(len(ATOMsrcs)):
        ATOMsrc=ATOMsrcs[i]
        print 'Searching for the ATOM source '+ATOMsrc+' in the master list of sources...'
        tmpsrc,tmpra,tmpdec,tmpz,tmpfglName=auto.readSourceList(ATOMsrc)
        # Check if the ATOM source is in the master list of sources.
        if tmpsrc is None:
            print "Your source "+ATOMsrc+" can't be found in the master list of sources, skip it."
            print
            continue
        FermiSrcInATOMSchedule.append(tmpsrc)
    
    # Return a list of sources observed with ATOM last night, which are also known to the "master" list of sources
    return FermiSrcInATOMSchedule


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
                      help='only process those sources which are in the current ATOM schedule, forcing daily-binned data to be processed (this will also automatically create weekly-binned data)')
    parser.add_option("-d", "--daily", action="store_true", dest="d", default=False,
                      help='use daily bins for the light curves (defaulted to weekly)')
    parser.add_option("-c", "--custom-threshold", action="store_true", dest="c", default=False,
                      help='use custom trigger thresholds from the master list of sources (defaulted to 1.e-6 ph cm^-2 s^-1)')
    parser.add_option("-l", "--long-term", action="store_true", dest="l", default=False,
                      help='generate a long term light curve, using the whole mission time (defaulted to False)')
    parser.add_option("-n", "--no-mail", action="store_true", dest="n", default=False,
                      help='do not send alert mails')
    parser.add_option("-t", "--test", action="store_true", dest="t", default=False,
                      help='for test purposes. Do not send the alert mail to everybody if a source is above the trigger threshold, but only to J.-P. Lenain (by default, mail alerts are sent to everybody)')
    parser.add_option("--dry-run", action="store_true", dest="dryRun", default=False,
                      help='simulate only what the pipeline would do, forcing the use of ATOM sources, without actually processing any Fermi/LAT event. This is useful to see if the master list of sources is up-to-date with the ATOM sources in the current schedule.')

    (opt, args) = parser.parse_args()


    # If dry run
    if opt.dryRun:
        DRYRUN=True
        # Force opt.a to be True
        opt.a=True
    else:
        DRYRUN=False
        
    # If daily bins
    if opt.d:
        DAILY=True
    else:
        DAILY=False

    # If custom thresholds
    if opt.c:
        USECUSTOMTHRESHOLD=True
    else:
        USECUSTOMTHRESHOLD=False

    # If test mode
    if opt.t:
        TEST=True
    else:
        TEST=False

    # If no mail is sent
    if opt.n:
        MAIL=False
    else:
        MAIL=True

    # If long term
    if opt.l:
        LONGTERM=True
    else:
        LONGTERM=False
    
    
    if(len(args)!=0):
        file=args[0]
        print "Overriding default list of source: using "+file
        auto=autoLC(file,customThreshold=USECUSTOMTHRESHOLD,daily=DAILY,longTerm=LONGTERM)
    else:
        auto=autoLC(customThreshold=USECUSTOMTHRESHOLD,daily=DAILY,longTerm=LONGTERM)

    ATOMsrcsInSchedule=readATOMschedule()
    
    # If process only sources which are in ATOM schedule
    if opt.a:
        src=[]
        ra=[]
        dec=[]
        z=[]
        fglName=[]
        for i in range(len(ATOMsrcsInSchedule)):
            ATOMsrc=ATOMsrcsInSchedule[i]
            tmpsrc,tmpra,tmpdec,tmpz,tmpfglName=auto.readSourceList(ATOMsrc)
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
    
    # If dry run, we exit here
    if DRYRUN:
        return True


    ## OBSOLETE:
    ## Do it the dirty way, invoking os.system
    #for i in range(nbSrc):
    #    cmd='echo "./automaticLightCurve.py "'+str(src[i]+' | batch')
    #    os.system(cmd)
    


    if MULTITHREAD:
        print "The MULTITHREAD flag is deprecated. Aborting..."
        return False

        ## Use the multiprocessing Python module
        #from multiprocessing import Process, Queue
        #q=Queue()
        #proc=[]
        #for i in range(nbSrc):
        #    print 'Starting proc ',i,' for source ',src[i]
        #    proc.append(Process(target=processSrc, args=(src[i],q)))
        #    proc[i].start()

    else:

        # Or use the shell command "parallel"
        if PARALLEL:
            options=[]
            
            # Set the options for automaticLightCurve
            autoOptions=[]
            if USECUSTOMTHRESHOLD:
                autoOptions.append("-c")
            if DAILY:
                autoOptions.append("-d")
            if MAIL is False:
                autoOptions.append("-n")
            if TEST:
                autoOptions.append("-t")
            if LONGTERM:
                autoOptions.append("-l")

            # Loop on sources
            for i in range(nbSrc):
                # If source is in ATOM schedule and DAILY is False, force the creation of a daily-binned light curve
                if src[i] in ATOMsrcsInSchedule and DAILY is False:
                    autoOptions.append("-d")
                options.append('\"./automaticLightCurve.py '+' '.join(autoOptions)+' '+str(src[i])+'\"')
            cmd="parallel --jobs 8 ::: "+" ".join(options)
            # use --dry-run just to test the parallel command
            os.system(cmd)

        else:
            # Or directly process everything sequentially
            # Loop on sources
            for i in range(nbSrc):
                print 'Starting process ',i,' for source ',src[i]
                # If source is in the ATOM schedule and DAILY is False, force the creation of a daily-binned light curve. The corresponding weekly-binned data will be automatically recreated, if missing, by automaticLightCurve.py
                if src[i] in ATOMsrcsInSchedule and DAILY is False:
                    DAILY=True
                processSrc(mysrc=src[i],useThresh=USECUSTOMTHRESHOLD,daily=DAILY,mail=MAIL,longTerm=LONGTERM,test=TEST)
    
    return True



if __name__ == '__main__':
    """
    Execute main()
    """

    main()
