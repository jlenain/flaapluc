#!/bin/env python

## SGE parameters (this should work at CCIN2P3):
#$ -N processAllSources
#$ -e /sps/hess/lpnhe/jlenain/fermi/log
#$ -o /sps/hess/lpnhe/jlenain/fermi/log
#$ -j yes
#$ -l ct=02:00:00
#$ -l mem_free=2048M
#$ -l fsize=2048M
#$ -l xrootd=0
#$ -l sps=1
#$ -l os=sl6
#$ -p 0
#$ -P P_hess -cwd


"""
Process all sources for automatic aperture photometry of interesting high energy sources, with parametric batch jobs.

@author Jean-Philippe Lenain <mailto:jlenain@in2p3.fr>
"""

import sys, os, datetime
from optparse import OptionParser
from ConfigParser import ConfigParser

# Import custom module
try:
    from automaticLightCurve import *
except ImportError:
    print "\033[91mERROR\033[0m Can't import automaticLightCurve"
    sys.exit(1)


def getConfig(configfile='./default.cfg'):
    config = ConfigParser()
    config.readfp(open(configfile))
    return config


def readATOMschedule(configfile='default.cfg'):
    """
    Read the ATOM schedule file for observations of tonight, automatically put on hess-lsw@lsw.uni-heidelberg.de by "copy_schedule_to_attel.sh" in ATOM pipeline.

    By default, it takes as argument the schedule file of today.
    """

    # Read configuration file
    config = getConfig(configfile=configfile)

    ATOMSchedulesDir = config.get('InputDirs','ATOMSchedulesDir')
    ATOMScheduleFile = datetime.date.today().strftime('%y%m%d')+'.sched'

    infile=ATOMSchedulesDir+'/'+ATOMScheduleFile

    print 'First, we look for the ATOM schedule for today: %s' % infile
    if not os.path.isfile(infile):
        print "\033[93mWARNING\033[0m readATOMschedule: %s does not exist. I will try to see if there is any ATOM schedule file for the last 10 days." % infile
    
        found=False
        # scan of ATOM schedues over the last 10 days
        for i in range(1,10):
            infile=ATOMSchedulesDir+'/'+(datetime.date.today()-datetime.timedelta(i)).strftime('%y%m%d')+'.sched'
            if os.path.isfile(infile):
                found=True
                print 'INFO readATOMschedule: I found an ATOM schedule ! I will use the file: %s' % infile
                print
            if found:
                break

        if not found:
            print '\033[93mWARNING\033[0m readATOMschedule: I could not find any ATOM schedule file for the last 10 days. I will not process any daily-binned data...'
            # Return an empty list
            return []
        

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
    auto=autoLC(configfile=configfile)
    for i in range(len(ATOMsrcs)):
        ATOMsrc=ATOMsrcs[i]
        print 'Searching for the ATOM source %s in the master list of sources...' % ATOMsrc
        tmpsrc,tmpra,tmpdec,tmpz,tmpfglName=auto.readSourceList(ATOMsrc)
        # Check if the ATOM source is in the master list of sources.
        if tmpsrc is None:
            print "\033[93mYour source can't be found in the master list of sources, skipping it.\033[0m" % ATOMsrc
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
                      help='only process those sources which are in the current ATOM schedule, forcing daily-binned data to be processed (this will also automatically create long time-binned data).')
    parser.add_option("-d", "--daily", action="store_true", dest="d", default=False,
                      help='use daily bins for the light curves (defaulted to weekly).')
    parser.add_option("-c", "--custom-threshold", action="store_true", dest="c", default=False,
                      help='use custom trigger thresholds provided in the master list of sources (defaulted to 1.e-6 ph cm^-2 s^-1).')
    parser.add_option("-l", "--long-term", action="store_true", dest="l", default=False,
                      help='generate a long term light curve, using the whole mission time (defaulted to False).')
    parser.add_option("-m", "--merge-long-term", action="store_true", dest="m", default=False,
                      help='merge the long-term month-by-month light curves together.')
    parser.add_option("-u","--update",action="store_true",dest="u",default=False,
                      help='update with new data for last month/year when used in conjunction of --merge-long-term. Otherwise, has no effect.')
    parser.add_option("-w", "--with-history", action="store_true", dest="history", default=False,
                      help='use the long-term history of a source to dynamically determine a flux trigger threshold, instead of using a fixed flux trigger threshold as done by default. This option makes use of the long-term data on a source, and assumes that these have been previously generated with the --merge-long-term option.')
    parser.add_option("-n", "--no-mail", action="store_true", dest="n", default=False,
                      help='do not send alert mails')
    parser.add_option("-t", "--test", action="store_true", dest="t", default=False,
                      help='for test purposes. Do not send the alert mail to everybody if a source is above the trigger threshold, but only to test recipients (by default, mail alerts are sent to everybody, cf. the configuration files).')
    parser.add_option("--dry-run", action="store_true", dest="dryRun", default=False,
                      help='only simulate what the pipeline would do, forcing the use of ATOM sources, without actually processing any Fermi/LAT event. This is useful to see if the master list of sources is up-to-date with the ATOM sources in the current schedule.')
    parser.add_option("--use-parallel", action="store_true", default=False, dest="parallel",
                      help="use the linux program 'parallel' to use multiple threads on a local server.")
    parser.add_option("--max-cpu", default=1, dest="MAXCPU", metavar="<MAXCPU>",
                      help="in conjunction with --use-parallel, defines the number of CPU to use. Using '%default' by default.")
    parser.add_option("-f", "--config-file", default='default.cfg', dest="CONFIGFILE", metavar="CONFIGFILE",
                      help="provide a configuration file. Using '%default' by default.")
    parser.add_option("-v", "--verbose", action="store_true", dest="v", default=False,
                      help='verbose output.')

    (opt, args) = parser.parse_args()

    CONFIGFILE=opt.CONFIGFILE

    global VERBOSE
    if opt.v:
        VERBOSE=True
    else:
        VERBOSE=False

    # If dry run
    if opt.dryRun:
        DRYRUN=True
    else:
        DRYRUN=False

    # If parallel
    if opt.parallel:
        PARALLEL=True
        if int(opt.MAXCPU) != 1:
            MAXCPU=int(opt.MAXCPU)
        else:
            # Force the script to process only one source at a time
            MAXCPU=1
    else:
        PARALLEL=False
        
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

    if TEST and not MAIL:
        print "\033[91mFATAL You asked for both the --test and --no-mail options."
        print "      These are mutually exclusive options.\033[0m"
        sys.exit(1)

    # If long term
    if opt.l:
        LONGTERM=True
    else:
        LONGTERM=False
    
    if LONGTERM is True and DAILY is True:
        print "\033[91mFATAL You asked for both the --long-term and --daily options."
        print "      Since this is too CPU intensive, this combination was disabled.\033[0m"
        sys.exit(1)


    # If merge long term light curves
    if opt.m:
        MERGELONGTERM=True
        # If update long term light curves with brand new data
        if opt.u:
            UPDATE=True
        else:
            UPDATE=False        
    else:
        MERGELONGTERM=False
        UPDATE=False        

    # If dynamical flux trigger threshold based on source history
    if opt.history:
        WITHHISTORY=True
    else:
        WITHHISTORY=False


    if(len(args)!=0):
        file=args[0]
        print "\033[93mINFO\033[0m Overriding default list of source: using %s" % file
        auto=autoLC(file=file,customThreshold=USECUSTOMTHRESHOLD,daily=DAILY,longTerm=LONGTERM,mergelongterm=MERGELONGTERM,withhistory=WITHHISTORY,configfile=CONFIGFILE)
    else:
        auto=autoLC(customThreshold=USECUSTOMTHRESHOLD,daily=DAILY,longTerm=LONGTERM,mergelongterm=MERGELONGTERM,withhistory=WITHHISTORY,update=UPDATE,configfile=CONFIGFILE)

    ATOMsrcsInSchedule=readATOMschedule(configfile=CONFIGFILE)
    
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

    print
    print "I will process %i sources." % nbSrc
    print
    


    # Use the shell command "parallel"
    if PARALLEL:
        options=[]
        
        # Set the options for automaticLightCurve
        autoOptions=[]
        autoOptions.append("--config-file="+CONFIGFILE)
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
        if MERGELONGTERM:
            autoOptions.append("-m")
            if UPDATE:
                autoOptions.append("-u")
        if WITHHISTORY:
            autoOptions.append("--with-history")

        # Loop on sources
        for i in range(nbSrc):
            # If source is in ATOM schedule and DAILY is False, force the creation of a daily-binned light curve
            if src[i] in ATOMsrcsInSchedule and DAILY is False and LONGTERM is False and MERGELONGTERM is False:
                # Put the -d option only for this source
                options.append('\"/bin/nice -n 10 ./automaticLightCurve.py -d '+' '.join(autoOptions)+' '+str(src[i])+'\"')
            else:
                options.append('\"/bin/nice -n 10 ./automaticLightCurve.py '+' '.join(autoOptions)+' '+str(src[i])+'\"')
        cmd="parallel --jobs "+str(MAXCPU)+" ::: "+" ".join(options)
        # use --dry-run just to test the parallel command
        if not DRYRUN:
            os.system(cmd)
        else:
            print cmd

    else:
        # Or directly process everything sequentially, using only 1 CPU
        if MERGELONGTERM:
            LONGTERM=True
        
        # Loop on sources
        for i in range(nbSrc):
            print 'Starting process %i for source %s' % (i,src[i])
            # If source is in the ATOM schedule and DAILY is False, force the creation of a daily-binned light curve.
            if src[i] in ATOMsrcsInSchedule and not DAILY and not LONGTERM and not MERGELONGTERM:
                tmpDAILY=True
                # We have to make sure that the corresponding weekly-binned data are created first (needed for daily PNG figure)
                if not DRYRUN:
                    processSrc(mysrc=src[i],useThresh=USECUSTOMTHRESHOLD,daily=False,mail=False,longTerm=LONGTERM,mergelongterm=MERGELONGTERM,withhistory=WITHHISTORY,update=UPDATE,configfile=CONFIGFILE)
                else:
                    print "processSrc(mysrc="+src[i]+",useThresh="+str(USECUSTOMTHRESHOLD)+",daily=False,mail=False,longTerm="+str(LONGTERM)+",mergelongterm="+str(MERGELONGTERM)+",withhistory="+str(WITHHISTORY)+",update="+str(UPDATE)+",configfile="+str(CONFIGFILE)+")"
            else:
                tmpDAILY=DAILY

            if not DRYRUN:
                processSrc(mysrc=src[i],useThresh=USECUSTOMTHRESHOLD,daily=tmpDAILY,mail=MAIL,longTerm=LONGTERM,test=TEST,mergelongterm=MERGELONGTERM,withhistory=WITHHISTORY,update=UPDATE,configfile=CONFIGFILE)
            else:
                print "processSrc(mysrc="+src[i]+",useThresh="+str(USECUSTOMTHRESHOLD)+",daily="+str(tmpDAILY)+",mail="+str(MAIL)+",longTerm="+str(LONGTERM)+",test="+str(TEST)+",mergelongterm="+str(MERGELONGTERM)+",withhistory="+str(WITHHISTORY)+",update="+str(UPDATE)+",configfile="+str(CONFIGFILE)+")"
            print

    return True



if __name__ == '__main__':
    """
    Execute main()
    """

    main()
