#!/bin/env python

"""
FLaapLUC: Fermi/LAT automatic aperture photometry Light C<->Urve

This module creates the web site presenting the results.

Highly inspired from Marcus Hauser's script 'create-result_website.py' for ATOM ADRAS.

@author Jean-Philippe Lenain
@version $Id$
@date $Date$
"""

from optparse import OptionParser
import os, sys, datetime, shutil

# Import custom module
try:
    from automaticLightCurve import automaticLightCurve as alc
except ImportError:
    print "ERROR Can't import automaticLightCurve"
    sys.exit(1)


def createResultWebsite(configfile='default.cfg'):
    
    # Retrieving information from the list of sources and from the configuration cfg file.
    auto=alc(configfile=configfile)
    autoLT=alc(longTerm=True,mergelongterm=True,configfile=configfile)

    # where all the graphs and data files are
    PATH_TO_FILES=auto.webpageDir+'/'

    WORKDIR=auto.workDir+'/'
    WORKDIRLT=autoLT.workDir+'/'
    src,ra,dec,z,fglName=auto.readSourceList()

    # Change umask, to make sure all files created hereafter in the webpage directory will be world-readable, in odrer that the web server can distribute the web page to anybody
    os.umask(0022)

    # create the directory if does not exist
    if not os.path.isdir(PATH_TO_FILES):
        os.makedirs(PATH_TO_FILES)

    # are we able to read and write there?
    if not (os.access(PATH_TO_FILES, os.W_OK | os.R_OK) ):
        sys.stderr.write("ERROR: path >%s< not read/writable!\n" % PATH_TO_FILES)
        sys.exit(1)

    # Remove everything from PATH_TO_FILES, for a fresh start
    for thisfile in os.listdir(PATH_TO_FILES):
        thisfilepath = os.path.join(PATH_TO_FILES,thisfile)
        os.remove(thisfilepath)



    # open file for writing
    f1 = open(PATH_TO_FILES+'index.html', 'w')

    f1.write("""
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
          "http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<title>FLaapLUC results</title>
  <meta name="ROBOTS" content="NOINDEX, NOFOLLOW">
  <link rel="stylesheet" href="../sty.css" type="text/css">
</head>
<body>


<div id="main">

<h1>The FLaapLUC (Fermi/LAT automatic aperture photometry Light C&harr;Urve) result page</h1>

<h2> Processing status</h2>

<p>
<strong style="color:green">system running</strong>, last data update was done at %s (UTC).
</p>


<h2> Results</h2>

Light curves:
<ul>
  <li>Horizontal red dotted line in the lightcurves shows the thresholds fixed for each source. By default, this is set to 10<sup>-6</sup> ph cm<sup>-2</sup> s<sup>-1</sup>, or is defined as long-term flux average over the whole mission.
  <li>Horizontal blue solid line in some plots represents the flux long-term average of a source.
  <li>Horizontal blue dotted lines in some plots represent the RMS around the flux long-term average.
</ul>

All data products are available from two sources: 

<ul> 
<li>the most recent results for the last 30 days of Fermi/LAT data, with weekly-binned light curves, and daily-binned as well if activity is seen.</li>

<li>a weekly-binned long-term light curve for the whole Fermi mission.</li>
</ul>

<!--
<p>
Moreover, if a source was on the ATOM schedule for last night, the most recent light curve will show both a weekly-binned light curve (blue data points) and a daily-binned light curve (red data points). If ATOM did not operate last night, the Fermi/LAT pipeline looks for ATOM schedule files up to 10 days in the past, in order to present here Fermi/LAT daily-binned light curves for the most recently observed ATOM sources.
</p>
-->
<p>
The ASCII data files
include all data points.
</p>

""" % (datetime.datetime.utcnow()) )
    
    
    ### print table header
    f1.write( """
<hr>

<div class="table-container">
<table border="1">
<thead border="1">
<tr>
	<th rowspan="2">Object name<br>(link to SIMBAD)</th>
	<th rowspan="2">RA, Dec (degrees)</th>
	<th colspan="2" style='background-color:#ccccff;'>Most recent data</th>

	<th colspan="2">Long-term data</th>
</tr>
<tr>
	<th style='background-color:#ccccff;'>Light Curve</th>
	<th style='background-color:#ccccff;'>ASCII data</th>

	<th>Light Curve</th>
	<th>ASCII data</th>
</tr>
</thead>
<tbody border="1">
""")
    
    
    
    
    ########################################################################
    # create main table:
    
    for i in range(len(src)):
        # Discard GPS sources
        if 'GPS' in src[i]:
            continue

        plotname         = src[i] + '_lc.png'
        dailyplotname    = src[i] + '_daily_lc.png'
        asciiname        = src[i] + '_lc.dat'
        dailyasciiname   = src[i] + '_daily_lc.dat'
        ltplotname       = 'lt_' + src[i] + '_lc.png'
        ltasciiname      = 'lt_' + src[i] + '_lc.dat'
        
        # copy files from WORKDIR to PATH_TO_FILES
        if os.path.isfile(WORKDIR+dailyplotname):
            shutil.copyfile(WORKDIR+dailyplotname,
                            PATH_TO_FILES + dailyplotname)
        else:
            if os.path.isfile(WORKDIR+plotname):
                shutil.copyfile(WORKDIR+plotname,
                                PATH_TO_FILES + plotname)
        if os.path.isfile(WORKDIR+asciiname):
            shutil.copyfile(WORKDIR+asciiname,
                            PATH_TO_FILES + asciiname)
        if os.path.isfile(WORKDIR+dailyasciiname):
            shutil.copyfile(WORKDIR+dailyasciiname,
                            PATH_TO_FILES + dailyasciiname)
        if os.path.isfile(WORKDIRLT+plotname):
            shutil.copyfile(WORKDIRLT+plotname,
                            PATH_TO_FILES + ltplotname)
        if os.path.isfile(WORKDIRLT+asciiname):
            shutil.copyfile(WORKDIRLT+asciiname,
                            PATH_TO_FILES + ltasciiname)

        linkSIMBAD = "http://simbad.u-strasbg.fr/simbad/sim-id?Ident=%s&NbIdent=1&submit=submit" % src[i]
        linkSIMBAD =  linkSIMBAD.replace('+', '%2b') # for URL conversion
        
        ### create one row for the HTML output table
        f1.write( '<tr> <!-- %s -->\n' % (src[i]) )
        f1.write( '    <td><a href="%s">%s</a></td>\n' % (linkSIMBAD, src[i]) )
        f1.write( '    <td style="font-size:smaller">%s %s</td>\n' %( ra[i], dec[i] ) )

        
        # most recent data
        # link the daily/weekly plot, if we have one (i.e. if the source was observed with ATOM last night)
        thisplot = plotname
        if os.path.isfile(PATH_TO_FILES + dailyplotname):
            thisplot = dailyplotname
        f1.write( '    <td style="background-color:#ccccff;"> <a href="%s"><img src="%s" width="128" alt=""></a></td>\n' % (thisplot, thisplot))

        txt = '    <td style="background-color:#ccccff;"> <a href="%s"> data </a>' % (asciiname)
        if os.path.isfile(PATH_TO_FILES + dailyasciiname):
            txt += '<br><a href="%s">daily data</a>' % (dailyasciiname)
        txt += '</td>\n'
        f1.write(txt)
        
        # long-term data
        f1.write( '    <td> <a href="%s"><img src="%s" width="128" alt=""></a></td>\n' % (ltplotname, ltplotname) )
        f1.write( '    <td> <a href="%s"> data </a></td>\n' % (ltasciiname) )
        f1.write( '</tr>\n\n' )

    
    ### static stuff at the end
    f1.write( """
</tbody>
</table>
</div>

<h3>Notes:</h3>
<ul>
  <li>This webpage is updated on a daily basis, around 8:00 (UTC).</li>

  <li>If you need more information, please send a mail to: <a href="mailto:jlenain@in2p3.fr">jlenain@in2p3.fr</a>
  </li>

</ul>

<div style="font-size:xx-small">
<hr>
last update:<br>
  <script type="text/javascript" >
     document.write(document.lastModified)
  </script>
</div>


</body>
</html>
""")
    
    f1.close()
    return True


def main(argv=None):
    """
    Main procedure
    """

    helpmsg="""create FLaapLUC result webpage.
"""
    
    parser = OptionParser(version="%prog:  $Id$",
                      usage=helpmsg)
    parser.add_option("--config-file", default='default.cfg', dest="CONFIGFILE", metavar="CONFIGFILE",
                                    help="provide a configuration file. Using '%default' by default.")
    (opt, args) = parser.parse_args()

    CONFIGFILE=opt.CONFIGFILE

    createResultWebsite(configfile=CONFIGFILE)

    return True


if __name__ == '__main__':
    """
    Execute main()
    """

    main()
