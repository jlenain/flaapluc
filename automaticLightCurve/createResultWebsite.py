#! /usr/bin/python

"""
Fermi/LAT light curve pipeline: create website for presenting results.

Highly inspired from Marcus Hauser's script 'create-result_website.py' for ATOM ADRAS.

@author Jean-Philippe Lenain
@version $Id$
@date $Date$
"""

from optparse import OptionParser
import os, sys, datetime, shutil

# Import custom module
try:
    from automaticLightCurve import *
except ImportError:
    print "ERROR Can't import automaticLightCurve"
    sys.exit(1)



# where all the graphs and data files are
PATH_TO_FILES="/var/www/html/results/"

# create the directory if does not exist
if not os.path.isdir(PATH_TO_FILES):
    os.makedirs(PATH_TO_FILES)

# are we able to read and write there?
if not (os.access(PATH_TO_FILES, os.W_OK | os.R_OK) ):
    sys.stderr.write("ERROR: path >%s< not read/writable!\n" % PATH_TO_FILES)
    sys.exit(1)


def createResultWebsite():
    
    # open file for writing
    f1 = open(PATH_TO_FILES+'index.html', 'w')

    f1.write("""
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
          "http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<title>Fermi/LAT light curve results</title>
  <meta name="ROBOTS" content="NOINDEX, NOFOLLOW">
  <link rel="stylesheet" href="../sty.css" type="text/css">
</head>
<body>


<div id="main">

<h1>The Fermi/LAT light curve result page at the Landessternwarte</h1>

<h2> Processing status</h2>

<p>
<strong style="color:green">system running</strong>, last data update was done at %s (UTC).
</p>


<h2> Results</h2>

Light curves:
<ul>
  <li>Horizontal red dotted lines in the lightcurves show the thresholds fixed for each source. By default, this is set to 1.e-6 ph cm^-2 s^-1.
</ul>

All data products are available from two sources: 

<ul> 
<li>First, the most recent results for the last 70 days of Fermi/LAT data.</li>

<li>Secondly, the long-term light curve for the whole Fermi mission, for each source.</li>
</ul>

<p>
The ASCII data files
include all data points.
</p>

""" % (datetime.datetime.utcnow()) )
    
    
    ### print table header
    f1.write( """
<hr>
<table border="1" style="font-size:70%">
<tr>
	<th rowspan="2">object name<br>(link to SIMBAD)</th>
	<th rowspan="2">RA / Dec (degrees)</th>
	<th colspan="2" style='background-color:#e7e5bc;'>most recent data</th>

	<th colspan="2">long-term data</th>
</tr>
<tr>
	<th style='background-color:#e7e5bc;'>LC</th>
	<th style='background-color:#e7e5bc;'>ASCII data</th>

	<th>LC</th>
	<th>ASCII data</th>
</tr>
""")
    
    
    # Retrieving information from the list of sources
    auto=autoLC()
    autoLT=autoLC(longTerm=True,mergelongterm=True)
    WORKDIR=auto.workDir+'/'
    WORKDIRLT=autoLT.workDir+'/'
    src,ra,dec,z,fglName=auto.readSourceList()
    
    
    ########################################################################
    # create main table:
    
    for i in range(len(src)):
        plotname      = src[i] + '_lc.png'
        dailyplotname = src[i] + '_daily_lc.png'
        asciiname     = src[i] + '_lc.dat'
        ltplotname    = 'lt_'+src[i] + '_lc.png'
        ltasciiname   = 'lt_'+src[i] + '_lc.dat'
        
        # copy files from WORKDIR to PATH_TO_FILES
        if os.path.isfile(WORKDIR+dailyplotname):
            shutil.copyfile( WORKDIR+dailyplotname,
                             PATH_TO_FILES + dailyplotname)
        else:
            if os.path.isfile(WORKDIR+plotname):
                shutil.copyfile( WORKDIR+plotname,
                                 PATH_TO_FILES + plotname)
        if os.path.isfile(WORKDIR+asciiname):
            shutil.copyfile( WORKDIR+asciiname,
                             PATH_TO_FILES + asciiname)
        if os.path.isfile(WORKDIRLT+plotname):
            shutil.copyfile( WORKDIRLT+plotname,
                             PATH_TO_FILES + ltplotname)
        if os.path.isfile(WORKDIRLT+asciiname):
            shutil.copyfile( WORKDIRLT+asciiname,
                             PATH_TO_FILES + ltasciiname)

        linkSIMBAD = "http://simbad.u-strasbg.fr/simbad/sim-id?Ident=%s&NbIdent=1&submit=submit" % src[i]
        linkSIMBAD =  linkSIMBAD.replace('+', '%2b') # for URL conversion
        
        ### create one row for the HTML output table
        f1.write( '<tr> <!-- %s -->\n' % (src[i]) )
        f1.write( '    <td><a href="%s">%s</a></td>\n' % (linkSIMBAD, src[i]) )
        f1.write( '    <td style="font-size:smaller">%s %s</td>\n' %( ra[i], dec[i] ) )
        
        # most recent data
        # link the daily/weekly plot, if we have one (i.e. if the source was observed with ATOM last night)
        if os.path.isfile(PATH_TO_FILES + dailyplotname):
            f1.write( '    <td style="background-color:#e7e5bc;"> <a href="%s">PNG</a></td>\n' % (dailyplotname) )
        else:
            f1.write( '    <td style="background-color:#e7e5bc;"> <a href="%s">PNG</a></td>\n' % (plotname) )
        f1.write( '    <td style="background-color:#e7e5bc;"> <a href="%s"> data </a></td>\n' % (asciiname) )
        
        # long-term data
        f1.write( '    <td> <a href="%s">PNG</a></td>\n' % (ltplotname) )
        f1.write( '    <td> <a href="%s"> data </a></td>\n' % (ltasciiname) )
        f1.write( '</tr>\n\n' )
        
    
    ### static stuff at the end
    f1.write( """
</table>

<h3>Notes:</h3>
<ul>
  <li>This webpage is updated daily.</li>

  <li>If you need more information, please send a mail to: <a href="mailto:jp.lenain@lsw.uni-heidelberg.de">jp.lenain@lsw.uni-heidelberg.de</a>
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

    helpmsg="""create Fermi/LAT light curve result webpage.

Note: due to different reasons, the PATHNAME for all created files is FIXED to %s !
""" % PATH_TO_FILES
    
    parser = OptionParser(version="%prog:  $Id$",
                      usage=helpmsg)

    (opt, args) = parser.parse_args()

    createResultWebsite()

    return True


if __name__ == '__main__':
    """
    Execute main()
    """

    main()
