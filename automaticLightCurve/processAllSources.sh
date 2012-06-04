#!/bin/bash
#
#Process all sources for automatic aperture photometry of interesting 2FGL sources, with parametric batch jobs.
#
# @author Jean-Philippe Lenain <mailto:jplenain@lsw.uni-heidelberg.de>
# @date $Date$
# @version $Id$

listSrc="listSources.txt"
nbSrc=$(cat $listSrc | grep -v "#" | wc -l)

. ${FERMI_DIR}/fermi-init.sh
export PYTHONPATH=/usr/lib/python2.6/site-packages:$PYTHONPATH
#export PATH=/usr/local/heasoft-6.11/x86_64-unknown-linux-gnu-libc2.12.2/bin:/usr/kerberos/sbin:/usr/kerberos/bin:/usr/local/root/5.28.00/bin:/usr/local/bin:/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin:/home/jplenain/hess/bin:/usr/local/hess/bin:/opt/dell/srvadmin/bin:/home/jplenain/local/matlab/bin:/home/jplenain/bin

for src in `cat $listSrc | grep -v "#" | awk '{print $1}'`; do
    #nice ./automaticLightCurve.py $src &
    echo "./automaticLightCurve.py $src" | batch
    #disown
done
