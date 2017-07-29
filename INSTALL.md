# Installation

It is recommended to install FLaapLUC under a python environment using e.g. [conda](https://www.continuum.io/downloads). The procedure is then the following:

```conda create -n flaapluc python=2.7 astropy matplotlib ephem
source activate flaapluc
git clone https://github.com/jlenain/flaapluc
cd flaapluc
python setup.py install
```

# Configuration

FLaapLUC needs to use the Fermi-LAT ScienceTools as well as [Enrico](http://enrico.readthedocs.io), which should be pre-installed. The configuration is then the following:

```
. $FERMI_DIR/fermi-init.sh
. $ENRICO_DIR/enrico-init.sh
source activate flaapluc
```
