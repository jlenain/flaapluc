# Installation

## Using conda

It is recommended, and by far the easiest way, to install FLaapLUC using e.g. [conda](https://www.continuum.io/downloads):

```
conda install -c jlenain flaapluc
```

If you want to create a new environment for FLaapLUC, the following will do so, handling all required dependencies:

```
conda create -n flaapluc -c jlenain flaapluc
```

## Build from source

It is of course still possible to build FLaapLUC from the source code:

```
git clone https://github.com/jlenain/flaapluc
cd flaapluc
python setup.py install

```


# Configuration

FLaapLUC needs to use the [Fermi-LAT ScienceTools](https://fermi.gsfc.nasa.gov/ssc/data/analysis/software/) as well as [Enrico](http://enrico.readthedocs.io), which should be pre-installed on your system. The configuration is then the following:

```
. $FERMI_DIR/fermi-init.sh
. $ENRICO_DIR/enrico-init.sh
source activate flaapluc
```
