# Installation

## Using conda

It is recommended, and by far the easiest way, to install FLaapLUC using e.g. [conda](https://www.continuum.io/downloads):

```
conda install -c jlenain -c conda-forge -c fermi flaapluc
```

If you want to create a new environment for FLaapLUC, the following will do so, handling almost all required dependencies:

```
conda create -n flaapluc -c jlenain -c conda-forge -c fermi flaapluc
```

FLaapLUC depends on the [Fermi-LAT ScienceTools](https://fermi.gsfc.nasa.gov/ssc/data/analysis/software/) as well as [Enrico](http://enrico.readthedocs.io). The [Fermi-LAT ScienceTools](https://fermi.gsfc.nasa.gov/ssc/data/analysis/software/) should automatically come with a `conda` installation of FLaapLUC. [Enrico](http://enrico.readthedocs.io) should be manually installed, though.

## Build from source

It is of course still possible to build FLaapLUC from the source code to enjoy the bleeding-edge developments (and bugs!):

```
git clone https://github.com/jlenain/flaapluc
cd flaapluc
python setup.py install
```


# Configuration

To activate your configuration for FLaapLUC, just do the following:

```
conda activate flaapluc
```
