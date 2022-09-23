[![ascl](https://img.shields.io/badge/ascl-1709.011-blue.svg?colorB=262255)](http://ascl.net/1709.011)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.906991.svg)](https://doi.org/10.5281/zenodo.906991)
[![Anaconda-Server Badge](https://anaconda.org/jlenain/flaapluc/badges/version.svg)](https://anaconda.org/jlenain/flaapluc)

# FLaapLUC

Welcome to the FLaapLUC (Fermi-LAT automatic aperture photometry Light C<->Urve) pipeline !

This pipeline can be used to quickly generate short- or long-term Fermi-LAT light curves. It relies on the Fermi-LAT aperture photometric analysis (cf. e.g. http://fermi.gsfc.nasa.gov/ssc/data/analysis/scitools/aperture_photometry.html), which basically assumes that the data set, within 1 degree around the source, is background-free.

Thus, it is meant to have a rough feeling about the variability of a source, but is absolutely not reliable enough to obtain publication-ready results.

The pipeline can be used to assess the variability of a particular source, using the script `flaapluc`, or to process all the sources included in a "master" list of sources `listSources.txt` using the script `flaapluc-allsources`.


## Principle

The pipeline generates light curves from Fermi-LAT data, using an unusually small "region of interest" of 1 degree, in which the data set is assumed to be background-free. No likelihood fit is performed with the aperture photometry method. If you want a reliable result, please consider performing a standard Fermi-LAT data analysis (cf. http://fermi.gsfc.nasa.gov/ssc/data/analysis/scitools/).

For a given source, the pipeline first looks whether the source of interest has a counterpart in the 4FGL catalogue. This is done by scanning the "master" list of sources `listSources.txt`, which includes a column with the 4FGL name, if any, for a given source.

- If the source is included in the 4FGL catalogue, the scripts use the catalogue FITS file as given in the input configuration file to generate an XML model file using the script `make4FGLxml.py`. [This script](http://fermi.gsfc.nasa.gov/ssc/data/analysis/user/) was contributed to the *Fermi* Science Tools by T. Johnson. It means that we assume that the spectral parameters of our source of interest is as per the 4FGL catalogue. This is of course not valid for a source with a strong spectral variability.
- If the source is not included in the 4FGL catalogue, the pipeline assumes that the spectrum of the source is a power-law with a photon index Gamma=-2.5. This value is encoded in the file `automaticLightCurve.py` in the function `processSrc`. Feel free to change this value.

On a local machine, the script `flaapluc-allsources` can call the shell command `parallel`, to process different sources in parallel (1 thread per source), in an Unix "nice" way. It means that if other processes are running on the machine, the pipeline will not overload the machine.

If the last flux point of a source is above the trigger threshold, a mail will be generated and sent, with a PNG figure of the light curve in attachement.

The trigger threshold is set to 1.e-6 ph cm^-2 s^-1 by default. However, this threshold can be individually set for each source in the file `listSources.txt`. Moreover, a long term data set can be computed, and a dynamical trigger threshold can be used to rely on long term data to assess if a source is active or not. Look at the configuration files to set your custom sigma level above the long term averaged light curve to generate a trigger.


## Installation and configuration

Information on installation and setup can be found in `INSTALL.md`.


## Implementation

The pipeline includes:

- `flaapluc`: to process a single source.
- `flaapluc-allsources`: to process a bunch of sources, which actually calls the class `automaticLightCurve` for each individual source.
- `listSources.txt`: an example of a list of sources.
- `dummy.cfg`: a configuration file where the user enters the configuration to make the pipeline work, such as file paths for inputs and outputs. You can of course have use several configuration files in one site.
- `README.md`: this file.

Help can be found on the two individual scripts, by typing:

```sh
flaapluc -h
```

```sh
flaapluc-allsources -h
```

The pipeline uses an allsky file as input, which could easily be created using [`enrico`](http://enrico.readthedocs.org/en/latest/index.html). `FLaapLUC` is actually using `enrico` and can use it to generate an allsky file on the fly for you.

`FLaapLUC` depends on the *Fermi* Science Tools. The installation process should take care of this dependency, but see https://github.com/fermi-lat/Fermitools-conda/wiki for more details. More specifically, the `automaticLightCurve` class should be able to find the `gt_apps` Python module provided in the *Fermi* Science Tools.
