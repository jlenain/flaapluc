# FLaapLUC

Welcome to the FLaapLUC (Fermi/LAT automatic aperture photometry Light C<->Urve) pipeline !

This pipeline can be used to quickly generate short- or long-term Fermi/LAT light curves. It relies on the Fermi/LAT aperture photometric analysis (cf. e.g. http://fermi.gsfc.nasa.gov/ssc/data/analysis/scitools/aperture_photometry.html), which basically assumes that the data set, within 1 degree around the source, is background-free.

Thus, it is meant to have a rough feeling about the variability of a source, but is absolutely not reliable enough to obtain publication-ready data.

The pipeline can be used to assess the variability of a particular source, using the script `automaticLightCurve.py`, or to process all the sources included in a "master" list of sources `listSources.txt` using the script `processAllSources.py`. Basically, `processAllSources.py` just calls an instance of the `autoLC` class in `automaticLightCurve.py` for each entry in the list of sources.


## Principle

The pipeline generates light curves from Fermi/LAT data, using an unusually small "region of interest" of 1 degree, in which the data set is assumed to be background-free. No likelihood fit is performed with the aperture photometry method. If you want a reliable result, please consider performing a standard Fermi/LAT data analysis (cf. http://fermi.gsfc.nasa.gov/ssc/data/analysis/scitools/).

For a given source, the pipeline first looks whether the source of interest has a counterpart in the 3FGL catalogue. This is done by scanning the "master" list of sources `listSources.txt`, which includes a column with the 3FGL name, if any, for a given source.

- If the source is included in the 3FGL catalogue, the scripts use the catalogue FITS file `gll_psc_v16.fit` to generate an XML model file using the script `make3FGLxml.py`. [This script](http://fermi.gsfc.nasa.gov/ssc/data/analysis/user/) was contributed to the *Fermi* Science Tools by T. Johnson. It means that we assume that the spectral parameters of our source of interest is as per the 3FGL catalogue. This is of course not valid for a source with a strong spectral variability.
- If the source is not included in the 3FGL catalogue, the pipeline assumes that the spectrum of the source is a power-law with a photon index Gamma=-2.5. This value is encoded in the file `automaticLightCurve.py` in the function `processSrc`. Feel free to change this value.

On a local machine, the script `processAllSources.py` can call the shell command `parallel`, to process different sources in parallel (1 thread per source), in an Unix "nice" way. It means that if other processes are running on the machine, the pipeline will not overload the machine.

If the last flux point of a source is above the trigger threshold, a mail will be generated and sent, with a PNG figure of the light curve in attachement.

The trigger threshold is set to 1.e-6 ph cm^-2 s^-1 by default. However, this threshold can be individually set for each source in the file `listSources.txt`. Moreover, a long term data set can be computed, and a dynamical trigger threshold can be used to rely on long term data to assess if a source is active or not. Look at the configuration files to set your custom sigma level above the long term averaged light curve to generate a trigger.


## Download

The pipeline is hosted in a git public repository on GitHub.

To download the pipeline, type in a terminal:

```
git clone https://github.com/jlenain/flaapluc.git
```

## Implementation

The pipeline includes:

- automaticLightCurve.py: to process a single source.
- processAllSources.py: to process a bunch of sources, which actually calls the script `automaticLightCurve.py` for each individual source.
- listSources.txt: an example of a list of sources.
- dummy.cfg: a configuration file where the user enters the configuration to make the pipeline work, such as file paths for inputs and outputs. You can of course have use several configuration files in one site.
- README.md: this file.

Help can be found on the two individual scripts, by typing:

```sh
./automaticLightCurve.py -h
```

```sh
./processAllSources.py -h
```

The pipeline depends on the *Fermi* Science Tools, and it is assumed that you have a local installation of this software on your machine (for more details, see http://fermi.gsfc.nasa.gov/ssc/data/analysis/software/). In more details, the two Python scripts 'automaticLightCurve.py' and 'processAllSources.py' should be able to find the `gt_apps` Python modules provided in the *Fermi* Science Tools (be careful to have your shell variable `$PYTHONPATH` up-to-date !).
