from setuptools import setup, find_packages

import sys
if sys.version_info >= (3,0):
    sys.exit('Sorry, Python >= 3.0 is not supported')

setup(name='flaapluc',
      version='1.2.3',
      description='Quick generation of alerts from Fermi-LAT data',
      url='',
      author='Jean-Philippe Lenain',
      author_email='jlenain@in2p3.fr',
      license='BSD',
      packages=find_packages(),
      install_requires=[
          'astropy',
          'matplotlib',
          'ephem'
      ],
      scripts=[
          'bin/flaapluc',
          'bin/flaapluc-allsources'
      ],
      zip_safe=False
)
