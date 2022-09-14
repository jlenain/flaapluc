from setuptools import setup, find_packages

setup(name='flaapluc',
      version='2.0.0',
      description='Quick generation of alerts from Fermi-LAT data',
      url='',
      author='Jean-Philippe Lenain',
      author_email='jlenain@in2p3.fr',
      license='BSD',
      packages=find_packages(),
      install_requires=[
          'astropy',
          'matplotlib',
          'ephem',
          'urllib3'
      ],
      scripts=[
          'bin/flaapluc',
          'bin/flaapluc-allsources'
      ],
      zip_safe=False
)
