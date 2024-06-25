#!/usr/bin/env python
from setuptools import setup

version = re.findall(r"__version__ = \"(.*?)\"", open(fname).read())[0]

setup(name='RedPipe',
      version='0.1.1',
      description='Optical Photometric & Spectroscopic Reduction Pipeline',
      long_description=open('README.rst').read(),
      classifiers=['Development Status :: 0 - Alpha',
                   'Intended Audience :: Science/Research',
                   'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
                   'Programming Language :: Python :: 3.9'
                   'Topic :: Scientific/Engineering :: Astronomy'],
      url='https://github.com/sPaMFouR/RedPipe',
      author='Avinash Singh',
      author_email='avinash21292@gmail.com',
      install_requires=['numpy', 'scipy', 'pyraf', 'astropy', 'matplotlib', 'ephem', 'imreg_dft', 'specutils'],
      packages=['photometry', 'spectrosscopy', 'filters', 'plotting', 'analysis',
                'preprocess'],
      include_package_data=True,
      zip_safe=False)ß