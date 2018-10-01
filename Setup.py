from setuptools import setup

setup(name='RedPipe',
      version='1.0.0a1',
      description='Optical Photometric & Spectroscopic Reduction Pipeline',
      long_description=open('README.rst').read(),
      classifiers=['Development Status :: 3 - Alpha',
                   'Intended Audience :: Science/Research',
                   'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
                   'Programming Language :: Python :: 2.7'],
      url='https://github.com/sPaMFouR',
      author='Avinash Singh',
      author_email='avinash.singh@iiap.res.in',
      install_requires=['numpy', 'scipy', 'pyraf', 'astropy', 'matplotlib', 'ephem', 'imreg_dft', 'specutils'],
      packages=['Photometry', 'Spectroscopy', 'ObsPlan', 'PhotPreProcess', 'FluxCalib',
                'SpecPreProcess', 'AlignFrames', 'CalcAirmass'],
      include_package_data=True,
      zip_safe=False)
