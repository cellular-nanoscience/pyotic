#!/usr/bin/env python3

from setuptools import setup, find_packages, findall
#from distutils.core import setup
import os
import re
import importlib


def find_data_files(ldir):
    """
    Find files in the given directory in the list of directories.

    Returns a list of (directory, list_of_files) tuples as it is
    necessary for the data_files arguments of setuptools.setup()
    function.
    """
    dfiles = []
    for d in ldir:
        dfiles.append((d, [os.path.join(d, fl)
                           for fl in os.listdir(d)
                           if os.path.isfile(os.path.join(d, fl))
                           ]))

    return dfiles


# read the version from VERSION.txt
with open('VERSION.txt') as version_file:
    __version__ = version_file.read().strip()

# write VERSION.txt files in both packages
for pkg in ['pyoti', 'pyotc']:
    ptf = os.path.join(*['.', pkg, 'VERSION.txt'])
    with open(ptf, 'w') as vf:
        vf.write('{}\n'.format(__version__))
    print('created version file: {}'.format(ptf))

# get the required packages from requirements.txt and check if they are already
# installed. Probably obsolete for newer setuptools versions when using required
# packages arguments of setup() function.
with open('requirements.txt', 'r') as rfile:
    requs = [re.split('([\s=<>]+)', l.replace('\n', '').strip(' '))
             for l in rfile.readlines()
             if not(l.startswith('#')) and l != '\n']

abort = False

for requ in requs:
    try:
        importlib.import_module(requ[0])
    except:
        abort = True
        print("Module {} not available.".format(requ[0]))
        if len(requ) > 1:
           print("\tPyOTIC requires version {}".format(''.join(requ[1:]).strip()))

if abort:
    print("Try to install the missing package(s) first.")
        
if not abort:
    long_desc = """
                Time-dependent signal analysis and calibration of optical
                tweezers by power spectral denstity analysis.
                """
    license = 'Apache-2.0'

    setup(name='pyotic',
          version=__version__,
          author='Tobias Jachowski and Steve Simmert',
          author_email='tobias.jachowski@uni-tuebingen.de; '
                       'steve.simmert@uni-tuebingen.de',
          url='https://github.com/cellular-nanoscience/pyotic',
          license=license,
          description="Time-dependent signal analysis and calibration of optical tweezers.",
          long_description = long_desc,
          platforms = ['Windows', 'Linux', 'Mac OS X'],
          classifiers=['Intended Audience :: Science/Research',
                       'Intended Audience :: Education',
                       'Framework :: Jupyter',
                       'License :: OSI Approved :: Apache Software License',
                       'Operating System :: OS Independent',
                       'Programming Language :: Python',
                       'Topic :: Scientific/Engineering :: Physics',
                       'Topic :: Scientific/Engineering :: Visualization', ],
          python_requires='>=3, <4',
          packages=find_packages(),
          package_dir = {'pyotc': 'pyotc',
                         'pyoti': 'pyoti'},
          package_data={'pyotc': ['VERSION.txt'],
                        'pyoti': ['VERSION.txt', '*.cfg', '*.py', '*.pyx',]},
          data_files=find_data_files(['pyoti/etc/',
                                      'pyoti/etc/record/',
                                      ])
          )
