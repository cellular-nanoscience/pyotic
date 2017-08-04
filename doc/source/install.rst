Prerequisites
=============

PyOTIC was written in `Python 3.4`_ and uses the following packages:

  - `scipy`_ >= 0.14.0
  - `numpy`_ >= 1.8.2
  - `matplotlib`_ >= 1.5.1
  - `ZODB`_ >= 4.3.1
  - `cloudpickle`_ >= 0.1.1
  - `ipywidgets`_ >= 4.0.3
  - `lmfit`_ >= 0.9.3
  - `pint`_ == 0.5.2

.. _Python 3.4: https://www.python.org/download/releases/3.4.5/
.. _scipy: https://www.scipy.org/
.. _numpy: http://www.numpy.org/
.. _matplotlib: http://matplotlib.org/
.. _ZODB: http://www.zodb.org/en/latest/
.. _cloudpickle: https://pypi.python.org/pypi/cloudpickle/
.. _ipywidgets: https://pypi.python.org/pypi/ipywidgets/
.. _lmfit: https://lmfit.github.io/lmfit-py/
.. _pint: http://pint.readthedocs.org


Download and Installation
=========================

Installation via the Python site packages folder
------------------------------------------------

Windows
^^^^^^^

The PyOTIC software package was tested with `WinPython`_ version
3.4.4.2 64bit.  Download and install the corresponding WinPython
distribution from the website.  Copy the *pyotic* package folder to::

    WINPYTHONDIR\python-3.4.4.amd64\Lib\site-packages

Open a "WinPython command prompt" (can be found in the WINPYTHONDIR) and
install additionally required packages::

    pip3 install ZODB

If you want to read in calibration files created by the PyOTC software package,
also install pint:
::

    pip3 install pint

.. _WinPython: http://winpython.github.io/
.. WINPYTHONDIR = WinPython-64bit-3.4.4.2Qt5

pip package manager
^^^^^^^^^^^^^^^^^^^

.. todo::
   add pip installation how to

via git
^^^^^^^
.. todo::
   git clone ##url##

   python3 setup.py install


via system path addition
^^^^^^^^^^^^^^^^^^^^^^^^

Within your python source code add the path to the *path variable* by::
      
   import sys
   sys.path.append('/path/to/pyoti/')
   sys.path.append('/path/to/pyotc/')

      
After the installation, you can use the pyoti package by::
 
   import pyoti


Acknowledgements
================
.. todo::
   add acknowledgements

License
=======
.. todo::
   add license text
