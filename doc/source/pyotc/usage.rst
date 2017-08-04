=================================
 What is PyOTC and what does it?
=================================

PyOTC is a python optical tweezers calibration software package. Its
aim is to provide means to precisely calibrate a trapped particle by
analyzing the power spectral density of its thermal motion.

In particular, it provides means for a height-dependent calibration of
a trapped bead. This calibration approach is especially interesting
for experiments that are carried out close to the sample-chamber
surface, where conventional PSD analysis is difficult or even
impossible. It is also useful in cases where the size of the bead
**or** the viscosity of the sample medium and/or the actual height of
the bead is unknown.

The package provides the following features:

* calculation of power spectral density from time-dependent signals
* saving and loading of PSD measurments
* active and passive power spectral density analysis
* fitting of PSDs to different models
* consideration of low-pass filtering and aliasing of the signals
* height-dependent analysis of active PSD measurements

Height-dependent analysis of active PSD calibrations enables you to
determine the

* particle size or 
* viscosity of the medium,
* location of the sample surface and
* height-dependence of the trap stiffness and displacement sensitivity
 
PyOTC is written in `Python <www.python.org>`_ 3.4 and was developed
for the use with `Jupyter <http://jupyter.org/>`_.

=================
 Getting started
=================

.. todo::
   
   getting started

===========
Calibration
===========

.. toctree::
   :maxdepth: 1

   ./psd_calibration.rst

========
Examples
========

.. todo::
   add examples

=========
TODO list
=========

.. todolist::
