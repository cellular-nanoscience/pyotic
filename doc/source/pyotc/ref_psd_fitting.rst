==================
PSD fitting module
==================

.. topic:: Manage the least squares fitting of power spectra.

.. contents::
   :depth: 4

PSDFit class
------------

.. autoclass:: psd_fitting.PSDFit
   :members: __init__, get_freq, get_psd, get_psd_err, set_bounds, exclude_freq, analytical_lorentzian_fit, plot_anal_fits, collective_psd_fit, setup_fit, fit_psd, fit_psds, plot_fits, get_basepower, calc_pc_results, calc_ac_results, print_results

Functions related to fitting
----------------------------

Create fitting parameters
^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: psd_fitting
   :members: gen_psd_fit_pars


Create a fitting model
^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: psd_fitting
   :members: gen_model_fun


Fit to a model to a psd
^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: psd_fitting
   :members: fit_psd


Fit a model collectively
^^^^^^^^^^^^^^^^^^^^^^^^

Call function
"""""""""""""
Fit a model to a set of power spectra simultaneously. This can be useful to find the parameters of a QPDs low-pass filter characteristics. The idea is, to use a set of different power spectra (e.g. with different beads at a fixed trapping-laser power) and collectively find the best fit value for the cut-off frequency and the low-pass filtering efficiency $\alpha$.

.. automodule:: psd_fitting
   :members: make_collective_psd_fit

Residual function
"""""""""""""""""

.. automodule:: psd_fitting
   :members: collective_psd_fit_fun


Outlier determination
^^^^^^^^^^^^^^^^^^^^^

.. automodule:: psd_fitting
   :members: is_outlier


Analytically solve the least-squares problem of for a Lorentzian
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: psd_fitting
   :members: analytical_lsq_lorentzian

Find optimal boundaries
^^^^^^^^^^^^^^^^^^^^^^^

Try to find better boundaries (f_min, f_max) for the analytical least-squares solution of a Lorentzian-shaped power spectrum.

.. automodule:: psd_fitting
   :members: calc_anal_lsq_opt


