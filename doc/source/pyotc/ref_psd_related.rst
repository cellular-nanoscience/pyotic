.. 
   .. sectnum::

==========
PSD module
==========

.. topic:: Manage power spectral densities measurements.

.. contents::
   :depth: 4


PSD class
---------

.. autoclass:: psd.PSD
   :members: __init__, get_freq, get_psd, get_err, get_f_sample,
             set_f_sample, get_temp, set_temp, get_temp_err,
             set_temp_err, get_radius, set_radius, get_radius_err,
             set_radius_err, get_height, set_height, df, is_lateral,
             reset_mask, add_mask, exclude_freq, exclude_freq_outside,
             get_density, get_viscosity

PSDMeasurement class
--------------------

.. autoclass:: psd.PSDMeasurement
   :members: get_names, get_laterality, set_f_sample, set_temp,
             set_temp_err, set_radius, set_radius_err, set_height,
             get_freq, get_psd, get_psd_err, get_stokes_drag,
             get_stokes_drag_err, add_psd, load, set_ac_params, save,
             get_ex_freq, get_ex_amplitude, get_ex_amplitude_err,
             get_ex_power, get_ex_power_err, reset_masks,
             exclude_freq, exclude_freq_outside, plot_psds


PSD related functions
---------------------

These function are directly related to a power spectrum.

Generate a power spectrum
^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: psd
   :members: gen_PSD_from_time_series

Lorentzian-shaped power spectrum
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: psd
   :members: lorentzian_psd

Hydrodynamically correct power spectrum
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: psd
   :members: hydro_psd

PSD modification
----------------

Low-pass filter
^^^^^^^^^^^^^^^

You can either generate a low-pass filtered PSD by multiplying the
output of the ``low_pass_filter`` function to the power spectral
density values or use the decorator ``apply_low_pass_filter`` to
generate a function with the respective input parameters.

.. automodule:: psd
   :members: low_pass_filter, apply_low_pass_filter

Aliased power spectrum
^^^^^^^^^^^^^^^^^^^^^^

The aliased version of a power spectral density function can be
created via this decorator.

.. automodule:: psd
   :members: apply_aliasing

Example
"""""""

Some example code that calculates aliasing of a Lorentzian. ::

 from pyotc.psd import lorentzian_psd, apply_aliasing
 from pyotc import add_plot_to_figure

 from scipy import arange

 fs = 30000 # Hz
 f = arange(1, fs / 2, 1)
 D = 1e-12 # in µm²/s
 f_c = 1500 # in Hz

 aliased_lorentzian = apply_aliasing(lorentzian_psd, fs)

 fig = plt.figure()
 ax = add_plot_to_figure(fig, f, lorentzian_psd(f, D, f_c))
 add_plot_to_figure(fig, f, aliased_lorentzian(f, D, f_c),
                    axis=ax, logplot=True,
                   title='Aliasing a PSD',
                   xlabel='Frequency (Hz)',
                   ylabel='PSD ($\mathrm{V^2/Hz}$)')

.. image:: ./pix/aliasing_a_psd.png
   :align: center
   :scale: 30%

