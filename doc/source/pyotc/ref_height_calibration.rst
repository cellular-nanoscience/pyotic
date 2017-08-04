=========================
Height Calibration module
=========================

.. topic:: Manage height-dependent power spectral density measurements
           and their fits.

.. contents::
   :depth: 4

HeightCalibration class
-----------------------

general methods
^^^^^^^^^^^^^^^

.. autoclass:: height_calibration.HeightCalibration
   :members: get_heights, add_height_offset, set_height_unit,
             get_radius, set_radius, get_radius_er, set_radius_err,
             get_wavelength, set_wavelength, reset_mask, get_mask,
             exclude_height, exclude_by_max_rel_drag,
             exclude_heights_outside, generate_outliers_mask,
             reset_psd_masks


Manage psd fits
^^^^^^^^^^^^^^^

.. automethod:: height_calibration.HeightCalibration.add_psdfit
.. automethod:: height_calibration.HeightCalibration.get_psd_files
.. automethod:: height_calibration.HeightCalibration.gen_psd_fits
.. automethod:: height_calibration.HeightCalibration.setup_psd_fits
.. automethod:: height_calibration.HeightCalibration.fit_psds


Height-dependent data
^^^^^^^^^^^^^^^^^^^^^

.. automethod:: height_calibration.HeightCalibration.get_dissens
.. automethod:: height_calibration.HeightCalibration.set_dissens_unit
.. automethod:: height_calibration.HeightCalibration.get_trap_stiffness
.. automethod:: height_calibration.HeightCalibration.set_trap_stiffness_unit
.. automethod:: height_calibration.HeightCalibration.get_drag
.. automethod:: height_calibration.HeightCalibration.set_drag_unit
.. automethod:: height_calibration.HeightCalibration.get_red_chi2
.. automethod:: height_calibration.HeightCalibration.save_hc_data
.. automethod:: height_calibration.HeightCalibration.load_hc_data


Fit height-dependent data
^^^^^^^^^^^^^^^^^^^^^^^^^

.. automethod:: height_calibration.HeightCalibration.fit_rel_drag
.. automethod:: height_calibration.HeightCalibration.fit_height_data
.. automethod:: height_calibration.HeightCalibration.write_results_to_file

Plotting
^^^^^^^^

.. automethod:: height_calibration.HeightCalibration.plot_drag
.. automethod:: height_calibration.HeightCalibration.plot_rel_drag_fit
.. automethod:: height_calibration.HeightCalibration.plot_rel_draf_fit_result
.. automethod:: height_calibration.HeightCalibration.plot_dissens
.. automethod:: height_calibration.HeightCalibration.plot_dissens_fit
.. automethod:: height_calibration.HeightCalibration.plot_trap_stiffness
.. automethod:: height_calibration.HeightCalibration.plot_trap_stiffness_fit
.. automethod:: height_calibration.HeightCalibration.plot_redchi2
.. automethod:: height_calibration.HeightCalibration.plot_pc_results
.. automethod:: height_calibration.HeightCalibration.plot_results


Functions
---------

.. autofunction:: height_calibration.gen_height_fit_pars
.. autofunction:: height_calibration.fit_height_data
