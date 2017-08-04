.. contents::

Calibration?
------------

Let's say, you do an experiment with an optical tweezers. For example,
a molecular motor, e.g. kinesin-1, is bound to a microsphere (=bead),
which is made of polystyrene with a certain diameter. You trap bead,
move it to a microtuble and it happens - because you are lucky - that
the bead consumes the present ATP and starts walking along the
microtuble. Because the bead is still held by the optical trap, the
kinesin needs to pull the bead out of center of the trap and each step
will lead to a little higher force that the trap excerts to the
bead. You measure the displament of the bead interferometrically with
a position sensitive device like a quadrant photodiode (QPD). The
signal that you record will be given in Volts and will depend on many
different physical parameters that are defined by the optical tweezers
setup itself, the bead, the sample, *etc.* However, converting the
signal into a physical quantity, such as the displacement of the bead
in nanometers or the force that pulls the bead out of the center of
the trap, is called calibration.

There are different approaches to calibrate an optical trap. An
overview is given in [Jun et al.]

.. [Jun et al.] 1. Jun, Y., Tripathy, S. K., Narayanareddy, B. R. J.,
                Mattson-Hoss, M. K. & Gross, S. P. Calibration of
                Optical Tweezers for In Vivo Force Measurements: How
                do Different Approaches Compare? Biophys. J. 107,
                1474–1484 (2014).

Common approaches calibrate by:

   +----------------------------------------+----------------------------+----------------------------+-------------+
   | Method                                 | Prior known parameters     | Calibrating parameter      | Ref.        |
   +========================================+============================+============================+=============+
   | equipartition theorem                  | temperature,               | trap stiffness             |             |
   |                                        | displacement sensitivity   |                            |             |
   +----------------------------------------+----------------------------+----------------------------+-------------+
   | scanning through stuck bead            | --                         | displacement sensitivity   |             |
   +----------------------------------------+----------------------------+----------------------------+-------------+
   | Boltzmann statistics                   | temperature,               | trap stiffness             |             |
   |                                        | displacement sensitivity   |                            |             |
   +----------------------------------------+----------------------------+----------------------------+-------------+
   | viscous drag                           | flow speed of medium       | displacement sensitivity   |             |
   +----------------------------------------+----------------------------+----------------------------+-------------+
   | passive power spectral density         | temperature,               | trap stiffness             | [#ber2004]_ |
   |                                        | displacement sensitivity   |                            |             |
   |                                        | (via radius)               |                            |             |
   +----------------------------------------+----------------------------+----------------------------+-------------+
   | active power spectral density          | temperature                | trap stiffness             | [#tol2006]_ |
   |                                        |                            | displacement sensitivity   |             |
   +----------------------------------------+----------------------------+----------------------------+-------------+

.. [#ber2004] K. Berg-Sørensen & H. Flyvbjerg 2004 - `doi
              <http://dx.doi.org/10.1063/1.1645654>`_
.. [#tol2006] S. F. Tolić-Nørrelykke, et al. 2006 - `doi
              <ttp://scitation.aip.org/content/aip/journal/rsi/77/10/10.1063/1.2356852>`_


Passive PSD analysis
--------------------

In passive PSD analysis, the brownian motion of a trapped particle is
monitored and its power spectral density is used to calibrate for the
apperent trap stiffness.  This is done by fitting a model to the PSD.

.. figure:: ./pix/PSD_passive.png
   :scale: 70 %
   :alt: power spectral density
   :align: center

   Power spectral density of a trapped microsphere.


Active PSD analysis
-------------------

In active calibration you simultaneously calibrate for both, the trap
stiffness and the displacement senitivity. While the former is done
the same way as it is done in passive PSD analysis, the latter
directly measures the displacement sensitivity. This is done by
exciting the bead sinusoidially with a known amplitude and frequency
and measuring the displacement at that frequency and comparing it to
the theoretical value.

.. figure:: ./pix/PSD_active.png
   :scale: 70 %
   :alt: excited power spectral density
   :align: center

   Power spectral density of a trapped microsphere with an additional
   sinusoidal movement at 32 Hz.


Height dependent calibration
----------------------------

A height-dependent calibration uses the height-dependence of the
viscous drag close to a plane surface, which is well described by
`Faxen's law <https://en.wikipedia.org/wiki/Fax%C3%A9n%27s_law>`_.
This can be used to determine:

 - the exact height of the sphere
 - the size of the sphere **or**
 - the viscosity of the medium
 
.. figure:: ./pix/drag_fit.png
   :scale: 70 %
   :alt: Faxén's law fit to relative drag.
   :align: center

   Fit of the relative drag to Faxén's law.
