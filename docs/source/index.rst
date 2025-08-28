Welcome to The PPOLs Model's documentation!
===========================================
Intro
-----
This Python package is a planet formation model that explores the "pebble snow" mechanism. Users can easily set and run planet formation models for a variety of disk models. Major features include:

* a disk model with explicit stellar mass and disk mass
   * disk mass can be fraction of stellar mass, or have the solid disk mass set exactly
* dust mass that depletes as it converts into pebbles
* temperature profile that is a simple power law or account for viscous/irradiation heating (e.g. Ida et al. 2016)
* a snow line (water-ice) based on temperature, set explicitly, or self-consistent with disk parameters
   * the snow line that can evole (or not) based on disk parameters
* any number of planetesimal seed masses, with any initial mass, introduced at any time
* planetesimal seed masses that grow simultaneously via pebble accretion and remove ("filter") pebble mass as it drifts inward
* accounting for dry vs water-rich solid mass that is accreted with time

The model was built to test outcomes of pebble accretion, specifically during first ~10 million years of the protoplanetary disk when gas is still present. It may serve well as inputs for longer n-body simulations.

.. include:: intro
.. toctree::
   :maxdepth: 2
   :caption: Contents:

   intro
   thePPOLsCode

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`