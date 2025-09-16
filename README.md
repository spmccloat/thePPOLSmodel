# The PPOLs Model
A planet formation model using the "pebble snow" mechanism. Create a disk and arrangement of planetesimal seed masses to grow protoplanets during the gas phase of the protoplanetry disk, track the location of an evolving snow line, and account for the accretion of dry vs. water-rich pebbles.

The name credits the two models that are combined and built upon: pebble-predictor (PP) by Drazkowska et al. (2021) and "epsilon" by Ormel & Liu (OL). See below for links, credits, and license information.

An introductory paper is McCloat, Mulders, & Fieber-Beyer (2025).

HOW TO USE THE CODE

Users are best served by [reading the docs](https://spmccloat.github.io/thePPOLSmodel/), along with the jupyter notebook Tutorial.ipynb.

LICENSE

This code is licensed under the GNU General Public License v3. Among other things, this means that if you modified this code for a publication, the original code needs to be cited (see below), the changes need to be stated, and the new code needs to be open sourced under the same license.

MODIFICATIONS TO ORIGINAL CODE

As per the license, the code presented in "The PPOLs Model" uses pebble-predictor [(Drazkowska et al., 2021)](https://www.aanda.org/articles/aa/abs/2021/03/aa39925-20/aa39925-20.html), [![DOI](https://zenodo.org/badge/300679267.svg)](https://zenodo.org/badge/latestdoi/300679267). The original code and calculations of pebble-predictor remain intact, with the following modifications:

- function now accepts parameter for snow line location, which defaults to no snow line
- disk dust surface density (SigmaDust) is reduced by 1/2 inside snow line when present
- returns parameter for Toomre stability criterion (Q), calculated from values of the disk model

Users interested in a disk model that calculates flux averaged Stokes number and pebble flux alone may be better served by using the original pebble-predictor; these modifications are intended for use within the PPOLs Model.

The PPOLs Mdel uses the pebble accretion efficiency code "epsilon" from [Liu & Ormel (2018)](https://ui.adsabs.harvard.edu/abs/2018A%26A...615A.138L/abstract), and [Ormel & Liu (2018)](https://ui.adsabs.harvard.edu/abs/2018A%26A...615A.178O/abstract). Their code has not been modified from its original version (except for documentation). You can find the code on [Chris Ormel's website](https://staff.fnwi.uva.nl/c.w.ormel/software.html)

CREDITS

If you use this code for your work, please cite the corresponding paper, this github repository and/or the corresponding Zenodo DOI.