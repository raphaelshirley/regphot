#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 15 16:41:07 2016

The package REGPHOT is designed for the analysis of resolved galaxy images.

It can run chi squared minimisation using Galfit or pyprofit on single band images as well
as GalfitM on multiband data.

It can run these in batch on a catalogue of objects either by cutting out postage stamps 
from a large image or downloading them from SDSS. Or from a folder of already cutout postage 
stamps.

I am currently working on developing tools for Bayesian investigation of galaxy profile
parameter space. Using MultiNest i hope to compute the Bayesian evidence for competing 
models.

e.g. one Sersic vs two Sersic, Sersic vs Gaussian vs two Gaussians.

Example:
	examples/ contains some Jupyter notebooks investigating the relevant problems to 
	building this code including an investigation of the Simard Sersic and B/D fits.
	
Todo:
	* Implement priors and liklihood for simple case of uniform or seperable normal.
	* Run tests on NGC 450
	* Run tests on Simard B/D vs Sersic F-test comparisons

Dr Raphael Shirley
Daphne Jackson Research Fellow
Astronomy Centre
University of Sussex

www.raphaelshirley.co.uk
@raphaelshirley

@author: rs548
"""

