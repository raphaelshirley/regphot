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

import subprocess
from os.path import dirname


def git_version():
    """Returns the git version of the module
        This function returns a string composed of the abbreviated Git hash of the
        module source, followed by the date of the last commit.  If the source has
        some local modifications, “ [with local modifications]” is added to the
        string.
        This is used to print the exact version of the source code that was used
        inside a Jupiter notebook.
        """
    module_dir = dirname(__file__)
    
    command_hash = "cd {} && git rev-list --max-count=1 " \
        "--abbrev-commit HEAD".format(module_dir)
    command_date = "cd {} && git log -1 --format=%cd" \
        .format(module_dir)
    command_modif = "cd {} && git diff-index --name-only HEAD" \
        .format(module_dir)

    try:
        commit_hash = subprocess.check_output(command_hash, shell=True)\
            .decode('ascii').strip()
        commit_date = subprocess.check_output(command_date, shell=True)\
            .decode('ascii').strip()
        commit_modif = subprocess.check_output(command_modif, shell=True)\
            .decode('ascii').strip()
        
        version = "{} ({})".format(commit_hash, commit_date)
        if commit_modif:
            version += " [with local modifications]"
    except subprocess.CalledProcessError:
        version = "Unable to determine version."
    
    return version

