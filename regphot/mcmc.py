#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 13 16:17:59 2017

As a first pass I wanty to try and fit one sersic profile to an image with 
stars and another galaxy using chisquared from galfit

params to fit in a sersic:
    RA                  (0 to image size)
    Dec                 (0 to image size)
    Magnitude           (depth to max in sky?)
    half light radius   (psf to max beyond local group?)
    sersic index        (1 to 4?)
    axis ratio          (1 to 10?)
    position angle      (0 to 180)
    
In a dream world I'll get two peaks one for each galaxy in such a way that
I can use it to suggest trying two sersic profiles

The sersic profiles then then be sued as weights for assigning flux - note that
this would be implicitly assuming that morphology was small in comparison to
surface brightness. As in variation across surface due to eg spiral arms was 
lower than total brightness. This may well be wrong meaning that we are 'using
the data twice'.

@author: rs548
"""

#import emcee
#import galfit
import numpy as np
from scipy.stats import multivariate_normal
import pymultinest
import math, os
if not os.path.exists("chains"): os.mkdir("chains")


#I may want to define a seperate model class that can produce a general model
def model1():
    """
    This model is composed of one Sersic profile
    """
    
def model2():
    """
    This model is composed of two Sersic profiles
    """

def prior(cube, ndim, nparams):
    """
    Returns the prior probabilty by stretching the unit cube. This can be built 
    on a multivariate Gaussian over the parameter space for a given number of 
    Sersic profiles. The required form for multinest takes the cube of 
    probabilities and returns the inverse cumulative distribution. The value 
    that has a probability less than the input value p such that  P(x<X) = p
    
    inputs:
        cube,       1-d array,          input a uniform prior on unit hypercube
        ndim,       int,                covariance matrix for parameter
        nparams,    int,                number of parameters in model
        
    outputs (this is not returned by the function but modified by it):
        cube,       1-d array,          the unit hypercube is mapped to the 
                                        parameter hypercube
    """
#    for i in range(nparams):
#        cube[i] *=  scipy.stats.norm(mu, std).ppf(cube[i])
    cube[0] = 
    cube[1] = 
    
    return prior
    
    
    
def lnlike(cube, ndim, nparams):
    """
    Returns the log liklihood for a given data point (image) and parameter 
    vector w.
    
    $ln(\mathcal{L}) = -1/2 \chi^2 + const $
    
    
    """
    #First the cube must be transformed into the parameter vector
    w = [None] * ndim
    for i in range(ndim):
        w[i] = 
    chi2 = chi2('pyprofit',cube)
    lnlike = -0.5* chi2
    return lnlike

def multinestevidence(mean,cov,model):
    """
    Input prior over paramters and black box function for evaluating ln liklihood
    based on chi2 from image/model and calcualte evidence for given model
    in preferebly the exact same way to emcee to allow comparison
    """
    prior1 = prior(mean,cov)
    pymultinest.run(lnlike, prior1.pdf, n_params, 
                    resume = True, verbose = True)
    
def emceeevidence():
    """
    Input prior over paramters and black box function for evaluating ln liklihood
    based on chi2 from image/model and calcualte evidence for given model
    in preferebly the exact same way to multinest to allow comparison
    """


