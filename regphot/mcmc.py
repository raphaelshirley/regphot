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

def prior(mean,cov):
    """
    Returns the prior probabilty over the parameter vector w. This can be built 
    on a multivariate Gaussian over the parameter space for a given number of 
    Sersic profiles
    
    inputs:
        w,      1-d array,      parameter vector
        cov,    array,          covariance matrix for parameter
        
    outputs:
        prior,  scipy object,   normalised probability density for gaussian 
                                prior with specified covariance matrix the
                                output object can generate a pdf for a given 
                                parameter vector by calling "prior.pdf(w)"
    """
    prior  = multivariate_normal(mean,cov)
    
    return prior
    
    
    
def lnlike(image,w):
    """
    Returns the log liklihood for a given data point (image) and parameter 
    vector w.
    
    $ln(\mathcal{L}) = \chi^2 $
    
    
    """
    chi2 = chi2('pyprofit',w)
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


