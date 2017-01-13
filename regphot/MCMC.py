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

def fakeTestChi2(ra,dec,mag,re,n,axis,angle,image):
    print(mag)
    x = np.array((ra,dec,mag,re,n,axis,angle))
    ivar = np.array((10,10,10,10,1,0.5,10))
    #the following gives the ln of the gaussian
    return -0.5 * np.sum(ivar * x ** 2)