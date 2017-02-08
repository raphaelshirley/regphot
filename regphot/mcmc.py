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

testModel = [['sersic',
              1791., 
              1222.,
              19.,
              107.,
              0.8,
              0.0000,
              0.0000,
              0.0000,
              0.7,
              -8.5],
             ['sersic',
              1902., 
              1387.,
              22.,
              37.,
              0.,
              0.0000,
              0.0000,
              0.0000,
              0.,
              -53.],
             ['sky',
              -1.651e-03,
              0.000e+00,
              0.000e+00]]

def fakeTestChi2(image,model):
    print('fake test chisquared fucntion doesnt use input image ',image)
    
    trueModel =  [['sersic',
              1790.7958, 
              1222.2855,
              19.4770,
              107.7908,
              0.8121,
              0.0000,
              0.0000,
              0.0000,
              0.7140,
              -8.5677],
             ['sersic',
              1902.9004,
              1387.3171,
              22.2516,
              37.0798,
              0.7359,
              0.0000,
              0.0000,
              0.0000,
              0.4275,
              -53.2329],
             ['sky',
              -1.651e-03,
              0.000e+00,
              0.000e+00]]
    chiSquared = 0
    for componentN in range(1,len(model)):
        component = model[componentN]
        testcomponent = trueModel[componentN]
        for m in range(1,len(component)):
            chiSquared = chiSquared + (component[m] - testcomponent[m])**2.0
        
    #the following gives the ln of the gaussian
    return chiSquared

def fakeTestLogChi2(image,model):
    return np.log(fakeTestChi2(image,model))

d = fakeTestLogChi2('/path/to/fakeimage.fits', testModel)