#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 13 12:23:18 2016

@author: rs548
"""


from __future__ import division, print_function

from astropy.io import fits
from astropy.wcs import WCS
from astropy.nddata import Cutout2D
#from astropy import units as u
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

#Open FITS catalogue
hdulist = fits.open('cosmos-hyperleda-sdss-dr13-2massXSC-4reg-flagged_HELP.fits')

#Open all the huge FITS images
imageY = fits.open('./data/ADP.2016-03-17T08:31:54.127.fits')
imageJ = fits.open('./data/ADP.2016-03-17T08:31:54.147.fits')
imageH = fits.open('./data/ADP.2016-03-17T08:31:54.060.fits')
imageKs = fits.open('./data/ADP.2016-03-17T08:31:54.113.fits')
imageNB118 =fits.open('./data/ADP.2016-03-17T08:31:54.107.fits')


w = WCS('./data/ADP.2016-03-17T08:31:54.127.fits')
nisin = 0
nnotin = 0 
for source in hdulist[1].data:
    #print(source[0])
    xpix, ypix = w.all_world2pix(source[1],source[2],0)
    #top left is 150.77177 +02.80870
    #bottom right is 149.29941 +01.59448
    dimension = 150
    size = (dimension,dimension)
    position = (xpix , ypix )
    if (source[1] <= 150.77177 
        and source[1] >= 149.29941 
        and source[2] <= 2.80870 
        and source[2] >= 1.59448):
        print(source[0], 'is in COSMOS field')
        nisin = nisin + 1
        cutout = Cutout2D(imageY[0].data, position, size)
        fig = plt.figure()
        plt.imshow(cutout.data, cmap='gray', norm=LogNorm(), interpolation='none')
        fig.savefig('./output/' + source[0] + '.png')
    else:
        print(source[0], 'is NOT in COSMOS field')
        nnotin = nnotin + 1
        
        
print('there are ', nisin, 'objects in the field, and ', nnotin, ' out of the field.')