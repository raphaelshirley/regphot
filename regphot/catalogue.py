#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 13 12:23:18 2016

@author: rs548
"""

from astropy.io import fits

#Open FITS catalogue
hdulist = fits.open('cosmos-hyperleda-sdss-dr13-2massXSC-4reg-flagged_HELP.fits')

#Open all the huge FITS images
imageY = fits.open('./data/ADP.2016-03-17T08:31:54.127.fits')
imageJ = fits.open('./data/ADP.2016-03-17T08:31:54.127.fits')
imageH = fits.open('./data/ADP.2016-03-17T08:31:54.127.fits')
imageKs = fits.open('./data/ADP.2016-03-17T08:31:54.127.fits')
imageNB118 =fits.open('./data/ADP.2016-03-17T08:31:54.127.fits')

for source in hdulist[1].data:
    print(source[0])