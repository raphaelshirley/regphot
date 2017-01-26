#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 11 13:20:20 2017

@author: rs548
"""

from astroquery.sdss import SDSS
from astropy import coordinates
from astropy.io import fits
import os

pos = coordinates.SkyCoord('0h8m05.63s +14d50m23.3s', frame='icrs')
xid = SDSS.query_region(pos, spectro=True)
sp = SDSS.get_spectra(matches=xid)
im = SDSS.get_images(matches=xid, band='g')
plate = xid[0]['plate']

print(plate)

for f in os.listdir('/Users/rs548/Documents/Science/PeteHurley/SDSS/'):
    print(f)

def getPlateFits(position,bandName):
    #First check if it is already there by looping through files
    pos = coordinates.SkyCoord(position, frame='icrs')
    xid = SDSS.query_region(pos)
#    bandName = 'g'
#    plate = xid[0]['plate']
    images = SDSS.get_images(matches=xid, band=bandName)
#    filename = str(plate) + '-' + bandName + '.fits'
#    alreadyDownloaded = False
#    for f in os.listdir('/Users/rs548/Documents/Science/PeteHurley/SDSS/'):
#        print(f)
#        if f == filename:
#            alreadyDownloaded = True
#    
#    if alreadyDownloaded == False:
#        im = SDSS.get_images(matches=xid, band=bandName)
#        im.writeto('/Users/rs548/Documents/Science/PeteHurley/SDSS/' + filename)
#    else:
#        im = fits.open('/Users/rs548/Documents/Science/PeteHurley/SDSS/' + filename)
        
    return images
