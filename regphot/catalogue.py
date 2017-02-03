#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 13 12:23:18 2016

@author: rs548
"""


from __future__ import division, print_function

from os import chdir
from astropy.io import fits
from astropy.wcs import WCS
from astropy.nddata import Cutout2D
#from astropy import units as u
import matplotlib.pyplot as plt
#from matplotlib.colors import LogNorm
import time
t0 = time.time()


# MY MODULES
#from regphot import getsdss
from getsdss import getPlateFits
from galfit import optimise
from galfitm import optimiseM




#############################INPUTS############################################
workDir = '/Users/rs548/Documents/Science/PeteHurley/'
chdir(workDir + 'TRASH/')

#Open FITS catalogue
hdulist = fits.open(workDir + 
                    'cosmos-hyperleda-sdss-dr13-2massXSC-4reg-flagged_HELP.fits')



#Options
field = 'UltraVISTA'
runGalfit = True
runGalfitM = False
runUltraVISTA = False
runSDSS = True
plotfigs = False
downloadSDSS = False


#Open all the huge FITS images
imageY = fits.open(workDir + 'data/ADP.2016-03-17T08:31:54.127.fits')
imageJ = fits.open(workDir + 'data/ADP.2016-03-17T08:31:54.147.fits')
imageH = fits.open(workDir + 'data/ADP.2016-03-17T08:31:54.060.fits')
imageKs = fits.open(workDir + 'data/ADP.2016-03-17T08:31:54.113.fits')
imageNB118 =fits.open(workDir + 'data/ADP.2016-03-17T08:31:54.107.fits')

images = (('Y', imageY, 1020),
          ('J', imageJ, 1250),
          ('H', imageH, 1650),
          ('Ks', imageKs, 2200),
          ('NB118', imageNB118, 1180))

#Open Y band for coordinate transformations
w = WCS(workDir + 'data/ADP.2016-03-17T08:31:54.127.fits')
#############################INPUTS############################################
nisin = 0
nnotin = 0 
nSDSS = 0

#Go through every object, check it is in field and then do stuff with it
#I deally this will eventually cutout an appropriate fits file and run
#Galfit or Galfitm on it 
for source in hdulist[1].data:
    print('Now running on ' + source[0])
    #print(source[0])
    xpix, ypix = w.all_world2pix(source[1],source[2],0)
    #top left is 150.77177 +02.80870
    #bottom right is 149.29941 +01.59448
    dimension = 150
    size = (dimension,dimension)
    position = (xpix , ypix )
    #Should replace this bit with MOC test
    if (runUltraVISTA
        and source[1] <= 150.77177 
        and source[1] >= 149.29941 
        and source[2] <= 2.80870 
        and source[2] >= 1.59448):
        print(source[0], ' at', source[1],source[2], 'is in ', field)
        nisin = nisin + 1
        cutout = Cutout2D(imageY[0].data, position, size)

        fig = plt.figure()
        #norm=LogNorm(),
        plt.imshow(cutout.data, cmap='gray',  interpolation='none')
        plt.title(source[0])
        fig.savefig(workDir + 'output/' + source[0] + '.png')
        try:
            fits.writeto(workDir + 'galfit/' + source[0] + '.fits', cutout.data)
        except:
            os.remove(workDir + 'galfit/' + source[0] + '.fits')
            fits.writeto(workDir + 'galfit/' + source[0] + '.fits', cutout.data)  

        for band in images:
            try:
                fits.writeto(workDir + 'galfitm/' + source[0] +'-'+ band[0] + '.fits', 
                             Cutout2D(band[1][0].data, position, size).data)
            except:
                os.remove(workDir + 'galfitm/' + source[0] +'-'+  band[0] + '.fits')
                fits.writeto(workDir + 'galfitm/' + source[0] +'-'+  band[0] + '.fits', 
                             Cutout2D(band[1][0].data, position, size).data) 
                

            

        if runGalfit:
            optimise(source[0] + '.fits')
        if runGalfitM:
            optimiseM(images,source[0])
    else:
        #print(source[0], 'is NOT in ', field)
        nnotin = nnotin + 1
                
    
    
    if runSDSS:
        #For a given object position go through all the SDSS bands
        SDSSimages = []
        #3543,4770,6231,7625,9134
        SDSSwavelengths = {'u':3543,'g':4770,'r':6231,'i':7625,'z':9134}
        for band in ('u','g','r','i','z'):
            
            if downloadSDSS:
                try:
                    print('downloading from astroquery SDSS')
                    SDSSplate = getPlateFits((str(source[1]) + 'd ' + str(source[2]) + 'd' ), band)
                    
                except:
                    print(source[0], ' is not in SDSS')
                    continue
                SDSSwcs = WCS(SDSSplate[0][0].header)
                SDSSxpix, SDSSypix = SDSSwcs.all_world2pix(source[1],source[2],0)
                #print('SDSS pixelnumbers are', SDSSxpix, SDSSypix)
                SDSSposition = (SDSSxpix, SDSSypix)
                SDSSsize = 150
                
                SDSSCutout = Cutout2D(SDSSplate[0][0].data, SDSSposition, SDSSsize) #WCS????
                try:
                    fits.writeto(workDir + 'SDSS/' + source[0] +'-'+ band + '.fits', 
                                 SDSSCutout.data)
                                 #Cutout2D(SDSSplate[0][0].data, SDSSposition, SDSSsize).data)                 
                except:

                    os.remove(workDir + 'SDSS/' + source[0] +'-'+  band + '.fits')
                    fits.writeto(workDir + 'SDSS/' + source[0] +'-'+ band + '.fits', 
                                 SDSSCutout.data)            
            else:
                try:                
                    SDSSCutout = fits.open(workDir + 'SDSS/' + source[0] + '-' + band + '.fits')
                except:
                    print('no ', band, ' band fits cutout for ', source[0])
                    continue
                     
            if band == 'g' and plotfigs:
                #CHECK IMAGE
                fig = plt.figure()
                #norm=LogNorm(),
                
                plt.imshow(SDSSCutout.data, cmap='gray',  interpolation='none')
                plt.title(source[0])
                #fig.savefig(workDir + 'output/' + source[0] + '.png')
                plt.close()
            

                
            
            SDSSimages = SDSSimages + [[band,SDSSCutout,SDSSwavelengths[band] ]]

            
            if band == 'g' and runGalfit:
                #run Galfit on g band
                optimise(source[0] + '-' + band + '.fits', workDir + 'SDSS/')
                nSDSS = nSDSS + 1
                
        #print(len(SDSSimages))
        #runGalfitM = False
        if runGalfitM:
            #run GalfitM on all bands
            optimiseM(SDSSimages,source[0])
        

                
imageY.close()
imageJ.close()
imageH.close()
imageKs.close()
imageNB118.close()

print('Code ran on ', nisin, 'objects in UltraVISTA and ignored ', nnotin, ' out of the field.')
print('Galfit ran on ',nSDSS, ' SDSSobjects')
print('GalfitM ran on ????????? objects')

t1= time.time()
print('Code completed in a total of ', t1-t0, ' seconds.')