#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 13 12:23:18 2016

A principle module of the regphot python package.

This script goes though a supplied catalogue and for every object it can:
    
Cutout a postage stamp from a supplied fits image.
Download a fits plate and cutout the required postage stamp from SDSS using astroquery

It can then run either Galfit on one band or GalfitM on a set of bands.

Eventually this will manage a batch run for my deblending algorithms.

Raphael Shirley
Daphne Jackson Research Fellow
Astronomy Centre, University of Sussex
@author: rs548
@raphaelshirley
www.raphaelshirley.co.uk
"""
from __future__ import division, print_function

from os import chdir
from os import remove
from astropy.io import fits
from astropy.wcs import WCS
from astropy.nddata import Cutout2D #this produces an inheritor to the numpy
# array class and can not be saved to a fits file
#from astropy import units as u
import matplotlib.pyplot as plt
#from matplotlib.colors import LogNorm
import time
t0 = time.time()


# MY MODULES
#from regphot import getsdss
from utils import getPlateFits
from galfit import optimise
from galfitm import optimiseM




#############################INPUTS############################################
if __name__ == '__main__':
    workDir = '/Users/rs548/Documents/Science/PeteHurley/'
    chdir(workDir + 'TRASH/')
    UVfolder = 'UV/'
    SDSSfolder = 'SDSS/'
    
    #Open FITS catalogue
    COSMOSlist = 'cosmos-hyperleda-sdss-dr13-2massXSC-4reg-flagged_HELP.fits'
    XMMlist = 'xmm-hyperleda-sdss-dr13-2massXSC-4reg-flagged.fits'
    hdulist = fits.open(workDir + COSMOSlist)
    
    
    #Options
    runUltraVISTA = False
    makeUVcutouts = False
    
    runSDSS = True
    downloadSDSS = False
    
    runGalfit = False
    runGalfitM = True
    plotfigs = False
    options = {'runUltraVISTA':False,
               'makeUVcutouts':False,
               'runSDSS':True,
               'downloadSDSS':False,
               'runGalfit':False,
               'runGalfitM':True,
               'plotfigs':False}
    
    
    
    dimension = 150   #Size of cutout to fit (pixels)
    
    
    #Location of the UltraVISTA FITS images
    UVYfits = '/Volumes/Raph500/data/ADP.2016-03-17T08-31-54.127.fits'
    UVJfits = '/Volumes/Raph500/data/ADP.2016-03-17T08-31-54.147.fits'
    UVHfits = '/Volumes/Raph500/data/ADP.2016-03-17T08-31-54.060.fits'
    UVKsfits = '/Volumes/Raph500/data/ADP.2016-03-17T08-31-54.113.fits'
    UVNB118fits = '/Volumes/Raph500/data/ADP.2016-03-17T08-31-54.107.fits'
    
    UVfits = (('Y', UVYfits, 1020),
              ('J', UVJfits, 1250),
              ('H', UVHfits, 1650),
              ('Ks', UVKsfits, 2200),
              ('NB118', UVNB118fits, 1180))




#############################INPUTS############################################

def makeUVcutouts(source,dimension,UVfits,workdir,UVfolder):
    """
    Take in a fits image and cutout a postage stamp of a given size at a given 
    location
    """
    UVwcs = WCS(UVfits[0][1])
    UVxpix, UVypix = UVwcs.all_world2pix(source['ra_sdss'],source['dec_sdss'],0)
    #top left is 150.77177 +02.80870
    #bottom right is 149.29941 +01.59448
    
    UVsize = (dimension,dimension)
    UVposition = (UVxpix , UVypix )

    UVimages=[]
    for band in UVfits:
         

        UVCutoutFits = fits.open(band[1])       
        cutout = Cutout2D(UVCutoutFits[0].data, UVposition, UVsize,
                      mode='partial', fill_value=0.0)
        UVCutoutFits[0].data = cutout.data
        xpixoffset = - (UVxpix -(UVsize[0]/2) - UVCutoutFits[0].header['CRPIX1'])
        ypixoffset = - (UVypix -(UVsize[1]/2) - UVCutoutFits[0].header['CRPIX2'])
        UVCutoutFits[0].header['CRPIX1'] = (xpixoffset,
        'X of reference pixel (modified rs548)')
        UVCutoutFits[0].header['CRPIX2'] = (ypixoffset,
        'Y of reference pixel (modified rs548)')
        UVCutoutFits[0].header['CRVAL1'] = (source['ra_sdss'],
        'RA of reference pixel (deg) (modified rs548)')
        UVCutoutFits[0].header['CRVAL2'] = (source['dec_sdss'],
        'Dec of reference pixel (deg) (modified rs548)')
        UVCutoutFits[0].header['NAXIS1'] = (150,
        'xpixels (modified by catalogue.py)')                
        UVCutoutFits[0].header['NAXIS2'] = (150,
        'ypixels (modified by catalogue.py)') 
        try:
            UVCutoutFits.writeto(workDir + UVfolder 
                                 + source[0] +'-'+ band[0] + '.fits')
        except:
            remove(workDir + UVfolder + source[0] +'-'+  band[0] + '.fits')
            UVCutoutFits.writeto(workDir + UVfolder 
                                 + source[0] +'-'+  band[0] + '.fits') 
         
        UVimages = UVimages +[[band[0],UVCutoutFits,band[2]]]
    return UVimages
    
def openimages(source,bandinfo,workdir,folder):
    """
    Open an array of fits images, one for each band
    """
    images=[]
    for band in bandinfo:

        try:
            CutoutFits = fits.open(workDir + folder + 
                                       source[0] + '-' + band[0] + '.fits')
        except:
            print('No UltraVISTA fits cutout for ',source[0])
         
        images = images +[[band[0],CutoutFits,band[2]]]
    return images
    
def downloadSDSScutouts(source,size,SDSSbandinfo):
    """
    Download the fits plate for a given location for each band and cutout a 
    postage stamp and put it in the appropriate folder
    """

def UltraVISTA(source, 
               UVfits,
               size = dimension,
               makeUVcutouts = False,
               runGalfit=False, 
               runGalfitM=False,
               infolder = None,
               outfolder = None):
    """
    for a given source in a catalogue open/create a fits image and run Galfit
    if requested
    """
    if makeUVcutouts:
        images = makeUVcutouts(source,size,UVfits)
    else:
        images = openimages(source,UVfits)
        
    if runGalfit:
        optimise(images[0],source[0],field = 'UltraVISTA')
    if runGalfitM:
        optimiseM(images,source[0],)
        
    
    
def SDSS(source,
         size = dimension,
         downloadSDSS = False,
         runGalfit=False, 
         runGalfitM=False,
         infolder = None,
         outfolder = None):
    """
    for a given source download and cutout a fits if required and run Galfit
    or GalfitM (or eventually pyprofit)
    """
    SDSSbandinfo = (('u', '', 354.3),
                    ('g', '', 477.0),
                    ('r', '', 623.1),
                    ('i', '', 762.5),
                    ('z', '', 913.4))
    if downloadSDSS:
        images = downloadSDSScutouts(source,size,SDSSbandinfo)
    else:
        images = openimages(source,infolder)
        
    if runGalfit:
        optimise(images[0],source[0],field = 'UltraVISTA')
    if runGalfitM:
        optimiseM(images,source[0],)
    
def runcatalogue(catfile,UVfits= None,runUltraVISTA = False, runSDSS = False):
    """
    Loop over a catalogue and run the UltraVISTA/SDSS image cutting and Galfit
    optimisation code. Reimplement the scripting below as functions to make 
    this operate like a proper Python package.
    """
    nrc = 0
    hdulist = fits.open(catfile)
    for source in hdulist:
        if runUltraVista:
            UltraVISTA(source,UVfits)
        if runSDSS:
            SDSS(source)
        nrc += 1
    print( nrc)
        
    
nisin = 0
nnotin = 0 
nSDSS = 0
nGalfitM = 0
n = 0 
#Go through every object, check it is in field and then do stuff with it
#I deally this will eventually cutout an appropriate fits file and run
#Galfit or Galfitm on it 
for source in hdulist[1].data:
    objectt0 = time.time()
    n += 1
    print('Now running on object number ',n,': ', source[0])
    #print(source[0])
    
    #Should replace this bit with MOC test
    if (runUltraVISTA
        and source['ra_sdss'] <= 150.77177 
        and source['ra_sdss'] >= 149.29941 
        and source['dec_sdss'] <= 2.80870 
        and source['dec_sdss'] >= 1.59448):
        print(source[0], ' at', source['ra_sdss'],source['dec_sdss'], 'is in UltraVISTA')
        nisin = nisin + 1

        #Open Y band for coordinate transformations
        UVwcs = WCS(UVfits[0][1])
        UVxpix, UVypix = UVwcs.all_world2pix(source['ra_sdss'],source['dec_sdss'],0)
        #top left is 150.77177 +02.80870
        #bottom right is 149.29941 +01.59448
        
        UVsize = (dimension,dimension)
        UVposition = (UVxpix , UVypix )

        UVimages=[]
        for band in UVfits:
                 
            if makeUVcutouts:
                UVCutoutFits = fits.open(band[1])       
                cutout = Cutout2D(UVCutoutFits[0].data, UVposition, UVsize,
                              mode='partial', fill_value=0.0)
                UVCutoutFits[0].data = cutout.data
                xpixoffset = - (UVxpix -(UVsize[0]/2) - UVCutoutFits[0].header['CRPIX1'])
                ypixoffset = - (UVypix -(UVsize[1]/2) - UVCutoutFits[0].header['CRPIX2'])
                UVCutoutFits[0].header['CRPIX1'] = (xpixoffset,
                'X of reference pixel (modified rs548)')
                UVCutoutFits[0].header['CRPIX2'] = (ypixoffset,
                'Y of reference pixel (modified rs548)')
                UVCutoutFits[0].header['CRVAL1'] = (source['ra_sdss'],
                'RA of reference pixel (deg) (modified rs548)')
                UVCutoutFits[0].header['CRVAL2'] = (source['dec_sdss'],
                'Dec of reference pixel (deg) (modified rs548)')
                UVCutoutFits[0].header['NAXIS1'] = (150,
                'xpixels (modified by catalogue.py)')                
                UVCutoutFits[0].header['NAXIS2'] = (150,
                'ypixels (modified by catalogue.py)') 
                try:
                    UVCutoutFits.writeto(workDir + UVfolder 
                                         + source[0] +'-'+ band[0] + '.fits')
                except:
                    remove(workDir + UVfolder + source[0] +'-'+  band[0] + '.fits')
                    UVCutoutFits.writeto(workDir + UVfolder 
                                         + source[0] +'-'+  band[0] + '.fits') 
                
            else:
                try:
                    UVCutoutFits = fits.open(workDir + UVfolder + 
                                               source[0] + '-' + band[0] + '.fits')
                except:
                    print('No UltraVISTA fits cutout for ',source[0])
                 
                UVimages = UVimages +[[band[0],UVCutoutFits,band[2]]]
                
            if plotfigs:
                fig = plt.figure()
                #norm=LogNorm(),
                plt.imshow(UVCutoutFits[0].data, cmap='gray',  interpolation='none')
                plt.title(source[0])
                fig.savefig(workDir + 'pics/' + source[0] + '.png')
                plt.close('all')                    

            
        galfitBand = 'Y'                
        if runGalfit:
            optimise(source[0] + '-' + galfitBand + '.fits', workDir + UVfolder,
                     source, field='UltraVISTA')
        if runGalfitM:
            optimiseM(UVimages,source[0],field='UltraVISTA')
    else:
        #print(source[0], 'is NOT in ', field)
        nnotin = nnotin + 1
                
    
    
    if runSDSS:
        #For a given object position go through all the SDSS bands
        SDSSimages = []

        #354.3,477.0,623.1,762.5,913.4 (nm)
        SDSSwavelengths = {'u':354.3,'g':477.0,'r':623.1,'i':762.5,'z':913.4}
        for band in ('u','g','r','i','z'):
            
            if downloadSDSS:
                try:
                    #print('downloading ', band, ' band from astroquery SDSS')
                    tdownload0 = time.time()
                    SDSSplate = getPlateFits((str(source['ra_sdss']) + 'd ' 
                                              + str(source['dec_sdss']) + 'd' ), band)
                    print('downloaded ', band, 
                    ' band from astroquery SDSS in ', 
                    round(time.time()-tdownload0,1),
                    's.')
                except:
                    print(source[0], ' is not in SDSS')
                    continue
                SDSSwcs = WCS(SDSSplate[0][0].header)
                SDSSxpix, SDSSypix = SDSSwcs.all_world2pix(source['ra_sdss'],
                                                           source['dec_sdss'],0)
                #print('SDSS pixelnumbers are', SDSSxpix, SDSSypix)
                SDSSposition = (SDSSxpix, SDSSypix)
                SDSSsize = dimension 
                
                #The following cuts out a 150 pix array then creates a new fits ith 
                #that data but the original header file for exptime etc but
                #update reference pixels so coordinates are hopefully correct
                SDSSCutout = Cutout2D(SDSSplate[0][0].data, SDSSposition, SDSSsize,
                                      mode='partial', fill_value=0.0) #WCS????
                #SDSSCutout.header = SDSSplate[0][0].header
                SDSSCutoutFits = SDSSplate[0]
                SDSSCutoutFits[0].data = SDSSCutout.data
                SDSSCutoutFits[0].header = SDSSplate[0][0].header
                xpixoffset = - (SDSSxpix -(SDSSsize/2) - SDSSplate[0][0].header['CRPIX1'])
                ypixoffset = - (SDSSypix -(SDSSsize/2) - SDSSplate[0][0].header['CRPIX2'])
                SDSSCutoutFits[0].header['CRPIX1'] = (xpixoffset,
                'X of reference pixel (modified by catalogue.py)')
                SDSSCutoutFits[0].header['CRPIX2'] = (ypixoffset,
                'Y of reference pixel (modified by catalogue.py)')
                SDSSCutoutFits[0].header['NAXIS1'] = (150,
                'xpixels (modified by catalogue.py)')                
                SDSSCutoutFits[0].header['NAXIS2'] = (150,
                'ypixels (modified by catalogue.py)')                
                
                try:
                    SDSSCutoutFits.writeto(workDir + SDSSfolder + source[0] +'-'
                                       + band + '.fits') #.data?
                                 #Cutout2D(SDSSplate[0][0].data, SDSSposition, SDSSsize).data)                 
                except:

                    remove(workDir + SDSSfolder + source[0] +'-'+  band + '.fits')
                    SDSSCutoutFits.writeto(workDir + SDSSfolder + source[0] +'-'
                                       + band + '.fits' ) #.data?   


            else:
                try:                
                    SDSSCutoutFits = fits.open(workDir + SDSSfolder + 
                                               source[0] + '-' + band + '.fits')
                except:
                    print('no ', band, ' band fits cutout for ', source[0])
                    continue
            if band == 'g' and plotfigs:
                #CHECK IMAGE
                fig = plt.figure()
                #norm=LogNorm(),
                
                plt.imshow(SDSSCutoutFits.data, cmap='gray',  interpolation='none')
                plt.title(source[0])
                #fig.savefig(workDir + 'output/' + source[0] + '.png')
                plt.close()
                     
  

                
            
            SDSSimages = SDSSimages + [[band,SDSSCutoutFits,SDSSwavelengths[band] ]]

            
            if band == 'r' and runGalfit:
                #run Galfit on g band
                optimise(source[0] + '-' + band + '.fits', workDir + SDSSfolder, source)
                

        nSDSS = nSDSS + 1        
        #print(len(SDSSimages))
        #runGalfitM = False

        if runGalfitM:
            #run GalfitM on all bands
            optimiseM(SDSSimages,source[0])
            nGalfitM = nGalfitM + 1
        
        
        
        print('Total for ', source[0],' was ',round(time.time()-objectt0,1),' s')
                


print('Code ran on ', nisin, 'objects in UltraVISTA')
print('Code ran on ', nSDSS, 'objects in SDSS')

#print('GalfitM ran on ',nGalfitM,' SDSS objects')

t1= time.time()
print('Code completed in a total of ', round((t1-t0)/60,2), ' minutes.')

