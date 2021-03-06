#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 11 13:20:20 2017

Some astroquery functions for SDSS

Also a couple fo functions for dealing with Simard's dat files and querying SDSS
for objID info to out put Simard B/D fits for a given object.

@author: rs548
"""

from astroquery.sdss import SDSS
from astropy import coordinates
from astropy.io import fits
from astropy.wcs import WCS
from astropy.nddata import Cutout2D

import os
import csv

#pos = coordinates.SkyCoord('0h8m05.63s +14d50m23.3s', frame='icrs')
#xid = SDSS.query_region(pos, spectro=True)
#sp = SDSS.get_spectra(matches=xid)
#im = SDSS.get_images(matches=xid, band='g')
#plate = xid[0]['plate']
#
#print(plate)
#
#for f in os.listdir('/Users/rs548/Documents/Science/PeteHurley/SDSS/'):
#    print(f)

def getPlateFits(position,bandName):
    """Get SDSS plate in fits format
    
    This function takes a position and band and returns the available fits images from 
    SDSS using Astroquery
    
    Parameters
    ----------
    position: string
        Ra/Dec position with units e.g. ' 35.2345342d 5.3452346d'
    bandName: string
        Name of band e.g. i from ugrizy


    Returns
    -------
    images astropy.table.Table
        i think the return is an HDUlist of all the available fits for a position/band.
    """
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

def cutoutFits(image,position,size):
    """"Cutout an image from within a fits file and return a fits image with correct wcs
    
    The Cutout2D object does not have the normal fits methods of associated header data
    So this function exists do cutout the array and then put it back in the image 
    with corrected wcs.
    
    Parameters
    ----------
    image: ndarray
        The fits image
    position: string
        Name of band e.g. i from ugrizy
    size: int, array-like, astropy.units.Quantity
    	Size to cut in pixels or can be angle if use astropy Quantity


    Returns
    -------
    images astropy.io.fits
        This will be the same as the input image but a cutout
    
    """"
    cutout = Cutout2D(image[0].data, 
                      position, 
                      size,
                      mode='partial', 
                      fill_value=0.0,
                      image[0].wcs)
    image.data = cutout.data

    image[0].header['CRPIX1'] = (xpixoffset,
    'X of reference pixel (modified rs548)')
    image[0].header['CRPIX2'] = (ypixoffset,
    'Y of reference pixel (modified rs548)')
    image[0].header['CRVAL1'] = (position[0],
    'RA of reference pixel (deg) (modified rs548)')
    image[0].header['CRVAL2'] = (position[0],
    'Dec of reference pixel (deg) (modified rs548)')
    image[0].header['NAXIS1'] = (len(cutout.data),
    'xpixels (modified by catalogue.py)')                
    image[0].header['NAXIS2'] = (len(cutout.data),
    'ypixels (modified by catalogue.py)') 
    return image
    #try:
    #    UVCutoutFits.writeto(workDir + UVfolder 
    #                             + source[0] +'-'+ band[0] + '.fits')
    #except:
    #    remove(workDir + UVfolder + source[0] +'-'+  band[0] + '.fits')
    #    UVCutoutFits.writeto(workDir + UVfolder 
    #                             + source[0] +'-'+  band[0] + '.fits')
    



def objid2dr7id(objID):
    """
    This function takes a SDSS DR13 objID and returns a DR7 objID based on
    RA/DEC
    """
    
    sql = """
SELECT
 p.ra,p.dec
FROM PhotoObj AS p
WHERE p.objid = {objID}
    """.format(**vars())
    result = SDSS.query_sql_async(sql)
    radecstring = result.content.split('\n')
    try:    
        ra,dec = radecstring[2].split(',')
    except:
        return False, 'empty'
    ra,dec = float(ra),float(dec)
    print( ra,dec)
    sql2 = """
SELECT
 p.objid,p.ra,p.dec
FROM PhotoObj AS p
JOIN dbo.fGetNearestObjEq({ra},{dec}, 1) AS pN
ON p.objID = pN.objID
    """.format(**vars())
    result2 = SDSS.query_sql_async(sql2, data_release=7)
    datastring = result2.content
    try:    
        data = datastring.split('\n')[1]
        dr7objid,ra2,dec2 = data.split(',')
        ra2,dec2 = float(ra2),float(dec2)
        print(dr7objid,ra2,dec2 )
    except:
        print('No neighbour found for ',objID)
        dr7objid = 'empty'
        

    try:
        assert abs(ra - ra2) < 0.0003 #1 arcsec = 0.000277 deg
        assert abs(dec - dec2) < 0.0003
        matchfound = True
    except:
        print('No match found for ',objID)
        matchfound = False
    return matchfound,dr7objid
    
def getsimardrow(dr7id,simarddata):

    datfile = open(simarddata, 'rb') 
    reader = csv.reader(datfile, delimiter = ' ')
        #outputary = [row for row in datfile if row[0]==dr7id]
    outputary  = ['empty']
    for row in reader:
        print( row)
        if row[0] == dr7id:
            outputary = row
        
    return outputary
       
    
def returncsv(inputcatalogue,simarddata,outputfile):
    #for every source in catalogue check Simard List 
    sourcelist = fits.open(inputcatalogue)
    outcsvfile = open(outputfile, 'w')
    simardwriter = csv.writer(outcsvfile, delimiter=',')
    simardwriter.writerow(['sourcename','objID_sdss', 'objID_sdss(dr7)'])
    for source in sourcelist[1].data:
        matchfound,dr7id = objid2dr7id(source['objID_sdss'])
        if matchfound:
            simardwriter.writerow([source[0],source['objID_sdss']] 
            + getsimardrow(dr7id,simarddata))
        else:
            print('Writing blank line for ', source[0])
            simardwriter.writerow([source[0],source['objID_sdss'],dr7id,'problem with this source'])
    outcsvfile.close()        
    #If source in Simard list then write the Simard data to CSV    
    
    
    


if __name__ == '__main__':
    wd = '/Users/rs548/Documents/Science/PeteHurley/'
    #dr7id = objid2dr7id(str(1237654600489566525)) #587722952230175138  1237654600489566525
    #xmm-hyperleda-sdss-dr13-2massXSC-4reg-flagged.fits
    #cosmos-hyperleda-sdss-dr13-2massXSC-4reg-flagged_HELP.fits
    returncsv(wd + 'xmm-hyperleda-sdss-dr13-2massXSC-4reg-flagged.fits',
             wd + 'Simard/fulltable1.dat',
             wd + 'Simard/Simardcross-xmm.csv')
    #sim = getsimardrow('587739630633418838','/Users/rs548/Documents/Science/PeteHurley/Simard/fulltable1.dat')
    
