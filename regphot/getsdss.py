#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 11 13:20:20 2017

Some astroquery functions for SDSS

@author: rs548
"""

from astroquery.sdss import SDSS
from astropy import coordinates
from astropy.io import fits
import os

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

def objid2dr7id(objID):
    sql = """
SELECT
 p.ra,p.dec
FROM PhotoObj AS p
WHERE p.objid = {objID}
    """.format(**vars())
    result = SDSS.query_sql_async(sql)
    radecstring = result.content.split('\n')
    ra,dec = radecstring[2].split(',')
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
    data = datastring.split('\n')[1]
    dr7objid = data.split(',')[0]
    return dr7objid
    
    
def returncsv(inputcatalogue,simarddata):
    #for every source in catalogue check Simard List 
    sourcelist = fits.open(inputcatalogue)
    for source in sourcelist:
        dr7id = objid2dr7id(source['objID_sdss'])

    #If source in Simard list then write the Simard data to CSV    
    
    
    


if __name__ == '__main__':
    dr7id = objid2dr7id(str(1237654600489566525)) #587722952230175138  1237654600489566525
    
    
