#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 13:42:05 2017

This script is for post processing data produced by other parts of the codebase.

It contains function definitions which may also be useful to other modules.

eventually this should write the models to a csv for use by XID+

@author: rs548
"""

from astropy.io import fits

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

import numpy as np
import pandas as pd
from scipy import stats, integrate
import seaborn as sns
sns.set(color_codes=True)

from os import listdir
from os import getcwd
from os import remove

import csv


def checkcatalogue(sdssid1,cat2):
    #First get the SDSSid for the provided objID
    matchfound = False
    for source in cat2:
        if sdssid1 == sdssid2:
            matchfound = True
    return matchfound
    
def comparecatalogues(cat1,cat2):
    print('There are ', len(cat1), ' objects in catalogue 1')
    print('There are ', len(cat2), ' objects in catalogue 2')
    nmatch=0
    for source1 in cat1:
        if checkcatalogue(source1['SDSSid'],cat2):
            nmatch += 1
            
    print('There are ', nmatch, ' objects from catalogue 1 in catalogue 2')

def printGraphs(folder):
    numberOutputs = 0
    for filename in listdir(folder):
        if filename[-11:] == 'output.fits':
            images = fits.open(folder + filename)
            fig = plt.figure()
            fig.suptitle(filename)
            #norm=LogNorm(),
            plt.subplot(131)
            plt.imshow(images[1].data, cmap='gray',  interpolation='none')
            plt.title('Image')
            plt.subplot(132)
            plt.imshow(images[2].data, cmap='gray',  interpolation='none')
            plt.title('Model')
            plt.subplot(133)
            plt.imshow(images[3].data, cmap='gray',  interpolation='none')
            plt.title('Residual')
        
            
            plt.show()
            plt.close()
            images.close()
            numberOutputs = numberOutputs + 1
            #remove('/Users/rs548/Documents/Science/PeteHurley/SDSS/' + filename)
            

def print5bandGraphs(folder,band):
    numberOutputs = 0
    for filename in listdir(folder):
        if filename[-11:] == 'output.fits':
            images = fits.open(folder + filename)
            fig = plt.figure()
            fig.suptitle(filename)
            #norm=LogNorm(),
            plt.subplot(131)
            plt.imshow(images[band + 1].data, cmap='gray',  interpolation='none')
            plt.title('Image')
            plt.subplot(132)
            plt.imshow(images[band + 6].data, cmap='gray',  interpolation='none')
            plt.title('Model')
            plt.subplot(133)
            plt.imshow(images[band + 11].data, cmap='gray',  interpolation='none')
            plt.title('Residual')
        
            
            plt.show()
            plt.close()
            images.close()
            numberOutputs = numberOutputs + 1
            #remove('/Users/rs548/Documents/Science/PeteHurley/SDSS/' + filename)
            


def oneModel(output):
    lognorm = True
    image = fits.open(output)
    fig = plt.figure()
    fig.suptitle(output)
            #norm=LogNorm(),
    
    
    plt.imshow(image[1].data, cmap='gray',  interpolation='none',norm=LogNorm())
    plt.title('Image')
   
    fig = plt.figure()
    plt.imshow(image[2].data, cmap='gray',  interpolation='none',norm=LogNorm())
    plt.title('Model')
   
    fig = plt.figure()
    plt.imshow(image[3].data, cmap='gray',  interpolation='none',norm=LogNorm())
    plt.title('Residual')
            

    image.close()


        


def generateTables(folder,bandnames=['u','g','r','i','z']):
    numberObjects = 0
    writer = csv.writer(open( folder + 'out.csv', 'wb'))
    paramnames = ['OBJID',
                'CHISQ',
                'RA',
                'DEC',
                'R_e']
    writer.writerow(paramnames)
    for filename in listdir(folder):
        if filename[-11:] == 'output.fits':
            output = fits.open(folder + filename)
            numBands = ((len(output) -1)/3) -1
            for band in range(0,numBands):
                allbandparams = []
                for param in paramnames:                       
                    allbandparams += [output[band+numBands].header[param]]
                writer.writerow(allbandparams)
    return writer
            
           
            
            
def printAllBandGraphs(folder):
    """
    Go though a folder and print all the passband images/models/residuals
    for every Galfit output file. Will have to be modified for pyprofit.
    """
    numberOutputs = 0
    for filename in listdir(folder):
        if filename[-11:] == 'output.fits':
            images = fits.open(folder + filename)
            numBands = ((len(images) -1)/3) -1
            print(numBands)
            
            fig,axarr = plt.subplots(nrows=numBands, ncols=3, sharex=True, 
                                     sharey=True, figsize=(10,10))
            plt.suptitle(filename)
            #norm=LogNorm(),
            axarr[0,0].set_title('Image')
            axarr[0,1].set_title('Model')
            axarr[0,2].set_title('Residual')
            
            for band in range(0,numBands):
                axarr[band,0].imshow(images[band].data, 
                     cmap='gray',  interpolation='none')
                axarr[band,1].imshow(images[band + numBands].data, 
                     cmap='gray',  interpolation='none')
                axarr[band,2].imshow(images[band + 2*numBands].data, 
                     cmap='gray',  interpolation='none')

        
            
            plt.show()
            plt.close()
            images.close()
            numberOutputs = numberOutputs + 1
            #remove('/Users/rs548/Documents/Science/PeteHurley/SDSS/' + filename)
            print('done a file')
    plt.close('all')
    
def generateTables(folder,bandnames=['u','g','r','i','z']):
    """
    A function to go through a folder of GalfitM output files (fits) and
    print all the Sersic or Sersic/bulge parameters to a CSV
    """
    numberObjects = 0
    outfile = open( folder + 'out.csv', 'w')
    writer = csv.writer(outfile)
    #Define a non general set of params to pull out for a SDSS UGRIZ fit
    paramnames = ['DATAIN_U',
                'CHISQ',
                '1_XC_U',
                '1_YC_U',
                '1_MAG_U', '1_MAG_G', '1_MAG_R','1_MAG_I','1_MAG_Z',
                '1_RE_U',
                '1_N_U',
                '1_AR_U',
                '1_PA_U']
    writer.writerow(paramnames)
    for filename in listdir(folder):
        if filename[-11:] == 'output.fits':
            output = fits.open(folder + filename)
            numBands = ((len(output) -1)/3) -1
            #for band in range(0,numBands):
            allbandparams = []
            for param in paramnames:    
                #print(band,numBands,param)                   
                allbandparams += [output[6].header[param]]
            writer.writerow(allbandparams)
    return writer



    
        
    

    
    
        
if __name__ == '__main__':

    #printGraphs('/Users/rs548/Documents/Science/PeteHurley/UVG/')
    #printAllBandGraphs('/Users/rs548/Documents/Science/PeteHurley/SDSS-M-BD/')
    #print5bandGraphs('/Users/rs548/Documents/Science/PeteHurley/SM/',3)
    #oneModel('/Users/rs548/Documents/Science/Blended/g-output.fits')
    #generateTables('/Users/rs548/Documents/Science/PeteHurley/SDSS-XM/')
    
    
