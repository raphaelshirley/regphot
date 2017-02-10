#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 13:42:05 2017

this scirpt is for post processing the outputs from batch runs of galfit or galfity m

eventually this should write the models to a csv for use by XID+

@author: rs548
"""

from astropy.io import fits

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from os import listdir
from os import getcwd
from os import remove

import csv


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




def generateTables(folder):
    numberObjects = 0
    writer = csv.writer(open( folder + 'out.csv', 'wb'))
    headings = ['CHISQ',
                'RA',
                'DEC',
                'R_e']
    writer.writerow(headings)
    for filename in listdir(folder):
        if filename[-11:] == 'output.fits':
            output = fits.open(folder + filename)
            allParams = [output[1]['CHISQ'],
                         output[1]['RA'],
                         output[1]['DEC'],
                         output[1]['R_e']]
            writer.writerow(allParams)
            
def printAllBandGraphs(folder):
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
        
if __name__ == '__main__':
    printAllBandGraphs('/Users/rs548/Documents/Science/PeteHurley/SM/')
    #print5bandGraphs('/Users/rs548/Documents/Science/PeteHurley/SM/',3)
    #oneModel('/Users/rs548/Documents/Science/Blended/g-output.fits')
    #generateTables('/Users/rs548/Documents/Science/PeteHurley/SM/')
