#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 13:42:05 2017

@author: rs548
"""

from astropy.io import fits

import matplotlib.pyplot as plt

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
            images.close()
            numberOutputs = numberOutputs + 1
            #remove('/Users/rs548/Documents/Science/PeteHurley/SDSS/' + filename)
            
printGraphs('/Users/rs548/Documents/Science/PeteHurley/SDSS/')

def generateTables(folder):
    numberObjects = 0
    writer = csv.writer(open( folder + 'out.csv', 'wb'))
    headings = ['CHISQ',
                'RA',
                'DEC',
                'R_e']
    for filename in listdir(folder):
        if filename[-11:] == 'output.fits':
            output = fits.open(folder + filename)
            allParams = output[0]['CHISQ']
            writer.writerow(allParams)
            
        

