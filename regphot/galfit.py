#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 14 16:21:46 2016

Defines a function to call Galfit for a single fits image sersic optimisation

Also a funtion to calculate chisquared for a given sersic parametrerisation

@author: rs548
"""

from __future__ import division, print_function

import subprocess
#import astropy
from astropy.io import fits
from os import listdir
from os import getcwd
from os import rename
from os import chdir
from os import remove
from textwrap import dedent
import time

trashDir = '/Users/rs548/Documents/Science/PeteHurley/TRASH/'
#chdir(workDir + 'TRASH/')

wkdir = getcwd()
#print('at galfit import wkdir is' + wkdir)

def optimise(fitsfile, folder, source, field=None,**kwargs):
    #print('galfit.py has been called')
    if field == None:
        print('no field assigned for Galfit run, assuming SDSS')
        field = 'SDSS'
    
    #field = 'SDSS'
    if field == 'UltraVISTA':
        plateScale = '0.14997825'
        magphotzero = '30.0'
        inputDir = '/Users/rs548/Documents/Science/PeteHurley/UV/'
        outputDir = '/Users/rs548/Documents/Science/PeteHurley/UVG/'

        
    elif field == 'SDSS':
        plateScale = '0.396127'
        # SDSS from http://classic.sdss.org/dr7/algorithms/fluxcal.html
        #magphotzero = 24.63 #u band 
        magphotzero = '25.11' #g band
        #magphotzero = 24.80 #r band
        #magphotzero = 24.36 #i band
        #magphotzero = 22.83 #z band
        inputDir = '/Users/rs548/Documents/Science/PeteHurley/SDSS/'
        outputDir = '/Users/rs548/Documents/Science/PeteHurley/SG/' 


        
    filename=str(fitsfile)
    fileRoot = filename[:-5]
    restart = False
    
    rGuess = source['expRad_r_sdss'] / float(plateScale)
    axisRatioGuess = source['expAB_r_sdss']
    method = 0
    psf = 'none' #folder + 'PSF-g.fits' #'PSF-g.fits' psfNOTOK.fits
    inputText = """
    ===============================================================================
    # IMAGE and GALFIT CONTROL PARAMETERS
    A) {folder}{filename}            # Input data image (FITS file)
    B) {outputDir}{fileRoot}-output.fits        # Output data image block
    C) none                # Sigma image name (made from data if blank or "none") 
    D) {psf}   #        # Input PSF image and (optional) diffusion kernel
    E) 1                   # PSF fine sampling factor relative to data 
    F) none                # Bad pixel mask (FITS image or ASCII coord list)
    G) /Users/rs548/Documents/Science/PeteHurley/12xyequal.constraints                # File with parameter constraints (ASCII file) 
    H) 1   150   1    150   # Image region to fit (xmin xmax ymin ymax)
    I) 100    100          # Size of the convolution box (x y)
    J) {magphotzero}              # Magnitude photometric zeropoint (UltraVISTA = 30.0) ???
    K) {plateScale} {plateScale}        # Plate scale (dx dy)    [arcsec per pixel] (UltraVISTA= 0.14997825  0.14997825 )
    O) regular             # Display type (regular, curses, both)
    P) {method}                   # Choose: 0=optimize, 1=model, 2=imgblock, 3=subcomps

    
    # INITIAL FITTING PARAMETERS
    #
    #   For object type, the allowed functions are: 
    #       nuker, sersic, expdisk, devauc, king, psf, gaussian, moffat, 
    #       ferrer, powsersic, sky, and isophote. 
    #  
    #   Hidden parameters will only appear when they're specified:
    #       C0 (diskyness/boxyness), 
    #       Fn (n=integer, Azimuthal Fourier Modes),
    #       R0-R10 (PA rotation, for creating spiral structures).
    # 
    # -----------------------------------------------------------------------------
    #   par)    par value(s)    fit toggle(s)    # parameter description 
    # -----------------------------------------------------------------------------
    
    # Object number: 1
     0) sersic                 #  object type
     1) 75.0  75.0  1 1        #  position x, y
     3) 20.0    1              #  Integrated magnitude	
     4) 4      1             #  R_e (half-light radius)   [pix]
     5) 4.0         1             #  Sersic index n (de Vaucouleurs n=4) 
     6) 0.0000      0          #     ----- 
     7) 0.0000      0          #     ----- 
     8) 0.0000      0          #     ----- 
     9) 1      1             #  axis ratio (b/a)  
    10) 0.0    1               #  position angle (PA) [deg: Up=0, Left=90]
     Z) 0                      #  output option (0 = resid., 1 = Don't subtract) 

    # Object number: 2
     0) expdisk                #  object type
     1) 75.0  75.0  1 1        #  position x, y
     3) 17.0        1          #  Integrated magnitude	
     4) {rGuess}        1          #  R_s    [pix]
     9) {axisRatioGuess}      1             #  axis ratio (b/a)  
    10) 0.0    1               #  position angle (PA) [deg: Up=0, Left=90]
     Z) 0                      #  output option (0 = resid., 1 = Don't subtract) 
    
    # Object number: 3
     0) sky                    #  object type
     1) 0.001      1          #  sky background at center of fitting region [ADUs]
     2) 0.0000      0          #  dsky/dx (sky gradient in x)
     3) 0.0000      0          #  dsky/dy (sky gradient in y)
     Z) 0                      #  output option (0 = resid., 1 = Don't subtract) 
    
    ================================================================================
    """.format(**vars())
    
    inputFile = open(outputDir + fileRoot + '.feedme', "w")
    inputFile.write(inputText)
    inputFile.close()
    
    log_file = open(outputDir + fileRoot + '-log.txt', 'a')
    

    
    #print('at galfit run wkdir is' + wkdir)
    if restart:
        print('restarting with last run parameters')
    else:
        #print( listdir(trashDir))
        for restartFile in listdir(trashDir):
            if restartFile[0:7] == 'galfit.':
                rename(trashDir + restartFile, trashDir + 'wasnotused'  + restartFile)
                #print(trashDir + restartFile, trashDir + 'wasnotused'  + restartFile)
                
            

    #'-imax', '99', 
    print('Galfit.py is about to call galfit on ' + filename)
    #print(inputText)
    t0 = time.time()
    subprocess.call(['/usr/local/bin/galfit', 
                    outputDir 
                    + fileRoot + '.feedme'], 
                    stdout=log_file)
    trun = time.time() - t0
    print('Galfit took ',  round(trun,1), 's to run.')
    
    log_file.close()
    
    #fitsOutput = fits.open('imgblock.fits')
    #fitsOutput.close()
    #chiSq = fitsOutput[2].header['CHISQ']
    #print('ChiSq is ',chiSq )
    
###############################PROFIT INPUT TO COPY##########################
#    profit_model = {'width':  data.image.shape[1],
#                    'height': data.image.shape[0],
#                    'magzero': data.magzero,
#                    'psf': data.psf,
#                    'profiles': {'sersic': sparams}
#                   }
#    print(profit_model)
#    if use_mask:
#        profit_model['calcmask'] = data.calcregion
#    return allparams, np.array(pyprofit.make_model(profit_model))
    
#get chi squared for use by mcmc code
def chisquared(image,model,*args,**kwargs):
    """
    This function should take an input that can be taken by pyprofit and return 
    the same object that pyprofit would return. For rigour any code should
    therefore be able to be run using both pyprofit and Galfit. This will allow 
    comparison and help with publishing given the long standing of Galfit.
    
    
    """
    temp = 'temp/'
    folder = '/Users/rs548/Documents/Science/Blended/'
    imageName = 'tempOut'
    #Call galfit for param values model only
    
    nObjects = len(model)
    
    
    allObjects = ''
    n = 1
    for component in model:
        if component[0] == 'sersic':
            thisObject = """
                         # Object number: {n}
                         0) {component[0]}                 #  object type
                         1) {component[1]}  {component[2]}  1 1  #  position x, y
                         3) {component[3]}    1          #  Integrated magnitude	
                         4) {component[4]}     1          #  R_e (half-light radius)   [pix]
                         5) {component[5]}      1          #  Sersic index n (de Vaucouleurs n=4) 
                         6) 0.0000      0          #     ----- 
                         7) 0.0000      0          #     ----- 
                         8) 0.0000      0          #     ----- 
                         9) {component[9]}      1          #  axis ratio (b/a)  
                         10) {component[10]}    1          #  position angle (PA) [deg: Up=0, Left=90]
                         Z) 0                      #  output option (0 = resid., 1 = Don't subtract) 
                         """.format(**vars())
        elif component[0] == 'sky':
            thisObject = """
                         # Object number: {n}
                         0) {component[0]}                 #  object type
                         1) {component[1]}   1  #  Sky background
                         2) 0.000e+00
                         3) 0.000e+00
                         Z) 0                      #  output option (0 = resid., 1 = Don't subtract) 
                         """.format(**vars())
            
        allObjects = allObjects + dedent(thisObject)     
                    
        n = n + 1
                    
                    
    #return chi squared measure of liklihood
    method = '2'
    inputText = """
===============================================================================
# IMAGE and GALFIT CONTROL PARAMETERS
A) {image}            # Input data image (FITS file)
B) {folder}{temp}{imageName}-output.fits        # Output data image block
C) none                # Sigma image name (made from data if blank or "none") 
D) none #/Users/rs548/Documents/Science/Blended/PSF-g.fits   #        # Input PSF image and (optional) diffusion kernel
E) 1                   # PSF fine sampling factor relative to data 
F) none                # Bad pixel mask (FITS image or ASCII coord list)
G) none                # File with parameter constraints (ASCII file) 
H) 1448   2048   889    1489   # Image region to fit (xmin xmax ymin ymax)
I) 100    100          # Size of the convolution box (x y)
J) 25.110              # Magnitude photometric zeropoint 
K) 0.396127  0.396127       # Plate scale (dx dy)    [arcsec per pixel]
O) regular             # Display type (regular, curses, both)
P) {method}                   # Choose: 0=optimize, 1=model, 2=imgblock, 3=subcomps

# INITIAL FITTING PARAMETERS
#
#   For object type, the allowed functions are: 
#       nuker, sersic, expdisk, devauc, king, psf, gaussian, moffat, 
#       ferrer, powsersic, sky, and isophote. 
#  
#   Hidden parameters will only appear when they're specified:
#       C0 (diskyness/boxyness), 
#       Fn (n=integer, Azimuthal Fourier Modes),
#       R0-R10 (PA rotation, for creating spiral structures).
# 
# -----------------------------------------------------------------------------
#   par)    par value(s)    fit toggle(s)    # parameter description 
# -----------------------------------------------------------------------------

{allObjects}



================================================================================
    """.format(**vars())
    
    inputFile = open(folder + temp 
                     + imageName + '.feedme', "w")
    inputFile.write(inputText)
    inputFile.close()

    log_file = open(imageName + "-log.txt", "a")
    
    #'-imax', '99', 
    print('Galfit.py is about to call galfit for chiSquared')
    print(inputText)
    t0 = time.time()
    subprocess.call(['/usr/local/bin/galfit', 
                     folder + temp 
                     + imageName + '.feedme'], 
                    stdout=log_file)
    print('galfit.chisquared took ',time.time() - t0, 's to run galfit')
    log_file.close()
    
    
    output = fits.open(folder + temp 
              + imageName + '-output.fits')
    #Get chi squared from fits file
    chiSquared = output[2].header['CHISQ']
    #remove(folder + temp 
    #          + imageName + '-output.fits')
    print('galfit.chisquared gives a value of ', chiSquared)
    return chiSquared

    
    