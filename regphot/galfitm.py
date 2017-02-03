#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 14 16:21:46 2016

@author: rs548
"""

from __future__ import division, print_function

import subprocess
#import astropy
#from astropy.io import fits

def optimiseM(images,source):
    print('galfitm.py has been called')
    inputDir = '/Users/rs548/Documents/Science/PeteHurley/SDSS/'
    outputDir = '/Users/rs548/Documents/Science/PeteHurley/SDSSMOutput/'
    filenames=''
    sourcename = source
    bandlabels=''
    n=0
    for band in images:
        if not n == 0:
            filenames = filenames + ','
            bandlabels = bandlabels + ','
        n = n + 1

        filenames = (filenames + inputDir
                                + source + '-' + band[0] + '.fits')
        bandlabels = bandlabels + band[0]


    method = 0
    inputText = """
    ================================================================================
    # GUIDE TO INPUT FILE FOR GALFITM (a product of the MegaMorph project)
    Including multi-band fitting, non-parametric component and MultiNest sampling.
    CSL = comma separated list (must not contain any whitespace)
    Where several lines are given for the same letter code, these are alternative
    examples. The behaviour for multiple lines with the same letter is undefined.
    ================================================================================
    # IMAGE and GALFIT CONTROL PARAMETER
    
    # Input data images (CSL of FITS filenames)
    # the number of input data images defines <nbands>
    # the order of the bands must be maintained in all multi-band options
    # the first band in the list is the 'reference band'
    A) {filenames}
    
    # Band labels (CSL of <nbands> labels containing no whitespace)
    # (these must be unique in a case-insensitive manner)
    # (can be omitted if fitting a single band)
    A1) {bandlabels}
    
    # Band wavelengths (CSL of values)
    # (choice of wavelength units is arbitrary, as long as consistent,
    #  but affects the resulting wavelength-dependence parameters)
    #Ultravista: 1020,1250,1650,2200,1180
    #SDSS: 3543,4770,6231,7625,9134
    A2) 3543,4770,6231,7625,9134
    
    # Output data image block (FITS filename)
    B) {outputDir}{sourcename}-output.fits
    
    # Sigma image name (CSL of <nbands> FITS filenames or "none")
    # (if an individual filename is specified as "none", then that sigma
    #  image will be made from data; if the whole entry consists of just a
    #  single "none", then all sigma images will be made from data.)
    # One can also add a minimum sigma value, such that any galfit-created
    # sigma image will have a minimum of that value times the sky-subtracted
    # input data.
    C) none
 
    # Input PSF image (CSL of <nbands> FITS filenames) 
    # and a single diffusion kernel (FITS filename, # or omitted)
    D) {inputDir}PSF-u.fits,{inputDir}PSF-g.fits,{inputDir}PSF-r.fits,{inputDir}PSF-i.fits,{inputDir}PSF-z.fits
    #D) psf-r.fits,psf-g.fits,psf-i.fits  kernel.fits
    
    # PSF fine sampling factor relative to data 
    E) 1                   
    
    # Bad pixel mask (CSL of <nbands> FITS image or ASCII coord list)
    # (if an individual filename is specified as "none", then a blank
    #  mask will be used; if the whole entry consists of just a single
    #  "none", then all masks will be blank.)
    F) none
    #F) mask-r.fits,mask-g.fits,mask-i.fits
    #F) none,none,none
    #F) mask-r.fits,none,mask-i.fits
    
    # File with parameter constraints (ASCII file)
    G) none # {inputDir}bulgedisklocation.constraints              
    #G) constraints_filename 12xyequal.constraints
    
    # Image region to fit (xmin xmax ymin ymax)
    H) 1    150   1    150
    
    # Size of the convolution box (x y)
    I) 100    100
    
    # Magnitude photometric zeropoint (CSL of <nbands> values)
    #UltraVISTA from fits header 30.0,30.0,30.0,30.0,30.0
    # SDSS from http://classic.sdss.org/dr7/algorithms/fluxcal.html
    J) 24.63,25.11,24.80,24.36,22.83
    
    # Plate scale (dx dy)   [arcsec per pixel]
    #UltraVISTA 0.14997825  0.14997825
    #SDSS 0.396127  0.396127 
    K) 0.396127  0.396127 
    
    # Display type (regular, curses, both)
    O) regular             
    
    # Options: 0=normal run; 1,2=make model/imgblock & quit
    P) 0                   
    
    # Non-parametric component
    U) 0     # Standard parametric fitting
    #U) 1     # Turn on non-parametric component with SED homogenisation
    #U) -1    # Turn on non-parametric component without SED homogenisation
    #U) 10 0.75 25 4 40  # Customise the non-parametric schedule (n,a,b,c,d)
                        # Every n iterations, the nonparametric image is updated
    		    # using a fraction of the filtered residuals, npf:
    		    # npf = a * b^c / (b^c + |i - d|^c),
    		    # where i is the iteration number.
    
    # MultiNest
    V) 0     # Use standard Levenburg-Marquardt algorithm
    #V) 1	 # Use MultiNest sampling algorithm (experimental and slow!)
    #V) 1 0 500 0.8 0.5 100000   # Customise MultiNest options:
                                # ceff,nlive,efr,tol,maxiter
    
    # Output options
    # Valid options: blank,input,model,residual -- as usual
    #                component -- individual model components
    #                psf -- input psf image
    #                sigma -- sigma image (input or created)
    #                mask -- mask image
    #                nonparam -- nonparam image (if appropriate)
    #                datasub -- input minus nonparam image (if appropriate)
    #                itertable -- table of parameters at each iteration
    #W) default  # == blank,input,model,residual and assumed if omitted
    #W) none     # (or any other invalid option) no images output
    #W) model    # or any other valid option
    W) input,model,residual    # or any comma separated list of valid options
    
       
    # Sersic function --------------------------------------------------------------
    
    # Only this first function includes multi-band examples, but the same approach
    # should work for all these functions.
    
    0) sersic     # Object type
    1) 75.  1    # position x [pixel]  (constant with wavelength)
    2) 75.  1    # position y [pixel]
    3) 20.0,20.0,20.0,20.0,20.0  3     # total magnitude in each band
    4) 5.0,5.0,5.0,5.0,5.0   2     # R_e in each band
    5) 2.0,2.0,2.0,2.0,2.0  1     # Sersic exponent in each band
    9) 1.0,1.0,1.0,1.0,1.0  1     # axis ratio (b/a) in each band
    10) 0.0            1     # position angle (PA), same value in each band
    Z) 0                  #  Skip this model in output image?  (yes=1, no=0)
    
    #0) sersic     # Object type
    #1) 300.  1    # position x [pixel]  (constant with wavelength)
    #2) 350.  1    # position y [pixel]
    #3) 20.00,0.2,0.0  1,1,1 cheb # total magnitude and coeffs of wavelength dependence
    #4) 4.30,0.1,0.0   1,1,0 cheb # R_e and coeffs of wavelength dependence
    #5) 5.20,0.0,0.0   1,0,0 cheb # Sersic exponent and coeffs of wavelength dependence
    #9) 0.30,0.0,0.0   1     cheb # axis ratio (b/a) and coeffs of wavelength dependence
    #10) 10.0           1     cheb # position angle (PA) and coeffs of wavelength dependence
    #Z) 0                  #  Skip this model in output image?  (yes=1, no=0)
    
    
  
          #  Skip this model in output image?  (yes=1, no=0)
    
     # sky --------------------------------------------------------------------------
    
     0) sky
     1) 0.005       0       # sky background       [ADU counts]
     2) 0.000      0       # dsky/dx (sky gradient in x) 
     3) 0.000      0       # dsky/dy (sky gradient in y) 
     Z) 0                  #  Skip this model in output image?  (yes=1, no=0)
    
    
    
    """.format(**vars())
    
    inputFile = open(outputDir 
                     + sourcename + '.feedme', "w")
    inputFile.write(inputText)
    inputFile.close()
    
    log_file = open(outputDir + sourcename + "-log.txt", "a")
    
    #'-imax', '99', 
    print('Galfitm.py is about to call galfitm')
    print(inputText)
    subprocess.call(['/Users/rs548/Documents/Code/galfitm/galfitm-1.2.1-osx', 
                    outputDir 
                    + sourcename + '.feedme'], 
                    stdout=log_file)
    
    log_file.close()
    
    #fitsOutput = fits.open('imgblock.fits')
    #fitsOutput.close()
    #chiSq = fitsOutput[2].header['CHISQ']
    #print('ChiSq is ',chiSq )
    