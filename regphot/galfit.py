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

def optimise(fitsfile):
    print('galfit.py has been called')
    filename=str(fitsfile)
    method = 0
    inputText = """
    ===============================================================================
    # IMAGE and GALFIT CONTROL PARAMETERS
    A) /Users/rs548/Documents/Science/PeteHurley/galfit/{filename}            # Input data image (FITS file)
    B) /Users/rs548/Documents/Science/PeteHurley/galfit/{filename}-output.fits        # Output data image block
    C) none                # Sigma image name (made from data if blank or "none") 
    D) ./galfit/psfNOTOK.fits   #        # Input PSF image and (optional) diffusion kernel
    E) 1                   # PSF fine sampling factor relative to data 
    F) none                # Bad pixel mask (FITS image or ASCII coord list)
    G) none                # File with parameter constraints (ASCII file) 
    H) 1   150   1    150   # Image region to fit (xmin xmax ymin ymax)
    I) 100    100          # Size of the convolution box (x y)
    J) 26.563              # Magnitude photometric zeropoint 
    K) 0.14997825  0.14997825        # Plate scale (dx dy)    [arcsec per pixel]
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
     1) 75.0  75.0  1 1  #  position x, y
     3) 20.0    1          #  Integrated magnitude	
     4) 10.0      1          #  R_e (half-light radius)   [pix]
     5) 2.0      1          #  Sersic index n (de Vaucouleurs n=4) 
     6) 0.0000      0          #     ----- 
     7) 0.0000      0          #     ----- 
     8) 0.0000      0          #     ----- 
     9) 1.0      1          #  axis ratio (b/a)  
    10) 0.0    1          #  position angle (PA) [deg: Up=0, Left=90]
     Z) 0                      #  output option (0 = resid., 1 = Don't subtract) 
    
    # Object number: 2
     0) sky                    #  object type
     1) 1.3920      1          #  sky background at center of fitting region [ADUs]
     2) 0.0000      0          #  dsky/dx (sky gradient in x)
     3) 0.0000      0          #  dsky/dy (sky gradient in y)
     Z) 0                      #  output option (0 = resid., 1 = Don't subtract) 
    
    ================================================================================
    """.format(**vars())
    
    inputFile = open('/Users/rs548/Documents/Science/PeteHurley/galfit/' + filename + '.feedme', "w")
    inputFile.write(inputText)
    inputFile.close()
    
    log_file = open("a_log.txt", "a")
    
    #'-imax', '99', 
    print('Galfit.py is about to call galfit')
    print(inputText)
    subprocess.call(['/usr/local/bin/galfit', '/Users/rs548/Documents/Science/PeteHurley/galfit/' + filename + '.feedme'], 
                    stdout=log_file)
    
    log_file.close()
    
    #fitsOutput = fits.open('imgblock.fits')
    #fitsOutput.close()
    #chiSq = fitsOutput[2].header['CHISQ']
    #print('ChiSq is ',chiSq )
    