# -*- coding: utf-8 -*-
"""
REGPHOT python software package to investigate parameter space for multi-galaxy
profiles on densely populated highly resolved images with foreground stars

Written by Dr Raphael Shirley
Astronomy Centre
University of Sussex
raphael.shirley@sussex.ac.uk
"""

import subprocess
#import astropy
from astropy.io import fits


method = 2
inputText = """
===============================================================================
# IMAGE and GALFIT CONTROL PARAMETERS
A) frame-g-004858-1-0480.fits            # Input data image (FITS file)
B) sdss.fits       # Output data image block
C) none                # Sigma image name (made from data if blank or "none") 
D) psf.fits   #        # Input PSF image and (optional) diffusion kernel
E) 1                   # PSF fine sampling factor relative to data 
F) none                # Bad pixel mask (FITS image or ASCII coord list)
G) none                # File with parameter constraints (ASCII file) 
H) 1602    2043   992    1489   # Image region to fit (xmin xmax ymin ymax)
I) 100    100          # Size of the convolution box (x y)
J) 26.563              # Magnitude photometric zeropoint 
K) 0.038  0.038        # Plate scale (dx dy)    [arcsec per pixel]
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
 1) 1770.93  1184.06  1 1  #  position x, y
 3) 28.68     1          #  Integrated magnitude	
 4) 230.1160      1          #  R_e (half-light radius)   [pix]
 5) 1.90      1          #  Sersic index n (de Vaucouleurs n=4) 
 6) 0.0000      0          #     ----- 
 7) 0.0000      0          #     ----- 
 8) 0.0000      0          #     ----- 
 9) 0.68      1          #  axis ratio (b/a)  
10) -13.3690    1          #  position angle (PA) [deg: Up=0, Left=90]
 Z) 0                      #  output option (0 = resid., 1 = Don't subtract) 

# Object number: 2
 0) sky                    #  object type
 1) 1.3920      1          #  sky background at center of fitting region [ADUs]
 2) 0.0000      0          #  dsky/dx (sky gradient in x)
 3) 0.0000      0          #  dsky/dy (sky gradient in y)
 Z) 0                      #  output option (0 = resid., 1 = Don't subtract) 

================================================================================
""".format(**vars())

inputFile = open("thisRun.feedme", "w")
inputFile.write(inputText)
inputFile.close()

log_file = open("a_log.txt", "a")

#'-imax', '99', 
subprocess.call(['/usr/local/bin/galfit', 'galfitSDSS.feedme'], 
                stdout=log_file)

log_file.close()

fitsOutput = fits.open('imgblock.fits')
fitsOutput.close()
chiSq = fitsOutput[2].header['CHISQ']
print('ChiSq is ',chiSq )