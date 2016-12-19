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
    print('galfitm.py has been called')
    filename=str(fitsfile)
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
    A) gal-r.fits,gal-g.fits,gal-i.fits
    
    # Band labels (CSL of <nbands> labels containing no whitespace)
    # (these must be unique in a case-insensitive manner)
    # (can be omitted if fitting a single band)
    A1) r,g,i
    
    # Band wavelengths (CSL of values)
    # (choice of wavelength units is arbitrary, as long as consistent,
    #  but affects the resulting wavelength-dependence parameters)
    A2) 6220,4750,7630
    
    # Output data image block (FITS filename)
    B) /Users/rs548/Documents/Science/PeteHurley/galfitm/{filename}-output.fits
    
    # Sigma image name (CSL of <nbands> FITS filenames or "none")
    # (if an individual filename is specified as "none", then that sigma
    #  image will be made from data; if the whole entry consists of just a
    #  single "none", then all sigma images will be made from data.)
    # One can also add a minimum sigma value, such that any galfit-created
    # sigma image will have a minimum of that value times the sky-subtracted
    # input data.
    C) none
    #C) sig-r.fits,sig-g.fits,sig-i.fits
    #C) none,none,none
    #C) sig-r.fits,none,sig-i.fits         # perhaps unwise to do this in practice
    #C) none    0.1			      # min. sigma is 10% of the data flux
    
    # Input PSF image (CSL of <nbands> FITS filenames) 
    # and a single diffusion kernel (FITS filename, # or omitted)
    D) psf-r.fits,psf-g.fits,psf-i.fits  #  
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
    G) none              
    #G) constraints_filename
    
    # Image region to fit (xmin xmax ymin ymax)
    H) 1    93   1    93
    
    # Size of the convolution box (x y)
    I) 100    100
    
    # Magnitude photometric zeropoint (CSL of <nbands> values)
    J) 26.563,25.123,27.987
    
    # Plate scale (dx dy)   [arcsec per pixel]
    K) 0.038  0.038
    
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
    
    
    # THE OBJECT LIST BELOW can be however long or short as the complexity
    # requires.  The user has complete freedom to mix and match the components
    # by duplicating each object block.
    
    # INITIAL FITTING PARAMETERS
    #
    # column 1: Parameter number
    #
    # column 2:
    # -- Parameter 0: the allowed functions are: sersic, nuker, expdisk
    #    	       	  edgedisk, devauc, king, moffat, gaussian, ferrer, psf, sky
    # -- Parameter 1-10: value of the initial parameters
    # -- and Parameter C0: For diskiness/boxiness (<0 = disky, >0 = boxy)
    #  	     	  By default this is a CSL of the values of the parameter
    #  	     	  in each of the input bands.  This may be optionally indicated
    #		  by putting the word 'band' at the end of the line (before any
    #  	     	  comment).
    #  	     	  Only specifying a single value with multiple input bands
    #  	     	  assumes the same value for all bands.
    #                 This can also optionally be specified in terms of Chebyshev
    #	       	  coefficients by adding the word 'cheb' at the end of the line
    #  	     	  (before any comment).
    #	       	  In this case one should give a CSL of at most <nbands> values
    #		  corresponding to coefficients of a Chebyshev series.
    # 	     	  First value of the CSL specifies the parameter value at the
    #	     	  average wavelength of the input bands.
    #	     	  Additional values in the CSL specify the variation in that
    #	     	  parameter value with wavelength, from linear (1st-order),
    #	     	  quadratic (2nd-order), up to <nbands>-order (which should be
    #	     	  equivalent to fitting the value independently for each band
    #		  Values omitted from the end of the CSL are assumed to be zero.
    # -- Parameter Z: Outputting image options, the options are:
    #              	  0 = normal, i.e. subtract final model from the data to create
    #		      the residual image
    #	      	  1 = Leave in the model -- do not subtract from the data
    #
    # column 3: This may be specified in one of two ways:
    #           An integer giving the order of the Chebyshev series , e.g.,
    #             0 = fixed to input value(s)
    #             1 = fit a constant offset from the input value(s)
    #             2 = fit a linear function of wavelength
    #             3 = fit a quadratic function of wavelength, etc.
    #           Note that, for >2, the input values are fit by a polynomial function
    #           of the specified order before fitting begins.
    #           Alternatively, one may give a CSL of at most <nbands> integers
    #           indicating whether or not that coefficient is allowed to vary
    #           (yes = 1, no = 0).  Values omitted from the end of the CSL are
    #           assumed to be zero.
    #
    # column 4: comment
    
    # Sersic function --------------------------------------------------------------
    
    # Only this first function includes multi-band examples, but the same approach
    # should work for all these functions.
    
     0) sersic     # Object type
     1) 300.  1    # position x [pixel]  (constant with wavelength)
     2) 350.  1    # position y [pixel]
     3) 20.0,21.0,22.0  3     # total magnitude in each band
     4) 4.30,4.40,4.5   2     # R_e in each band
     5) 5.20,5.20,5.20  1     # Sersic exponent in each band
     9) 0.30,0.30,0.30  1     # axis ratio (b/a) in each band
    10) 10.0            1     # position angle (PA), same value in each band
     Z) 0                  #  Skip this model in output image?  (yes=1, no=0)
    
     0) sersic     # Object type
     1) 300.  1    # position x [pixel]  (constant with wavelength)
     2) 350.  1    # position y [pixel]
     3) 20.00,0.2,0.0  1,1,1 cheb # total magnitude and coeffs of wavelength dependence
     4) 4.30,0.1,0.0   1,1,0 cheb # R_e and coeffs of wavelength dependence
     5) 5.20,0.0,0.0   1,0,0 cheb # Sersic exponent and coeffs of wavelength dependence
     9) 0.30,0.0,0.0   1     cheb # axis ratio (b/a) and coeffs of wavelength dependence
    10) 10.0           1     cheb # position angle (PA) and coeffs of wavelength dependence
     Z) 0                  #  Skip this model in output image?  (yes=1, no=0)
    
    
  
    # PSF fit ----------------------------------------------------------------------
    
     0) psf                # object type
     1) 300.       1       # position x [pixel]
     2) 357.4      1       # position y [pixel]
     3) 18.5       1       # total magnitude     
     Z) 0                  #  Skip this model in output image?  (yes=1, no=0)
    
    # sky --------------------------------------------------------------------------
    
     0) sky
     1) 0.77       0       # sky background       [ADU counts]
     2) 0.000      0       # dsky/dx (sky gradient in x) 
     3) 0.000      0       # dsky/dy (sky gradient in y) 
     Z) 0                  #  Skip this model in output image?  (yes=1, no=0)
    
    
    
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
    