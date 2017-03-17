#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 15 16:40:10 2016

REGPHOT: Resolved Galactic PHOTometry

This code is designed for investigating models of galaxy profiles in a very 
general/robust way for images of highly blended resolved galaxies.

It will calculate the Bayesian evidence for a model with a given number of 
Sersic/Gaussian profiles in order to determine how many to use and will
then provide a maximum liklihood (chi squared minimum) output parameters
for the favoured model.

Dr Raphael Shirley
Daphne Jackson Research Fellow
Astronomy Centre
University of Sussex
www.raphaelshirley.co.uk
@raphaelshirley
raphael.shirley@sussex.ac.uk

@author: rs548
"""



from setuptools import setup

setup(
      name='regphot',
      version='1.0',
      description='Resolved Galactic PHOTometry',
      author='Raphael Shirley',
      author_email='raphael.shirley@sussex.ac.uk',
      url='https://github.com/raphaelshirley/regphot',
      classifiers=[
          "License :: MIT Public License",
          "Operating System :: OS Independent",
          "Programming Language :: Python :: 2.7",
          "Topic :: Scientific/Engineering :: Astronomy"
      ],
      
)
