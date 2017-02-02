#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  1 16:10:07 2017

This should define a class Model which is an internal representation of a model.

It should allow arbitrary definition of functions and should be able to call 
chisquared from either galfit or profit.


@author: rs548
"""

class Model:
    
    def __init__(self, name, image):
        self.name = name
        self.components = []
        
        
    def addComponent(self, component):    #Dictionaries?
        self.components = self.components + [component]
        
        
    def getGalfitChiSq(self):
        print('send components to galfit module to produce feedme and return chisq')
        
    def getProFitChiSq(self):
        print('send components to ProFit return chisq')