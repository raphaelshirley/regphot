#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 17 15:49:33 2017

I have a long term plan to experiment with using an RDF database for galaxy 
catalogues.

Making a distinction between objects and sources; each object will link to any 
sources associated with it from different instruments.

thus the job of the code is to maintain an RDF database of objects updated based
on catalogues of sources from different surveys.

A sample query could for instance return all objects extracted from both SDSS
and UltraVISTA.

@author: rs548
"""

from rdflib import Namespace, Literal, Graph, XSD
from astropy.io import fits



typesdict ={'L':, # Logical 
            'X':, # bit
            'B':XSD.byte, # Unsigned byte
            'I':, # 16-bit int
            'J':XSD.int, # 32-bit int
            'K':XSD.long, # 64-bit int
            'A':, # character
            'E':, # single prec float
            'D':, # double prec float
            'C':, # single prec compl
            'M':, # double precis complex
            'P':, # array descriptor
            'Q':} # arraydescriptor

def catalogue2rdf(cataloguefits, folder):
    
    n = Namespace('http://raphaelshirley.co.uk/galaxies/')
    cat = fits.open(cataloguefits)
    

    columns = cat[1].data.names
    predicates = columns
    datatypes = cat[1].data.formats
    xsdtypes = datatypes
    for n in range(0,len(xsdtypes)):
        xsdtypes[n] = typesdict[str(xsdtypes[n])]
               
    for source in cat:
        
        sourcenode = BNode()
        sourcename=source[0] #this should be a generated URI that we can then maintain
        rdffile = open(folder + sourcename + '.rdf', 'w')
        g= Graph()
        for datum in range(0,len(source)):
            ]
            g.add((n.sourcename ,predicates[datum] ,Literal(source[datum], datatype=xsdtypes[datum])))
        print(g.serialise(format='turtle'))
        rdffile.write(g.serialize(format='xml'))
        rdffile.close()
            
if __name__ == '__main__':
    catalogue2rdf('/Users/rs548/Documents/Science/PeteHurley/cosmos-hyperleda-sdss-dr13-2massXSC-4reg-flagged_HELP.fits', 
                  '/Users/rs548/Documents/Science/PeteHurley/rdf/')
    
