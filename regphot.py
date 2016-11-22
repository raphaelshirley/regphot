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


log_file = open("a_log.txt", "a")

#'-imax', '99', 
subprocess.call(['/usr/local/bin/galfit', 'galfitSDSS.feedme'], 
                stdout=log_file)

log_file.close()