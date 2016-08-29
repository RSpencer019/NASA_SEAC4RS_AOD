# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 16:21:24 2016

@author: rsspenc3
"""

import os

for folder in os.listdir(os.getcwd()):
    n = 0    
    for i in range(len(list(names))):             
        if list(names)[i] in folder:
            n += 1
    if n != 1:
        print ('Removed ', folder)
        os.remove(folder)