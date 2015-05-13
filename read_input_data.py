# -*- coding: utf-8 -*-
"""
Created on Sun May 03 15:31:37 2015

@author: af26
"""

input_data_file = open('input.dat', 'r')

input_dict = {}

for line in input_data_file:
    wordlist = line.split()
    input_dict[wordlist[0]] = wordlist[-1]
    
print input_dict
    
