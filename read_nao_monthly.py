# -*- coding: utf-8 -*-
"""
Created on Mon Nov 02 11:07:02 2015

@author: af26
"""

def read_nao_monthly(filename):
    year = []
    month = []
    nao = []
    for line in infile:
        wordlist = line.split()
        year.append(int(wordlist[0]))
        month.append(int(wordlist[1]))
        nao.append(float(wordlist[2]))
    return year, month, nao
    
filename = ('C:/Users/af26/NAO/nao_monthly.txt')
infile = open(filename,'r')

year, month, nao = read_nao_monthly(infile)

print [nao[i] for i in range(len(nao)) if month[i] == 4 and year[i] == 1970]