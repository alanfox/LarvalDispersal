# -*- coding: utf-8 -*-
"""
Created on Thu May 07 10:42:26 2015

@author: af26
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as syop

x = -0.5


depth = 200.0

def f_builder_1(x,hij):
    N = 40
    B = 0.05
    theta = 8.0
    hc = 150.0
    def f2(Sk):
        h = (hij - hc)/hij
        Ck1 = (1.0 - B) * np.sinh(theta * Sk) / np.sinh(theta)
        Ck2 = (B * (np.tanh(theta * (Sk + 0.5)) - np.tanh(theta * 0.5)) / 
                   (2 * np.tanh(0.5* theta)))
                   
        Ck = Ck1 + Ck2
        return (Sk + h * (Ck - Sk)) - x
    return f2

def f_builder_2(hij):
    N = 40
    B = 0.05
    theta = 8.0
    hc = 150.0
    def f2(Sk):
        h = (hij - hc)/hij
        Ck1 = (1.0 - B) * np.sinh(theta * Sk) / np.sinh(theta)
        Ck2 = (B * (np.tanh(theta * (Sk + 0.5)) - np.tanh(theta * 0.5)) / 
                   (2 * np.tanh(0.5* theta)))
                   
        Ck = Ck1 + Ck2
        return (Sk + h * (Ck - Sk))
    return f2
    
f1 = f_builder_1(x, depth) 
f2 = f_builder_2(depth)

z = syop.brentq(f1,-1.1,0.1)

print x, z, f2(z)

print -40*z, int((-40 * z)//1), (-40 * z)% 1
        
   
