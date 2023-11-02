# -*- coding: utf-8 -*-
"""
Created on Thu Feb  2 15:16:21 2023

@author: AndrewTilman
"""

" MODEL PARAMETERS"
paramSet = '1'
###! cost of treatment
c = 300

###! survalence params 'h' 'i' and 'd' are the underlying states of the trees
pr_h = 7/10 #frac healthy
pr_i = 1/4 #frac infested
pr_d = 1/20 #frac dead

###! assessment params: pr_Hh is the prob of assessing a tree as 'H' healthy given that it is truly healthy 'h'
pr_Hh = 89/100
pr_Ih = 10/100
pr_Dh = 1/100

pr_Hi = 49/100
pr_Ii = 50/100
pr_Di = 1/100

pr_Hd = 1/100
pr_Id = 19/100
pr_Dd = 80/100

###! outcomes of treatment
hth = 1 # prob of a tree staying healty 'h', given that it is treated 't' and healthy 'h'
huh = 7/10  # prob of a tree staying healty 'h', given that it is untreated 'u' and healthy 'h'
hti = 9/10 # prob of a tree becoming healty 'h', given that it is treated 't' and infested 'i'
hui = 0 # prob of a tree becoming healty 'h', given that it is untreated 'u' and infested 'i'
htd = 0   # prob of saving a dead tree -> 0
hud = 0   # prob of saving a dead tree -> 0

###! Tree values
# V_o = value of a healthy tree to its owner 
# W_o = cost of a dead tree to its owner
# Delta_o = V_o + W_o = net benefit of saving a tree to its owner

a = 600 # min value of Delta_o
b = 900 # max value of Delta_o

V_m = 500 # value of a healthy tree to society 
W_m = 50 # Cost of a dead tree to society 
W_m_alt = 500 #cost of a dead public tree to society, including removal cost
Delta_m = V_m + W_m # value of saving a private tree to society
Delta_m_alt = V_m + W_m_alt # value of saving a publictree to society
vFrac = 7/10
wFrac = 3/10
delt_min = 0
delt_max = 1000

#PLOT options
uMax = 600
uMin=-150
sMax=150
sMin = -2
sMax2=125

params_dict = {'c':c,
        'a':a,
        'b': b, 
        'pr_h': pr_h, 
        'pr_i': pr_i,
        'pr_d': pr_d,
        'pr_Hh': pr_Hh,
        'pr_Ih':pr_Ih,
        'pr_Dh': pr_Dh,
        'pr_Hi': pr_Hi,
        'pr_Ii': pr_Ii,
        'pr_Di':pr_Di,
        'pr_Hd':pr_Hd,
        'pr_Id':pr_Id,
        'pr_Dd':pr_Dd,
        'hth':hth,
        'huh': huh,
        'hti': hti,
        'hui': hui,
        'htd':htd,
        'hud':hud,
        'V_m':V_m,
        'W_m':W_m,
        'W_m_alt':W_m_alt,
        'Delta_m':Delta_m,
        'Delta_m_alt':Delta_m_alt,
        'vFrac':vFrac,
        'wFrac':wFrac,
        'delt_min':delt_min,
        'delt_max':delt_max,
        'paramSet':paramSet,
        'uMax':uMax,
        'uMin':uMin,
        'sMax':sMax,
        'sMax2':sMax2,
        'sMin':sMin
         }