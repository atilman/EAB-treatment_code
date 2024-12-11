# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 14:12:55 2024

@author: AndrewTilman
"""


" MODEL PARAMETERS"
paramSet = 'midNoAssess'

#DYNAMIC SIM params
beta = 1       # infest rate
gamma = 3/10   # death rate
alpha = 1      # recovery rate

###! cost of treatment
c = 250

###! survalence params 'h' 'i' and 'd' are the underlying states of the trees
pr_h = 16/20 #frac healthy
pr_i = 3/20 #frac infested
pr_d = 1/20 #frac dead

###! assessment params: pr_Hh is the prob of assessing a tree as 'H' healthy given that it is truly healthy 'h'
pr_Hh = 1/3
pr_Ih = 1/3
pr_Dh = 1/3

pr_Hi = 1/3
pr_Ii = 1/3
pr_Di = 1/3

pr_Hd = 1/3
pr_Id = 1/3
pr_Dd = 1/3

###! Expected outcomes of treatment
t_horiz = 3
eps_h = .97
eps_i = .5
spill_weight = 1

# htd = 0   # prob of saving a dead tree -> 0
# hud = 0   # prob of saving a dead tree -> 0
# hui = 0 # prob of a tree becoming healty 'h', given that it is untreated 'u' and infested 'i'
# huh = hui + (1-hui)*np.exp(-pr_i*beta*t_horiz)  # prob of a tree staying healty 'h', given that it is untreated 'u' and healthy 'h'
# hth = 1 - tau_h*tau_i*(1-huh) # prob of a tree staying healty 'h', given that it is treated 't' and healthy 'h'
# hti = 1-tau_i*(1-hui) # prob of a tree becoming healty 'h', given that it is treated 't' and infested 'i'



###! Tree values
# V_o = value of a healthy tree to its owner 
# W_o = cost of a dead tree to its owner
# Delta_o = V_o + W_o = net benefit of saving a tree to its owner

a = 675 # min value of Delta_o
b = 1100 # max value of Delta_o

V_m = 1000 # value of a healthy tree to society 
W_m = 150 # Cost of a dead tree to society 
W_m_alt = 850 #cost of a dead public tree to society, including removal cost
Delta_m = V_m + W_m # value of saving a private tree to society
Delta_m_alt = V_m + W_m_alt # value of saving a publictree to society
vFrac = 5/10
wFrac = 5/10
delt_min = 0
delt_max = 1000

#PLOT options
uMax = 600
uMin=-150
sMax=250
sMin = -2
sMax2 = 300




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
        # 'hui': hui,
        # 'htd':htd,
        # 'hud':hud,
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
        'sMin':sMin,
        'beta':beta,
        'gamma':gamma,
        'alpha':alpha,
        't_horiz': t_horiz,
        'eps_h': eps_h,
        'eps_i': eps_i,
        'spill_weight':spill_weight
         }