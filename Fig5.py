# -*- coding: utf-8 -*-

"""
Created on Fri May  5 15:44:15 2023

@author: AndrewTilman
"""
import numpy as np
import matplotlib.pyplot as plt
" MODEL PARAMETERS"
"Choose a ParamSet 1 2 3 or 4"
from ParamSet1 import params_dict
import functions as fn


"SAVE FIGS???"
saveFigs = False

"COLORS"
healthy = '#6DB33F'
infested = '#CFCA7A'
dead = '#8e7028'
opac = .8

"PLOT OPTIONS"
sMax=params_dict['sMax']
sMin = params_dict['sMin']

"OTHER PARAMS"

###! cost of treatment
c = params_dict['c']

###! survalence params 'h' 'i' and 'd' are the underlying states of the trees
pr_h = params_dict['pr_h'] #frac healthy
pr_i = params_dict['pr_i'] #frac infested
pr_d = params_dict['pr_d'] #frac dead

###! assessment params: pr_Hh is the prob of assessing a tree as 'H' healthy given that it is truly healthy 'h'
pr_Hh = params_dict['pr_Hh']
pr_Ih = params_dict['pr_Ih']
pr_Dh = params_dict['pr_Dh']

pr_Hi = params_dict['pr_Hi']
pr_Ii = params_dict['pr_Ii']
pr_Di = params_dict['pr_Di']

pr_Hd = params_dict['pr_Hd']
pr_Id = params_dict['pr_Id']
pr_Dd = params_dict['pr_Dd']

###! outcomes of treatment
hth = params_dict['hth'] # prob of a tree staying healty 'h', given that it is treated 't' and healthy 'h'
huh = params_dict['huh']  # prob of a tree staying healty 'h', given that it is untreated 'u' and healthy 'h'
hti = params_dict['hti'] # prob of a tree becoming healty 'h', given that it is treated 't' and infested 'i'
hui = params_dict['hui'] # prob of a tree becoming healty 'h', given that it is untreated 'u' and infested 'i'
htd = params_dict['htd']   # prob of saving a dead tree -> 0
hud = params_dict['hud']   # prob of saving a dead tree -> 0

###! Tree values
# V_o = value of a healthy tree to its owner 
# W_o = cost of a dead tree to its owner
# Delta_o = V_o + W_o = net benefit of saving a tree to its owner

a = params_dict['a'] # min value of Delta_o
b = params_dict['b'] # max value of Delta_o

V_m = params_dict['V_m'] # value of a healthy tree to society 
W_m = params_dict['W_m'] # Cost of a dead tree to society

vFrac = params_dict['vFrac'] # fraction of social benefit of saving a tree derived from ecosystem service benefit
wFrac = params_dict['wFrac'] # fraction of social beenfit of saving a dree derived from avoiding dealing with dead tree
delt_min = params_dict['delt_min'] # min social value of saving a tree 
delt_max = params_dict['delt_max'] # max social value of saving a tree

Delta_m = V_m + W_m # value of saving a tree to society
Delta_m_alt = params_dict['Delta_m_alt']

figName = params_dict['paramSet']

if vFrac + wFrac != 1:
    print('fix your value fractions! vFrac / wFrac so they add to one!')
if pr_h + pr_i +pr_d != 1:
    print('fix your survailence probs!')
if pr_Hh + pr_Ih + pr_Dh != 1 or pr_Hi + pr_Ii + pr_Di != 1 or pr_Hd + pr_Id + pr_Dd != 1 :
    print('fix your assessment accuracy probs so they sum to 1!')
if np.max((hth,huh,hti,hui,htd,hud))>1 or np.min((hth,huh,hti,hui,htd,hud))<0:
    print('sort out your treatment sucess probablities!')    
    
###! Bayes Theorem: pr_hH is the prob that a tree assessed as 'H' is truly healthy 'h'

pr_hH = fn.prob_xH('h',pr_h,pr_i,pr_d,pr_Hh,pr_Hi,pr_Hd)
pr_iH = fn.prob_xH('i',pr_h,pr_i,pr_d,pr_Hh,pr_Hi,pr_Hd)
pr_dH = fn.prob_xH('d',pr_h,pr_i,pr_d,pr_Hh,pr_Hi,pr_Hd)


pr_hI = fn.prob_xI('h',pr_h,pr_i,pr_d,pr_Ih,pr_Ii,pr_Id)
pr_iI = fn.prob_xI('i',pr_h,pr_i,pr_d,pr_Ih,pr_Ii,pr_Id)
pr_dI = fn.prob_xI('d',pr_h,pr_i,pr_d,pr_Ih,pr_Ii,pr_Id)


pr_hD = fn.prob_xD('h',pr_h,pr_i,pr_d,pr_Dh,pr_Di,pr_Dd)
pr_iD = fn.prob_xD('i',pr_h,pr_i,pr_d,pr_Dh,pr_Di,pr_Dd)
pr_dD = fn.prob_xD('d',pr_h,pr_i,pr_d,pr_Dh,pr_Di,pr_Dd)

###! fraction of trees with each classification
pr_H = pr_Hh * pr_h + pr_Hi * pr_i + pr_Hd * pr_d
pr_I = pr_Ih * pr_h + pr_Ii * pr_i + pr_Id * pr_d
pr_D = pr_Dh * pr_h + pr_Di * pr_i + pr_Dd * pr_d

###! marginal impact of treatment
t_h = hth - huh
t_i = hti - hui
t_d = htd - hud


###! exected impact of treatemt on saving a tree, given its assessed state
k_H = fn.k_X(t_h,t_i,t_d,pr_hH,pr_iH,pr_dH)
k_I = fn.k_X(t_h,t_i,t_d,pr_hI,pr_iI,pr_dI)
k_D = fn.k_X(t_h,t_i,t_d,pr_hD,pr_iD,pr_dD)


###! Plots


Delta_range = np.linspace(delt_min, delt_max, 500)
sH = np.zeros(len(Delta_range))
sI = np.zeros(len(Delta_range)) 
sD = np.zeros(len(Delta_range))
pr_tH = np.zeros(len(Delta_range))
pr_tI = np.zeros(len(Delta_range))
pr_tD = np.zeros(len(Delta_range))
EUmH = np.zeros(len(Delta_range))
EUmI = np.zeros(len(Delta_range))
EUmD = np.zeros(len(Delta_range))
pr_survH = np.zeros(len(Delta_range))
pr_survI = np.zeros(len(Delta_range))
pr_survD = np.zeros(len(Delta_range))
pr_survTot = np.zeros(len(Delta_range))

for i in range( len(Delta_range)):
    sH[i] = fn.s_opt(c,k_H,Delta_range[i],a,b)
    sI[i] = fn.s_opt(c,k_I,Delta_range[i],a,b)
    sD[i] = fn.s_opt(c,k_D,Delta_range[i],a,b)
    pr_tH[i] = fn.pr_treated(c,k_H,Delta_range[i],a,b)
    pr_tI[i] = fn.pr_treated(c,k_I,Delta_range[i],a,b)
    pr_tD[i] = fn.pr_treated(c,k_D,Delta_range[i],a,b)
    EUmH[i] = fn.EU_m(c,k_H,Delta_range[i],a,b,fn.Pi(vFrac*Delta_range[i],wFrac*Delta_range[i],huh,hui,pr_hH,pr_iH,pr_dH))
    EUmI[i] = fn.EU_m(c,k_I,Delta_range[i],a,b,fn.Pi(vFrac*Delta_range[i],wFrac*Delta_range[i],huh,hui,pr_hI,pr_iI,pr_dI))
    EUmD[i] = fn.EU_m(c,k_D,Delta_range[i],a,b,fn.Pi(vFrac*Delta_range[i],wFrac*Delta_range[i],huh,hui,pr_hD,pr_iD,pr_dD))
    pr_survH[i] = fn.pr_survX(fn.pr_treated(c,k_H,Delta_range[i],a,b),pr_hH,pr_iH,pr_dH,hth,huh,hti,hui,htd,hud)
    pr_survI[i] = fn.pr_survX(fn.pr_treated(c,k_I,Delta_range[i],a,b),pr_hI,pr_iI,pr_dI,hth,huh,hti,hui,htd,hud)
    pr_survD[i] = fn.pr_survX(fn.pr_treated(c,k_D,Delta_range[i],a,b),pr_hD,pr_iD,pr_dD,hth,huh,hti,hui,htd,hud)
    pr_survTot[i] =  pr_survH[i]*pr_H + pr_survI[i]*pr_I + pr_survD[i]*pr_D
    
fig1 = plt.figure(figsize=(8,2.75))
ax1 = fig1.add_subplot(131)
ax2 = fig1.add_subplot(132)
ax4 = fig1.add_subplot(133)


plt.tight_layout()
ax1.plot(Delta_range,sH, color=healthy,linewidth=3,alpha=opac,label='$\hat h$')
ax1.plot(Delta_range,sI, color=infested,linewidth=3,alpha=opac,label='$\hat i$')
ax1.plot(Delta_range,sD, color=dead,linewidth=3,alpha=opac,label='$\hat d$')

ax2.plot(Delta_range,pr_tH, color=healthy,linewidth=3,alpha=opac,label='$\hat h$')
ax2.plot(Delta_range,pr_tI, color=infested,linewidth=3,alpha=opac,label='$\hat i$')
ax2.plot(Delta_range,pr_tD, color=dead,linewidth=3,alpha=opac,label='$\hat d$')

ax4.plot(Delta_range,pr_survH, color=healthy,linewidth=3,alpha=opac,label='$\hat h$')
ax4.plot(Delta_range,pr_survI, color=infested,linewidth=3,alpha=opac,label='$\hat i$')
ax4.plot(Delta_range,pr_survD, color=dead,linewidth=3,alpha=opac,label='$\hat d$')


rows=['$a\mapsto \Delta_o$min','$b\mapsto\Delta_o$max','$c$',
      '$P(h)$','$P(i)$','$P(d)$',
      '$P(\hat h\mid h)$','$P(\hat i\mid h)$','$P(\hat d\mid h)$',
      '$P(\hat h\mid i)$','$P(\hat i\mid i)$','$P(\hat d\mid i)$',
      '$P(\hat h\mid d)$','$P(\hat i\mid d)$','$P(\hat d\mid d)$',
      '$h_{th}$','$h_{uh}$',
      '$h_{ti}$','$h_{ui}$',
      '$h_{td}$','$h_{ud}$',
      '$V_m$','$W_m$']
columns = ['value']
cell_text = [[a],[b],[c],
             [pr_h],[pr_i],[pr_d],
             [pr_Hh],[pr_Ih],[pr_Dh],
             [pr_Hi],[pr_Ii],[pr_Di],
             [pr_Hd],[pr_Id],[pr_Dd],
             [hth],[huh],
             [hti],[hui],
             [htd],[hud],
             [str(vFrac)+' $*\Delta_m$'],[str(wFrac)+' $*\Delta_m$']]


ax1.set_ylabel('Optimal Subsidy, $s^*_{\hat\phi}$')
ax1.set_xlabel('$\Delta_m$');
ax1.legend(loc='upper left')
ax2.set_ylabel('Treatment Probability')
ax2.set_xlabel('$\Delta_m$');
ax4.set_ylabel('Survival Probability')
ax4.set_xlabel('$\Delta_m$');
ax1.set_ylim([sMin, sMax])
ax1.set_xlim([delt_min, delt_max])
ax2.set_ylim([-.01, 1.05])
ax2.set_xlim([delt_min, delt_max])
ax4.set_ylim([-.01, 1.05])
ax4.set_xlim([delt_min, delt_max])
plt.tight_layout()



"SAVE FIGS"
if saveFigs == True:
    fig1.savefig("FIGS/optSub_"+figName+".pdf",bbox_inches='tight', dpi=150)
    fig1.savefig("FIGS/optSub_"+figName+".png", bbox_inches='tight',dpi=250)