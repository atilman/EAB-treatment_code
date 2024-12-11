# -*- coding: utf-8 -*-
"""
Created on Fri May  5 15:44:15 2023

@author: AndrewTilman
"""
import numpy as np
import matplotlib.pyplot as plt
" MODEL PARAMETERS"
"Choose a ParamSet to import"
from ParamSet_high import params_dict
import functions as fn

"SAVE FIGS???"
saveFigs = True

"SIMTYPE"
#! SimType defines the Subsidy and treatment probs for simulation
    #! 1: no private treatment or public treatment
    #! 2: no private treatment but optimal public treatment
    #! 3: no private subsidies and no public treatment
    #! 4: no private subsidies but optimal public treatment
    #! 5: optimal private subsidies but no public treatment
    #! 6: optimal private subsidies and optimal public treatment
simType = 6
params_dict["simType"] = simType
#COLORS
healthy = '#6DB33F'
infested = '#CFCA7A'
dead = '#8e7028'
opac = .75
linWidPub = 3
linWidPriv = 4
linWid = 4

t_max = 50
stepsize = .05
# timeGap=3.5
# sims= 8
privFrac = 6/10
pubFrac = 1- privFrac
y0 = [99/100,1/100,0]
y0_mix = [y0[0]*privFrac,y0[1]*privFrac,y0[2]*privFrac,y0[0]*pubFrac,y0[1]*pubFrac,y0[2]*pubFrac]
#y0_mix = [0,0,0,75/100,15/100,10/100]
beta = params_dict['beta']
gamma = params_dict['gamma']
alpha = params_dict['alpha']

t_eval = np.arange(0, t_max, stepsize)
args = (beta,gamma,alpha,params_dict)


"PLOT OPTIONS"
uMax = params_dict['uMax']
uMin=params_dict['uMin']
sMax=params_dict['sMax']
sMin = params_dict['sMin']
sMax2=params_dict['sMax2']
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



t_horiz = params_dict['t_horiz']
eps_h = params_dict['eps_h']
eps_i = params_dict['eps_i']
spilloverCare = params_dict['spill_weight']


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
# if np.max((hth,huh,hti,hui,htd,hud))>1 or np.min((hth,huh,hti,hui,htd,hud))<0:
#     print('sort out your treatment sucess probablities!')    
    
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

###! marginal impact of treatment on focal tree survival
muh = fn.m_uh(t_horiz,pr_i,beta,gamma)
mth = fn.m_th(t_horiz,pr_i,beta,gamma,eps_h,eps_i)
mui = fn.m_ui(t_horiz,gamma)
mti = fn.m_ti(t_horiz,gamma,eps_i)
mud = 1
mtd = 1

t_h = muh - mth 
t_i = mui - mti
t_d = mud - mtd

###! marginal impact of treatment on focal tree surrounding tree survival
smuh = fn.sm_uh(t_horiz,pr_i,pr_h,beta,gamma)
smth = fn.sm_th(t_horiz,pr_i,pr_h,beta,gamma,alpha,eps_h,eps_i)
smui = fn.sm_ui(pr_h,beta,gamma)
smti = fn.sm_ti(pr_h,beta,gamma,alpha,eps_i)
smud = 0
smtd = 0

s_h = (smuh - smth) * spilloverCare
s_i = (smui - smti) * spilloverCare
s_d = (smud - smtd) * spilloverCare


###! exected impact of treatemt on focal tree survival, given its assessed state
k_H = fn.k_X(t_h,t_i,t_d,pr_hH,pr_iH,pr_dH)
k_I = fn.k_X(t_h,t_i,t_d,pr_hI,pr_iI,pr_dI)
k_D = fn.k_X(t_h,t_i,t_d,pr_hD,pr_iD,pr_dD)

###! exected impact of treatemt of a focal tree on community tree survival, given its assessed state
L_H = fn.L_X(s_h,s_i,s_d,pr_hH,pr_iH,pr_dH)
L_I = fn.L_X(s_h,s_i,s_d,pr_hI,pr_iI,pr_dI)
L_D = fn.L_X(s_h,s_i,s_d,pr_hD,pr_iD,pr_dD)


###! Plots
opac = .8

Delta_range = np.linspace(delt_min, delt_max, 500)
sH = np.zeros(len(Delta_range)) # Subsidy level for H aassessed tree
sI = np.zeros(len(Delta_range)) # Subsidy level for I aassessed tree
sD = np.zeros(len(Delta_range)) # Subsidy level for D aassessed tree
pr_tH = np.zeros(len(Delta_range)) # Treatment probability for H assessed tree
pr_tI = np.zeros(len(Delta_range)) # Treatment probability for I assessed tree
pr_tD = np.zeros(len(Delta_range)) # Treatment probability for D assessed tree

pr_survH = np.zeros(len(Delta_range))
pr_survI = np.zeros(len(Delta_range))
pr_survD = np.zeros(len(Delta_range))


for i in range( len(Delta_range)):
    sH[i] = fn.s_opt(c,k_H,L_H,Delta_range[i],a,b)
    sI[i] = fn.s_opt(c,k_I,L_I,Delta_range[i],a,b)
    sD[i] = fn.s_opt(c,k_D,L_D,Delta_range[i],a,b)
    pr_tH[i] = fn.pr_treated(c,k_H,L_H,Delta_range[i],a,b)
    pr_tI[i] = fn.pr_treated(c,k_I,L_I,Delta_range[i],a,b)
    pr_tD[i] = fn.pr_treated(c,k_D,L_D,Delta_range[i],a,b)

    pr_survH[i] = fn.pr_survX(fn.pr_treated(c,k_H,L_H,Delta_range[i],a,b),pr_hH,pr_iH,pr_dH,1-mth,1-muh,1-mti,1-mui,1-mtd,1-mud)
    pr_survI[i] = fn.pr_survX(fn.pr_treated(c,k_I,L_I,Delta_range[i],a,b),pr_hI,pr_iI,pr_dI,1-mth,1-muh,1-mti,1-mui,1-mtd,1-mud)
    pr_survD[i] = fn.pr_survX(fn.pr_treated(c,k_D,L_D,Delta_range[i],a,b),pr_hD,pr_iD,pr_dD,1-mth,1-muh,1-mti,1-mui,1-mtd,1-mud)
    
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


ax1.set_ylabel('Optimal Subsidy, $s^*_{\hat\phi}$')
ax1.set_xlabel('$\Delta_m$');
ax2.legend(loc='lower right')
#ax3.legend(loc='lower right',bbox_to_anchor=(1.34, .29))
ax2.set_ylabel('Treatment Probability')
ax2.set_xlabel('$\Delta_m$');
# ax3.set_ylabel('Expected muni utility, $E[U^*_{m}]$')
# ax3.set_xlabel('$\Delta_m$');
# ax3.set_yticks([0])
ax4.set_ylabel(str(t_horiz)+'-year Survival Prob.')
ax4.set_xlabel('$\Delta_m$');
ax1.set_ylim([sMin, sMax])
ax1.set_xlim([delt_min, delt_max])
# ax3.set_ylim([uMin, uMax])
# ax3.set_xlim([delt_min, delt_max])
ax2.set_ylim([-.01, 1.05])
ax2.set_xlim([delt_min, delt_max])
ax4.set_ylim([-.01, 1.05])
ax4.set_xlim([delt_min, delt_max])

ax1.text(-0.05, 1.08, '(a)', transform=ax1.transAxes, va='top', fontsize=9, fontweight='bold')
ax2.text(-0.05, 1.08, '(b)', transform=ax2.transAxes, va='top', fontsize=9, fontweight='bold')
ax4.text(-0.05, 1.08, '(c)', transform=ax4.transAxes, va='top', fontsize=9, fontweight='bold')



plt.tight_layout()



"SAVE FIGS"
if saveFigs == True:
    fig1.savefig("FIGS/Fig2_"+figName+".png", bbox_inches='tight',dpi=250)


