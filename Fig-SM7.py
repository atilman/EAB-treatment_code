# -*- coding: utf-8 -*-
"""
Created on 11 / 5 / 2024

@author: AndrewTilman
"""
import numpy as np
import matplotlib.pyplot as plt
" MODEL PARAMETERS"

"Choose a ParamSets to import"
from ParamSet_middle_noAssess import params_dict as pd0
from ParamSet_middle import params_dict as pd1
from ParamSet_middle_perfect import params_dict as pd2

import functions as fn
from scipy.integrate import solve_ivp

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

"COLORS"
healthy = '#6DB33F'
infested = '#CFCA7A'
dead = '#8e7028'
opac = .75
linWidPub = 3
linWidPriv = 4
linWid = 4

t_max = 250
stepsize = .5

privFrac = 6/10
pubFrac = 1- privFrac
y0 = [99/100,1/100,0]
y0_mix = [y0[0]*privFrac,y0[1]*privFrac,y0[2]*privFrac,y0[0]*pubFrac,y0[1]*pubFrac,y0[2]*pubFrac]


figName = ''


#%%
sol_mix = [0,0,0]
i = 0
for params_dict in [pd0,pd1,pd2]:
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

    delt_min = params_dict['delt_min'] # min social value of saving a tree 
    delt_max = params_dict['delt_max'] # max social value of saving a tree

    Delta_m = V_m + W_m # value of saving a tree to society
    Delta_m_alt = params_dict['Delta_m_alt']
    params_dict["simType"] = 6
    fun_mix = lambda t,y: fn.hid_pub_priv(t,y,*args)
    sol_mix[i] = solve_ivp(fun_mix, [0, t_max], y0_mix ,method='LSODA', t_eval=t_eval)
    i = i+1

#%%
# Labels for figs, plotting options, etc.

titles = ['No assessment \n','Baseline assessment \naccuracy','Perfect assessmet \naccuracy']
labels = ['Healthy','Infested','Dying / dead']
labels1 = ['Healthy-Private','Infested-Private','Dying-Private','Healthy-Public','Infested-Public','Dying-Public']
panelAlphas = ['(a)','(b)','(c)','(d)','(e)','(f)']
linWid = [linWidPriv,linWidPriv,linWidPriv,linWidPub,linWidPub,linWidPub]
col = [healthy,infested,dead,healthy,infested,dead]
linSty = ['solid','solid','solid','dashed','dashed','dashed']
opacity = [opac,opac,opac,1,1,1]
textSize = 14


# Dynamics Figure 
fig, axs2 = plt.subplots(2, 3, figsize=(12, 7.66))

for i, ax in enumerate(axs2.flat):
    if i <3:
        for j in range(len(labels)):
            ax.plot(sol_mix[i].t, sol_mix[i].y[j]+sol_mix[i].y[j+3],label=labels[j], linewidth=linWid[j], color=col[j])
        ax.set_title(titles[i],fontsize=textSize)  
        ax.set_xlabel('Time (years)',fontsize=textSize)
        ax.set_ylabel('Fraction of trees',fontsize=textSize)
        ax.text(0.02, 1.125, panelAlphas[i], transform=ax.transAxes, va='top', fontsize=textSize, fontweight='bold')
        if i==0:
            ax.legend(bbox_to_anchor=(1, .65),loc='upper right',fontsize=12)
    if i> 2:
        for j in range(len(labels1)):
            ax.plot(sol_mix[i-3].t, sol_mix[i-3].y[j], label=labels1[j],linestyle=linSty[j], linewidth=linWid[j], color=col[j], alpha=opacity[j])
        ax.set_title(titles[i-3],fontsize=textSize)  
        ax.set_xlabel('Time (years)',fontsize=textSize)
        ax.set_ylabel('Fraction of trees',fontsize=textSize)
        ax.text(0.02, 1.125, panelAlphas[i], transform=ax.transAxes, va='top', fontsize=textSize, fontweight='bold')
        if i==3:
            ax.legend(bbox_to_anchor=(1, .68),loc='upper right',fontsize=12)


fig.tight_layout()

"SAVE FIGS"
if saveFigs == True:
    fig.savefig("FIGS/FigSM7"+figName+".png", bbox_inches='tight',dpi=250)


