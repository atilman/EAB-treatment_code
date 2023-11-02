# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
" MODEL PARAMETERS"
"Choose a ParamSet 1 2 3 or 4"
from ParamSet1 import params_dict
import functions as fn
from scipy.integrate import solve_ivp

"SAVE FIGS???"
saveFigs = False


"COLORS"
healthy = '#6DB33F'
infested = '#CFCA7A'
dead = '#8e7028'
opac = .75

"PLOT OPTIONS"
linWidPub = 6
linWidPriv = 3
linWid = 4

t_max=50
stepsize=.05
privFrac = 6/10
pubFrac = 1- privFrac
y0 = [99/100,1/100,0]
y0_mix = [y0[0]*privFrac,y0[1]*privFrac,y0[2]*privFrac,y0[0]*pubFrac,y0[1]*pubFrac,y0[2]*pubFrac]

"Epidemiological parameters"
beta = 1
gamma = 3/10
alpha = 1
t_eval = np.arange(0, t_max, stepsize)
args = (beta,gamma,alpha,params_dict)



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

Delta_m = V_m + W_m # value of saving a tree to society


figName = params_dict['paramSet']

# if vFrac + wFrac != 1:
#     print('fix your value fractions! vFrac / wFrac so they add to one!')
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


#%%
"SIMTYPE"
#! SimType defines the Subsidy and treatment probs for simulation
    #! 1: no private treatment or public treatment
    #! 2: no private treatment but optimal public treatment
    #! 3: no private subsidies and no public treatment
    #! 4: no private subsidies but optimal public treatment
    #! 5: optimal private subsidies but no public treatment
    #! 6: optimal private subsidies and optimal public treatment
    
sol_mix = [0,0,0,0,0,0]
for i in range(6):
    params_dict["simType"] = i+1
    fun_mix = lambda t,y: fn.hid_pub_priv(t,y,*args)
    sol_mix[i] = solve_ivp(fun_mix, [0, t_max], y0_mix ,method='LSODA', t_eval=t_eval)


#%%
fig3 = plt.figure(figsize=(12,6.66))
dx1 = fig3.add_subplot(231)
dx2 = fig3.add_subplot(234)
dx3 = fig3.add_subplot(232)
dx4 = fig3.add_subplot(235)
dx5 = fig3.add_subplot(233)
dx6 = fig3.add_subplot(236)

fig2 = plt.figure(figsize=(12,6.66))
cx1 = fig2.add_subplot(231)
cx2 = fig2.add_subplot(234)
cx3 = fig2.add_subplot(232)
cx4 = fig2.add_subplot(235)
cx5 = fig2.add_subplot(233)
cx6 = fig2.add_subplot(236)

dx1.plot(sol_mix[0].t, sol_mix[0].y[0],label='h(t)-Private',linewidth=linWidPriv,color=healthy)
dx1.plot(sol_mix[0].t, sol_mix[0].y[1],label='i(t)-Private',linewidth=linWidPriv,color=infested)
dx1.plot(sol_mix[0].t, sol_mix[0].y[2],label='d(t)-Private',linewidth=linWidPriv,color=dead)
dx1.plot(sol_mix[0].t, sol_mix[0].y[3],label='h(t)-Public',linewidth=linWidPub,color=healthy,alpha=opac)
dx1.plot(sol_mix[0].t, sol_mix[0].y[4],label='i(t)-Public',linewidth=linWidPub,color=infested,alpha=opac)
dx1.plot(sol_mix[0].t, sol_mix[0].y[5],label='d(t)-Public',linewidth=linWidPub,color=dead,alpha=opac)

dx2.plot(sol_mix[1].t, sol_mix[1].y[0],label='h(t)-Private',linewidth=linWidPriv,color=healthy)
dx2.plot(sol_mix[1].t, sol_mix[1].y[1],label='i(t)-Private',linewidth=linWidPriv,color=infested)
dx2.plot(sol_mix[1].t, sol_mix[1].y[2],label='d(t)-Private',linewidth=linWidPriv,color=dead)
dx2.plot(sol_mix[1].t, sol_mix[1].y[3],label='h(t)-Public',linewidth=linWidPub,color=healthy,alpha=opac)
dx2.plot(sol_mix[1].t, sol_mix[1].y[4],label='i(t)-Public',linewidth=linWidPub,color=infested,alpha=opac)
dx2.plot(sol_mix[1].t, sol_mix[1].y[5],label='d(t)-Public',linewidth=linWidPub,color=dead,alpha=opac)

dx3.plot(sol_mix[2].t, sol_mix[2].y[0],label='h(t)-Private',linewidth=linWidPriv,color=healthy)
dx3.plot(sol_mix[2].t, sol_mix[2].y[1],label='i(t)-Private',linewidth=linWidPriv,color=infested)
dx3.plot(sol_mix[2].t, sol_mix[2].y[2],label='d(t)-Private',linewidth=linWidPriv,color=dead)
dx3.plot(sol_mix[2].t, sol_mix[2].y[3],label='h(t)-Public',linewidth=linWidPub,color=healthy,alpha=opac)
dx3.plot(sol_mix[2].t, sol_mix[2].y[4],label='i(t)-Public',linewidth=linWidPub,color=infested,alpha=opac)
dx3.plot(sol_mix[2].t, sol_mix[2].y[5],label='d(t)-Public',linewidth=linWidPub,color=dead,alpha=opac)

dx4.plot(sol_mix[3].t, sol_mix[3].y[0],label='h(t)-Private',linewidth=linWidPriv,color=healthy)
dx4.plot(sol_mix[3].t, sol_mix[3].y[1],label='i(t)-Private',linewidth=linWidPriv,color=infested)
dx4.plot(sol_mix[3].t, sol_mix[3].y[2],label='d(t)-Private',linewidth=linWidPriv,color=dead)
dx4.plot(sol_mix[3].t, sol_mix[3].y[3],label='h(t)-Public',linewidth=linWidPub,color=healthy,alpha=opac)
dx4.plot(sol_mix[3].t, sol_mix[3].y[4],label='i(t)-Public',linewidth=linWidPub,color=infested,alpha=opac)
dx4.plot(sol_mix[3].t, sol_mix[3].y[5],label='d(t)-Public',linewidth=linWidPub,color=dead,alpha=opac)

dx5.plot(sol_mix[4].t, sol_mix[4].y[0],label='h(t)-Private',linewidth=linWidPriv,color=healthy)
dx5.plot(sol_mix[4].t, sol_mix[4].y[1],label='i(t)-Private',linewidth=linWidPriv,color=infested)
dx5.plot(sol_mix[4].t, sol_mix[4].y[2],label='d(t)-Private',linewidth=linWidPriv,color=dead)
dx5.plot(sol_mix[4].t, sol_mix[4].y[3],label='h(t)-Public',linewidth=linWidPub,color=healthy,alpha=opac)
dx5.plot(sol_mix[4].t, sol_mix[4].y[4],label='i(t)-Public',linewidth=linWidPub,color=infested,alpha=opac)
dx5.plot(sol_mix[4].t, sol_mix[4].y[5],label='d(t)-Public',linewidth=linWidPub,color=dead,alpha=opac)

dx6.plot(sol_mix[5].t, sol_mix[5].y[0],label='H -Private',linewidth=linWidPriv,color=healthy)
dx6.plot(sol_mix[5].t, sol_mix[5].y[1],label='I -Private',linewidth=linWidPriv,color=infested)
dx6.plot(sol_mix[5].t, sol_mix[5].y[2],label='D -Private',linewidth=linWidPriv,color=dead)
dx6.plot(sol_mix[5].t, sol_mix[5].y[3],label='H -Public',linewidth=linWidPub,color=healthy,alpha=opac)
dx6.plot(sol_mix[5].t, sol_mix[5].y[4],label='I -Public',linewidth=linWidPub,color=infested,alpha=opac)
dx6.plot(sol_mix[5].t, sol_mix[5].y[5],label='D -Public',linewidth=linWidPub,color=dead,alpha=opac)

dx2.set_xlabel('time (years)',fontsize=14)
dx4.set_xlabel('time (years)',fontsize=14)
dx6.set_xlabel('time (years)',fontsize=14)
dx1.set_title('No private treatment',fontsize=14)
dx3.set_title('No private subsidies',fontsize=14)
dx5.set_title('Optimal private subsidies',fontsize=14)
dx1.set_ylabel('No public treatment',fontsize=14)
dx2.set_ylabel('Optimal public treatment',fontsize=14)
dx6.legend(loc='right',fontsize=12)

cx1.plot(sol_mix[0].t, sol_mix[0].y[0]+sol_mix[0].y[3],label='h(t)',linewidth=linWid,color=healthy)
cx1.plot(sol_mix[0].t, sol_mix[0].y[1]+sol_mix[0].y[4],label='i(t)',linewidth=linWid,color=infested)
cx1.plot(sol_mix[0].t, sol_mix[0].y[2]+sol_mix[0].y[5],label='d(t)',linewidth=linWid,color=dead)

cx2.plot(sol_mix[1].t, sol_mix[1].y[0]+sol_mix[1].y[3],label='h(t)',linewidth=linWid,color=healthy)
cx2.plot(sol_mix[1].t, sol_mix[1].y[1]+sol_mix[1].y[4],label='i(t)',linewidth=linWid,color=infested)
cx2.plot(sol_mix[1].t, sol_mix[1].y[2]+sol_mix[1].y[5],label='d(t)',linewidth=linWid,color=dead)

cx3.plot(sol_mix[2].t, sol_mix[2].y[0]+sol_mix[2].y[3],label='h(t)',linewidth=linWid,color=healthy)
cx3.plot(sol_mix[2].t, sol_mix[2].y[1]+sol_mix[2].y[4],label='i(t)',linewidth=linWid,color=infested)
cx3.plot(sol_mix[2].t, sol_mix[2].y[2]+sol_mix[2].y[5],label='d(t)',linewidth=linWid,color=dead)

cx4.plot(sol_mix[3].t, sol_mix[3].y[0]+sol_mix[3].y[3],label='h(t)',linewidth=linWid,color=healthy)
cx4.plot(sol_mix[3].t, sol_mix[3].y[1]+sol_mix[3].y[4],label='i(t)',linewidth=linWid,color=infested)
cx4.plot(sol_mix[3].t, sol_mix[3].y[2]+sol_mix[3].y[5],label='d(t)',linewidth=linWid,color=dead)

cx5.plot(sol_mix[4].t, sol_mix[4].y[0]+sol_mix[4].y[3],label='h(t)',linewidth=linWid,color=healthy)
cx5.plot(sol_mix[4].t, sol_mix[4].y[1]+sol_mix[4].y[4],label='i(t)',linewidth=linWid,color=infested)
cx5.plot(sol_mix[4].t, sol_mix[4].y[2]+sol_mix[4].y[5],label='d(t)',linewidth=linWid,color=dead)

cx6.plot(sol_mix[5].t, sol_mix[5].y[0]+sol_mix[5].y[3],label='Healthy',linewidth=linWid,color=healthy)
cx6.plot(sol_mix[5].t, sol_mix[5].y[1]+sol_mix[5].y[4],label='Infested',linewidth=linWid,color=infested)
cx6.plot(sol_mix[5].t, sol_mix[5].y[2]+sol_mix[5].y[5],label='Dying / dead',linewidth=linWid,color=dead)



cx2.set_xlabel('time (years)',fontsize=14)
cx4.set_xlabel('time (years)',fontsize=14)
cx6.set_xlabel('time (years)',fontsize=14)
cx1.set_title('No private treatment',fontsize=14)
cx3.set_title('No private subsidies',fontsize=14)
cx5.set_title('Optimal private subsidies',fontsize=14)
cx1.set_ylabel('No public treatment',fontsize=14)
cx2.set_ylabel('Optimal public treatment',fontsize=14)
cx6.legend(loc='right',fontsize=12)

"SAVE FIGS"
if saveFigs == True:
    fig2.savefig("FIGS/dynamMixedAvgd_"+figName+".pdf",bbox_inches='tight', dpi=150)
    fig2.savefig("FIGS/dynamMixedAvgd_"+figName+".png", bbox_inches='tight',dpi=250)
    fig3.savefig("FIGS/dynamMixedSep_"+figName+".pdf",bbox_inches='tight', dpi=150)
    fig3.savefig("FIGS/dynamMixedSep_"+figName+".png", bbox_inches='tight',dpi=250)