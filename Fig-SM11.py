# -*- coding: utf-8 -*-
"""
Created on Fri May  5 15:44:15 2023

@author: AndrewTilman
"""
import numpy as np
import matplotlib.pyplot as plt
" MODEL PARAMETERS"
"Choose a ParamSet to import"
from ParamSet_middle_a0bHigh import params_dict
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

privFrac = 6/10
pubFrac = 1- privFrac
y0 = [99/100,1/100,0]
y0_mix = [y0[0]*privFrac,y0[1]*privFrac,y0[2]*privFrac,y0[0]*pubFrac,y0[1]*pubFrac,y0[2]*pubFrac]

beta = params_dict['beta']
gamma = params_dict['gamma']
alpha = params_dict['alpha']

t_eval = np.arange(0, t_max, stepsize)
args = (beta,gamma,alpha,params_dict)


"PLOT OPTIONS"

sMax2=params_dict['sMax2']
"OTHER PARAMS"

###! cost of treatment
c = params_dict['c']



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

figName = ''


#%%
"CALCS OF OPTIMAL SUBS for HEATMAPS"
# from matplotlib.cm import ScalarMappable

delta = 0.0021
levs = np.linspace(0, 250, 300)
levsPR = np.linspace(0, 1, 300)
colorMap = 'viridis'



dFrac = np.linspace(0.001, .99, 300) #frac dying / dead trees
iFrac = np.linspace(0.001, .99, 300) # frac infested trees
dFRAC, iFRAC = np.meshgrid(dFrac, iFrac)

sH = np.zeros((len(dFrac),len(iFrac))) #subsidy levels
sI = np.zeros((len(dFrac),len(iFrac)))
sD = np.zeros((len(dFrac),len(iFrac)))
ptH = np.zeros((len(dFrac),len(iFrac))) # Treatment probability for private tree
ptI = np.zeros((len(dFrac),len(iFrac)))
ptD = np.zeros((len(dFrac),len(iFrac)))
ptH_noSub = np.zeros((len(dFrac),len(iFrac))) # Treatment probability W/O subsidies for private tree
ptI_noSub = np.zeros((len(dFrac),len(iFrac)))
ptD_noSub = np.zeros((len(dFrac),len(iFrac)))
ptH_pub = np.zeros((len(dFrac),len(iFrac))) # treatment probability for public tree
ptI_pub = np.zeros((len(dFrac),len(iFrac)))
ptD_pub = np.zeros((len(dFrac),len(iFrac)))
dfSHOW = np.zeros((len(dFrac),len(iFrac))) # a rescaling to make equilateral triangles rather than right triangles for simplex


for i in range(len(dFrac)):
    for j in range(len(iFrac)):
        if dFrac[i] + iFrac[j] > 1:

            sH[i,j]=np.nan
            sI[i,j]=np.nan
            sD[i,j]=np.nan
            ptH[i,j]=np.nan
            ptI[i,j]=np.nan
            ptD[i,j]=np.nan
            ptH_noSub[i,j]=np.nan
            ptI_noSub[i,j]=np.nan
            ptD_noSub[i,j]=np.nan
            ptH_pub[i,j]=np.nan
            ptI_pub[i,j]=np.nan
            ptD_pub[i,j]=np.nan

        else:
            
            muh = fn.m_uh(t_horiz,iFRAC[i,j],beta,gamma)
            mth = fn.m_th(t_horiz,iFRAC[i,j],beta,gamma,eps_h,eps_i)
            mui = fn.m_ui(t_horiz,gamma)
            mti = fn.m_ti(t_horiz,gamma,eps_i)
            mud = 1
            mtd = 1
            
            t_h = muh - mth 
            t_i = mui - mti
            t_d = mud - mtd
            
            hFrac = 1-(dFrac[i] + iFrac[j])
            smuh = fn.sm_uh(t_horiz,iFRAC[i,j],hFrac,beta,gamma)
            smth = fn.sm_th(t_horiz,iFRAC[i,j],hFrac,beta,gamma,alpha,eps_h,eps_i)
            smui = fn.sm_ui(hFrac,beta,gamma)
            smti = fn.sm_ti(hFrac,beta,gamma,alpha,eps_i)
            smud = 0
            smtd = 0

            s_h = (smuh - smth) * spilloverCare
            s_i = (smui - smti) * spilloverCare
            s_d = (smud - smtd) * spilloverCare
            
            
            kHtemp = fn.k_X(t_h,t_i,t_d,\
                                fn.prob_xH('h',1-dFRAC[i,j]-iFRAC[i,j],iFRAC[i,j],dFRAC[i,j],pr_Hh,pr_Hi,pr_Hd),\
                                fn.prob_xH('i',1-dFRAC[i,j]-iFRAC[i,j],iFRAC[i,j],dFRAC[i,j],pr_Hh,pr_Hi,pr_Hd),\
                                fn.prob_xH('d',1-dFRAC[i,j]-iFRAC[i,j],iFRAC[i,j],dFRAC[i,j],pr_Hh,pr_Hi,pr_Hd))
            kItemp = fn.k_X(t_h,t_i,t_d,\
                                fn.prob_xI('h',1-dFRAC[i,j]-iFRAC[i,j],iFRAC[i,j],dFRAC[i,j],pr_Ih,pr_Ii,pr_Id),\
                                fn.prob_xI('i',1-dFRAC[i,j]-iFRAC[i,j],iFRAC[i,j],dFRAC[i,j],pr_Ih,pr_Ii,pr_Id),\
                                fn.prob_xI('d',1-dFRAC[i,j]-iFRAC[i,j],iFRAC[i,j],dFRAC[i,j],pr_Ih,pr_Ii,pr_Id))
            kDtemp = fn.k_X(t_h,t_i,t_d,\
                                fn.prob_xD('h',1-dFRAC[i,j]-iFRAC[i,j],iFRAC[i,j],dFRAC[i,j],pr_Dh,pr_Di,pr_Dd),\
                                fn.prob_xD('i',1-dFRAC[i,j]-iFRAC[i,j],iFRAC[i,j],dFRAC[i,j],pr_Dh,pr_Di,pr_Dd),\
                                fn.prob_xD('d',1-dFRAC[i,j]-iFRAC[i,j],iFRAC[i,j],dFRAC[i,j],pr_Dh,pr_Di,pr_Dd))
                
            LHtemp = fn.L_X(s_h,s_i,s_d,\
                                fn.prob_xH('h',1-dFRAC[i,j]-iFRAC[i,j],iFRAC[i,j],dFRAC[i,j],pr_Hh,pr_Hi,pr_Hd),\
                                fn.prob_xH('i',1-dFRAC[i,j]-iFRAC[i,j],iFRAC[i,j],dFRAC[i,j],pr_Hh,pr_Hi,pr_Hd),\
                                fn.prob_xH('d',1-dFRAC[i,j]-iFRAC[i,j],iFRAC[i,j],dFRAC[i,j],pr_Hh,pr_Hi,pr_Hd))
            LItemp = fn.L_X(s_h,s_i,s_d,\
                                fn.prob_xI('h',1-dFRAC[i,j]-iFRAC[i,j],iFRAC[i,j],dFRAC[i,j],pr_Ih,pr_Ii,pr_Id),\
                                fn.prob_xI('i',1-dFRAC[i,j]-iFRAC[i,j],iFRAC[i,j],dFRAC[i,j],pr_Ih,pr_Ii,pr_Id),\
                                fn.prob_xI('d',1-dFRAC[i,j]-iFRAC[i,j],iFRAC[i,j],dFRAC[i,j],pr_Ih,pr_Ii,pr_Id))
            LDtemp = fn.L_X(s_h,s_i,s_d,\
                                fn.prob_xD('h',1-dFRAC[i,j]-iFRAC[i,j],iFRAC[i,j],dFRAC[i,j],pr_Dh,pr_Di,pr_Dd),\
                                fn.prob_xD('i',1-dFRAC[i,j]-iFRAC[i,j],iFRAC[i,j],dFRAC[i,j],pr_Dh,pr_Di,pr_Dd),\
                                fn.prob_xD('d',1-dFRAC[i,j]-iFRAC[i,j],iFRAC[i,j],dFRAC[i,j],pr_Dh,pr_Di,pr_Dd))
                                

            sH[i,j] = fn.s_opt(c,kHtemp,LHtemp,Delta_m,a,b)
            sI[i,j] = fn.s_opt(c,kItemp,LItemp,Delta_m,a,b)
            sD[i,j] = fn.s_opt(c,kDtemp,LDtemp,Delta_m,a,b)

            ptH[i,j] = fn.pr_treated(c,kHtemp,LHtemp,Delta_m,a,b)
            ptI[i,j] = fn.pr_treated(c,kItemp,LItemp,Delta_m,a,b)
            ptD[i,j] = fn.pr_treated(c,kDtemp,LDtemp,Delta_m,a,b)
            
            ptH_noSub[i,j] = fn.pr_treated(c,kHtemp,LHtemp,0,a,b)
            ptI_noSub[i,j] = fn.pr_treated(c,kItemp,LItemp,0,a,b)
            ptD_noSub[i,j] = fn.pr_treated(c,kDtemp,LDtemp,0,a,b)
            
            ptH_pub[i,j] = fn.pr_treated_pub(c,kHtemp,LHtemp,Delta_m_alt,a,b)
            ptI_pub[i,j] = fn.pr_treated_pub(c,kItemp,LItemp,Delta_m_alt,a,b)
            ptD_pub[i,j] = fn.pr_treated_pub(c,kDtemp,LDtemp,Delta_m_alt,a,b)
            

            dfSHOW[i,j]= dFRAC[i,j] + 1/2*iFRAC[i,j]
#%%
"PLOTS OF OPTIMAL SUBS WITH HEATMAPS"

"LABELS AND SUCH"
rowTitles = ['Assessed \nhealthy','','','','Assessed \ninfested','','','','Assessed \ndying / dead','','','']
colTitles = [' Optimal subsidy \n','Treatment probability \n   w/ optimal subsidy', 'Treatment probability \n   w/o subsidy','Treatment probability \n   for public trees','','','','','','','','','','','','','','','','','','']
heatMapsData = [sH,ptH,ptH_noSub,ptH_pub,sI,ptI,ptI_noSub,ptI_pub,sD,ptD,ptD_noSub,ptD_pub]
cbMax = [sMax2+1,2,2,2,sMax2+1,2,2,2,sMax2+1,2,2,2]
cbSteps = [50,1,1,1,50,1,1,1,50,1,1,1]
levsList = [levs,levsPR,levsPR,levsPR,levs,levsPR,levsPR,levsPR,levs,levsPR,levsPR,levsPR]
panelAlphas = ['(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)','(j)','(k)','(l)']

"PLOTTING"
fig4, bxs = plt.subplots(3, 4, figsize=(15.3, 9))

for i, bx in enumerate(bxs.flat):
    bx.text(.35,1.05,'All infested',size='12')
    bx.text(-.05,-.1,'All healthy',size='12')
    bx.text(.76,-.1,'All dying',size='12')
    bx.text(-.5,.5,rowTitles[i],size='16')
    bx.text(.12,1.25,colTitles[i],size='16')
    bx.set_frame_on(None)
    BX = bx.contourf( dfSHOW,iFRAC,heatMapsData[i],levels=levsList[i],cmap=colorMap)
    bx.set_ylim([0,1])
    bx.set_xlim([0,1])
    fig4.colorbar(BX,ticks=range(0,cbMax[i], cbSteps[i]),ax=bx)
    bx.set_xticks(())
    bx.set_yticks(())
    bx.text(-0.05, 1.08, panelAlphas[i], transform=bx.transAxes, va='top', fontsize=12, fontweight='bold')

fig4.tight_layout()

"SAVE PLOTS?"
if saveFigs == True:
    fig4.savefig("FIGS/FigSM11"+figName+".png", bbox_inches='tight',dpi=250)

