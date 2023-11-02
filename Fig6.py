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
opac = .75

"PLOT OPTIONS"
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
Delta_m_alt = params_dict['Delta_m_alt']

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
"PLOTS OF OPTIMAL SUBS WITH HEATMAPS"
# from matplotlib.cm import ScalarMappable

delta = 0.0021
levs = np.linspace(0, 125, 300)
levsPR = np.linspace(0, 1, 300)
colorMap = 'viridis'

dFrac = np.linspace(0.001, .99, 300)
iFrac = np.linspace(0.001, .99, 300)
dFRAC, iFRAC = np.meshgrid(dFrac, iFrac)

sH = np.zeros((len(dFrac),len(iFrac)))
sI = np.zeros((len(dFrac),len(iFrac)))
sD = np.zeros((len(dFrac),len(iFrac)))
ptH = np.zeros((len(dFrac),len(iFrac)))
ptI = np.zeros((len(dFrac),len(iFrac)))
ptD = np.zeros((len(dFrac),len(iFrac)))
ptH_noSub = np.zeros((len(dFrac),len(iFrac)))
ptI_noSub = np.zeros((len(dFrac),len(iFrac)))
ptD_noSub = np.zeros((len(dFrac),len(iFrac)))
ptH_pub = np.zeros((len(dFrac),len(iFrac)))
ptI_pub = np.zeros((len(dFrac),len(iFrac)))
ptD_pub = np.zeros((len(dFrac),len(iFrac)))
dfSHOW = np.zeros((len(dFrac),len(iFrac)))

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
                                
            sH[i,j] = fn.s_opt(c,kHtemp,Delta_m,a,b)
            sI[i,j] = fn.s_opt(c,kItemp,Delta_m,a,b)
            sD[i,j] = fn.s_opt(c,kDtemp,Delta_m,a,b)
            
            ptH[i,j] = fn.pr_treated(c,kHtemp,Delta_m,a,b)
            ptI[i,j] = fn.pr_treated(c,kItemp,Delta_m,a,b)
            ptD[i,j] = fn.pr_treated(c,kDtemp,Delta_m,a,b)
            
            ptH_noSub[i,j] = fn.pr_treated(c,kHtemp,0,a,b)
            ptI_noSub[i,j] = fn.pr_treated(c,kItemp,0,a,b)
            ptD_noSub[i,j] = fn.pr_treated(c,kDtemp,0,a,b)
            
            ptH_pub[i,j] = fn.pr_treated_pub(c,kHtemp,Delta_m_alt,a,b)
            ptI_pub[i,j] = fn.pr_treated_pub(c,kItemp,Delta_m_alt,a,b)
            ptD_pub[i,j] = fn.pr_treated_pub(c,kDtemp,Delta_m_alt,a,b)
            
            dfSHOW[i,j]= dFRAC[i,j] + 1/2*iFRAC[i,j]

fig2 = plt.figure(figsize=(15.3,9))

bx1 = fig2.add_subplot(341)
plt.text(.35,1.05,'$P(i)=1$',size='12')
plt.text(-.05,-.1,'$P(h)=1$',size='12')
plt.text(.75,-.1,'$P(d)=1$',size='12')
plt.text(-.25,.5,'$\hat h$',size='20')
plt.text(.12,1.25,'Optimal Subsidy',size='16')
plt.box(on=None)

cx1 = fig2.add_subplot(342)
plt.text(.35,1.05,'$P(i)=1$',size='12')
plt.text(-.05,-.1,'$P(h)=1$',size='12')
plt.text(.75,-.1,'$P(d)=1$',size='12')
plt.text(.10,1.25,'Treatment probability',size='16')
plt.text(.2,1.15,'w/ optimal subsidy',size='16')
plt.box(on=None)

dx1 = fig2.add_subplot(343)
plt.text(.35,1.05,'$P(i)=1$',size='12')
plt.text(-.05,-.1,'$P(h)=1$',size='12')
plt.text(.75,-.1,'$P(d)=1$',size='12')
plt.text(.10,1.25,'Treatment probability',size='16')
plt.text(.2,1.15,'w/o subsidy',size='16')
plt.box(on=None)

ex1 = fig2.add_subplot(344)
plt.text(.35,1.05,'$P(i)=1$',size='12')
plt.text(-.05,-.1,'$P(h)=1$',size='12')
plt.text(.75,-.1,'$P(d)=1$',size='12')
plt.text(.10,1.25,'Treatment probability',size='16')
plt.text(.2,1.15,'for public trees',size='16')
plt.box(on=None)

bx2 = fig2.add_subplot(345)
plt.text(.35,1.05,'$P(i)=1$',size='12')
plt.text(-.05,-.1,'$P(h)=1$',size='12')
plt.text(.75,-.1,'$P(d)=1$',size='12')
plt.text(-.25,.5,'$\hat i$',size='20')
plt.box(on=None)

cx2 = fig2.add_subplot(346)
plt.text(.35,1.05,'$P(i)=1$',size='12')
plt.text(-.05,-.1,'$P(h)=1$',size='12')
plt.text(.75,-.1,'$P(d)=1$',size='12')
plt.box(on=None)

dx2 = fig2.add_subplot(347)
plt.text(.35,1.05,'$P(i)=1$',size='12')
plt.text(-.05,-.1,'$P(h)=1$',size='12')
plt.text(.75,-.1,'$P(d)=1$',size='12')
plt.box(on=None)

ex2 = fig2.add_subplot(348)
plt.text(.35,1.05,'$P(i)=1$',size='12')
plt.text(-.05,-.1,'$P(h)=1$',size='12')
plt.text(.75,-.1,'$P(d)=1$',size='12')
plt.box(on=None)

bx3 = fig2.add_subplot(3,4,9)
plt.text(.35,1.05,'$P(i)=1$',size='12')
plt.text(-.05,-.1,'$P(h)=1$',size='12')
plt.text(.75,-.1,'$P(d)=1$',size='12')
plt.text(-.25,.5,'$\hat d$',size='20')
plt.box(on=None)

cx3 = fig2.add_subplot(3,4,10)
plt.text(.35,1.05,'$P(i)=1$',size='12')
plt.text(-.05,-.1,'$P(h)=1$',size='12')
plt.text(.75,-.1,'$P(d)=1$',size='12')
plt.box(on=None)

dx3 = fig2.add_subplot(3,4,11)
plt.text(.35,1.05,'$P(i)=1$',size='12')
plt.text(-.05,-.1,'$P(h)=1$',size='12')
plt.text(.75,-.1,'$P(d)=1$',size='12')
plt.box(on=None)

ex3 = fig2.add_subplot(3,4,12)
plt.text(.35,1.05,'$P(i)=1$',size='12')
plt.text(-.05,-.1,'$P(h)=1$',size='12')
plt.text(.75,-.1,'$P(d)=1$',size='12')
plt.box(on=None)

BX1 = bx1.contourf( dfSHOW,iFRAC,sH,levels=levs,cmap=colorMap)
bx1.set_ylim([0,1])
bx1.set_xlim([0,1])
fig2.colorbar(BX1,ticks=range(0,sMax2+1, 50),ax=bx1)
bx1.set_xticks([])
bx1.set_yticks([])

CX1 = cx1.contourf( dfSHOW,iFRAC,ptH,levels=levsPR,cmap=colorMap)
cx1.set_ylim([0,1])
cx1.set_xlim([0,1])
fig2.colorbar(CX1,ticks=range(0,2,1),ax=cx1)
cx1.set_xticks(())
cx1.set_yticks(())

DX1 = dx1.contourf( dfSHOW,iFRAC,ptH_noSub,levels=levsPR,cmap=colorMap)
dx1.set_ylim([0,1])
dx1.set_xlim([0,1])
fig2.colorbar(DX1,ticks=range(0,2,1),ax=dx1)
dx1.set_xticks(())
dx1.set_yticks(())

EX1 = ex1.contourf( dfSHOW,iFRAC,ptH_pub,levels=levsPR,cmap=colorMap)
ex1.set_ylim([0,1])
ex1.set_xlim([0,1])
fig2.colorbar(EX1,ticks=range(0,2,1),ax=ex1)
ex1.set_xticks(())
ex1.set_yticks(())

BX2 = bx2.contourf( dfSHOW,iFRAC,sI,levels=levs,cmap=colorMap)
bx2.set_ylim([0,1])
bx2.set_xlim([0,1])
fig2.colorbar(BX2,ticks=range(0,sMax2+1, 50),ax=bx2)
bx2.set_xticks(())
bx2.set_yticks(())

CX2 = cx2.contourf( dfSHOW,iFRAC,ptI,levels=levsPR,cmap=colorMap)
cx2.set_ylim([0,1])
cx2.set_xlim([0,1])
fig2.colorbar(CX2,ticks=range(0,2,1),ax=cx2)
cx2.set_xticks(())
cx2.set_yticks(())

DX2 = dx2.contourf( dfSHOW,iFRAC,ptI_noSub,levels=levsPR,cmap=colorMap)
dx2.set_ylim([0,1])
dx2.set_xlim([0,1])
fig2.colorbar(DX2,ticks=range(0,2,1),ax=dx2)
dx2.set_xticks(())
dx2.set_yticks(())

EX2 = ex2.contourf( dfSHOW,iFRAC,ptI_pub,levels=levsPR,cmap=colorMap)
ex2.set_ylim([0,1])
ex2.set_xlim([0,1])
fig2.colorbar(EX2,ticks=range(0,2,1),ax=ex2)
ex2.set_xticks(())
ex2.set_yticks(())

BX3 = bx3.contourf( dfSHOW,iFRAC,sD,levels=levs,cmap=colorMap)
bx3.set_ylim([0,1])
bx3.set_xlim([0,1])
fig2.colorbar(BX3,ticks=range(0,sMax2+1, 50),ax=bx3);
bx3.set_xticks(())
bx3.set_yticks(())


CX3 = cx3.contourf( dfSHOW,iFRAC,ptD,levels=levsPR,cmap=colorMap)
cx3.set_ylim([0,1])
cx3.set_xlim([0,1])
fig2.colorbar(CX3,ticks=range(0,2,1),ax=cx3)
cx3.set_xticks(())
cx3.set_yticks(())

DX3 = dx3.contourf( dfSHOW,iFRAC,ptD_noSub,levels=levsPR,cmap=colorMap)
dx3.set_ylim([0,1])
dx3.set_xlim([0,1])
fig2.colorbar(DX3,ticks=range(0,2,1),ax=dx3)
dx3.set_xticks(())
dx3.set_yticks(())

EX3 = ex3.contourf( dfSHOW,iFRAC,ptD_pub,levels=levsPR,cmap=colorMap)
ex3.set_ylim([0,1])
ex3.set_xlim([0,1])
fig2.colorbar(EX3,ticks=range(0,2,1),ax=ex3)
ex3.set_xticks(())
ex3.set_yticks(())

plt.tight_layout()

if saveFigs == True:
    fig2.savefig("FIGS/simplex_"+figName+".pdf",bbox_inches='tight', dpi=150)
    fig2.savefig("FIGS/simplex_"+figName+".png", bbox_inches='tight',dpi=250)