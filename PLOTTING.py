# -*- coding: utf-8 -*-
"""
Created on Tue Dec 18 11:15:37 2018

@author: toppgabr
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Dec 17 16:00:25 2018

@author: toppgabr
"""

import matplotlib.pyplot as plt  
import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib as mpl
from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition, mark_inset)
from matplotlib.colors import ListedColormap
#from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition, mark_inset)
    
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['lines.linewidth'] = 2

mpl.rcParams['lines.markersize'] = 10
mpl.rcParams['font.size'] = 20  # <-- change fonsize globally
mpl.rcParams['legend.fontsize'] = 14
mpl.rcParams['axes.titlesize'] = 20
mpl.rcParams['axes.labelsize'] = 20
mpl.rcParams['axes.linewidth'] = 1.5
mpl.rcParams['xtick.major.size'] = 10
mpl.rcParams['ytick.major.size'] = 10
mpl.rcParams['xtick.major.width'] = 1.5
mpl.rcParams['ytick.major.width'] = 1.5
mpl.rcParams['xtick.direction'] = 'out'
mpl.rcParams['ytick.direction'] = 'out'
mpl.rcParams['figure.titlesize'] = 20
mpl.rcParams['figure.figsize'] = [12,5]
mpl.rcParams['text.usetex'] = True

C1 = "crimson"
C2 = "orange"
C3 = "#117A65"
C4 = "#48C9B0"
C5 = "#34495E"

m_max = 1
nband=(2*m_max+1)
omega = 1.5
kk = np.linspace(0,1,2)
nokp = 1
RANGE=1 
 
BW=np.zeros(RANGE)
BW_EFF=np.zeros(RANGE)
BW_EFF0=np.zeros(RANGE)
BW_EFF1=np.zeros(RANGE)
BW_EFF2=np.zeros(RANGE)
BW_EFF3=np.zeros(RANGE)
BW_EFF4=np.zeros(RANGE)
BW_EFF5=np.zeros(RANGE)
BW_EFF6=np.zeros(RANGE)
BW_EFF7=np.zeros(RANGE)
BW_EFF8=np.zeros(RANGE)
BW_EFF9=np.zeros(RANGE)

FGAPS=np.zeros(RANGE)
FGAPS_EFF=np.zeros(RANGE)
FGAPS_EFF0=np.zeros(RANGE)
FGAPS_EFF1=np.zeros(RANGE)
FGAPS_EFF2=np.zeros(RANGE)
FGAPS_EFF3=np.zeros(RANGE)
FGAPS_EFF4=np.zeros(RANGE)
FGAPS_EFF5=np.zeros(RANGE)
FGAPS_EFF6=np.zeros(RANGE)
FGAPS_EFF7=np.zeros(RANGE)
FGAPS_EFF8=np.zeros(RANGE)
FGAPS_EFF9=np.zeros(RANGE)

AA = np.array([0.0, 0.005, 0.01, 0.015, 0.02, 0.025])

Nmax = 512
dim_new= Nmax

case = 0
for mm in range(RANGE):
    
    print(mm)
    case = mm
    
    file_BANDS = open('DATA/bands_DOWN_'+str(case)+'.dat','r')
    MAT_BANDS_DOWN = np.loadtxt(file_BANDS)
    file_BANDS.close()
    
    mu = (MAT_BANDS_DOWN[0,int(Nmax/2)]+MAT_BANDS_DOWN[0,int(Nmax/2)-1])*0.5
    MAT_BANDS_DOWN = MAT_BANDS_DOWN - mu
    
    file_BANDS = open('DATA/bands_floquet_DOWN_'+str(case)+'.dat','r')
    FLOQUET_BANDS_DOWN = np.loadtxt(file_BANDS)-mu
    file_BANDS.close()
    
    file = open('DATA/FLOQUET_EFF_'+str(case)+'.dat','r')
    MAT_BANDS_FLOQUET_EFF = np.loadtxt(file)-mu
    file.close()

    file = open('DATA/FLOQUET_EFF0_'+str(case)+'.dat','r')
    MAT_BANDS_FLOQUET_EFF0 = np.loadtxt(file)
    file.close()
            
    file = open('DATA/FLOQUET_EFF1_'+str(case)+'.dat','r')
    MAT_BANDS_FLOQUET_EFF1 = np.loadtxt(file)
    file.close()
        
    file = open('DATA/FLOQUET_EFF2_'+str(case)+'.dat','r')
    MAT_BANDS_FLOQUET_EFF2 = np.loadtxt(file)
    file.close()
       
    file = open('DATA/FLOQUET_EFF3_'+str(case)+'.dat','r')
    MAT_BANDS_FLOQUET_EFF3 = np.loadtxt(file)
    file.close()
    
    file = open('DATA/FLOQUET_EFF4_'+str(case)+'.dat','r')
    MAT_BANDS_FLOQUET_EFF4 = np.loadtxt(file)
    file.close()
    
    file = open('DATA/FLOQUET_EFF5_'+str(case)+'.dat','r')
    MAT_BANDS_FLOQUET_EFF5 = np.loadtxt(file)
    file.close()
    
    file = open('DATA/FLOQUET_EFF6_'+str(case)+'.dat','r')
    MAT_BANDS_FLOQUET_EFF6 = np.loadtxt(file)
    file.close()
    
    file = open('DATA/FLOQUET_EFF7_'+str(case)+'.dat','r')
    MAT_BANDS_FLOQUET_EFF7 = np.loadtxt(file)
    file.close()

    file = open('DATA/FLOQUET_EFF8_'+str(case)+'.dat','r')
    MAT_BANDS_FLOQUET_EFF8 = np.loadtxt(file)
    file.close()    
    
    file = open('DATA/FLOQUET_EFF9_'+str(case)+'.dat','r')
    MAT_BANDS_FLOQUET_EFF9 = np.loadtxt(file)
    file.close()    

    ## FULLHamiltonian
    file_BANDS = open('DATA/bands_DOWN_'+str(case)+'.dat','r')
    MAT_BANDS = np.loadtxt(file_BANDS)-mu
    file_BANDS.close()
    
    ## Measure bandwidth
    BW_GS = (MAT_BANDS_DOWN[int(nokp),int(Nmax/2)]-MAT_BANDS_DOWN[int(nokp),int(Nmax/2)-1])*1000
    BW[mm] = (FLOQUET_BANDS_DOWN[int(nokp),int(m_max*Nmax+Nmax/2)]-FLOQUET_BANDS_DOWN[int(nokp),int(m_max*Nmax+Nmax/2)-1])*1000
    BW_EFF[mm] = (MAT_BANDS_FLOQUET_EFF[int(1.0*nokp),int(Nmax/2)]-MAT_BANDS_FLOQUET_EFF[int(1.0*nokp),int(Nmax/2)-1])*1000
    BW_EFF0[mm] = (MAT_BANDS_FLOQUET_EFF0[int(1.0*nokp),int(Nmax/2)]-MAT_BANDS_FLOQUET_EFF0[int(1.0*nokp),int(Nmax/2)-1])*1000    
    BW_EFF1[mm] = (MAT_BANDS_FLOQUET_EFF1[int(1.0*nokp),int(Nmax/2)]-MAT_BANDS_FLOQUET_EFF1[int(1.0*nokp),int(Nmax/2)-1])*1000
    BW_EFF2[mm] = (MAT_BANDS_FLOQUET_EFF2[int(1.0*nokp),int(Nmax/2)]-MAT_BANDS_FLOQUET_EFF2[int(1.0*nokp),int(Nmax/2)-1])*1000
    BW_EFF3[mm] = (MAT_BANDS_FLOQUET_EFF3[int(1.0*nokp),int(Nmax/2)]-MAT_BANDS_FLOQUET_EFF3[int(1.0*nokp),int(Nmax/2)-1])*1000
    BW_EFF4[mm] = (MAT_BANDS_FLOQUET_EFF4[int(1.0*nokp),int(Nmax/2)]-MAT_BANDS_FLOQUET_EFF4[int(1.0*nokp),int(Nmax/2)-1])*1000
    BW_EFF5[mm] = (MAT_BANDS_FLOQUET_EFF5[int(1.0*nokp),int(Nmax/2)]-MAT_BANDS_FLOQUET_EFF5[int(1.0*nokp),int(Nmax/2)-1])*1000
    BW_EFF6[mm] = (MAT_BANDS_FLOQUET_EFF6[int(1.0*nokp),int(Nmax/2)]-MAT_BANDS_FLOQUET_EFF6[int(1.0*nokp),int(Nmax/2)-1])*1000
    BW_EFF7[mm] = (MAT_BANDS_FLOQUET_EFF7[int(1.0*nokp),int(Nmax/2)]-MAT_BANDS_FLOQUET_EFF7[int(1.0*nokp),int(Nmax/2)-1])*1000
    BW_EFF8[mm] = (MAT_BANDS_FLOQUET_EFF8[int(1.0*nokp),int(Nmax/2)]-MAT_BANDS_FLOQUET_EFF8[int(1.0*nokp),int(Nmax/2)-1])*1000
    BW_EFF9[mm] = (MAT_BANDS_FLOQUET_EFF9[int(1.0*nokp),int(Nmax/2)]-MAT_BANDS_FLOQUET_EFF9[int(1.0*nokp),int(Nmax/2)-1])*1000
    
    print('Number of bands considered: '+str(np.size(MAT_BANDS[0,:])))    
    print("Effective bandwidth: "+str(BW[mm]-BW_GS))
    
    ## Measure gap
    FGAPS[mm] = (FLOQUET_BANDS_DOWN[0,int(m_max*Nmax+Nmax/2)]-FLOQUET_BANDS_DOWN[0,int(m_max*Nmax+Nmax/2-1)])*1000
    FGAPS_EFF[mm] = (MAT_BANDS_FLOQUET_EFF[0,int(Nmax/2)]-MAT_BANDS_FLOQUET_EFF[0,int(Nmax/2-1)])*1000
    FGAPS_EFF0[mm] = (MAT_BANDS_FLOQUET_EFF0[0,int(Nmax/2)]-MAT_BANDS_FLOQUET_EFF0[0,int(Nmax/2-1)])*1000
    FGAPS_EFF1[mm] = (MAT_BANDS_FLOQUET_EFF1[0,int(Nmax/2)]-MAT_BANDS_FLOQUET_EFF1[0,int(Nmax/2-1)])*1000
    FGAPS_EFF2[mm] = (MAT_BANDS_FLOQUET_EFF2[0,int(Nmax/2)]-MAT_BANDS_FLOQUET_EFF2[0,int(Nmax/2-1)])*1000
    FGAPS_EFF3[mm] = (MAT_BANDS_FLOQUET_EFF3[0,int(Nmax/2)]-MAT_BANDS_FLOQUET_EFF3[0,int(Nmax/2-1)])*1000
    FGAPS_EFF4[mm] = (MAT_BANDS_FLOQUET_EFF4[0,int(Nmax/2)]-MAT_BANDS_FLOQUET_EFF4[0,int(Nmax/2-1)])*1000
    FGAPS_EFF5[mm] = (MAT_BANDS_FLOQUET_EFF5[0,int(Nmax/2)]-MAT_BANDS_FLOQUET_EFF5[0,int(Nmax/2-1)])*1000
    FGAPS_EFF6[mm] = (MAT_BANDS_FLOQUET_EFF6[0,int(Nmax/2)]-MAT_BANDS_FLOQUET_EFF6[0,int(Nmax/2-1)])*1000
    FGAPS_EFF7[mm] = (MAT_BANDS_FLOQUET_EFF7[0,int(Nmax/2)]-MAT_BANDS_FLOQUET_EFF7[0,int(Nmax/2-1)])*1000
    FGAPS_EFF8[mm] = (MAT_BANDS_FLOQUET_EFF8[0,int(Nmax/2)]-MAT_BANDS_FLOQUET_EFF8[0,int(Nmax/2-1)])*1000   
    FGAPS_EFF9[mm] = (MAT_BANDS_FLOQUET_EFF9[0,int(Nmax/2)]-MAT_BANDS_FLOQUET_EFF9[0,int(Nmax/2-1)])*1000  
    
    print("Floquet band gap (exact): "+str(FGAPS[mm]))
    print("Floquet band gap (H_EFF): "+str(FGAPS_EFF[mm]))
    

DELTA_TRIVIAL = 2.*(0.2293780611802569*AA)**2/1.5*1000
## see EQ_BANDS for Fermi velocity
print(DELTA_TRIVIAL)

fig1 = plt.figure(101)
gs1 = gridspec.GridSpec(1, 2)    

############################################################################### 00

ax1 = fig1.add_subplot(gs1[0,0])
ax1.set_title('(a)', loc='left', y=1.1, pad=-14)
ax1.plot(AA, FGAPS_EFF7, color=C1, label=r'$\mathrm{\langle flat,m | H_{eff}^{A}| flat,n \rangle}$')
ax1.plot(AA, FGAPS, color='k', linestyle='--', linewidth=2.0, label=r'$\mathrm{H_{FLOQ}^{TRUNC}}$')
ax1.plot(AA,DELTA_TRIVIAL, color=C5, linewidth=2.0, label=r'$\mathrm{2(v_F \cdot A)^2/\Omega}$')
ax1.set_xlabel(r'$\mathrm{A_0}$ $\mathrm{(\AA^{-1})}$')
ax1.set_xticks([0.0, 0.005, 0.010, 0.015, 0.020, 0.025])
ax1.set_ylabel(r'$\mathrm{\Delta_{K_1}}$ $\mathrm{(meV)}$')
plt.legend(loc='upper left',frameon=False)

############################################################################### 01
 
# =============================================================================
# fit_lin1t= np.polyfit(np.log(AA), np.log(DELTA_TRIVIAL), 1)
# p_lin1t = np.poly1d(fit_lin1t)
# fit_lin1t = lambda x: np.exp(p_lin1t(np.log(x)))
# pdx_lin1t = np.polyder(p_lin1t,1)
# print("gap (trivial): "+str(p_lin1t))
# 
# fit_lin17= np.polyfit(np.log(AA), np.log(np.abs(FGAPS_EFF7)), 1)
# p_lin17 = np.poly1d(fit_lin17)
# fit_lin17 = lambda x: np.exp(p_lin17(np.log(x)))
# pdx_lin17 = np.polyder(p_lin17,1)
# print("gap (17): "+str(p_lin17))
# 
# ax12= fig1.add_subplot(gs1[0,1])
# ax12.loglog(AA, np.abs(DELTA_TRIVIAL), 'k')
# ax12.loglog(AA, np.abs(FGAPS_EFF7), 'r')
# ax12.loglog(AA, fit_lin1t(AA), "x", color='k',  linewidth=2.0, label='$\mathrm{fitted}$ $\mathrm{slope}:$ $\mathrm{2.0}$')
# ax12.loglog(AA, fit_lin17(AA), "x", color='red',  linewidth=2.0, label='$\mathrm{fitted}$ $\mathrm{slope}:$ $\mathrm{2.0}$')
# 
# ax12.set_ylabel(r'$|\mathrm{\Delta_{K_1}}|$ $\mathrm{(meV)}$')
# plt.legend(loc='lower right',frameon=False)
# =============================================================================
############################################################################### 10

ax2 = fig1.add_subplot(gs1[0,1])
ax2.set_title('(b)', loc='left', y=1.1, pad=-14)
ax2.plot(AA, BW_EFF-BW_GS, color=C1, label=r'$\mathrm{H_{eff}}$')
ax2.plot(AA, BW_EFF0-BW_GS, color=C2, label=r'$\mathrm{\mathrm{\langle m | H^{0}_{eff} | m \rangle}}$')
ax2.plot(AA, BW_EFF3, color=C3, label=r'$\mathrm{\mathrm{\langle m | H^{AA}_{eff} | m \rangle}}$')
ax2.plot(AA, BW_EFF7, color=C4, label=r'$\mathrm{\mathrm{\langle flat,m | H^{A}_{eff} |  flat,n \rangle}}$')
ax2.plot(AA, BW_EFF0+BW_EFF3+BW_EFF7-BW_GS, color=C5, label='SUM')
ax2.plot(AA, BW-BW_GS, color='k', linestyle='--', label=r'$\mathrm{H_{FLOQ}^{TRUNC}}$')
ax2.set_xlabel(r'$\mathrm{A_0}$ $\mathrm{(\AA^{-1})}$')
ax2.set_xticks([0.0, 0.005, 0.010, 0.015, 0.020, 0.025])
ax2.set_ylabel(r'$\mathrm{\Delta_\Gamma-\Delta_\Gamma^0}$ $\mathrm{(meV)}$')
plt.legend(loc='upper left', frameon=False)

############################################################################### 11

LIMIT = RANGE
# =============================================================================
# fit_lin0= np.polyfit(np.log(AA[0:LIMIT]), np.log(np.abs(BW[0:LIMIT]-BW_GS)), 1)
# p_lin0 = np.poly1d(fit_lin0)
# fit_lin0 = lambda x: np.exp(p_lin0(np.log(x)))
# pdx_lin0 = np.polyder(p_lin0,1)
# print("bandwidth (H0): "+str(p_lin0))
#    
# =============================================================================
# =============================================================================
# fit_lin3= np.polyfit(np.log(AA[0:LIMIT]), np.log(np.abs(BW_EFF3[0:LIMIT])), 1)
# p_lin3 = np.poly1d(fit_lin3[0:LIMIT])
# fit_lin3 = lambda x: np.exp(p_lin3(np.log(x)))
# pdx_lin3 = np.polyder(p_lin3,1)
# print("bandwidth (AA intra): "+str(p_lin3))
# 
# ax22= fig1.add_subplot(gs1[1,1])
# 
# #ax22.loglog(AA, np.abs(BW-BW_GS), 'b')
# ax22.loglog(AA, np.abs(BW_EFF3), 'b')
# #ax22.loglog(AA,fit_lin0(AA), "k", linestyle='dashed', linewidth=1.5, label='$\mathrm{fitted} \mathrm{slope}: \mathrm{2}$')
# ax22.loglog(AA,fit_lin3(AA), "x", color='b', linewidth=2.0, label='$\mathrm{fitted}$ $\mathrm{slope}:$ $\mathrm{2.0}$')
# 
# ax22.set_xlabel(r'$\mathrm{driving}$ $\mathrm{amplitude}$ $\mathrm{(\AA^{-1})}$')
# ax22.set_ylabel(r'$|\mathrm{\Delta_{\Gamma}}|$ $\mathrm{(meV)}$')
# plt.legend(loc='lower right',frameon=False)
# =============================================================================

plt.tight_layout()
plt.subplots_adjust(hspace=0.0, wspace=0.40)
plt.plot()
