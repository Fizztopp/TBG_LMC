# -*- coding: utf-8 -*-
"""
This file plots the single-particle dispersions for chosen rescaling
and twist angle along the high symmetry path.
@author: toppgabr
"""

import matplotlib.pyplot as plt  
import numpy as np
import spglib
import matplotlib.gridspec as gridspec
import matplotlib as mpl
from mpl_toolkits.mplot3d import axes3d, Axes3D

mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['lines.linewidth'] = 3
mpl.rcParams['lines.markersize'] = 15
mpl.rcParams['font.size'] = 14  # <-- change fonsize globally
mpl.rcParams['legend.fontsize'] = 20
mpl.rcParams['axes.titlesize'] = 20
mpl.rcParams['axes.labelsize'] = 20
mpl.rcParams['xtick.major.size'] = 10
mpl.rcParams['ytick.major.size'] = 10
mpl.rcParams['xtick.major.width'] = 1
mpl.rcParams['ytick.major.width'] = 1
mpl.rcParams['xtick.direction'] = 'inout'
mpl.rcParams['ytick.direction'] = 'inout'
mpl.rcParams['figure.titlesize'] = 28
mpl.rcParams['figure.figsize'] = [10.,10]

### m and n define the scaled twist angle $\theta'$ used in the computations
m = 31
n = m+1 
## Number of kpoints per subpath
nokp = 96
Nmesh = 32

# number of lattice sites 
LS_real = 4*(m**2 + n**2 + m*n)
print("Number of atomic sites: "+str(LS_real))

# twist angle
theta = np.arccos(0.5*(m**2 + n**2 + 4*m*n)/(m**2 + n**2 + m*n))
print("Theta in rad: "+str(theta))                                                                    
print("Theta in degree: "+str(theta*360./(2.*np.pi)))

# unscaled twist angle $\theta$
desired_theta = 1.0501208797943464
print("Theta in degree (rescaled): "+str(desired_theta))   
desired_theta = 1.0501208797943464*np.pi/180.
print("Theta in rad (rescaled): "+str(desired_theta))    

# rescaling coefficient
Lambda = np.sin(theta/2.)/(np.sin(desired_theta/2.)) 
print("Rescaling coefficient: "+str(Lambda))

# lattice constant of monolayer graphene (m)
a0 = 2.46##*10.**-10.
print("Monolayer lattice constant in AA: "+str(a0))
a = a0*Lambda 
print("Monolayer lattice constant in AA (rescaled): "+str(a))

# NN bond length
b0 = a0/np.sqrt(3.)
print("Intralayer NN distance in AA: "+str(b0))
b = a/np.sqrt(3.)
print("Intralayer NN distance in AA (rescaled): "+str(b))
# determines whether one uses NN intralayer hopping only (i.e. NNN, NNNN etc intralayer hoppings are not included if NN_bool=1)
NN_bool = 1 

# interlayer distance (m)
d0 = 1.3618*a0 
print("Interlayer distance in AA: "+str(d0))
d = 1.3618*a 
print("Interlayer distance in AA (rescaled): "+str(d))

# NN hopping strength
tt = -2.7
print("Intralayer hopping in eV: "+str(tt))
t0 = tt/Lambda
print("Intralayer hopping in eV (rescaled): "+str(t0))

# interlayer hopping term
t1 = -0.11*tt
print("Interlayer hopping in eV: "+str(t1))
t1 = -0.11*tt*Lambda
print("Interlayer hopping in eV (rescaled): "+str(t1))

# determines the decay of the hopping terms
betat = 7.2 



###############################################################################

## monolayer lattice vectors:
a1 = a*np.array([1.,0.])
a2 = a*np.array([0.5,np.sqrt(3.)/2.])

a1_orig = a1
a2_orig = a2

## determine the Bravais vectors of superlattice
def R_mat_func(theta):
    return np.array([[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]])

a1_l1g = np.dot(R_mat_func(-theta/2.),a1)
a2_l1g = np.dot(R_mat_func(-theta/2.),a2)

a1_l2g = np.dot(R_mat_func(theta/2.),a1)
a2_l2g = np.dot(R_mat_func(theta/2.),a2)

L1 = m*a1_l1g + n*a2_l1g
L2 = np.dot(R_mat_func(np.pi/3.),L1)

L_mat1 = np.zeros((3))
L_mat2 = np.zeros((3))

L_mat1[0] = L1[0]
L_mat1[1] = L1[1]
L_mat1[2] = 3.

L_mat2[0] = L2[0]
L_mat2[1] = L2[1]
L_mat2[2] = 3.

# Set basis vectors of TBG
a1 = L1
a2 = L2

print(a1)
print(a2)

np.savetxt('L_VECS.dat', np.array([a1,a2]))                      

###############################################################################
## create coordinates of the atomic lattice
## saves them to variable "r_points_final":

Labs = np.linalg.norm(L1)
Nx = 100
Ny = Nx

nn = 50

r_points = np.zeros((2*Nx*Ny,4))
r_points2 = d*np.ones((2*Nx*Ny,4))

ii = 0
for i in range(1,Nx+1):
    for j in range(1,Ny+1):
        
        r_points[ii,0] = -(nn-1)*a1_orig[0] + i*a1_orig[0] -(nn-1)*a2_orig[0] + j*a2_orig[0]
        r_points[ii,1] = -(nn-1)*a1_orig[1] + i*a1_orig[1] -(nn-1)*a2_orig[1] + j*a2_orig[1]
        
        r_orig = np.array([r_points[ii,0],r_points[ii,1]])
        r_new2 = np.dot(R_mat_func(theta/2.),r_orig)
        r_new = np.dot(R_mat_func(-theta/2.),r_orig)
        
        r_points2[ii,0] = r_new2[0]
        r_points2[ii,1] = r_new2[1]
        r_points2[ii,3] = 0
        
        r_points[ii,0] = r_new[0]
        r_points[ii,1] = r_new[1]
        r_points[ii,3] = 0
        
        ii = ii + 1
        
        tau_b = (a1_orig - 2.*a2_orig)/3.
        
        r_points[ii,0] = -(nn-1)*a1_orig[0] + i*a1_orig[0] -(nn-1)*a2_orig[0] + j*a2_orig[0] + tau_b[0]
        r_points[ii,1] = -(nn-1)*a1_orig[1] + i*a1_orig[1] -(nn-1)*a2_orig[1] + j*a2_orig[1] + tau_b[1]
        
        r_orig = np.array([r_points[ii,0],r_points[ii,1]])
        r_new2 = np.dot(R_mat_func(theta/2.),r_orig)
        r_new = np.dot(R_mat_func(-theta/2.),r_orig)
        
        r_points2[ii,0] = r_new2[0]
        r_points2[ii,1] = r_new2[1]
        r_points2[ii,3] = 1
        
        r_points[ii,0] = r_new[0]
        r_points[ii,1] = r_new[1]
        r_points[ii,3] = 1
        
        ii = ii + 1

r_points_both = np.append(r_points,r_points2,axis=0)
r_points_final = np.array([])

for i in range(np.size(r_points_both[:,0])):
    temp = r_points_both[i,3]
    
    r_tmp = np.array([r_points_both[i,0],r_points_both[i,1],r_points_both[i,2]])
   
    n_z = np.array([0.,0.,1.])

    b1 = np.array([L1[0],L1[1],0.])
    b2 = np.array([L2[0],L2[1],0.])

    alpha = np.divide(np.dot(n_z,np.cross(b2,r_tmp)),np.dot(n_z,np.cross(b2,b1)))
    beta = np.divide(np.dot(n_z,np.cross(b1,r_tmp)),np.dot(n_z,np.cross(b1,b2)))
    
    r_tmp=np.append(r_tmp,temp)
    
    if(alpha <=0.5 and alpha>-0.5 and beta<=0.5 and beta>-0.5):
        r_points_final = np.append(r_points_final,r_tmp)
        
r_points_final = r_points_final.reshape(int(np.size(r_points_final)/4),4)         
N_FINAL = np.size(r_points_final[:,0])
print("Number of atomic sites (rescaled): "+str(N_FINAL))

## Save atomic positions to file
file = open('Unit_Cell.dat','w')
for i in range(N_FINAL):
    for j in range(4):
        file.write("%s " % r_points_final[i,j])
    file.write("\n")    
file.close() 

## PLOT Points
fig1 = plt.figure(1)
ax1 = fig1.add_subplot(111, projection='3d')
for i in range(int(np.size(r_points_final[:,0])/2)):
     if(r_points_final[i,3]>0.5):
         ax1.scatter(r_points_final[i,0], r_points_final[i,1], c="r", marker=".")
     else:
         ax1.scatter(r_points_final[i,0], r_points_final[i,1], c="b", marker=".")    
 
ax1.scatter(a1[0], a1[1], c="r", marker="x")
ax1.scatter(a2[0], a2[1], c="r", marker="x")
 
ax1.set_xlabel("x-axis")
ax1.set_ylabel("y-axis")
ax1.set_zlabel("z-axis")

plt.show()
 
###############################################################################

# Calculate spatial differences
R1 = np.zeros(N_FINAL)
R2 = np.zeros(N_FINAL)

for i in range(N_FINAL):
    R1[i] = r_points_final[i,0] #x
    R2[i] = r_points_final[i,1] #y


B = np.array([a1, a2])
B = 2.*np.pi*np.linalg.inv(B)

GM1 = B[:,0]
GM2 = B[:,1]

R1_diff = np.zeros((N_FINAL,N_FINAL))
R2_diff = np.zeros((N_FINAL,N_FINAL))

for i in range(N_FINAL):
    for j in range(N_FINAL):
        R1_diff[i,j] = -(R1[i] - R1[j])
        R2_diff[i,j] = -(R2[i] - R2[j])

r_cells = np.zeros((3,3,N_FINAL,3))

# Calulate atomic sites of neigbouring unit cells
for i in range(3):
    for j in range(3):
        r_vecs_tmp = np.zeros((N_FINAL,3))    
        L_shift = np.zeros((N_FINAL,2))
        L_shift[:,0] = L_shift[:,0] + (i-1)*L1[0] + (j-1)*L2[0]
        L_shift[:,1] = L_shift[:,1] + (i-1)*L1[1] + (j-1)*L2[1]
        
        ## Shift x,y compinent by Bravais basis vector
        r_vecs_tmp[:,:2] = r_points_final[:,:2] + L_shift
        r_vecs_tmp[:,2] = r_points_final[:,2]

        r_cells[i,j,:,:] = r_vecs_tmp

## PLOT atomic sites of super cell + NN cells
fig2 = plt.figure(2)
ax2 = fig2.add_subplot(111, projection='3d')

for i in range(3):
    for j in range(3):
        ax2.scatter(r_cells[i,j,:,0], r_cells[i,j,:,1], r_cells[i,j,:,2], c="r", marker=".")
        ax2.scatter(r_cells[1,1,:,0], r_cells[1,1,:,1], r_cells[1,1,:,2], c="b", marker=".")
ax2.set_zlim(-1*10**-8, 1*10**-8)
ax2.set_xlabel("x-axis")
ax2.set_ylabel("y-axis")
ax2.set_zlabel("z-axis")

plt.show()

#Calculate hopping matrix elements
###############################################################################

# spatial cut off for the hopping terms 
spatial_cut_off = 4.0*b

#M_inter = np.zeros((3,3,N_FINAL,N_FINAL))

# =============================================================================
# for i in range(N_FINAL):
#     
#     r_i = r_points_final[i,:]
#     for i1 in range(3):
#        for i2 in range(3):
#             r_diff = np.sqrt((r_cells[i1,i2,:,0]-r_i[0])**2.+(r_cells[i1,i2,:,1]-r_i[1])**2.+(r_cells[i1,i2,:,2]-r_i[2])**2.)
#             
#             r_inds = np.where(r_diff <= spatial_cut_off)
# 
#             for j in range(np.size(r_inds)):
#                 rx = r_cells[i1,i2,r_inds[0][j],0] - r_i[0]
#                 ry = r_cells[i1,i2,r_inds[0][j],1] - r_i[1]
#                 rz = r_cells[i1,i2,r_inds[0][j],2] - r_i[2]
#                 R_vec = np.array([rx,ry,rz])
#                 R = np.linalg.norm(R_vec)
# 
#                 if(NN_bool == 1 and rz == 0 and np.sqrt(rx**2. + ry**2.) >= a-a/10.):       
#                     continue
#                             
#                 if (R<a*10.**-5.):
#                     continue 
#         
#                 l = np.dot(R_vec,np.array([1.,0.,0.]))/R
#                 m = np.dot(R_vec,np.array([0.,1.,0.]))/R
#                 n = np.dot(R_vec,np.array([0.,0.,1.]))/R
#                 
#                 M_inter[i1,i2,i,r_inds[0][j]] = t0*np.exp(-betat*np.abs(R-b)/b)*((rx**2.+ry**2.)/R**2.) + t1*np.exp(-betat*np.abs(R-d)/b)*(rz**2./R**2.)
#     
#                 if np.isnan(M_inter[i1,i2,i,r_inds[0][j]]):
#                     print('PROBLEM: Created nan onject!!!')
# =============================================================================


#TList=np.array([[0 0],[1 0],[-1,0],[0 1],[0,-1],[1,1],[-1,1],[-1,-1],[1,-1]],[M_inter[1,1,:,:],M_inter[2,1,:,:],M_inter[0,1,:,:],M_inter[1,2,:,:],M_inter[1,0,:,:],M_inter[2,2,:,:],M_inter[0,2,:,:],M_inter[0,0,:,:],M_inter[2,0,:,:]])

### Create k-points of high-symmetry path and BZ
###############################################################################

# Here the k-vectors are set:
K1_plus = np.dot(R_mat_func(theta/2.),np.array([4.*np.pi/a,0])/3.)
K2_plus = np.dot(R_mat_func(-theta/2.),np.array([4.*np.pi/a,0])/3.)

Kdiff_plus = -K2_plus + K1_plus
Kdiff_plus2 = Kdiff_plus/2.

q0 = K2_plus + Kdiff_plus2

k_gamma_to_Kbar = np.linalg.norm(Kdiff_plus)*np.array([np.cos(30.*np.pi/180.),np.sin(30.*np.pi/180.)])
k_gamma_to_M = np.linalg.norm(Kdiff_plus)*np.array([np.cos(30.*np.pi/180.),0])
k_M_to_K = -Kdiff_plus/2.;

Delta_k1 = -k_gamma_to_Kbar/nokp
Delta_k2 = k_gamma_to_M/nokp
Delta_k3 = k_M_to_K/nokp

Gamma_point = K1_plus - k_gamma_to_Kbar

K_PATH = np.zeros((3*nokp,3))

# =============================================================================
# for i in range(nokp):
#    K_PATH[i,0] = i*2*Delta_k1[0]/10 + K1_plus[0] - Gamma_point[0] - nokp*Delta_k1[0]/10
#    K_PATH[i,1] = i*2*Delta_k1[1]/10 + K1_plus[1] - Gamma_point[1] - nokp*Delta_k1[1]/10
#    
#    K_PATH[i+nokp,0] = -(i*2*Delta_k1[0]/10 + K1_plus[0] - Gamma_point[0] - nokp*Delta_k1[0]/10)
#    K_PATH[i+nokp,1] = -(i*2*Delta_k1[1]/10 + K1_plus[1] - Gamma_point[1] - nokp*Delta_k1[1]/10)
# =============================================================================
       

for i in range(nokp):
   K_PATH[i,0] = i*Delta_k1[0] + K1_plus[0] - Gamma_point[0]
   K_PATH[i,1] = i*Delta_k1[1] + K1_plus[1] - Gamma_point[1]
   
   K_PATH[i+nokp,0] = i*Delta_k2[0] + Gamma_point[0] - Gamma_point[0]
   K_PATH[i+nokp,1] = i*Delta_k2[1] + Gamma_point[1] - Gamma_point[1]
       
   K_PATH[i+2*nokp,0] = i*Delta_k3[0] + q0[0] - Gamma_point[0]
   K_PATH[i+2*nokp,1] = i*Delta_k3[1] + q0[1] - Gamma_point[1]
   
   
#K_PATH = np.flipud(K_PATH);
N_PATH = np.size(K_PATH[:,0])
print("Number of k-points of PATH: "+str(N_PATH))

file = open('k_path_full.dat','w')
for i in range(N_PATH):
    for j in range(3):
        file.write("%s " % K_PATH[i,j])
    file.write("\n")    
file.close()

K_POINTS = np.zeros((2,3))
K_POINTS[0,0] = K_PATH[0,0]
K_POINTS[0,1] = K_PATH[0,1]

K_POINTS[1,0] = K_PATH[nokp,0]
K_POINTS[1,1] = K_PATH[nokp,1]

file = open('k_path.dat','w')
for i in range(2):
    for j in range(3):
        file.write("%s " % K_POINTS[i,j])
    file.write("\n")    
file.close()
#############################################################Set k-mesh by hand
# Real space basis vectors of super cell in m
A1 = np.array([a1[0], a1[1], 0.])                                              
A2 = np.array([a2[0], a2[1], 0.])                                               
A3 = np.array([0.0, 0.0, 1.0])*d                                               

print(np.linalg.norm(A1))
print(np.linalg.norm(A2))

# Reciprocal space basis vectors of super cell in 1/m
B1 = 2.*np.pi*np.cross(A2,A3)/np.dot(A1,np.cross(A2,A3))                     
B2 = 2.*np.pi*np.cross(A3,A1)/np.dot(A2,np.cross(A3,A1))                     
B3 = 2.*np.pi*np.cross(A1,A2)/np.dot(A3,np.cross(A1,A2))

A_BZ = np.linalg.norm(np.cross(B1, B2))
N_B1 = np.linalg.norm(B1)
N_B2 = np.linalg.norm(B2)
dk = N_B1/(Nmesh-1)

## Calculate PATCH
MAT_BZ = np.zeros((Nmesh*Nmesh,3)) 
k0 = np.array([0.,0.,0.])-(B1+B2)/2.

for i in range(Nmesh):
    for j in range(Nmesh):
        MAT_BZ[i+Nmesh*j,:] = k0+i*dk*B1/N_B1+j*dk*B2/N_B2

file = open('k_BZ.dat','w')
file1 = open('k_weights.dat','w')    
for i in range(Nmesh*Nmesh):
    file1.write("%s " % 1.0)
    file1.write("\n")  
    for j in range(3):
        file.write("%s " % MAT_BZ[i,j])
    file.write("\n")    
file.close()
file1.close()

N_BZ = np.size(MAT_BZ[:,0])
print("Number of k-points of BZ: "+str(N_BZ))

# =============================================================================
# ## PLOT Points
# fig3 = plt.figure(3)
# ax3 = fig3.add_subplot(111, projection='3d')
# ax3.scatter(K_PATH[:,0], K_PATH[:,1], K_PATH[:,2], c="r", marker="o")
# ax3.scatter(MAT_BZ[:,0], MAT_BZ[:,1], MAT_BZ[:,2], c="k", marker=".", label="BZ")
# ax3.set_xlabel("x-axis")
# ax3.set_ylabel("y-axis")
# 
# plt.show()
# =============================================================================


#### Automatic with weights

A1_AA = np.array([a1[0], a1[1], 0.])#*10.**10.                                            
A2_AA = np.array([a2[0], a2[1], 0.])#*10.**10.                                                
A3_AA = np.array([0.0, 0.0, 1.0])*d#*10.**10.   

def Calc_BZ():    
    '''
    Calculates the k-vectors of the irreducable/full BZ
    '''
    
    mesh = [Nmesh, Nmesh, 1]    
    lattice = np.array([A1_AA, A2_AA, A3_AA])                                                    
    positions = [[0.0, 0.0, 0.0]]
    numbers= [1]             
             
    cell = (lattice, positions, numbers)
    
    print('spacegroup: ' +str(spglib.get_spacegroup(cell, symprec=1e-10)))
    
    # caclulatio of irr. BZ vectors + weights
    mapping, grid = spglib.get_ir_reciprocal_mesh(mesh, cell, is_shift=[0.,0.,0.])
    MAT_help = grid[np.unique(mapping)]/np.array(mesh, dtype=float)
    MAT_irr_BZ = np.zeros((np.size(MAT_help[:,0]),3))       
    for k in range(1,np.size(MAT_help[:,0])):
        MAT_irr_BZ[k,:] = B1*MAT_help[k,0] + B2*MAT_help[k,1] + B3*MAT_help[k,2] # transform from fractional to cartesian coordinates

    print("Number of kpoints: %d (irr BZ)" % len(np.unique(mapping)))
    num_kpoints = np.size(MAT_irr_BZ[:,0])

    weights = (np.unique(mapping,return_counts=True)[1])
    print("Number of kpoints: %d (full BZ, check of weights)" % weights.sum())                    
    
    # caclulatio of full BZ vectors (weights = 1) 
    MAT_BZ_full = np.array(grid, dtype=float)
    for k in range(1,np.size(MAT_BZ_full[:,0])):
        MAT_BZ_full[k,:] = B1*MAT_BZ_full[k,0] + B2*MAT_BZ_full[k,1] + B3*MAT_BZ_full[k,2]
    print("Number of kpoints: %d (full BZ)" % np.size(MAT_BZ_full[:,0]))
    
    file = open('k_BZ_irr.dat','w')
    for i in range(num_kpoints):
        for j in range(3):
            file.write("%s " % MAT_irr_BZ[i][j])
        file.write("\n")    
    file.close()
    
    file = open('k_weights_irr.dat','w')
    for i in range(num_kpoints):
        file.write("%s " % (weights[i]*1.0))
        file.write("\n")    
    file.close()
    file = open('k_BZ_full.dat','w')
    for i in range(np.size(MAT_BZ_full[:,0])):
        for j in range(3):
            file.write("%s " % (MAT_BZ_full[i][j]/mesh[0]))
        file.write("\n")    
    file.close()
    
    file = open('k_weights_full.dat','w')
    for i in range(np.size(MAT_BZ_full[:,0])):
        file.write("%s " % 1.0)
        file.write("\n")    
    file.close()
    return MAT_irr_BZ, MAT_BZ_full/(mesh[0]*1.0)                 # in 1/(AA)   

# =============================================================================
# MAT_irr_BZ, MAT_BZ_full = Calc_BZ()
# 
# ## PLOT Points
# fig3 = plt.figure(3)
# ax3 = fig3.add_subplot(111, projection='3d')
# ax3.scatter(K_POINTS[:,0], K_POINTS[:,1], K_POINTS[:,2], c="r", marker="x", s=50)
# ax3.scatter(K_PATH[:,0], K_PATH[:,1], K_PATH[:,2], c="r", marker="o")
# ax3.scatter(K_PATH[0,0], K_PATH[0,1], K_PATH[0,2], c="b", marker="x")
# ax3.scatter(K_PATH[nokp-1 ,0], K_PATH[nokp-1 ,1], K_PATH[nokp-1 ,2], c="b", marker="x")
# ax3.scatter(MAT_BZ_full[:,0], MAT_BZ_full[:,1], MAT_BZ_full[:,2], c="k", marker="o", label="BZ")
# ax3.scatter(MAT_irr_BZ[:,0], MAT_irr_BZ[:,1], MAT_irr_BZ[:,2], c="g", marker="o", label="irr. BZ")
# ax3.set_xlabel("x-axis")
# ax3.set_ylabel("y-axis")
# 
# plt.show()
# 
# print("Deltak: "+str(np.linalg.norm(K_POINTS[1,:]-K_PATH[0,:]))) 
# =============================================================================
