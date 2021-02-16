'''
Student Name: K.R.Srinivas
Roll Number : EE18B136
EE2703 - Applied Programming and Lab Final Examination
'''
from pylab import *
import scipy
import matplotlib.pyplot as plt
import math
import time

# Defining Parameters

L_y   = 20     # The lengh of the capacitor along y direction
L_x   = 10     # The lengh of the capacitor along x direction
delta = 0.5    # the distance between nodes (assumed same along x and along y)
ratio = 0.5    # h/Ly ratio
No = 1500      # No. of iterations
accuracy = 1e-7# desired accuracy

epsilon_l = 2        # dielectric constant of liquid
epsilon_a = 1        # dielectric constant of air
epsilon   = 8.85e-12 #permitivity of freespace


M  = int(L_y/delta)   # No. of rows (the height given as the index k corresponding to h)
N  = int(L_x/delta)   # No. of coloumns (the number of nodes along y, including the boundary nodes)
k  = int(M-ratio*M)   # the height given as the index k corresponding to h.
#####################################################################################################################################
# Creating Meshgrid and Initialising Potential phi
m = np.linspace( 10, -10, num = M, dtype=float)
n = np.linspace( 5, -5, num = N, dtype=float)
phi = np.zeros( (M,N), dtype = float)           # matrix to hold the potential at different points in capacitor mesh
Y,X = meshgrid(n,m)                             # modeling capacitor as a mesh

#function to update potential
def phi_new(phi,phiold,t):
    phi[1:t,1:-1] = 0.25*(phiold[1:t,0:-2]+ phiold[1:t,2:]+ phiold[0:t-1,1:-1] + phiold[2:t+1,1:-1])           # to handle potential update in air region (using transformed difference equation)
    phi[t,1:-1] = (epsilon_a*phiold[t-1,1:-1] +epsilon_l*phiold[t+1,1:-1])/(epsilon_a+epsilon_l)               # to handle continuity of Dn at the surface (junction of air & liquid)
    phi[t+1:-1,1:-1] = 0.25*(phiold[t+1:-1,0:-2]+ phiold[t+1:-1,2:]+ phiold[t:-2,1:-1] + phiold[t+2:,1:-1])    # to handle potential update in liquid region (using transformed difference equation)
    return phi

# function to update potential at boundaries
def boundary(phi):
    phi[0,:]    = 1 # Top boundary
    phi[1:,N-1] = 0 # Right boundary
    phi[1:,0]   = 0 # Left boundary
    phi[M-1,:]  = 0 # Bottom boundary
    return phi

err = np.zeros(No)# initialise error array

for i in range(No):
    phiold = phi.copy()                     # copy the old phi
    phi = phi_new(phi,phiold,k)             # Updating the potential
    phi = boundary(phi)                     # applying boundary conditions
    err[i] = np.max(np.abs(phi-phiold))     # Appending errors for each iterations
    if(err[i] == 0):                        # command to break the loop if the condition is met
        print("Reached steady state at ",i," Iterations")
        break

# Electric field ;Electric field defined as the difference eqn of potential (Note:The difference eqn directly transforms to -dV/dx)
Ex = (1/2*(phi[1:-1,0:-2]-phi[1:-1,2:]))/(delta*1e-2) # electric field in x direction
Ey = (1/2*(phi[0:-2,1:-1]-phi[2:,1:-1]))/(delta*1e-2) # electric field in y direction

#plotting 2d contour of final potential & Vector plot of current flow
title("2D Contour plot of Potential & Vector plot of Electric Field")
xlabel("X")
ylabel("Y")
contourf(Y,X,phi)
quiver(Y[1:-1,1:-1],X[1:-1,1:-1],-Ex,-Ey)
colorbar()
show()

#plotting Error on semilog
title("Error - semilog plot")
xlabel("No of iterations")
ylabel("Error")
semilogy(range(No),err)
show()
#########################################################################################################################
#Q(e)
# Given the ratio h/Ly, we are supposed to calculate the Qtop and Qliquid (Q on the wall surface in contact
# with dielectric.) We approximate Gauss law to discrete domain to obtain the charge.
# Q_top = delta*N*epsilon_a*sum(Ey[0,:])
# Q_liq = delta*epsilon_l*(N*sum(Ey[M-1,:])+(M-k-1)*sum(Ex[k+1:,0])+(M-k-1)*sum(Ex[k+1:,N-1]))
ratio_list = np.array([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9])
Q_top = []
Q_liq = []

for i in ratio_list:                             # for every value of h/Ly ratio we compute the potential of the cap mesh and evaluate Qtop and Qliq
    height_index = int(M-i*M)

    m_ = np.linspace( 10, -10, num = M, dtype=float)
    n_ = np.linspace( 5, -5, num = N, dtype=float)
    phi_ = np.zeros( (M,N), dtype = float)
    Y1,X1 = meshgrid(n_,m_)

    err_ = np.zeros(No)

    for i in range(No):
        phiold_ = phi_.copy()                     # copy the old phi
        phi_ = boundary(phi_)                     # applying boundary conditions
        phi_ = phi_new(phi_,phiold_,height_index) # Updating the potential
        err_[i] = np.max(np.abs(phi_-phiold_))    # Appending errors for each iterations
        if(err_[i] < accuracy):
            break

    Ex_ = (1/2*(phi_[1:-1,0:-2]-phi_[1:-1,2:]))/(delta*1e-2)  # electric field in x direction
    Ey_ = (1/2*(phi_[0:-2,1:-1]-phi_[2:,1:-1]))/(delta*1e-2)  # electric field in y direction

    Q_top.append(delta*epsilon*epsilon_a*sum(np.abs(Ey_[0,:]))*1e-2)                                                               # Qtop defined by approximate Gauss law, computing the charge on top metal plate.
    Q_liq.append(-(delta*1e-2*epsilon*epsilon_l*(sum(np.abs(Ey_[M-3,:]))+sum(np.abs(Ex_[height_index:,0]))+sum(np.abs(Ex_[height_index:,N-3]))))) # Qliq defined by approximate Gauss law, specifically for bottom surface and sides
    # sum(np.abs(Ey_[M-3,:])): electric field components normal to bottom surface of liquid
    # sum(np.abs(Ex_[height_index:,0])): electric field component normal to left wall in contact with liquid
    # sum(np.abs(Ex_[height_index:,N-3])): electric field component normal to right wall in contact with liquid

plot(ratio_list, Q_top, 'go--', linewidth=2, markersize=8)        # plotting charge on top metal plate vs h/Ly ratio
xlabel('h/Ly ratio')
ylabel('Q top (C/m)')
title('Plot of Q top vs h/Ly ratio')
show()

plot(ratio_list, Q_liq, 'bo--', linewidth=2, markersize=8)        # plotting charge on the surface of wall in contact with dielectric vs h/Ly ratio
xlabel('h/Ly ratio')
ylabel('Q liquid (C/m)')
title('Plot of Q liquid vs h/Ly ratio')
show()
###################################################################################################################
#Q (f)

# For h = 0.5Ly.To compute Ex and Ey on the at (m+0:5;n+0:5), i.e., at the centre of mesh cells. We attempt to obtain
# electric field at (m,n) differently and then try approximating the electric field at the center of the mesh cell (m+0.5,n+0.5).

E_x = (phi[:,:-1]-phi[:,1:])/(delta*1e-2) # calculation of E field modified in this problem to find E field at mesh cell center
E_y = (phi[:-1,:]-phi[1:,:])/(delta*1e-2) # calculation of E field modified in this problem to find E field at mesh cell center

E_x1 = 1/2*(E_x[:-1,:]+E_x[1:,:]) # averaging out the electric field along the sides of mesh cell
E_y1 = 1/2*(E_y[:,:-1]+E_y[:,1:]) # averaging out the electric field along the sides of mesh cell

# Creating a mesh with centers of phi mesh
m = np.linspace( 10, -10, num = M-1, dtype=float)
n = np.linspace( 5, -5, num = N-1, dtype=float)
Y_,X_ = meshgrid(n,m)

#plotting Electric field vector at each mesh cell center.
title("Vector plot of Electric Field at Mesh-Cell center")
xlabel("X")
ylabel("Y")
quiver(Y_,X_,-E_x1,-E_y1)
show()

# Displacement vector D is given as D = epsilon*E; At the interface the normal component of D must be ideally
# conserved in electrostatic condition.(D_liq - D_air).n = 0 => epsilon_air*E_air = epsilon_liq*E_liq.
# It is quite accurate to compute the the D vector at the mesh cell center above and below k th level
# of phi matrix and verify their difference nad hence continuity. The electric field component at mesh cell centers
# above height h correspond to index k-1 and those below height h correspont to index k.

D_air = epsilon_a*E_y1[k-1,:] # Displacement vector in air medium above interface
D_liq = epsilon_l*E_y1[k,:]   # Displacement vector in liquid medium above interface

D_diff = (np.abs(D_air-D_liq)/D_air)*100

print("The percentage change in value of Displacement vector (normal component) across the interface (m=k) is: ")
print(D_diff)

print(" ")

########################################################################################################################
# Q (g)
# Verify Snell's law i this case; Verify if epsilon_a*sin(theta_i) = epsilon_l*sin(theta_r)
incidence = np.zeros(N-2)
refraction = np.zeros(N-2)
angle_change = np.zeros(N-2)

for i in range(N-2):
    incidence[i]= epsilon_a*np.sin(math.atan2(Ey[k-2,i],Ex[k-2,i]))                       # epsilon_a*sin(theta_i)
    refraction[i] = epsilon_l*np.sin(math.atan2(Ey[k,i],Ex[k,i]))                         # epsilon_l*sin(theta_r)
    angle_change[i] = (math.atan2(Ey[k-2,i],Ex[k-2,i])-math.atan2(Ey[k,i],Ex[k,i]))       # change in direction of electric field interms of angle_change

print("The change in direction of electric field in terms of angle deviation (radians):")
print(angle_change)
print(" ")

snell_difference_percentage = ((np.abs(incidence-refraction))/incidence)*100

print("Snell's Law validity (Percentage):")
print(snell_difference_percentage)  # looking at epsilon_a*sin(theta_i) and epsilon_l*sin(theta_r) its quite evident that snells law doesnt hold.
print(" ")
#########################################################################################################################
# Q(c)
# to verify the effectiveness of parallel computation using vectorization over iteration through loops, an example is presented
A = np.linspace( 100, -100, num = 100000, dtype=float)
B = np.zeros(100000)

t1 = time.time()
for i in range(100000):
    B[i] = 0.5*A[i]
t2 = time.time()

print("Time taken by loop implementation:",1e6*(t2-t1),"ms")

C = np.linspace( 100, -100, num = 100000, dtype=float)
D = np.zeros(100000)

t3 = time.time()
D[:] = (0.5)*C[:]
t4 = time.time()

print("Time taken by vectorised implementation:",1e6*(t4-t3),"ms")
