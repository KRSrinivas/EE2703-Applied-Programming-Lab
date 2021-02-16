######################
# K.R.Srinivas EE18B136
# Assignment-3
######################
import sys
from pylab import *
#import numpy as np
import scipy.special as sp
import math

A = 1.05
B = -0.105
# True Data Generation
def g (t, A = 1.05, B = -0.105):
   return(A*sp.jv(2, t)+B*t)
#data loading
try:
    data = loadtxt("fitting.dat")
except Exception:
    print("Data file unavailable!")
    exit()

sigma = np.logspace(-1,-3,9) # Standard Deviation of Noise
#plotting corrupt data
for i in range(1, 10):
        plot(data[:, 0], data[:, i], label=r'$\sigma_'+str(i)+'$ = ' + str(np.around(sigma[i-1], 4)))
        legend()
        title('Data corresponding to different added noise ')

# Plotting true data
plot(data[:, 0], g(data[:, 0]), '#000000', label="True Value")
legend()
grid()
xlabel("t")
ylabel("f(t)+ noise")
title("True Data Vs Corrupt Data ")
show()

#Error Bar Plot at every 5th element
errorbar(data[::5, 0], data[::5, 1], sigma[0], fmt="ro", label="Error Bar")
legend()
grid()
plot(data[:, 0], g(data[:, 0]), '#000000', label="f(t)")
legend()

# Annotation
annotate("Corrupt Data",(data[15, 0], data[15, 1]), xytext=(40, -40), textcoords="offset points", arrowprops={"arrowstyle": "-|>"})
annotate("True Data",(data[3, 0], g(data[3, 0])-0.01), xytext=(-20, 35), textcoords="offset points", arrowprops={"arrowstyle": "-|>"})
xlabel("t")
ylabel("Corrupt Data")
title("Error vs True Value")
show()

M = c_[sp.jv(2, data[:, 0]), data[:, 0]]# matrix holding J(t) & t
coefficient= [A, B] #coefficient matrix
G = dot(M, coefficient)
G1 = np.array(g(data[:, 0],A,B))
print(np.array_equal(G,G1))

A_ = linspace(0, 2, 21) #declaring different values of A
B_ = linspace(-0.2, 0, 21)
error = zeros((len(A_), len(B_)))

# Error calculation to plot contours
for i in range(len(A_)):
    for j in range(len(B_)):
        error[i, j] = mean(square(g(data[:, 0], A_[i], B_[j])-g(data[:, 0])))

# Plotting contours
CS = contour(A_, B_, error, levels=20)
xlabel("A")
ylabel("B")
title(r"Contour Plot of $\epsilon_{ij}$")
clabel(CS, inline=1, fontsize=10)
a = np.unravel_index(np.argmin(error[:,:]),error.shape)# extracting index corresponding to min error
plot(A_[a[0]],B_[a[1]],"ro") # plotting the coefficient corresponding to min error
grid()
annotate(" Minima", (A_[a[0]],B_[a[1]]), xytext=(-50, -30), textcoords="offset points", arrowprops={"arrowstyle": "-|>"})
show()
#estimating A and B using least square method
E = [lstsq(M, data[:, i], rcond = -1) [0] for i in range(1,10)] # estimate of A and B
#print(E)
E = np.asarray(E)
#print(E)
A_E = square(abs(E[:,0]-np.ones(9)*A)) #squared error in A
B_E = square(abs(E[:,1]-np.ones(9)*B)) #squared error in B
#plot error vs sigma
plot(sigma, A_E, '+--',linewidth=0.4, label='A_E', dashes=(8, 10))
plot(sigma, B_E, '+--',linewidth=0.4, label = 'B_E' ,dashes=(8, 10))
ylabel('MS Error',fontsize=15)
xlabel(r'$\sigma_{n}$',fontsize=15)
legend()
title('Variation of Error with Noise')
grid()
show()
#plot loglog scale of error with sigma
loglog(sigma,A_E,'bo',label = 'A_E')
loglog(sigma,B_E,'ro',label = 'B_E')
legend()
errorbar(sigma, A_E, std(A_E), fmt='bo')
errorbar(sigma, B_E, std(B_E), fmt='ro')
ylabel('MS Error',fontsize=15)
xlabel(r'$\sigma_{n}$',fontsize=15)
title('Variation of Error with Noise')
grid()
show()
