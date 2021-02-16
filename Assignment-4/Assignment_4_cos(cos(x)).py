#Importing Libraries
import scipy.integrate as integrate
import numpy as np
import math
import matplotlib.pyplot as plt
import scipy.special as special
from pylab import *

def func(x): # defining the function
    return np.cos(np.cos(x))

x = linspace(-2*pi, 4*pi, 400)
# Period of function created using fourier coefficients will be 2pi
period = 2*pi
# Plotting original function vs expected function for cos(cos(x)) in semilog
semilogy(x, func(x), 'k', label="Original Function")
semilogy(x, func(x % period), '--',label="Expected Function from fourier series")
legend(loc='upper right')
ylabel(r'$\cos(\cos(x))\rightarrow$',fontsize=15)
xlabel(r'x$\rightarrow$',fontsize=15)
plt.show()

# defining the integrands in fourier coefficient calculation
def a_dc():
    return (1/(2*math.pi))*(integrate.quad( func, 0, 2*math.pi )[0])

def integrand_a(n,x):
    return (func(x))*math.cos(n*x)

def integrand_b(n,x):
    return (func(x))*math.sin(n*x)

def A_n(n):
    return (1/math.pi)*(integrate.quad(lambda x: integrand_a(n,x), 0, 2*math.pi)[0])

def B_n(n):
    return (1/math.pi)*(integrate.quad(lambda x: integrand_b(n,x), 0, 2*math.pi)[0])
#Calculating & defining vector to store the fourier coefficients
def coef_vec():
    coefficient = []
    coefficient.append(a_dc())
    for i in range(1,26):
        coefficient.append(A_n(i))
        coefficient.append(B_n(i))
    return coefficient
# defining the matrix of harmonics
def argument_matrix(nrow, ncol, x):
    A = zeros((nrow, ncol))
    A[:, 0] = 1
    for k in range(1, int((ncol+1)/2)):
        A[:, 2*k-1] = cos(k*x)
        A[:, 2*k] = sin(k*x)
    return A

x_ = linspace(-2*pi, 4*pi, 400)
# Function to compute function from coefficients with argument as coefficient vector 'c'
def compute_fn(coef):
    A = argument_matrix(400, 51, x_)
    f_fourier = np.dot(A,coef)
    return f_fourier

func_coef = np.array(coef_vec())
func_fourier = compute_fn(func_coef)
# Plotting the Function computed using Fourier Coefficients
plot(x_, func_fourier, 'ro', label="Function using Fourier Coefficients")
legend(loc='upper right')
xlabel("x")
ylabel("y=cos(cos(x))")
grid()
show()

index = [i for i in range(51)]
#Plotting the Fourier Coefficients
plt.figure()
semilogy((func_coef[1::2]), 'ro', label=r"$a_{n}$ using Integration")
semilogy((func_coef[2::2]), 'go', label=r"$b_{n}$ using Integration")
legend(loc='upper right')
plt.xlabel("x")
plt.ylabel("Coefficient")
plt.title("Fourier coefficients of $cos(cos(x))$ (semi-log)")

grid()
show()

plt.figure()
loglog((func_coef[1::2]), 'ro', label=r"$a_{n}$ using Integration")
loglog((func_coef[2::2]), 'go', label=r"$b_{n}$ using Integration")
legend(loc='upper right')
plt.xlabel("x")
plt.ylabel("Coefficient")
plt.title("Fourier coefficients of $cos(cos(x))$ (log-log)")

grid()
show()
# Using lstsq method to evaluate coefficients
b = [func(i) for i in x_]
A_ = argument_matrix(400, 51, x_)
c = lstsq(A_, b, rcond=None)[0]
index = [i for i in range(51)]

semilogy((c[1::2]), 'go', label=r"$a_{n}$ using Least Squares")
semilogy((c[2::2]), 'ro', label=r"$b_{n}$ using Least Squares")
xlabel("x")
ylabel("Coefficient")
title("SemiLog Plot of Coefficients")
legend(loc='upper right')
grid()
show()

# Compute the deviation of lstsq coefficients from integral coefficient
deviation = np.abs(coef_vec() - c)
max_deviation = np.amax(deviation)
print(max_deviation)

A = argument_matrix(400, 51, x_)
f_approx = np.dot(A,c)
# Reconstruct the function from coefficients calculated by lstsq method
plot(x_,f_approx,'go',label="lstsq")
legend(loc='upper right')

x_ = linspace(-2*pi, 4*pi, 400)
ylist = list(map(lambda y: func(y), x_))
plt.plot(x_, ylist,label="True Value")
xlabel(r'n$\rightarrow$',fontsize=15)
ylabel(r'$f(x)\rightarrow$',fontsize=15)
legend(loc='upper right')
plt.show()
