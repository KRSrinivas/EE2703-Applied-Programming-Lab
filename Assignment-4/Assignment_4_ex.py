# Importing libraries
from pylab import *
from scipy.integrate import quad
import numpy as np

def func(x):# defining the function
    return np.exp(x)

x = linspace(-2*pi, 4*pi, 400)
# Period of function created using fourier coefficients will be 2pi
period = 2*pi

# Plotting original function vs expected function for exp(x) in semilog
semilogy(x, func(x), 'k', label="Original Function")
semilogy(x, func(x % period), '--',label="Expected Function from fourier series")
legend(loc='upper right')
ylabel(r'$e^{x}\rightarrow$',fontsize=15)
xlabel(r'x$\rightarrow$',fontsize=15)
plt.show()
# defining the integrands in fourier coefficient calculation
def integrand_a(x, k, f):
    return f(x)*cos(k*x)

def integrand_b(x, k, f):
    return f(x)*sin(k*x)

#Calculating & defining vector to store the fourier coefficients
def coef_vec(f):
    coeff = []
    coeff.append((quad(f, 0, 2*pi)[0])/(2*pi))
    for i in range(1, 26):
        coeff.append((quad(integrand_a, 0, 2*pi, args=(i, f))[0])/pi)
        coeff.append((quad(integrand_b, 0, 2*pi, args=(i, f))[0])/pi)

    return coeff
# defining the matrix of harmonics
def argument_matrix(nrow, ncol, x):
    A = zeros((nrow, ncol))  # allocate space for A
    A[:, 0] = 1  # col 1 is all ones
    for k in range(1, int((ncol+1)/2)):
        A[:, 2*k-1] = cos(k*x)  # cos(kx) column
        A[:, 2*k] = sin(k*x)  # sin(kx) column
    return A

# Function to compute function from coefficients with argument as coefficient vector 'c'
def compute_fn(c):
    A = argument_matrix(400, 51, x)
    f_fourier =  np.dot(A,c)
    return f_fourier


# Initialising empty lists to store coefficients for both functions
exp_coeff = []
exp_coeff1 = []

exp_coeff1 = coef_vec(func)
# to store absolute value of coefficients
exp_coeff = np.abs(exp_coeff1)

# Computing function using fourier coeff
exp_fn_fourier = compute_fn(exp_coeff1)

# Plotting the Function computed using Fourier Coefficients
semilogy(x, exp_fn_fourier, 'ro', label="Function using Fourier Coefficients")
legend()
xlabel("x")
ylabel("y= e^x")
grid()
show()
#Plotting the Fourier Coefficients
semilogy((exp_coeff[1::2]), 'ro', label=r"$a_{n}$ using Integration")
semilogy((exp_coeff[2::2]), 'go', label=r"$b_{n}$ using Integration")
legend()
title(" Fourier coefficients of $e^{x}$ (semi-log)")
xlabel("n")
ylabel("Magnitude of coeffients")
grid()
show()


loglog((exp_coeff[1::2]), 'ro', label=r"$a_{n}$ using Integration")
loglog((exp_coeff[2::2]), 'ko', label=r"$b_{n}$ using Integration")
legend(loc='upper right')
title(" Fourier coefficients of $e^{x}$ (Log-Log)")
xlabel("n")
ylabel("Magnitude of coeffients")
grid()
show()

# Using lstsq method to evaluate coefficients
x_ = linspace(0, 2*pi, 400)
b = [func(i) for i in x_]
A = argument_matrix(400, 51, x_)
c = lstsq(A, b, rcond=-1)[0]
# To plot magnitude of coefficients
coeff_exp = np.abs(c)

semilogy((coeff_exp[1::2]), 'go', label=r"$a_{n}$ using Least Squares")
semilogy((coeff_exp[2::2]), 'ro', label=r"$b_{n}$ using Least Squares")
xlabel("x")
ylabel("Coefficient")
title("SemiLog Plot of Coefficients")
legend(loc='upper right')
grid()
show()

# Compute the deviation of lstsq coefficients from integral coefficient
def coeff_compare():
    deviations = []
    max_dev = 0
    deviations = np.abs(exp_coeff1 - c)
    max_dev = np.amax(deviations)
    return deviations, max_dev

dev1, maxdev1 = coeff_compare()

print("Maximum deviation in exp coefficients : ", maxdev1)

x1 = linspace(0, 2*pi, 400)

# Reconstruct the function from coefficients calculated by lstsq method
def fn_create_lstsq(c):
    f_lstsq = []
    A = argument_matrix(400, 51, x1)
    f_lstsq = np.dot(A,c)
    return f_lstsq

exp_fn_lstsq = fn_create_lstsq(c)

semilogy(x1, exp_fn_lstsq, 'go',label="Inverse Fourier Transform From Least Squares")
legend()
xlabel(r'n$\rightarrow$',fontsize=15)
ylabel(r'$f(x)\rightarrow$',fontsize=15)
grid()
show()
