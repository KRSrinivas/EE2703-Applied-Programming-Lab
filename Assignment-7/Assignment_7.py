import numpy as np
from sympy import *
import scipy.signal as sp
init_session
import pylab as p

def display(i,x,y,xlabel='t',ylabel='x',title='response'):

    p.figure(i)
    p.plot(x,y,'-r',label=r'$V_{o}$')
    p.xlabel(xlabel,fontsize=15)
    p.ylabel(ylabel,fontsize=15)
    p.title(title)
    p.legend(loc ='upper right')
    p.grid(True)
    p.show()

def lti_conversion(xpr, s=symbols('s')): # function to convert sympy polynomial to sp.lti polynomial

    num, denom = simplify(xpr).as_numer_denom()
    H = Poly(num, s).all_coeffs() ,Poly(denom, s).all_coeffs()
    l_num, l_den = [lambdify((), c)() for c in H]  # convert to floats
    return sp.lti(l_num, l_den)


def lowpass(R1,R2,C1,C2,G,Vi): # function to solve for output voltage in lowpass filter

    s = symbols('s')
    A = Matrix([[0,0,1,-1/G],[-1/(1+s*R2*C2),1,0,0],[0,-G,G,1],[-1/R1-1/R2-s*C1,1/R2,0,s*C1]])
    b = Matrix([0,0,0,-Vi/R1])
    V = A.inv()*b
    return (A,b,V)

def highpass(R1,R3,C1,C2,G,Vi): # function to solve for output voltage in highpass filter

    s = symbols('s')
    A = Matrix([[0,0,1,-1/G],[-1/(1+1/(s*R3*C2)),1,0,0],[0,-G,G,1],[-s*C1-s*C2-1/R1,s*C2,0,1/R1]])
    b = Matrix([0,0,0,-Vi*s*C1])
    V = A.inv()*b
    return (A,b,V)

# Obtaining the Magnitude response for the lowpass circuit
s = symbols('s')
A,b,V = lowpass(10000,10000,1e-9,1e-9,1.586,1)
Vo = V[3]
H = lti_conversion(Vo)
w = p.logspace(0,8,801)
ss = 1j*w
hf = lambdify(s,Vo,'numpy')
v = hf(ss)
p.loglog(w,abs(v),lw=2)
p.xlabel(r'$w\rightarrow$')
p.ylabel(r'$|H(jw)|\rightarrow$')
p.title('Magnitude response for Circuit 1')
p.grid(True)
p.show()

# Obtaining the step response for the lowpass circuit
t = np.linspace(0,0.001,1000)
Vo = sp.step(H,T=t)
display( 0, Vo[0], Vo[1], r't$\rightarrow$', r'$V_{o}\rightarrow$', 'Step Response for Circuit 1')

# Obtaining the response for mixed frequency sinusoid in lowpass filter
t = np.linspace(0,0.01,100000)
Vi = np.multiply((np.sin(2000*np.pi*t)+np.cos(2000000*np.pi*t)),np.heaviside(t,0.5))
Vo = sp.lsim(H,Vi,T=t)
p.figure(1)
p.plot(Vo[0],Vi,label=r'$V_{in}$')
display( 1, Vo[0], Vo[1], r't$\rightarrow$', r'$V\rightarrow$', 'Mixed Frequency Sinusoid Response for Circuit 1')

# Magnitude response of the high pass filter
A,b,V = highpass(10000,10000,1e-9,1e-9,1.586,1)
Vo = V[3]
H = lti_conversion(Vo)
w = p.logspace(0,8,801)
ss = 1j*w
hf = lambdify(s,Vo,'numpy')
v = hf(ss)
p.loglog(w,abs(v),lw=2)
p.xlabel(r'$w\rightarrow$')
p.ylabel(r'$|H(jw)|\rightarrow$')
p.title('Magnitude response of the Circuit 2')
p.grid(True)
p.show()

# Obtaining the response for damped sinusoids in highpass filters
t = np.linspace(0,10,1000)
Vi = np.multiply(np.multiply(np.exp(-0.5*t),np.sin(2*np.pi*t)),np.heaviside(t,0.5))
Vo = sp.lsim(H,Vi,T=t)
p.figure(2)
p.plot(Vo[0],Vi,label=r'$V_{in}$')
display(2,Vo[0],Vo[1],r't$\rightarrow$',r'$V\rightarrow$', 'Damped Sinusoids Response for Circuit 2')

#Obtaining the Step response for highpass filters
t = np.linspace(0,0.001,1000)
Vo = sp.step(H,T=t)
display(0,Vo[0],Vo[1],r't$\rightarrow$',r'$V_{o}\rightarrow$', 'Step Response for Circuit 2')
