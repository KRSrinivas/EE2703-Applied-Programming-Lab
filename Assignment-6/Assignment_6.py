from pylab import *
import scipy.signal as sp

def func(t, w, alpha):
    return cos(w*t)*exp(-alpha*t)

def input(t):
    return cos(1e3*t)-cos(1e6*t)

#function to solve laplace domain equation and generate transient response
def time_response( alpha, wo, k):
    Dr = polymul( [1,0,pow(k,2)],[1,2*alpha,(pow(wo,2)+pow(alpha,2))])
    Nr = poly1d( [1,alpha] )
    H = sp.lti(Nr, Dr)
    t, x = sp.impulse(H, None, linspace(0, 100, 10000))
    return H, t, x

X, t1, x1 = time_response(0.5, 1.5, 1.5)
X, t2, x2 = time_response(0.05, 1.5, 1.5)

# plot of x(t) with decay of 0.5
plot(t1, x1, label="decay = 0.5")
legend()
title("Time Response of System-1")
xlabel(r"$t \to $")
ylabel(r"$x(t) \to $")
grid()
show()

# plot of x(t) with decay of 0.05
plot(t2, x2, label="decay = 0.05")
legend()
title("Time Response of System-2")
xlabel(r"$t \to $")
ylabel(r"$x(t) \to $")
grid()
show()

# frequency response (Magnitude and Phase) of the system for the given input
fig = plt.figure()

ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)

w,S,phi=X.bode()

ax1.semilogx(w,S)
ax1.set_ylabel('YdB-Magnitude')
ax1.set_xlabel(r'Frequency(rad\sec)')
ax1.set_title('Bode Plot of the Transfer Function')

ax2.semilogx(w,phi)
ax2.set_ylabel('Phase')
ax2.set_xlabel(r'Frequency(rad\sec)')

show()

# Plot of x(t) with different input frequencies
for wo in arange(1.4, 1.6, 0.05):
    X, t1, x1 = time_response(0.5, wo, 1.5)
    t = linspace(0, 200, 20000)
    t, y, svec = sp.lsim(X, func(t, wo, 0.05), t)
    plot(t, y, label="$w$ = %g rad/s" % (wo))
    legend()

xlabel(r"$t \to $")
ylabel(r"$x(t) \to $")
title("Varying frequency response")
grid()
show()

#coupled Equation solver
X_s  = sp.lti([1,0,2],[1,0,3,0])# Computes the impulse response of the transfer function
t, x = sp.impulse(X_s, None, linspace(0, 100, 10000))
Y_s = sp.lti([2],[1,0,3,0])# Computes the impulse response of the transfer function
t1, x1 = sp.impulse(Y_s, None, linspace(0, 100, 10000))


plot(t, x, label='x(t)')
legend()
plot(t1, x1, label='y(t)')
legend()
title("Time Response of X(s) & Y(s)")
xlabel(r"$t \to $")
ylabel(r"$x(t) \to $")
grid()
show()

# RLC 2-Port network solver
R = 100
C = 1e-6
L = 1e-6

H_s = sp.lti( [1], [L*C,R*C,1] )#transfer funtion of the given RLC network

fig = plt.figure()

ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)

w,S,phi=H_s.bode()# Calculates magnitude and phase response

ax1.semilogx(w,S)
ax1.set_ylabel('HdB-Magnitude')
ax1.set_xlabel(r'Frequency(rad\sec)')
ax1.set_title('Bode Plot of the Transfer Function')

ax2.semilogx(w,phi)
ax2.set_ylabel('Phase')
ax2.set_xlabel(r'Frequency(rad\sec)')

show()
# Plot of Vo(t) for 0<t<30usec
t = linspace( 0, 30e-6, 10000)
t, output, svec = sp.lsim(H_s, input(t), t)
plot( t, output)
xlabel(r"$t \to $")
ylabel(r"$y(t) \to $")
grid()
show()
# plot of Vo(t) for large time i.e at steady state
t1 = linspace( 0, 30e-3, 10000)
t1, output, svec = sp.lsim(H_s, input(t), t1)
plot( t1, output)
xlabel(r"$t \to $")
ylabel(r"$y(t) \to $")
grid()
show()
