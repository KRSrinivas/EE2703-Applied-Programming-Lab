from pylab import *
import numpy as np

#*******************************************************************************
# Finding FFT of the function sin(5x)
x = linspace(0,4*pi/5,1025); x=x[:-1]
y = sin(5*x)
Y = fftshift(fft(y))/1024.0
w = linspace(-512,511,1024)

figure()

subplot(2,1,1)
plot(w,abs(Y),lw=2)
xlim([-10,10])
ylabel(r"$|Y|$",size=16)
title(r"Spectrum of $\sin(5t)$")
grid(True)

subplot(2,1,2)

ii=where(abs(Y)>1e-3)
plot(w[ii],angle(Y[ii]),'go',lw=2)
xlim([-10,10])
ylabel(r"Phase of $Y$",size=16)
xlabel(r"$k$",size=16)
grid(True)
show()

#*******************************************************************************
# Finding FFT of the function (1+0.1cos(t))cos(10t)
t=linspace(-4*pi,4*pi,513);t=t[:-1]
y=(1+0.1*cos(t))*cos(10*t)
Y=fftshift(fft(y))/512.0
w=linspace(-64,64,513);w=w[:-1]

figure()

subplot(2,1,1)
plot(w,abs(Y),lw=2)
xlim([-15,15])
ylabel(r"$|Y|$",size=16)
title(r"Spectrum of $\left(1+0.1\cos\left(t\right)\right)\cos\left(10t\right)$")
grid(True)

subplot(2,1,2)

ii=where(abs(Y)>1e-3)
plot(w[ii],angle(Y[ii]),'go',lw=2)
xlim([-15,15])
ylabel(r"Phase of $Y$",size=16)
xlabel(r"$\omega$",size=16)
grid(True)
show()

#*******************************************************************************
# Finding FFT of the function sin^3(t)
t=linspace(-4*pi,4*pi,513);t=t[:-1]
y= pow(sin(t),3)
Y=fftshift(fft(y))/512.0
w=linspace(-64,64,513);w=w[:-1]

figure()

subplot(2,1,1)
plot(w,abs(Y),lw=2)
xlim([-5,5])
ylabel(r"$|Y|$",size=16)
title(r"Spectrum of $\sin^{3}\left(t\right)$")
grid(True)

subplot(2,1,2)

ii=where(abs(Y)>1e-3)
plot(w[ii],angle(Y[ii]),'go',lw=2)
xlim([-5,5])
ylabel(r"Phase of $Y$",size=16)
xlabel(r"$\omega$",size=16)
grid(True)
show()

#*******************************************************************************
# Finding FFT of the function cos^3(t)
t=linspace(-4*pi,4*pi,513);t=t[:-1]
y= pow(cos(t),3)
Y=fftshift(fft(y))/512.0
w=linspace(-64,64,513);w=w[:-1]

figure()

subplot(2,1,1)
plot(w,abs(Y),lw=2)
xlim([-5,5])
ylabel(r"$|Y|$",size=16)
title(r"Spectrum of $\cos^{3}\left(t\right)$")
grid(True)

subplot(2,1,2)

ii=where(abs(Y)>1e-3)
plot(w[ii],angle(Y[ii]),'go',lw=2)
xlim([-5,5])
ylabel(r"Phase of $Y$",size=16)
xlabel(r"$\omega$",size=16)
grid(True)
show()

#*******************************************************************************
# Finding FFT of the function cos(20t+5cos(t))
t=linspace(-4*pi,4*pi,513);t=t[:-1]
y= cos(20*t+5*cos(t))
Y=fftshift(fft(y))/512.0
w=linspace(-64,64,513);w=w[:-1]

figure()

subplot(2,1,1)
plot(w,abs(Y),lw=2)
xlim([-40,40])
ylabel(r"$|Y|$",size=16)
title(r"Spectrum of $\cos(20t+5\cos(t))\left(t\right)$")
grid(True)

subplot(2,1,2)

ii=where(abs(Y)>1e-3)
plot(w[ii],angle(Y[ii]),'go',lw=2)
xlim([-40,40])
ylabel(r"Phase of $Y$",size=16)
xlabel(r"$\omega$",size=16)
grid(True)
show()

#*******************************************************************************

# initial window_size and sampling rate defined
window_size = 2*pi
sampling_rate = 128
# tolerance for error
tolerance = 1e-15
# normalisation factor derived
norm_factor = (window_size)/(2*pi*(sampling_rate))

# In order to decrease the error in IFFT we increase the Window-Size iteratively.
# Also to overcome the issue of aliasing it is necessary to increase the sampling
# rate when window_size is increased.
for i in range(1, 10):

    t = linspace(-window_size/2, window_size/2,sampling_rate+1)[:-1]
    y = exp(-pow(t, 2)/2)
    N = sampling_rate
    Y = fftshift((fft(ifftshift(y)))*norm_factor)
    w_lim = (2*pi*N/((window_size)))
    w = linspace(-(w_lim/2), (w_lim/2), (sampling_rate+1))[:-1]

    # actual Y
    actual_Y = (1/sqrt(2*pi))*exp(-pow(w, 2)/2)
    error = (np.mean(np.abs(np.abs(Y)-actual_Y))) #Calculating error
    print("Absolute error at Iteration - %g is : %g" % ((i, error)))

    if(error < tolerance):
        print("\nAccuracy of the DFT is: %g and Iterations took: %g" %((error, i)))
        print("Best Window_size: %g , Sampling_rate: %g" %((window_size, sampling_rate)))
        break
    else:
        window_size = window_size*2
        sampling_rate = (sampling_rate)*2
        norm_factor = (window_size)/(2*pi*(sampling_rate))

figure()

subplot(2,1,1)
plot(w,abs(Y),lw=2)
xlim([-40,40])
ylabel(r"$|Y|$",size=16)
title(r"Spectrum of $exp(-x^{2}/2\left(t\right)$")
grid(True)

subplot(2,1,2)

ii=where(abs(Y)>1e-3)
plot(w[ii],angle(Y[ii]),'go',lw=2)
xlim([-40,40])
ylabel(r"Phase of $Y$",size=16)
xlabel(r"$\omega$",size=16)
grid(True)
show()

plot(w, abs(actual_Y),label=r"$\frac{1}{\sqrt{2}\pi} e^{\frac{\ - \omega ^{2}}{2}}$")
title("Exact Fourier Transform of Gaussian")
xlim([-10, 10])
ylabel(r"$Y(\omega) \to$")
xlabel(r"$\omega \to$")
grid()
legend()
show()
