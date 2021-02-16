from pylab import *
import mpl_toolkits.mplot3d.axes3d as p3
import os,sys
import scipy
import scipy.linalg as s
import matplotlib.pyplot as plt

# Obtain Parameters via commandline arguments
if(len(sys.argv)==5):
    Nx = int(sys.argv[1])# size along x
    Ny = int(sys.argv[2])# size along y
    radius = int(sys.argv[3]) # radius of central lead
    Niter = int(sys.argv[4]) # number of iterations to perform
else:
    Nx = 25; # size along x
    Ny = 25; # size along y
    radius = 8; # radius of central lead
    Niter = 1500; # number of iterations to perform

# Creating Meshgrid and Initialising Potential phi
x = np.linspace(-0.5,0.5,num=Nx,dtype=float)
y = np.linspace(-0.5,0.5,num=Ny,dtype=float)
phi = np.zeros((Nx,Ny),dtype = float)
Y,X = meshgrid(y,x)
phi[np.where(X**2+Y**2<(0.35)**2)] = 1.0

#plot potential
xlabel("X")
ylabel("Y")
contourf(X,Y,phi)
colorbar()
show()

#function to update potential
def phi_new(phi,phiold):
    phi[1:-1,1:-1]=0.25*(phiold[1:-1,0:-2]+ phiold[1:-1,2:]+ phiold[0:-2,1:-1] + phiold[2:,1:-1])
    return phi

# function to update potential at boundaries
def boundary(phi,central_portion = np.where(X**2+Y**2<(0.35)**2)):
    phi[:,0] = phi[:,1] # Left Boundary
    phi[:,Nx-1] = phi[:,Nx-2] # Right Boundary
    phi[0,:] = phi[1,:] # Top Boundary
    phi[Ny-1,:] = 0 # Bottom boundary connected to ground
    phi[central_portion] = 1.0 #recreating boundary values at electrode
    return phi

err = np.zeros(Niter)# initialise error array

for k in range(Niter):
    phiold = phi.copy()# copy the old phi
    phi = phi_new(phi,phiold)# Updating the potential
    phi = boundary(phi)# applying boundary conditions
    err[k] = np.max(np.abs(phi-phiold))# Appending errors for each iterations
    if(err[k] == 0):
        print("Reached steady state at ",k," Iterations")
        break

#plotting Error on semilog
title("Error - semilog plot")
xlabel("No of iterations")
ylabel("Error")
semilogy(range(Niter),err)
show()
#plotting Error on loglog
title("Error - loglog plot")
xlabel("No of iterations")
ylabel("Error")
loglog((np.asarray(range(Niter))+1),err)
loglog((np.asarray(range(Niter))+1)[::50],err[::50],'ro')
legend(["real","every 50th value"])
show()

#function for getting best fit
def best_fit(y,Niter,lastn=0):
    log_err = np.log(err)[-lastn:]
    X = np.vstack([(np.arange(Niter)+1)[-lastn:],np.ones(log_err.shape)]).T
    log_err = np.reshape(log_err,(1,log_err.shape[0])).T
    return s.lstsq(X, log_err)[0]

def plot_error(err,Niter,a,a_,b,b_):
    title("Best fit for error on a loglog scale")
    xlabel("No of iterations")
    ylabel("Error")
    x = np.asarray(range(Niter))+1
    loglog(x,err)
    loglog(x[::100],np.exp(a+b*np.asarray(range(Niter)))[::100],'ro')
    loglog(x[::100],np.exp(a_+b_*np.asarray(range(Niter)))[::100],'go')
    legend(["errors","fit1","fit2"])
    show()

    title("Best fit for error on a semilog scale")
    xlabel("No of iterations")
    ylabel("Error")
    semilogy(x,err)
    semilogy(x[::100],np.exp(a+b*np.asarray(range(Niter)))[::100],'ro')
    semilogy(x[::100],np.exp(a_+b_*np.asarray(range(Niter)))[::100],'go')
    legend(["errors","fit1","fit2"])
    show()

def find_net_error(a,b,Niter):
    return -a/b*np.exp(b*(Niter+0.5))

b,a = best_fit(err,Niter)
b_,a_ = best_fit(err,Niter,500)
plot_error(err,Niter,a,a_,b,b_)
#plotting cumulative error
iter=np.arange(100,1501,100)
grid(True)
title(r'Plot of Cumulative Error values On a loglog scale')
loglog(iter,np.abs(find_net_error(a_,b_,iter)),'ro')
xlabel("iterations")
ylabel("Net  maximum error")
show()

#plotting 3d contour of final potential
fig1=plt.figure(4)     # open a new figure
ax=p3.Axes3D(fig1) # Axes3D is the means to do a surface plot
title('The 3-D surface plot of the potential')
surf = ax.plot_surface(Y, X, phi.T, rstride=1, cstride=1, cmap=plt.cm.jet)
fig1.colorbar(surf, shrink=0.5, aspect=5)
show()

#plotting 2d contour of final potential
title("2D Contour plot of potential")
xlabel("X")
ylabel("Y")
x_c,y_c=np.where(X**2+Y**2<(0.35)**2)
plot((x_c-Nx/2)/Nx,(y_c-Ny/2)/Ny,'ro')
contourf(Y,X[::-1],phi)
colorbar()
show()

Jx = (1/2*(phi[1:-1,0:-2]-phi[1:-1,2:]))
Jy = (1/2*(phi[:-2,1:-1]-phi[2:,1:-1]))

title("Vector plot of current flow")
quiver(Y[1:-1,1:-1],-X[1:-1,1:-1],-Jx[:,::-1],-Jy)
x_c,y_c=np.where(X**2+Y**2<(0.35)**2)
plot((x_c-Nx/2)/Nx,(y_c-Ny/2)/Ny,'ro')
show()
