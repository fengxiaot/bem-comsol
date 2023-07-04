import numpy as np
from scipy.optimize import newton
from globalvar import *

def Eqposition(N,dVfun,ddVfun,inipos):
    def func(zeta): # the equations
        y=np.empty(shape=(N))
        for n in range(N):
            y[n] = dVfun(zeta[n]) + e/(4*pi*epsilon0*r0)*(-np.sum(1/(zeta[n]-zeta[:n])**2)+np.sum(1/(zeta[n]-zeta[n+1:])**2))
        return y
    
    def dfunc(zeta): # first order derivative
        y=np.empty(shape=(N))
        for n in range(N):
            y[n]= ddVfun(zeta[n]) + e/(4*pi*epsilon0*r0)*2*(np.sum(1/(zeta[n]-zeta[:n])**3)-np.sum(1/(zeta[n]-zeta[n+1:])**3))
        return y
    zeta0 = newton(func,inipos,fprime=dfunc,maxiter=100000) # Newton method
    return(zeta0)

def AxialFull(N,ddVfun,zeta0):
    ddU = np.empty((N,N))
    for n in range(N):
        for m in range(N):
            if m==n:
                ddU[n,m]= 4*pi*epsilon0*r0/e * ddVfun(zeta0[n]) +\
                    2*(np.sum(1/(zeta0[n]-zeta0[:n])**3) - np.sum(1/(zeta0[n]-zeta0[n+1:])**3))  # m == n
            else:
                ddU[n,m]=-2/np.absolute(zeta0[n]-zeta0[m])**3 # m!=n
    w,v = np.linalg.eig(ddU)
    idx = np.argsort(w)
    w = w[idx]
    v = v[:,idx]
    return ([w,v])

def AxialMode(N,ddVfun,zeta0):
    [w,v] = AxialFull(N,ddVfun,zeta0)
    return np.sqrt(w/m*1/r0**2 * e**2/(4*pi*epsilon0*r0))
