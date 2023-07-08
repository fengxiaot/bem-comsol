import numpy as np
from scipy.optimize import root_scalar,root
from globalvar import *

def SymEqposition(dVfun,N=2):
    '''
    SymEqpostion will find the equilibrium postion and eigenmodes based on the assumption that the electric potential is symetric.
    It will iterate through positive $z$ coordinate, and compute the derivative of potential function $V$ at $[-z,z]$ until it finds root.
    SymEqpostion only supports 2 ions. 
    '''
    def dfunc(N,zeta):
        y=np.empty(shape=(N))
        for n in range(N):
            y[n] = dVfun(zeta[n]) + e/(4*pi*epsilon0*r0)*(-np.sum(1/(zeta[n]-zeta[:n])**2)+np.sum(1/(zeta[n]-zeta[n+1:])**2))
        return y
    def solfunc(pos):
        zeta = np.array([-pos,pos])
        return (dfunc(2,zeta)[0]-dfunc(2,zeta)[1])/2
    # ==================================================
    # This code can find the root interval
    pos0 = 1e-3
    step = 0.02
    while solfunc(pos0)*solfunc(pos0+step)>0:
        pos0 = pos0 + step
    # ==================================================
    eqpos = root_scalar(solfunc,method='bisect',bracket=[pos0,pos0+step]).root
    zeta0 = np.array([-eqpos,eqpos])
    return zeta0

def Eqposition(dVfun,ddVfun,inipos,N=2,method='hybr'):
    '''
    Eqposition will find the equilibrium position of ions. Note: All functions should be vectorized pyfuncs.
    - dVfun: the first order derivative of 1D electric potential function
    - ddVfun: the second order derivative of 1D electric potential function
    - inipos: initial guess position that will be passed into scipy.optimize.root
    - N: ion number
    - method: methods provided by scipy.optimize.root
    '''
    def func(zeta): # the equations
        y=np.empty(shape=(N))
        for n in range(N):
            y[n] = dVfun(zeta[n]) + e/(4*pi*epsilon0*r0)*(-np.sum(1/(zeta[n]-zeta[:n])**2)+np.sum(1/(zeta[n]-zeta[n+1:])**2))
        return y
    
    def jac(zeta): # Jacobi matrix for scipy.optimize.root
        J = np.empty((N,N))
        for n in range(N):
            for m in range(N):
                if m==n:
                    J[n,m]= ddVfun(zeta[n]) + e/(4*pi*epsilon0*r0)*2*(np.sum(1/(zeta[n]-zeta[:n])**3)-np.sum(1/(zeta[n]-zeta[n+1:])**3))  # m == n
                else:
                    J[n,m]= -2*e/(4*pi*epsilon0*r0)/np.absolute(zeta[n]-zeta[m])**3 # m!=n
        return J
    zeta0 = root(func,inipos,method=method,jac=jac).x
    return zeta0

def AxialFull(ddVfun,zeta0,N=2):
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

def AxialMode(ddVfun,zeta0,N=2):
    '''
    AxialMode() will output the frequency f, NOT angular frequency
    '''
    [w,v] = AxialFull(ddVfun,zeta0,N)
    return np.sqrt(w/m*1/r0**2 * e**2/(4*pi*epsilon0*r0))/(2*pi)
