import numpy as np
from scipy.optimize import root_scalar,root
from scipy.interpolate import UnivariateSpline
import scipy.constants as const
import re

e = const.e
pi = const.pi
epsilon0 = const.epsilon_0
mCa = 40.078 * const.u

def unit_to_num(unit):
    '''
    unit_to_num() transfers string that indicates unit to order of magnitude.
    '''
    if unit == '[mm]' or unit == 'mm':
        return 1e-3
    elif unit == '[um]' or unit == 'um':
        return 1e-6
    elif unit == '[nm]' or unit == 'nm':
        return 1e-9
    else:
        raise TypeError('No type matching!')

def eqposition_sym(z:np.ndarray, esbeV:np.ndarray, unit='[um]', k=5, smooth=1e-6, pos0=1e-3, step=0.02):
    '''
    eqposition_sym() will find the equilibrium postion and eigenmodes based on the assumption that the electric potential function is symmetric.
    
    Input
    --------
    z: np.ndarray
        Coordinate.
    - esbeV: Electric potential.
    - unit: Unit of z. Default value is '[um]'.
    - k: Degree of the smoothing spline. Must be `1 <= k <= 5`.
    - s: Smoothing factor of spline.
    - pos0: The starting point for the root search (in the unit of rescaled coordinate).
    - step: The step length for the root search (in the unit of rescaled coordinate).

    Output
    --------
    - z:

    It will iterate through positive $z$ coordinates, and compute the derivative of potential function $V$ at $[-z,z]$ until it finds root.
    eqposition_sym() only supports 2 ions. 
    '''

    # Rescale the coordinate to [-1,1]
    zlim = np.max(z)
    zeta = z/zlim
    r0 = zlim * unit_to_num(unit)

    # Interpolation
    Vspl = UnivariateSpline(zeta,esbeV,k=k,s=smooth) # spl.get_residual() <= s
    dVspl = Vspl.derivative(1)

    # Define the function to be solved
    def dfunc(zeta,N=2):
        y=np.empty(N)
        for n in range(N):
            y[n] = dVspl(zeta[n]) + e/(4*pi*epsilon0*r0)*(-np.sum(1/(zeta[n]-zeta[:n])**2)+np.sum(1/(zeta[n]-zeta[n+1:])**2))
        return y
    def solfunc(pos):
        zeta = np.array([-pos,pos])
        return (dfunc(zeta)[0]-dfunc(zeta)[1])/2

    # Find the root interval
    while solfunc(pos0)*solfunc(pos0+step)>0:
        pos0 = pos0 + step
    
    # Solve for the equilibrium position
    eqpos = root_scalar(solfunc,method='bisect',bracket=[pos0,pos0+step]).root
    zeta0 = np.array([-eqpos,eqpos])
    z0 = zeta0 * zlim
    return z0

def eqposition(z:np.ndarray, esbeV:np.ndarray, inipos, N=2, unit='[um]', k=5, smooth=1e-6, method='hybr'):
    '''
    eqposition will find the equilibrium position of ions.

    Input:
    - z: Coordinate. The unit of inipos should be the parameter `unit`.
    - esbeV: Electric potential.
    - inipos: Initial guess for ion positions that will be passed into `scipy.optimize.root`. The unit of inipos should be the parameter `unit`.
    - N: ion number. Default is 2 ions.
    - unit: Unit of z.
    - k: Degree of the smoothing spline. Must be 1 <= k <= 5.
    - s: Smoothing factor.
    - m: Mass of ion. Default value is Ca+.
    - method: methods provided by scipy.optimize.root
    '''

    # Rescale the coordinate to [-1,1]
    zlim = max(np.abs(np.max(z)),np.abs(np.min(z)))
    zeta = z/zlim
    inipos_zeta = inipos/zlim
    r0 = zlim * unit_to_num(unit)

    # Interpolation
    Vspl = UnivariateSpline(zeta,esbeV,k=k,s=smooth) # spl.get_residual() <= s
    dVspl = Vspl.derivative(1)
    ddVspl = Vspl.derivative(2)

    def solfun(zeta): # the equations
        y=np.empty(shape=(N))
        for n in range(N):
            y[n] = dVspl(zeta[n]) + e/(4*pi*epsilon0*r0)*(-np.sum(1/(zeta[n]-zeta[:n])**2)+np.sum(1/(zeta[n]-zeta[n+1:])**2))
        return y
    
    def jac(zeta): # Jacobi matrix for scipy.optimize.root
        J = np.empty((N,N))
        for i in range(N):
            for j in range(N):
                if j==i:
                    J[i,j]= ddVspl(zeta[i]) + e/(4*pi*epsilon0*r0)*2*(np.sum(1/(zeta[i]-zeta[:i])**3)-np.sum(1/(zeta[i]-zeta[i+1:])**3))  # m == n
                else:
                    J[i,j]= -2*e/(4*pi*epsilon0*r0)/np.absolute(zeta[i]-zeta[j])**3 # m!=n
        return J
    
    zeta0 = root(solfun,inipos_zeta,method=method,jac=jac).x
    z0 = zeta0 * zlim
    return z0

def mode1d(z:np.ndarray, esbeV:np.ndarray, z0:np.ndarray, N=2, unit='[um]', k=5, smooth=1e-6, m=mCa):
    '''
    mode1d() will solve for frequency f, NOT angular frequency.

    Input:
    - z: Coordinate.
    - esbeV: Electric potential.
    - z0: Equilibrium postion of ions.
    - N: Number of ions. Default value is 2.
    '''

    # Rescale the coordinate to [-1,1]
    zlim = max(np.abs(np.max(z)),np.abs(np.min(z)))
    zeta = z/zlim
    zeta0 = z0/zlim
    r0 = zlim * unit_to_num(unit)

    # Interpolation
    Vspl = UnivariateSpline(zeta,esbeV,k=k,s=smooth) # spl.get_residual() <= s
    ddVspl = Vspl.derivative(2)

    # Solve for the eigenmodes
    ddU = np.empty((N,N))
    for i in range(N):
        for j in range(N):
            if j==i:
                ddU[i,j]= 4*pi*epsilon0*r0/e * ddVspl(zeta0[i]) +\
                    2*(np.sum(1/(zeta0[i]-zeta0[:i])**3) - np.sum(1/(zeta0[i]-zeta0[i+1:])**3))  # m == n
            else:
                ddU[i,j]=-2/np.absolute(zeta0[i]-zeta0[j])**3 # m!=n
    w,v = np.linalg.eig(ddU)
    idx = np.argsort(w)
    w = w[idx]
    # v = v[:,idx]
    
    return np.sqrt(w/m*1/r0**2 * e**2/(4*pi*epsilon0*r0))/(2*pi)

def eqposition_analytic(Vfunc, dVfunc, ddVfunc, inipos, N=2, unit='[um]', method='hybr'):
    '''
    eqposition_analytic will find the equilibrium position of ions.

    Input:
    - z: Coordinate.
    - esbeV: Electric potential.
    - inipos: Initial guess for ion positions that will be passed into `scipy.optimize.root`
    - N: ion number. Default is 2 ions.
    - unit: Unit of z and inipos. Default is micrometer.
    - k: Degree of the smoothing spline. Must be 1 <= k <= 5.
    - s: Smoothing factor.
    - m: Mass of ion. Default value is Ca+.
    - method: methods provided by scipy.optimize.root
    '''

    # Rescale the coordinate to [-1,1]
    r0 = unit_to_num(unit)

    def solfun(zeta): # the equations
        y=np.empty(shape=(N))
        for n in range(N):
            y[n] = dVfunc(zeta[n]) + e/(4*pi*epsilon0*r0)*(-np.sum(1/(zeta[n]-zeta[:n])**2)+np.sum(1/(zeta[n]-zeta[n+1:])**2))
        return y
    
    def jac(zeta): # Jacobi matrix for scipy.optimize.root
        J = np.empty((N,N))
        for i in range(N):
            for j in range(N):
                if j==i:
                    J[i,j]= ddVfunc(zeta[i]) + e/(4*pi*epsilon0*r0)*2*(np.sum(1/(zeta[i]-zeta[:i])**3)-np.sum(1/(zeta[i]-zeta[i+1:])**3))  # i==j
                else:
                    J[i,j]= -2*e/(4*pi*epsilon0*r0)/np.absolute(zeta[i]-zeta[j])**3 # i!=j
        return J
    
    z0 = root(solfun,inipos,method=method,jac=jac).x
    return z0

def mode1d_analytic(ddVfunc, z0, N=2, unit='[um]', method='hybr', m=mCa):
    '''
    mode1d() will solve for frequency f, NOT angular frequency.

    Input:
    - z: Coordinate.
    - esbeV: Electric potential.
    - z0: Equilibrium postion of ions.
    - N: Number of ions. Default value is 2.
    '''

    # Rescale the coordinate to [-1,1]
    r0 = unit_to_num(unit)
    # Solve for the eigenmodes
    ddU = np.empty((N,N))
    for i in range(N):
        for j in range(N):
            if j==i:
                ddU[i,j]= 4*pi*epsilon0*r0/e * ddVfunc(z0[i]) +\
                    2*(np.sum(1/(z0[i]-z0[:i])**3) - np.sum(1/(z0[i]-z0[i+1:])**3))  # m == n
            else:
                ddU[i,j]=-2/np.absolute(z0[i]-z0[j])**3 # m!=n
    w,v = np.linalg.eig(ddU)
    idx = np.argsort(w)
    w = w[idx]
    # v = v[:,idx]
    
    return np.sqrt(w/m*1/r0**2 * e**2/(4*pi*epsilon0*r0))/(2*pi)
