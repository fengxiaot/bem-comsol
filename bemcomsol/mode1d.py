import numpy as np
from scipy.optimize import root_scalar,root
from scipy.interpolate import UnivariateSpline
import scipy.constants as const

e = const.e
pi = const.pi
epsilon0 = const.epsilon_0
mCa = 40.078 * const.u
m_e = const.m_e

def unit_to_num(unit):
    '''
    unit_to_num() transfers string that indicates unit to order of magnitude.
    '''
    if unit == '[m]' or unit == 'm':
        return 1
    elif unit == '[mm]' or unit == 'mm':
        return 1e-3
    elif unit == '[um]' or unit == 'um':
        return 1e-6
    elif unit == '[nm]' or unit == 'nm':
        return 1e-9
    else:
        raise TypeError('No type matching!')

def eqposition_sym(z:np.ndarray, esbeV:np.ndarray, unit='[um]', k=5, s=1e-6, pos0=1e-3, step=0.02):
    '''
    eqposition_sym() will find the equilibrium postion and eigenmodes
    based on the assumption that the electric potential function is symmetric.
    
    NOTE: eqposition_sym() only support 2 ions.

    Parameters
    ----------
    z : np.ndarray
        Coordinate.
    esbeV : np.ndarray
        Corresponding electric potential V(z).
    unit : str
        Unit of z. Default value is '[um]'.
    k : int
        Degree of the smoothing spline. Must be `1 <= k <= 5`.
    s : float
        Smoothing factor of spline.
    pos0 : float
        The starting point for the root search (in the unit of rescaled coordinate).
    step : float
        The step length for the root search (in the unit of rescaled coordinate).

    Returns
    -------
    z0 : np.ndarray
        Equilibrium positions.
    '''

    # Rescale the coordinate to [-1,1]
    zlim = np.max(z)
    zeta = z/zlim
    r0 = zlim * unit_to_num(unit)

    # Interpolation
    Vspl = UnivariateSpline(zeta,esbeV,k=k,s=s) # spl.get_residual() <= s
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

def eqposition(z:np.ndarray, esbeV:np.ndarray, inipos:np.ndarray, N=2, unit='[um]', k=5, s=1e-6, method='hybr'):
    '''
    Find the equilibrium position of ions.

    Parameters
    ----------
    z : np.ndarray
        Coordinate.
    esbeV : np.ndarray
        Corresponding electric potential V(z).
    inipos : np.ndarray
        Initial guess of equilibrium position of ions.
    N : int
        Number of ions. Default value is 2.
    unit : str
        Unit of z. Default value is '[um]'.
    k : int
        Degree of the smoothing spline. Must be `1 <= k <= 5`.
    s : float
        Smoothing factor of spline.
    method : str
        Type of solver provided by `scipy.optimize.root`

    Returns
    -------
    z0 : np.ndarray
        Equilibrium positions.
    '''

    # Rescale the coordinate to [-1,1]
    zlim = max(np.abs(np.max(z)),np.abs(np.min(z)))
    zeta = z/zlim
    inipos_zeta = inipos/zlim
    r0 = zlim * unit_to_num(unit)

    # Interpolation
    Vspl = UnivariateSpline(zeta,esbeV,k=k,s=s) # spl.get_residual() <= s
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

def mode1d(z:np.ndarray, esbeV:np.ndarray, z0:np.ndarray, N=2, unit='[um]', k=5, s=1e-6, m=mCa):
    '''
    Solve for trap frequencies.

    NOTE: The result is not angular frequency but ordinary frequency.

    Parameters
    ----------
    z : np.ndarray
        Coordinate.
    esbeV : np.ndarray
        Corresponding electric potential V(z).
    z0 : np.ndarray
        Equilibrium position of ions.
    N : int
        Number of ions. Default value is 2.
    unit : str
        Unit of z. Default value is '[um]'.
    k : int
        Degree of the smoothing spline. Must be `1 <= k <= 5`.
    s : float
        Smoothing factor of spline.
    m : float
        Mass of ions.
    method : str
        Type of solver provided by `scipy.optimize.root`
      
    Returns
    -------
    w : np.ndarray
        Eigenmodes.
    '''

    # Rescale the coordinate to [-1,1]
    zlim = max(np.abs(np.max(z)),np.abs(np.min(z)))
    zeta = z/zlim
    zeta0 = z0/zlim
    r0 = zlim * unit_to_num(unit)

    # Interpolation
    Vspl = UnivariateSpline(zeta,esbeV,k=k,s=s) # spl.get_residual() <= s
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
    Find the equilibrium position of ions.

    Parameters
    ----------
    Vfunc : function
        Electric potential function.
    dVfunc : function
        First order derivative of electric potential function.
    ddVfunc : function
        Second order derivative of electric potential function.
    inipos : np.ndarray
        Initial guess of equilibrium position of ions.
    N : int
        Number of ions. Default value is 2.
    unit : str
        Unit of z. Default value is '[um]'.
    method : str
        Type of solver provided by `scipy.optimize.root`

    Returns
    -------
    z0 : np.ndarray
        Equilibrium positions.
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
    Solve for trap frequencies.

    NOTE: The result is not angular frequency but ordinary frequency.

    Parameters
    ----------
    ddVfunc : function
        Second order derivative of the electrostatic potential function.
    z0 : np.ndarray
        Equilibrium position of ions.
    N : int
        Number of ions. Default value is 2.
    unit : str
        Unit of z. Default value is '[um]'.

    m : float
        Mass of ions.
    
    Returns
    -------
    w : np.ndarray
        Eigenmodes.
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
