import scipy as sp
import numpy as np
import matplotlib
import matplotlib.pyplot as py



#Some global settings for matplotlib plots
matplotlib.rcParams['font.size'] = 12
matplotlib.rcParams['font.weight'] = 'bold'



#Equation whose root defines the homogeneous equilibrium solutions:
def func(beta,mu):
    return lambda rho : rho-(1.0-rho)*np.exp(beta*(mu+5.0*rho))



# Helpful wrapper for the root finder in scipy which handles its annoying errors:
def find_roots(f):
  roots = []
  #Compute a root over the entire range [0,1]
  sol = sp.optimize.brentq(f,0,1) #Use Brent's method - will find smallest root
  roots.append(sol)

  #Find if there are any other roots by using intervals above this solution
  tol = 1e-2 #Small shift above found root
  start = sol + tol
  step = 1e-1 #Small step used to find bounding interval for next root
  finish = sol + step
  #Continue looping and adding roots until finish becomes unphysical (>1)

  while finish <= 1.0:
    #Brent method will fail if f(start) and f(finish) have the same sign,
    #so keep trying until it works:

    try:
      sol = sp.optimize.brentq(f,start,finish)
      roots.append(sol)
      start = sol + tol #Pick new interval start beyond root just found
      finish = sol + step #Assign new finish beyond it also and continue looping

    except:
      finish += step  #Enlarge interval until root can be found (if any)
  return roots




#Main script:


N = 100 #Density samples to use
rho = np.linspace(0,1,N) #Density values to evaluate



#Evaluate the density profiles for the various parameter regimes
muVals = [-3.0]   #Chemical potentials to test



# --- Low temperature case ---
beta = 2/3
funcs = [func(beta,mu) for mu in muVals] #Define a list of functions with the beta and mu specified



#Evaluate the functions for plotting
lowT = [f(rho) for f in funcs] #Returns a list of Nx1 float arrays

roots = [find_roots(f) for f in funcs] #Find the roots for each function




# #Plot
# fig1, ax1 = py.subplots()
# for curve in lowT:
#   ax1.plot(rho,curve)

# ax1.plot(roots, np.zeros(len(roots)), 'o', color = "red")

# ax1.set_xlabel('Density Rho')
# ax1.set_ylabel('Function')

kB = 1E-23
a = 1E-23


def rho_plus_minus(sign, yarr):
    # sign is +1 or -1
    # return (1/2)*(1 + sign * np.sqrt(1 - 4*kB*yarr/(5*a)))
    return (1/2)*(1 + sign * np.sqrt(1 - 4*kB*yarr/(5*a)))

def chem_pot_func():
    fig2, ax2 = py.subplots()
    ax2.set_ylim(0, 1.6)

    yarr = np.linspace(0, 1.25, N)

    rho1 = rho_plus_minus(1, yarr)
    rho2 = rho_plus_minus(-1, yarr)


    xarr1 = kB*yarr*np.log(rho1/(rho2)) - 5*a*rho1
    xarr2 = kB*yarr*np.log(rho2/(rho1)) - 5*a*rho2
    xarr3 = np.ones(100)*(-5/2 * a)


    ax2.plot(xarr1,yarr, "blue")
    ax2.plot(xarr2,yarr, "blue")
    ax2.plot(xarr3,yarr, "red")

    ax2.plot(np.array([-2,-2.5,-3,-2,-2.5,-3])*1E-23,[1,1,1,3/2,3/2,3/2], 'o', color = "orange")

betaarr = np.linspace(0,2,N)
muarr = np.linspace(-5,0,N)

funcs = []

for beta in betaarr:
    funcs += [func(beta, mu) for mu in muarr]

roots = [find_roots(f) for f in funcs] #Find the roots for each function

lowdensities = [min(i) if min(i) < 0.5 else 0 for i in roots]
highdensities = [max(i) if max(i) > 0.5 else 0 for i in roots]

def density_colormap():
    fig2, ax2 = py.subplots()
    ax2.set_ylim(0, 1.6)

    yarr = np.linspace(0, 1.25, N)

    rho1 = rho_plus_minus(1, yarr)
    rho2 = rho_plus_minus(-1, yarr)


    xarr1 = kB*yarr*np.log(rho1/(rho2)) - 5*a*rho1
    xarr2 = kB*yarr*np.log(rho2/(rho1)) - 5*a*rho2
    xarr3 = np.ones(100)*(-5/2 * a)


    ax2.plot(xarr1,yarr, "blue")
    ax2.plot(xarr2,yarr, "blue")
    ax2.plot(xarr3,yarr, "red")

    ax2.plot(xarr1, yarr, cmap=lowdensities)
    py.show()

density_colormap()
