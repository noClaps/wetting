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

# for i in range[curve]:
#     if 


#Plot
fig1, ax1 = py.subplots()
for curve in lowT:
  ax1.plot(rho,curve)
  
ax1.plot(roots, np.zeros(len(roots)), 'o', color = "red")

ax1.set_xlabel('Density Rho')
ax1.set_ylabel('Function')

