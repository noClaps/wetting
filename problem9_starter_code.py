import scipy.optimize as sp
import numpy as np
import matplotlib.pyplot as plt

#Equation whose root defines the homogeneous equilibrium solutions:
def func(beta,mu):
    return lambda rho : rho-(1.0-rho)*np.exp(beta*(mu+5.0*rho))

# Helpful wrapper for the root finder in scipy which handles its annoying errors:
def find_roots(f):
  roots = []
  #Compute a root over the entire range [0,1]
  sol = sp.brentq(f,0,1) #Use Brent's method - will find smallest root
  roots.append(sol)
  #Find if there are any other roots by using intervals above this solution
  tol = 1e-3 #Small shift above found root
  start = sol + tol
  step = 1e-2 #Small step used to find bounding interval for next root
  finish = sol + step
  #Continue looping and adding roots until finish becomes unphysical (>1)
  while finish <= 1.0:
    #Brent method will fail if f(start) and f(finish) have the same sign,
    #so keep trying until it works:
    try:
      sol = sp.brentq(f,start,finish)
      roots.append(sol)
      start = sol + tol #Pick new interval start beyond root just found
      finish = sol + step #Assign new finish beyond it also and continue looping
    except:
      finish += step  #Enlarge interval until root can be found (if any)
  return roots

#Main script:

Lx = 4 #Number of sites along x-axis
Ly = 4 #Number of sites along y-axis
beta = 1 #Inverse temperature beta*epsilon
mu = -2.5 #Chemical potential mu/epsilon

rho_0 = 0.51 #Initial density
tol = 1e-12 #Convergence tolerance
count = 30000 #Upper limit for iterations
alpha  = 0.01 #Mixing parameter

#Solve equations iteratively:
conv = 1
cnt = 1
rho = rho_0*np.ones([Lx,Ly])
rho_new = np.zeros([Lx,Ly]);
while conv >= tol and cnt<count:
  cnt = cnt + 1
  for i in range(Lx):
    for j in range(Ly):
      #Handle the periodic boundaries for x and y:
      left = np.mod((i-1),Lx) #i-1, maps -1 to Lx-1
      right = np.mod((i+1),Lx) #i+1, maps Lx to 0
      down = np.mod((j-1),Ly) #j-1, maps -1 to Ly-1
      up = np.mod((j+1),Ly) #j+1, maps Ly to 0
      rho_new[i,j] = (1 - rho[i,j])*np.exp(beta*(rho[i,down] + rho[i,up] + rho[left,j] + rho[right,j] + (1/4)*(rho[left,down] + rho[right,down] + rho[left,up] + rho[right,up]) + mu))

  conv = sum(sum((rho - rho_new)**2)); #Compute the convergence parameter.
  rho = alpha*rho_new + (1 - alpha)*rho #Mix the new and old solutions.

print(rho[0,0])
sol = find_roots(func(beta,mu))
print(sol)
