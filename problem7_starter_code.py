import scipy as sp
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

# Some global settings for matplotlib plots
matplotlib.rcParams['font.size'] = 12
matplotlib.rcParams['font.weight'] = 'bold'

# Helper functions
def func(beta,mu):
    return lambda rho : rho-(1.0-rho)*np.exp(beta*(mu+5.0*rho))

def find_roots(f):
  roots = []
  # Compute a root over the entire range [0,1]
  sol = sp.optimize.brentq(f,0,1) #Use Brent's method - will find smallest root
  roots.append(sol)

  # Find if there are any other roots by using intervals above this solution
  tol = 1e-2 # Small shift above found root
  start = sol + tol
  step = 1e-2 # Small step used to find bounding interval for next root
  finish = sol + step
  # Continue looping and adding roots until finish becomes unphysical (>1)

  while finish <= 1.0:
    # Brent method will fail if f(start) and f(finish) have the same sign,
    # so keep trying until it works:

    try:
      sol = sp.optimize.brentq(f,start,finish)
      roots.append(sol)
      start = sol + tol # Pick new interval start beyond root just found
      finish = sol + step # Assign new finish beyond it also and continue looping

    except:
      finish += step  # Enlarge interval until root can be found (if any)
  return roots

def rho_plus_minus(sign, kBTe):
    # sign is +1 or -1
    return (1/2)*(1 + sign * np.sqrt(1 - 4*kBTe/5))


# Plotting functions
def plotroots(num_samples = 100):
    rho = np.linspace(0,1,num_samples) # Density values to evaluate

    # Evaluate the density profiles for the various parameter regimes
    muVals = [-3.0]   # Chemical potentials to test

    # --- Low temperature case ---
    beta = 2/3
    funcs = [func(beta,mu) for mu in muVals] # Define a list of functions with the beta and mu specified

    # Evaluate the functions for plotting
    lowT = [f(rho) for f in funcs] # Returns a list of Nx1 float arrays
    roots = [find_roots(f) for f in funcs] # Find the roots for each function

    fig1, ax1 = plt.subplots()
    for curve in lowT:
        ax1.plot(rho,curve)

    ax1.plot(roots, np.zeros(len(roots)), 'o', color = "red")

    ax1.set_xlabel('Density Rho')
    ax1.set_ylabel('Function')

def chem_pot_func(num_samples = 100):
    fig2, ax2 = plt.subplots()
    ax2.set_ylim(0, 1.6)

    kBTe = np.linspace(0, 1.25, num_samples)

    rho1 = rho_plus_minus(1, kBTe)
    rho2 = rho_plus_minus(-1, kBTe)


    leftline = kBTe*np.log(rho1/rho2) - 5*rho1
    rightline = kBTe*np.log(rho2/rho1) - 5*rho2
    centreline = np.ones(100)*(-5/2)

    ax2.set_xlabel(r"$\mu / \epsilon$")
    ax2.set_ylabel(r"$k_B T / \epsilon$")

    ax2.plot(leftline,kBTe, "blue")
    ax2.plot(rightline,kBTe, "blue")
    ax2.plot(centreline,kBTe, "red")

    ax2.plot(np.array([-2,-2.5,-3,-2,-2.5,-3]),[1,1,1,3/2,3/2,3/2], 'o', color = "orange")


def phase_density(num_samples = 100):
    kBTarr = np.linspace(0.01,2,num_samples)
    muarr = np.linspace(-5,0,num_samples)
    densities = np.ones((num_samples,num_samples))

    for i in range(num_samples):
        for j in range(num_samples):
            f = func(1/kBTarr[i], muarr[j])
            roots = find_roots(f)

            if len(roots) == 1:
                densities[i,j] = roots[0]
            else:
                densities[i,j] = -1

    kBTe = np.linspace(0.01, 1.25, num_samples)

    rho1 = rho_plus_minus(1, kBTe)
    rho2 = rho_plus_minus(-1, kBTe)

    leftline = kBTe*np.log(rho1/rho2) - 5*rho1
    rightline = kBTe*np.log(rho2/rho1) - 5*rho2
    centreline = np.ones(num_samples)*(-5/2)
    criticaltempline = np.ones(num_samples)*1.25

    fig3, ax3 = plt.subplots(figsize=(8,6))

    plt.pcolormesh(muarr, kBTarr, densities)

    plt.plot(leftline, kBTe, "blue", linewidth=4)
    plt.plot(rightline, kBTe, "blue", linewidth=4)
    plt.plot(centreline, kBTe, "red", linewidth=3)
    plt.plot(np.linspace(-5, 0, num_samples), criticaltempline, color="black", linestyle="--", linewidth=2)

    plt.xlabel(r"$\mu / \epsilon$")
    plt.ylabel(r"$k_BT/\epsilon$")
    plt.title("Phase Diagram with Density Gradient")
    cbar = plt.colorbar()
    cbar.set_label(r"Density $\rho$", rotation=270, labelpad=20)

def main():
    # plotroots()
    # chem_pot_func()
    phase_density()
    plt.show()

if __name__ == "__main__":
    main()
