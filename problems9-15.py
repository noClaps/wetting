import scipy.optimize as sp
import numpy as np
import matplotlib.pyplot as plt
import math
from matplotlib import cm
from matplotlib.colors import ListedColormap


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


def problem9():
    Lx = 4 #Number of sites along x-axis
    Ly = 4 #Number of sites along y-axis
    beta = 1 #Inverse temperature beta*epsilon
    mu = -2.5 #Chemical potential mu/epsilon

    rho_0 = 0.49 #Initial density
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

    # plt.imshow(rho, extent=(0, Lx, 0, Ly), vmin=-1, vmax=1)
    # cbar = plt.colorbar()
    # cbar.set_label(r"Density $\rho$", rotation=270, labelpad=20)

    # plt.xlabel("Lattice points")
    # plt.ylabel("Lattice points")

    # plt.grid(True, linestyle="--")
    # plt.xticks(np.linspace(0, Lx, Lx+1))
    # plt.yticks(np.linspace(0, Ly, Ly+1))
    # for (i, j), z in np.ndenumerate(rho):
    #     plt.text(j+0.5, i+0.5, "{:0.3f}".format(z), ha="center", va="center")

    # plt.title(rf"Equilibrium potential of {Lx}x{Ly} 2D lattice with $\beta={beta}$ and $\mu={mu}$")

    # plt.show()

    return rho


def problem10(mu):
    Lx = 4 #Number of sites along x-axis
    Ly = 20 #Number of sites along y-axis

    beta = 1.2 #Inverse temperature beta*epsilon
    beta_epsilon_wall = 4/3

    rho_0 = 0.51 #Initial density
    tol = 1e-12 #Convergence tolerance
    count = 30000 #Upper limit for iterations
    alpha  = 0.01 #Mixing parameter

    #Solve equations iteratively:
    conv = 1
    cnt = 1
    global rho
    rho = rho_0*np.ones([Lx,Ly])
    rho_new = np.zeros([Lx,Ly])

    sol = find_roots(func(beta,mu))

    rho[:,0] = 0
    rho[:,-1] = min(sol)

    while conv >= tol and cnt<count:
        cnt = cnt + 1
        for i in range(Lx):
            for j in range(1,Ly-1):
                #Handle the periodic boundaries for x and y:
                left = np.mod((i-1),Lx)
                right = np.mod((i+1),Lx)
                down = np.mod((j-1),Ly)
                up = np.mod((j+1),Ly)

                v_j = -beta_epsilon_wall*j**(-3)
                rho_new[i,j] = (1 - rho[i,j])*np.exp(beta*(rho[i,down] + rho[i,up] + rho[left,j] + rho[right,j] + (1/4)*(rho[left,down] + rho[right,down] + rho[left,up] + rho[right,up]) + mu - v_j))

        conv = sum(sum((rho - rho_new)**2)); #Compute the convergence parameter.
        rho = alpha*rho_new + (1 - alpha)*rho #Mix the new and old solutions.

        rho[:,0] = 0
        rho[:,-1] = min(sol)



    rho_t = rho.T
    rho_t[0, :] = math.inf


    plt.pcolor(rho_t, vmin=-1, vmax=1)
    cbar = plt.colorbar()
    cbar.set_label(r"Density $\rho$", rotation=270, labelpad=20)

    plt.xlabel("Lattice points")
    plt.ylabel("Lattice points")

    plt.grid(True, linestyle="--")
    plt.xticks(np.linspace(0, Lx, Lx+1))
    plt.yticks(np.linspace(0, Ly, Ly+1))
    for (i, j), z in np.ndenumerate(rho):
        plt.text(i+0.5, j+0.5, "{:0.1f}".format(z), ha="center", va="center", fontsize=6)

    plt.title(rf"Equilibrium potential of {Lx}x{Ly} 2D lattice with $\beta={beta}$ and $\mu={mu}$")

    plt.plot()



def problem11(mu, Ly, beta_epsilon_wall = 1.6):

    beta = 1.2 #Inverse temperature beta*epsilon
    epsilon_wall = beta_epsilon_wall / beta

    rho_0 = 0.51 #Initial density
    tol = 1e-12 #Convergence tolerance
    count = 30000 #Upper limit for iterations
    alpha  = 0.01 #Mixing parameter


    conv = 1
    cnt = 1

    gas_bulk = min(find_roots(func(beta,mu)))

    rho = rho_0*np.ones(Ly)*gas_bulk
    rho[0] = 0
    rho_new = np.zeros(Ly);


    while conv >= tol and cnt<count:
        for j in range(Ly):
            #Handle the periodic boundaries for x and y:
            down = np.mod((j-1),Ly) #j-1, maps -1 to Ly-1
            up = np.mod((j+1),Ly) #j+1, maps Ly to 0
            if j == Ly - 1:
                rho_new[j] = gas_bulk
            elif j == 0:
                rho_new[j] = 0
            else:
                rho_new[j] = (1 - rho[j])*np.exp(beta*((3/2)*(rho[down] + rho[up]) + 2*rho[j] + mu + epsilon_wall*(j)**(-3)))

        conv = sum((rho - rho_new)**2); #Compute the convergence parameter.
        rho = alpha*rho_new + (1 - alpha)*rho #Mix the new and old solutions.
        cnt = cnt + 1

    return rho, gas_bulk


def plotproblem11():
    Ly = 10

    rho1,_ = problem11(-2.67, Ly)
    rho2,_ = problem11(-2.53, Ly)

    plt.plot(np.linspace(0, Ly, Ly), rho1, "blue", marker="o")
    plt.plot(np.linspace(0, Ly, Ly), rho2, "red", marker="x")

    bulkrhopg = np.ones(3)*min(rho1[1:])
    bulkrhopl = np.ones(3)*max(rho2)
    plt.plot([0,1,2], bulkrhopg, color="black", linestyle="--")
    plt.plot([0,1,2], bulkrhopl, color="black", linestyle="--")


    plt.legend([r"$\mu/\epsilon$ = -2.67", r"$\mu/\epsilon$ = -2.53"])

    plt.xlabel(r"Lattice points $y/\sigma$")
    plt.ylabel(r"Density $\rho\sigma^2$")

    plt.title("Density profiles close to the wall")
    plt.show()


def problem12(beta_eps_wall):

    num_vals = 20


    mu_coex = -2.5
    mu_vals = np.linspace(-2.8, mu_coex - 0.01, num_vals)

    Gamma = np.zeros(num_vals)

    for i, mu in enumerate(mu_vals):
        rhoi, rhog = problem11(mu, num_vals, beta_eps_wall)
        Gamma[i] = sum(rhoi - rhog)


    return Gamma


def plotproblem12(beta_vals):

    num_vals = 20
    mu_coex = -2.5

    legend = []

    for beta in beta_vals:
        Gamma = problem12(beta)
        plt.plot(np.linspace(-2.8 - mu_coex, 0, num_vals), Gamma)
        legend.append(rf"$\beta\epsilon_w$ = {beta}")

    plt.legend(legend)

    plt.xlabel(r"$\beta (\mu - \mu_{coex})$")
    plt.ylabel(r"Adsoption $\Gamma$")

    plt.title("Adsorption at the wall")
    plt.show()






def problem14(beta_epsilon_wall):
    Lx = 100 #Number of sites along x-axis
    Ly = 40 #Number of sites along y-axis

    beta = 1.2 #Inverse temperature beta*epsilon
    epsilon_wall = beta_epsilon_wall / beta

    mu_coex = -5/2

    tol = 1e-12 #Convergence tolerance
    count = 30000 #Upper limit for iterations
    alpha  = 0.01 #Mixing parameter

    #Solve equations iteratively:
    conv = 1
    cnt = 1


    minsol = min(find_roots(func(beta,mu_coex)))
    maxsol = max(find_roots(func(beta,mu_coex)))

    global rho
    rho = np.ones([Lx,Ly])*minsol


    shape = np.ones([30,20])*maxsol
    num_particles = 30*20
    rho[34:64,1:21] = shape
    rho[:,0] = maxsol


    rho_initial = rho
    rho_new = np.zeros([Lx,Ly])

    fancy_N = num_particles * 1/np.sum(rho)




    while conv >= tol and cnt<count:
        cnt = cnt + 1
        for i in range(Lx):
            for j in range(1,Ly-1):
                #Handle the periodic boundaries for x and y:
                left = np.mod((i-1),Lx)
                right = np.mod((i+1),Lx)
                down = np.mod((j-1),Ly)
                up = np.mod((j+1),Ly)

                v_j = -epsilon_wall*j**(-3)
                rho_new[i,j] = (1 - rho[i,j])*np.exp(beta*(rho[i,down] + rho[i,up] + rho[left,j] + rho[right,j] + (1/4)*(rho[left,down] + rho[right,down] + rho[left,up] + rho[right,up]) + mu_coex - v_j))

                new_fancy_N = num_particles * 1/np.sum(rho_new)
                if new_fancy_N != fancy_N:
                    print(new_fancy_N)
                    break


        conv = sum(sum((rho - rho_new)**2)); #Compute the convergence parameter.
        rho = alpha*rho_new + (1 - alpha)*rho #Mix the new and old solutions.

        rho[:,0] = maxsol
        rho[:,-1] = minsol


    rho_t = rho.T
    rho_t[0, :] = 0

    jet_cmap = cm.get_cmap("jet", 512)
    piss_cmap = ListedColormap(jet_cmap(np.linspace(0.2, 0.65, 256)))

    plt.pcolor(rho_t, vmin=0, vmax=1, cmap = piss_cmap)
    cbar = plt.colorbar()
    cbar.set_label(r"Density $\rho\sigma^2$", rotation=270, labelpad=20)

    plt.xlabel(r"Lattice points $x/\sigma$")
    plt.ylabel(r"Lattice points $y/\sigma$")

    plt.title(rf"Droplet shape on {Lx}x{Ly} lattice, $\beta\epsilon_w={beta_epsilon_wall}$")

    plt.show()


def main():
    # problem9()


    # problem10(-2.53)


    # plotproblem11()

    # problem12()

    # P12_beta_vals = [0.5, 2.0]
    # plotproblem12(P12_beta_vals)

    # P13_beta_vals = [1.6, 1.7, 1.8]
    # P13_beta_vals16 = [1.6]
    # plotproblem12(P13_beta_vals)
    # plotproblem12(P13_beta_vals16)

    problem14(0.5)




main()
