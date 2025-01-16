import numpy as np

BOLTZMANN = 1.38e-23 # J/K
PLANCK = 6.63e-34 # J/Hz

class Vec2D:
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def is_nearest_neighbor(self, vec):
        return abs(self.x - vec.x) + abs(self.y - vec.y) == 1

    def is_next_nearest(self, vec):


def density(beta, potential, phi, chem_pot):
    expon = np.exp(-beta * (potential + phi - chem_pot))
    return expon / (1+expon)

def ref_sys_part_fn(beta, potential, phi, chem_pot):
    expon = np.exp(-beta * (potential + phi - chem_pot))
    return np.prod(1 + expon)

def interaction_sum(eps, num_sites):
    eps_nn = eps
    eps_nnn = eps/4

    sum = 0

    for i in range(num_sites):
        for j in range(num_sites):
            if abs(i - j) == 1: # nearest neighbor
                sum += eps_nn
            elif abs(i - j) == 2: # next-nearest neighbor
                sum += eps_nnn

    return sum

def main():
    num_sites = 10
    chem_pot = 0 # unimpl
    eps = 1

    temp = 293 # K
    beta = 1/(BOLTZMANN * temp) # 1/J

    potential = np.zeros(num_sites)
    phi = np.zeros(num_sites)

    ref_sys = ref_sys_part_fn(beta, potential, phi, chem_pot)
    print("Ξ:", ref_sys)

    rho = density(beta, potential, phi, chem_pot)
    print("⍴:", rho)

    zeta_LJ = interaction_sum(eps, num_sites)
    print("ζ_LJ:", zeta_LJ)

main()
