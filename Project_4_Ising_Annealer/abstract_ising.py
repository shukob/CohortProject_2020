import numpy as np
from numba import jit
from numba.experimental import jitclass


class AbstractIsing:

    def __init__(self, method='metropolis'):
        self.method = method

    def energy(self):
        pass
        """Returns the energy of the current spin configuration"""

    def energy_diff(self, *coords):
        pass
        """Returns the energy difference resulting from flipping the site at the given coordinates"""

    def rand_site(self):
        pass
        """Selects a site in the lattice at random"""

    def local_energy_difference(self, *site):
        pass
        """Local energy difference E(1) - E(0) """


    def mc_step(self, T):
        if self.method == 'metropolis':
            return self.mc_step_metropolis(T)
        elif self.method == 'heat_bath':
            return self.mc_step_heat_bath(T)
        elif self.method == 'swendsen_wang':
            pass

    def mc_step_metropolis(self, T):
        """Performs a full update of the given Ising model using the Metropolis-Hastings algorithm"""
        current_energy = self.energy()
        for _ in range(self.num_spins):
            site = self.rand_site()
            dE = self.energy_diff(*site)

            if (dE < 0) or (np.random.rand() < np.exp(-dE / T)):
                current_energy += dE
                self.spins[site] *= -1

        return current_energy

    def mc_step_swendsen_wang(self, T):

        pass

    def mc_step_heat_bath(self, T):
        current_energy = self.energy()
        for _ in range(self.num_spins):
            site = self.rand_site()
            dE = self.local_energy_difference(*site)
            if (np.random.rand() < np.exp(-dE / T)):
                if self.spins[site] == 1:
                    current_energy -= dE
                self.spins[site] = -1
            else:
                if self.spins[site] == -1:
                    current_energy += dE
                self.spins[site] = 1
        return current_energy


## For Swendsen-Wang algorithm
@jit(nopython=True)
def get_cluster_num(num, cluster_num):
    if num == cluster_num[num]:
        return num
    else:
        return get_cluster_num(cluster_num[num], cluster_num)


@jit(nopython=True)
def update_cluster_num(ni, nj, cluster_num):
    ci = get_cluster_num(ni, cluster_num)
    cj = get_cluster_num(nj, cluster_num)
    if ci < cj:
        cluster_num[cj] = ci
    else:
        cluster_num[ci] = cj


@jit(nopython=True)
def make_bond(S, prob, L):
    N = L ** 2
    bond = np.zeros((L, L, 2), dtype=np.int64)
    ran = np.random.rand(2 * N).reshape(L, L, 2)
    for ix in range(L):
        for iy in range(L):
            if S[ix, iy] * S[(ix + 1) % L, iy] > 0:  ## spins are parallel
                if ran[ix, iy, 0] < prob:
                    bond[ix, iy, 0] = 1
            if S[ix, iy] * S[ix, (iy + 1) % L] > 0:  ## spins are parallel
                if ran[ix, iy, 1] < prob:
                    bond[ix, iy, 1] = 1
    return bond


@jit(nopython=True)
def make_cluster(bond, L):
    N = L ** 2
    cluster_num = np.arange(N)
    for ix in range(L):
        for iy in range(L):
            if bond[ix, iy, 0] > 0:  ## connected
                ni = ix + iy * L
                nj = (ix + 1) % L + iy * L
                update_cluster_num(ni, nj, cluster_num)
            if bond[ix, iy, 1] > 0:  ## connected
                ni = ix + iy * L
                nj = ix + (iy + 1) % L * L
                update_cluster_num(ni, nj, cluster_num)
        ## count total cluster number

    cluster_num_count = np.zeros(N, dtype=np.int64)
    for i in range(N):
        nc = get_cluster_num(i, cluster_num)
        cluster_num[i] = nc
        cluster_num_count[nc] += 1

    total_cluster_num = 0

    true_cluster_num = np.zeros(N, dtype=np.int64)
    true_cluster_num_count = np.zeros(N, dtype=np.int64)
    for nc in range(N):
        if cluster_num_count[nc] > 0:
            true_cluster_num[nc] = total_cluster_num
            true_cluster_num_count[total_cluster_num] = cluster_num_count[nc]
            total_cluster_num += 1

    for i in range(N):
        cluster_num[i] = true_cluster_num[cluster_num[i]]
    return cluster_num.reshape(L, L), true_cluster_num_count[:total_cluster_num]


@jit(nopython=True)
def flip_spin(S, cluster_num, cluster_num_count, flip, L):
    total_cluster_num = cluster_num_count.shape[0]
    ran = np.random.rand(total_cluster_num)
    spin_direction = np.zeros(total_cluster_num, dtype=np.int64)

    for i in range(total_cluster_num):
        if ran[i] < 1.0 / (1.0 + np.exp(-2.0 * flip * cluster_num_count[i])):
            spin_direction[i] = 1
        else:
            spin_direction[i] = -1

    for ix in range(L):
        for iy in range(L):
            S[ix, iy] = spin_direction[cluster_num[ix, iy]]


def Swendsen_Wang(S, prob, flip, L):
    N = L ** 2

    ## make bond configulations
    bond = make_bond(S, prob, L)

    ## make clusters
    cluster_num, cluster_num_count = make_cluster(bond, L)
    ## update spin

    flip_spin(S, cluster_num, cluster_num_count, flip, L)
    ## for imporoved estimator
    Nc2 = np.sum(cluster_num_count.astype(float) ** 2)

    return S, Nc2 / N ** 2