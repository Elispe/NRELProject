# Parameters of the power system
# File built out of the Julia code
# 3 areas, 1 conventional gen per area
# Gives non-Hurwitz A, I think it is wrong
from math import pi
import numpy as np
from pipenv.cli.command import update

SBASE = 1000  # MW
VBASE = 380  # kV
ZBASE = (VBASE * 1e3) ** 2 / (SBASE * 1e6)
ohm_per_km = 0.3
omega_s = 2 * pi * 60

#### Topology of the multi-area power system ####

# 3 area system
G = np.array([[1, 1, 1], [1, 1, 0], [1, 0, 1]])

# Number of areas
N = len(G)

# Number of defense matrices
m_v = N

# Number of attack matrices
m_w = N - 1

# Parameters
Sgen = 30000 * np.ones(N) + np.random.choice([1000, -1000, 500, -500], size=N)
H = 5 * np.ones(N) + np.random.choice([1, -1], size=N)
D = 5 * np.ones(N) + np.random.choice([0.1, -0.1, 0.2, -0.2], size=N)
M = np.array([2 * h * sgen / SBASE for (h, sgen) in zip(H, Sgen)])
tau = 2 * np.ones(N) + np.random.choice([0.2, -0.2, 0.4, -0.4], size=N)
Kgovs = 0.1 * np.ones(N) + np.random.choice([0.01, -0.01, 0.02, -0.02], size=N)
Kgov = np.array([1 / kgov * sgen / SBASE for (kgov, sgen) in zip(Kgovs, Sgen)])
X = ohm_per_km * 100 / ZBASE  # X - admittance value
beta = np.array([1 / kgov * sgen / SBASE + d for (kgov, sgen, d) in zip(Kgovs, Sgen, D)])
tau_z = 5 * np.ones(N) + np.random.choice([0.1, -0.1, 0.2, -0.2], size=N)


# Function for forming a diagonal block
def form_diag_block(i):
    # get the number of neighbors of node i, i.e. its degree:
    degree = sum(G[i, :]) - 1
    mat = np.array([[-D[i] / M[i], 1 / M[i], -1 / M[i], 0], [-Kgov[i] / tau[i], -1 / tau[i], 0, 1 / tau[i]],
                    [omega_s / (degree * X), 0, 0, 0], [-beta[i] / tau_z[i], 0, -1 / tau_z[i], 0]])
    return mat


# System matrix block
sys_mat_block = {}
for i in range(N):
    for j in range(N):
        if i == j:
            key = (i, j)
            mat = form_diag_block(i)
        elif i != j and G[i, j] == 0:
            key = (i, j)
            mat = np.zeros((4, 4))
        else:
            key = (i, j)
            mat = np.zeros((4, 4))
            mat[2, 0] = -omega_s / X

        sys_mat_block.update({key: mat})

# Defense matrices
Ds = {}
last_i = 0
for i in range(m_v):
    mat = np.zeros((4 * N, 4 * N))
    mat[last_i, last_i] = -1 / M[i]
    last_i += 4
    Ds.update({i: mat})

# Attack matrices
Ks = {}
last_i = 0
for i in range(m_w):
    mat = np.zeros((4 * N, 4 * N))
    mat[last_i, last_i] = -1 / M[i]
    last_i += 4
    Ks.update({i: mat})

eps = 0.0001

# Ai matrices for distributed optimization
area_mats = {}
for n in range(N):
    mat = np.zeros((4 * N, 4 * N))
    last_i = 0
    for i in range(N):
        last_j = 0
        for j in range(N):
            if i == n and j == n:
                mat[last_i:last_i+4,last_j:last_j+4] = sys_mat_block[(i, j)]
            elif i == n or j == n and G[i, j] == 1:
                mat[last_i:last_i + 4, last_j:last_j + 4] = sys_mat_block[(i, j)]/2
            else:
                mat[last_i:last_i + 4, last_j:last_j + 4] = np.zeros((4,4))
            last_j += 4
        last_i += 4
    area_mats.update({n: mat})

# The A matrix
A = np.zeros((4 * N, 4 * N))
last_i = 0
for i in range(N):
    last_j = 0
    for j in range(N):
        if G[i,j]==1:
            A[last_i:last_i+4,last_j:last_j+4] = sys_mat_block[(i, j)]
        else:
            mat[last_i:last_i + 4, last_j:last_j + 4] = np.zeros((4,4))
        last_j+=4
    last_i+=4

# Make A Hurwitz
# A -= 5*np.identity(4*N)

eigval, eigvec = np.linalg.eig(A)

print("Eigenvalues")
for e in eigval:
    print(e)

print("Matrix A rank", np.linalg.matrix_rank(A))

