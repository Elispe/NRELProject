# Generate Hurwitz matrix A for 2 interconnected areas, with 1 gen per area.
# Full rank

from constant_4 import *

SBASE = 1000  # MW
VBASE = 380  # kV
ZBASE = (VBASE * 1e3) ** 2 / (SBASE * 1e6)
ohm_per_km = 0.3
omega_s = 2 * np.pi * 60


class LTI():
    def __init__(self):
        Sgen1 = 29000  # MW
        H1 = 5  # s
        Kgov1 = 0.1
        Sgen2 = 31000  # MW
        H2 = 5  # s
        Kgov2 = 0.1
        self.D = {1: 5, 2: 5}
        self.M = {1: 2 * H1 * Sgen1 / SBASE, 2: 2 * H2 * Sgen2 / SBASE}
        self.tau = {1: 2, 2: 2}
        self.Kgov = {1: 1 / Kgov1 * Sgen1 / SBASE, 2: 1 / Kgov2 * Sgen2 / SBASE}
        self.X12 = ohm_per_km * 100 / ZBASE
        print('Line impedance in pu: %.3f' % self.X12)
        self.beta = {1: 1 / Kgov1 * Sgen1 / SBASE + self.D[1], 2: 1 / Kgov2 * Sgen2 / SBASE + self.D[2]}
        self.tau_z = {1: 5, 2: 5}
        self.x_history = list()
        self.u_history = list()
        self.load1 = 30000 / SBASE
        self.load2 = 30000 / SBASE
        self.Psch12 = -1000 / SBASE

    def build_AB(self):
        # x = [w1,Pm1,w2,Pm2,Ptie12,z1,z2]
        self.A = np.array([[-self.D[1] / self.M[1], 1 / self.M[1], 0, 0, -1 / self.M[1], 0, 0],
                           [-self.Kgov[1] / self.tau[1], -1 / self.tau[1], 0, 0, 0, 1 / self.tau[1], 0],
                           [0, 0, -self.D[2] / self.M[2], 1 / self.M[2], 1 / self.M[2], 0, 0],
                           [0, 0, -self.Kgov[2] / self.tau[2], -1 / self.tau[2], 0, 0, 1 / self.tau[2]],
                           [omega_s / self.X12, 0, -omega_s / self.X12, 0, 0, 0, 0],
                           [-self.beta[1] / self.tau_z[1], 0, 0, 0, -1 / self.tau_z[1], 0, 0],
                           [0, 0, -self.beta[2] / self.tau_z[2], 0, 1 / self.tau_z[2], 0, 0],
                           ])

        self.B = np.array([[-1 / self.M[1], 0, 0],
                           [0, 0, 0],
                           [0, -1 / self.M[2], 0],
                           [0, 0, 0],
                           [0, 0, 0],
                           [0, 0, 1 / self.tau_z[1]],
                           [0, 0, -1 / self.tau_z[2]]
                           ])
        eigval = np.linalg.eigvals(self.A)

        print("rank(A): ", np.linalg.matrix_rank(self.A))

        print("Eigenvalues, real part")
        for e in eigval:
            print(np.real(e))
            if np.real(e) > 0:
                print("non-Hurwitz")

    def build_u(self):
        self.u = np.ones((3, 1))
        self.u[0] = self.load1
        self.u[1] = self.load2
        self.u[2] = self.Psch12

    def c2d(self):
        n = self.A.shape[0]  # Assuming square matrices

        M1 = np.eye(n) + 0.5 * dt * self.A
        M2 = np.eye(n) - 0.5 * dt * self.A

        self.Ad = M1 @ np.linalg.inv(M2)
        self.Bd = 0.5 * dt * (np.eye(n) + self.Ad) @ self.B

    def run(self):
        max_iter = int(t_max / dt)
        iter = 0
        self.build_AB()
        self.c2d()
        self.build_u()
        self.x = -np.linalg.inv(self.A) @ self.B @ self.u
        #print(self.x)
        while iter <= max_iter:
            # store results
            self.x_history.append(self.x)
            self.u_history.append(self.u)
            # update u
            if abs(iter * dt - 1) < 0.5 * dt:
                self.load1 += 1000 / SBASE
                self.build_u()
            self.x = self.Ad @ self.x + self.Bd @ self.u
            iter += 1


if __name__ == '__main__':
    LTIObj = LTI()
    LTIObj.run()
