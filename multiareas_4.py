# Generate Hurwitz matrix A for 4 interconnected areas, with 1 gen per area.
# Not full rank

from constant_4 import *

SBASE = 1000  # MW
VBASE = 380  # kV
ZBASE = (VBASE * 1e3) ** 2 / (SBASE * 1e6)
ohm_per_km = 0.3
omega_s = 2 * np.pi * 60


class LTI():
    def __init__(self):
        Sgen1 = 15000  # MW
        H1 = 5  # s
        Kgov1 = 0.05
        Sgen2 = 16000  # MW
        H2 = 3  # s
        Kgov2 = 0.05
        Sgen3 = 15000  # MW
        H3 = 6  # s
        Kgov3 = 0.05
        Sgen4 = 16000  # MW
        H4 = 4  # s
        Kgov4 = 0.1
        self.D = {1: 5, 2: 4, 3: 3, 4: 6}
        self.M = {1: 2 * H1 * Sgen1 / SBASE, 2: 2 * H2 * Sgen2 / SBASE, 3: 2 * H3 * Sgen3 / SBASE,
                  4: 2 * H4 * Sgen4 / SBASE}
        self.tau = {1: 4, 2: 2, 3: 2.5, 4: 3.5}
        self.Kgov = {1: 1 / Kgov1 * Sgen1 / SBASE, 2: 1 / Kgov2 * Sgen2 / SBASE, 3: 1 / Kgov3 * Sgen3 / SBASE,
                     4: 1 / Kgov4 * Sgen4 / SBASE}
        self.X12 = ohm_per_km * 20 / ZBASE
        self.X13 = ohm_per_km * 100 / ZBASE
        self.X34 = ohm_per_km * 20 / ZBASE
        print('Line impedance in pu: %.3f' % self.X12)
        self.beta = {1: 1 / Kgov1 * Sgen1 / SBASE + self.D[1], 2: 1 / Kgov2 * Sgen2 / SBASE + self.D[2],
                     3: 1 / Kgov3 * Sgen3 / SBASE + self.D[3], 4: 1 / Kgov4 * Sgen4 / SBASE + self.D[4]}
        self.tau_z = {1: 5, 2: 5, 3: 5, 4: 5}
        self.x_history = list()
        self.u_history = list()
        self.load1 = 17000 / SBASE
        self.load2 = 15000 / SBASE
        self.load3 = 14000 / SBASE
        self.load4 = 16000 / SBASE
        self.Psch12 = 0 / SBASE
        self.Psch13 = -1000 / SBASE
        self.Psch34 = 0 / SBASE

    def build_AB(self):
        # x =  [w1,Pm1,w2,Pm2,w3,Pm3,w4,Pm4,Ptie1,Ptie2,Ptie3,Ptie4,z1,z2,z3,z4]
        self.A = np.array(
            [[-self.D[1] / self.M[1], 1 / self.M[1], 0, 0, 0, 0, 0, 0, -1 / self.M[1], 0, 0, 0, 0, 0, 0, 0],
             [-self.Kgov[1] / self.tau[1], -1 / self.tau[1], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 / self.tau[1], 0, 0, 0],
             [0, 0, -self.D[2] / self.M[2], 1 / self.M[2], 0, 0, 0, 0, 0, -1 / self.M[2], 0, 0, 0, 0, 0, 0],
             [0, 0, -self.Kgov[2] / self.tau[2], -1 / self.tau[2], 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 / self.tau[2], 0, 0],
             [0, 0, 0, 0, -self.D[3] / self.M[3], 1 / self.M[3], 0, 0, 0, 0, -1 / self.M[3], 0, 0, 0, 0, 0],
             [0, 0, 0, 0, -self.Kgov[3] / self.tau[3], -1 / self.tau[3], 0, 0, 0, 0, 0, 0, 0, 0, 1 / self.tau[3], 0],
             [0, 0, 0, 0, 0, 0, -self.D[4] / self.M[4], 1 / self.M[4], 0, 0, 0, -1 / self.M[4], 0, 0, 0, 0],
             [0, 0, 0, 0, 0, 0, -self.Kgov[4] / self.tau[4], -1 / self.tau[4], 0, 0, 0, 0, 0, 0, 0, 1 / self.tau[4]],
             [omega_s * (1 / self.X12 + 1 / self.X13), 0, -omega_s / self.X12, 0, -omega_s / self.X13, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0],
             [-omega_s / self.X12, 0, omega_s * (1 / self.X12), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
             [-omega_s / self.X13, 0, 0, 0, omega_s * (1 / self.X13 + 1 / self.X34), 0, -omega_s / self.X34, 0, 0, 0, 0,
              0, 0, 0, 0, 0],
             [0, 0, 0, 0, -omega_s / self.X34, 0, omega_s * (1 / self.X34), 0, 0, 0, 0, 0, 0, 0, 0, 0],
             [-self.beta[1] / self.tau_z[1], 0, 0, 0, 0, 0, 0, 0, -1 / self.tau_z[1], 0, 0, 0, 0, 0, 0, 0],
             [0, 0, -self.beta[2] / self.tau_z[2], 0, 0, 0, 0, 0, 0, -1 / self.tau_z[2], 0, 0, 0, 0, 0, 0],
             [0, 0, 0, 0, -self.beta[3] / self.tau_z[3], 0, 0, 0, 0, 0, -1 / self.tau_z[3], 0, 0, 0, 0, 0],
             [0, 0, 0, 0, 0, 0, -self.beta[4] / self.tau_z[4], 0, 0, 0, 0, -1 / self.tau_z[4], 0, 0, 0, 0],
             ])
        # u = [Pload1, Pload2, Pload3, Pload4, Psch12, Psch13, Psch34]
        self.B = np.array([[-1 / self.M[1], 0, 0, 0, 0, 0, 0],
                           [0, 0, 0, 0, 0, 0, 0],
                           [0, -1 / self.M[2], 0, 0, 0, 0, 0],
                           [0, 0, 0, 0, 0, 0, 0],
                           [0, 0, -1 / self.M[3], 0, 0, 0, 0],
                           [0, 0, 0, 0, 0, 0, 0],
                           [0, 0, 0, -1 / self.M[4], 0, 0, 0],
                           [0, 0, 0, 0, 0, 0, 0],
                           [0, 0, 0, 0, 0, 0, 0],
                           [0, 0, 0, 0, 0, 0, 0],
                           [0, 0, 0, 0, 0, 0, 0],
                           [0, 0, 0, 0, 0, 0, 0],
                           [0, 0, 0, 0, 1 / self.tau_z[1], 1 / self.tau_z[1], 0],
                           [0, 0, 0, 0, 0, -1 / self.tau_z[2], 0],
                           [0, 0, 0, 0, 0, -1 / self.tau_z[3], 1 / self.tau_z[3]],
                           [0, 0, 0, 0, 0, 0, -1 / self.tau_z[4]],
                           ])

        print("rank(A): ", np.linalg.matrix_rank(self.A))

        eigval = np.linalg.eigvals(self.A)

        print("Eigenvalues, real part")
        for e in eigval:
            print(np.real(e))
            if np.real(e) > 0:
                print("non-Hurwitz")

        print("Matrix A rank", np.linalg.matrix_rank(self.A))

    def build_u(self):
        self.u = np.ones((7, 1))
        self.u[0] = self.load1
        self.u[1] = self.load2
        self.u[2] = self.load3
        self.u[3] = self.load4
        self.u[4] = self.Psch12
        self.u[5] = self.Psch13
        self.u[6] = self.Psch34

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
