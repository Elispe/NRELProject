import numpy as np
from network import loadIEEE39bus
from constant import beta, omega_s, ohm_per_km, ZBASE
import os



class LTI_system():
    def __init__(self, M_DERS, D_DERS, Tz):
        self.gen, self.line, self.trafos, self.gen_areaA, self.gen_areaB, self.gen_areaC, self.Sbase, self.ng, self.nb, self.nb_single_area = loadIEEE39bus()
        self.M_DERS = M_DERS
        self.D_DERS = D_DERS
        self.Tz = Tz
        self.beta = beta
        self.X = ohm_per_km * 100 / ZBASE  # X - admittance value
        self.ngA = len(self.gen_areaA)
        self.ngB = len(self.gen_areaB)
        self.ngC = len(self.gen_areaC)
        self.build_params_mat(self.gen_areaA)
        self.build_params_mat(self.gen_areaB)
        self.build_params_mat(self.gen_areaC)
        self.build_A_mat()
        #self.build_B_mat()

    def build_params_mat(self, gen_area):
        ng_single_area = len(gen_area)
        M = self.M_DERS * np.ones((self.nb_single_area, 1))
        D = self.D_DERS * np.ones((self.nb_single_area, 1))
        Tg = np.zeros((ng_single_area,ng_single_area))
        Kp = np.zeros((ng_single_area, 1))
        alpha = 1 / ng_single_area * np.ones((ng_single_area, 1)) # EquivalentSharing
        idx = 0
        for _, value in gen_area.items():
            M[idx] = 2 * value.H * value.Sgen / self.Sbase
            D[idx] = value.D
            Tg[idx][idx] = value.Tsm
            Kp[idx] = value.Kg * value.Sgen / self.Sbase
            idx += 1
        if gen_area == self.gen_areaA:
            self.M_effA = sum(M)
            self.D_netA = sum(D)
            self.TgA=Tg
            self.KpA=Kp
            self.alphaA = alpha
        elif gen_area == self.gen_areaB:
            self.M_effB = sum(M)
            self.D_netB = sum(D)
            self.TgB=Tg
            self.KpB=Kp
            self.alphaB = alpha
        else:
            self.M_effC = sum(M)
            self.D_netC = sum(D)
            self.TgC=Tg
            self.KpC=Kp
            self.alphaC = alpha

    def build_A_mat(self):
        self.A = np.block([
            [-self.D_netA / self.M_effA, 1 / self.M_effA * np.ones((1, self.ngA)), np.zeros((1,9)), 1 / self.M_effA, np.zeros((1,5))],
            [-np.linalg.inv(self.TgA) @ self.KpA, -np.linalg.inv(self.TgA), np.zeros((3, 12)), np.linalg.inv(self.TgA) @ self.alphaA,np.zeros((3, 2))],
            [np.zeros((1,4)), -self.D_netB / self.M_effB, 1 / self.M_effB * np.ones((1, self.ngB)), np.zeros((1,6)),
             1 / self.M_effB, np.zeros((1,4))],
            [np.zeros((3, 4)),-np.linalg.inv(self.TgB) @ self.KpB, -np.linalg.inv(self.TgB), np.zeros((3, 9)),
             np.linalg.inv(self.TgB) @ self.alphaB, np.zeros((3, 1))],
            [np.zeros((1,8)), -self.D_netC / self.M_effC, 1 / self.M_effC * np.ones((1, self.ngC)), np.zeros((1,2)), 1 / self.M_effC, np.zeros((1,3))],
            [np.zeros((4, 8)), -np.linalg.inv(self.TgC) @ self.KpC, -np.linalg.inv(self.TgC), np.zeros((4, 5)),
             np.linalg.inv(self.TgC) @ self.alphaC],
            [omega_s / (2 * self.X), np.zeros((1,self.ngA)),-omega_s / self.X, np.zeros((1,self.ngB)), -omega_s / self.X, np.zeros((1,self.ngC+6))],
            [-omega_s / self.X, np.zeros((1,self.ngA)),omega_s / (2 * self.X), np.zeros((1,self.ngB)),-omega_s / self.X, np.zeros((1,self.ngC+6))],
            [-omega_s / self.X, np.zeros((1,self.ngA)),-omega_s / self.X, np.zeros((1,self.ngB)),omega_s / (2 * self.X), np.zeros((1,self.ngC+6))],
            [-1 / self.Tz * self.beta, np.zeros((1, self.ng+2)), -1 / self.Tz,np.zeros((1, 5))],
            [np.zeros((1, self.ngA + 1)),-1 / self.Tz * self.beta, np.zeros((1, self.ng + 2-self.ngA)), -1 / self.Tz, np.zeros((1, 4))],
            [np.zeros((1, self.ngA +  self.ngB + 2)),-1 / self.Tz * self.beta, np.zeros((1, self.ngC + 2)), -1 / self.Tz, np.zeros((1, 3))]
        ])

        self.eig, _ = np.linalg.eig(self.A)
        if np.array([e < 0 for e in np.real(self.eig)]).all():
            self.Hurwitz = True
        else:
            self.Hurwitz = False

        if self.Hurwitz:
            print("Hurwitz")
        else:
            print("non Hurwitz")

        # save matrix A
        if not os.path.exists('matrixA'):
            os.makedirs('matrixA')
            np.savetxt('matrixA//A.txt', self.A)

    # def build_B_mat(self):
    #     self.B = np.block([[np.block([-1 / self.M_eff, np.zeros((1, self.ng))])],
    #                        [np.block([np.zeros((self.ng, 1)), np.linalg.inv(self.Tg) @ (
    #                                    np.eye(self.ng) - self.alpha * np.ones((1, self.ng)))])],
    #                        [np.block([np.zeros((1, 1 + self.ng))])]])

