import numpy as np
from network import loadIEEE39bus
from constant import beta


class LTI_system():
    def __init__(self, M_DERS, D_DERS, Tz):
        self.gen, self.line, self.trafos, self.gen_areaA, self.gen_areaB, self.gen_areaC, self.Sbase, self.ng, self.nb, self.ngA, self.ngB, self.ngC, self.nbA, self.nbB, self.nbC = loadIEEE39bus()
        self.M_DERS = M_DERS
        self.D_DERS = D_DERS
        self.Tz = Tz
        self.build_params_mat(self.gen_areaA.items())
        self.build_params_mat(self.gen_areaB.items())
        self.build_params_mat(self.gen_areaC.items())
        self.beta = beta
        self.build_A_mat()
        #self.build_B_mat()

        self.alpha = 1 / self.ngA * np.ones((self.ngA, 1)) #EquivalentSharing


    def build_params_mat(self, gen_area):
        M = self.M_DERS * np.ones((self.nbA, 1))
        D = self.D_DERS * np.ones((self.nbA, 1))
        # self.Tg = np.zeros((self.ng, self.ng))
        # self.Kp = np.zeros((self.ng, 1))
        idx = 0
        for _, value in gen_area:
            M[idx] = 2 * value.H * value.Sgen / self.Sbase
            D[idx] = value.D
            # self.Tg[idx][idx] = value.Tsm
            # self.Kp[idx] = value.Kg * value.Sgen / self.Sbase
            idx += 1
        if gen_area == self.gen_areaA.items():
            self.M_effA = sum(M)
            self.D_netA = sum(D)
        elif gen_area == self.gen_areaB.items():
            self.M_effB = sum(M)
            self.D_netB = sum(D)
        else:
            self.M_effC = sum(M)
            self.D_netC = sum(D)


    def build_A_mat(self):
        self.A = np.block([
            [-self.D_netA / self.M_effA, 1 / self.M_effA * np.ones((1, self.ngA)), np.zeros((1,9)), 1 / self.M_effA, np.zeros((1,5))],
            [0,0,0,0,-self.D_netB / self.M_effB, 1 / self.M_effB * np.ones((1, self.ngB)),0,0,0,0,0, 0, 1 / self.M_effB, 0, 0, 0, 0],
            [0, 0, 0, 0, 0,0,0,0,-self.D_netC / self.M_effC, 1 / self.M_effC * np.ones((1, self.ngC)), 0, 0, 1 / self.M_effC, 0, 0, 0],
            [np.zeros((16,19))]
            #[-np.linalg.inv(self.Tg) @ self.Kp, -np.linalg.inv(self.Tg), np.linalg.inv(self.Tg) @ self.alpha],
            #[1 / self.Tz * self.beta, 1 / self.Tz * np.ones((1, self.ng)), -1 / self.Tz]
        ])


        self.eig, _ = np.linalg.eig(self.A)
        if np.array([e < 0 for e in np.real(self.eig)]).all():
            self.Hurwitz = True
        else:
            self.Hurwitz = False

        if self.Hurwitz:
            print("Hurwitz")

    # def build_B_mat(self):
    #     self.B = np.block([[np.block([-1 / self.M_eff, np.zeros((1, self.ng))])],
    #                        [np.block([np.zeros((self.ng, 1)), np.linalg.inv(self.Tg) @ (
    #                                    np.eye(self.ng) - self.alpha * np.ones((1, self.ng)))])],
    #                        [np.block([np.zeros((1, 1 + self.ng))])]])

