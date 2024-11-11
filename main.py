from buildSystem import LTI_system
from constant import tau_z, M_DERS, D_DERS
import numpy as np

lti = LTI_system(M_DERS, D_DERS, tau_z)

print("M_effA,lti.D_netA:", lti.M_effA,lti.D_netA)
print("M_effB,lti.D_netB:", lti.M_effB,lti.D_netB)
print("M_effC,lti.D_netC:", lti.M_effC,lti.D_netC)

print(lti.A)

eigval, eigvec = np.linalg.eig(lti.A)

print("Positive eigenvalues")
for e in eigval:
    if np.real(e)>0:
        print(e)

print("Matrix A rank", np.linalg.matrix_rank(lti.A))


