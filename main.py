from buildSystem import LTI_system
from constant import tau_z, M_DERS, D_DERS
import numpy as np

lti = LTI_system(M_DERS, D_DERS, tau_z)

eigval, eigvec = np.linalg.eig(lti.A)

print("rank(A)=", np.linalg.matrix_rank(lti.A))

print("Eigenvalues, real part:")
for e in eigval:
    print(np.real(e))