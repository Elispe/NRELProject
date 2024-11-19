import numpy as np
#Simulation parameters
dt = 0.001
t_max = 100
#Power system parameters
tau_z = 10
Part_Factors = "EquivalentSharing" #EquivalentSharing or ProportionalSharing or RandomDispatch
M_DERS = [40,40,40,40,40,40]
D_DERS = [1.5,-100,-200,-150,-120,-170]
beta = -0.1
Kg = 1/0.05
Tsm = 2
#Attack parameters
tau_d = 2
N0 = 5
tau_a = 4
T0 = 2
#Type of simulation
sim_config = "Transient" #TimeVaryingLoad or PerturbIC or Transient
if sim_config=="TimeVaryingLoad":
      #Time varying Load
      t_cond = np.arange(0,t_max,0.01)
      u_cond = 5*np.sin(0.01*t_cond)
      eps1 = 0
      eps2 = 0
elif sim_config=="PerturbIC":
      #Perturb IC
      t_cond = [0]
      u_cond = [0]
      eps1 = 0.05
      eps2 = 0.0005
else:
      #Transient
      t_cond = [10]
      u_cond = [5]
      eps1 = 0
      eps2 = 0
#Saving results
filename = 'test'