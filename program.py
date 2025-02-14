# Program for calculating the time-dependent stress and strain
# of a viscoelastic material under deformation.
# Second-order time integration using Runge-Kutta algorithm (Heun's method).
# The input parameters are provided in input.json.

##############################################################
# Differential equation for the visco-elastic stress tensor:
#              ▽  
# T + lambda * T + f(T) = 2*eta*D

# with T the stress tensor, lambda the relaxation time,
# f a function of T depending on the visco-elastic model,
# eta the viscosity and D = (L + L^T)/2 the rate-of-deformation
# tensor, where L_ij = dv_i/dx_j is the velocity gradient tensor,
# with v the velocity vector and x position vector.
# Furthermore, ▽ is the upper-convected time-derivative defined by:
#  ▽    .
# ( ) = () - L . T - T . L^T,
#        .
# where ( ) the material time derivative.
# This equation can be rewritten in the form:
# .
# T = L.T + T.L^T- 1 / lambda * [ T + f(T) ] + 2*G*D,
# where G = eta / lambda is the shear modulus.
# This can in short be written as:
# .
# T = rhs(L,T)
###############################################################

import json, numpy as np
from rhs_viscoelastic_model import rhs_viscoelastic_model
from modelc import ModelC
import matplotlib.pyplot as plt

# Read parameters from input file
with open("input.json") as file:
    data = json.load(file)  # Load JSON into a dictionary

# Model number
modelnr = data["modelnr"]
# Relaxation time
lamb = data["lamb"]
# Viscosity
Gmod = data["Gmod"]
# Non-linear parameter (for Giesekus)
alpha = data["alpha"]
# Flow number
flownr = data["flownr"]
# Strain rate
rate = data["rate"]
# Time step for time integration
deltat = data["deltat"]
# Number of time steps for time integration
nsteps = data["nsteps"]


# Start simulation
print("\n=========Simulation started=========\n\n"
      "Deforming a visco-elastic material. "
      "Calculating stress and deformation.\n")

# Create material model
print("=======Creating material model======\n")
model = ModelC(modelnr,lamb,Gmod,alpha)
print(model)

# Apply deformation, start time integration
print("========Applying deformation========\n")

# Initialize stress tensor to zero stress
stress_n = np.zeros((3,3))

# Set strain rate tensor
match flownr:
      case 1:
        # Shear flow
        gradv = np.array([[ 0, rate, 0 ],
                          [ 0,    0, 0 ],
                          [ 0,    0, 0 ]]);
      case 2:
        # Extensional flow
        gradv = np.array([[ rate,       0,       0 ],
                          [ 0,    -rate/2,       0 ],
                          [ 0,          0, -rate/2 ]]);

print("Rate-of-deformation tensor:\n", gradv,'\n')


# Start time integration using 2nd order Runge-Kutta algorithm (Heun's method)
print("========Start time integration========\n")

stress = [0] * nsteps
strain = [0] * nsteps
for step in range(0,nsteps,1):
    # Update stress tensor
    # Prediction step
    k1 = rhs_viscoelastic_model(gradv,model,stress_n)
    stress_np1 = stress_n + deltat * k1
    # Actual step
    k2 = rhs_viscoelastic_model(gradv,model,stress_np1)
    stress_np1 = stress_n + deltat / 2 * ( k1 + k2 )

    stress_n = stress_np1
    stress[step] = np.linalg.norm(stress_np1)
    strain[step] = (step+1)*rate*deltat
    
    # Print results
    if ((step+1)%(nsteps/10)==0):
       print("=======================",str(int((step+1)/nsteps*100))+"% completed=====================")
    print("Strain:",str(round(strain[step],3)) +
          ", Stress: " + str(round(stress[step],3)))
    
# Plot results

plt.plot(strain,stress)
plt.ylabel("Stress")
plt.xlabel("Strain")
plt.xlim(left=0)              # Start x-axis from 0
plt.xlim(right=20)              # End x-axis at 10
plt.ylim(bottom=0)            # Start y-axis from 0
plt.show()

