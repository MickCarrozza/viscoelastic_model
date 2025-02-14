import sys, numpy as np

def rhs_viscoelastic_model (gradv,model,stress):

    # Testing
    if not model.is_created():
        print("Error in rhs_viscoelastic_model: model not created.\n"
              "First create model with model = ModelC(modelnr,lamb,Gmod).\n")
        sys.exit()
        
    if ( isinstance(gradv,  np.ndarray) and gradv.shape != (3, 3) and \
         isinstance(stress, np.ndarray) and gradv.shape != (3, 3) ):
        print("Error in rhs_viscoelastic_model: gradv and stress should be 3x3 numpy matrices.")
        sys.exit()

#   Retrieve material parameters from model
    lamb = model.lamb
    Gmod = model.Gmod
    alpha = model.alpha

# Calculate f-function (see program.py)
    match model.nr:
        case 1:
            f = 0
        case 2:
            f = alpha / Gmod * np.dot(stress,stress)

#   Calculate rhs (see program.py)
    gradvT = np.transpose(gradv)
    rhs = np.dot(gradv,stress) + np.dot(stress,gradvT) - \
          1 / lamb * ( stress + f ) + Gmod * ( gradv + gradvT )

    return rhs