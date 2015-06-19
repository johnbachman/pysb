from pysb import *
from pysb.integrate import Solver, get_jacobian_matrix, get_model_odes_as_str,\
                           eqn_substitutions
import numpy as np
from matplotlib import pyplot as plt
from pysb.bng import generate_equations
import sympy

class Sensitivity(object):
    """Calculate the sensitivity timecourses, dy(t)/dp.

    Parameters
    ----------
    model : pysb.core.Model
        The model for calculating the sensitivities.
    t : np.array
        Array of timepoints.
    observables : list of pysb.core.Observable
        The observables to calculate the sensitivities of.
    parameters : list of pysb.core.Parameter
        The parameters to perform sensitivity analysis against.
    """

    def __init__(self, model, t, observables=None, parameters=None):
        # For each column s_i(t) in the sensitivity matrix,
        # s_i'(t) = J * s_i(t) + dy'(t)/dp
        # where J is the Jacobian matrix and dy'(t)/dp_i is the derivative
        # of the ODE RHS equations with respect to the parameters.
        #
        # First, we generate the sensitivity matrix for the whole system.
        # We may be able to reduce this later on.

        # Generate the equations for the model
        if not model.odes:
            generate_equations(model)
        # Get the model ODE RHS as a string
        code_eqs = get_model_odes_as_str(model)
        code_eqs += '\n'

        # SENSITIVITY MATRIX -----------------------------------------
        # Get the Jacobian matrix for the model
        jac_matrix = get_jacobian_matrix(model)
        # Calculate the derivatives of each equation with respect to the
        # parameters
        dfdp_matrix = sympy.Matrix(len(model.species), len(model.parameters),
                              lambda i, j: sympy.diff(model.odes[i],
                                                      model.parameters[j].name))
        s_matrix = sympy.Matrix(len(model.species), len(model.parameters),
                              lambda i, j: sympy.Symbol('sens_%d_%d' % (i, j)))
        # The sensitivity matrix: s_i'(t) = J * s_i(t) + dy'(t)/dp
        sdot_matrix = jac_matrix * s_matrix + dfdp_matrix
        # Prepare the sensitivity entries as a list of strings
        sdot_eqs_list = []
        (nrows, ncols) = sdot_matrix.shape
        for i in range(nrows):
            for j in range(ncols):
                entry = sdot_matrix[i, j]
                # Skip zero entries in the sensitivity matrix
                if entry == 0:
                    continue
                sdot_eq_str = 'sdot[%d, %d] = %s;' % (i, j, sympy.ccode(entry))
                sdot_eqs_list.append(sdot_eq_str)
        code_eqs += eqn_substitutions(model, '\n'.join(sdot_eqs_list))

        print code_eqs
        # Next:
        # 3. Create an inline function with macros for YDOT1 and SDOT2
        # 4. Create run method, analogous to Solver
        #    (alternatively, create Solver object and fill in functions, fields)
        # 5. Integrate and return results, plot
        #
        # Q. Make faster by calculating the Jacobian of the extended system?
        #    Probably not, as this would be huge, and CVODES wouldn't use it.
        # 6. Simplify sensitivity equations by only including 
        # 5. 

        # Print TODO
        print model.species
        print model.parameters
        print sympy.Matrix(model.odes).__repr__()
        print jac_matrix.__repr__()
        print dfdp_matrix.__repr__()
        print s_matrix.__repr__()
        print sdot_matrix.__repr__()


if __name__ == '__main__':

    # Create an mRNA + Protein model
    Model()
    """
    Monomer('m', []) # mRNA
    Monomer('p', []) # protein
    Rule('mrna_synth', None >> m(), Parameter('k1', 1.))
    Rule('mrna_deg', m() >> None, Parameter('k2', 0.1))
    Rule('protein_synth', m() >> m() + p(), Parameter('k3', 0.02))
    Rule('protein_deg', p() >> None, Parameter('k4', 0.01))
    # Initial conditions
    Observable('m_', m())
    Observable('p_', p())
    """
    Monomer('A', [])
    Rule('A_deg', A() >> None, Parameter('k', 0.01))
    Initial(A(), Parameter('A_0'))
    Observable('A_', A())

    t = np.linspace(0, 500, 100)
    #sens = Sensitivity(model, t, observables=[m_], parameters=[k1, k2])
    sens = Sensitivity(model, t)

    #sol = Solver(model, t)
    #sol.run()
    #plt.ion()
    #plt.figure()
    #plt.plot(t, sol.yobs['m_'], label='mRNA')
    #plt.plot(t, sol.yobs['p_'], label='Protein')
    #plt.legend(loc='lower right')
