from pysb import *
from pysb.integrate import Solver, get_jacobian_matrix
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
        #
        # Generate the equations for the model
        if not model.odes:
            generate_equations(model)
        print model.species # TODO
        print sympy.Matrix(model.odes).__repr__() # TODO

        # Get the Jacobian matrix for the model
        jac_matrix = get_jacobian_matrix(model)
        print jac_matrix.__repr__() # TODO

        # Calculate the derivatives of each equation with respect to the
        # parameters
        dfdp_matrix = sympy.Matrix(len(model.species), len(model.parameters),
                              lambda i, j: sympy.diff(model.odes[i],
                                                      model.parameters[j].name))
        print dfdp_matrix.__repr__()

        s_matrix = sympy.Matrix(len(model.species), len(model.parameters),
                                lambda i, j: sympy.Symbol('sens%d' % i))
        print s_matrix.__repr__()

        sdot_matrix = jac_matrix * s_matrix + dfdp_matrix
        print sdot_matrix.__repr__()
        return

if __name__ == '__main__':

    # Create an mRNA + Protein model
    Model()
    Monomer('m', []) # mRNA
    Monomer('p', []) # protein
    Rule('mrna_synth', None >> m(), Parameter('k1', 1.))
    Rule('mrna_deg', m() >> None, Parameter('k2', 0.1))
    Rule('protein_synth', m() >> m() + p(), Parameter('k3', 0.02))
    Rule('protein_deg', p() >> None, Parameter('k4', 0.01))
    # Initial conditions
    Observable('m_', m())
    Observable('p_', p())

    t = np.linspace(0, 500, 100)
    sens = Sensitivity(model, t, observables=[m_], parameters=[k1, k2])

    #sol = Solver(model, t)
    #sol.run()
    #plt.ion()
    #plt.figure()
    #plt.plot(t, sol.yobs['m_'], label='mRNA')
    #plt.plot(t, sol.yobs['p_'], label='Protein')
    #plt.legend(loc='lower right')
