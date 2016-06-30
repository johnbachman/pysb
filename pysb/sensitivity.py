# TODO TODO TODO
# 4. Add back in calculations of sens.yobs and sens.yexpr!!!
# 4.5 Calculate derivs of initial conditions for exprs
# 6. Simplify sensitivity equations by only including 
# TODO TODO TODO

from pysb import *
from pysb.integrate import Solver, get_jacobian_matrix, get_model_odes_as_str,\
                           eqn_substitutions, default_integrator_options
import numpy as np
from matplotlib import pyplot as plt
from pysb.bng import generate_equations
import sympy
import re
import itertools
from scipy.integrate import ode
from scipy.weave import inline

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

    def __init__(self, model, tspan, observables=None, parameters=None,
                 integrator='vode', **integrator_options):
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
        # A dict to keep track of the ydot vector index for the sensitivity
        # matrix entry  at i, j
        self.sdot_ix_map = {}
        # We will start indexing the extended set of ODEs where we left off
        # with the RHS odes:
        sdot_cur_ix = len(model.odes)
        for i in range(nrows):
            for j in range(ncols):
                entry = sdot_matrix[i, j]
                # TODO TODO TODO
                # Can't skip zero entries easily, because other entries
                # in the sensitivity matrix may refer to s_ij that are zero.
                # Ideally, the fix would be to figure out which s_ij are zero,
                # and then substitute those 0s into the equations so that some
                # of the terms would disappear.
                # Skip zero entries in the sensitivity matrix
                #if entry == 0:
                #    continue
                # TODO TODO TODO
                sdot_eq_str = 'ydot[%d] = %s;' % (sdot_cur_ix,
                                                  sympy.ccode(entry))
                #sdot_eq_str = 'sdot[%d, %d] = %s;' % (i, j, sympy.ccode(entry))
                sdot_eqs_list.append(sdot_eq_str)
                # Keep track of the vector index for this (i, j) entry
                self.sdot_ix_map[(i, j)] = sdot_cur_ix
                # Increment the sdot index
                sdot_cur_ix += 1
        # We'll need to know how many equations there are in total to
        # initialize the y and ydot vectors correctly
        num_eqns = sdot_cur_ix
        # Replace sens_xxx_xxx with y[i], where i is the appropriate
        # index for this ydot equation:
        sdot_eqs = '\n'.join(sdot_eqs_list)
        sdot_eqs = re.sub(r'sens_(\d+)_(\d+)',
                          lambda m: 'y[%s]' % \
                             self.sdot_ix_map[(int(m.group(1)), int(m.group(2)))],
                          sdot_eqs)
        code_eqs += eqn_substitutions(model, sdot_eqs)

        # Test inline
        Solver._test_inline()
        # If we can't use weave.inline to run the C code, compile it as Python
        # code instead for use with exec. Note: C code with array indexing,
        # basic math operations, and pow() just happens to also be valid
        # Python.  If the equations ever have more complex things in them, this
        # might fail.
        if not Solver._use_inline:
            code_eqs_py = compile(code_eqs, '<%s odes>' % model.name, 'exec')
        else:
            for arr_name in ('ydot', 'y', 'p'):
                macro = arr_name.upper() + '1'
                code_eqs = re.sub(r'\b%s\[(\d+)\]' % arr_name,
                                  '%s(\\1)' % macro, code_eqs)
        print "---- code_eqs"
        print code_eqs
        def rhs(t, y, p):
            ydot = self.ydot
            # note that the evaluated code sets ydot as a side effect
            if Solver._use_inline:
                inline(code_eqs, ['ydot', 't', 'y', 'p']);
            else:
                exec code_eqs_py in locals()
            return ydot

        # FIXME FIXME FIXME
        # From here down is nearly a complete copy paste of Solver;
        # The difference is that the y and ydot arrays are of size
        # num_eqns rather than len(model.species)
        # FIXME FIXME FIXME
        # build integrator options list from our defaults and any kwargs passed
        # to this function
        options = {}
        try:
            options.update(default_integrator_options[integrator])
        except KeyError as e:
            pass
        options.update(integrator_options)

        # Initialize variables
        self.model = model
        self.tspan = tspan
        self.y = np.ndarray((len(tspan), num_eqns))
        self.ydot = np.ndarray(num_eqns)
        # Initialize record array for observable timecourses
        if len(model.observables):
            self.yobs = np.ndarray(len(tspan), zip(model.observables.keys(),
                                                      itertools.repeat(float)))
            # FIXME FIXME FIXME Could never get the view to work correctly
            self.yobs_sens = np.ndarray((len(model.parameters), len(tspan)),
                           zip(model.observables.keys(), itertools.repeat(float)))
        # TODO TODO TODO What is purpose of this section?
        else:
            self.yobs = np.ndarray((len(tspan), 0))
            self.yobs_sens = np.ndarray((0, len(model.parameters), len(tspan)))
        # Initialize view of observables record array
        self.yobs_view = self.yobs.view(float).reshape(len(self.yobs), -1)
        # FIXME FIXME FIXME Could never get the view to work correctly
        #self.yobs_sens_view = self.yobs_sens.view(float).reshape(
        #              len(model.observables), len(model.parameters), len(tspan))
        # Initialize array for expression timecourses
        exprs = model.expressions_dynamic()
        if len(exprs):
            self.yexpr = np.ndarray(len(tspan), zip(exprs.keys(),
                                                       itertools.repeat(float)))
        else:
            self.yexpr = np.ndarray((len(tspan), 0))
        # Initialize an instance of scipy.integrate.ode
        self.integrator = ode(rhs).set_integrator(integrator, **options)

        # Print TODO
        print "---- model.species"
        print model.species
        print "---- model.parameters"
        print model.parameters
        print "---- model.odes"
        print sympy.Matrix(model.odes).__repr__()
        print "---- jac_matrix"
        print jac_matrix.__repr__()
        print "---- dfdp_matrix"
        print dfdp_matrix.__repr__()
        print "---- S_matrix"
        print s_matrix.__repr__()
        print "---- Sdot_matrix"
        print sdot_matrix.__repr__()


    # FIXME FIXME Copy-pasted from Solver
    def run(self, param_values=None, y0=None):
        """Perform an integration.

        Returns nothing; access the Solver object's ``y``, ``yobs``, or
        ``yobs_view`` attributes to retrieve the results.

        Parameters
        ----------
        param_values : vector-like, optional
            Values to use for every parameter in the model. Ordering is
            determined by the order of model.parameters. If not specified,
            parameter values will be taken directly from model.parameters.
        y0 : vector-like, optional
            Values to use for the initial condition of all species. Ordering is
            determined by the order of model.species. If not specified, initial
            conditions will be taken from model.initial_conditions (with initial
            condition parameter values taken from `param_values` if specified).
        """

        if param_values is not None:
            # accept vector of parameter values as an argument
            if len(param_values) != len(self.model.parameters):
                raise Exception("param_values must be the same length as "
                                "model.parameters")
            if not isinstance(param_values, np.ndarray):
                param_values = np.array(param_values)
        else:
            # create parameter vector from the values in the model
            param_values = np.array([p.value for p in self.model.parameters])

        subs = dict((p, param_values[i])
                    for i, p in enumerate(self.model.parameters))
        if y0 is not None:
            # accept vector of initial values as an argument
            if len(y0) != self.y.shape[1]:
                raise Exception("y0 must be the same length as y")
            if not isinstance(y0, np.ndarray):
                y0 = np.array(y0)
        else:
            y0 = np.zeros((self.y.shape[1],))
            for cp, value_obj in self.model.initial_conditions:
                if value_obj in self.model.parameters:
                    pi = self.model.parameters.index(value_obj)
                    value = param_values[pi]
                elif value_obj in self.model.expressions:
                    raise Exception("Can't handle initial conditions from "
                                    "expressions yet") # TODO TODO TODO
                    value = value_obj.expand_expr().evalf(subs=subs)
                else:
                    raise ValueError("Unexpected initial condition value type")
                si = self.model.get_species_index(cp)
                if si is None:
                    raise Exception("Species not found in model: %s" % repr(cp))
                y0[si] = value
                # TODO TODO TODO Take actual derivative
                y0[self.sdot_ix_map[(si, pi)]] = 1.

        # perform the actual integration
        self.integrator.set_initial_value(y0, self.tspan[0])
        # Set parameter vectors for RHS func and Jacobian
        self.integrator.set_f_params(param_values)
        self.y[0] = y0
        i = 1
        while (self.integrator.successful() and
               self.integrator.t < self.tspan[-1]):
            self.y[i] = self.integrator.integrate(self.tspan[i])
            i += 1
        if self.integrator.t < self.tspan[-1]:
            self.y[i:, :] = 'nan'

        for i, obs in enumerate(self.model.observables):
            self.yobs_view[:, i] = \
                (self.y[:, obs.species] * obs.coefficients).sum(1)
        obs_names = self.model.observables.keys()
        obs_dict = dict((k, self.yobs[k]) for k in obs_names)
        for expr in self.model.expressions_dynamic():
            expr_subs = expr.expand_expr().subs(subs)
            func = sympy.lambdify(obs_names, expr_subs, "np")
            self.yexpr[expr.name] = func(**obs_dict)

        # Create the sensitivity matrix as a view on the original ydot array
        # Get only the sensitivity entries and transpose so tspan is last
        # dimension:
        self.ysens = self.y[:, len(model.odes):].T.reshape(
                       (len(model.odes), len(model.parameters), len(self.tspan)))
        # Create the sensitivity matrix of the observables
        for i, obs in enumerate(self.model.observables):
            # We don't need to bother multiplying by the observables coefficients
            # if they are all ones!
            if np.all(obs.coefficients == 1):
                #self.yobs_sens_view[i, :, :] = \
                #                self.ysens[obs.species, :, :].sum(axis=0)
                self.yobs_sens[obs_name] = \
                                self.ysens[obs.species, :, :].sum(axis=0)
            else:
                # FIXME FIXME FIXME
                # There's got to be a better way to do this multiplication
                # (i.e., without a for loop!)
                for p_ix in range(len(model.parameters)):
                    self.yobs_sens[obs.name][p_ix] = \
                      (self.ysens[obs.species, p_ix, :].T *
                                    obs.coefficients).sum(axis=1).T
                    #self.yobs_sens_view[i, p_ix, :] = \
                    #  (self.ysens[obs.species, p_ix, :].T *
                    #                obs.coefficients).sum(axis=1).T

def exp_decay_model():
    Model()
    Monomer('A', ['s'], {'s': ['on', 'off']})
    Rule('A_deg', A(s='on') >> A(s='off'), Parameter('k', 0.01))
    Initial(A(s='on'), Parameter('A_0', 100))
    Observable('A_on', A(s='on'))
    Observable('A_off', A(s='off'))
    Observable('A_tot', A())
    t = np.linspace(0, 500, 100)
    sens = Sensitivity(model, t)
    sens.run()

    # Analytical solutions for A(s='on')
    y = A_0.value * np.exp(-k.value * t)
    dydk = A_0.value * -t * np.exp(-k.value * t)
    dydA0 = np.exp(-k.value * t)
    plt.figure()
    # A_on
    plt.subplot(1, 2, 1)
    plt.plot(t, sens.yobs['A_on'], linewidth=2, label='integr')
    plt.plot(t, y, linewidth=2, color='r', linestyle='--', label='analyt')
    plt.legend(loc='upper right')
    plt.title('A_on')
    # A_off
    plt.subplot(1, 2, 2)
    plt.plot(t, sens.yobs['A_off'], linewidth=2, label='integr')
    plt.legend(loc='lower right')
    plt.title('A_off')

    # Sensitivity Plots
    plt.figure()
    # dA_on/dk
    plt.subplot(2, 2, 1)
    plt.plot(t, sens.ysens[0, 0], linewidth=2, label='integr')
    plt.plot(t, dydk, linewidth=2, color='r', linestyle='--', label='analyt')
    plt.legend(loc='lower right')
    plt.title('Sens of A_on to k')
    # dA_on/dA0
    plt.subplot(2, 2, 2)
    plt.plot(t, sens.ysens[0, 1], linewidth=2, label='integr')
    plt.plot(t, dydA0, linewidth=2, color='r', linestyle='--', label='analyt')
    plt.legend(loc='upper right')
    plt.title('Sens of A_on to A0')
    # dA_off/dk
    plt.subplot(2, 2, 3)
    plt.plot(t, sens.ysens[1, 0], linewidth=2, label='integr')
    plt.legend(loc='upper right')
    plt.title('Sens of A_off to k')
    # dA_on/dA0
    plt.subplot(2, 2, 4)
    plt.plot(t, sens.ysens[1, 1], linewidth=2, label='integr')
    plt.legend(loc='lower right')
    plt.title('Sens of A_off to A0')
    return sens

def mrna_protein_model():
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
    sens = Sensitivity(model, t)
    sens.run()

    """
    plt.figure()
    plt.subplot(1, 2, 1)
    plt.plot(t, sens.yobs['m_'], label='mRNA', linewidth=2)
    plt.title('mRNA')
    plt.subplot(1, 2, 2)
    plt.plot(t, sol.yobs['p_'], label='Protein', linewidth=2)
    plt.title('mRNA')
    """

    plt.figure()
    subplot_ix = 1
    for obs_ix, obs_name in enumerate(['m_', 'p_']):
        for p_ix in range(4):
            plt.subplot(2, 4, subplot_ix)
            plt.plot(t, sens.yobs_sens[obs_name][p_ix])
            plt.title('%s, p%d' % (obs_name, p_ix))
            subplot_ix += 1

if __name__ == '__main__':

    plt.ion()
    sens = exp_decay_model()
    #sens = mrna_protein_model()

