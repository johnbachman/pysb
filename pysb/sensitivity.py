from pysb import *
from pysb.integrate import Solver
import numpy as np
from matplotlib import pyplot as plt

class Sensitivity(object):
    def __init__(self, model):
        pass

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
    sol = Solver(model, t)
    sol.run()

    plt.ion()
    plt.figure()
    plt.plot(t, sol.yobs['m_'], label='mRNA')
    plt.plot(t, sol.yobs['p_'], label='Protein')
    plt.legend(loc='lower right')
