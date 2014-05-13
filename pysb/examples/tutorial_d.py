from pysb import *

Model()
Monomer('Raf', ['s', 'k'], {'s': ['u', 'p']})
Monomer('MEK', ['s218', 's222', 'k'], {'s218': ['u', 'p'], 's222': ['u', 'p']})
Parameter('kf', 1e-5)
Parameter('Raf_0', 7e4)
Parameter('MEK_0', 3e6)

Rule('Raf_binds_MEK_at_s218',
     Raf(k=None) + MEK(s218=None) >> Raf(k=1) % MEK(s218=1),
     kf)
