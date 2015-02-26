#!/bin/env pytriqs

import pytriqs.utility.mpi as mpi
from pytriqs.gf.local import *
from pytriqs.operators import *
from pytriqs.archive import HDFArchive
from pytriqs.applications.impurity_solvers.cthyb import *

def print_master(msg):
    if mpi.rank==0: print msg

print_master("Welcome to asymm_bath test (1 band with a small asymmetric hybridization function).")
print_master("This test helps to detect sampling problems.")

# H_loc parameters
beta = 40.0
ed = -1.0
U = 2.0

use_prop = True

#epsilon = [0.0,0.1,0.2,0.3,0.4,0.5]
epsilon = [0.4]
V = 0.2

# Parameters
n_iw = 1025
n_tau = 10001

p = {}
p["max_time"] = -1
p["random_name"] = ""
p["random_seed"] = 123 * mpi.rank + 567
p["length_cycle"] = 50
p["n_warmup_cycles"] = 20000
p["n_cycles"] = 1000000
p["make_histograms"] = True
p["use_proposed"] = use_prop

# Block structure of GF
gf_struct = {'up':[0], 'dn' : [0]}

# Hamiltonian
H = U*n("up",0)*n("dn",0)

# Histograms to be saved
histos_to_save = [('histo_opcount0',int),
                  ('histo_opcount1',int),
                  ('histo_opcount_total',int),
                  ('histo_insert_length_proposed',float),
                  ('histo_insert_length_accepted',float),
                  ('histo_remove_length_proposed',float),
                  ('histo_remove_length_accepted',float),
                  ('histo_shift_length_proposed',float),
                  ('histo_shift_length_accepted',float)]

# Quantum numbers
qn = []
qn.append(n("up",0))
qn.append(n("dn",0))
p["partition_method"] = "quantum_numbers"
p["quantum_numbers"] = qn

# Construct solver
S = SolverCore(beta=beta, gf_struct=gf_struct, n_tau=n_tau, n_iw=n_iw)

def read_histo(f,type_of_col_1):
    histo = []
    for line in f:
        cols = filter(lambda s: s, line.split(' '))
        histo.append((type_of_col_1(cols[0]),float(cols[1]),float(cols[2])))
    return histo

if use_prop:
    filename = "asymm_bath_with.h5"
else:
    filename = "asymm_bath_without.h5"

if mpi.rank==0:
    arch = HDFArchive(filename,'w')

# Set hybridization function
for e in epsilon:
    delta_w = GfImFreq(indices = [0], beta=beta)
    delta_w << (V**2) * inverse(iOmega_n - e)

    S.G0_iw["up"] << inverse(iOmega_n - ed - delta_w)
    S.G0_iw["dn"] << inverse(iOmega_n - ed - delta_w)

    S.solve(h_loc=H, **p)

    if mpi.rank==0:
        d = {'G_tau':S.G_tau, 'beta':beta, 'U':U, 'ed':ed, 'V':V, 'e':e}

        for histo_name, type_of_col_1 in histos_to_save:
            d[histo_name] = read_histo(open(histo_name+'.dat','r'),type_of_col_1)

        #arch['epsilon_' + str(e)] = d
        arch['G'] = S.G_tau
