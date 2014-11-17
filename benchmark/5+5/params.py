#!/bin/env pytriqs

import pytriqs.utility.mpi as mpi

results_file_name = "5_plus_5.h5"
matrix_results_file_name = "5_plus_5.matrix.h5"
ref_file_name = "5_plus_5.ref.dat"

# Block structure of GF
L = 2 # d-orbital
cubic_names = ("xy","yz","z^2","xz","x^2-y^2")
mkind = lambda sn, cn: (sn+'_'+cn, 0)

# Input parameters
beta = 40.
mu = 26

use_interaction = True
U = 4.0
J = 0.7
F0 = U
# Sasha uses coefficient 0.63 instead of 0.625
F2 = J*(14.0/(1.0 + 0.63))
F4 = F2*0.63

# Dump the local Hamiltonian to a text file (set to None to disable dumping)
H_dump = "H.txt"
# Dump Delta parameters to a text file (set to None to disable dumping)
Delta_dump = "Delta_params.txt"

# Hybridization function parameters
# Delta(\tau) is diagonal in the basis of cubic harmonics
# Each component of Delta(\tau) is represented as a list of single-particle
# terms parametrized by pairs (V_k,\epsilon_k).
delta_params={"xy" : {'V':0.5,'e':-0.2},
              "yz" : {'V':0.5,'e':-0.15},
              "z^2" : {'V':0.5,'e':-0.1},
              "xz" : {'V':0.5,'e':0.05},
              "x^2-y^2" : {'V':0.5,'e':0.4}}

atomic_levels = {('up_xy',0) :        -0.2,
                 ('dn_xy',0) :        -0.2,
                 ('up_yz',0) :        -0.15,
                 ('dn_yz',0) :        -0.15,
                 ('up_z^2',0) :       -0.1,
                 ('dn_z^2',0) :       -0.1,
                 ('up_xz',0) :        0.05,
                 ('dn_xz',0) :        0.05,
                 ('up_x^2-y^2',0) :   0.4,
                 ('dn_x^2-y^2',0) :   0.4}

n_iw = 1025
n_tau = 10001

p = {}
p["max_time"] = -1
p["random_name"] = ""
p["random_seed"] = 123 * 83 + 567
p["verbosity"] = 3
p["length_cycle"] = 50
p["n_warmup_cycles"] = 50000
p["n_cycles"] = 7500000
p["partition_method"] = "autopartition"
p["measure_g_tau"] = False
p["move_shift"] = False
p["measure_pert_order"] = False
p["make_histograms"] = False
p["use_trace_estimator"] = False
