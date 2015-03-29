# Generated automatically using the command :
# wrapper_desc_generator.py ../c++/solver_core.hpp -p -mpytriqs.applications.impurity_solvers.cthyb -o cthyb --moduledoc "The cthyb solver"
from wrap_generator import *

# The module
module = module_(full_name = "pytriqs.applications.impurity_solvers.cthyb", doc = "The cthyb solver")

# All the triqs C++/Python modules
module.use_module('gf')
module.use_module('operators')

# Add here all includes beyond what is automatically included by the triqs modules
module.add_include("../c++/solver_core.hpp")

# Add here anything to add in the C++ code at the start, e.g. namespace using
module.add_preamble("""
#include <triqs/python_tools/converters/pair.hpp>
#include <triqs/python_tools/converters/map.hpp>
#include <triqs/python_tools/converters/vector.hpp>
using namespace triqs::gfs;
using triqs::utility::many_body_operator;
using namespace cthyb;
#include "./converters.hxx"
""")

# The class solver_core
c = class_(
        py_type = "SolverCore",  # name of the python class
        c_type = "solver_core",   # name of the C++ class
)

c.add_constructor("""(double beta, std::map<std::string,indices_type> gf_struct, int n_iw = 1025, int n_tau = 10001, int n_l = 50)""", 
                  doc = """ """)

c.add_method("""void solve (**cthyb::solve_parameters_t)""", 
             doc = """  Parameter Name      Type                                   Default                                        Documentation

  h_loc               real_operator_t                        --                                             Atomic Hamiltonian
  n_cycles            int                                    --                                             Number of QMC cycles
  partition_method    std::string                            "autopartition"                                Partition method
  quantum_numbers     std::vector<real_operator_t>           std::vector<real_operator_t>{}                 Quantum numbers
  length_cycle        int                                    50                                             Length of a single QMC cycle
  n_warmup_cycles     int                                    5000                                           Number of cycles for thermalization
  random_seed         int                                    34788+928374*boost::mpi::communicator().rank() Seed for random number generator
  random_name         std::string                            ""                                             Name of random number generator
  max_time            int                                    -1                                             Maximum runtime in seconds, use -1 to set infinite
  verbosity           int                                    ((boost::mpi::communicator().rank()==0)?3:0)   Verbosity level
  move_shift          bool                                   true                                           Add move_shift as a move?
  use_trace_estimator bool                                   false                                          Calculate the full trace or use an estimate?
  measure_g_tau       bool                                   true                                           Measure G(tau)?
  measure_g_l         bool                                   false                                          Measure G_l (Legendre)?
  measure_pert_order  bool                                   false                                          Measure perturbation order?
  make_histograms     bool                                   false                                          Make histograms of the trace computation?
  proposal_prob       std::map<std::string, double>          (std::map<std::string,double>{})               Operator insertion/removal probabilities for different blocks
  static_observables  std::map<std::string, real_operator_t> (std::map<std::string,real_operator_t>{})      List of static observables to be measured (with their names)  """)

c.add_property(name = "last_solve_parameters",
               getter = cfunction("cthyb::solve_parameters_t get_last_solve_parameters ()"),
               doc = """Set of parameters used in the last call to solve """)

c.add_property(name = "G0_iw",
               getter = cfunction("block_gf_view<imfreq> G0_iw ()"),
               doc = """G0(iw) in imaginary frequencies """)

c.add_property(name = "Delta_tau",
               getter = cfunction("block_gf_view<imtime> Delta_tau ()"),
               doc = """Delta(tau) in imaginary time """)

c.add_property(name = "G_tau",
               getter = cfunction("block_gf_view<imtime> G_tau ()"),
               doc = """G(tau) in imaginary time """)

c.add_property(name = "G_l",
               getter = cfunction("block_gf_view<legendre> G_l ()"),
               doc = """G_l in Legendre polynomials representation """)

c.add_property(name = "static_observables",
               getter = cfunction("std::map<std::string,double> static_observables ()"),
               doc = """Accumulated static observable with a given name """)

c.add_property(name = "atomic_gf",
               getter = cfunction("block_gf_view<imtime> atomic_gf ()"),
               doc = """Atomic G(tau) in imaginary time """)

c.add_property(name = "average_sign",
               getter = cfunction("mc_sign_type average_sign ()"),
               doc = """Monte Carlo average sign """)

c.add_property(name = "eigensystems",
               getter = cfunction("std::vector<std::pair<vector<double>,matrix<double>>> get_eigensystems ()"),
               doc = """Eigensystems of the atomic problem
  Returns a list of pairs (E,U), where H = U * diag(E) * U^+ (for each subspace) """)

module.add_class(c)

module.generate_code()