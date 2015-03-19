from cthyb import SolverCore
from pytriqs.gf.local import *
import pytriqs.utility.mpi as mpi

class Solver(SolverCore):

    def __init__(self, beta, gf_struct, n_iw=1025, n_tau=10001, n_l=30):
        """
        Initialise the solver.

        Parameters
        ----------
        beta : scalar
               Inverse temperature.
        gf_struct : dict{str:list}
                    Structure of the Green's functions. It must be a
                    dictionary which maps the name of each block of the
                    Green's function as a string to a list of integer
                    indices.
                    For example: { 'up': [1,2,3], 'down', [1,2,3] }.
        n_iw : integer, optional
               Number of Matsubara frequencies used for the Green's functions.
        n_tau: integer, optional
               Number of imaginary time points used for the Green's functions.
        n_l: integer, optional
             Number of legendre polynomials to use in accumulations of the Green's functions.
        
        """        
        # Initialise the core solver
        SolverCore.__init__(self, beta, gf_struct, n_iw=n_iw, n_tau=n_tau, n_l=n_l)

        self.Sigma_iw = self.G0_iw.copy()
        self.Sigma_iw.zero()
        self.G_iw = self.G0_iw.copy()
        self.G_iw.zero()
        self.gf_struct = gf_struct
        self.n_iw = n_iw
        self.n_tau = n_tau

    def solve(self, h_loc, **params_kw):
        """
        Solve the impurity problem.
        If measure_g_tau (default = True), G_iw and Sigma_iw will be calculated and their tails fitted. 
        In addition to the solver parameters, parameters to control the tail fitting can be provided.

        Parameters
        ----------
        h_loc : Operator object
                The local Hamiltonian of the impurity problem to be solved.
        fit_n_moments : integer, optional, default = 3
                        Number of moments to fit in the tail of Sigma_iw.
        fit_known_moments : dict{str:TailGf object}, optional, default = {block_name: TailGf(dim1, dim2, n_moments, order_min)}
                            Known moments of Sigma_iw, given as a TailGf object.
        fit_min_n : integer, optional, default = int(0.8 * self.n_iw)
                    Index of iw from which to start fitting.
        fit_max_n : integer, optional, default = n_iw
                    Index of iw to fit until.

        """

        # If tail parameters provided for Sigma_iw fitting, use them, otherwise use defaults
        if not (("fit_min_n" in params_kw) or ("fit_max_n" in params_kw)):
	    if mpi.is_master_node(): 
                warning = ("!------------------------------------------------------------------------------------!\n"
                           "! WARNING: Using default high-frequency tail fitting parameters in the CTHYB solver. !\n"
                           "! You should check that the fitting range is suitable for your calculation!          !\n"
                           "!------------------------------------------------------------------------------------!")
                print warning

        fit_n_moments = params_kw.pop("fit_n_moments", 3)
        if "fit_known_moments" in params_kw:
            fit_known_moments = params_kw.pop("fit_known_moments")
        else:
            fit_known_moments = {}
            for name, indices in self.gf_struct.items():
                dim = len(indices)
                fit_known_moments[name] = TailGf(dim,dim,0,0) # TailGf(dim1, dim2, n_moments, order_min)
        fit_min_n = params_kw.pop("fit_min_n", int(0.8 * self.n_iw)) # Fit last 20% of frequencies
        fit_max_n = params_kw.pop("fit_max_n", self.n_iw)

        # Call the core solver's solve routine
        SolverCore.solve(self, h_loc = h_loc, **params_kw)

        # Post-processing: 
        # (only supported for G_tau, to permit compatibility with dft_tools)
        if self.last_solve_parameters["measure_g_tau"] == True:
            # Fourier transform G_tau to obtain G_iw
            for name, g in self.G_tau:
                self.G_iw[name] << Fourier(g)

            # Solve Dyson's eq to obtain Sigma_iw and fit the tail
            self.Sigma_iw = inverse(self.G0_iw) - inverse(self.G_iw)
            for name, sig in self.Sigma_iw:
                 self.Sigma_iw[name].fit_tail(fit_known_moments[name], fit_n_moments, fit_min_n, fit_max_n)
                 # Now fix the tail of G_iw
                 self.G_iw[name].tail = inverse( inverse(self.G0_iw[name].tail) - self.Sigma_iw[name].tail )
