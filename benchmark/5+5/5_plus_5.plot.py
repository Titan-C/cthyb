#!/bin/env pytriqs

from pytriqs.archive import *
from pytriqs.gf.local import *
from pytriqs.plot.mpl_interface import *
from matplotlib.backends.backend_pdf import PdfPages
from pytriqs.applications.impurity_solvers.cthyb import change_mesh

from multiorbital import *
from params import *

# Read the reference table file
tau = []
data = []
for line in open(ref_file_name,'r'):
    cols = line.split()
    tau.append(float(cols[0]))
    data.append([-float(c) for c in cols[1:]])

g_ref = GfImTime(indices = range(len(data[0])), beta=tau[-1], n_points=len(tau), kind='full_bins')
for nt, d in enumerate(data):
    for nc, val in enumerate(d):
        g_ref.data[nt,nc,nc] = val

# Calculate theoretical curves
if not use_interaction:
    g_theor = GfImTime(indices = range(len(cubic_names)), beta=beta, n_points=n_tau)
    for nc, cn in enumerate(cubic_names):
        V = delta_params[cn]['V']
        e = delta_params[cn]['e']

        e1 = e - V
        e2 = e + V

        g_theor_w = GfImFreq(indices = [0], beta=beta)
        g_theor_w <<= 0.5*inverse(iOmega_n - e1) + 0.5*inverse(iOmega_n - e2)
        g_theor[nc,nc] <<= InverseFourier(g_theor_w)

# Read the results
arch = HDFArchive(results_file_name,'r')

pp = PdfPages('G.pdf')

for sn in spin_names:
    for nc, cn in enumerate(cubic_names):
        plt.clf()

        gf = change_mesh(arch['G_tau'][mkind(sn,cn)[0]],200)

        # Plot the results
        oplot(gf, name="cthyb")
        # Plot the reference curve
        if use_interaction:
            oplot(g_ref[nc,nc], name="ED")
        else:
            oplot(g_theor[nc,nc], name="ED")

        axes = plt.gca()
        axes.set_title('$' + sn + '$, $' + cn + '$')
        axes.set_xlabel('$\\tau$')
        axes.set_ylabel('$G(\\tau)$')
        axes.set_ylim(-1,0)
        axes.legend(loc='lower center',prop={'size':14})

        pp.savefig(plt.gcf())

pp.close()
