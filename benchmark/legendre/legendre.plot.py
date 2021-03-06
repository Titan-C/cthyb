#!/bin/env pytriqs

from pytriqs.archive import *
from pytriqs.gf.local import *
from pytriqs.plot.mpl_interface import plt

arch = HDFArchive('legendre.h5','r')
ed_arch = HDFArchive('legendre.ed.h5','r')

for spin in ("up","dn"):
    plt.plot(arch['G_l'][spin].data.flatten(),label="cthyb," + {'up':"$\uparrow\uparrow$",'dn':"$\downarrow\downarrow$"}[spin])
    plt.plot(ed_arch[spin].data.flatten(),label="ED," + {'up':"$\uparrow\uparrow$",'dn':"$\downarrow\downarrow$"}[spin])

axes = plt.gca()
axes.set_xlabel('$l$')
axes.set_ylabel('$G_l$')
axes.set_ylim((-2.0,1.0))
axes.legend(loc='lower center',prop={'size':10})

plt.savefig('G_l.pdf')