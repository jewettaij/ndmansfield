#!/usr/bin/env python

g_program_name = __file__.split('/')[-1]
g_version_str  = '0.1'
g_date_str     = '2016-3-27'

import sys


class InputError(Exception):
    """ A generic exception object containing a string for error reporting.

    """

    def __init__(self, err_msg):
        self.err_msg = err_msg
    def __str__(self):
        return self.err_msg
    def __repr__(self):
        return str(self)



n_atoms = int(sys.argv[1])
n_bonds = n_atoms-1
Lx = n_atoms
Ly = n_atoms
Lz = n_atoms
if len(sys.argv) > 2:
    n_bonds = int(sys.argv[2])
if len(sys.argv) > 5:
    Lx = float(sys.argv[3])
    Ly = float(sys.argv[4])
    Lz = float(sys.argv[5])
n_atom_types = 1
n_bond_types = 1

sys.stdout.write("\n" +
                 str(n_atoms) + " atoms\n" +
                 str(n_bonds) + " bonds\n" +
                 str(n_atom_types) + " atom types\n" +
                 str(n_bond_types) + " bond types\n" +
                 str(-0.5*Lx) + " " + str(0.5*Lx) + " xlo xhi\n" +
                 str(-0.5*Ly) + " " + str(0.5*Ly) + " ylo yhi\n" +
                 str(-0.5*Lz) + " " + str(0.5*Lz) + " zlo zhi\n" +
                 "\n")

sys.stdout.write("\n"
                 "Atoms\n"
                 "\n")

for i in range(0, n_atoms):
    x = y = z = i
    a_id = i+1
    m_id = 1
    a_type = 1
    charge = 0.0
    sys.stdout.write(str(a_id) + " " +
                     str(m_id) + " " +
                     str(a_type) + " " +
                     str(charge) + " " +
                     str(x) + " " + str(y) + " " + str(z) + "\n")

sys.stdout.write("\n"
                 "Bonds\n"
                 "\n")

for i in range(0, n_bonds):
    b_id = i+1
    b_type = 1
    a1_id = i+1
    a2_id = i+2
    if a2_id > n_atoms:
        a2_id -= n_atoms
    sys.stdout.write(str(b_id) + " " +
                     str(b_type) + " " +
                     str(a1_id) + " " + str(a2_id) + "\n")

