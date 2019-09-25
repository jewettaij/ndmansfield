ndmansfield
===========

Usage:
```
ndmansfield -box xsize ysize zsize -tsave tsave [options] > traj.raw
```
Example:
```
ndmansfield -box 10 10 10 -tsave 2000 -tstop 20000 | tail -n 1001 > traj.raw
```

The "ndmansfield" program generates a file containing the coordinates of a
random path(s) which visits every site in a lattice inside a rectangular box.
(These paths have been used to model the conformations of coarse-grained
biopolymers, for example.)

The "ndmansfield" program is an implementation of
"Unbiased sampling of lattice Hamilton path ensembles"
Marc L. Mansfield, J. Chem Phys, 2006
The
[move-set used in each Monte-Carlo iteration](images/Mansfield_monte-carlo_move_JCP2006_Fig1.png)
is the same as the
["BACKBITE" move from the Mansfield J. Chem Phys 1982 paper (figure 1c).](images/Mansfield_monte-carlo_move_JCP1982_Fig1.png)

Each step in the path corresponds to a location in the lattice, and the
the coordinates of each step in the path is saved as 3 integers on a separate
line of the file. (...Assuming your lattice is 3 dimensional.  See below.)
This file is written to the standard output.

This program uses a Monte-Carlo procedure to to generate many such paths.
Each of these paths is saved to the output file.  Different paths are delimited
by blank lines in the output file.

This code works in an arbitrary number of dimensions, and with arbitrary box
sizes. All of this code successfully compiles using g++ 4.8.4 (tested 2016-1-14)


# Usage Syntax:
```
ndmansfield -box xsize ysize zsize -tsave tsave [options] > trajectory.raw
```

# Required Arguments

## -box xsize ysize zsize

How big is the box that the path lives in?
(xsize ysize and zsize are positive integers)

*(Note: The number of numeric arguments following
the "-box" command determines the number of
dimensions the polymer will live in.  For a 2-D
simulation, you would omit the "zsize" argument,
for example.)*

## -tsave tsave

After a certain number of iterations, save the
coordinates of the path in a file using
.RAW format (see below).
### Cyclic paths:
If cyclic paths are requested (see below), then
the simulation will iterate until a cyclic path
conformation discovered.  In that case, cyclic
paths which were discovered less than tsave
iterations ago are discarded


# Optional Arguments:

## -tstop tstop

Specify how many iterations of Monte-Carlo to run
before halting the simulation.  If specified, then
the simulation is halted when the iteration count
reaches tstop (an integer).
(By default tstop=infinity)


## -tstart tstart

Start running the simulation with the iteration
count set to tstart. This is only potentially
useful if you are continuing a simulation you
stopped earlier and you want to keep track of
how many total iterations have ellapsed.
(By default, tstart=0.)


## -startcrd init_crd.raw

Supply a .RAW file xyz(...) coordinate file
with the starting coordinates of your path.
Each line in the file should contains a list of
integers corresponding to the coordinates in
step in the path.
The number of lines in the file should match
the number of steps in the path, which should
match the number of lattice sites in the box.
By default, the initial path is a simple "zig-zag"
shape which fills the box.
Each line in the file should have g_dim integers
(...where "g_dim" is the number of dimensions
of the lattice and is usually 3.)


## -seed n

Select the seed for the random number generator.


## -cyclic yes/no

Limit results to cyclic paths?
*(Warning: Cyclic paths which do not cross periodic
boundaries are only possible in boxes
which contain an even number of
lattice sites.)*


## -cyclic-direction d

Limit yourself to paths whose starting and ending
points lie on opposite sides of the box along
direction d.  For example, if d=0, then the path's
starting and ending points must lie on opposite
sides of the box along the x direction. IE they
have the lowest and highest possible x coordinate,
and the remaining coordinates must be identical.
If d=1, then the two ends must lie in opposite
sides of the box along the y direction
d=2 ... z direction, etc...



# Optional Arguments that effect physical properties of the chain:

(By default all Monte-Carlo moves will be accepted.
If you start adding energetic bias to some shapes over others,
then Monte-Carlo moves that increase chain energy are more likely to
be rejected. All energies are in units of reduced temperature, kB*T)

## -bend-energy E

If specified, then 90-degree bends in the chain
will increase the energy of the chain by E.

## -twist-energy E

If specified, then right-handed torsions in
the chain backbone will increase the energy
of the chain by E.  Similarly, left-handed
torsions in the chain will reduce the energy by E.
Formally, a right-handed torsion occurs any time
the vector triple product ((v1 x v2) . v3) between
three successive bonds in the chain backbone
(vectors v1, v2, and v3) is positive (+1).  
A left-handed torsion occurs whenever the
tripple-product ((v1 x v2) . v3) is negative (-1).
(Possible values for the triple product are -1,0,1)

*(Warning: This command line argument only makes
sense for simulations in 3-dimensions.)*
