This is a hastilly written note to explain a hack I use
to generate images of snapshots from ndmansfield trajectories.
This uses the tools I am familiar with (VMD, topotools, moltemplate).

(Comment/Recommendation
 When you display the polymer in VMD, you can use the internal-tachyon 
 renderer with ambient occlusion: "Display"->"Display Settings"->"Amb. Occl.")

I think this proceedure works, but some of these steps may 
be wrong.  Even so, hopefully it gives you an idea how to 
proceed.  I'm sure there are better ways of doing this.

--- Instructions ----

1) If you have not already, install VMD.
   Also install moltemplate (www.moltemplate.org)

2) 
Use the "make_LAMMPS_data_file.py" script to create a file
(usually named "system.data") this way:
 EXAMPLE: For a polymer trapped in a box with 16 x 16 x 16 lattice sites:

N=$((16*16*16))

./make_LAMMPS_data_file.py $N $((N-1)) > system.data.tmp

# For a circular chain use:  ./make_LAMMPS_data_file.py $N $N > system.data.tmp

3) Create a file containing the coordinates you want:

# Example assuming "traj_ndmansfield.raw" was created using ndmansfield:
TIME=100 # <-- which snapshot do you want to extract from your trajectory file?
./select_structure $TIME $N  < traj_ndmansfield.raw > coords.raw

#(coords.raw is a 3-column ASCII txt file containing 3 coordinates on each line)

4) Now copy the coordinates from coords.raw into the data file.
   One way is to use the "raw2data.py" utility that comes with moltemplate:

raw2data.py system.data.tmp < coords.raw > system.data


5) Finally, load this file "system.data" into VMD using the instructions here:

README_using_VMD_to_view_a_LAMMPS_data_file.txt

