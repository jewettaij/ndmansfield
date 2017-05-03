#ifndef _NDMANSFIELD_H   //This check insures that C-preprocessor only includes
#define _NDMANSFIELD_H   //the following text once (not multiple times)





#include<iostream>
#include<fstream>
#include<string>
#include<sstream>
#include<vector>

//#define NDEBUG         //<--dissables assert()
#include<cassert>      // This defines the "assert()" function.
                       // You can ask me what this does if you want.
                       // Strictly speaking you don't need to use assert()s
                       // in the code, but they are helpful for debugging.
#include "err_report.h"

#include "globals.h"




class Hamiltonian; 



class NDmansfield
{
  long *aSize;          //How big is the box that the polymer lives in?

  long *aIseqFromIloc;  //A 1-dimensional array which stores the sequential
                        //location of the monomer (from a polymer) occupying 
                        //a particular site on the lattice, which is determined
                        //from the index into this array.
                        //
                        //For example, in 3-dimensions, the monomer located
                        //at position x,y,z is (probably) stored at:
                        //
                        // aIseqFromIloc[ aSize[1]*aSize[0]*z + aSize[0]*y + x ]
                        //
                        //(where x,y,z are integers, and
                        //   x ranges from 0 ... aSize[0]-1
                        //   y ranges from 0 ... aSize[1]-1
                        //   z ranges from 0 ... aSize[2]-1)
                        //
                        // This is just an example.  The actual mapping from
                        // memory-space to spatial coordinates is implementation
                        // dependent.  To keep things abstract, you should
                        // navigate this array using these functions:
                        //   CoordsFromIloc()  and  IlocFromCoords()

  long *aIlocFromIseq; //An inverse lookup-table for aIseqFromIloc

  long num_cells;    // Total number of cells in the lattice
                     // = aSize[0] * aSize[1] * aSize[2] ...
  long poly_length;  // Length of the polymer.
                     // (poly_length <= num_cells, because some of the sites
                     //  are unavailable, such as the sites on the boundaries.)

  // Reverse the order of the chain in the interval
  // from [iseqA, iseqB], inclusive
  void ReverseInc(long iseqA, long iseqB);

  friend istream& operator >> (istream& in,  NDmansfield& ndmansfield);
  friend ostream& operator << (ostream& out, NDmansfield const& ndmansfield);
  friend Hamiltonian;  // The "Hamiltonian" object calculates energy
                       // used for Boltzmann weighting during MonteCarlo.

  void Clear(); // fill the contents of aIseqFromIloc[] and aIlocFromIseq[]
                // with blanks (which are asumed to be preallocated)

  // Fill the available lattice with a simple default polymer conformation
  void DefaultPolyShape();

  // Some frequently accessed arrays (of size g_dim) used by MonteCarloStep()
  long *aEndLoc;
  long *aNeighLoc;

public:

  inline long NumCells() const { return num_cells; }
  inline long PolyLength() const { return poly_length; }

  // Return the lattice size in direction "d":
  inline long GetSize(int d) const {
    return aSize[d];
  }



  NDmansfield(const long *set_size);
  ~NDmansfield();


  static const long NOT_FOUND = -1; // Return this when asking for contents in
                                     // out-of-bounds or unoccupied cells


  // IlocFromCoords(aCoords)
  // Returns the index position into the array aIseqFromIloc[iloc] 
  // corresponding to the explicit x,y,z coordinates stored in the 
  // aCoords[] array.

  // In this version, I do not do any bounds-checking.
  // I assume the polymer lies at least 1 lattice site away from the boundary.
  // (That way, we can safely check the nearest-neighbor sites without worry.)

  inline long 
  IlocFromCoords(long const *aCoords) const
  {
    long i = aCoords[g_dim-1];
    for(int d=g_dim-2; d >= 0; --d) {
      i *= aSize[d];
      i += aCoords[d];
    }

    // ALTERNATE METHOD:
    // long i;
    // int d = g_dim-1;
    // i = aCoords[d];
    // --d;
    // while (d >= 0) {
    //   i *= aSize[d];
    //   i += aCoords[d];
    //   --d;
    // }

    assert((0 <= i) && (i < num_cells));
    return i;
  }


  // The next version deals with periodic boundary conditions (untested)
  inline long 
  IlocFromCoordsPeriodic(long const *aCoords) const
  {
    long i = aCoords[g_dim-1];
    // loop coordinates when out of bounds
    if (i > aSize[g_dim-1])
      i -= aSize[g_dim-1];
    else if (i < 0)
      i += aSize[g_dim-1];

    for(int d=g_dim-2; d >= 0; --d) {
      i *= aSize[d];
      long iincr = aCoords[d];

      // loop coordinates when out of bounds
      if (iincr > aSize[d])
        iincr -= aSize[d];
      else if (i < 0)
        iincr += aSize[d];

      i += iincr;
    }

    // ALTERNATE METHOD:
    // long i;
    // int d = g_dim-1;
    // i = aCoords[d];
    // --d;
    // while (d >= 0) {
    //   i *= aSize[d];
    //   i += aCoords[d];
    //   --d;
    // }

    assert((0 <= i) && (i < num_cells));
    return i;
  }



  // CoordsFromIloc(aloc, iloc)
  // converts the index "iloc" corresponding to the place in the 1-D table where
  // a particular lattice is stored, into an array of coordinates (eg x,y,z)

  inline void CoordsFromIloc(long *aCoords, long iloc) const
  {
    assert((0 <= iloc) && (iloc < num_cells));
    long next_i;
    long size;
    for(int d=0; d < g_dim-1; ++d)
    {
      size = aSize[d];
      next_i = iloc / size;
      //aCoords[d] = iloc % size;  
      aCoords[d] = iloc-size*next_i; //this probably computes the remainder faster
      assert((0 <= aCoords[d]) && (aCoords[d] < aSize[d]));
      iloc = next_i;
    }
    aCoords[g_dim-1] = iloc;
    assert((0 <= aCoords[g_dim-1]) && (aCoords[g_dim-1] < aSize[g_dim-1]));
  }


  //ALTERNATE VERSION
  //inline void CoordsFromIloc(long *aCoords, long iloc) const
  //{
  //  long num_cells_in_dimension_d;
  //  long coordinate;
  //  for(int d=g_dim-1; d > 0; --d) {
  //    num_cells_in_dimension_d = aNumCellsInDimension[d];
  //    coordinate = iloc / num_cells_in_dimension_d;
  //    iloc -= coordinate * num_cells_in_dimension_d;  //iloc is the remainder
  //    aCoords[d] = coordinate;
  //    assert((0 <= aCoords[d]) && (aCoords[d] < aSize[d]));
  //  }
  //  aCoords[0] = iloc;
  //  assert((0 <= aCoords[0]) && (aCoords[0] < aSize[0]));
  //}


  inline void CoordsFromIseq(long *aCoords, long iseq) const {
    assert((0 <= iseq) && (iseq < poly_length));
    CoordsFromIloc(aCoords, aIlocFromIseq[iseq]);
  }

  inline long IseqFromCoords(long const *aCoords) const {
    long iloc = IlocFromCoords(aCoords);
    return aIseqFromIloc[iloc];
  }


  // The following function executes one iteration of Monte Carlo on the lattice
  // It returns whether or not the move was accepted (true or false)
  bool
  MonteCarloStep(Hamiltonian& hamiltonian, //calculates energies
                 double& total_energy); // keep track of total energy here


  bool IsCyclic(); //This function returns true of the end of the path is
                   //on a lattice site which is adjacent to the beginning.

}; //struct NDmansfield





//Write the current state of the system to a file.
istream& operator >> (istream& in,  NDmansfield& ndmansfield);

//Read in the current state of the system from a file.
ostream& operator << (ostream& out, NDmansfield const& ndmansfield);






#endif //#ifndef _NDMANSFIELD_H
