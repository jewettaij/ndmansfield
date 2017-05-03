#include <iostream>
#include <fstream>
#include <string>
#include <sstream>     // defines stringstream
#include <cstdlib>     // defines atod()
#include <cassert>
using namespace std;

#include "random_gen.h"

#include "ndmansfield.h"   // Contains a description of the NDmansfield data structure
                       // we will be using to store the state of our system.

#include "hamiltonian.h" // The object that will calculate the energy.


// NDmansfield::MonteCarloStep() 
//    Perform a single iteration in the simulation.


bool NDmansfield::
MonteCarloStep(Hamiltonian& hamiltonian,
               double& energy)
{
  // Choose either end-point of the chain
  long iseq_end = RANDOM_INT(2)*(poly_length-1); //--> 0 or poly_length-1
  long iseq_end_predecessor = 1;
  if (iseq_end == poly_length-1)
    iseq_end_predecessor = poly_length-2;  

  // Find the spatial location of this end-point
  CoordsFromIseq(aEndLoc, iseq_end);

  // Pick one of the neighbors of this endpoint at random.
  // Choose from one of the 2*g_dim neighboring cells
  for (long d=0; d < g_dim; d++)
    aNeighLoc[d] = aEndLoc[d];
  long direction = RANDOM_INT(2*g_dim) - g_dim;
  if (direction >= 0) 
    aNeighLoc[direction] += 1;
  else
    aNeighLoc[-1-direction] -= 1;

  long iseq_neigh = IseqFromCoords(aNeighLoc);
  if (iseq_neigh == NOT_FOUND)
    return true; //do nothing. move does not modify the chain (move is accepted)
  if (iseq_neigh == iseq_end_predecessor)
    return true; //do nothing. move does not modify the chain (move is accepted)

  // Reverse the order of the chain between iseq_end and iseq_neigh-1
  // (or between iseq_neigh+1 and iseq_end)
  long iseqA, iseqB;
  if (iseq_end < iseq_neigh) {
    iseqA = iseq_end;
    iseqB = iseq_neigh-1;
  }
  else {
    iseqA = iseq_neigh+1;
    iseqB = iseq_end;
  }

  // #################### COMMENTING OUT: ########################
  //// The Monte-Carlo move reverses the order of the chain in the interval
  //// from [iseqA, iseqB], inclusive
  //ReverseInc(iseqA, iseqB);
  //double deltaE = hamiltonian.CalcEnergyChange(*this, iseqA, iseqB);
  //double x = RANDOM_REAL_0_1();
  //bool accept_move = true;
  //if (x > (1.0 / (exp(deltaE) + 1.0))) // same as (x > 0.5(1-tanh(deltaE/2)))
  //  accept_move = false;
  //if (accept_move)
  //  energy += deltaE;
  //else
  //  // Undo the move.
  //  // Reverse the order of the chain in the interval from [iseqA, iseqB]
  //  ReverseInc(iseqA, iseqB);



  // #################### NEW STRATEGY: ##########################
  // First, calculate the change in energy if the move were accepted.
  // Do this BEFORE actually implementing the move.
  // (The MC moves can be expensive (O(N)), so implement the move only if it is
  //  accepted. This saves time, especially if the acceptance proability is low)

  double deltaE = hamiltonian.CalcEnergyChange(*this, iseqA, iseqB);
  bool accept_move = true;


  //Metropolis acceptance criteria  <--  COMMENTING OUT
  //if (deltaE > 0.0) {
  //  double x=RANDOM_REAL_0_1();
  //  if (x > exp(- deltaE))
  //    accept_move = false;
  //}


  // Glauber acceptance criteria
  double x = RANDOM_REAL_0_1();
  if (x > (1.0 / (exp(deltaE) + 1.0))) // same as (x > 0.5(1-tanh(deltaE/2)))
    accept_move = false;

  if (accept_move) {
    // The Monte-Carlo move reverses the order of the chain in the interval
    // from [iseqA, iseqB], inclusive
    ReverseInc(iseqA, iseqB); // <-- takes O(N) time
    energy += deltaE;
  }

  return accept_move;

} //MoveGenerator::Move()




// Reverse the order of the chain in the interval
// from [iseqA, iseqB], inclusive

void NDmansfield::
ReverseInc(long iseqA, long iseqB) {
  assert(iseqA <= iseqB);

  long num_swaps = (1 + iseqB - iseqA) / 2;  // integer div rounds down
  for (long ii=0; ii < num_swaps; ii++) {
    long iseq_a = iseqA + ii;           // location in the sequence
    long iseq_b = iseqB - ii;           // location in the sequence

    // Each lattice site has a unique index (independent of the polymer shape)
    // which represents the location in the 1-D table where that site is stored.
    long iloc_a = aIlocFromIseq[iseq_a]; // which lattice cell has iseq_a?
    long iloc_b = aIlocFromIseq[iseq_b]; // which lattice cell has iseq_b?
    // (You can use the GetLocation() function to look up the coordinates of
    //  these lattice sites.)
    //
    // Check the integrity of the inverse-lookup table:
    assert((0<=iseq_a) and (iseq_a <= poly_length));
    assert((0<=iseq_b) and (iseq_b <= poly_length));
    assert((0<=iloc_a) and (iloc_a <= num_cells));
    assert((0<=iloc_b) and (iloc_b <= num_cells));
    assert(aIseqFromIloc[iloc_a] == iseq_a);
    assert(aIseqFromIloc[iloc_b] == iseq_b);
    // Now swap the monomers located at iseq_a, iseq_b along the chain.
    aIlocFromIseq[iseq_a] = iloc_b;
    aIlocFromIseq[iseq_b] = iloc_a;
    aIseqFromIloc[iloc_a] = iseq_b;
    aIseqFromIloc[iloc_b] = iseq_a;
  }
} // NDmansfield::ReverseInc(long iseqA, long iseqB)




// Check to see whether the chain makes a closed loop
bool NDmansfield::IsCyclic() {
  long *aCoordsEndA = new long [g_dim];
  long *aCoordsEndB = new long [g_dim];
  CoordsFromIseq(aCoordsEndA, 0);
  CoordsFromIseq(aCoordsEndB, poly_length-1);
  long distsq = 0;
  for (int d=0; d < g_dim; d++) {
    long delta = aCoordsEndA[d] - aCoordsEndB[d];
    distsq += delta*delta;
  }
  return (distsq == 1);
  delete [] aCoordsEndA;
  delete [] aCoordsEndB;
}



//  -------- File IO --------




//Read the current state of the system from a file.
istream& operator >> (istream& in, NDmansfield& ndmansfield)
{
  long *aCoords = new long [g_dim];

  ndmansfield.Clear(); // Start with an empty lattice (polymer length 0)
  long iseq;
  for(iseq=0; iseq < ndmansfield.PolyLength(); ++iseq)
  {
    string line;
    getline(in, line);
    if (! in)
      break;
    stringstream line_ss(line);
    long d;
    for (d=0; d < g_dim; d++) {
      line_ss >> aCoords[d];
      if (! line_ss)
	break;
    }

    if (d != g_dim) {
      if (d>0) {
        stringstream errmsg;
        errmsg << "Error near line "<<iseq+1<<":\n"
               << "Expected " << g_dim 
               << " numbers on each line.\n";
        throw InputErr(errmsg.str());
      }
      else {    // skip over blank lines
	iseq--;
        continue;
      }
    }

    // Check for out-of-bounds errors:
    for (int d=0; d < g_dim; d++)
    {
      if (aCoords[d] <= 0) {
        stringstream errmsg;
        errmsg << "Error near line "<<iseq+1<<":\n"
               << "  Coordinates must be strictly positive integers.\n";
        throw InputErr(errmsg.str());
      }
      else if (ndmansfield.aSize[2]-1 <= aCoords[d]) {
        stringstream errmsg;
        errmsg << "Error near line "<<iseq+1<<":\n"
               << "  These coordinates lie outside the dimensions of the lattice\n"
               << "  (defined by the \"-box\" arguments), which are currently:\n  ";
        for (d = 0; d < g_dim-1; d++)
          errmsg << ndmansfield.aSize[d]-2 << " ";
        errmsg << ndmansfield.aSize[g_dim-1]-2 << "\n";
        throw InputErr(errmsg.str());
      }
    }

    //Now (finally) store the data:
    long iloc = ndmansfield.IlocFromCoords(aCoords);
    assert(iloc != NDmansfield::NOT_FOUND);
    ndmansfield.aIlocFromIseq[iseq] = iloc;
    ndmansfield.aIseqFromIloc[iloc] = iseq;
  } //for(long iseq=0; iseq < ndmansfield.num_cells; ++iseq)
  ndmansfield.poly_length = iseq;

  delete [] aCoords;

  return in;

} //istream& operator >> (istream& in,  NDmansfield& ndmansfield)



//Write the current state of the system to a file.
ostream& operator << (ostream& out, NDmansfield const& ndmansfield)
{
  long *aCoords = new long [g_dim];
  for(long iseq=0; iseq < ndmansfield.poly_length; ++iseq) {
    ndmansfield.CoordsFromIseq(aCoords, iseq);
    for (int d=0; d < g_dim-1; d++)
      out << aCoords[d] << " ";
    out << aCoords[g_dim-1] << "\n";
  }
  out << "\n";
  delete [] aCoords;
  return out;
}





//  -------- Memory allocation --------




NDmansfield::NDmansfield(const long *set_size) {
  assert(g_dim > 1);
  aSize = new long [g_dim];

  assert(set_size);
  for (int d=0; d < g_dim; d++)
    aSize[d] = set_size[d];

  // Create an array of boolean variables 
  // containing that number of entries
  num_cells=1;
  for(int d=0; d < g_dim; ++d)
    num_cells *= aSize[d]; // "*=" means multiply by, then store
  assert(num_cells > 0);

  aIseqFromIloc = new long [num_cells];
  aIlocFromIseq = new long [num_cells];

  if ((! aIseqFromIloc) || (! aIlocFromIseq)) {
    stringstream err_msg;
    err_msg << "Error in memory allocation: file \""
            <<__FILE__<<"\":"<<__LINE__<<"\n"
            << "  Unable to allocate an array of size "<<num_cells<<"\n"<<flush;
    throw InputErr(err_msg.str());
  }

  Clear();   // <--Fill the contents of aIseqFromIloc[] and aIlocFromIseq[]
             //    arrays with blanks

  DefaultPolyShape(); // <-- Then fill them with a polymer (in a simple conformation)

  aEndLoc       = new long [g_dim];
  aNeighLoc     = new long [g_dim];

} // NDmansfield::NDmansfield()


void NDmansfield::Clear() {
  for (long i=0; i < num_cells; ++i) {
    aIseqFromIloc[i] = NDmansfield::NOT_FOUND;
    aIlocFromIseq[i] = NDmansfield::NOT_FOUND;
  }
}



NDmansfield::~NDmansfield()
{
  assert(aIseqFromIloc);
  assert(aIlocFromIseq);
  assert(aEndLoc);
  assert(aNeighLoc);
  assert(aSize);
  delete [] aIseqFromIloc;
  delete [] aIlocFromIseq;
  delete [] aEndLoc;
  delete [] aNeighLoc;
  delete [] aSize;
}





void
NDmansfield::DefaultPolyShape()
{
  // Generate a zig-zag shape space-filling curve which fills the rectangle,
  // ...leaving extra layer of empty lattice sites on each face of the rectangle
  // (to insure that checking for neighbors never exceeds boundaries of the box)
  // (If you plan to use periodic boundary conditions, then fill the whole box.)
  //
  // 2-dimensional example
  //
  // *   *   *   *   *   *   *   <-- sites only the boundary remain vacant
  //                          
  // *   *---*---*---*---*   *   <-- only a 5x5 box is filled with polymer
  //     |                           (in this example, with aSize={7,7})
  // *   *---*---*---*---*   *
  //                     |
  // *   *---*---*---*---*   *
  //     |                    
  // *   *---*---*---*---*   *
  //                     |    
  // *   *---*---*---*---*   *
  //                          
  // *   *   *   *   *   *   *

  long *aIncrDirection = new long [g_dim];
  long *aCoords = new long [g_dim];
  long *aEffSize = new long [g_dim];
  for (int d=0; d<g_dim; d++) {
    aIncrDirection[d] = 1;
    aCoords[d] = 0;
    aEffSize[d] = aSize[d]-2;
  }


  long iseq = 0;
  bool continue_looping = true;
  while (continue_looping) {

    //Store the coordinates:

    // First, coordinates should be incremented by 1 to avoid lying on 
    // the boundaries. (Lattice sites on boundaries must remain empty.)
    long aCoordsEff[g_dim];
    for (int d=0; d<g_dim; d++) {
      aCoordsEff[d] = aCoords[d]+1;
      assert((1 <= aCoordsEff[d]) && (aCoordsEff[d] < aSize[d]-1));
    }

    long iloc = IlocFromCoords(aCoordsEff);
    assert(iloc != NDmansfield::NOT_FOUND);
    aIlocFromIseq[iseq] = iloc;
    aIseqFromIloc[iloc] = iseq;

    //Now, "increment" the position of the cell (aCoords[]).
    // For example:
    //    if we are at (0,0,0), we go to (1,0,0)
    //    if we are at (1,0,0), we go to (2,0,0)
    //    if we are at (2,0,0), we go to (3,0,0)
    //                            .                 .
    //                            .                 .
    //                            .                 .
    //    if we are at (aEffSize[0]-1,0,0), we go to (aEffSize[0]-1,1,0)
    //    if we are at (aEffSize[0]-1,1,0), we go to (aEffSize[0]-2,1,0)
    //    if we are at (aEffSize[0]-2,1,0), we go to (aEffSize[0]-3,1,0)
    //                            .                 .
    //                            .                 .
    //                            .                 .
    //    if we are at (2,aEffSize[0]-1,0), we go to (1,aEffSize[0]-1,0)
    //    if we are at (1,aEffSize[0]-1,0), we go to (0,aEffSize[0]-1,0)
    //    if we are at (0,aEffSize[0]-1,0), we go to (0,aEffSize[0]-2,0)
    //    if we are at (1,aEffSize[0]-2,0), we go to (2,aEffSize[0]-2,0)
    //    if we are at (2,aEffSize[0]-2,0), we go to (3,aEffSize[0]-2,0)
    //                            .                 .
    //                            .                 .
    //                            .                 .
    // This is just like incrementing a multi-digit counter (with g_dim digits),
    // except that we reverse direction each time we exceed the boundaries, 
    // instead of looping back to 0.
    {
      bool carry_the_1 = true;
      int d=0;
      while (carry_the_1 && (d < g_dim))
      {
        //aCoords[d]++;
        if (((aIncrDirection[d] > 0) && (aCoords[d] == aEffSize[d]-1)) ||
            ((aIncrDirection[d] < 0) && (aCoords[d] == 0)))
        {
          //aCoords[d] = 0;
          aIncrDirection[d] = -aIncrDirection[d];
          ++d;
        }
        else {
          carry_the_1 = false;
          aCoords[d] += aIncrDirection[d];
        }
      }
      if (d >= g_dim)
        continue_looping = false;

    } //finished incrementing

    iseq++;
  } //while (continue_looping)

  poly_length = iseq;

  delete [] aIncrDirection;
  delete [] aCoords;
  delete [] aEffSize;
} //NDmansfield::DefaultPolyShape()
