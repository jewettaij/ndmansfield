#include <cassert>
#include "ndmansfield.h"
#include "hamiltonian.h"

// Right now the hamiltonian does nothing (yet)





// This function calculates the (total) energy of the system
double 
Hamiltonian::CalcEnergy(NDmansfield &ndmansfield)
{
  long num_bends = 0;
  long num_twists = 0;
  for (long i=0; i+2 < ndmansfield.poly_length; i++) {
    ndmansfield.CoordsFromIseq(aCoords_i,   i);
    ndmansfield.CoordsFromIseq(aCoords_ip1, i+1);
    ndmansfield.CoordsFromIseq(aCoords_ip2, i+2);
    num_bends += (1 + CosAngleABC(aCoords_i, aCoords_ip1, aCoords_ip2));

    // Torsional energy
    if ((twist_energy_coeff != 0.0) && 
	(g_dim == 3) &&
	(i+3 < ndmansfield.poly_length)) {
      ndmansfield.CoordsFromIseq(aCoords_ip3, i+3);
      num_twists += (1 - CosTorsionABCD(aCoords_i,
					aCoords_ip1,
					aCoords_ip2, 
					aCoords_ip3));
    }
  }
  return ((bend_energy_coeff * num_bends) + 
	  (twist_energy_coeff * num_twists));
} //Hamiltonian::CalcEnergy()






// This function calculates the change in energy due to a Mansfield move.
// A Mansfield Monte-Carlo move consists of reversing the order of all of
// the monomers in an interval containing one of the end-points of the chain.
// (In this function, the interval is denoted by [iseqA, iseqB] inclusive.)
// This function reports the difference in the bend-energy and twist-energies
// before and after the move.
double 
Hamiltonian::CalcEnergyChange(NDmansfield &ndmansfield,
                              // Interval where the chain was reversed 
                              // during the move:
                              long iseqA,
                              long iseqB)
{
  // Chains containing bends must pay an energetic penalty.
  // After the monte carlo move, bonds will be broken and formed near
  // the two interval endpoints.  
  // It would be slow to recalculate the entire energy of the chain.
  // Instead we only look for bends in the path in the viscinity of the two
  // ends that would have been distroyed or created by the monte-carlo move.
  // (The bonds and position of the beads in the chain between these 
  //  two endpoints will not be effected by the move.)
  long num_bends_before = 0;
  long num_bends_after = 0;
  long num_twists_before = 0;
  long num_twists_after = 0;

  // Coordinates of the two end-points and their successors & predecessors:

  assert(iseqA < iseqB);
  assert(iseqB-2 >= 0);                      //(all loops have length>=4)
  assert(iseqA+2 < ndmansfield.poly_length); //(all loops have length>=4)

  ndmansfield.CoordsFromIseq(aCoords_iseqA,   iseqA);
  ndmansfield.CoordsFromIseq(aCoords_iseqAp1, iseqA+1);
  ndmansfield.CoordsFromIseq(aCoords_iseqAp2, iseqA+2);
  ndmansfield.CoordsFromIseq(aCoords_iseqBm2, iseqB-2);
  ndmansfield.CoordsFromIseq(aCoords_iseqBm1, iseqB-1);
  ndmansfield.CoordsFromIseq(aCoords_iseqB,   iseqB);

  // Either iseqA or iseqB must lie at either end of the polymer.
  // This is to insure that the polymer does not have any discontinuous jumps
  // after the interval [iseqA,iseqB] is reversed during the Monte-Carlo step.


  if (iseqA == 0) {
    // NOTATION:
    // Comparison with figure1 from the Mansfield JCP2006 paper:
    // What I am calling iseqA, and iseqB+1 here correspond to "A" and "B"
    // in figure 1a of the Mansfield J.Chem.Phys. 2006 paper.

    // After the proposed move, the bond between iSeqB and iSeqB+1 
    // (ie. "B-1" and "B" from figure 1a of the JCP2006 paper) would be broken,
    // and replaced with a bond between iSeqA and iSeqB+1.
    // We must calculate the energy of these angles (which would be deleted)
    //    iseqB-1, iseqB,   iseqB+1
    //    iseqB,   iseqB+1, iseqB+2
    //and subtract them from the energy of the new angles which would be created
    //    iseqB+2, iseqB+1, iseqA
    //    iseqB+1, iseqA,   iseqA+1
    assert(iseqB+1 < ndmansfield.poly_length);  //iseqB+1 is iseqA's neighbor
    ndmansfield.CoordsFromIseq(aCoords_iseqBp1, iseqB+1);
    num_bends_before += (1 + CosAngleABC(aCoords_iseqBp1,
                                         aCoords_iseqB, 
                                         aCoords_iseqBm1));
    num_bends_after  += (1 + CosAngleABC(aCoords_iseqBp1, 
                                         aCoords_iseqA, 
                                         aCoords_iseqAp1));
   
    if ((twist_energy_coeff != 0.0) && (g_dim == 3)) {
      // Do the same thing for the torsion angles:
      num_twists_before += (1 - CosTorsionABCD(aCoords_iseqBp1, 
					       aCoords_iseqB,
					       aCoords_iseqBm1, 
					       aCoords_iseqBm2));
      num_twists_after  += (1 - CosTorsionABCD(aCoords_iseqBp1, 
					       aCoords_iseqB,
					       aCoords_iseqA, 
					       aCoords_iseqAp1));
    }

    if (iseqB+2 < ndmansfield.poly_length) {
      ndmansfield.CoordsFromIseq(aCoords_iseqBp2, iseqB+2);
      num_bends_before += (1 + CosAngleABC(aCoords_iseqBp2, 
                                           aCoords_iseqBp1,
                                           aCoords_iseqB));
      num_bends_after  += (1 + CosAngleABC(aCoords_iseqBp2, 
                                           aCoords_iseqBp1,
                                           aCoords_iseqA));
      if ((twist_energy_coeff != 0.0) && (g_dim == 3)) {
	// torsion:
	num_twists_before += (1 - CosTorsionABCD(aCoords_iseqBp2, 
						 aCoords_iseqBp1,
						 aCoords_iseqB, 
						 aCoords_iseqBm1));
	num_twists_after  += (1 - CosTorsionABCD(aCoords_iseqBp2, 
						 aCoords_iseqBp1,
						 aCoords_iseqA, 
						 aCoords_iseqAp1));
	if (iseqB+3 < ndmansfield.poly_length) {
	  ndmansfield.CoordsFromIseq(aCoords_iseqBp3, iseqB+3);
	  num_twists_before += (1 - CosTorsionABCD(aCoords_iseqBp3,
						   aCoords_iseqBp2,
						   aCoords_iseqBp1,
						   aCoords_iseqB));
	  num_twists_after  += (1 - CosTorsionABCD(aCoords_iseqBp3, 
						   aCoords_iseqBp2,
						   aCoords_iseqBp1,
						   aCoords_iseqA));
	}
      }
    }
  } //if (iseqA == 0)




  else if (iseqB == ndmansfield.poly_length-1) {
    // NOTATION:
    // Comparison with figure1 from the Mansfield JCP2006 paper:
    // What I am calling iseqA-1, and iseqB here correspond to "D" and "C"
    // in figure 1b of the Mansfield J.Chem.Phys. 2006 paper.

    // After the proposed move, the bond between iSeqA and iSeqA-1
    // (ie. "D+1" and "D" from figure 1b of the JCP2006 paper) would be broken,
    // and replaced with a bond between iSeqB and iSeqA-1.
    // We must calculate the energy of these angles (which would be deleted)
    //    iseqA-2, iseqA-1, iseqA
    //    iseqA-1, iseqA,   iseqA+1
    //and subtract them from the energy of the new angles which would be created
    //    iseqA-2, iseqA-1, iseqB
    //    iseqA-1, iseqB,   iseqB-1
    assert(iseqA-1 >= 0);    // because iseqA-1 is a neigbor of iseqB
    ndmansfield.CoordsFromIseq(aCoords_iseqAm1, iseqA-1);
    num_bends_before += (1 + CosAngleABC(aCoords_iseqAm1, 
                                         aCoords_iseqA, 
                                         aCoords_iseqAp1));
    num_bends_after  += (1 + CosAngleABC(aCoords_iseqAm1, 
                                         aCoords_iseqB, 
                                         aCoords_iseqBm1));

    if ((twist_energy_coeff != 0.0) && (g_dim == 3)) {
      // Do the same thing for the torsion angles:
      num_twists_before += (1 - CosTorsionABCD(aCoords_iseqAm1, 
					       aCoords_iseqA,
					       aCoords_iseqAp1, 
					       aCoords_iseqAp2));
      num_twists_after +=  (1 - CosTorsionABCD(aCoords_iseqAm1, 
					       aCoords_iseqB, 
					       aCoords_iseqBm1, 
					       aCoords_iseqBm2));
    }
    if (iseqA-2 >= 0) {
      ndmansfield.CoordsFromIseq(aCoords_iseqAm2, iseqA-2);
      num_bends_before += (1 + CosAngleABC(aCoords_iseqAm2, 
                                           aCoords_iseqAm1, 
                                           aCoords_iseqA));
      num_bends_after  += (1 + CosAngleABC(aCoords_iseqAm2, 
                                           aCoords_iseqAm1,
                                           aCoords_iseqB));

      if ((twist_energy_coeff != 0.0) && (g_dim == 3)) {
	// torsion:
	num_twists_before += (1 - CosTorsionABCD(aCoords_iseqAm2, 
						 aCoords_iseqAm1,
						 aCoords_iseqA, 
						 aCoords_iseqAp1));
	num_twists_after  += (1 - CosTorsionABCD(aCoords_iseqAm2, 
						 aCoords_iseqAm1,
						 aCoords_iseqB, 
						 aCoords_iseqBm1));
	if (iseqA-3 >= 0) {
	  ndmansfield.CoordsFromIseq(aCoords_iseqAm3, iseqA-3);
	  num_twists_before += (1 - CosTorsionABCD(aCoords_iseqAm3, 
						   aCoords_iseqAm2,
						   aCoords_iseqAm1, 
						   aCoords_iseqA));
	  num_twists_after  += (1 - CosTorsionABCD(aCoords_iseqAm3, 
						   aCoords_iseqAm2,
						   aCoords_iseqAm1, 
						   aCoords_iseqB));
	}
      }
    }
  } //else if (iseqB == ndmansfield.poly_length)

  else 
    // The Monte-Carlo moves should only reverse interval from iseqA to iseqB 
    // if one of those two end points lies at one of the ends of the chain.
    // (Otherwise, the resulting chain would likely have discontinuos jumps 
    //  of length > 1 lattice site.)
    assert(false);


  return (bend_energy_coeff * (num_bends_after - num_bends_before) 
	  +
	  twist_energy_coeff * (num_twists_after - num_twists_before));

} //Hamiltonian::CalcEnergy()




Hamiltonian::Hamiltonian(NDmansfield& ndmansfield, 
                         double set_bend_energy,
                         double set_twist_energy)
{
  bend_energy_coeff = set_bend_energy;
  twist_energy_coeff = set_twist_energy;

  aCoords_im3 = new long [g_dim];
  aCoords_im2 = new long [g_dim];
  aCoords_im1 = new long [g_dim];
  aCoords_i = new long [g_dim];
  aCoords_ip1 = new long [g_dim];
  aCoords_ip2 = new long [g_dim];
  aCoords_ip3 = new long [g_dim];

  aCoords_iseqAm3 = new long [g_dim];
  aCoords_iseqAm2 = new long [g_dim];
  aCoords_iseqAm1 = new long [g_dim];
  aCoords_iseqA = new long [g_dim];
  aCoords_iseqAp1 = new long [g_dim];
  aCoords_iseqAp2 = new long [g_dim];

  aCoords_iseqBm2 = new long [g_dim];
  aCoords_iseqBm1 = new long [g_dim];
  aCoords_iseqB = new long [g_dim];
  aCoords_iseqBp1 = new long [g_dim];
  aCoords_iseqBp2 = new long [g_dim];
  aCoords_iseqBp3 = new long [g_dim];

  v1 = new long [g_dim];
  v2 = new long [g_dim];
  v3 = new long [g_dim];
  v1x2 = new long [g_dim];
}




Hamiltonian::~Hamiltonian()
{
  delete [] aCoords_im3;
  delete [] aCoords_im2;
  delete [] aCoords_im1;
  delete [] aCoords_i;
  delete [] aCoords_ip1;
  delete [] aCoords_ip2;
  delete [] aCoords_ip3;

  delete [] aCoords_iseqAm3;
  delete [] aCoords_iseqAm2;
  delete [] aCoords_iseqAm1;
  delete [] aCoords_iseqA;
  delete [] aCoords_iseqAp1;
  delete [] aCoords_iseqAp2;

  delete [] aCoords_iseqBm2;
  delete [] aCoords_iseqBm1;
  delete [] aCoords_iseqB;
  delete [] aCoords_iseqBp1;
  delete [] aCoords_iseqBp2;
  delete [] aCoords_iseqBp3;

  delete [] v1;
  delete [] v2;
  delete [] v3;
  delete [] v1x2;
}



