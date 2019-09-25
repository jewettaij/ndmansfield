#ifndef _HAMILTONIAN_H
#define _HAMILTONIAN_H


#include "ndmansfield.hpp" // Contains a description of the NDmansfield data 
                           // structure we will be using to store the state
                           // of our system as the simulation evolves.




class Hamiltonian
{
  double bend_energy_coeff;
  double twist_energy_coeff;
public:
  Hamiltonian(NDmansfield& ndmansfield, 
              double set_bend_energy,
              double set_twist_energy);
  ~Hamiltonian();

  double CalcEnergy(NDmansfield& ndmansfield); // Calculate total energy of system

  // Calculate the change in energy of the system upon a move:  
  double CalcEnergyChange(NDmansfield& ndmansfield,
                          long iseqA,   // Interval of chain (inclusive)
                          long iseqB);  // which is reversed during this move

 private:
  long *aCoords_im3;
  long *aCoords_im2;
  long *aCoords_im1;
  long *aCoords_i;
  long *aCoords_ip1;
  long *aCoords_ip2;
  long *aCoords_ip3;

  long *aCoords_iseqAm3;
  long *aCoords_iseqAm2;
  long *aCoords_iseqAm1;  // <--- location of monomer iseqA-1
  long *aCoords_iseqA;    // <--- location of monomer iseqA
  long *aCoords_iseqAp1;  // <--- location of monomer iseqA+1
  long *aCoords_iseqAp2;

  long *aCoords_iseqBm2;
  long *aCoords_iseqBm1;  // <--- location of monomer iseqB-1
  long *aCoords_iseqB;    // <--- location of monomer iseqB
  long *aCoords_iseqBp1;  // <--- location of monomer iseqB+1
  long *aCoords_iseqBp2;
  long *aCoords_iseqBp3;

  long *v1;
  long *v2;
  long *v3;
  long *v1x2;



  inline long CosAngleABC(const long *aCoordsA, 
                          const long *aCoordsB,
                          const long *aCoordsC) {
    // Return the dot-product between aCoordsA-aCoordsB and aCoordsC-aCoordsB
    long total = 0;
    for(int d=0; d<g_dim; ++d)
      total += (aCoordsA[d] - aCoordsB[d]) * (aCoordsC[d] - aCoordsB[d]);
    return total;
  }




  inline long CosTorsionABCD(const long *aCoordsA, 
                             const long *aCoordsB,
                             const long *aCoordsC,
                             const long *aCoordsD) {
    assert(g_dim == 3);
    // Returns the triple-product between these three vectors:
    //   v1 = (aCoordsB - aCoordsA)
    //   v2 = (aCoordsC - aCoordsB)
    //   v3 = (aCoordsD - aCoordsC)
    // (The "triple product" is the volume of a parallelepiped
    //  formed by the three vectors.  For vectors v1,v2,v3, the formula for 
    //  the triple product is : (v1 x v2).v3.  It is invariant with respect
    //  to cyclic permutations, so it is equivalent to (v2 x v3).v1
    //  If positive, it means that v1, v2, v3 connected together in succession
    //  form a right-handed helix.  If negative, then a left-handed helix.)
    long total = 0;
    for(int d=0; d<g_dim; ++d) {
      v1[d] = aCoordsB[d] - aCoordsA[d];
      v2[d] = aCoordsC[d] - aCoordsB[d];
      v3[d] = aCoordsD[d] - aCoordsC[d];
    }
    v1x2[0] = v1[1]*v2[2] - v1[2]*v2[1];
    v1x2[1] = v1[2]*v2[0] - v1[0]*v2[2];
    v1x2[2] = v1[0]*v2[1] - v1[1]*v2[0];
    for(int d=0; d<g_dim; ++d)
      total += v1x2[d] * v3[d];
    return total;
  }
};


#endif //#ifndef _HAMILTONIAN_H
