#include <cmath>
#include <cassert>
#include "ndmansfield.hpp"
#include "analysis.h"
using namespace std;



// Calculate how spatial position correlates with distance along the chain

void SpatialCorrelation(double *aSpatialCorrelation,
                        long N,
                        long **aaXYZ)
{
  assert(aSpatialCorrelation);
  // Now calculate the centroid (the average x,y,z position of the polymer)
  double xyz_ave[g_dim];
  for (int d = 0; d < g_dim; d++)
    xyz_ave[d] = 0.0;
  for (long i = 0; i < N; i++)
    for (int d = 0; d < g_dim; d++)
      xyz_ave[d] += aaXYZ[i][d];
  for (int d = 0; d < g_dim; d++)
    xyz_ave[d] /= N;
  // Then calculate the standard deviations in the x,y,z positions:
  double sigma_xyz[g_dim];
  for (int d = 0; d < g_dim; d++)
    sigma_xyz[d] = 0.0;
  for (long i = 0; i < N; i++)
    for (int d = 0; d < g_dim; d++)
      sigma_xyz[d] += (aaXYZ[i][d] - xyz_ave[d])*(aaXYZ[i][d] - xyz_ave[d]);
  for (int d = 0; d < g_dim; d++)
    sigma_xyz[d] = sqrt(sigma_xyz[d] / N);
  double i_ave = 0.5*N;
  double sigma_i = 0.0;
  for (long i = 0; i < N; i++)
    sigma_i += (i-i_ave)*(i-i_ave);
  sigma_i = sqrt(sigma_i / N);

  // Finally calculate the spatial correlation of x,y,z with i
  for (long i = 0; i < N; i++) {
    for (int d = 0; d < g_dim; d++) {
      aSpatialCorrelation[d] += (i-i_ave) * (aaXYZ[i][d] - xyz_ave[d]);
    }
  }
  for (int d = 0; d < g_dim; d++)
    aSpatialCorrelation[d] /= (N*sigma_xyz[d]*sigma_i);
}




void SpatialCorrelation(double *aSpatialCorrelation,
                        const NDmansfield &ndmansfield)
{
  // Load the coordinates into a 2-D array
  long N = ndmansfield.PolyLength();
  long **aaXYZ = new long* [N];
  for (long i = 0; i < N; i++) {
    aaXYZ[i] = new long [g_dim];
    //load coordinates for ith atom into aaXYZ[i]
    ndmansfield.CoordsFromIseq(aaXYZ[i], i);
  }

  SpatialCorrelation(aSpatialCorrelation, N, aaXYZ);

  // Cleanup
  for (long i = 0; i < N; i++)
    delete [] aaXYZ[i];
  delete [] aaXYZ;
}




void CountBondDirections(long *aBondCountPlusXYZ,
                         long *aBondCountMinusXYZ,
                         const NDmansfield &ndmansfield)
{
  // Load the coordinates into a 2-D array.
  // (This is not really necessary, but it makes the code easier to read.)
  long N = ndmansfield.PolyLength();
  long **aaXYZ = new long* [N];  
  for (long i = 0; i < N; i++) {
    aaXYZ[i] = new long [g_dim];
    //load coordinates for ith atom into aaXYZ[i]
    ndmansfield.CoordsFromIseq(aaXYZ[i], i);
  }

  CountBondDirections(aBondCountPlusXYZ,
                      aBondCountMinusXYZ,
                      N,
                      aaXYZ);

  // Cleanup
  for (long i = 0; i < N; i++)
    delete [] aaXYZ[i];
  delete [] aaXYZ;
}


// Calculate how spatial position correlates with distance along the chain

void CountBondDirections(long *aBondCountPlusXYZ,
                         long *aBondCountMinusXYZ,
                         long N,
                         long **aaXYZ)
{
  assert(aBondCountPlusXYZ);
  assert(aBondCountMinusXYZ);

  // Count the number of bonds in the +X,+Y,+Z and -X,-Y,-Z directions:
  for (int d = 0; d < g_dim; d++) {
    aBondCountPlusXYZ[d] = 0;
    aBondCountMinusXYZ[d] = 0;
  }
  for (long i = 0; i+1 < N; i++) {
    int d;
    for (d = 0; d < g_dim; d++) {
      if ((aaXYZ[i+1][d] - aaXYZ[i][d]) == 1) {
        aBondCountPlusXYZ[d] += 1;
        break;
      }
      else if ((aaXYZ[i+1][d] - aaXYZ[i][d]) == -1) {
        aBondCountMinusXYZ[d] += 1;
        break;
      }
    }
    assert(d < g_dim);
  }
}

