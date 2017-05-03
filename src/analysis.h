#ifndef _ANALYSIS_H
#define _ANALYSIS_H

#include "ndmansfield.h"


void CountBondDirections(long *aBondCountPlusXYZ,
			 long *aBondCountMinusXYZ,
			 long N,
                         long **aaXYZ);

void CountBondDirections(long *aBondCountPlusXYZ,
			 long *aBondCountMinusXYZ,
			 const NDmansfield &ndmansfield);

void SpatialCorrelation(double *aSpatialCorrelation,
                        long N,
                        long **aaXYZ);

void SpatialCorrelation(double *aSpatialCorrelation,
			const NDmansfield &ndmansfield);


#endif //#ifndef _ANALYSIS_H
