#ifndef _RANDOM_GEN_H
#define _RANDOM_GEN_H

#include <cassert>            // defines assert(), optional debugging utility.

#include <cstdlib>            // needed for random-number generation
                              // (specifically getenv() and seed48())
#include <ctime>              // required for "struct timezone" and
                              // "gettimeofday()" used to set the randomseed
#include <cmath>              // defines exp(), log(), sqrt(), cos(), sin()

#include <iostream>
using namespace std;



//RANDOM_INT(n) returns uniformly distributed integers ranging from 0 to n-1.
//The maximum allowed value of "n" is 2^31

inline long RANDOM_INT(unsigned long n)
{
  //assert(n < 2147483648-1); // Unnecessary? Probably...Check that n < 2^31-1
  return lrand48() % n; //Note: lrand48 returns a number between 0 and 2^31
}







//RANDOM_REAL_0_1() returns uniformly distributed "Real" numbers over [0, 1)
inline double RANDOM_REAL_0_1()
{
  return drand48(); //for some reason this is better than lrand48() / RAND_MAX
}







//RANDOM_GAUSSIAN() returns a random number in the range
//                            (-infinity, +infinity)
//with distrubition of P(x) proportional to exp( - x^2 / 2 )
//   (the "standard normal distribution".)
//The number should be Gaussian distributed, and have a variance, <x^2> = 1.

inline double RANDOM_GAUSSIAN()
{
  double u = RANDOM_REAL_0_1(); // generates a random number in the range [0,1)
  u = 1.0 - u;                  // After this command, u is in the range (0,1]
                                // We want u>0 because log(u) not defined at 0

  // The following line should generate r with a probability density of
  // P(r) proportional to  r * exp(- r^2/2)

  double r =  sqrt( - log (u) * 2.0 );
  //cerr << "u=" << u << ", log(u)=" << log(u) << ", sqrt(-log(u))="
  //     << sqrt(-log(u)) << endl;


  // Now generate another random number between 0 and 2*pi
  double theta = (2.0*M_PI) * RANDOM_REAL_0_1();

  // Convert from polar (r,theta) to cartesian coordinates (x,y)
  // The "x" and "y" values should both be Gaussian distributed.
  // (Probability density P(x) should be proportional to exp(-lambda x^2))

  //double return_val = r*cos(theta);
  //cerr << "    theta=" << theta << ", cos(theta)=" << cos(theta)
  //     << ", r=" << r << ", r*cos(theta)=" << return_val << endl;

  return r * cos( theta );

} //RANDOM_GAUSSIAN()







// RANDOM_EXPONENTIAL() returns a number in the range [0,infinity)
// with an exponential probability distribution  P(x) = exp(-x)
//    (Scale the final answer to get different decay rates.
//     P(x) = (1/K) exp(-K x))


inline double RANDOM_EXPONENTIAL()
{
  double u = RANDOM_REAL_0_1(); // generates a random number in the range [0,1)
  u = 1.0 - u;                  // After this command, u is in the range (0,1]
                                // We want u>0 because log(u) not defined at 0
  return -log(u);
}




// RANDOM_INIT()
// Invoke this function before calling anything else:
// This function chooses the random-seed for the random
// number generator based on the time.
// This insures that the program will never reproduce exactly the same
// results when run twice (unless run within 1 second of eachother).

inline void RANDOM_INIT(long seed=-1)
{
  time_t timer;
  struct tm y2k;
  double seconds;
  y2k.tm_hour = 0;   y2k.tm_min = 0; y2k.tm_sec = 0;
  y2k.tm_year = 100; y2k.tm_mon = 0; y2k.tm_mday = 1;
  time(&timer);  /* get current time; same as: timer = time(NULL)  */
  seconds = difftime(timer,mktime(&y2k));
  if (seed < 0) {
    seed = static_cast<long>(floor(seconds));
  }
  srand48( seed );
}



#endif //#ifndef _RANDOM_GEN_H
