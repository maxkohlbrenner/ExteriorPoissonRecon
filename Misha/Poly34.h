// poly34.h : solution of cubic and quartic equation
// (c) Khashin S.I. http://math.ivanovo.ac.ru/dalgebra/Khashin/index.html
// khash2 (at) gmail.com

#ifndef POLY34_INCLUDED
#define POLY34_INCLUDED

namespace Poly34
{
																							// x - array of size 2
																							// return 2: 2 real roots x[0], x[1]
																							// return 0: pair of complex roots: x[0]�i*x[1]
	static int SolveP2( double *x , double a , double b );									// solve equation x^2 + a*x + b = 0

																							// x - array of size 3
																							// return 3: 3 real roots x[0], x[1], x[2]
																							// return 1: 1 real root x[0] and pair of complex roots: x[1]�i*x[2]
	static int SolveP3( double *x , double a , double b , double c );						// solve cubic equation x^3 + a*x^2 + b*x + c = 0

																							// x - array of size 4
																							// return 4: 4 real roots x[0], x[1], x[2], x[3], possible multiple roots
																							// return 2: 2 real roots x[0], x[1] and complex x[2]�i*x[3], 
																							// return 0: two pair of complex roots: x[0]�i*x[1],  x[2]�i*x[3], 
	static int SolveP4( double *x , double a , double b , double c , double d );			// solve equation x^4 + a*x^3 + b*x^2 + c*x + d = 0 by Dekart-Euler method

																							// x - array of size 5
																							// return 5: 5 real roots x[0], x[1], x[2], x[3], x[4], possible multiple roots
																							// return 3: 3 real roots x[0], x[1], x[2] and complex x[3]�i*x[4], 
																							// return 1: 1 real root x[0] and two pair of complex roots: x[1]�i*x[2],  x[3]�i*x[4], 
	static int SolveP5( double *x , double a , double b , double c , double d , double e );	// solve equation x^5 + a*x^4 + b*x^3 + c*x^2 + d*x + e = 0
}

#include "Poly34.inl"
#endif // POLY34_INCLUDED