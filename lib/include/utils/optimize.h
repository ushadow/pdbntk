/*
 * optimize.h --- Brent and Powell function minimizers
 *
 *  Copyright (C) 2008, Wouter Boomsma, The Bioinformatics Centre, University of Copenhagen.
 *
 *  This file is part of Mocapy++.
 *
 *  Mocapy++ is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Mocapy++ is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with Mocapy++.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef OPTIMIZE_H
#define OPTIMIZE_H

#include "vector_nD.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <valarray>
#include <cstdarg>
#include <limits>

namespace mocapy {
	inline double NR_SIGN(double a, double b) { return ((b>=0.0) ? fabs(a) : -fabs(a)); }
	inline double NR_MAX(double a, double b) { return ((a>b) ? a : b ); } 
	
	// Global functions used by Bracket and Brent
	inline double len(double v) {return v;}
	inline double len(Vector_nD v) {return v.norm();}     
	
	
	//////////////////// Function ////////////////////
	// Supports various ways of function calls from optimize functions
	template <typename TYPE, typename FUNCTIONTYPE>
	class Function {
	public:
		 virtual double operator()(TYPE x) =0;
	
		 virtual ~Function(){};
	};
	
	template <typename TYPE, typename FUNCTIONTYPE>
	class ValueFunction: public Function<TYPE,FUNCTIONTYPE> {
	private:
		 FUNCTIONTYPE func;
	public:
		 ValueFunction(FUNCTIONTYPE func) {
			  this->func = func;          
		 }
		 
		 double operator()(TYPE x) {
			  return (*func)(x);
		 }
	};
	
	template <typename TYPE, typename FUNCTIONTYPE>
	class PointerFunction: public Function<TYPE,FUNCTIONTYPE> {
	private:
		 FUNCTIONTYPE func;
		 double **xPointer;
		 int size;
	public:
		 PointerFunction(FUNCTIONTYPE func, int size, double **xPointer) {
			  this->func = func;
			  this->xPointer = xPointer;
		 }
		 
		 double operator()(double x) {
			  *xPointer = x;
			  return (*func)();
		 }
	
		 double operator()(Vector_nD x) {
			  for (int i=0; i<size; i++) {
				   *xPointer[i] = x[i]; 
			  }
			  return (*func)();
		 }
	};
	
		 
	//////////////////// Bracket ////////////////////
	// Find points a, b, c that bracket the minimum of a given function
	template <typename TYPE>
	class Bracket {
	private:
		 // Constants
		 double goldStep;
		 double growLimit;
		 double tinyNumber;
		 int maxIterations;
		 
	public:
		 TYPE a, b, c;              // Three bracketing x values a<b<c or c<b<a
		 double fa, fb, fc;         // Corresponding function values
	
		 bool status;               // error status
		 int evaluations;           // number of function evaluations
	
		 Bracket(double(*f)(double), double a=0.0, double b=1.0) {
			  init(a, b);
			  createBracket(f);
		 }
	
		 template <typename FUNCTIONTYPE>
		 Bracket(Function<TYPE, FUNCTIONTYPE> f, double a, double b, double c) {
			  init(a, b);
			  this->a = a; this->fa = f(a);
			  this->b = b; this->fb = f(b);
			  this->c = c; this->fc = f(c);
			  this->status = true;
			  this->evaluations = 3;
		 }
	
		 template <typename FUNCTIONTYPE>
		 Bracket(FUNCTIONTYPE func, double a, double b, double c) {
			  ValueFunction<double, FUNCTIONTYPE> f(func);          
			  init(a, b);
			  this->a = a; this->fa = f(a);
			  this->b = b; this->fb = f(b);
			  this->c = c; this->fc = f(c);
			  this->status = true;
			  this->evaluations = 3;
		 }
		 
		 template <typename FUNCTIONTYPE>
		 Bracket(Function<double, FUNCTIONTYPE> &f, TYPE a=-1.0, TYPE b=1.0) {
			  init(a, b);
			  createBracket(f);
		 }
	
		 template <typename FUNCTIONTYPE>
		 Bracket(Function<Vector_nD, FUNCTIONTYPE> &f, Vector_nD a=Vector_nD(), Vector_nD b=Vector_nD()) {
			  init(a, b);
			  createBracket(f);
		 }
	
	//      template <typename FUNCTIONTYPE>     
		 void init(TYPE &a=TYPE(), TYPE &b=TYPE()) {
			  // Constants
			  this->goldStep = 1.618034;
			  this->growLimit = 110.0;
			  this->tinyNumber = 1e-21;
			  this->maxIterations = 50;
			  
	
			  this->a = a;
			  this->b = b;
		 
			  this->status = true;
			  this->evaluations = 0;
		 }
	
		 template <typename FUNCTIONTYPE>     
		 void createBracket(FUNCTIONTYPE &func) {
	//           std::cout << "********************************Bracketing\n";
			  TYPE unitVector = ((b-a)/(len(b-a)));
			  TYPE origin = a;
		 
			  fa = func(a); evaluations++;
			  fb = func(b); evaluations++;
		 
			  if (fb > fa) {
				   // Switch roles so f(a)>f(b) - downhill from a to b
				   TYPE tmp = a;
				   a = b;
				   b = tmp;
				   double ftmp = fa;
				   fa = fb;
				   fb = ftmp;
			  }
	
			  double ax = len(a-origin);
			  double bx = len(b-origin);
			  double cx;
	
			  double xmin = ax;
			  double fmin = fa;
			  
			  // Guess a value
			  cx = bx + (bx-ax)*goldStep;
			  fc = func(origin + unitVector * cx); evaluations++;
			  
			  int iterations=0;
			  while (fb >= fc) {
				   /* std::cout << ax << " " << bx << " " << cx << "\n"; */
				   /* std::cout << fa << " " << fb << " " << fc << "\n"; */
				   if (iterations > maxIterations) {
						status = false;
						a = origin + unitVector * xmin;
						fa = fmin;
						return;
				   }
				   iterations++;
	
				   double r = (bx-ax)*(fb-fc);
				   double q = (bx-cx)*(fb-fa);
				   double u = bx-((bx-cx)*q-(bx-ax)*r)/(2.0*NR_SIGN(NR_MAX(fabs(q-r), tinyNumber), q-r));
				   double ulim = bx+(cx-bx)*growLimit;
				   double fu;
		  
				   if ((bx-u)*(u-cx) > 0.0) {
						// Parabolic u is between b and c.
						// std::cout << "Parabolic u is between b and c.\n";
						fu=func(origin + unitVector * u); evaluations++;
						if (fu < fc) {
							 // std::cout << "Minimum is between b and c.\n";
							 // Minimum between b and c
							 ax=bx;
							 bx=u;
							 fa=fb;
							 fb=fu;
							 break;
						} else if (fu > fb) {
							 // std::cout << "Minimum is between a and u.\n";
							 // Minimum between a and u
							 cx=u;
							 fc=fu;
							 break;
						} 
						// std::cout << "Parabolic fit of no use - use default magnification\n";
						// Parabolic fit of no use - use default magnification
						u=cx+(cx-bx)*goldStep;
						fu=func(origin + unitVector * u); evaluations++;
				   } else if ((cx-u)*(u-ulim) > 0.0) {
						// Parabolic fit between c and allowed limit
						// std::cout << "Parabolic fit between c and allowed limit\n";
						fu=func(origin + unitVector * u); evaluations++;
						if (fu < fc) { 
							 bx = cx;
							 cx = u;
							 u = cx + (cx-bx)*goldStep;
							 fb = fc;
							 fc = fu;
							 fu = func(origin + unitVector * u); evaluations++;
	
						}
				   } else if ((u-ulim)*(ulim-cx) >= 0.0) {
						// set u to limit
						// std::cout << "set u to limit\n";
						u=ulim;
						fu=func(origin + unitVector * u); evaluations++;
				   } else {
						// Reject parabolic fit - use default magnification
						// std::cout << "Reject parabolic fit - use default magnification\n";
						u=cx+(cx-bx)*goldStep;
						fu=func(origin + unitVector * u); evaluations++;
				   }
	
				   if (fu < fmin) {
						fmin = fu;
						xmin = u;
				   }
				   
				   if (fu != fc) {   // Only update a if f(u) is different thatn f(c)
						ax = bx;
						fa = fb;
				   }
				   
				   bx = cx;
				   fb = fc;
	
				   cx = u;
				   fc = fu;
			  }
	
			  // Swap so a<b<c
			  if (ax > cx) {
				   c = origin + unitVector * ax;
				   b = origin + unitVector * bx;
				   a = origin + unitVector * cx;
				   double tmp = fa;
				   fa = fc;
				   fc = tmp;
			  } else {
				   a = origin + unitVector * ax;
				   b = origin + unitVector * bx;
				   c = origin + unitVector * cx;
			  }
	
			  if (fabs(fa - fb) < tinyNumber) { // If fa==fb call create bracket to get other boundary right
				   TYPE tmpx = a;
				   double tmpf = fa;
				   a = c;
				   fa = fc;
				   c = tmpx;
				   fc = tmpf;
				   createBracket(func);
			  }
	
			  return;
		 }
	
		 friend std::ostream& operator<<(std::ostream& Out, const Bracket<TYPE> &b) {
			  Out << "bracket = (" << b.a << ", " << b.b << ", " << b.c << ")" << "\t(";
			  Out << "evaluations=" << b.evaluations << ", ";
			  Out << "Status=";
			  if (b.status) {
				   Out << "Completed";
			  } else {
				   Out << "Terminated abnormally";
			  }
			  Out << ")\n";
			  return Out;
		 }
	};
	
	
	
	//////////////////// Brent ////////////////////
	
	// Optimize function of one variable
	// Finds minimum of function given a bracketing triplet (a, b, c) (see numerical recipes)
	// Returns minimum x value and optionally the corresponding function value
	// NOTE: The tolerance should be no smaller than the square root of the machines precision of double
	template <typename TYPE>
	class Brent {
	private:
		 // Constants
		 double minTol;
		 double goldenRatio;
		 double tinyNumber;
	
		 int maxIterations;
		 double xtol, ftol;
	public:
		 TYPE xmin;                 // minimal x-value 
		 double fmin;               // corresponding function value
		 int evaluations;
	
		 bool status;               // error status
	
		 // Default constructor
		 Brent() {};
	
		 // Initialization
		 void init(double xtol, double ftol, int maxIterations) {
			  // Constants
			  this->minTol = 1.0e-11;
			  this->goldenRatio = 0.3819660;
			  this->tinyNumber = 1e-25;
			  
	
			  evaluations = 0;
			  this->maxIterations = maxIterations;
			  this->xtol = xtol;
			  this->ftol = ftol;
			  
			  this->status = true;
			  
			  std::numeric_limits<double> l;
			  if (this->xtol<0) {
				   this->xtol = sqrt(l.epsilon());
			  }
			  if (this->ftol<0) {
				   this->ftol = sqrt(l.epsilon());
			  }
		 }
	
		 // Optimize - given point and direction - for functions already wrapped
		 template <typename FUNCTIONTYPE>
		 void optimize(Function<TYPE, FUNCTIONTYPE> &func, TYPE point, TYPE direction, double xtol=-1, double ftol=-1, int maxIterations = 500) {
			  compute(func, Bracket<TYPE>(func, point, point+(direction/len(direction))), xtol, maxIterations);
		 }
	
		 template <typename FUNCTIONTYPE>
		 void optimize(FUNCTIONTYPE func, Bracket<TYPE> bracket, double xtol=-1, double ftol=-1, int maxIterations = 500) {
			  ValueFunction<TYPE, FUNCTIONTYPE> f(func);
			  compute(f, bracket, xtol, ftol, maxIterations);          
		 }
	
		 // Optimize - default bracketing
		 template <typename FUNCTIONTYPE>
		 void optimize(FUNCTIONTYPE func, double xtol=-1, double ftol=-1, int maxIterations = 500) {
			  ValueFunction<TYPE, FUNCTIONTYPE> f(func);
			  compute(f, Bracket<TYPE>(f), xtol, ftol, maxIterations);
		 }
	
		 // Compute - main method
		 template <typename FUNCTIONTYPE>
		 void compute(Function<TYPE, FUNCTIONTYPE> &func, Bracket<TYPE> bracket, double xtol=-1, double ftol=-1, int maxIterations = 500) {
	//           std::cout << "********************************Brent\n";
			  double a, b;                // Minimum is bracketed between a and b
			  double x;                   // x-value corresponding to least function value so far
			  double w;                   // x-value corresponding to next-to-least function value so far
			  double v;                   // Previous value of w
			  double u;                   // x-value at which the function has been evaluated most recently
			  double xm;                  // midpoint between a and b (function is not evaluated here)
			  double fx, fw, fv, fu;      // function values
			  double xtol1, xtol2;        // Tolerance values
			  double d = 0.0;             // Step-size
			  double e = 0.0;             // Step-size of step i-2
	
			  init(xtol, ftol, maxIterations);
	
			  if (bracket.status == false) {
				   status=false;
				   // std::cerr << "Unconstrained bracket - no optimization done\n";
				   return;
			  }
	
			  TYPE unitVector = ((bracket.b-bracket.a)/(len(bracket.b-bracket.a)));
			  TYPE origin = bracket.a;
			  
			  // Ensure that a and b are in ascending order
			  if ((bracket.a-origin) < (bracket.c-origin)) {
				   a = len(bracket.a-origin);
				   b = len(bracket.c-origin);
			  } else {
				   a = len(bracket.c-origin);
				   b = len(bracket.a-origin);
			  }
	
	//           std::cout << bracket;
			  
			  // Initialize
			  x=w=v=len(bracket.b-origin);
			  fw=fv=fx=func(origin + unitVector*x);  evaluations++;
	
	
			  int i;
			  for (i=0; i<maxIterations; i++) {
	
				   xtol1 = this->xtol*fabs(x) + minTol;
				   xtol2 = 2*xtol1;
	
				   xm = 0.5 * (a+b);
	
				   // Termination test (convergence)
				   if (fabs(x-xm) <= (xtol2 - 0.5*(b-a))) {
						xmin = origin + unitVector * x;
						fmin = fx;
						return;
				   }
	
				   
				   if (fabs(e) <= xtol1) {
						// Golden step
						if (x >= xm) {
							 e = a - x;
						} else {
							 e = b - x;
						}
						d = goldenRatio * e;
				   } else {
						// Parabolic fit
						double r = (x-w)*(fx-fv);
						double q = (x-v)*(fx-fw);
						double p = (x-v)*q - (x-w)*r;
						q = 2.0 * (q-r);
						if (q > 0.0) {
							 p = -p;
						}
						q = fabs(q);
						double etemp = e;
						e = d;
						if ((fabs(p) < fabs(0.5*q*etemp)) && (p > q*(a-x)) && (p < q*(b-x))) {
							 // Use parabolic step
							 d = p/q;
							 u = x + d;
							 if ((u-a) < xtol2 or (b-u) < xtol2) {
								  if ((xm-x) >= 0) {
									   d = xtol1;
								  } else {
									   d = -xtol1;
								  }
							 }
						} else {
							 // Use golden step 
							 if (x >= xm) {
								  e = (a-x);
							 } else {
								  e = (b-x);
							 }
							 d = goldenRatio * e;
						}
				   }
	
				   // Use xtol1 if d is too small
				   if (fabs(d) >= this->xtol) {
						u = x + d;
				   } else {
						if (d >= 0) {
							 u = x + xtol1;
						} else {
							 u = x - xtol1;
						}
				   }
			  
				   fu = func(origin + unitVector * u); evaluations++;         // The function evaluation
	
				   if (fu <= fx) {
						if (u >= x) {
							 a = x;
						} else {
							 b = x;
						}
						v = w;
						w = x;
						x = u;
						fv = fw;
						fw = fx;
						fx = fu;
				   } else {
						if (u < x) {
							 a = u;
						} else {
							 b = u;
						}
						if ((fu <= fw) || (w==x)) {
							 v = w;
							 w = u;
							 fv = fw;
							 fw = fu;
						} else if (fu <= fv || v==x || v==w ) {
							 v = u;
							 fv = fu;
						}
				   }
			  }
	
			  xmin = origin + unitVector * x;
			  fmin = fx;
			  return;
		 }
			  
		 
		 // Output results
		 friend std::ostream& operator<<(std::ostream& Out, const Brent &b) {
		   Out << "(";
			  Out << "x=" << b.xmin << ", ";
			  Out << "f(x)=" << b.fmin << ", ";
			  Out << "evaluations=" << b.evaluations << ", ";
			  Out << "Status=";
			  if (b.status) {
				   Out << "Completed";
			  } else {
				   Out << "Terminated abnormally";
			  }
			  Out << ")\n";
			  return Out;
		 }
	};
	
	
	
	//////////////////// DirectionSet ////////////////////
	
	// Direction set used by powell
	class DirectionSet {
	public:
		 Vector_nD *data;
		 int size;
	
		 // Create direction set from unit vectors
		 DirectionSet(int size) {
			  init(size);
			  initAsUnitVector_nDs(size);
		 }
	
		 // Create direction set from vector list
		 DirectionSet(int size, Vector_nD *vectorList) {
			  init(size);
			  for (int i=0; i<size; i++) {
				   data[i] = vectorList[i];
			  }
		 }
	
		 // Copy constructor
		 DirectionSet(const DirectionSet &d) {
			  init(d.size);
			  for (int i=0; i<this->size; i++) {
				   this->data[i] = d.data[i];
			  }
		 }
	
		 // Destructor
		 ~DirectionSet() {
			  delete[] this->data;
		 }
	
		 // Initialization
		 void init(int size) {
			  this->size = size;
			  data = new Vector_nD[size];
		 }
			  
		 // Create unit vectors
		 void initAsUnitVector_nDs(int size) {
			  for (int i=0; i<size; i++) {
				   Vector_nD v(size);
				   for (int j=0; j<size; j++) {
						if (i==j) {
							 v[j] = 1;
						} else {
							 v[j] = 0;
						}
				   }
				   data[i] = v;
			  }
		 }
		 
		 //Overload indexing operator
		 Vector_nD operator[](const int index) const {
			  return data[index];
		 }
	
		 // Overload indexing operator (non-const)
		 Vector_nD& operator[](const int index) {
			  return data[index];
		 }
		 
		 // Output
		 friend std::ostream& operator<<(std::ostream& Out, const DirectionSet &d) {
			  for (int i=0; i<d.size; i++){
				   std::cout << d.data[i] << "\n";
			  }
			  return Out;
		 }
	};
	
	
	
	//////////////////// Powell ////////////////////
	
	// Modified Powell's method
	// Minimizes function of n variables given a starting point and an initial direction matrix (defaults to unit vectors)
	class Powell {
	private:
		 Brent<Vector_nD> brentOptimizer;
		 double tinyNumber;
	
	public:
		 Vector_nD xmin;                  // minimal x-value 
		 double fmin;                 // corresponding function value
		 double deltaf;             // Change in function value due to optimization
		 int evaluations;
		 int iterations;
	
		 double xtol;               // x-value tolerance
		 double ftol;               // Fractional tolerance in the function value (used for termination)
		 
		 bool status;               // error status
	
		 // Default constructor
		 Powell() {
			  this->init();
		 };
	
		 void init() {
			  this->tinyNumber = 1e-25;
		 }     
	
		 // Optimize - Main method
		 template <typename FUNCTIONTYPE>
		 void compute(FUNCTIONTYPE &func, Vector_nD startingPoint, DirectionSet directions,
					  int maxIterations=1000, double xtol=-1, double ftol=-1) {
	//           std::cout << "********************************Powell\n";
	//           std::cout.flush();
			  evaluations = 0;
			  this->status = true;
	
			  xmin.clear();
			  
			  // Initialize tolerances
			  this->xtol = xtol;
			  this->ftol = ftol;
			  std::numeric_limits<double> l;
			  if (this->xtol<0) {
				   this->xtol = sqrt(l.epsilon());
			  }
			  if (this->ftol<0) {
				   this->ftol = sqrt(l.epsilon());
			  }
	
			  assert(directions.size == startingPoint.size);
	
			  xmin = startingPoint;
			  fmin = func(xmin); evaluations++;
			  deltaf = fmin;
	
			  for (iterations=1; ; iterations++) {
				   Vector_nD xOld = xmin;                  // Keep track of previous values 
				   double fxOld = fmin;
	
				   int iLargest = -1;                // index corresponding to direction with largest decrease
				   double deltaFLargest = 0.0;       // function decrease corresponding to iLargest
	
				   // Call brent for all directions
				   for (int i=0; i<directions.size; i++) {
						Vector_nD direction = directions[i];
						brentOptimizer.optimize(func, xmin, direction, xtol);
						if (brentOptimizer.status == true) {
							 if ((fmin - brentOptimizer.fmin) > deltaFLargest) {
								  iLargest = i;
								  deltaFLargest = fmin - brentOptimizer.fmin;
							 }
							 xmin = brentOptimizer.xmin;     // Update x, fx and evaluations
							 fmin = brentOptimizer.fmin;
						}
						evaluations += brentOptimizer.evaluations;
				   }
	
				   // Termination criterium
				   if (2.0 * (fxOld-fmin) < this->ftol*(fabs(fxOld)+fabs(fmin))+tinyNumber) {
						break;
				   }
	
				   // Maximum number of iterations exceeded
				   if (iterations >= maxIterations) {
						status = false;
						break;
				   }
	
				   // Find extrapolated point based on Pn-P0 direction
				   Vector_nD extrapolated = xmin*2.0 - xOld;
				   double fExtrapolated = func(extrapolated); evaluations++;
				   Vector_nD direction = xmin-xOld;
	
				   // Check whether it is better to keep old direction set
				   if (fExtrapolated < fxOld) {
						if ((2.0*(fxOld-2.0*fmin+fExtrapolated) *
							 (fxOld - fmin - deltaFLargest) * (fxOld - fmin - deltaFLargest)) < 
							(deltaFLargest * (fxOld - fExtrapolated)*(fxOld - fExtrapolated))) {
							 
							 brentOptimizer.optimize(func, xmin, direction, xtol); // Call brent on new direction
							 if (brentOptimizer.status == true) {
								  xmin = brentOptimizer.xmin; // This is different than in NR - They save the old value of x before doing this step (why?)
								  fmin = brentOptimizer.fmin;
								  directions[iLargest] = directions[directions.size-1]; // Add new direction to direction set
								  directions[directions.size-1] = direction;
							 }
							 evaluations += brentOptimizer.evaluations;
						}
				   }
			  }
			  deltaf = deltaf - fmin;
		 }
	
		 
		 template <typename FUNCTIONTYPE>
		 void optimize(FUNCTIONTYPE func, Vector_nD startingPoint,
					   int maxIterations=1000, double xtol=-1, double ftol=-1) {
			  ValueFunction<Vector_nD, FUNCTIONTYPE> f(func);
			  compute(f, startingPoint, DirectionSet(startingPoint.size), maxIterations, xtol, ftol);
		 }
		 
		 template <typename FUNCTIONTYPE>
		 void optimize(FUNCTIONTYPE func, std::vector<double *> startingPoint, DirectionSet directions,
					   int maxIterations=1000, double xtol=-1, double ftol=-1) {
			  std::cout << "Starting Powell optimization\n";
			  std::cout.flush();
			  double **xPointers = new double*[startingPoint.size()];
			  for (unsigned int i=0; i<startingPoint.size(); i++) {
				   xPointers[i] = startingPoint[i];
			  }
			  PointerFunction<Vector_nD, FUNCTIONTYPE> f(func, startingPoint.size(), xPointers);
			  compute(f, Vector_nD(startingPoint), directions, maxIterations, xtol, ftol);
			  for (unsigned int i=0; i<startingPoint.size(); i++) {
				   *xPointers[i] = xmin[i];
			  }
			  delete[] xPointers;
		 }
	
		 template <typename FUNCTIONTYPE>
		 void optimize(FUNCTIONTYPE func, std::vector<double *> startingPoint,
					   int maxIterations=1000, double xtol=-1, double ftol=-1) {
		  optimize(func, startingPoint, DirectionSet(2), maxIterations, xtol, ftol);
		 }
		 
		 
		 // Output results
		 friend std::ostream& operator<<(std::ostream& Out, const Powell &p) {
			  Out << "(";
			  Out << "xmin=" << p.xmin << ", ";
			  Out << "f(xmin)=" << p.fmin << ", ";
			  Out << "delta_f=" << p.deltaf << ", ";
			  Out << "iterations=" << p.iterations << ", ";
			  Out << "evaluations=" << p.evaluations << ", ";
			  Out << "Status=";
			  if (p.status) {
				   Out << "Completed";
			  } else {
				   Out << "Terminated abnormally";
			  }
			  Out << ")\n";
			  return Out;
		 }
	};
}

#endif
