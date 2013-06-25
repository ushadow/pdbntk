// integrate.h --- Wrapper for fortran integration code
// Copyright (C) 2006-2008 Wouter Boomsma
//
// This file is part of BackboneDBN
//
// BackboneDBN is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// BackboneDBN is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public License
// along with BackboneDBN.  If not, see <http://www.gnu.org/licenses/>.


#ifndef INTEGRATE_H
#define INTEGRATE_H

#include <cstdlib>

extern "C" int dqags_(double (*f)(double*), double *a, double *b, double *epsabs, double *epsrel, double *result, double *abserr, int *neval, int *ier, int *limit, int *lenw, int *last, int *iwork, double *work);

extern "C" double dgamma_(double *x);

// Computes a definite integral
// Integrate func from a to b (possibly infinite interval) using a technique
// from the Fortran library QUADPACK
double integrateQuad(double(*func)(double *), double a, double b, int *evaluations=NULL,
            double epsabs=1.49e-8, double epsrel=1.49e-8, int limit=50);


// Version of quad taking additional arguments of type double
double integrateQuad(double(*func)(double *, double *), double extraArguments[], double a, double b, 
                     int *evaluations=NULL, double epsabs=1.49e-8, double epsrel=1.49e-8, int limit=50);


#endif
