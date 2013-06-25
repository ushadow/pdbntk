// functions.h --- interface to special functions
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


#ifndef BESSEL_H
#define BESSEL_H

/* extern "C" int dvi0_(int *M, double *X, double *F, double *WORK, int *IWORK, int *INFO); */
extern "C" double i0(double x);
extern "C" double i0e(double x);
extern "C" double i1(double x);
extern "C" double i1e(double x);
extern "C" double dgamma(double x);
extern "C" double lgam(double x);
extern "C" double zeta(double x, double q);

#endif
