/*
 * ShrinkageHelpers.h
 *
 *  Copyright (C) 2008, Pedro, The Bioinformatics Centre,
 *  University of Copenhagen.
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

#ifndef SHRINKAGEHELPERS_H_
#define SHRINKAGEHELPERS_H_

#include "../../framework/essbase.h"

namespace mocapy
{
double get_mean(MDArray<double>);
double get_var(MDArray<double>);
double unb_emp_covar(MDArray<double>,int,int);
double unb_emp_var(MDArray<double>,int);
double get_median(MDArray<double>);
MDArray<double> Sort(MDArray<double>);
void Swap(double&, double&);
int index_of_smallest(MDArray<double>, int, int);
double get_tot_mean(MDArray<double>);
MDArray<double> get_diag(MDArray<double>);
MDArray<double> get_off_diag(MDArray<double>);
MDArray<double> get_column(MDArray<double>, int);
double round_lambda(double);
MDArray<double> get_sub_means(MDArray<double>);
double w_var(MDArray<double>,int,int);
double w_covar(MDArray<double>,int,int,int,int);
MDArray<double> get_SDM(MDArray<double>);
MDArray<double> get_partial_correlation(MDArray<double>); 
MDArray<double> add_matrices(MDArray<double>,MDArray<double>);
MDArray<double> mult_num_to_matrix(MDArray<double>,double);
}

#endif // SHRINKAGEHELPERS_H_
