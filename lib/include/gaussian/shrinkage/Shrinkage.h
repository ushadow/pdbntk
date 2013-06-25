/*
 * Shrinkage.h
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

#ifndef SHRINKAGE_H_
#define SHRINKAGE_H_

#include "ShrinkageHelpers.h"

// GENERATE_W TIENE QUE CALCULARSE EN FUNCION DE S_MLE, NO DE LOS DATOS!!!!

namespace mocapy
{
class Shrinkage{

        protected:
                int N;
                int P;
                MDArray<double> data;
                MDArray<double> S_MLE;
                MDArray<double> wk;
                MDArray<double> StDM;
                MDArray<double> wk_StDM;
		MDArray<double> compute_f();
	        void generate_wk(MDArray<double>,int);
                void generate_MLE();
        public:
                Shrinkage(MDArray<double>);
		Shrinkage(MDArray<double>, MDArray<double>);
		MDArray<double> combine(MDArray<double>,MDArray<double>,double);
                MDArray<double> get_data();
                MDArray<double> get_S_MLE();
		void set_S_MLE(MDArray<double>);
                MDArray<double> get_wk();
                MDArray<double> get_StDM();
                MDArray<double> get_wk_StDM();
                double get_lambda_cor();
                double get_lambda_var();
                double get_lambda_A();
                double get_lambda_B();
                double get_lambda_C();
                double get_lambda_D();
                MDArray<double> get_target_A();
                MDArray<double> get_target_B();
                MDArray<double> get_target_C();
                MDArray<double> get_target_D();
};

}

#endif // SHRINKAGE_H_
