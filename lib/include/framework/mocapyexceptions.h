/*
 * MocapyExceptions.h
 *
 *  Copyright (C) 2008, Martin Paluszewski, The Bioinformatics Centre, University of Copenhagen.
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

#ifndef MOCAPYEXCEPTIONS_H_
#define MOCAPYEXCEPTIONS_H_

#include <exception>

namespace mocapy {

class MocapyExceptions: public std::exception {
public:
	MocapyExceptions(const char* new_msg) {
		msg = new_msg;
	}
	const char* what() const throw() {
		return msg;
	}

private:
	const char* msg;
};

}
#endif /* MOCAPYEXCEPTIONS_H_ */
