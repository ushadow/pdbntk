#include "cond_prob_dist.h"

namespace pdbntk {
void CondProbDist::setRandomGen(mocapy::RandomGen* rg) {
  density_->setRandomGen(rg);
}

void CondProbDist::update_ess(std::vector<double> & ptv) {
	ess_->add_ptv(ptv);
}

void CondProbDist::save_ess() {
	ess_->save_ess(weight_);
}

std::vector<mocapy::MDArray<double> > CondProbDist::get_ess() {
	return ess_->get_array();
}

void CondProbDist::do_M_step(std::vector<mocapy::MDArray<double> > & new_ess) {
  // Update the parameters
  density_->estimate(new_ess);
  // Get ready for new cycle
  ess_->clear();
}

uint CondProbDist::get_output_size() {
  return density_->get_output_size();
}

std::vector<mocapy::MDArray<double> > CondProbDist::get_parameters() {
  return density_->get_parameters();
}

}
