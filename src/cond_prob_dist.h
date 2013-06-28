#ifndef PDBNTK_COND_PROB_DIST_H_
#define PDBNTK_COND_PROB_DIST_H_

#include "framework/densitiesbase.h"

#include <boost/scoped_ptr.hpp>

namespace pdbntk {

enum eCPDType {
     DISCRETE, GAUSSIAN, DIRICHLET, KENT, VONMISES, VONMISES2D, POISSON, MULTINOMIAL, BIPPO
};

/// Represents a conditional probability distribution (CPD). CPD consists two
/// parts: the density function and the estimated sufficient statistics.
class CondProbDist {
public:
	CondProbDist(eCPDType type) : cpd_type_(type) {};

	// Return a sample of the i'th slice
	void sample(uint i);

	// Add a new sample to the ESS
	void update_ess(std::vector<double> & ptv);

	// Flag end of sequence sampling
	void save_ess();

	std::vector<mocapy::MDArray<double> > get_ess();

	void do_M_step(std::vector<mocapy::MDArray<double> > & ess);

	void set_densities(mocapy::DensitiesBase *d) { densities_.reset(d); }
  mocapy::DensitiesBase* get_densities() { return densities_.get(); }
	std::vector<mocapy::MDArray<double> > get_parameters();

	void setRandomGen(mocapy::RandomGen* rg);

	// Dimension
	uint get_output_size();
	eCPDType cpd_type() { return cpd_type_; }

  mocapy::ESSBase& get_ess_ref() { return *ess_; }

private:
  eCPDType cpd_type_;
  boost::scoped_ptr<mocapy::ESSBase> ess_;
  boost::scoped_ptr<mocapy::DensitiesBase> densities_;
  double weight_;
};

}

#endif // PDBNTK_COND_PROB_DIST_H_ 
