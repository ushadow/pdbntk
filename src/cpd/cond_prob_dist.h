#ifndef PDBNTK_COND_PROB_DIST_H_
#define PDBNTK_COND_PROB_DIST_H_

#include "framework/densitiesbase.h"

#include <boost/scoped_ptr.hpp>

namespace pdbntk {

enum eCPDType {
  DISCRETE, GAUSSIAN
};

/// Represents a conditional probability distribution (CPD). CPD consists two
/// parts: the density function and the estimated sufficient statistics.
class CondProbDist {
public:
  /// This object takes the ownership of \a ess and \a density.
	CondProbDist(eCPDType type, mocapy::ESSBase *ess, mocapy::DensitiesBase *density) 
      : cpd_type_(type), ess_(ess), density_(density) {};

/// Accessors
//@{
  uint node_size() const { return density_->get_node_size(); }
  mocapy::DensitiesBase* get_densities() { return density_.get(); }
//@}

	// Return a sample of the i'th slice
	void sample(uint i);

	// Add a new sample to the ESS
	void update_ess(std::vector<double> & ptv);

	// Flag end of sequence sampling
	void save_ess();

	std::vector<mocapy::MDArray<double> > get_ess();

	void do_M_step(std::vector<mocapy::MDArray<double> > & ess);

	void set_densities(mocapy::DensitiesBase *d) { density_.reset(d); }
	std::vector<mocapy::MDArray<double> > get_parameters();

	void setRandomGen(mocapy::RandomGen* rg);

	// Dimension
	uint get_output_size();
	eCPDType cpd_type() { return cpd_type_; }

  mocapy::ESSBase& get_ess_ref() { return *ess_; }

private:
  eCPDType cpd_type_;
  boost::scoped_ptr<mocapy::ESSBase> ess_;
  boost::scoped_ptr<mocapy::DensitiesBase> density_;
  double weight_;
};

}

#endif // PDBNTK_COND_PROB_DIST_H_ 
