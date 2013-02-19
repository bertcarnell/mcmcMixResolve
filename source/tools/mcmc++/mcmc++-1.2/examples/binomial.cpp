// standard includes
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
// Boost includes
#include <boost/any.hpp>
#include <boost/format.hpp>
#include <boost/spirit/core.hpp>
#include <boost/spirit/utility/lists.hpp> 
// MCMC++ includes
#include <mcmc++/Density.h>
#include <mcmc++/MCMC.h>
#include <mcmc++/intervals.h>
#include <mcmc++/statistics.h>
#include <mcmc++/util.h>

using std::cerr;
using std::cout;
using std::endl;
using std::vector;
using boost::any;
using boost::any_cast;
using namespace boost;
using namespace boost::spirit;
using Density::dbeta;
using Density::dbinom;

namespace {
  lot rng_(lot::RAN_MT, lot::ZERO);
}

#if defined(METRO)
typedef MetroStep StepType;
#else
#if defined(DIRECT)
typedef FunctionStep StepType;
#else
typedef SliceStep StepType;
#endif  // DIRECT
#endif  // METRO

class BinomialModel;

class p : public Parameter {
public:
  p(BinomialModel* bin)
    : Parameter("p"), bin_(bin)
  {
    Assign(0.5);
  }

  double llike(const double p0) const;

  double propose(const double current) const {
    return proposeBeta(current, 5, Util::dbl_eps);
  }

  double lQ(const double x, const double y) const {
    return logQBeta(x, y, 5, Util::dbl_eps);
  }

  const double Function(const bool doCalc = true) const;

private:
  BinomialModel* bin_;

};


class BinomialModel : public Model {
public:
  BinomialModel(const int nBurnin, const int nSample, const int thin, 
                const int n, const int k)
    : Model(nBurnin, nSample, thin, true), n_(n), k_(k)
  {
    step_.push_back(new StepType(new p(this)));
    step_[0]->SetBounds(Util::dbl_min, 1.0 - Util::dbl_eps);
  }

  void Report(std::ostream& outf) {
    double mean = static_cast<double>(k_+1)/static_cast<double>(n_+2);
    double variance = mean*(1.0 - mean)/(n_+2+1);
#if defined(METRO)
    outf << "Metropolis-Hastings sampler...";
#else
#if defined(DIRECT)
    outf << "Direct sampler...";
#else
    outf << "Slice sampler...";
#endif  // DIRECT
#endif  // METRO
    outf << endl << format("Predicted: %|10| %|10|")
      % mean % sqrt(variance) << endl;
    Model::Report(outf);
  }

  inline int K(void) const {
    return k_;
  }

  inline int N(void) const {
    return n_;
  }

  double Llike(const SampleVector& p0) const {
    double p = any_cast<double>(p0[0]);
    return dbinom(k_, n_, p, true);
  }

private:
  int n_, k_;
};

double 
p::llike(const double p) const {
  return dbinom(bin_->K(), bin_->N(), p, true);
}

const double
p::Function(const bool doCalc) const {
  return rng_.beta(bin_->K() + 1, bin_->N() - bin_->K() + 1);
}

int main(int ac, char** av) {
  int k;
  parse(av[1], int_p[assign(k)], space_p);
  int n;
  parse(av[2], int_p[assign(n)], space_p);

  int nBurnin;
  parse(av[3], int_p[assign(nBurnin)], space_p);
  int nSample;
  parse(av[4], int_p[assign(nSample)], space_p);
  int thin;
  parse(av[5], int_p[assign(thin)], space_p);

  BinomialModel model(nBurnin, nSample, thin, n, k);
  model.Simulation(cerr, false);
  model.Report(cout);
}
