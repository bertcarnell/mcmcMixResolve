// standard includes
#include <iostream>
#include <cmath>
#include <sstream>
#include <vector>
// Boost includes
#include <boost/format.hpp>
#include <boost/spirit/core.hpp>
#include <boost/spirit/utility/lists.hpp> 
// MCMC++ includes
#include <mcmc++/DataTable.h>
#include <mcmc++/Density.h>
#include <mcmc++/MCMC.h>
#include <mcmc++/intervals.h>
#include <mcmc++/statistics.h>
#include <mcmc++/util.h>

// using declarations
using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::vector;
// namespaces
using namespace boost;
using namespace boost::spirit;

namespace {
  boost::format formatter_("%|12| %|12| %|12| %|12| %|12| %|12|");
}

class MultinomialModel;

class lambda : public Parameter {
public:
  lambda(MultinomialModel* multi, const int index)
    : Parameter("lambda"), multi_(multi), idx_(index)
  {
    std::ostringstream ost;
    ost << "lambda[" << idx_ + 1 << "]";
    SetLabel(ost.str());
    Assign(1);
  }

  double llike(const double lambda0) const;
  double lPrior(const double lambda0) const;
  
private:
  MultinomialModel* multi_;
  const int idx_;

};

class p : public Parameter {
public:
  p(MultinomialModel* multi, const int idx)
    : Parameter("p"), multi_(multi), idx_(idx)
  {
    std::ostringstream ost;
    ost << "p[" << idx_ + 1 << "]";
    SetLabel(ost.str());
  }
  const double Function(const bool doCalc = true) const;

private:
  MultinomialModel* multi_;
  const int idx_;

};

class MultinomialModel : public Model {
public:
  MultinomialModel(const int nBurnin, const int nSample, const int thin,
                   vector<int>& x)
    : Model(nBurnin, nSample, thin), n_(x), size_(n_.size())
  {
    for (int i = 0; i < size_; ++i) {
      step_.push_back(new SliceStep(new lambda(this, i)));
      step_[i]->SetBounds(Util::dbl_min, Util::dbl_max);
    }
    for (int i = 0; i < size_; ++i) {
      step_.push_back(new FunctionStep(new p(this, i)));
    }
  }

  int N(const int idx) const {
    return n_[idx];
  }

  double Lambda(const int idx) const {
    return step_[idx]->Value();
  }

  int Size(void) const {
    return size_;
  }

  double CumLambda(void) const {
    double sum = 0.0;
    for (int i = 0; i < size_; ++i) {
      sum += Lambda(i);
    }
    return sum;
  }

  SampleVector Parameters(void) const {
    SampleVector x;
    // only interested in p, not lambda
    ModelSteps::const_iterator i = step_.begin() + size_;
    ModelSteps::const_iterator iEnd = step_.end();
    for (; i != iEnd; ++i) {
      x.push_back((*i)->Value());
    }
    return x;
  }

  void Report(std::ostream& outf) {
    std::vector<double> a(size_);
    double a0 = 0.0;
    for (int i = 0; i < size_; ++i) {
      a[i] = n_[i] + 1;
      a0  += n_[i] + 1;
    }
    std::vector<double> mu(size_);
    std::vector<double> sd(size_);
    for (int i = 0; i < size_; ++i) {
      mu[i] = a[i]/a0;
      sd[i] = sqrt(a[i]*(a0 - a[i])/(Util::sqr(a0)*(a0+1)));
    }
    outf.precision(6);
    outf << endl << "Predicted: "
         << endl << "  Mean:    " << mu
         << endl << "  s.d.:    " << sd << endl;
    Model::Report(outf);
  }

  void Summarize(const int i, std::ostream& outf) {
    int n = nSample_/nThin_;
    vector<double> x(n);
    for (int k = 0; k < n; ++k) {
      x[k] = any_cast<double>(results_[k][i]);
    }
    SimpleStatistic xStat(x);
    // add size_ because of offset for p[]
    outf << formatter_ % Label(i+size_)
      % xStat.Mean() % xStat.StdDev() % quantile(x, 0.025) 
      % quantile(x, 0.5) % quantile(x, 0.975) << endl;
  }

private:
  vector<int> n_;
  int size_;

};

double 
lambda::llike(const double lambda0) const {
  return Density::dpois(multi_->N(idx_), lambda0, true);
}

double 
lambda::lPrior(const double lambda0) const {
  return Density::dgamma(lambda0, 1.0, 1.0, true);
}

const double
p::Function(const bool doCalc) const {
  double k = multi_->Lambda(idx_);
  double n = multi_->CumLambda();
  return k/n;
}

int main(int ac, char** av) {
  string fileName = av[1];
  DataTable<int> data(false, false);
  data.Read(fileName);
  vector<int> x = data.ColumnVector(0);

  int nBurnin;
  parse(av[2], int_p[assign(nBurnin)], space_p);
  int nSample;
  parse(av[3], int_p[assign(nSample)], space_p);
  int thin;
  parse(av[4], int_p[assign(thin)], space_p);

  MultinomialModel model(nBurnin, nSample, thin, x);
  model.Simulation(cerr, true);
  model.Report(cout);
}
