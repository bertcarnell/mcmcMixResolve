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
#include <mcmc++/DataTable.h>
#include <mcmc++/Density.h>
#include <mcmc++/MCMC.h>
#include <mcmc++/intervals.h>
#include <mcmc++/statistics.h>
#include <mcmc++/util.h>

using std::cerr;
using std::cin;
using std::cout;
using std::endl;
using std::string;
using std::vector;
using boost::any;
using boost::any_cast;
using namespace boost;
using namespace boost::spirit;
using Density::dgamma;
using Density::dinvgamma;
using Density::dnorm;

typedef SliceStep StepType;

class NormalModel;

class mean : public Parameter {
public:
  mean(NormalModel* norm)
    : Parameter("mean"), norm_(norm)
  {
    Assign(0.0);
  }
  double llike(const double mu) const;
  double lPrior(const double mu0) const;

private:
  NormalModel* norm_;
};


class variance : public Parameter {
public:
  variance(NormalModel* norm) 
    : Parameter("variance"), norm_(norm)
  {
    Assign(1.0);
  }
  double llike(const double var) const;
  double lPrior(const double var) const;

private:
  NormalModel* norm_;
};

class precision : public Parameter {
public:
  precision(NormalModel* norm)
    : Parameter("precision"), norm_(norm)
  {}
  const double Function(const bool doCalc = true) const;

private:
  NormalModel* norm_;
};


class NormalModel : public Model {
public:
  NormalModel(const int nBurnin, const int nSample, const int thin, 
              const vector<double>& x, const bool usePrecision,
              const bool useMedian)
    : Model(nBurnin, nSample, thin, true, useMedian), x_(x), 
      usePrecision_(usePrecision)
  {
    step_.push_back(new StepType(new mean(this)));
    step_.push_back(new StepType(new variance(this)));
    step_.push_back(new FunctionStep(new precision(this)));
    // set lower bound on variance
    step_[1]->SetBounds(Util::dbl_min, Util::dbl_max);
  }

  double Mean() const {
    return step_[0]->Value();
  }

  double Var() const {
    return step_[1]->Value();
  }

  double X(const int index) const {
    return x_[index];
  }

  double llike(const double mu, const double sd) const {
    double llike = 0;
    unsigned nElem = x_.size();
    for (unsigned i = 0; i < nElem; ++i) {
      llike += dnorm(x_[i], mu, sd, true);
    }
    return llike;
  }

  double Llike(const SampleVector& p) const {
    double a;
    if (usePrecision_) {
      a = 1/any_cast<double>(p[2]);  //  p[2] == precision
    } else {
      a = any_cast<double>(p[1]);    //  p[1] == variance
    }
    double mu = any_cast<double>(p[0]);
    double sd = sqrt(a);
    return llike(mu, sd);
  }

private:
  const vector<double>& x_;
  bool usePrecision_;

};

double 
mean::llike(const double mu) const {
  double sd = sqrt(norm_->Var());
  return norm_->llike(mu, sd);
}

double
mean::lPrior(const double mu) const {
  return dnorm(mu, 0.0, sqrt(1000.0), true);
}

double
variance::llike(const double var) const {
  double mu = norm_->Mean();
  double sd = sqrt(var);
  return norm_->llike(mu, sd);
}

double
variance::lPrior(const double var) const {
  return dnorm(var, 0, 10, true);
}

const double
precision::Function(const bool doCalc) const {
  return 1.0/norm_->Var();
}



int main(int ac, char** av) {
  bool usePrecision = true;
  bool useMedian = false;
  while (av[1][0] == '-') {
    string arg = av[1];
    if (arg == "--prec") {
      usePrecision = true;
    } else if (arg == "--var") {
      usePrecision = false;
    } else if (arg == "--med") {
      useMedian = true;
    } else if (arg == "--mean") {
      useMedian = false;
    } else {
      cerr << "Unrecognized option: " << arg << endl;
      exit(1);
    }
    --ac;
    ++av;
  }

  DataTable<double> data(false, false);
  string fileName = av[1];
  data.Read(fileName);
  vector<double> x = data.ColumnVector(0);

  int nBurnin;
  parse(av[2], int_p[assign(nBurnin)], space_p);
  int nSample;
  parse(av[3], int_p[assign(nSample)], space_p);
  int thin;
  parse(av[4], int_p[assign(thin)], space_p);

  NormalModel model(nBurnin, nSample, thin, x, usePrecision, useMedian);
  model.Simulation(cerr, true);
  cout << "DIC based on " << (usePrecision ? "precision" : "variance")
       << " and " << (useMedian ? "median" : "mean") << endl;
  model.Report(cout);
}
