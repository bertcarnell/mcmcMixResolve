// standard headers
#include <sstream>
#include <string>
#include <vector>
// boost headers
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
// local headers
#include "mcmc++/MCMC.h"
#include "mcmc++/lot.h"

using namespace boost::unit_test_framework;

namespace {

  std::string result = "       Label         5.12         0.23         4.66         5.12         5.58";
  double precision = 1.0e-13;

}

class DummyModel : public Model {
public:
  DummyModel(void) 
    : Model(1, 1, 1)
  {}
  std::string Format(void) {
    std::ostringstream ost;
    ost << *summaryFormat_ % "Label"
      % 5.12 % 0.23 % 4.66 % 5.12 % 5.58;
    return ost.str();
  }
};

void
formatTest(void) {
  DummyModel model;
  std::string ost = model.Format();
  BOOST_CHECK_EQUAL(ost, result);
}

void
dirchTest(void) {
  int nElem = 5;
  // N.B.: 1.0e-1 and 1.0e4 come from MCMC.cpp, be sure to keep in sync
  double small = 0.5*(1.0e-1);
  double big = 2.0*1.0e4;
  std::vector<double> ok(nElem, 10.0);
  std::vector<double> oneSmall(ok);
  oneSmall[2] = small;
  std::vector<double> oneBig(ok);
  oneBig[3] = big;
  std::vector<double> oneBigoneSmall(oneSmall);
  oneBigoneSmall[3] = oneBig[3];
  for (int i = 0; i < nElem; ++i) {
    BOOST_CHECK_EQUAL(ok[i], 10.0);
  }
  for (int i = 0; i < nElem; ++i) {
    BOOST_CHECK_EQUAL(oneSmall[i], (i == 2) ? small : 10.0);
  }
  for (int i = 0; i < nElem; ++i) {
    BOOST_CHECK_EQUAL(oneBig[i], (i == 3) ? big : 10.0);
  }
  for (int i = 0; i < nElem; ++i) {
    BOOST_CHECK_EQUAL(oneBigoneSmall[i], (i == 2) ? small : 
                      (i == 3) ? big : 10.0);
  }
}

void
safeFreqTest(void) {
  using std::vector;
  double x = safeFreq(0.0);
  BOOST_CHECK_EQUAL(x, MCMC_ZERO_FREQ);
  x = safeFreq(1.0);
  BOOST_CHECK_EQUAL(x, 1.0-MCMC_ZERO_FREQ);
  x = safeFreq(0.5);
  BOOST_CHECK_EQUAL(x, 0.5);
  x = safeFreq(0.3, 0.001);
  BOOST_CHECK_CLOSE(x, 0.3004, precision);
  BOOST_CHECK_CLOSE(invSafeFreq(x, 0.001), 0.3, precision);
  vector<double> y(3);
  y[0] = 0.1;
  y[1] = 0.2;
  y[2] = 0.7;
  y = safeFreq(y, 0.0001);
  BOOST_CHECK_CLOSE(y[0], 0.10007, precision);
  BOOST_CHECK_CLOSE(y[1], 0.20004, precision);
  BOOST_CHECK_CLOSE(y[2], 0.69989, precision);
  y = invSafeFreq(y, 0.0001);
  BOOST_CHECK_CLOSE(y[0], 0.1, precision);
  BOOST_CHECK_CLOSE(y[1], 0.2, precision);
  BOOST_CHECK_CLOSE(y[2], 0.7, precision);
}

void
safeVectorTest(void) {
  int nElem = 4; // must be a power of 2
  std::vector<double> x(nElem, 0.0);
  x = safeFreq(x, MCMC_ZERO_FREQ);
  for (int i = 0; i < nElem; ++i) {
    BOOST_CHECK_EQUAL(x[i], MCMC_ZERO_FREQ);
  }
  std::vector<double> y(nElem, 1.0);
  y = safeFreq(y, MCMC_ZERO_FREQ);
  for (int i = 0; i < nElem; ++i) {
    BOOST_CHECK_EQUAL(y[i], 1.0 - (nElem-1)*MCMC_ZERO_FREQ);
  }
  std::vector<double> z(nElem, 1.0/nElem);
  z = safeFreq(z, MCMC_ZERO_FREQ);
  for (int i = 0; i < nElem; ++i) {
    BOOST_CHECK_EQUAL(z[i], 1.0/nElem);
  }
}

void
proposeTest(void) {
  lot& rng = GetRNG();
  rng.set_seed(1234L);
  BOOST_CHECK_CLOSE(proposeBeta(0.5, 10.0, Util::dbl_eps),
                    0.3443288124095045, precision);
  rng.set_seed(1234L);
  BOOST_CHECK_CLOSE(proposeBeta(0.5, 35.0, 1/(1+MCMC::MaxBetaPar)),
                    0.4147906755691068, precision);
  std::vector<double> half(2, 1.0/2.0);
  rng.set_seed(1234L);
  std::vector<double> p = proposeDirch(half, 15.0, Util::dbl_eps);
  BOOST_CHECK_CLOSE(p[0], 0.3687066232607837, precision);
  std::vector<double> third(3, 1.0/3.0);
  rng.set_seed(1234L);
  p = proposeDirch(third, 15.0, Util::dbl_eps);
  BOOST_CHECK_CLOSE(p[0], 0.1956208290497083, precision);
  BOOST_CHECK_CLOSE(p[1], 0.3824043236891421, precision);
  BOOST_CHECK_CLOSE(p[2], 0.4219748472611496, precision);
  rng.set_seed(1234L);
  p = proposeDirch(third, 10.0, Util::dbl_eps);
  BOOST_CHECK_CLOSE(p[0], 0.1680965225690985, precision);
  BOOST_CHECK_CLOSE(p[1], 0.3913248147021047, precision);
  BOOST_CHECK_CLOSE(p[2], 0.4405786627287968, precision);
}

void
lQTest(void) {
  BOOST_CHECK_CLOSE(logQBeta(0.3, 0.5, 10.0, 0.001), 
                    0.2000708713918532, precision);
  BOOST_CHECK_CLOSE(logQBeta(0.5, 0.3, 10.0, 0.001), 
                    -0.01574835696813909, 1.0e-11);
  BOOST_CHECK_CLOSE(logQBeta(0.2, 0.3, 10.0, 0.001), 
                    0.9701782010663446, precision);
  BOOST_CHECK_CLOSE(logQBeta(0.3, 0.2, 10.0, 0.001), 
                    0.5786386775517017, precision);
}

test_suite*
init_unit_test_suite(int ac, char** av) {
  test_suite* test = BOOST_TEST_SUITE("MCMC test suite");
  test->add( BOOST_TEST_CASE(&formatTest) );
  test->add( BOOST_TEST_CASE(&dirchTest) );
  test->add( BOOST_TEST_CASE(&safeFreqTest) );
  test->add( BOOST_TEST_CASE(&safeVectorTest) );
  test->add( BOOST_TEST_CASE(&proposeTest) );
  test->add( BOOST_TEST_CASE(&lQTest) ) ;
  return test;
}



