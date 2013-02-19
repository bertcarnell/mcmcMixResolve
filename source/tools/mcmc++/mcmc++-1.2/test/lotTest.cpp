// standard headers
#include <vector>
// boost headers
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
// MCMC++ headers
#include <mcmc++/lot.h>

using namespace boost::unit_test_framework;

namespace {

  const double rngPrecision = 1.0e-6;

}

void
ConstructorTest(void) {
  lot* rng = new lot;
  BOOST_CHECK(rng != 0); 
  delete rng;
}

void
KNUIntTest(void) {
  lot rng(lot::RAN_KNU);
  std::vector<long> a(2009);
  rng.ran_start(310952L);
  for (int m = 0; m <= 2009; ++m) {
    rng.ran_array(a, 1009);
  }
  BOOST_CHECK_EQUAL(a[0] , 995235265L);
  rng.ran_start(310952L);
  for (int m = 0; m <= 1009; ++m) {
    rng.ran_array(a,2009);
  }
  BOOST_CHECK_EQUAL(a[0] , 995235265L);
}

void
KNUDoubleTest(void) {
  lot rng(lot::RAN_KNU);
  std::vector<double> a(2009);
  rng.ranf_start(310952L);
  for (int m = 0; m < 2009; ++m) {
    rng.ranf_array(a,1009);
  }
  // N.B.: not full comparison of precision.
  BOOST_CHECK_CLOSE(rng.get_ran_u_0(), 0.36410514377569680455, 1.0e-16);
  rng.ranf_start(310952L);
  for (int m = 0; m < 1009; ++m) {
    rng.ranf_array(a,2009);
  }
  BOOST_CHECK_CLOSE(rng.get_ran_u_0(), 0.36410514377569680455, 1.0e-16);
}

void
MTTest(void) {
  // comparison values come directly from the output of 
  // http://www.math.keio.ac.jp/matumoto/CODES/MT2002/mt19937ar.c
  // (floating point output modified for additional accuracy)
  lot rng(lot::RAN_MT);
  unsigned long init[4]={0x123, 0x234, 0x345, 0x456};
  unsigned long length = 4;
  rng.MT_init_by_array(init, length);
  for (int i = 0; i < 999; ++i) {
    rng.MT_genrand_int();
  }
  BOOST_CHECK_EQUAL(rng.MT_genrand_int(), 3460025646UL);
  for (int i = 0; i < 999; ++i) {
    rng.uniform();
  }
  // accuracy = 1/(2^32) (ca. 2.328306e-10)
  BOOST_CHECK_CLOSE(rng.MT_genrand_with_zero(), 0.7215715157799423, 
                    1.0/4294967296.0);
  // types of MT
  BOOST_CHECK(rng.Set_MT(lot::OPEN));
  BOOST_CHECK(rng.Set_MT(lot::ZERO));
  BOOST_CHECK(rng.Set_MT(lot::ZERO_ONE));
  BOOST_CHECK(!rng.Set_MT(0));
  // make sure copy construction and assignment don't work
#if defined(CHECK_COPY_ASSIGNMENT)
    lot rng2(rng);
    lot rng3;
    rng3 = rng;
#endif
}

void
RTest(void) {
  lot rng;
  rng.MT_R_initialize(1);
  double x;
  x = rng.uniform();
  for (int j = 0; j < 50; ++j) {
    x = rng.uniform();
  }
  BOOST_CHECK_CLOSE(rng.uniform(), 0.47761962213553, rngPrecision);
  BOOST_CHECK_CLOSE(rng.beta(2,3), 0.68772855076513, rngPrecision);
  BOOST_CHECK_CLOSE(rng.expon(), 1.4352853433731, rngPrecision);
  BOOST_CHECK_CLOSE(rng.gamma(0.5, 4.0), 0.6194799961546, rngPrecision);
  BOOST_CHECK_EQUAL(rng.binom(25, 0.2), 4);
  BOOST_CHECK_CLOSE(rng.snorm(), 0.68167111439523, rngPrecision);
  BOOST_CHECK_EQUAL(rng.poisson(4.3), 2);
  BOOST_CHECK_CLOSE(rng.beta(0.5,3), 0.003370016262863173, rngPrecision);
  BOOST_CHECK_EQUAL(rng.binom(100, 0.45), 53);
  BOOST_CHECK_EQUAL(rng.poisson(22.3), 20);
  BOOST_CHECK_CLOSE(rng.cauchy(-1, 4), -2.816071291549914, rngPrecision);
  BOOST_CHECK_CLOSE(rng.cauchy(1.3, 0.1), 1.577732058772176, rngPrecision);
  BOOST_CHECK_CLOSE(rng.chisq(7), 14.17932338443496, rngPrecision);
  BOOST_CHECK_CLOSE(rng.chisq(13), 14.31890931779990, rngPrecision);
  BOOST_CHECK_CLOSE(rng.chisq(37.2), 32.13143324293112, rngPrecision);
  BOOST_CHECK_CLOSE(rng.f(1.3, 7.2), 0.1281989382307877, rngPrecision);
  BOOST_CHECK_CLOSE(rng.f(1, 7), 0.4669727131896132, rngPrecision);
  BOOST_CHECK_CLOSE(rng.f(17, 23), 1.565862756168847, rngPrecision);
  BOOST_CHECK_CLOSE(rng.f(47, 13), 1.270805865921048, rngPrecision);
  BOOST_CHECK_CLOSE(rng.f(47, 13), 3.258986248177208, rngPrecision);
  BOOST_CHECK_CLOSE(rng.f(470, 1), 5.798332616147139, rngPrecision);
  BOOST_CHECK_EQUAL(rng.geom(0.1), 4);
  BOOST_CHECK_EQUAL(rng.geom(0.81), 1);
  BOOST_CHECK_EQUAL(rng.geom(0.177), 17);
  BOOST_CHECK_EQUAL(rng.hypergeom(10, 10, 2), 1);
  BOOST_CHECK_EQUAL(rng.hypergeom(17, 87, 10), 1);
  BOOST_CHECK_EQUAL(rng.hypergeom(17, 87, 27), 4);
  BOOST_CHECK_CLOSE(rng.lnorm(0, 1), 0.57074655359235, rngPrecision);
  BOOST_CHECK_CLOSE(rng.lnorm(-1, 1), 0.1787736013404846, rngPrecision);
  BOOST_CHECK_CLOSE(rng.lnorm(-3.1, 1), 0.01578878771139496, rngPrecision);
  BOOST_CHECK_CLOSE(rng.lnorm(3.1, 0.1), 23.29267338532605, rngPrecision);
  BOOST_CHECK_CLOSE(rng.logis(0, 4), -0.8635735260679713, rngPrecision);
  BOOST_CHECK_CLOSE(rng.logis(4, 4), 9.062313282540877, rngPrecision);
  BOOST_CHECK_CLOSE(rng.logis(-8, 1), -6.001694658347600, rngPrecision);
  std::vector<double> p(5);
  double sum = 0.0;
  for (int i = 0; i < 5; ++i) {
    p[i] = i+1;
    sum += p[i];
  }
  for (int i = 0; i < 5; ++i) {
    p[i] /= sum;
  }
  std::vector<int> n = rng.multinom(25, p);
  std::vector<int> k1(5);
  k1[0] = 1;
  k1[1] = 1;
  k1[2] = 5;
  k1[3] = 9;
  k1[4] = 9;
  for (int i = 0; i < 5; ++i) {
    BOOST_CHECK_EQUAL(n[i], k1[i]);
  }
  n = rng.multinom(77, p);
  std::vector<int> k2(5);
  k2[0] = 4;
  k2[1] = 11;
  k2[2] = 19;
  k2[3] = 23;
  k2[4] = 20;
  for (int i = 0; i < 5; ++i) {
    BOOST_CHECK_EQUAL(n[i], k2[i]);
  }
  BOOST_CHECK_EQUAL(rng.nbinom(7.7, 0.3), 18);
  BOOST_CHECK_EQUAL(rng.nbinom(7, 0.2), 34);
  BOOST_CHECK_EQUAL(rng.nbinom(1, 0.2), 4);
  BOOST_CHECK_CLOSE(rng.t(1), 1.082284217247224, rngPrecision);
  BOOST_CHECK_CLOSE(rng.t(3.1), 0.1233535814493406, rngPrecision);
  BOOST_CHECK_CLOSE(rng.t(17), 0.579037415210614, rngPrecision);
  BOOST_CHECK_CLOSE(rng.t(17), -0.6748728388191424, rngPrecision);
  BOOST_CHECK_CLOSE(rng.weibull(0.3, 1), 0.1579586855386296, rngPrecision);
  BOOST_CHECK_CLOSE(rng.weibull(3.0, 1), 1.269657809132781, rngPrecision);
  BOOST_CHECK_CLOSE(rng.weibull(0.3, 11), 30.74244732942167, rngPrecision);
  BOOST_CHECK_CLOSE(rng.weibull(3.0, 11), 7.612015090903996, rngPrecision);
}

test_suite*
init_unit_test_suite(int ac, char** av) {
  test_suite* test = BOOST_TEST_SUITE("Lot test suite");
  test->add( BOOST_TEST_CASE(&ConstructorTest) );
  test->add( BOOST_TEST_CASE(&KNUIntTest) );
  test->add( BOOST_TEST_CASE(&KNUDoubleTest) );
  test->add( BOOST_TEST_CASE(&MTTest) );
  test->add( BOOST_TEST_CASE(&RTest) );
  return test;
}
