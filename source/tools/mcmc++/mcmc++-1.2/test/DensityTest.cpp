// boost headers
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
// MCMC++ headers
#include <mcmc++/Density.h>
#include <mcmc++/util.h>

using namespace boost::unit_test_framework;

namespace {

  double dbetaPrecision = 1.0e-13;
  double dbinomPrecision = 1.0e-13;
  double dcauchyPrecision = 1.0e-11;
  double dchisqPrecision = 1.0e-11;
  double dexpPrecision = 1.0e-13;
  double dfPrecision = 1.0e-14;
  double dgammaPrecision = 1.0e-13;
  double dgeomPrecision = 1.0e-14;
  double dhyperPrecision = 1.0e-13;
  double dinvgammaPrecision = 1.0e-12;
  double dlnormPrecision = 1.0e-13;
  double dnormPrecision = 1.0e-13;
  double dlogisPrecision = 1.0e-13;
  double dnbinomPrecision = 1.0e-13;
  double dpoisPrecision = 1.0e-13;
  double dtPrecision = 1.0e-12;
  double dweibullPrecision = 1.0e-13;
  double gamlnPrecision = 1.0e-12;

}

void
dbetaTest(void) {
  using Density::dbeta;
  BOOST_CHECK_CLOSE(dbeta(0.5, 1.0, 1.0, false), 1.0,
                    dbetaPrecision);
  BOOST_CHECK_CLOSE(dbeta(0.5, 1.0, 2.0, false), 1.0,
                    dbetaPrecision);
  BOOST_CHECK_CLOSE(dbeta(0.5, 1.0, 4.0, false), 0.5,
                    dbetaPrecision);
  BOOST_CHECK_CLOSE(dbeta(0.5, 1.0, 8.0, false), 0.0625,
                    dbetaPrecision);
  BOOST_CHECK_CLOSE(dbeta(0.2, 1.0, 1.0, false), 1.0,
                    dbetaPrecision);
  BOOST_CHECK_CLOSE(dbeta(0.2, 1.0, 2.0, false), 1.6,
                    dbetaPrecision);
  BOOST_CHECK_CLOSE(dbeta(0.2, 1.0, 4.0, false), 2.048,
                    dbetaPrecision);
  BOOST_CHECK_CLOSE(dbeta(0.2, 1.0, 8.0, false), 1.677721600000001,
                    dbetaPrecision);
  BOOST_CHECK_CLOSE(dbeta(0.2, 1.0, 4.0, true), 0.7168637071772614,
                    dbetaPrecision);
  BOOST_CHECK_CLOSE(dbeta(0.2, 1.0, 8.0, true), 0.5174366824803678,
                    dbetaPrecision);
}

void
dbinomTest(void) {
  using Density::dbinom;
  BOOST_CHECK_CLOSE(dbinom(5, 25, 0.5, false), 0.001583397388458249,
                    dbinomPrecision);
  BOOST_CHECK_CLOSE(dbinom(12, 25, 0.5, false), 0.1549810171127319,
                    dbinomPrecision);
  BOOST_CHECK_CLOSE(dbinom(1, 25, 0.5, false), 7.450580596923841e-07,
                    dbinomPrecision);
  BOOST_CHECK_CLOSE(dbinom(24, 25, 0.5, false), 7.450580596923841e-07,
                    dbinomPrecision);
  BOOST_CHECK_CLOSE(dbinom(12, 25, 0.5, true), -1.864452639803422,
                    dbinomPrecision);
  BOOST_CHECK_CLOSE(dbinom(1, 25, 0.5, true), -14.10980368913043,
                    dbinomPrecision);
}

void
dgammaTest(void) {
  using Density::dgamma;
  BOOST_CHECK_CLOSE(dgamma(1.0, 1.0, 1.0, false), 0.3678794411714423, 
                    dgammaPrecision);
  BOOST_CHECK_CLOSE(dgamma(1.0, 1.0, 0.5, false), 0.2706705664732254, 
                    dgammaPrecision);
  BOOST_CHECK_CLOSE(dgamma(1.0, 2.0, 10.0, false), 0.009048374180359595, 
                    dgammaPrecision);
  BOOST_CHECK_CLOSE(dgamma(1.0, 0.1, 0.5, false), 0.01524661247061709, 
                    dgammaPrecision);
  BOOST_CHECK_CLOSE(dgamma(1.0, 0.1, 0.1, false), 6.0077867261999e-06,
                    dgammaPrecision);
  BOOST_CHECK_CLOSE(dgamma(0.25, 0.1, 0.1, false), 0.03782484120647526,
                    dgammaPrecision);
  BOOST_CHECK_CLOSE(dgamma(1.0, 0.1, 10.0, false),  0.07554920138253074,
                    dgammaPrecision);
  BOOST_CHECK_CLOSE(dgamma(0.25, 0.1, 10.0, false),  0.2835671747172699,
                    dgammaPrecision);
  BOOST_CHECK_CLOSE(dgamma(0.25, 0.1, 0.1, true), -3.2747892174269,
                    dgammaPrecision);
  BOOST_CHECK_CLOSE(dgamma(1.0, 0.1, 10.0, true), -2.58297116103361,
                    dgammaPrecision);
}

void
dhyperTest(void) {
  using Density::dhyper;
  BOOST_CHECK_CLOSE(dhyper(1, 5, 5, 2, false), 0.5555555555555556,
                    dhyperPrecision);
  BOOST_CHECK_CLOSE(dhyper(10, 50, 500, 20, false), 1.358226980349045e-06,
                    dhyperPrecision);
  BOOST_CHECK_CLOSE(dhyper(10, 10, 15, 24, false), 0.6,
                    dhyperPrecision);
  BOOST_CHECK_CLOSE(dhyper(10, 50, 500, 20, true), -13.50933039968045,
                    dhyperPrecision);
  BOOST_CHECK_CLOSE(dhyper(1, 5, 5, 2, true), -0.5877866649021188,
                    dhyperPrecision);
}

void
dnormTest(void) {
  using Density::dnorm;
  BOOST_CHECK_CLOSE(dnorm(0.0, 0.0, 1.0, false), 0.3989422804014327,
                    dnormPrecision);
  BOOST_CHECK_CLOSE(dnorm(10.0, 10.0, 1.0, false), 0.3989422804014327,
                    dnormPrecision);
  BOOST_CHECK_CLOSE(dnorm(0.0, 1.0, 1.0, false), 0.2419707245191434,
                    dnormPrecision);
  BOOST_CHECK_CLOSE(dnorm(10.0, 11.0, 1.0, false), 0.2419707245191434,
                    dnormPrecision);
  BOOST_CHECK_CLOSE(dnorm(0.0, 1.0, 0.25, false), 0.0005353209030595415,
                    dnormPrecision);
  BOOST_CHECK_CLOSE(dnorm(10.0, 11.0, 0.25, false), 0.0005353209030595415,
                    dnormPrecision);
  BOOST_CHECK_CLOSE(dnorm(10.0, 11.0, 1.0, true), -1.418938533204673,
                    dnormPrecision);
  BOOST_CHECK_CLOSE(dnorm(0.0, 1.0, 0.25, true), -7.53264417208478,
                    dnormPrecision);
}

void
dinvgammaTest(void) {
  using Density::dinvgamma;
  BOOST_CHECK_CLOSE(dinvgamma(1.0, 1.0, 1.0, false), 0.3678794411714423, 
                    dinvgammaPrecision);
  BOOST_CHECK_CLOSE(dinvgamma(1.0, 1.0, 0.5, false), 0.3032653298563167,
                    dinvgammaPrecision);
  BOOST_CHECK_CLOSE(dinvgamma(1.0, 2.0, 10.0, false), 0.004539992976248485,
                    dinvgammaPrecision);
  BOOST_CHECK_CLOSE(dinvgamma(1.0, 0.1, 0.5, false), 0.0594852218356498, 
                    dinvgammaPrecision);
  BOOST_CHECK_CLOSE(dinvgamma(1.0, 0.1, 0.1, false), 0.07554920138253074,
                    dinvgammaPrecision);
  BOOST_CHECK_CLOSE(dinvgamma(0.25, 0.1, 0.1, false), 0.2571624316925206,
                    dinvgammaPrecision);
  BOOST_CHECK_CLOSE(dinvgamma(4.0, 0.1, 0.1, false),  0.01772294841982936,
                    dinvgammaPrecision);
  BOOST_CHECK_CLOSE(dinvgamma(0.25, 0.1, 10.0, false),  2.583128674255218e-18,
                    dinvgammaPrecision);
}

void
dcauchyTest(void) {
  using Density::dcauchy;
  BOOST_CHECK_CLOSE(dcauchy(0.0, 0.0, 1.0, false), 0.31830988618379,
                    dcauchyPrecision);
  BOOST_CHECK_CLOSE(dcauchy(2.0, 0.0, 1.0, false), 0.063661977236758,
                    dcauchyPrecision);
  BOOST_CHECK_CLOSE(dcauchy(4.0, 0.0, 1.0, false), 0.018724110951988,
                    dcauchyPrecision);
  BOOST_CHECK_CLOSE(dcauchy(0.0, 2.0, 1.0, false), 0.063661977236758,
                    dcauchyPrecision);
  BOOST_CHECK_CLOSE(dcauchy(0.0, 4.0, 1.0, false), 0.018724110951988,
                    dcauchyPrecision);
  BOOST_CHECK_CLOSE(dcauchy(0.0, 4.0, 3.0, false), 0.038197186342055,
                    dcauchyPrecision);
}

void
dchisqTest(void) {
  using Density::dchisq;
  BOOST_CHECK_CLOSE(dchisq(2.0, 1.0, false), 0.10377687435515,
                    dchisqPrecision);
  BOOST_CHECK_CLOSE(dchisq(4.0, 1.0, false), 0.026995483256594,
                    dchisqPrecision);
  BOOST_CHECK_CLOSE(dchisq(20.0, 10.0, false), 0.0094583187005177,
                    dchisqPrecision);
  BOOST_CHECK_CLOSE(dchisq(20.0, 40.0, false), 0.0018660813139988,
                    dchisqPrecision);
}

void
dexpTest(void) {
  using Density::dexp;
  BOOST_CHECK_CLOSE(dexp(-1.0, 1.0, true), Util::log_dbl_min,
                    dexpPrecision);
  BOOST_CHECK_CLOSE(dexp(0.0, 1.0, false), 1.0,
                    dexpPrecision);
  BOOST_CHECK_CLOSE(dexp(0.0, 1.0/4.0, false), 4.0,
                    dexpPrecision);
  BOOST_CHECK_CLOSE(dexp(2.0, 1.0, false), 0.1353352832366127,
                    dexpPrecision);
  BOOST_CHECK_CLOSE(dexp(4.0, 1.0, false), 0.01831563888873418,
                    dexpPrecision);
  BOOST_CHECK_CLOSE(dexp(4.0, 1.0/4.0, false), 4.501406988770365e-07,
                    dexpPrecision);
}

void
dfTest(void) {
  using Density::df;
  BOOST_CHECK_CLOSE(df(0.01, 2.0, 7.0, false), 0.9872432555205566,
                    dfPrecision);
  BOOST_CHECK_CLOSE(df(2.0, 2.0, 7.0, false), 0.1308199855581678,
                    dfPrecision);
  BOOST_CHECK_CLOSE(df(10.0, 2.0, 7.0, false), 0.002300404674081928,
                    dfPrecision);
  BOOST_CHECK_CLOSE(df(0.01, 20.0, 7.0, false), 3.515190882096264e-11,
                    dfPrecision);
  BOOST_CHECK_CLOSE(df(2.0, 20.0, 7.0, false), 0.1802546569667958,
                    dfPrecision);
  BOOST_CHECK_CLOSE(df(10.0, 20.0, 7.0, false), 0.0007150525359714545,
                    dfPrecision);
}
void
dgeomTest(void) {
  using Density::dgeom;
  BOOST_CHECK_CLOSE(dgeom(2, 0.5, false), 0.125,
                    dgeomPrecision);
  BOOST_CHECK_CLOSE(dgeom(2, 0.75, false), 0.046875,
                    dgeomPrecision);
  BOOST_CHECK_CLOSE(dgeom(2, 0.9, false), 0.008999999999999996,
                    dgeomPrecision);
  BOOST_CHECK_CLOSE(dgeom(4, 0.9, false), 8.999999999999992e-05,
                    dgeomPrecision);
  BOOST_CHECK_CLOSE(dgeom(2, 0.1, false), 0.08100000000000002,
                    dgeomPrecision);
  BOOST_CHECK_CLOSE(dgeom(4, 0.1, false), 0.06561,
                    dgeomPrecision);
}

void
dlnormTest(void) {
  using Density::dlnorm;
  BOOST_CHECK_CLOSE(dlnorm(2, 1, 0.2, false), 0.3073921572711245,
                    dlnormPrecision);
  BOOST_CHECK_CLOSE(dlnorm(4, 3.2, 2.1, false), 0.03270809710914876,
                    dlnormPrecision);
  BOOST_CHECK_CLOSE(dlnorm(8, 3.2, 2.1, false), 0.02059552400964706,
                    dlnormPrecision);
}

void
dlogisTest(void) {
  using Density::dlogis;
  BOOST_CHECK_CLOSE(dlogis(-4, 0, 1, false), 0.01766270621329112,
                    dlogisPrecision);
  BOOST_CHECK_CLOSE(dlogis(0, 0, 1, false), 0.25,
                    dlogisPrecision);
  BOOST_CHECK_CLOSE(dlogis(1, 0, 1, false), 0.1966119332414819,
                    dlogisPrecision);
  BOOST_CHECK_CLOSE(dlogis(4.3, 1.2, 3.4, false), 0.06014030111045845,
                    dlogisPrecision);
  BOOST_CHECK_CLOSE(dlogis(20.4, 1.2, 3.4, false), 0.001030328501897008,
                    dlogisPrecision);
}

void
dnbinomTest(void) {
  using Density::dnbinom;
  BOOST_CHECK_CLOSE(dnbinom(0, 1, 0.2, false), 0.2,
                    dnbinomPrecision);
  BOOST_CHECK_CLOSE(dnbinom(0, 2, 0.2, false), 0.04000000000000001,
                    dnbinomPrecision);
  BOOST_CHECK_CLOSE(dnbinom(3, 2.3, 0.2, false), 0.06873816465951731,
                    dnbinomPrecision);
  BOOST_CHECK_CLOSE(dnbinom(13, 2.3, 0.2, false), 0.03643153869578965,
                    dnbinomPrecision);
  BOOST_CHECK_CLOSE(dnbinom(27, 2.3, 0.2, false), 0.003918649731744054,
                    dnbinomPrecision);
}

void
dpoisTest(void) {
  using Density::dpois;
  BOOST_CHECK_CLOSE(dpois(0, 0.1, false), 0.9048374180359595,
                    dpoisPrecision);
  BOOST_CHECK_CLOSE(dpois(0, 7.1, false), 0.0008251049232659046,
                    dpoisPrecision);
  BOOST_CHECK_CLOSE(dpois(0, 17.1, false), 3.745970556295245e-08,
                    dpoisPrecision);
  BOOST_CHECK_CLOSE(dpois(7, 17.1, false), 0.003177653916729565,
                    dpoisPrecision);
  BOOST_CHECK_CLOSE(dpois(17, 17.1, false), 0.09625642348924972,
                    dpoisPrecision);
}

void
dweibullTest(void) {
  using Density::dweibull;
  BOOST_CHECK_CLOSE(dweibull(0.1, 2.1, 1, false), 0.1654891674285188,
                    dweibullPrecision);
  BOOST_CHECK_CLOSE(dweibull(11.1, 2.1, 1, false), 2.517712990224233e-67,
                    dweibullPrecision);
  BOOST_CHECK_CLOSE(dweibull(11.1, 2.1, 4.3, false), 0.0009119002829841303,
                    dweibullPrecision);
}

void dtTest(void) {
  using Density::dt;
  BOOST_CHECK_CLOSE(dt(1.2, 3.2, false), 0.1692811541918251,
                    dtPrecision);
  BOOST_CHECK_CLOSE(dt(3.56, 7, false), 0.006170322975786277,
                    dtPrecision);
  BOOST_CHECK_CLOSE(dt(7, 7, false), 9.399205342584660e-05,
                    dtPrecision);
  BOOST_CHECK_CLOSE(dt(1.2, 3.2, true), -1.776194312094514,
                    dtPrecision);
  BOOST_CHECK_CLOSE(dt(3.56, 7, true), -5.088004096278244,
                    dtPrecision);
  BOOST_CHECK_CLOSE(dt(7, 7, true), -9.272300317290719,
                    dtPrecision);
}

void gamlnTest(void) {
  using Density::gamln;
  using Density::MinGammaPar;
  using Density::MaxGammaPar;
  BOOST_CHECK_CLOSE(gamln(sqrt(MinGammaPar)), 44.3614195558365, gamlnPrecision);
  BOOST_CHECK_CLOSE(gamln(MinGammaPar), 88.722839111673, gamlnPrecision);
  BOOST_CHECK_CLOSE(gamln(1.0), 0.0, gamlnPrecision);
  BOOST_CHECK_CLOSE(gamln(2.0), 0.0, gamlnPrecision);
  BOOST_CHECK_CLOSE(gamln(4.0), 1.791759469228055, gamlnPrecision);
  BOOST_CHECK_CLOSE(gamln(16.0), 27.89927138384089, gamlnPrecision);
  BOOST_CHECK_CLOSE(gamln(256.0), 1161.712101118401, gamlnPrecision);
  BOOST_CHECK_CLOSE(gamln(sqrt(MaxGammaPar)), 799877009219260383232.0, gamlnPrecision);
  BOOST_CHECK_CLOSE(gamln(MaxGammaPar), 2.985053532594476e+40, gamlnPrecision);
}

test_suite*
init_unit_test_suite(int ac, char** av) {
  test_suite* test = BOOST_TEST_SUITE("Util test suite");
  test->add( BOOST_TEST_CASE(&dbetaTest) );
  test->add( BOOST_TEST_CASE(&dbinomTest) );
  test->add( BOOST_TEST_CASE(&dcauchyTest) );
  test->add( BOOST_TEST_CASE(&dchisqTest) );
  test->add( BOOST_TEST_CASE(&dexpTest) );
  test->add( BOOST_TEST_CASE(&dfTest) );
  test->add( BOOST_TEST_CASE(&dgammaTest) );
  test->add( BOOST_TEST_CASE(&dgeomTest) );
  test->add( BOOST_TEST_CASE(&dinvgammaTest) );
  test->add( BOOST_TEST_CASE(&dhyperTest) );
  test->add( BOOST_TEST_CASE(&dlnormTest) );
  test->add( BOOST_TEST_CASE(&dlogisTest) );
  test->add( BOOST_TEST_CASE(&dnbinomTest) );
  test->add( BOOST_TEST_CASE(&dnormTest) );
  test->add( BOOST_TEST_CASE(&dpoisTest) );
  test->add( BOOST_TEST_CASE(&dweibullTest) );
  test->add( BOOST_TEST_CASE(&dtTest) );
  test->add( BOOST_TEST_CASE(&gamlnTest) );
  return test;
}
