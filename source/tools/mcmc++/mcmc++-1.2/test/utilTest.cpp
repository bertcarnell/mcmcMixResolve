// standard headers
#include <iostream>
#include <vector>
// boost headers
#include <boost/test/unit_test.hpp>
// local headers
#include "mcmc++/util.h"

using std::vector;
using namespace boost::unit_test_framework;

namespace {

  class BadAssert {};

}

void
AssertTest(void) {
  using Util::Assert;
  BOOST_CHECK_THROW(Assert<BadAssert>(false), BadAssert);
  BOOST_CHECK_NO_THROW(Assert<BadAssert>(true));
}

void
vector_castTest(void) {
  using Util::vector_cast;
  vector<int> intVector(5);
  for (int i = 0; i < 5; ++i) {
    intVector[i] = i;
  }
  vector<double> toDoubleVector = vector_cast<double>(intVector);
  BOOST_CHECK_EQUAL(intVector[0], toDoubleVector[0]);
  BOOST_CHECK_EQUAL(intVector[1], toDoubleVector[1]);
  BOOST_CHECK_EQUAL(intVector[2], toDoubleVector[2]);
  BOOST_CHECK_EQUAL(intVector[3], toDoubleVector[3]);
  BOOST_CHECK_EQUAL(intVector[4], toDoubleVector[4]);
  vector<double> doubleVector(5);
  for (int i = 0; i < 5; ++i) {
    doubleVector[i] = intVector[i] + static_cast<double>(i)/5;
  }
  vector<int> toIntVector = vector_cast<int>(doubleVector);
  BOOST_CHECK_EQUAL(intVector[0], toIntVector[0]);
  BOOST_CHECK_EQUAL(intVector[1], toIntVector[1]);
  BOOST_CHECK_EQUAL(intVector[2], toIntVector[2]);
  BOOST_CHECK_EQUAL(intVector[3], toIntVector[3]);
  BOOST_CHECK_EQUAL(intVector[4], toIntVector[4]);
}

void FlushVectorTest(void) {
  using Util::FlushVector;
  using std::cout;
  using std::endl;
  const int vSize = 484;
  vector<double> x(vSize, 1.0);
  BOOST_CHECK_EQUAL(x.size(), vSize);
  FlushVector(x);
  BOOST_CHECK_EQUAL(x.size(), 0);
  BOOST_CHECK_EQUAL(x.capacity(), 0);
  vector<int> y(vSize, 1);
  BOOST_CHECK_EQUAL(y.size(), vSize);
  FlushVector(y);
  BOOST_CHECK_EQUAL(y.size(), 0);
  BOOST_CHECK_EQUAL(y.capacity(), 0);
  vector<unsigned> z(vSize, 1);
  BOOST_CHECK_EQUAL(z.size(), vSize);
  FlushVector(z);
  BOOST_CHECK_EQUAL(z.size(), 0);
  BOOST_CHECK_EQUAL(z.capacity(), 0);
}

test_suite*
init_unit_test_suite(int ac, char** av) {
  test_suite* test = BOOST_TEST_SUITE("Util test suite");
  test->add( BOOST_TEST_CASE(&vector_castTest) );
  test->add( BOOST_TEST_CASE(&AssertTest) );
  test->add( BOOST_TEST_CASE(&FlushVectorTest) );
  return test;
}
