#include <iostream>
// standard headers
#include <vector>
// boost headers
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
// local headers
#include <mcmc++/DataTable.h>
#include <mcmc++/intervals.h>

using namespace boost::unit_test_framework;

namespace {
  const double precision = 1.0e-7;
}

void
quantileTest(void) {
  DataTable<double> table;
  enum DataTableResult result = table.Read("mosquito-full.txt");
  BOOST_CHECK_EQUAL(result, readSuccess);
  BOOST_CHECK_EQUAL(table.Rows(), 5000U);
  BOOST_CHECK_EQUAL(table.Columns(), 2U);
  std::vector<double> x = table.ColumnVector(1);
  BOOST_CHECK_EQUAL(x[0], 0.167647);
  BOOST_CHECK_EQUAL(x[5], 0.0875195);
  BOOST_CHECK_EQUAL(x[6], 0.175942);
  BOOST_CHECK_CLOSE(quantile(x, 0.0),    0.01978140, precision);
  BOOST_CHECK_CLOSE(quantile(x, 0.0001), 0.01978140, precision);
  BOOST_CHECK_CLOSE(quantile(x, 0.0002), 0.0206255688, precision);
  BOOST_CHECK_CLOSE(quantile(x, 0.001),  0.023629172466667, precision);
  BOOST_CHECK_CLOSE(quantile(x, 0.01),   0.033941776333333, precision);
  BOOST_CHECK_CLOSE(quantile(x, 0.025),  0.041603739166667, precision);
  BOOST_CHECK_CLOSE(quantile(x, 0.05),   0.048515735, precision);
  BOOST_CHECK_CLOSE(quantile(x, 0.1),    0.057701953333333, precision);
  BOOST_CHECK_CLOSE(quantile(x, 0.25),   0.078359008333333, precision);
  BOOST_CHECK_CLOSE(quantile(x, 0.5),    0.1094055, precision);
  BOOST_CHECK_CLOSE(quantile(x, 0.75),   0.150324, precision);
  BOOST_CHECK_CLOSE(quantile(x, 0.90),   0.201796066666667, precision);
  BOOST_CHECK_CLOSE(quantile(x, 0.95),   0.23380630, precision);
  BOOST_CHECK_CLOSE(quantile(x, 0.975),  0.2678637, precision);
  BOOST_CHECK_CLOSE(quantile(x, 0.99),   0.30950809, precision);
  BOOST_CHECK_CLOSE(quantile(x, 0.999),  0.388703995, precision);
  BOOST_CHECK_CLOSE(quantile(x, 0.9998), 0.4285371647, precision);
  BOOST_CHECK_CLOSE(quantile(x, 0.9999), 0.431048, precision);
  BOOST_CHECK_CLOSE(quantile(x, 1.0),    0.431048, precision);
}

test_suite*
init_unit_test_suite(int ac, char** av) {
  test_suite* test = BOOST_TEST_SUITE("Intervals test suite");
  test->add( BOOST_TEST_CASE(&quantileTest) );
  return test;
}
