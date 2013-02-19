#undef ENABLE_DATA_TABLE_GRAMMAR

// standard headers
#include <vector>
// boost headers
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
// MCMC++ headers
#include <mcmc++/DataTable.h>

using namespace boost::unit_test_framework;

void 
ConstructorTest(void) {
  DataTable<double> dblTable;
  BOOST_CHECK_EQUAL(dblTable.Rows(), 0U);
  DataTable<int> intTable;
  BOOST_CHECK_EQUAL(intTable.Rows(), 0U);
}

void
ReadDataTest(void) {
  DataTable<double> table;
  enum DataTableResult result = table.Read("mosquito-full.txt");
  BOOST_CHECK_EQUAL(result, readSuccess);
  BOOST_CHECK_EQUAL(table.Rows(), 5000U);
  BOOST_CHECK_EQUAL(table.Columns(), 2U);

  result = table.Read("mosquito-fzero.txt");
  BOOST_CHECK_EQUAL(result, notEmptyError);

  table.Flush();
  BOOST_CHECK_EQUAL(table.Rows(), 0U);
  BOOST_CHECK_EQUAL(table.Columns(), 0U);

  result = table.Read("mosquito-fzero.txt");
  BOOST_CHECK_EQUAL(result, readSuccess);
  BOOST_CHECK_EQUAL(table.Rows(), 5000U);
  BOOST_CHECK_EQUAL(table.Columns(), 1U);

  table.Flush();
  result = table.Read("mosquito-thetazero.txt");
  BOOST_CHECK_EQUAL(result, openError);

}

void
AccessDataTest(void) {
  DataTable<double> table;
  enum DataTableResult result = table.Read("mosquito-full.txt");
  BOOST_CHECK_EQUAL(result, readSuccess);
  BOOST_CHECK_EQUAL(table.Columns(), 2);
  BOOST_CHECK(table.ColumnLabel(0) == "f");
  BOOST_CHECK(table.ColumnLabel(1) == "theta");
  BOOST_CHECK_CLOSE(table.Value(4, 0), 0.0267797, 1.0e-6);
}

void
RowLabelTest(void) {
  DataTable<int> table(false, true);
  enum DataTableResult result = table.Read("row-labels.txt");
  BOOST_CHECK_EQUAL(result, readSuccess);
  BOOST_CHECK_EQUAL(table.Rows(), 50U);
  BOOST_CHECK_EQUAL(table.Columns(), 21U);
  BOOST_CHECK_EQUAL(table.Value(7,2), 31);
  BOOST_CHECK_EQUAL(table.Value(10,4), 22);
  BOOST_CHECK_EQUAL(table.Value(2,20), 1);
  BOOST_CHECK_EQUAL(table.Value(3,20), 2);
  BOOST_CHECK(table.RowLabel(8) == "roW9");
  BOOST_CHECK(table.RowLabel(9) == "RoW-10");
  BOOST_CHECK(table.RowLabel(10)== "row-11");
  BOOST_CHECK(table.RowLabel(11) == "ROW_12");
  BOOST_CHECK(table.RowLabel(8) != "row9");
  BOOST_CHECK(table.RowLabel(8) == "roW9");
}

void
ColumnRowLabelTest(void) {
  DataTable<int> table(true, true);
  enum DataTableResult result = table.Read("test-count.txt");
  BOOST_CHECK_EQUAL(result, readSuccess);
  BOOST_CHECK_EQUAL(table.Rows(), 3);
  BOOST_CHECK_EQUAL(table.Columns(), 6);
  BOOST_CHECK_EQUAL(table.ColumnLabels(), 7);
  BOOST_CHECK(table.RowLabel(2) == "4");
  BOOST_CHECK(table.ColumnLabel(2) == "2");
  BOOST_CHECK(table.ColumnLabel(5) == "5");
  BOOST_CHECK(table.ColumnLabel(6) == "6");
  BOOST_CHECK_EQUAL(table.Value(1,3), 2);
}

void
ColumnVectorTest(void) {
  DataTable<double> table;
  enum DataTableResult result = table.Read("mosquito-full.txt");
  BOOST_CHECK_EQUAL(result, readSuccess);
  std::vector<double> x = table.ColumnVector(1);
  BOOST_CHECK_EQUAL(x[0], 0.167647);
  BOOST_CHECK_EQUAL(x[5], 0.0875195);
  BOOST_CHECK_EQUAL(x[6], 0.175942);
}

void
RowVectorTest(void) {
  DataTable<double> table;
  enum DataTableResult result = table.Read("mosquito-full.txt");
  BOOST_CHECK_EQUAL(result, readSuccess);
  std::vector<double> x = table.RowVector(5);
  BOOST_CHECK_EQUAL(x[0], 0.0477034);
  BOOST_CHECK_EQUAL(x[1], 0.0875195);
}

void
ExceptionTest(void) {
  DataTable<double> table;
  enum DataTableResult result = table.Read("mosquito-full.txt");
  BOOST_CHECK_EQUAL(result, readSuccess);
  BOOST_CHECK_THROW(table.Value(-1,0), BadRow);
  BOOST_CHECK_THROW(table.Value(0,-1), BadCol);
  BOOST_CHECK_THROW(table.Value(5000,0), BadRow);
  BOOST_CHECK_THROW(table.Value(0, 2), BadCol);
}

void
SetValueTest(void) {
  DataTable<int> table(true, true);
  enum DataTableResult result = table.Read("test-count.txt");
  BOOST_CHECK_EQUAL(result, readSuccess);
  BOOST_CHECK_EQUAL(table.Value(0, 0), 9);
  BOOST_CHECK_EQUAL(table.Value(0, 2), 5);
  BOOST_CHECK_EQUAL(table.Value(2, 1), 3);
  table.SetValue(0, 0, -1);
  table.SetValue(0, 2, -4);
  table.SetValue(2, 1, -9);
  BOOST_CHECK_EQUAL(table.Value(0, 0), -1);
  BOOST_CHECK_EQUAL(table.Value(0, 2), -4);
  BOOST_CHECK_EQUAL(table.Value(2, 1), -9);
}

void
TestMixedTest(void) {
  DataTable<double> table(true, false);
  enum DataTableResult result = table.Read("sla-data.txt");
  BOOST_CHECK_EQUAL(result, readSuccess);
}

test_suite*
init_unit_test_suite(int ac, char** av) {
  test_suite* test = BOOST_TEST_SUITE("DataTable test suite");
  test->add( BOOST_TEST_CASE(&ConstructorTest) );
  test->add( BOOST_TEST_CASE(&ReadDataTest) );
  test->add( BOOST_TEST_CASE(&AccessDataTest) );
  test->add( BOOST_TEST_CASE(&RowLabelTest) );
  test->add( BOOST_TEST_CASE(&ColumnRowLabelTest) );
  test->add( BOOST_TEST_CASE(&ColumnVectorTest) );
  test->add( BOOST_TEST_CASE(&RowVectorTest) );
  test->add( BOOST_TEST_CASE(&ExceptionTest) );
  test->add( BOOST_TEST_CASE(&SetValueTest) );
  test->add( BOOST_TEST_CASE(&TestMixedTest) );
  return test;
}

