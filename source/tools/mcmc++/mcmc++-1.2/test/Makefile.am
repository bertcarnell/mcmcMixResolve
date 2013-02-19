SOURCE_DIR = ../src

noinst_PROGRAMS = DataTableTest DensityTest lotTest utilTest MCMCTest \
intervalsTest

DataTableTest_SOURCES = DataTableTest.cpp
DensityTest_SOURCES   = DensityTest.cpp
lotTest_SOURCES       = lotTest.cpp
utilTest_SOURCES      = utilTest.cpp
MCMCTest_SOURCES      = MCMCTest.cpp
intervalsTest_SOURCES = intervalsTest.cpp

AM_LDFLAGS = -L$(SOURCE_DIR)

INCLUDES = -I$(SOURCE_DIR)
LDADD = -lmcmc++ -lstdc++ -l$(BOOST_UNIT_TEST)

TESTS = $(noinst_PROGRAMS)

EXTRA_DIST = mosquito-full.txt mosquito-fzero.txt row-labels.txt \
    test-count.txt
