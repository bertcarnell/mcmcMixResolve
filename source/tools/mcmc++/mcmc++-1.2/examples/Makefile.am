SOURCE_DIR = ../src

noinst_PROGRAMS = binomial normal multinomial normal-mixture

binomial_SOURCES       = binomial.cpp
normal_SOURCES         = normal.cpp
multinomial_SOURCES    = multinomial.cpp
normal_mixture_SOURCES = normal-mixture.cpp

AM_LDFLAGS = -L$(SOURCE_DIR)

INCLUDES = -I$(SOURCE_DIR)
LDADD = -lmcmc++ -lstdc++
AM_CPPFLAGS = ${RELEASE}

EXTRA_DIST = norm.txt multinomial.txt

