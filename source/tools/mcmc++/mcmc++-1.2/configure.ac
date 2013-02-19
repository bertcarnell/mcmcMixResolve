#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_INIT(mcmc++, 1.2, kent@darwin.eeb.uconn.edu, mcmc++)
AC_CONFIG_SRCDIR([src/DataTable.cpp])
AC_CONFIG_SRCDIR([test/DataTableTest.cpp])
AC_CONFIG_SRCDIR([examples/binomial.cpp])

# Handle arguments to configure.
config_opts="${ac_configure_args}"
AC_SUBST(config_opts)

# default is to define NDEBUG during compile
# test the argument
AC_ARG_ENABLE([debug],
[AC_HELP_STRING([--enable-debug],
                [include assertions for debugging @<:@no@:>@])],
[if test "${enableval}" = "yes"
  then
    RELEASE=""
  else
    RELEASE="-DNDEBUG"
  fi],
[RELEASE="-DNDEBUG"])
AC_SUBST([RELEASE])

# default SUFFIX is blank. Depending on how Boost is installed, a suffix
# may be necessary for the linker to find the appropriate libs for
# testing.
AC_ARG_WITH([boost-test],
[AC_HELP_STRING([--with-boost-test@<:@=suffix@:>@],
                [suffix of boost_unit_test_framework @<:@(none)@:>@])],
[if test "${withval}" = "no"
  then
    SUFFIX=""
  else
    SUFFIX="${withval}"
  fi
])
BOOST_UNIT_TEST=boost_unit_test_framework${SUFFIX}
AC_SUBST([BOOST_UNIT_TEST]) 

# sets build system type
# N.B.: This set of scripts assumes that the extension for library files
# is .a even under Windows. bjam gives windows libraries a .lib extension
# under MinGW, so it will be necessary to rename them for testing to work.
build_system=`uname`
case $build_system in
  MINGW32*)
    BUILD_SYSTEM="-DMINGW32=1"
    ;;
  *)
    BUILD_SYSTEM=
    ;;
esac
AC_SUBST([BUILD_SYSTEM])

# Automake
AM_INIT_AUTOMAKE

# Checks for programs.
AC_PROG_CXX
AC_PROG_MAKE_SET

# Makefiles
AC_CONFIG_FILES([Makefile src/Makefile test/Makefile examples/Makefile])

# Start C++ specific checks
AC_LANG_PUSH(C++)

# Checks for header files
AC_CHECK_HEADER([boost/format.hpp],[],[echo "Use CPPFLAGS to set include path for Boost"; exit])
AC_CHECK_HEADER([boost/progress.hpp],[],[echo "Use CPPFLAGS to set include path for Boost"; exit])
AC_CHECK_HEADER([boost/spirit/core.hpp],[],[echo "Use CPPFLAGS to set include path for Boost"; exit])
AC_CHECK_HEADER([boost/spirit/utility/lists.hpp],[],[echo "Use CPPFLAGS to set include path for Boost"; exit], [#include <boost/spirit/core.hpp>])
AC_CHECK_HEADER([boost/static_assert.hpp],[],[echo "Use CPPFLAGS to set include path for Boost"; exit])
AC_CHECK_HEADER([boost/tokenizer.hpp],[],[echo "Use CPPFLAGS to set include path for Boost"; exit])
AC_CHECK_HEADER([boost/type_traits.hpp],[],[echo "Use CPPFLAGS to set include path for Boost"; exit])

# Checks for library functions.
AC_CHECK_FUNCS([floor pow sqrt expm1 isnan])

# End C++ specific checks
AC_LANG_POP(C++)


# library
AC_PROG_RANLIB

AC_OUTPUT
