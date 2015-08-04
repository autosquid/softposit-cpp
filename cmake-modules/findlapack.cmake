# - Find LAPACK
# Intel Threading Building Blocks offers a rich and complete approach to expressing parallelism in a C++ program
# www.threadingbuildingblocks.org
#
# The module defines the following variables:
#  LAPACK_FOUND - the system has nlopt
#  LAPACK_INCLUDE_DIR - where to find nlopt.h
#  LAPACK_INCLUDE_DIRS - nlopt includes
#  LAPACK_LIBRARY - where to find the nlopt library
#  LAPACK_LIBRARIES - aditional libraries
#  LAPACK_ROOT_DIR - root dir (ex. /usr/local)

# set LAPACK_INCLUDE_DIR
find_path (LAPACK_INCLUDE_DIR
  NAMES
  lapacke.h
  PATHS
    "/usr/local/include"
    "/usr/local/opt/lapack/include"
    "${LAPACK_ROOT}/include"
  DOC
    "lapack include directory"
)

# set LAPACK_INCLUDE_DIRS
set (LAPACK_INCLUDE_DIRS ${LAPACK_INCLUDE_DIR} )

# set LAPACK_LIBRARY
find_library (LAPACK_LIBRARY
  NAMES blas lapack
  PATHS
    "/usr/local/lib"
    "/usr/local/opt/lapack/lib"
    "${LAPACK_ROOT}/lib"
  DOC "lapack library location"
)

# set LAPACK_LIBRARIES
set ( LAPACK_LIBRARIES ${LAPACK_LIBRARY} )

# root dir
# try to guess root dir from include dir
if ( LAPACK_INCLUDE_DIR )
  string ( REGEX REPLACE "(.*)/include.*" "\\1" LAPACK_ROOT_DIR ${LAPACK_INCLUDE_DIR} )

# try to guess root dir from library dir
elseif ( LAPACK_LIBRARY )
  string ( REGEX REPLACE "(.*)/lib[/|32|64].*" "\\1" LAPACK_ROOT_DIR ${LAPACK_LIBRARY} )
endif ()

# handle REQUIRED and QUIET options
include ( FindPackageHandleStandardArgs )

find_package_handle_standard_args (lapack  DEFAULT_MSG LAPACK_LIBRARY
  LAPACK_INCLUDE_DIR
  LAPACK_INCLUDE_DIRS
  LAPACK_LIBRARIES
  LAPACK_ROOT_DIR
)


mark_as_advanced (
  LAPACK_LIBRARY
  LAPACK_LIBRARIES
  LAPACK_INCLUDE_DIR
  LAPACK_INCLUDE_DIRS
  LAPACK_ROOT_DIR
  LAPACK_INTERFACE_VERSION
)
