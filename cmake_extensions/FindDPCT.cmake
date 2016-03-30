# This module finds an installed DPCT package.
#
# It sets the following variables:
#  DPCT_FOUND              - Set to false, or undefined, if lemon isn't found.
#  DPCT_INCLUDE_DIR        - Lemon include directory.
#  DPCT_LIBRARIES          - Lemon library files
FIND_PATH(DPCT_INCLUDE_DIR dpct/graph.h PATHS /usr/include /usr/local/include ${CMAKE_INCLUDE_PATH} ${CMAKE_PREFIX_PATH}/include $ENV{DPCT_ROOT}/include ENV CPLUS_INCLUDE_PATH)
FIND_LIBRARY(DPCT_LIBRARIES
  NAMES dpct
  PATHS $ENV{DPCT_ROOT} ENV LD_LIBRARY_PATH ENV LIBRARY_PATH
)

GET_FILENAME_COMPONENT(DPCT_LIBRARY_PATH ${DPCT_LIBRARIES} PATH)
SET( DPCT_LIBRARY_DIR ${DPCT_LIBRARY_PATH} CACHE PATH "Path to dpct library.")

# handle the QUIETLY and REQUIRED arguments and set DPCT_FOUND to TRUE if
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(DPCT DEFAULT_MSG DPCT_LIBRARIES DPCT_INCLUDE_DIR)

MARK_AS_ADVANCED( DPCT_INCLUDE_DIR DPCT_LIBRARIES DPCT_LIBRARY_DIR )
