CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
SET(CMAKE_LEGACY_CYGWIN_WIN32 0)
PROJECT(pdbntk C CXX Fortran)

SET(CMAKE_VERBOSE_MAKEFILE ON) 

SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/lib/include
    ${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/gtest-1.6.0/include)

SET(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)
SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

SET(MOCAPYLIB Mocapy)
SET(DAILIB dai)
SET(GTESLIB gtest)
SET(GTESTMAINLIB gtest_main)
SET(GLOG glog)

# Recurse into the "src" and "examples" subdirectories. This does not actually
# cause another cmake executable to run. The same process will walk through
# the project's entire directory structure.
ADD_SUBDIRECTORY(${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/gtest-1.6.0)
ADD_SUBDIRECTORY(${CMAKE_CURRENT_SOURCE_DIR}/src)
ADD_SUBDIRECTORY(${CMAKE_CURRENT_SOURCE_DIR}/app)
ADD_SUBDIRECTORY(${CMAKE_CURRENT_SOURCE_DIR}/test)

# settings for making distribution package
#set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Mocapy is a Dynamic Bayesian Network toolkit")
#set(CPACK_PACKAGE_VENDOR "Thomas Hamelryck and Martin Paluszewski")
##set(CPACK_PACKAGE_DESCRIPTION_FILE ${CMAKE_CURRENT_SOURCE_DIR}/README)
#set(CPACK_GENERATOR TGZ)
#set(
#        CPACK_SOURCE_PACKAGE_FILE_NAME
#        "Mocapy-${MOCAPY_VERSION}"
#        CACHE INTERNAL "tarball basename"
#        )
#set(CPACK_SOURCE_GENERATOR TGZ)
#include(CPack)

