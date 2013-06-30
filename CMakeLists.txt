cmake_minimum_required(VERSION 2.8)
set(CMAKE_LEGACY_CYGWIN_WIN32 0)
project(pdbntk C CXX Fortran)

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/include
    /usr/include
    ${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/gteset-1.6.0)

SET(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)
SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

SET(MOCAPYLIB Mocapy)
SET(DAILIB dai)
SET(GTESLIB gteset)

# Recurse into the "src" and "examples" subdirectories. This does not actually
# cause another cmake executable to run. The same process will walk through
# the project's entire directory structure.
add_subdirectory (${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/gtest-1.6.0)
add_subdirectory (${CMAKE_CURRENT_SOURCE_DIR}/src)
add_subdirectory (${CMAKE_CURRENT_SOURCE_DIR}/app)
add_subdirectory (${CMAKE_CURRENT_SOURCE_DIR}/test)

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

