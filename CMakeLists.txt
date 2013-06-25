cmake_minimum_required(VERSION 2.8)
set(CMAKE_LEGACY_CYGWIN_WIN32 0)
project(pdbntk)

include_directories(lib/include)
# Recurse into the "src" and "examples" subdirectories. This does not actually
# cause another cmake executable to run. The same process will walk through
# the project's entire directory structure.
add_subdirectory(src)
#add_subdirectory (examples)
#add_subdirectory (tests)

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

