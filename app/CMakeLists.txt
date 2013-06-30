SET(BUILD_RELEASE ON)
SET(BUILD_DEBUG OFF)
SET(BUILD_PROFILE OFF)

SET(CMAKE_VERBOSE_MAKEFILE ON) 

SET(Boost_USE_STATIC_LIBS ON)

FIND_PACKAGE(Boost COMPONENTS serialization program_options thread REQUIRED)
FIND_PACKAGE(OpenMP REQUIRED)

#FIND_LIBRARY(PTHREAD pthread)

# Make sure the compiler can find include files from pdbntk.
include_directories (${pdbntk_SOURCE_DIR}/src)

IF (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  SET (PLATFORM "linux")
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
MESSAGE(STATUS "Platform: ${PLATFORM}")
LINK_DIRECTORIES(${pdbntk_SOURCE_DIR}/lib/${PLATFORM})

# find lapack
FIND_LIBRARY(LAPACK_LIBRARY NAMES lapack )
IF(NOT LAPACK_LIBRARY)
message("error: could not find lapack library")
ENDIF(NOT LAPACK_LIBRARY)

# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

if (Fortran_COMPILER_NAME STREQUAL "gfortran")
  # gfortran
  set (CMAKE_FLIB "gfortran")
elseif (Fortran_COMPILER_NAME STREQUAL "g77")
  # g77
  set (CMAKE_FLIB "g2c")
else (Fortran_COMPILER_NAME STREQUAL "gfortran")
#  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
#  message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
#  message ("No optimized Fortran link flags are known, we just try -lgfortran")
  set (CMAKE_FLIB "gfortran")
endif (Fortran_COMPILER_NAME STREQUAL "gfortran")

INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})

IF(BUILD_RELEASE)
SET(CMAKE_CXX_FLAGS "-O3 -fPIC -Wall -fmessage-length=0 -fopenmp")
SET(PDBNTKLIB pdbntk)
ENDIF(BUILD_RELEASE)

IF(BUILD_DEBUG)
SET(CMAKE_CXX_FLAGS "-g -O0 -Wall -fmessage-length=0 -fopenmp")
SET(PDBNTKLIB pdbntk CACHE STRING "Mocapy library used when linking executables")
ENDIF(BUILD_DEBUG)

#IF(BUILD_PROFILE)
#SET(CMAKE_CXX_FLAGS "-pg -O3 -Wall -fmessage-length=0")
#SET(MOCAPYLIB MocapyP CACHE STRING "Mocapy library used when linking executables")
#ENDIF(BUILD_PROFILE)
#
#IF(BUILD_DUMMYRANDOM)
#SET(CMAKE_CXX_FLAGS "-DRANDDUMMY -O3 -Wall -fmessage-length=0")
#SET(MOCAPYLIB MocapyDR CACHE STRING "Mocapy library used when linking executables")
#ENDIF(BUILD_DUMMYRANDOM)
#
#
SET(PROGS parallel_em_hmm)
#SET(PROGS_BOOST_OPTIONS_THREADS mlr-uni)

FOREACH(p ${PROGS})
        add_executable(${p} ${p}.cpp)
        target_link_libraries (${p} ${PDBNTKLIB} ${MOCAPYLIB} ${DAILIB} ${Boost_SERIALIZATION_LIBRARY} ${LAPACK_LIBRARY} ${CMAKE_FLIB})
ENDFOREACH(p)
#
#FOREACH(p ${PROGS_BOOST_OPTIONS_THREADS})
#        add_executable(${p} ${p}.cpp)
#        target_link_libraries (${p} ${MOCAPYLIB} ${Boost_SERIALIZATION_LIBRARY} ${Boost_PROGRAM_OPTIONS_LIBRARY} ${Boost_THREAD_LIBRARY} ${PTHREAD} ${LAPACK_LIBRARY} ${CMAKE_FLIB})
#ENDFOREACH(p)
#
## copy data directory
#exec_program(
#        ${CMAKE_COMMAND}
#        ARGS
#        "-E create_symlink \"${CMAKE_CURRENT_SOURCE_DIR}/data\" \"${CMAKE_CURRENT_BINARY_DIR}/data\""
#        )
##exec_program(
##        ${CMAKE_COMMAND}
##        ARGS
##        "-E copy_directory \"${CMAKE_CURRENT_SOURCE_DIR}/data\" \"${CMAKE_CURRENT_BINARY_DIR}\""
##        )
