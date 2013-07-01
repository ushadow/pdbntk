SET(BUILD_RELEASE ON)
SET(BUILD_DEBUG OFF)
SET(BUILD_PROFILE OFF)
SET(BUILD_DUMMYRANDOM OFF)

SET(Boost_ADDITIONAL_VERSIONS "1.39" "1.39.0" "1.38" "1.38.0" "1.37" "1.37.0" "1.36" "1.36.0" )
SET(Boost_USE_STATIC_LIBS ON)

FIND_PACKAGE(Boost COMPONENTS serialization 1.39)
IF(NOT BOOST_serialization_FOUND )
       FIND_PACKAGE(Boost COMPONENTS serialization 1.39.0)
ENDIF(NOT BOOST_serialization_FOUND)
IF(NOT BOOST_serialization_FOUND )
       FIND_PACKAGE(Boost COMPONENTS serialization 1.38)
ENDIF(NOT BOOST_serialization_FOUND)
IF(NOT BOOST_serialization_FOUND )
                FIND_PACKAGE(Boost COMPONENTS serialization 1.38.0)
ENDIF(NOT BOOST_serialization_FOUND)
IF(NOT BOOST_serialization_FOUND )
		FIND_PACKAGE(Boost 1.37 COMPONENTS serialization )
ENDIF(NOT BOOST_serialization_FOUND)
IF(NOT BOOST_serialization_FOUND )
                FIND_PACKAGE(Boost COMPONENTS serialization 1.37.0)
ENDIF(NOT BOOST_serialization_FOUND)
IF(NOT BOOST_serialization_FOUND )
		FIND_PACKAGE(Boost 1.36 COMPONENTS serialization )
ENDIF(NOT BOOST_serialization_FOUND)
IF(NOT BOOST_serialization_FOUND )
                FIND_PACKAGE(Boost COMPONENTS serialization 1.36.0)
ENDIF(NOT BOOST_serialization_FOUND)

# Make sure the compiler can find include files from mocapy.
include_directories (${pdbntk_SOURCE_DIR}/src)

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
LINK_DIRECTORIES(${Boost_LIBRARY_DIRS} ${LIBRARY_OUTPUT_PATH})

IF (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  SET (PLATFORM "linux")
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
MESSAGE(STATUS "Platform: ${PLATFORM}")
LINK_DIRECTORIES(${pdbntk_SOURCE_DIR}/lib/${PLATFORM})
LINK_DIRECTORIES(/usr/local/lib)

IF(BUILD_RELEASE)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -fmessage-length=0")
SET(PDBNTKLIB pdbntk)
ENDIF(BUILD_RELEASE)

IF(BUILD_DEBUG)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 -Wall -fmessage-length=0")
SET(PDBNTKLIB pdbntkD)
ENDIF(BUILD_DEBUG)

IF(BUILD_PROFILE)
SET(CMAKE_CXX_FLAGS "-pg -O3 -Wall -fmessage-length=0")
SET(MOCAPYLIB MocapyP CACHE STRING "Mocapy library used when linking executables")
ENDIF(BUILD_PROFILE)

IF(BUILD_DUMMYRANDOM)
SET(CMAKE_CXX_FLAGS "-DRANDDUMMY -O3 -Wall -fmessage-length=0")
SET(MOCAPYLIB MocapyDR CACHE STRING "Mocapy library used when linking executables")
ENDIF(BUILD_DUMMYRANDOM)

SET(PROGS jtree_test)

FOREACH(p ${PROGS})
        add_executable(${p} ${p}.cpp)
        target_link_libraries (${p} ${PDBNTKLIB} ${MOCAPYLIB} ${DAILIB} ${GLOG}
            ${GTESLIB} ${GTESTMAINLIB} ${Boost_SERIALIZATION_LIBRARY}
            ${LAPACK_LIBRARY} ${CMAKE_FLIB})
ENDFOREACH(p)
