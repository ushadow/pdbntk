SET(BUILD_RELEASE ON)
#SET(BUILD_DEBUG OFF)
#SET(BUILD_PROFILE OFF)
#SET(BUILD_DUMMYRANDOM OFF)
#
#SET(Boost_USE_MULTITHREAD OFF)
       
SET(Boost_USE_STATIC_LIBS ON)
FIND_PACKAGE(Boost COMPONENTS program_options graph REQUIRED)

INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})

FIND_PACKAGE(OpenMP REQUIRED)

FIND_PACKAGE(GMPXX REQUIRED)
INCLUDE_DIRECTORIES(${GMPXX_INCLUDE_DIR})
INCLUDE_DIRECTORIES(${GMP_INCLUDE_DIR})

FILE(GLOB_RECURSE SRC RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp *.h)
SET(PDBNTK_SRC ${SRC})

# Build release library
IF(BUILD_RELEASE)
ADD_LIBRARY(pdbntk STATIC ${PDBNTK_SRC}) 
SET_TARGET_PROPERTIES(pdbntk PROPERTIES COMPILE_FLAGS "-O3 -fPIC -c -fmessage-length=0 -fopenmp")
ENDIF(BUILD_RELEASE)

#IF(BUILD_DEBUG)
#add_library (MocapyD 
#        ${Mocapy_src}
#        ${Integrate_src}
#        ${Functions_src}
#        )
#SET_TARGET_PROPERTIES(MocapyD PROPERTIES COMPILE_FLAGS "-O0 -fPIC -g -Wall -c -fmessage-length=0")
#ENDIF(BUILD_DEBUG)
#
#IF(BUILD_PROFILE)
#add_library (MocapyP 
#        ${Mocapy_src}
#        ${Integrate_src}
#        ${Functions_src}
#        )
#SET_TARGET_PROPERTIES(MocapyP PROPERTIES COMPILE_FLAGS "-pg -fPIC -O3 -Wall -c -fmessage-length=0")
#ENDIF(BUILD_PROFILE)
#
#IF(BUILD_DUMMYRANDOM)
#add_library (MocapyDR 
#        ${Mocapy_src}
#        ${Integrate_src}
#        ${Functions_src}
#        )
#SET_TARGET_PROPERTIES(MocapyDR PROPERTIES COMPILE_FLAGS "-DRANDDUMMY  -O3 -fPIC -Wall -c -fmessage-length=0")
#ENDIF(BUILD_DUMMYRANDOM)
#
#
## install commands
##INSTALL(TARGETS Mocapy 
##        LIBRARY DESTINATION lib
##        ARCHIVE DESTINATION lib
##        )
