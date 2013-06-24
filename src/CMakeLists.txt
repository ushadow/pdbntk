#SET(BUILD_RELEASE ON)
#SET(BUILD_DEBUG OFF)
#SET(BUILD_PROFILE OFF)
#SET(BUILD_DUMMYRANDOM OFF)
#
#SET(Boost_USE_MULTITHREAD OFF)
#       
#FIND_PACKAGE(Boost 1.39)
#
#INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
#LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})

SET(pdbntk_src
 src/parallel_em_engine.cpp 
)
#
#SET(LIBRARY_OUTPUT_PATH ${MOCAPY_SOURCE_DIR}/libs)
#
## build release library
#IF(BUILD_RELEASE)
add_library (pdbntk1 STATIC
  ${pdbntk_src}
        )
      #SET_TARGET_PROPERTIES(pdbntk PROPERTIES COMPILE_FLAGS " -O3 -fPIC -c
      #-fmessage-length=0 -I${pdbntk_SOURCE_DIR}/src")
#ENDIF(BUILD_RELEASE)
#
##IF(BUILD_RELEASE)
##add_library (Mocapy SHARED
##        ${Mocapy_src}
##        ${Integrate_src}
##        ${Functions_src}
##        )
#
##SET_TARGET_PROPERTIES(Mocapy PROPERTIES COMPILE_FLAGS " -O3 -fPIC -c -fmessage-length=0 -I${MOCAPY_SOURCE_DIR}/src")
##ENDIF(BUILD_RELEASE)
#
#
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
