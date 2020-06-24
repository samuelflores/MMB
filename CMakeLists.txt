# Create MMB build system using CMake.
#
# We check the BUILD_MMB_{SHARED,STATIC} variables to determine
# whether to build dynamically linked, statically linked, or both
# versions of the executable.

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(MMB)

SET( MMB_MAJOR_VERSION 3 )
SET( MMB_MINOR_VERSION 01)
SET( MMB_PATCH_VERSION 0 )

SET(CMAKE_BUILD_TYPE "Release" CACHE STRING "This is a Release build by default" )
#SET(CMAKE_GENERATOR "Xcode" CACHE STRING "Since you on a Mac, we assume you want to use Xcode" FORCE) 
#SET(CMAKE_GENERATOR "Xcode" CACHE STRING "Since you on a Mac, we assume you want to use Xcode" FORCE) 

SET(Build_NtC            FALSE CACHE BOOL 
        "Build NtC features of MMB. This greatly increases the size of the BiopolymerClassContainer object in memory")
SET(Build_CCP4           FALSE CACHE BOOL 
        "Build CCP4 maps support features for MMB?")
SET(BUILD_MMB_SHARED_LIB TRUE CACHE BOOL 
        "Build MMB dynamic library")
SET(BUILD_MMB_SHARED_EXEC TRUE CACHE BOOL 
        "Build MMB executable using shared libraries")
SET(BUILD_MMB_STATIC FALSE CACHE BOOL
        "Build MMB using static SimTK libraries.")
SET(BUILD_MMB_GUI_TEST FALSE CACHE BOOL
        "Build MMB_gui test program")
SET(SeqAn_INCLUDE_DIR "" CACHE PATH
        "Where to find Seqan headers")
SET(OpenMM_INCLUDE_DIR "" CACHE PATH
        "Where to find OpenMM main header")
#SET(OpenMM_INCLUDE_REFERENCE_DIR "" CACHE PATH
#        "Where to find OpenMM reference files")
#SET(OpenMM_INCLUDE_SUB_DIR "" CACHE PATH
#        "I don't know how to explain what this is except on my machine it's /usr/local/openmm/include/openmm")
SET(OpenMM_INSTALL_DIR "" CACHE PATH
        "OpenMM install directory. Used to find libs.")
#Add the OpenMM library
SET(LIBCCP4_INCLUDE_DIR "" CACHE PATH
        "Where to find libccp4 includes")
SET(LIBCCP4_LIB_DIR "" CACHE PATH
        "Where to find libccp4 library")

SET(OpenMM_SHARED_LIBS
        optimized OpenMM)

SET(OpenMM_SHARED_LIBS_D
        debug OpenMM_d)

SET(SimTK_SHARED_LIBS 
        optimized SimTKmolmodel 
        optimized SimTKsimbody 
        optimized SimTKmath
        optimized SimTKcommon
    )
SET(SimTK_SHARED_LIBS_D
        debug SimTKmolmodel_d 
        debug SimTKsimbody_d 
        debug SimTKmath_d
        debug SimTKcommon_d)

SET(SimTK_STATIC_LIBS 
        optimized SimTKmolmodel_static 
        optimized SimTKsimbody_static 
        optimized SimTKmath_static
        optimized SimTKcommon_static)
SET(SimTK_STATIC_LIBS_D
        debug SimTKmolmodel_static_d 
        debug SimTKsimbody_static_d 
        debug SimTKmath_static_d
        debug SimTKcommon_static_d)


# These extra libraries are only available as a shared, optimized.
IF(WIN32)
    SET(SimTK_GENERAL_LIBS SimTKlapack pthreadVC2 ws2_32)
ELSEIF(APPLE)
    SET(SimTK_GENERAL_LIBS lapack blas)
ELSE() #Linux
    SET(SimTK_GENERAL_LIBS lapack blas dl rt)
ENDIF()

#simbody
SET(SimTK_INSTALL_DIR $ENV{SimTK_INSTALL_DIR} CACHE PATH
    "Directory where SimTK simbody is installed, e.g. /usr/local/SimTK.")

IF(NOT SimTK_INSTALL_DIR)
    MESSAGE(FATAL_ERROR 
"Expected SimTK_INSTALL_DIR environment var to be set, or set SimTK in CMake")
ENDIF(NOT SimTK_INSTALL_DIR)

find_package(Simbody CONFIG REQUIRED HINTS "${SimTK_INSTALL_DIR}/lib/cmake")

#molmodel
SET(SimTKmolmodel ${SimTK_INSTALL_DIR} CACHE PATH
    "Directory where SimTK molmodel is installed, e.g. /usr/local/SimTK.")

IF(NOT SimTKmolmodel)
    MESSAGE(FATAL_ERROR 
"Expected SimTK_INSTALL_DIR environment var to be set, or set SimTKmolmodel in CMake")
ENDIF(NOT SimTKmolmodel)

# If we're building on a 64 bit Linux machine, look in the lib64 directory
# rather than just lib. (This might not be right on all Linux versions.)

# As of Ubuntu 12, 64 bit libraries are in lib, not lib64, so we are 
# not setting this on any platform now. If you are building on an
# older system you might have to set it manually.
SET(LIB64)
IF( UNIX AND NOT APPLE AND ${CMAKE_SIZEOF_VOID_P} EQUAL 8 )
    SET(LIB64 64) # LIB64 is just blank, in keeping with new molmodel practice
ENDIF()

# Make everything go in the same binary directory. (These are CMake-defined
# variables.)
SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
SET(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})

# for some reason it's not reading the ${OpenMM_INCLUDE_DIR}, so i hard coded the OpenMM directory: 
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/include ${Simbody_INCLUDE_DIR} ${SimTKmolmodel}/include ${OpenMM_INCLUDE_DIR} ${OpenMM_INCLUDE_DIR}/openmm ${OpenMM_INCLUDE_DIR}/openmm/reference ${SeqAn_INCLUDE_DIR}  )
LINK_DIRECTORIES(${Simbody_LIB_DIR}${LIB64} ${SimTKmolmodel}/lib ${CMAKE_BINARY_DIR} ${OpenMM_INSTALL_DIR}/lib   )

##################
### Add CCP4 Paths
##################
IF    ( Build_CCP4 )
	INCLUDE_DIRECTORIES ( ${LIBCCP4_INCLUDE_DIR} )
	LINK_DIRECTORIES    ( ${LIBCCP4_LIB_DIR}     )
	ADD_DEFINITIONS     ( -DCPP4_MAPS_USAGE )
ENDIF ( Build_CCP4 )

# On Mac, build 32 bit binaries.
IF(APPLE)
    SET(CMAKE_OSX_ARCHITECTURES "i386" CACHE STRING 
        "The processor architectures to build for")
ENDIF(APPLE)

# On Mac, use Xcode by default
IF(APPLE)
    SET(CMAKE_GENERATOR "Xcode" CACHE STRING "Since you on a Mac, we assume you want to use Xcode"  ) 
ELSEIF(UNIX)
    SET(CMAKE_GENERATOR "Unix Makefiles" CACHE STRING "Since you on Unix, we assume you want to use a Makefile" ) 
ENDIF()

# List source files explicitly here since there is junk in the src
# directory that we don't want to include here.
SET(MMB_LIBRARY_SOURCE_FILES
    ${PROJECT_SOURCE_DIR}/src/Utils.cpp  
    ${PROJECT_SOURCE_DIR}/src/BiopolymerClass.cpp  
    #//${PROJECT_SOURCE_DIR}/src/elliptic_integral.cpp
    ${PROJECT_SOURCE_DIR}/src/UnitCellParameters.cpp  
    # This should remain commented out:
    #${PROJECT_SOURCE_DIR}/src/NtCInteractionParameterReader.cpp 
    ${PROJECT_SOURCE_DIR}/src/AtomSpringContainer.cpp  
    ${PROJECT_SOURCE_DIR}/src/DisplacementContainer.cpp  
    ${PROJECT_SOURCE_DIR}/src/MobilizerContainer.cpp  
    ${PROJECT_SOURCE_DIR}/src/MoleculeContainer.cpp  
    ${PROJECT_SOURCE_DIR}/src/DensityContainer.cpp  
    ${PROJECT_SOURCE_DIR}/src/BasePairContainer.cpp  
    ${PROJECT_SOURCE_DIR}/src/ContactContainer.cpp  
    ${PROJECT_SOURCE_DIR}/src/ConstraintContainer.cpp  
    ${PROJECT_SOURCE_DIR}/src/Repel.cpp  
    ${PROJECT_SOURCE_DIR}/src/MonoAtoms.cpp  
    ${PROJECT_SOURCE_DIR}/src/ParameterReader.cpp  
    ${PROJECT_SOURCE_DIR}/src/BaseInteractionParameterReader.cpp  
    ${PROJECT_SOURCE_DIR}/src/TetherForce.cpp  
    ${PROJECT_SOURCE_DIR}/src/BiopolymerClassTwoTransformForces.cpp  
    ${PROJECT_SOURCE_DIR}/src/WaterDroplet.cpp  
    ${PROJECT_SOURCE_DIR}/src/DensityMap.cpp  
    ${PROJECT_SOURCE_DIR}/src/DensityForce.cpp
    ${PROJECT_SOURCE_DIR}/src/ElectrostaticPotentialGridForce.cpp
    ${PROJECT_SOURCE_DIR}/src/MoleculeContainer.cpp  
    ${PROJECT_SOURCE_DIR}/src/ErrorManager.cpp 
    ${PROJECT_SOURCE_DIR}/src/PeriodicPdbAndEnergyWriter.cpp
    )

IF (Build_NtC)
    SET(MMB_LIBRARY_SOURCE_FILES
        ${MMB_LIBRARY_SOURCE_FILES}
        # SCF manually commented these 3  out. Once Peter fixes the issue with the 2D arrays in BiopolymerClassContainer, we can put them back in. Or have a build varaible to determine whether to build these:
        ${PROJECT_SOURCE_DIR}/src/NTC_PARAMETER_READER.cpp
        ${PROJECT_SOURCE_DIR}/src/NtCClassContainer.cpp
        ${PROJECT_SOURCE_DIR}/src/NtCForces.cpp
        )
ENDIF(Build_NtC)

SET(MMB_SOURCE_FILES
    ${PROJECT_SOURCE_DIR}/src/RNABuilder.cpp  )

SET(MMB_GUI_TEST_SOURCE_FILES
    ${PROJECT_SOURCE_DIR}/src/GUI_test.cpp  )

SET(MMB_PYTHON_WRAPPER_SOURCE_FILES
    ${PROJECT_SOURCE_DIR}/src/python_wrapper.cpp)

# We'll just collect all the .h files from the include directory -- some
# might be junk but it doesn't matter since they won't actually get
# included.
FILE(GLOB MMB_HEADER_FILES "${PROJECT_SOURCE_DIR}/include/*.h")

SET(SHARED_MMB_TARGET MMBlib)
SET(STATIC_MMB_TARGET MMBlib_static)

SET(SHARED_MMB_PYTHON_WRAPPER_TARGET MMB_Python_wrapper)


# Ensure that debug libraries have "_d" appended to their names.
# CMake gets this right on Windows automatically with this definition.
IF (${CMAKE_GENERATOR} MATCHES "Visual Studio")
    SET(CMAKE_DEBUG_POSTFIX "_d" CACHE INTERNAL "" FORCE)
ENDIF (${CMAKE_GENERATOR} MATCHES "Visual Studio")

IF (UNIX AND CMAKE_BUILD_TYPE MATCHES Debug)
    SET(SHARED_MMB_TARGET ${SHARED_MMB_TARGET}_d)
    SET(STATIC_MMB_TARGET ${STATIC_MMB_TARGET}_d)
ENDIF (UNIX AND CMAKE_BUILD_TYPE MATCHES Debug)

IF (BUILD_MMB_SHARED_LIB)
    # Create shared library for most of MMB
    ADD_LIBRARY(${SHARED_MMB_TARGET} SHARED
        ${MMB_LIBRARY_SOURCE_FILES}
        ${MMB_HEADER_FILES})
    SET_TARGET_PROPERTIES(${SHARED_MMB_TARGET} 
        PROPERTIES
        COMPILE_FLAGS "-DMMB_BUILDING_SHARED_LIBRARY"
        PROJECT_LABEL "MMBlib (dynamic)"
        SOVERSION "${MMB_MAJOR_VERSION}.${MMB_MINOR_VERSION}")
    TARGET_LINK_LIBRARIES(${SHARED_MMB_TARGET}
              ${SimTK_SHARED_LIBS_D}
              ${SimTK_SHARED_LIBS}
              ${OpenMM_SHARED_LIBS_D}
              ${OpenMM_SHARED_LIBS}
              ${SimTK_GENERAL_LIBS})
    #################
	### Add CCP4 Libs
	#################
	IF    ( Build_CCP4 )
		TARGET_LINK_LIBRARIES ( ${SHARED_MMB_TARGET} ccp4c )
	ENDIF ( Build_CCP4 )

    # Install files
    INSTALL(TARGETS ${SHARED_MMB_TARGET}
        LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
        RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
    )
    set (CMAKE_VERBOSE_MAKEFILE on)


#IF(WIN32)
#    SET( ws2_32)
#ELSEIF(APPLE)
#    SET(LIB_SUFFIX ".dylib")
#ELSE() #Linux
#    SET(LIB_SUFFIX ".so")
#ENDIF()
    #MESSAGE(${LIB64})
    #MESSAGE(${SimTK_INSTALL_DIR}/lib${LIB64}/${CMAKE_SHARED_LIBRARY_PREFIX}*${CMAKE_SHARED_LIBRARY_SUFFIX})
    FILE(GLOB SimTK_INSTALL_DIR_FILES  "${SimTK_INSTALL_DIR}/lib${LIB64}/${CMAKE_SHARED_LIBRARY_PREFIX}*${CMAKE_SHARED_LIBRARY_SUFFIX}") 
    FILE(GLOB Molmodel_INSTALL_DIR_FILES  "${SimTKmolmodel}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}*${CMAKE_SHARED_LIBRARY_SUFFIX}") 
    FILE(GLOB CMAKE_BINARY_DIR_FILES  "${CMAKE_BINARY_DIR}/${CMAKE_SHARED_LIBRARY_PREFIX}*${CMAKE_SHARED_LIBRARY_SUFFIX}") 
    FILE(GLOB OpenMM_INSTALL_DIR_FILES  "${OpenMM_INSTALL_DIR}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}*${CMAKE_SHARED_LIBRARY_SUFFIX}") 
    #MESSAGE(${SimTK_INSTALL_DIR_FILES})
    INSTALL(FILES   ${SimTK_INSTALL_DIR_FILES}               
        ${Molmodel_INSTALL_DIR_FILES}
        ${CMAKE_BINARY_DIR_FILES}
        ${OpenMM_INSTALL_DIR_FILES}
	DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
    )

    foreach (file IN LISTS ${SimTK_SHARED_LIBS})
            #MESSAGE("My file variable is "${file})
	    #INSTALL(FILES   ${SimTK_INSTALL_DIR}/lib${file}.so               
	    #	DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
	    #)
    endforeach (file)
    # This is recursively copying entire ${Simbody_LIB_DIR}${LIB64} ! not optimal. Try using files with e.g. $<TARGET_FILE_NAME:tgt> http://www.cmake.org/cmake/help/v3.0/manual/cmake-generator-expressions.7.html#manual:cmake-generator-expressions(7)
    # ideally also only get release libs
    #INSTALL(DIRECTORY ${Simbody_LIB_DIR}${LIB64}/ 
    #    ${SimTKmolmodel}/lib${LIB64} 
    #    ${CMAKE_BINARY_DIR} 
    #    ${OpenMM_INSTALL_DIR}/lib 
    #    DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
    #)


    #FILE(GLOB MMB_HEADER_FILES "${PROJECT_SOURCE_DIR}/include/*.h")
    #FILE (GLOB SimTK_LIB_GLOB  "${SimTKmolmodel}/lib${LIB64}/*")
    #INSTALL (FILES ${SimTK_LIB_GLOB} 
    #     DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
    #)
#SimTK_SHARED_LIBS
    #INSTALL(TARGETS ${SimTK_SHARED_LIBS}
	#LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib) 
    #foreach(file ${Simbody_LIB_DIR}${LIB64}/lib${SimTK_SHARED_LIBS})
    #    install (FILES ${file} DESTINATION ${CMAKE_INSTALL_PREFIX})
    #endforeach()
    
#LINK_DIRECTORIES(${Simbody_LIB_DIR}${LIB64} ${SimTKmolmodel}/lib${LIB64} ${CMAKE_BINARY_DIR} ${OpenMM_INSTALL_DIR}/lib   )

ENDIF (BUILD_MMB_SHARED_LIB)

IF (BUILD_MMB_SHARED_EXEC)
    # Link with shared (dynamic) library
    ADD_EXECUTABLE(MMB ${MMB_HEADER_FILES} ${MMB_SOURCE_FILES})
    SET_TARGET_PROPERTIES(MMB
        PROPERTIES 
        PROJECT_LABEL "MMB - Dynamic")
    TARGET_LINK_LIBRARIES(MMB 
              debug     MMBlib_d
              optimized MMBlib
              ${SimTK_SHARED_LIBS_D}
              ${SimTK_SHARED_LIBS}
              ${OpenMM_SHARED_LIBS_D}
              ${OpenMM_SHARED_LIBS}
              ${SimTK_GENERAL_LIBS})
              
    #################
	### Add CCP4 Libs
	#################
	IF    ( Build_CCP4 )
		TARGET_LINK_LIBRARIES ( MMB ccp4c )
	ENDIF ( Build_CCP4 )

    INSTALL(TARGETS MMB
        RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
    )

    INSTALL(FILES ${PROJECT_SOURCE_DIR}/include/resources/parameters.csv
        DESTINATION ${CMAKE_INSTALL_PREFIX}/share/MMB
    )
    
    INSTALL(DIRECTORY ${PROJECT_SOURCE_DIR}/scripts/
        DESTINATION ${CMAKE_INSTALL_PREFIX}/scripts
    )
    INSTALL(DIRECTORY ${PROJECT_SOURCE_DIR}/examples/
        DESTINATION ${CMAKE_INSTALL_PREFIX}/examples
    )

ENDIF (BUILD_MMB_SHARED_EXEC)

IF (BUILD_MMB_STATIC)
    # Create static library for most of MMB
    ADD_LIBRARY(${STATIC_MMB_TARGET} STATIC
        ${MMB_LIBRARY_SOURCE_FILES}
        ${MMB_HEADER_FILES})
    SET_TARGET_PROPERTIES(${STATIC_MMB_TARGET} 
        PROPERTIES
        COMPILE_FLAGS "-DMMB_BUILDING_STATIC_LIBRARY -DSimTK_USE_STATIC_LIBRARIES"
        PROJECT_LABEL "MMBlib (static)")
    TARGET_LINK_LIBRARIES(${STATIC_MMB_TARGET}
              ${SimTK_STATIC_LIBS_D}
              ${SimTK_STATIC_LIBS}
              ${OpenMM_SHARED_LIBS_D}
              ${OpenMM_SHARED_LIBS}
              ${SimTK_GENERAL_LIBS})
              
    #################
	### Add CCP4 Libs
	#################
	IF    ( Build_CCP4 )
		TARGET_LINK_LIBRARIES ( ${STATIC_MMB_TARGET} ccp4c )
	ENDIF ( Build_CCP4 )

    # Link with static library
    ADD_EXECUTABLE(MMB_static ${MMB_HEADER_FILES} ${MMB_SOURCE_FILES})
    SET_TARGET_PROPERTIES(MMB_static
        PROPERTIES
        COMPILE_FLAGS "-DSimTK_USE_STATIC_LIBRARIES -DMMB_USE_STATIC_LIBRARIES"
        PROJECT_LABEL "MMB - Static")
    TARGET_LINK_LIBRARIES(MMB_static
              debug     MMBlib_static_d
              optimized MMBlib_static
              ${SimTK_STATIC_LIBS_D}
              ${SimTK_STATIC_LIBS}
              ${OpenMM_SHARED_LIBS_D}
              ${OpenMM_SHARED_LIBS}
              ${SimTK_GENERAL_LIBS})
              
    #################
	### Add CCP4 Libs
	#################
	IF    ( Build_CCP4 )
		TARGET_LINK_LIBRARIES ( MMB_static ccp4c )
	ENDIF ( Build_CCP4 )

    # Install files
    INSTALL(TARGETS MMB_static ${STATIC_MMB_TARGET}
        RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
        ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
    )

    INSTALL(FILES ${PROJECT_SOURCE_DIR}/include/resources/parameters.csv
        DESTINATION ${CMAKE_INSTALL_PREFIX}/share/MMB
    )

ENDIF (BUILD_MMB_STATIC)

IF (BUILD_MMB_GUI_TEST)

    # Create shared library for the Python wrapper
    ADD_LIBRARY(${SHARED_MMB_PYTHON_WRAPPER_TARGET} SHARED
        ${MMB_PYTHON_WRAPPER_SOURCE_FILES})
    SET_TARGET_PROPERTIES(${SHARED_MMB_PYTHON_WRAPPER_TARGET} 
        PROPERTIES
        COMPILE_FLAGS "-DMMB_BUILDING_SHARED_LIBRARY"
        PROJECT_LABEL "MMB_Python_wrapper (dynamic)")
    TARGET_LINK_LIBRARIES(${SHARED_MMB_PYTHON_WRAPPER_TARGET}
              ${SHARED_MMB_TARGET}
              ${SimTK_SHARED_LIBS_D}
              ${SimTK_SHARED_LIBS}
              ${OpenMM_SHARED_LIBS_D}
              ${OpenMM_SHARED_LIBS}
              ${SimTK_GENERAL_LIBS})
              
    #################
	### Add CCP4 Libs
	#################
	IF    ( Build_CCP4 )
		TARGET_LINK_LIBRARIES ( ${SHARED_MMB_PYTHON_WRAPPER_TARGET} ccp4c )
	ENDIF ( Build_CCP4 )

    # Install files
    INSTALL(TARGETS ${SHARED_MMB_PYTHON_WRAPPER_TARGET}
        RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
        LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
    )

    INSTALL(DIRECTORY ${CMAKE_SOURCE_DIR}/share/
        DESTINATION ${CMAKE_INSTALL_PREFIX}/share
        PATTERN .svn EXCLUDE
    )

ENDIF (BUILD_MMB_GUI_TEST)


