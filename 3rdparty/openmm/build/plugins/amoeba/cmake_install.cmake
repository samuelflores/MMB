# Install script for directory: /home/sam/github/openmm/plugins/amoeba

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local/openmm")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE FILE FILES "/home/sam/github/openmm/plugins/amoeba/openmmapi/include/OpenMMAmoeba.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/openmm" TYPE FILE FILES
    "/home/sam/github/openmm/plugins/amoeba/openmmapi/include/openmm/AmoebaAngleForce.h"
    "/home/sam/github/openmm/plugins/amoeba/openmmapi/include/openmm/AmoebaBondForce.h"
    "/home/sam/github/openmm/plugins/amoeba/openmmapi/include/openmm/AmoebaGeneralizedKirkwoodForce.h"
    "/home/sam/github/openmm/plugins/amoeba/openmmapi/include/openmm/AmoebaInPlaneAngleForce.h"
    "/home/sam/github/openmm/plugins/amoeba/openmmapi/include/openmm/AmoebaMultipoleForce.h"
    "/home/sam/github/openmm/plugins/amoeba/openmmapi/include/openmm/AmoebaOutOfPlaneBendForce.h"
    "/home/sam/github/openmm/plugins/amoeba/openmmapi/include/openmm/AmoebaPiTorsionForce.h"
    "/home/sam/github/openmm/plugins/amoeba/openmmapi/include/openmm/AmoebaStretchBendForce.h"
    "/home/sam/github/openmm/plugins/amoeba/openmmapi/include/openmm/AmoebaTorsionTorsionForce.h"
    "/home/sam/github/openmm/plugins/amoeba/openmmapi/include/openmm/AmoebaVdwForce.h"
    "/home/sam/github/openmm/plugins/amoeba/openmmapi/include/openmm/AmoebaWcaDispersionForce.h"
    "/home/sam/github/openmm/plugins/amoeba/openmmapi/include/openmm/HippoNonbondedForce.h"
    "/home/sam/github/openmm/plugins/amoeba/openmmapi/include/openmm/amoebaKernels.h"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/openmm/internal" TYPE FILE FILES
    "/home/sam/github/openmm/plugins/amoeba/openmmapi/include/openmm/internal/AmoebaAngleForceImpl.h"
    "/home/sam/github/openmm/plugins/amoeba/openmmapi/include/openmm/internal/AmoebaBondForceImpl.h"
    "/home/sam/github/openmm/plugins/amoeba/openmmapi/include/openmm/internal/AmoebaGeneralizedKirkwoodForceImpl.h"
    "/home/sam/github/openmm/plugins/amoeba/openmmapi/include/openmm/internal/AmoebaInPlaneAngleForceImpl.h"
    "/home/sam/github/openmm/plugins/amoeba/openmmapi/include/openmm/internal/AmoebaMultipoleForceImpl.h"
    "/home/sam/github/openmm/plugins/amoeba/openmmapi/include/openmm/internal/AmoebaOutOfPlaneBendForceImpl.h"
    "/home/sam/github/openmm/plugins/amoeba/openmmapi/include/openmm/internal/AmoebaPiTorsionForceImpl.h"
    "/home/sam/github/openmm/plugins/amoeba/openmmapi/include/openmm/internal/AmoebaStretchBendForceImpl.h"
    "/home/sam/github/openmm/plugins/amoeba/openmmapi/include/openmm/internal/AmoebaTorsionTorsionForceImpl.h"
    "/home/sam/github/openmm/plugins/amoeba/openmmapi/include/openmm/internal/AmoebaVdwForceImpl.h"
    "/home/sam/github/openmm/plugins/amoeba/openmmapi/include/openmm/internal/AmoebaWcaDispersionForceImpl.h"
    "/home/sam/github/openmm/plugins/amoeba/openmmapi/include/openmm/internal/HippoNonbondedForceImpl.h"
    "/home/sam/github/openmm/plugins/amoeba/openmmapi/include/openmm/internal/windowsExportAmoeba.h"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libOpenMMAmoeba.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libOpenMMAmoeba.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libOpenMMAmoeba.so"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/home/sam/github/openmm/build/libOpenMMAmoeba.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libOpenMMAmoeba.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libOpenMMAmoeba.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libOpenMMAmoeba.so"
         OLD_RPATH "/home/sam/github/openmm/build:"
         NEW_RPATH "")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libOpenMMAmoeba.so")
    endif()
  endif()
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/sam/github/openmm/build/plugins/amoeba/wrappers/cmake_install.cmake")
  include("/home/sam/github/openmm/build/plugins/amoeba/platforms/reference/cmake_install.cmake")
  include("/home/sam/github/openmm/build/plugins/amoeba/serialization/tests/cmake_install.cmake")

endif()

