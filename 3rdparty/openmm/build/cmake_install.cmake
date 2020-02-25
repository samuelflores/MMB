# Install script for directory: /home/sam/github/openmm

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
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE FILE FILES "/home/sam/github/openmm/openmmapi/include/OpenMM.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/openmm" TYPE FILE FILES
    "/home/sam/github/openmm/olla/include/openmm/Kernel.h"
    "/home/sam/github/openmm/olla/include/openmm/KernelFactory.h"
    "/home/sam/github/openmm/olla/include/openmm/KernelImpl.h"
    "/home/sam/github/openmm/olla/include/openmm/Platform.h"
    "/home/sam/github/openmm/olla/include/openmm/PluginInitializer.h"
    "/home/sam/github/openmm/olla/include/openmm/kernels.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/AndersenThermostat.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/BAOABLangevinIntegrator.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/BrownianIntegrator.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/CMAPTorsionForce.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/CMMotionRemover.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/CompoundIntegrator.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/Context.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/CustomAngleForce.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/CustomBondForce.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/CustomCVForce.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/CustomCentroidBondForce.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/CustomCompoundBondForce.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/CustomExternalForce.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/CustomGBForce.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/CustomHbondForce.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/CustomIntegrator.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/CustomManyParticleForce.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/CustomNonbondedForce.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/CustomTorsionForce.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/Force.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/GBSAOBCForce.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/GayBerneForce.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/HarmonicAngleForce.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/HarmonicBondForce.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/Integrator.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/LangevinIntegrator.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/LocalEnergyMinimizer.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/MonteCarloAnisotropicBarostat.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/MonteCarloBarostat.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/MonteCarloMembraneBarostat.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/NonbondedForce.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/NoseHooverChain.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/NoseHooverIntegrator.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/OpenMMException.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/PeriodicTorsionForce.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/RBTorsionForce.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/RMSDForce.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/State.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/System.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/TabulatedFunction.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/Units.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/VariableLangevinIntegrator.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/VariableVerletIntegrator.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/Vec3.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/VerletIntegrator.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/VirtualSite.h"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/openmm/internal" TYPE FILE FILES
    "/home/sam/github/openmm/openmmapi/include/openmm/internal/AndersenThermostatImpl.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/internal/AssertionUtilities.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/internal/CMAPTorsionForceImpl.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/internal/CMMotionRemoverImpl.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/internal/CompiledExpressionSet.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/internal/ContextImpl.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/internal/CustomAngleForceImpl.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/internal/CustomBondForceImpl.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/internal/CustomCVForceImpl.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/internal/CustomCentroidBondForceImpl.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/internal/CustomCompoundBondForceImpl.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/internal/CustomExternalForceImpl.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/internal/CustomGBForceImpl.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/internal/CustomHbondForceImpl.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/internal/CustomIntegratorUtilities.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/internal/CustomManyParticleForceImpl.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/internal/CustomNonbondedForceImpl.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/internal/CustomTorsionForceImpl.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/internal/ForceImpl.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/internal/GBSAOBCForceImpl.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/internal/GayBerneForceImpl.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/internal/HarmonicAngleForceImpl.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/internal/HarmonicBondForceImpl.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/internal/MSVC_erfc.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/internal/MonteCarloAnisotropicBarostatImpl.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/internal/MonteCarloBarostatImpl.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/internal/MonteCarloMembraneBarostatImpl.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/internal/NonbondedForceImpl.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/internal/OSRngSeed.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/internal/PeriodicTorsionForceImpl.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/internal/RBTorsionForceImpl.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/internal/RMSDForceImpl.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/internal/SplineFitter.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/internal/ThreadPool.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/internal/VectorExpression.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/internal/hardware.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/internal/timer.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/internal/vectorize.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/internal/vectorize8.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/internal/vectorize_neon.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/internal/vectorize_pnacl.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/internal/vectorize_ppc.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/internal/vectorize_sse.h"
    "/home/sam/github/openmm/openmmapi/include/openmm/internal/windowsExport.h"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/openmm/reference" TYPE FILE FILES
    "/home/sam/github/openmm/platforms/reference/include/ObcParameters.h"
    "/home/sam/github/openmm/platforms/reference/include/RealVec.h"
    "/home/sam/github/openmm/platforms/reference/include/ReferenceAndersenThermostat.h"
    "/home/sam/github/openmm/platforms/reference/include/ReferenceAngleBondIxn.h"
    "/home/sam/github/openmm/platforms/reference/include/ReferenceBAOABDynamics.h"
    "/home/sam/github/openmm/platforms/reference/include/ReferenceBondForce.h"
    "/home/sam/github/openmm/platforms/reference/include/ReferenceBondIxn.h"
    "/home/sam/github/openmm/platforms/reference/include/ReferenceBrownianDynamics.h"
    "/home/sam/github/openmm/platforms/reference/include/ReferenceCCMAAlgorithm.h"
    "/home/sam/github/openmm/platforms/reference/include/ReferenceCMAPTorsionIxn.h"
    "/home/sam/github/openmm/platforms/reference/include/ReferenceConstraintAlgorithm.h"
    "/home/sam/github/openmm/platforms/reference/include/ReferenceConstraints.h"
    "/home/sam/github/openmm/platforms/reference/include/ReferenceCustomAngleIxn.h"
    "/home/sam/github/openmm/platforms/reference/include/ReferenceCustomBondIxn.h"
    "/home/sam/github/openmm/platforms/reference/include/ReferenceCustomCVForce.h"
    "/home/sam/github/openmm/platforms/reference/include/ReferenceCustomCentroidBondIxn.h"
    "/home/sam/github/openmm/platforms/reference/include/ReferenceCustomCompoundBondIxn.h"
    "/home/sam/github/openmm/platforms/reference/include/ReferenceCustomDynamics.h"
    "/home/sam/github/openmm/platforms/reference/include/ReferenceCustomExternalIxn.h"
    "/home/sam/github/openmm/platforms/reference/include/ReferenceCustomGBIxn.h"
    "/home/sam/github/openmm/platforms/reference/include/ReferenceCustomHbondIxn.h"
    "/home/sam/github/openmm/platforms/reference/include/ReferenceCustomManyParticleIxn.h"
    "/home/sam/github/openmm/platforms/reference/include/ReferenceCustomNonbondedIxn.h"
    "/home/sam/github/openmm/platforms/reference/include/ReferenceCustomTorsionIxn.h"
    "/home/sam/github/openmm/platforms/reference/include/ReferenceDynamics.h"
    "/home/sam/github/openmm/platforms/reference/include/ReferenceForce.h"
    "/home/sam/github/openmm/platforms/reference/include/ReferenceGayBerneForce.h"
    "/home/sam/github/openmm/platforms/reference/include/ReferenceHarmonicBondIxn.h"
    "/home/sam/github/openmm/platforms/reference/include/ReferenceKernelFactory.h"
    "/home/sam/github/openmm/platforms/reference/include/ReferenceKernels.h"
    "/home/sam/github/openmm/platforms/reference/include/ReferenceLJCoulomb14.h"
    "/home/sam/github/openmm/platforms/reference/include/ReferenceLJCoulombIxn.h"
    "/home/sam/github/openmm/platforms/reference/include/ReferenceLincsAlgorithm.h"
    "/home/sam/github/openmm/platforms/reference/include/ReferenceMonteCarloBarostat.h"
    "/home/sam/github/openmm/platforms/reference/include/ReferenceNeighborList.h"
    "/home/sam/github/openmm/platforms/reference/include/ReferenceNoseHooverChain.h"
    "/home/sam/github/openmm/platforms/reference/include/ReferenceObc.h"
    "/home/sam/github/openmm/platforms/reference/include/ReferencePME.h"
    "/home/sam/github/openmm/platforms/reference/include/ReferencePairIxn.h"
    "/home/sam/github/openmm/platforms/reference/include/ReferencePlatform.h"
    "/home/sam/github/openmm/platforms/reference/include/ReferenceProperDihedralBond.h"
    "/home/sam/github/openmm/platforms/reference/include/ReferenceRMSDForce.h"
    "/home/sam/github/openmm/platforms/reference/include/ReferenceRbDihedralBond.h"
    "/home/sam/github/openmm/platforms/reference/include/ReferenceSETTLEAlgorithm.h"
    "/home/sam/github/openmm/platforms/reference/include/ReferenceStochasticDynamics.h"
    "/home/sam/github/openmm/platforms/reference/include/ReferenceTabulatedFunction.h"
    "/home/sam/github/openmm/platforms/reference/include/ReferenceVariableStochasticDynamics.h"
    "/home/sam/github/openmm/platforms/reference/include/ReferenceVariableVerletDynamics.h"
    "/home/sam/github/openmm/platforms/reference/include/ReferenceVelocityVerletDynamics.h"
    "/home/sam/github/openmm/platforms/reference/include/ReferenceVerletDynamics.h"
    "/home/sam/github/openmm/platforms/reference/include/ReferenceVirtualSites.h"
    "/home/sam/github/openmm/platforms/reference/include/SimTKOpenMMRealType.h"
    "/home/sam/github/openmm/platforms/reference/include/SimTKOpenMMUtilities.h"
    "/home/sam/github/openmm/platforms/reference/include/fftpack.h"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/lepton" TYPE FILE FILES
    "/home/sam/github/openmm/libraries/lepton/include/lepton/CompiledExpression.h"
    "/home/sam/github/openmm/libraries/lepton/include/lepton/CustomFunction.h"
    "/home/sam/github/openmm/libraries/lepton/include/lepton/Exception.h"
    "/home/sam/github/openmm/libraries/lepton/include/lepton/ExpressionProgram.h"
    "/home/sam/github/openmm/libraries/lepton/include/lepton/ExpressionTreeNode.h"
    "/home/sam/github/openmm/libraries/lepton/include/lepton/Operation.h"
    "/home/sam/github/openmm/libraries/lepton/include/lepton/ParsedExpression.h"
    "/home/sam/github/openmm/libraries/lepton/include/lepton/Parser.h"
    "/home/sam/github/openmm/libraries/lepton/include/lepton/windowsIncludes.h"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/sfmt" TYPE FILE FILES "/home/sam/github/openmm/libraries/sfmt/include/sfmt/SFMT.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libOpenMM.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libOpenMM.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libOpenMM.so"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/home/sam/github/openmm/build/libOpenMM.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libOpenMM.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libOpenMM.so")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libOpenMM.so")
    endif()
  endif()
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/sam/github/openmm/build/wrappers/cmake_install.cmake")
  include("/home/sam/github/openmm/build/platforms/reference/tests/cmake_install.cmake")
  include("/home/sam/github/openmm/build/platforms/cpu/cmake_install.cmake")
  include("/home/sam/github/openmm/build/plugins/amoeba/cmake_install.cmake")
  include("/home/sam/github/openmm/build/plugins/rpmd/cmake_install.cmake")
  include("/home/sam/github/openmm/build/plugins/drude/cmake_install.cmake")
  include("/home/sam/github/openmm/build/serialization/cmake_install.cmake")
  include("/home/sam/github/openmm/build/wrappers/python/cmake_install.cmake")
  include("/home/sam/github/openmm/build/docs-source/cmake_install.cmake")
  include("/home/sam/github/openmm/build/tests/cmake_install.cmake")
  include("/home/sam/github/openmm/build/examples/cmake_install.cmake")

endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/home/sam/github/openmm/build/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
