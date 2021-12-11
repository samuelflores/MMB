cmake_minimum_required(VERSION 3.2)

project(MMB LANGUAGES CXX)

set(MMB_MAJOR_VERSION 3)
set(MMB_MINOR_VERSION 4)
set(MMB_PATCH_VERSION 0)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(ENABLE_NTC "Enable NtC code" ON)
option(ENABLE_LEPTON "Enable SimTK Lepton math parser" ON)
option(BUILD_SHARED "Build shared binaries" ON)
option(BUILD_STATIC "Build static binaries" OFF)
option(BUILD_PYTHON_WRAPPER "Build Python wrapper around libMMB API" OFF)
option(ASSUME_VERSIONED_SIMBODY OFF)

set(
    MMB_VERSION
    "${MMB_MAJOR_VERSION}.${MMB_MINOR_VERSION}${MMB_PATCH_VERSION}"
)

set(
    MMB_SONAME_VERSION
    "${MMB_MAJOR_VERSION}.${MMB_MINOR_VERSION}"
)

### Locate dependencies ###

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake_modules")

find_package(Simbody REQUIRED)
find_package(Molmodel REQUIRED)
find_package(OpenMM REQUIRED)
find_package(ZLIB REQUIRED)

### Setup include directories ###

include_directories(SYSTEM ${Simbody_INCLUDE_DIR})
include_directories(SYSTEM ${Molmodel_INCLUDE_DIR})
include_directories(SYSTEM ${OpenMM_INCLUDE_DIR})
include_directories(SYSTEM ${ZLIB_INCLUDE_DIRS})

if (NOT "${GEMMI_DIR}" STREQUAL "")
    include_directories(SYSTEM "${GEMMI_DIR}/include")
endif ()

if (NOT "${SEQAN_DIR}" STREQUAL "")
    include_directories(SYSTEM "${SEQAN_DIR}/include")
endif ()

if (ENABLE_LEPTON)
    if ("${LEPTON_DIR}" STREQUAL "")
        set(LEPTON_DIR "3rdparty/Lepton1.3")
    endif ()
    include_directories(SYSTEM "${LEPTON_DIR}/include")
endif ()

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include")

### Setup compiler options ###

if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
    if (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 9.0.0)
        add_definitions(-DUSE_MMB_CONSTEXPR)
    endif ()
elseif (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
    if (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 3.4)
        add_definitions(-DUSE_MMB_CONSTEXPR)
    endif ()
endif ()

if (WIN32)
    add_definitions(-DNOMINMAX -D_USE_MATH_DEFINES -DWIN32_LEAN_AND_MEAN)

    if (MSVC)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj /MP")
    endif ()
elseif (UNIX)
    if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU" OR
        ${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fvisibility=hidden -Wall -pedantic")
    endif ()
endif ()

if (ENABLE_NTC)
    add_definitions(-DNTC_ENABLED)
endif ()

# Sanity checks
if ((NOT BUILD_SHARED) AND (NOT BUILD_STATIC))
    message(FATAL_ERROR "Neither shared not static build is enabled. There is nothing to build")
endif ()

if (BUILD_SHARED)
    if (NOT Molmodel_LIBRARIES)
        message(FATAL_ERROR "Shared build enabled but shared Molmodel libraries were not found")
    endif ()
    if (NOT OpenMM_LIBRARY)
        message(FATAL_ERROR "Shared build enabled but shared OpenMM libraries were not found")
    endif ()
endif ()

if (BUILD_STATIC)
    if (NOT Molmodel_STATIC_LIBRARIES)
        message(FATAL_ERROR "Static build enabled but static Molmodel libraries were not found")
    endif ()
    if (NOT OpenMM_STATIC_LIBRARY)
        message(FATAL_ERROR "Static build enabled but static OpenMM libraries were not found")
    endif ()
endif ()

if (ASSUME_VERSIONED_SIMBODY)
    set(SIMBODY_VERTAG _${Simbody_VERSION})
endif ()

set(CMAKE_DEBUG_POSTFIX "_d")

set(
    LIBMMB_SOURCE_FILES
    src/Utils.cpp
    src/BiopolymerClass.cpp
    src/UnitCellParameters.cpp
    src/AtomSpringContainer.cpp
    src/DisplacementContainer.cpp
    src/MobilizerContainer.cpp
    src/MoleculeContainer.cpp
    src/DensityContainer.cpp
    src/BasePairContainer.cpp
    src/ContactContainer.cpp
    src/ConstraintContainer.cpp
    src/Repel.cpp
    src/MonoAtoms.cpp
    src/Spiral.cpp
    src/ParameterReader.cpp
    src/BaseInteractionParameterReader.cpp
    src/TetherForce.cpp
    src/BiopolymerClassTwoTransformForces.cpp
    src/WaterDroplet.cpp
    src/DensityMap.cpp
    src/DensityForce.cpp
    src/ElectrostaticPotentialGridForce.cpp
    src/MoleculeContainer.cpp
    src/MMBLogger.cpp
    src/PeriodicPdbAndEnergyWriter.cpp
    src/CifOutput.cpp
    src/ProgressWriter.cpp
)

if (ENABLE_NTC)
    set(
        LIBMMB_SOURCE_FILES
        ${LIBMMB_SOURCE_FILES}
        src/NTC_PARAMETER_READER.cpp
        src/NtCClassContainer.cpp
        src/NtCForces.cpp
    )
endif ()

if (ENABLE_LEPTON)
    set(
        LIBMMB_SOURCE_FILES
        ${LIBMMB_SOURCE_FILES}
        ${LEPTON_DIR}/src/CompiledExpression.cpp
        ${LEPTON_DIR}/src/ExpressionProgram.cpp
        ${LEPTON_DIR}/src/ExpressionTreeNode.cpp
        ${LEPTON_DIR}/src/Operation.cpp
        ${LEPTON_DIR}/src/ParsedExpression.cpp
        ${LEPTON_DIR}/src/Parser.cpp
    )
endif ()

set(
    MMB_SOURCE_FILES
    src/RNABuilder.cpp
)

if (BUILD_PYTHON_WRAPPER)
    set(
        PYTHON_WRAPPER_SOURCE_FILES
        src/python_wrapper.cpp
    )
endif ()

### Build the binaries ###

set(LIBMMB_SHARED_TARGET MMB)
set(LIBMMB_STATIC_TARGET MMB_static)
set(MMB_SHARED_TARGET MMB_EXEC)
set(MMB_STATIC_TARGET MMB_EXEC_static)

# Fixup Simbody library names in case we are building against versioned Simbody libraries
if (ASSUME_VERSIONED_SIMBODY)
    foreach(SL ${Simbody_LIBRARIES})
        list(APPEND SIMBODY_FIXED_LIBS ${SL}${SIMBODY_VERTAG})
        list(APPEND SIMBODY_FIXED_STATIC_LIBS ${SL}${SIMBODY_VERTAG}_static)
    endforeach()
else ()
    set(SIMBODY_FIXED_LIBS ${Simbody_LIBRARIES})
    set(SIMBODY_FIXED_STATIC_LIBS ${Simbody_STATIC_LIBRARIES})
endif ()

if (BUILD_PYTHON_WRAPPER)
    set(LIBMMB_PYTHON_WRAPPER_SHARED_TARGET MMB_Python_wrapper)
    set(LIBMMB_PYTHON_WRAPPER_STATIC_TARGET MMB_Python_wrapper_static)
endif ()

## Build libMMB
if (BUILD_SHARED)
    add_library(
        ${LIBMMB_SHARED_TARGET} SHARED
        ${LIBMMB_SOURCE_FILES}
    )
    set_target_properties(
        ${LIBMMB_SHARED_TARGET}
        PROPERTIES
            PROJECT_LABEL "MMB (shared)"
            SOVERSION ${MMB_SONAME_VERSION}
    )
    target_compile_definitions(${LIBMMB_SHARED_TARGET} PRIVATE -DMMB_BUILDING_SHARED_LIBRARY)

    target_link_libraries(
        ${LIBMMB_SHARED_TARGET}
        PUBLIC  ${SIMBODY_FIXED_LIBS}
        PUBLIC  ${Molmodel_LIBRARIES}
        PRIVATE ${OpenMM_LIBRARY}
        PRIVATE ${ZLIB_LIBRARIES}
    )

    if (BUILD_PYTHON_WRAPPER)
        add_library(
            ${LIBMMB_PYTHON_WRAPPER_SHARED_TARGET}
            ${PYTHON_WRAPPER_SOURCE_FILES}
        )
        set_target_properties(
            ${LIBMMB_PYTHON_WRAPPER_SHARED_TARGET}
            PROPERTIES
                PROJECT_LABEL "libMMB Python wrapper (shared)"
                SOVERSION ${MMB_SONAME_VERSION}
        )

        target_link_libraries(
            ${LIBMMB_PYTHON_WRAPPER_SHARED_TARGET}
            PRIVATE ${LIBMMB_SHARED_TARGET}
        )
    endif ()

    add_executable(
        ${MMB_SHARED_TARGET}
        ${MMB_SOURCE_FILES}
    )
    set_target_properties(
        ${MMB_SHARED_TARGET}
        PROPERTIES
            OUTPUT_NAME MMB
    )

    target_link_libraries(
        ${MMB_SHARED_TARGET}
        PRIVATE ${LIBMMB_SHARED_TARGET}
    )
endif ()

if (BUILD_STATIC)
    add_library(
        ${LIBMMB_STATIC_TARGET} STATIC
        ${LIBMMB_SOURCE_FILES}
    )
    set_target_properties(
        ${LIBMMB_STATIC_TARGET}
        PROPERTIES
            PROJECT_LABEL "MMB (static)"
    )
    target_compile_definitions(${LIBMMB_STATIC_TARGET} PRIVATE -DMMB_BUILDING_STATIC_LIBRARY -DSimTK_USE_STATIC_LIBRARIES)

    target_link_libraries(
        ${LIBMMB_STATIC_TARGET}
        PUBLIC  ${SIMBODY_FIXED_STATIC_LIBS}
        PUBLIC  ${Molmodel_STATIC_LIBRARIES}
        PRIVATE ${OpenMM_STATIC_LIBRARY}
        PRIVATE ${ZLIB_LIBRARIES}
    )

    if (BUILD_PYTHON_WRAPPER)
        add_library(
            ${LIBMMB_PYTHON_WRAPPER_STATIC_TARGET}
            ${PYTHON_WRAPPER_SOURCE_FILES}
        )
        set_target_properties(
            ${LIBMMB_PYTHON_WRAPPER_STATIC_TARGET}
            PROPERTIES
                PROJECT_LABEL "libMMB Python wrapper (static)"
        )

        target_link_libraries(
            ${LIBMMB_PYTHON_WRAPPER_STATIC_TARGET}
            PRIVATE ${LIBMMB_STATIC_TARGET}
        )
    endif ()

    add_executable(
        ${MMB_STATIC_TARGET}
        ${MMB_SOURCE_FILES}
    )
    set_target_properties(
        ${MMB_STATIC_TARGET}
        PROPERTIES
            OUTPUT_NAME MMB_static
    )
    target_compile_definitions(${MMB_STATIC_TARGET} PRIVATE -DMMB_USE_STATIC_LIBRARY -DSimTK_USE_STATIC_LIBRARIES)

    target_link_libraries(
        ${MMB_STATIC_TARGET}
        PRIVATE ${LIBMMB_STATIC_TARGET}
    )

endif()

### Install the project ###

if (UNIX)
    include(GNUInstallDirs)
elseif (WIN32)
    set(${CMAKE_INSTALL_LIBDIR} "lib")
    set(${CMAKE_INSTALL_INCLUDEDIR} "include")
    set(${CMAKE_INSTALL_BINDIR} "bin")
endif ()

## Install libMMB
if (BUILD_SHARED)
    install(
        TARGETS ${LIBMMB_SHARED_TARGET} EXPORT MmbTargets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )
endif ()

if (BUILD_STATIC)
    install(
        TARGETS ${LIBMMB_STATIC_TARGET} EXPORT MmbTargets
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
endif ()

## Install MMB executable
if (BUILD_SHARED)
    install(
        TARGETS ${MMB_SHARED_TARGET}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )
endif ()

if (BUILD_STATIC)
    install(
        TARGETS ${MMB_STATIC_TARGET}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )
endif ()

if (WIN32)
    set(MMB_CMAKE_DIR cmake)
else ()
    set(MMB_CMAKE_DIR ${CMAKE_INSTALL_LIBDIR}/cmake/mmb/)
endif ()

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/cmake/MmbConfigVersion.cmake
    VERSION "${MMB_VERSION}"
    COMPATIBILITY SameMinorVersion
)
install(
    FILES ${CMAKE_CURRENT_BINARY_DIR}/cmake/MmbConfigVersion.cmake
    DESTINATION ${MMB_CMAKE_DIR}
)

install(EXPORT MmbTargets DESTINATION "${MMB_CMAKE_DIR}")

configure_file(
    ${CMAKE_SOURCE_DIR}/cmake/pkgconfig/mmb.pc.in
    ${CMAKE_CURRENT_BINARY_DIR}/cmake/pkgconfig/mmb.pc @ONLY
)
install(
    FILES ${CMAKE_CURRENT_BINARY_DIR}/cmake/pkgconfig/mmb.pc
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig/
)
