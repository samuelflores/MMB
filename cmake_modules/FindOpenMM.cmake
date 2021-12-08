# FindOpenMM.cmake
# Try to find OpenMM molecular dynamics API
#
# Created August 26, 2009 by Christopher Bruns
# Simbios National Center for Biomedical Computation
# Stanford University
#
# Once done this will define:
#
#   OpenMM_FOUND - Whether search for OpenMM libraries and headers succeeded.
#   OpenMM_INCLUDE_DIR - location of OpenMM.h
#   OpenMM_LIBRARY - location of shared libOpenMM.so or OpenMM.dll
#   OpenMM_DEBUG_LIBRARY - location of shared libOpenMM_d.so or OpenMM_d.dll
#   OpenMM_LIBRARIES - cmake variable with one or both of OpenMM_LIBRARY and OpenMM_DEBUG_LIBRARY
#   OpenMM_STATIC_LIBRARY - location of libOpenMM_static.a or OpenMM_static.lib
#
#
# == Using OpenMM libraries ==
#
#     find_package(OpenMM)
#     if(OpenMM_FOUND)
#         include_directories(${OpenMM_INCLUDE_DIR})
#         target_link_libraries(Foo ${OpenMM_LIBRARIES})
#     endif()
#

# Generic locations to search for OpenMM files
set(OpenMM_SEARCH_PATHS
    ${SimTK_SDK} # SimTK.org nightly builds location
    ${OpenMM_SDK}
    "/usr/local"
    $ENV{ProgramFiles}
    "C:/Program Files"
)

# Include files
find_path(OpenMM_INCLUDE_DIR "OpenMM.h"
    PATHS ${OpenMM_SEARCH_PATHS}
    PATH_SUFFIXES "include" "OpenMM/include" "openmm/include"
    DOC "Location of OpenMM include files"
)

# Shared Release library

# Look first for libraries next door to include files
find_library(OpenMM_LIBRARY NAMES OpenMM
    PATHS "${OpenMM_INCLUDE_DIR}/.."
    NO_DEFAULT_PATH
    DOC "Location of OpenMM shared library"
)
# If above search fails, search for OpenMM library in all locations
find_library(OpenMM_LIBRARY NAMES OpenMM
    PATHS ${OpenMM_SEARCH_PATHS}
    PATH_SUFFIXES "lib" "OpenMM/lib" "openmm/lib"
)

find_library(OpenMM_DEBUG_LIBRARY NAMES OpenMM_d
    PATHS "${OpenMM_INCLUDE_DIR}/.."
    DOC "Location of OpenMM debug library"
    NO_DEFAULT_PATH
)
find_library(OpenMM_DEBUG_LIBRARY NAMES OpenMM_d
    PATHS ${OpenMM_SEARCH_PATHS}
    PATH_SUFFIXES "lib" "OpenMM/lib" "openmm/lib"
)

# Static library

find_library(OpenMM_STATIC_LIBRARY NAMES OpenMM_static
    PATHS "${OpenMM_INCLUDE_DIR}/.."
    DOC "Location of OpenMM static library"
    NO_DEFAULT_PATH
)
find_library(OpenMM_STATIC_LIBRARY NAMES OpenMM_static
    PATHS ${OpenMM_SEARCH_PATHS}
    PATH_SUFFIXES "lib" "OpenMM/lib" "openmm/lib"
)

# Set composite OpenMM_LIBRARIES variable
set(LIBS "")
if(OpenMM_LIBRARY AND OpenMM_DEBUG_LIBRARY)
    set(LIBS optimized ${OpenMM_LIBRARY} debug ${OpenMM_DEBUG_LIBRARY})
elseif(OpenMM_LIBRARY)
    set(LIBS ${OpenMM_LIBRARY})
elseif(OpenMM_DEBUG_LIBRARY)
    set(LIBS ${OpenMM_DEBUG_LIBRARY})
elseif(OpenMM_STATIC_LIBRARY)
    set(LIBS ${OpenMM_STATIC_LIBRARY})
endif(OpenMM_LIBRARY AND OpenMM_DEBUG_LIBRARY)
set(OpenMM_LIBRARIES "${LIBS}" CACHE STRING "OpenMM Link Libraries" FORCE)

# CMAKE 2.4 does not have FindPackageHandleStandardArgs
# So to heck with CMAKE 2.4...
set(cmv "${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}")
if(cmv MATCHES "2.4")
else(cmv MATCHES "2.4")
    include(FindPackageHandleStandardArgs OPTIONAL)
    find_package_handle_standard_args(OpenMM DEFAULT_MSG OpenMM_INCLUDE_DIR OpenMM_LIBRARIES)
endif(cmv MATCHES "2.4")

mark_as_advanced(
    OpenMM_INCLUDE_DIR
    OpenMM_LIBRARY
    OpenMM_DEBUG_LIBRARY
    OpenMM_STATIC_LIBRARY
    OpenMM_LIBRARIES
)
