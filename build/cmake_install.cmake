# Install script for directory: /Users/Sam/svn/RNAToolbox/trunk

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
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

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/local/lib/libMMBlib.dylib")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/usr/local/lib" TYPE SHARED_LIBRARY FILES "/Users/Sam/svn/RNAToolbox/trunk/build/libMMBlib.dylib")
  if(EXISTS "$ENV{DESTDIR}/usr/local/lib/libMMBlib.dylib" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/usr/local/lib/libMMBlib.dylib")
    execute_process(COMMAND "/usr/bin/install_name_tool"
      -id "libMMBlib.dylib"
      "$ENV{DESTDIR}/usr/local/lib/libMMBlib.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/usr/local/SimTK/lib"
      "$ENV{DESTDIR}/usr/local/lib/libMMBlib.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/Sam/svn/RNAToolbox/trunk/build"
      "$ENV{DESTDIR}/usr/local/lib/libMMBlib.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/usr/local/openmm/lib"
      "$ENV{DESTDIR}/usr/local/lib/libMMBlib.dylib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/usr/local/lib/libMMBlib.dylib")
    endif()
  endif()
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/local/bin/MMB")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/usr/local/bin" TYPE EXECUTABLE FILES "/Users/Sam/svn/RNAToolbox/trunk/build/MMB")
  if(EXISTS "$ENV{DESTDIR}/usr/local/bin/MMB" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/usr/local/bin/MMB")
    execute_process(COMMAND "/usr/bin/install_name_tool"
      -change "/Users/Sam/svn/RNAToolbox/trunk/build/libMMBlib.dylib" "libMMBlib.dylib"
      "$ENV{DESTDIR}/usr/local/bin/MMB")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/usr/local/SimTK/lib"
      "$ENV{DESTDIR}/usr/local/bin/MMB")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/Sam/svn/RNAToolbox/trunk/build"
      "$ENV{DESTDIR}/usr/local/bin/MMB")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/usr/local/openmm/lib"
      "$ENV{DESTDIR}/usr/local/bin/MMB")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/usr/local/bin/MMB")
    endif()
  endif()
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/local/share/MMB/parameters.csv")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/usr/local/share/MMB" TYPE FILE FILES "/Users/Sam/svn/RNAToolbox/trunk/include/resources/parameters.csv")
endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

file(WRITE "/Users/Sam/svn/RNAToolbox/trunk/build/${CMAKE_INSTALL_MANIFEST}" "")
foreach(file ${CMAKE_INSTALL_MANIFEST_FILES})
  file(APPEND "/Users/Sam/svn/RNAToolbox/trunk/build/${CMAKE_INSTALL_MANIFEST}" "${file}\n")
endforeach()
