# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/sam/github/openmm

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/sam/github/openmm/build

# Include any dependencies generated for this target.
include plugins/amoeba/serialization/tests/CMakeFiles/TestSerializeAmoebaPiTorsionForce.dir/depend.make

# Include the progress variables for this target.
include plugins/amoeba/serialization/tests/CMakeFiles/TestSerializeAmoebaPiTorsionForce.dir/progress.make

# Include the compile flags for this target's objects.
include plugins/amoeba/serialization/tests/CMakeFiles/TestSerializeAmoebaPiTorsionForce.dir/flags.make

plugins/amoeba/serialization/tests/CMakeFiles/TestSerializeAmoebaPiTorsionForce.dir/TestSerializeAmoebaPiTorsionForce.cpp.o: plugins/amoeba/serialization/tests/CMakeFiles/TestSerializeAmoebaPiTorsionForce.dir/flags.make
plugins/amoeba/serialization/tests/CMakeFiles/TestSerializeAmoebaPiTorsionForce.dir/TestSerializeAmoebaPiTorsionForce.cpp.o: ../plugins/amoeba/serialization/tests/TestSerializeAmoebaPiTorsionForce.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/sam/github/openmm/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object plugins/amoeba/serialization/tests/CMakeFiles/TestSerializeAmoebaPiTorsionForce.dir/TestSerializeAmoebaPiTorsionForce.cpp.o"
	cd /home/sam/github/openmm/build/plugins/amoeba/serialization/tests && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/TestSerializeAmoebaPiTorsionForce.dir/TestSerializeAmoebaPiTorsionForce.cpp.o -c /home/sam/github/openmm/plugins/amoeba/serialization/tests/TestSerializeAmoebaPiTorsionForce.cpp

plugins/amoeba/serialization/tests/CMakeFiles/TestSerializeAmoebaPiTorsionForce.dir/TestSerializeAmoebaPiTorsionForce.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/TestSerializeAmoebaPiTorsionForce.dir/TestSerializeAmoebaPiTorsionForce.cpp.i"
	cd /home/sam/github/openmm/build/plugins/amoeba/serialization/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/sam/github/openmm/plugins/amoeba/serialization/tests/TestSerializeAmoebaPiTorsionForce.cpp > CMakeFiles/TestSerializeAmoebaPiTorsionForce.dir/TestSerializeAmoebaPiTorsionForce.cpp.i

plugins/amoeba/serialization/tests/CMakeFiles/TestSerializeAmoebaPiTorsionForce.dir/TestSerializeAmoebaPiTorsionForce.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/TestSerializeAmoebaPiTorsionForce.dir/TestSerializeAmoebaPiTorsionForce.cpp.s"
	cd /home/sam/github/openmm/build/plugins/amoeba/serialization/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/sam/github/openmm/plugins/amoeba/serialization/tests/TestSerializeAmoebaPiTorsionForce.cpp -o CMakeFiles/TestSerializeAmoebaPiTorsionForce.dir/TestSerializeAmoebaPiTorsionForce.cpp.s

plugins/amoeba/serialization/tests/CMakeFiles/TestSerializeAmoebaPiTorsionForce.dir/TestSerializeAmoebaPiTorsionForce.cpp.o.requires:

.PHONY : plugins/amoeba/serialization/tests/CMakeFiles/TestSerializeAmoebaPiTorsionForce.dir/TestSerializeAmoebaPiTorsionForce.cpp.o.requires

plugins/amoeba/serialization/tests/CMakeFiles/TestSerializeAmoebaPiTorsionForce.dir/TestSerializeAmoebaPiTorsionForce.cpp.o.provides: plugins/amoeba/serialization/tests/CMakeFiles/TestSerializeAmoebaPiTorsionForce.dir/TestSerializeAmoebaPiTorsionForce.cpp.o.requires
	$(MAKE) -f plugins/amoeba/serialization/tests/CMakeFiles/TestSerializeAmoebaPiTorsionForce.dir/build.make plugins/amoeba/serialization/tests/CMakeFiles/TestSerializeAmoebaPiTorsionForce.dir/TestSerializeAmoebaPiTorsionForce.cpp.o.provides.build
.PHONY : plugins/amoeba/serialization/tests/CMakeFiles/TestSerializeAmoebaPiTorsionForce.dir/TestSerializeAmoebaPiTorsionForce.cpp.o.provides

plugins/amoeba/serialization/tests/CMakeFiles/TestSerializeAmoebaPiTorsionForce.dir/TestSerializeAmoebaPiTorsionForce.cpp.o.provides.build: plugins/amoeba/serialization/tests/CMakeFiles/TestSerializeAmoebaPiTorsionForce.dir/TestSerializeAmoebaPiTorsionForce.cpp.o


# Object files for target TestSerializeAmoebaPiTorsionForce
TestSerializeAmoebaPiTorsionForce_OBJECTS = \
"CMakeFiles/TestSerializeAmoebaPiTorsionForce.dir/TestSerializeAmoebaPiTorsionForce.cpp.o"

# External object files for target TestSerializeAmoebaPiTorsionForce
TestSerializeAmoebaPiTorsionForce_EXTERNAL_OBJECTS =

TestSerializeAmoebaPiTorsionForce: plugins/amoeba/serialization/tests/CMakeFiles/TestSerializeAmoebaPiTorsionForce.dir/TestSerializeAmoebaPiTorsionForce.cpp.o
TestSerializeAmoebaPiTorsionForce: plugins/amoeba/serialization/tests/CMakeFiles/TestSerializeAmoebaPiTorsionForce.dir/build.make
TestSerializeAmoebaPiTorsionForce: libOpenMMAmoeba.so
TestSerializeAmoebaPiTorsionForce: libOpenMM.so
TestSerializeAmoebaPiTorsionForce: /usr/lib/x86_64-linux-gnu/libdl.so
TestSerializeAmoebaPiTorsionForce: plugins/amoeba/serialization/tests/CMakeFiles/TestSerializeAmoebaPiTorsionForce.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/sam/github/openmm/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../../../../TestSerializeAmoebaPiTorsionForce"
	cd /home/sam/github/openmm/build/plugins/amoeba/serialization/tests && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/TestSerializeAmoebaPiTorsionForce.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
plugins/amoeba/serialization/tests/CMakeFiles/TestSerializeAmoebaPiTorsionForce.dir/build: TestSerializeAmoebaPiTorsionForce

.PHONY : plugins/amoeba/serialization/tests/CMakeFiles/TestSerializeAmoebaPiTorsionForce.dir/build

plugins/amoeba/serialization/tests/CMakeFiles/TestSerializeAmoebaPiTorsionForce.dir/requires: plugins/amoeba/serialization/tests/CMakeFiles/TestSerializeAmoebaPiTorsionForce.dir/TestSerializeAmoebaPiTorsionForce.cpp.o.requires

.PHONY : plugins/amoeba/serialization/tests/CMakeFiles/TestSerializeAmoebaPiTorsionForce.dir/requires

plugins/amoeba/serialization/tests/CMakeFiles/TestSerializeAmoebaPiTorsionForce.dir/clean:
	cd /home/sam/github/openmm/build/plugins/amoeba/serialization/tests && $(CMAKE_COMMAND) -P CMakeFiles/TestSerializeAmoebaPiTorsionForce.dir/cmake_clean.cmake
.PHONY : plugins/amoeba/serialization/tests/CMakeFiles/TestSerializeAmoebaPiTorsionForce.dir/clean

plugins/amoeba/serialization/tests/CMakeFiles/TestSerializeAmoebaPiTorsionForce.dir/depend:
	cd /home/sam/github/openmm/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/sam/github/openmm /home/sam/github/openmm/plugins/amoeba/serialization/tests /home/sam/github/openmm/build /home/sam/github/openmm/build/plugins/amoeba/serialization/tests /home/sam/github/openmm/build/plugins/amoeba/serialization/tests/CMakeFiles/TestSerializeAmoebaPiTorsionForce.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : plugins/amoeba/serialization/tests/CMakeFiles/TestSerializeAmoebaPiTorsionForce.dir/depend

