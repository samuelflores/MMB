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
include platforms/reference/tests/CMakeFiles/TestReferenceBAOABLangevinIntegrator.dir/depend.make

# Include the progress variables for this target.
include platforms/reference/tests/CMakeFiles/TestReferenceBAOABLangevinIntegrator.dir/progress.make

# Include the compile flags for this target's objects.
include platforms/reference/tests/CMakeFiles/TestReferenceBAOABLangevinIntegrator.dir/flags.make

platforms/reference/tests/CMakeFiles/TestReferenceBAOABLangevinIntegrator.dir/TestReferenceBAOABLangevinIntegrator.cpp.o: platforms/reference/tests/CMakeFiles/TestReferenceBAOABLangevinIntegrator.dir/flags.make
platforms/reference/tests/CMakeFiles/TestReferenceBAOABLangevinIntegrator.dir/TestReferenceBAOABLangevinIntegrator.cpp.o: ../platforms/reference/tests/TestReferenceBAOABLangevinIntegrator.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/sam/github/openmm/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object platforms/reference/tests/CMakeFiles/TestReferenceBAOABLangevinIntegrator.dir/TestReferenceBAOABLangevinIntegrator.cpp.o"
	cd /home/sam/github/openmm/build/platforms/reference/tests && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/TestReferenceBAOABLangevinIntegrator.dir/TestReferenceBAOABLangevinIntegrator.cpp.o -c /home/sam/github/openmm/platforms/reference/tests/TestReferenceBAOABLangevinIntegrator.cpp

platforms/reference/tests/CMakeFiles/TestReferenceBAOABLangevinIntegrator.dir/TestReferenceBAOABLangevinIntegrator.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/TestReferenceBAOABLangevinIntegrator.dir/TestReferenceBAOABLangevinIntegrator.cpp.i"
	cd /home/sam/github/openmm/build/platforms/reference/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/sam/github/openmm/platforms/reference/tests/TestReferenceBAOABLangevinIntegrator.cpp > CMakeFiles/TestReferenceBAOABLangevinIntegrator.dir/TestReferenceBAOABLangevinIntegrator.cpp.i

platforms/reference/tests/CMakeFiles/TestReferenceBAOABLangevinIntegrator.dir/TestReferenceBAOABLangevinIntegrator.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/TestReferenceBAOABLangevinIntegrator.dir/TestReferenceBAOABLangevinIntegrator.cpp.s"
	cd /home/sam/github/openmm/build/platforms/reference/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/sam/github/openmm/platforms/reference/tests/TestReferenceBAOABLangevinIntegrator.cpp -o CMakeFiles/TestReferenceBAOABLangevinIntegrator.dir/TestReferenceBAOABLangevinIntegrator.cpp.s

platforms/reference/tests/CMakeFiles/TestReferenceBAOABLangevinIntegrator.dir/TestReferenceBAOABLangevinIntegrator.cpp.o.requires:

.PHONY : platforms/reference/tests/CMakeFiles/TestReferenceBAOABLangevinIntegrator.dir/TestReferenceBAOABLangevinIntegrator.cpp.o.requires

platforms/reference/tests/CMakeFiles/TestReferenceBAOABLangevinIntegrator.dir/TestReferenceBAOABLangevinIntegrator.cpp.o.provides: platforms/reference/tests/CMakeFiles/TestReferenceBAOABLangevinIntegrator.dir/TestReferenceBAOABLangevinIntegrator.cpp.o.requires
	$(MAKE) -f platforms/reference/tests/CMakeFiles/TestReferenceBAOABLangevinIntegrator.dir/build.make platforms/reference/tests/CMakeFiles/TestReferenceBAOABLangevinIntegrator.dir/TestReferenceBAOABLangevinIntegrator.cpp.o.provides.build
.PHONY : platforms/reference/tests/CMakeFiles/TestReferenceBAOABLangevinIntegrator.dir/TestReferenceBAOABLangevinIntegrator.cpp.o.provides

platforms/reference/tests/CMakeFiles/TestReferenceBAOABLangevinIntegrator.dir/TestReferenceBAOABLangevinIntegrator.cpp.o.provides.build: platforms/reference/tests/CMakeFiles/TestReferenceBAOABLangevinIntegrator.dir/TestReferenceBAOABLangevinIntegrator.cpp.o


# Object files for target TestReferenceBAOABLangevinIntegrator
TestReferenceBAOABLangevinIntegrator_OBJECTS = \
"CMakeFiles/TestReferenceBAOABLangevinIntegrator.dir/TestReferenceBAOABLangevinIntegrator.cpp.o"

# External object files for target TestReferenceBAOABLangevinIntegrator
TestReferenceBAOABLangevinIntegrator_EXTERNAL_OBJECTS =

TestReferenceBAOABLangevinIntegrator: platforms/reference/tests/CMakeFiles/TestReferenceBAOABLangevinIntegrator.dir/TestReferenceBAOABLangevinIntegrator.cpp.o
TestReferenceBAOABLangevinIntegrator: platforms/reference/tests/CMakeFiles/TestReferenceBAOABLangevinIntegrator.dir/build.make
TestReferenceBAOABLangevinIntegrator: libOpenMM.so
TestReferenceBAOABLangevinIntegrator: /usr/lib/x86_64-linux-gnu/libdl.so
TestReferenceBAOABLangevinIntegrator: platforms/reference/tests/CMakeFiles/TestReferenceBAOABLangevinIntegrator.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/sam/github/openmm/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../../../TestReferenceBAOABLangevinIntegrator"
	cd /home/sam/github/openmm/build/platforms/reference/tests && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/TestReferenceBAOABLangevinIntegrator.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
platforms/reference/tests/CMakeFiles/TestReferenceBAOABLangevinIntegrator.dir/build: TestReferenceBAOABLangevinIntegrator

.PHONY : platforms/reference/tests/CMakeFiles/TestReferenceBAOABLangevinIntegrator.dir/build

platforms/reference/tests/CMakeFiles/TestReferenceBAOABLangevinIntegrator.dir/requires: platforms/reference/tests/CMakeFiles/TestReferenceBAOABLangevinIntegrator.dir/TestReferenceBAOABLangevinIntegrator.cpp.o.requires

.PHONY : platforms/reference/tests/CMakeFiles/TestReferenceBAOABLangevinIntegrator.dir/requires

platforms/reference/tests/CMakeFiles/TestReferenceBAOABLangevinIntegrator.dir/clean:
	cd /home/sam/github/openmm/build/platforms/reference/tests && $(CMAKE_COMMAND) -P CMakeFiles/TestReferenceBAOABLangevinIntegrator.dir/cmake_clean.cmake
.PHONY : platforms/reference/tests/CMakeFiles/TestReferenceBAOABLangevinIntegrator.dir/clean

platforms/reference/tests/CMakeFiles/TestReferenceBAOABLangevinIntegrator.dir/depend:
	cd /home/sam/github/openmm/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/sam/github/openmm /home/sam/github/openmm/platforms/reference/tests /home/sam/github/openmm/build /home/sam/github/openmm/build/platforms/reference/tests /home/sam/github/openmm/build/platforms/reference/tests/CMakeFiles/TestReferenceBAOABLangevinIntegrator.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : platforms/reference/tests/CMakeFiles/TestReferenceBAOABLangevinIntegrator.dir/depend

