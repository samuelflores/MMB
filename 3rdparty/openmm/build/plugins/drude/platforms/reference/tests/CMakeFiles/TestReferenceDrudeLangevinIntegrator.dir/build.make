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
include plugins/drude/platforms/reference/tests/CMakeFiles/TestReferenceDrudeLangevinIntegrator.dir/depend.make

# Include the progress variables for this target.
include plugins/drude/platforms/reference/tests/CMakeFiles/TestReferenceDrudeLangevinIntegrator.dir/progress.make

# Include the compile flags for this target's objects.
include plugins/drude/platforms/reference/tests/CMakeFiles/TestReferenceDrudeLangevinIntegrator.dir/flags.make

plugins/drude/platforms/reference/tests/CMakeFiles/TestReferenceDrudeLangevinIntegrator.dir/TestReferenceDrudeLangevinIntegrator.cpp.o: plugins/drude/platforms/reference/tests/CMakeFiles/TestReferenceDrudeLangevinIntegrator.dir/flags.make
plugins/drude/platforms/reference/tests/CMakeFiles/TestReferenceDrudeLangevinIntegrator.dir/TestReferenceDrudeLangevinIntegrator.cpp.o: ../plugins/drude/platforms/reference/tests/TestReferenceDrudeLangevinIntegrator.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/sam/github/openmm/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object plugins/drude/platforms/reference/tests/CMakeFiles/TestReferenceDrudeLangevinIntegrator.dir/TestReferenceDrudeLangevinIntegrator.cpp.o"
	cd /home/sam/github/openmm/build/plugins/drude/platforms/reference/tests && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/TestReferenceDrudeLangevinIntegrator.dir/TestReferenceDrudeLangevinIntegrator.cpp.o -c /home/sam/github/openmm/plugins/drude/platforms/reference/tests/TestReferenceDrudeLangevinIntegrator.cpp

plugins/drude/platforms/reference/tests/CMakeFiles/TestReferenceDrudeLangevinIntegrator.dir/TestReferenceDrudeLangevinIntegrator.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/TestReferenceDrudeLangevinIntegrator.dir/TestReferenceDrudeLangevinIntegrator.cpp.i"
	cd /home/sam/github/openmm/build/plugins/drude/platforms/reference/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/sam/github/openmm/plugins/drude/platforms/reference/tests/TestReferenceDrudeLangevinIntegrator.cpp > CMakeFiles/TestReferenceDrudeLangevinIntegrator.dir/TestReferenceDrudeLangevinIntegrator.cpp.i

plugins/drude/platforms/reference/tests/CMakeFiles/TestReferenceDrudeLangevinIntegrator.dir/TestReferenceDrudeLangevinIntegrator.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/TestReferenceDrudeLangevinIntegrator.dir/TestReferenceDrudeLangevinIntegrator.cpp.s"
	cd /home/sam/github/openmm/build/plugins/drude/platforms/reference/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/sam/github/openmm/plugins/drude/platforms/reference/tests/TestReferenceDrudeLangevinIntegrator.cpp -o CMakeFiles/TestReferenceDrudeLangevinIntegrator.dir/TestReferenceDrudeLangevinIntegrator.cpp.s

plugins/drude/platforms/reference/tests/CMakeFiles/TestReferenceDrudeLangevinIntegrator.dir/TestReferenceDrudeLangevinIntegrator.cpp.o.requires:

.PHONY : plugins/drude/platforms/reference/tests/CMakeFiles/TestReferenceDrudeLangevinIntegrator.dir/TestReferenceDrudeLangevinIntegrator.cpp.o.requires

plugins/drude/platforms/reference/tests/CMakeFiles/TestReferenceDrudeLangevinIntegrator.dir/TestReferenceDrudeLangevinIntegrator.cpp.o.provides: plugins/drude/platforms/reference/tests/CMakeFiles/TestReferenceDrudeLangevinIntegrator.dir/TestReferenceDrudeLangevinIntegrator.cpp.o.requires
	$(MAKE) -f plugins/drude/platforms/reference/tests/CMakeFiles/TestReferenceDrudeLangevinIntegrator.dir/build.make plugins/drude/platforms/reference/tests/CMakeFiles/TestReferenceDrudeLangevinIntegrator.dir/TestReferenceDrudeLangevinIntegrator.cpp.o.provides.build
.PHONY : plugins/drude/platforms/reference/tests/CMakeFiles/TestReferenceDrudeLangevinIntegrator.dir/TestReferenceDrudeLangevinIntegrator.cpp.o.provides

plugins/drude/platforms/reference/tests/CMakeFiles/TestReferenceDrudeLangevinIntegrator.dir/TestReferenceDrudeLangevinIntegrator.cpp.o.provides.build: plugins/drude/platforms/reference/tests/CMakeFiles/TestReferenceDrudeLangevinIntegrator.dir/TestReferenceDrudeLangevinIntegrator.cpp.o


# Object files for target TestReferenceDrudeLangevinIntegrator
TestReferenceDrudeLangevinIntegrator_OBJECTS = \
"CMakeFiles/TestReferenceDrudeLangevinIntegrator.dir/TestReferenceDrudeLangevinIntegrator.cpp.o"

# External object files for target TestReferenceDrudeLangevinIntegrator
TestReferenceDrudeLangevinIntegrator_EXTERNAL_OBJECTS =

TestReferenceDrudeLangevinIntegrator: plugins/drude/platforms/reference/tests/CMakeFiles/TestReferenceDrudeLangevinIntegrator.dir/TestReferenceDrudeLangevinIntegrator.cpp.o
TestReferenceDrudeLangevinIntegrator: plugins/drude/platforms/reference/tests/CMakeFiles/TestReferenceDrudeLangevinIntegrator.dir/build.make
TestReferenceDrudeLangevinIntegrator: libOpenMMDrudeReference.so
TestReferenceDrudeLangevinIntegrator: libOpenMMDrude.so
TestReferenceDrudeLangevinIntegrator: libOpenMM.so
TestReferenceDrudeLangevinIntegrator: /usr/lib/x86_64-linux-gnu/libdl.so
TestReferenceDrudeLangevinIntegrator: plugins/drude/platforms/reference/tests/CMakeFiles/TestReferenceDrudeLangevinIntegrator.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/sam/github/openmm/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../../../../../TestReferenceDrudeLangevinIntegrator"
	cd /home/sam/github/openmm/build/plugins/drude/platforms/reference/tests && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/TestReferenceDrudeLangevinIntegrator.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
plugins/drude/platforms/reference/tests/CMakeFiles/TestReferenceDrudeLangevinIntegrator.dir/build: TestReferenceDrudeLangevinIntegrator

.PHONY : plugins/drude/platforms/reference/tests/CMakeFiles/TestReferenceDrudeLangevinIntegrator.dir/build

plugins/drude/platforms/reference/tests/CMakeFiles/TestReferenceDrudeLangevinIntegrator.dir/requires: plugins/drude/platforms/reference/tests/CMakeFiles/TestReferenceDrudeLangevinIntegrator.dir/TestReferenceDrudeLangevinIntegrator.cpp.o.requires

.PHONY : plugins/drude/platforms/reference/tests/CMakeFiles/TestReferenceDrudeLangevinIntegrator.dir/requires

plugins/drude/platforms/reference/tests/CMakeFiles/TestReferenceDrudeLangevinIntegrator.dir/clean:
	cd /home/sam/github/openmm/build/plugins/drude/platforms/reference/tests && $(CMAKE_COMMAND) -P CMakeFiles/TestReferenceDrudeLangevinIntegrator.dir/cmake_clean.cmake
.PHONY : plugins/drude/platforms/reference/tests/CMakeFiles/TestReferenceDrudeLangevinIntegrator.dir/clean

plugins/drude/platforms/reference/tests/CMakeFiles/TestReferenceDrudeLangevinIntegrator.dir/depend:
	cd /home/sam/github/openmm/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/sam/github/openmm /home/sam/github/openmm/plugins/drude/platforms/reference/tests /home/sam/github/openmm/build /home/sam/github/openmm/build/plugins/drude/platforms/reference/tests /home/sam/github/openmm/build/plugins/drude/platforms/reference/tests/CMakeFiles/TestReferenceDrudeLangevinIntegrator.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : plugins/drude/platforms/reference/tests/CMakeFiles/TestReferenceDrudeLangevinIntegrator.dir/depend

