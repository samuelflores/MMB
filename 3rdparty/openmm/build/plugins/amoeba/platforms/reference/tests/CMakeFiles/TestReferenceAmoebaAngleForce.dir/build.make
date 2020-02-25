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
include plugins/amoeba/platforms/reference/tests/CMakeFiles/TestReferenceAmoebaAngleForce.dir/depend.make

# Include the progress variables for this target.
include plugins/amoeba/platforms/reference/tests/CMakeFiles/TestReferenceAmoebaAngleForce.dir/progress.make

# Include the compile flags for this target's objects.
include plugins/amoeba/platforms/reference/tests/CMakeFiles/TestReferenceAmoebaAngleForce.dir/flags.make

plugins/amoeba/platforms/reference/tests/CMakeFiles/TestReferenceAmoebaAngleForce.dir/TestReferenceAmoebaAngleForce.cpp.o: plugins/amoeba/platforms/reference/tests/CMakeFiles/TestReferenceAmoebaAngleForce.dir/flags.make
plugins/amoeba/platforms/reference/tests/CMakeFiles/TestReferenceAmoebaAngleForce.dir/TestReferenceAmoebaAngleForce.cpp.o: ../plugins/amoeba/platforms/reference/tests/TestReferenceAmoebaAngleForce.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/sam/github/openmm/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object plugins/amoeba/platforms/reference/tests/CMakeFiles/TestReferenceAmoebaAngleForce.dir/TestReferenceAmoebaAngleForce.cpp.o"
	cd /home/sam/github/openmm/build/plugins/amoeba/platforms/reference/tests && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/TestReferenceAmoebaAngleForce.dir/TestReferenceAmoebaAngleForce.cpp.o -c /home/sam/github/openmm/plugins/amoeba/platforms/reference/tests/TestReferenceAmoebaAngleForce.cpp

plugins/amoeba/platforms/reference/tests/CMakeFiles/TestReferenceAmoebaAngleForce.dir/TestReferenceAmoebaAngleForce.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/TestReferenceAmoebaAngleForce.dir/TestReferenceAmoebaAngleForce.cpp.i"
	cd /home/sam/github/openmm/build/plugins/amoeba/platforms/reference/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/sam/github/openmm/plugins/amoeba/platforms/reference/tests/TestReferenceAmoebaAngleForce.cpp > CMakeFiles/TestReferenceAmoebaAngleForce.dir/TestReferenceAmoebaAngleForce.cpp.i

plugins/amoeba/platforms/reference/tests/CMakeFiles/TestReferenceAmoebaAngleForce.dir/TestReferenceAmoebaAngleForce.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/TestReferenceAmoebaAngleForce.dir/TestReferenceAmoebaAngleForce.cpp.s"
	cd /home/sam/github/openmm/build/plugins/amoeba/platforms/reference/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/sam/github/openmm/plugins/amoeba/platforms/reference/tests/TestReferenceAmoebaAngleForce.cpp -o CMakeFiles/TestReferenceAmoebaAngleForce.dir/TestReferenceAmoebaAngleForce.cpp.s

plugins/amoeba/platforms/reference/tests/CMakeFiles/TestReferenceAmoebaAngleForce.dir/TestReferenceAmoebaAngleForce.cpp.o.requires:

.PHONY : plugins/amoeba/platforms/reference/tests/CMakeFiles/TestReferenceAmoebaAngleForce.dir/TestReferenceAmoebaAngleForce.cpp.o.requires

plugins/amoeba/platforms/reference/tests/CMakeFiles/TestReferenceAmoebaAngleForce.dir/TestReferenceAmoebaAngleForce.cpp.o.provides: plugins/amoeba/platforms/reference/tests/CMakeFiles/TestReferenceAmoebaAngleForce.dir/TestReferenceAmoebaAngleForce.cpp.o.requires
	$(MAKE) -f plugins/amoeba/platforms/reference/tests/CMakeFiles/TestReferenceAmoebaAngleForce.dir/build.make plugins/amoeba/platforms/reference/tests/CMakeFiles/TestReferenceAmoebaAngleForce.dir/TestReferenceAmoebaAngleForce.cpp.o.provides.build
.PHONY : plugins/amoeba/platforms/reference/tests/CMakeFiles/TestReferenceAmoebaAngleForce.dir/TestReferenceAmoebaAngleForce.cpp.o.provides

plugins/amoeba/platforms/reference/tests/CMakeFiles/TestReferenceAmoebaAngleForce.dir/TestReferenceAmoebaAngleForce.cpp.o.provides.build: plugins/amoeba/platforms/reference/tests/CMakeFiles/TestReferenceAmoebaAngleForce.dir/TestReferenceAmoebaAngleForce.cpp.o


# Object files for target TestReferenceAmoebaAngleForce
TestReferenceAmoebaAngleForce_OBJECTS = \
"CMakeFiles/TestReferenceAmoebaAngleForce.dir/TestReferenceAmoebaAngleForce.cpp.o"

# External object files for target TestReferenceAmoebaAngleForce
TestReferenceAmoebaAngleForce_EXTERNAL_OBJECTS =

TestReferenceAmoebaAngleForce: plugins/amoeba/platforms/reference/tests/CMakeFiles/TestReferenceAmoebaAngleForce.dir/TestReferenceAmoebaAngleForce.cpp.o
TestReferenceAmoebaAngleForce: plugins/amoeba/platforms/reference/tests/CMakeFiles/TestReferenceAmoebaAngleForce.dir/build.make
TestReferenceAmoebaAngleForce: libOpenMMAmoebaReference.so
TestReferenceAmoebaAngleForce: libOpenMMAmoeba.so
TestReferenceAmoebaAngleForce: libOpenMM.so
TestReferenceAmoebaAngleForce: /usr/lib/x86_64-linux-gnu/libdl.so
TestReferenceAmoebaAngleForce: plugins/amoeba/platforms/reference/tests/CMakeFiles/TestReferenceAmoebaAngleForce.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/sam/github/openmm/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../../../../../TestReferenceAmoebaAngleForce"
	cd /home/sam/github/openmm/build/plugins/amoeba/platforms/reference/tests && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/TestReferenceAmoebaAngleForce.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
plugins/amoeba/platforms/reference/tests/CMakeFiles/TestReferenceAmoebaAngleForce.dir/build: TestReferenceAmoebaAngleForce

.PHONY : plugins/amoeba/platforms/reference/tests/CMakeFiles/TestReferenceAmoebaAngleForce.dir/build

plugins/amoeba/platforms/reference/tests/CMakeFiles/TestReferenceAmoebaAngleForce.dir/requires: plugins/amoeba/platforms/reference/tests/CMakeFiles/TestReferenceAmoebaAngleForce.dir/TestReferenceAmoebaAngleForce.cpp.o.requires

.PHONY : plugins/amoeba/platforms/reference/tests/CMakeFiles/TestReferenceAmoebaAngleForce.dir/requires

plugins/amoeba/platforms/reference/tests/CMakeFiles/TestReferenceAmoebaAngleForce.dir/clean:
	cd /home/sam/github/openmm/build/plugins/amoeba/platforms/reference/tests && $(CMAKE_COMMAND) -P CMakeFiles/TestReferenceAmoebaAngleForce.dir/cmake_clean.cmake
.PHONY : plugins/amoeba/platforms/reference/tests/CMakeFiles/TestReferenceAmoebaAngleForce.dir/clean

plugins/amoeba/platforms/reference/tests/CMakeFiles/TestReferenceAmoebaAngleForce.dir/depend:
	cd /home/sam/github/openmm/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/sam/github/openmm /home/sam/github/openmm/plugins/amoeba/platforms/reference/tests /home/sam/github/openmm/build /home/sam/github/openmm/build/plugins/amoeba/platforms/reference/tests /home/sam/github/openmm/build/plugins/amoeba/platforms/reference/tests/CMakeFiles/TestReferenceAmoebaAngleForce.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : plugins/amoeba/platforms/reference/tests/CMakeFiles/TestReferenceAmoebaAngleForce.dir/depend

