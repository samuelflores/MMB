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
include plugins/amoeba/platforms/reference/tests/CMakeFiles/TestReferenceAmoebaPiTorsionForce.dir/depend.make

# Include the progress variables for this target.
include plugins/amoeba/platforms/reference/tests/CMakeFiles/TestReferenceAmoebaPiTorsionForce.dir/progress.make

# Include the compile flags for this target's objects.
include plugins/amoeba/platforms/reference/tests/CMakeFiles/TestReferenceAmoebaPiTorsionForce.dir/flags.make

plugins/amoeba/platforms/reference/tests/CMakeFiles/TestReferenceAmoebaPiTorsionForce.dir/TestReferenceAmoebaPiTorsionForce.cpp.o: plugins/amoeba/platforms/reference/tests/CMakeFiles/TestReferenceAmoebaPiTorsionForce.dir/flags.make
plugins/amoeba/platforms/reference/tests/CMakeFiles/TestReferenceAmoebaPiTorsionForce.dir/TestReferenceAmoebaPiTorsionForce.cpp.o: ../plugins/amoeba/platforms/reference/tests/TestReferenceAmoebaPiTorsionForce.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/sam/github/openmm/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object plugins/amoeba/platforms/reference/tests/CMakeFiles/TestReferenceAmoebaPiTorsionForce.dir/TestReferenceAmoebaPiTorsionForce.cpp.o"
	cd /home/sam/github/openmm/build/plugins/amoeba/platforms/reference/tests && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/TestReferenceAmoebaPiTorsionForce.dir/TestReferenceAmoebaPiTorsionForce.cpp.o -c /home/sam/github/openmm/plugins/amoeba/platforms/reference/tests/TestReferenceAmoebaPiTorsionForce.cpp

plugins/amoeba/platforms/reference/tests/CMakeFiles/TestReferenceAmoebaPiTorsionForce.dir/TestReferenceAmoebaPiTorsionForce.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/TestReferenceAmoebaPiTorsionForce.dir/TestReferenceAmoebaPiTorsionForce.cpp.i"
	cd /home/sam/github/openmm/build/plugins/amoeba/platforms/reference/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/sam/github/openmm/plugins/amoeba/platforms/reference/tests/TestReferenceAmoebaPiTorsionForce.cpp > CMakeFiles/TestReferenceAmoebaPiTorsionForce.dir/TestReferenceAmoebaPiTorsionForce.cpp.i

plugins/amoeba/platforms/reference/tests/CMakeFiles/TestReferenceAmoebaPiTorsionForce.dir/TestReferenceAmoebaPiTorsionForce.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/TestReferenceAmoebaPiTorsionForce.dir/TestReferenceAmoebaPiTorsionForce.cpp.s"
	cd /home/sam/github/openmm/build/plugins/amoeba/platforms/reference/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/sam/github/openmm/plugins/amoeba/platforms/reference/tests/TestReferenceAmoebaPiTorsionForce.cpp -o CMakeFiles/TestReferenceAmoebaPiTorsionForce.dir/TestReferenceAmoebaPiTorsionForce.cpp.s

plugins/amoeba/platforms/reference/tests/CMakeFiles/TestReferenceAmoebaPiTorsionForce.dir/TestReferenceAmoebaPiTorsionForce.cpp.o.requires:

.PHONY : plugins/amoeba/platforms/reference/tests/CMakeFiles/TestReferenceAmoebaPiTorsionForce.dir/TestReferenceAmoebaPiTorsionForce.cpp.o.requires

plugins/amoeba/platforms/reference/tests/CMakeFiles/TestReferenceAmoebaPiTorsionForce.dir/TestReferenceAmoebaPiTorsionForce.cpp.o.provides: plugins/amoeba/platforms/reference/tests/CMakeFiles/TestReferenceAmoebaPiTorsionForce.dir/TestReferenceAmoebaPiTorsionForce.cpp.o.requires
	$(MAKE) -f plugins/amoeba/platforms/reference/tests/CMakeFiles/TestReferenceAmoebaPiTorsionForce.dir/build.make plugins/amoeba/platforms/reference/tests/CMakeFiles/TestReferenceAmoebaPiTorsionForce.dir/TestReferenceAmoebaPiTorsionForce.cpp.o.provides.build
.PHONY : plugins/amoeba/platforms/reference/tests/CMakeFiles/TestReferenceAmoebaPiTorsionForce.dir/TestReferenceAmoebaPiTorsionForce.cpp.o.provides

plugins/amoeba/platforms/reference/tests/CMakeFiles/TestReferenceAmoebaPiTorsionForce.dir/TestReferenceAmoebaPiTorsionForce.cpp.o.provides.build: plugins/amoeba/platforms/reference/tests/CMakeFiles/TestReferenceAmoebaPiTorsionForce.dir/TestReferenceAmoebaPiTorsionForce.cpp.o


# Object files for target TestReferenceAmoebaPiTorsionForce
TestReferenceAmoebaPiTorsionForce_OBJECTS = \
"CMakeFiles/TestReferenceAmoebaPiTorsionForce.dir/TestReferenceAmoebaPiTorsionForce.cpp.o"

# External object files for target TestReferenceAmoebaPiTorsionForce
TestReferenceAmoebaPiTorsionForce_EXTERNAL_OBJECTS =

TestReferenceAmoebaPiTorsionForce: plugins/amoeba/platforms/reference/tests/CMakeFiles/TestReferenceAmoebaPiTorsionForce.dir/TestReferenceAmoebaPiTorsionForce.cpp.o
TestReferenceAmoebaPiTorsionForce: plugins/amoeba/platforms/reference/tests/CMakeFiles/TestReferenceAmoebaPiTorsionForce.dir/build.make
TestReferenceAmoebaPiTorsionForce: libOpenMMAmoebaReference.so
TestReferenceAmoebaPiTorsionForce: libOpenMMAmoeba.so
TestReferenceAmoebaPiTorsionForce: libOpenMM.so
TestReferenceAmoebaPiTorsionForce: /usr/lib/x86_64-linux-gnu/libdl.so
TestReferenceAmoebaPiTorsionForce: plugins/amoeba/platforms/reference/tests/CMakeFiles/TestReferenceAmoebaPiTorsionForce.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/sam/github/openmm/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../../../../../TestReferenceAmoebaPiTorsionForce"
	cd /home/sam/github/openmm/build/plugins/amoeba/platforms/reference/tests && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/TestReferenceAmoebaPiTorsionForce.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
plugins/amoeba/platforms/reference/tests/CMakeFiles/TestReferenceAmoebaPiTorsionForce.dir/build: TestReferenceAmoebaPiTorsionForce

.PHONY : plugins/amoeba/platforms/reference/tests/CMakeFiles/TestReferenceAmoebaPiTorsionForce.dir/build

plugins/amoeba/platforms/reference/tests/CMakeFiles/TestReferenceAmoebaPiTorsionForce.dir/requires: plugins/amoeba/platforms/reference/tests/CMakeFiles/TestReferenceAmoebaPiTorsionForce.dir/TestReferenceAmoebaPiTorsionForce.cpp.o.requires

.PHONY : plugins/amoeba/platforms/reference/tests/CMakeFiles/TestReferenceAmoebaPiTorsionForce.dir/requires

plugins/amoeba/platforms/reference/tests/CMakeFiles/TestReferenceAmoebaPiTorsionForce.dir/clean:
	cd /home/sam/github/openmm/build/plugins/amoeba/platforms/reference/tests && $(CMAKE_COMMAND) -P CMakeFiles/TestReferenceAmoebaPiTorsionForce.dir/cmake_clean.cmake
.PHONY : plugins/amoeba/platforms/reference/tests/CMakeFiles/TestReferenceAmoebaPiTorsionForce.dir/clean

plugins/amoeba/platforms/reference/tests/CMakeFiles/TestReferenceAmoebaPiTorsionForce.dir/depend:
	cd /home/sam/github/openmm/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/sam/github/openmm /home/sam/github/openmm/plugins/amoeba/platforms/reference/tests /home/sam/github/openmm/build /home/sam/github/openmm/build/plugins/amoeba/platforms/reference/tests /home/sam/github/openmm/build/plugins/amoeba/platforms/reference/tests/CMakeFiles/TestReferenceAmoebaPiTorsionForce.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : plugins/amoeba/platforms/reference/tests/CMakeFiles/TestReferenceAmoebaPiTorsionForce.dir/depend

