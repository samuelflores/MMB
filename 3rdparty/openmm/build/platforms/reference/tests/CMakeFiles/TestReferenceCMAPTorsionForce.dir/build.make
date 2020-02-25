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
include platforms/reference/tests/CMakeFiles/TestReferenceCMAPTorsionForce.dir/depend.make

# Include the progress variables for this target.
include platforms/reference/tests/CMakeFiles/TestReferenceCMAPTorsionForce.dir/progress.make

# Include the compile flags for this target's objects.
include platforms/reference/tests/CMakeFiles/TestReferenceCMAPTorsionForce.dir/flags.make

platforms/reference/tests/CMakeFiles/TestReferenceCMAPTorsionForce.dir/TestReferenceCMAPTorsionForce.cpp.o: platforms/reference/tests/CMakeFiles/TestReferenceCMAPTorsionForce.dir/flags.make
platforms/reference/tests/CMakeFiles/TestReferenceCMAPTorsionForce.dir/TestReferenceCMAPTorsionForce.cpp.o: ../platforms/reference/tests/TestReferenceCMAPTorsionForce.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/sam/github/openmm/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object platforms/reference/tests/CMakeFiles/TestReferenceCMAPTorsionForce.dir/TestReferenceCMAPTorsionForce.cpp.o"
	cd /home/sam/github/openmm/build/platforms/reference/tests && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/TestReferenceCMAPTorsionForce.dir/TestReferenceCMAPTorsionForce.cpp.o -c /home/sam/github/openmm/platforms/reference/tests/TestReferenceCMAPTorsionForce.cpp

platforms/reference/tests/CMakeFiles/TestReferenceCMAPTorsionForce.dir/TestReferenceCMAPTorsionForce.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/TestReferenceCMAPTorsionForce.dir/TestReferenceCMAPTorsionForce.cpp.i"
	cd /home/sam/github/openmm/build/platforms/reference/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/sam/github/openmm/platforms/reference/tests/TestReferenceCMAPTorsionForce.cpp > CMakeFiles/TestReferenceCMAPTorsionForce.dir/TestReferenceCMAPTorsionForce.cpp.i

platforms/reference/tests/CMakeFiles/TestReferenceCMAPTorsionForce.dir/TestReferenceCMAPTorsionForce.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/TestReferenceCMAPTorsionForce.dir/TestReferenceCMAPTorsionForce.cpp.s"
	cd /home/sam/github/openmm/build/platforms/reference/tests && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/sam/github/openmm/platforms/reference/tests/TestReferenceCMAPTorsionForce.cpp -o CMakeFiles/TestReferenceCMAPTorsionForce.dir/TestReferenceCMAPTorsionForce.cpp.s

platforms/reference/tests/CMakeFiles/TestReferenceCMAPTorsionForce.dir/TestReferenceCMAPTorsionForce.cpp.o.requires:

.PHONY : platforms/reference/tests/CMakeFiles/TestReferenceCMAPTorsionForce.dir/TestReferenceCMAPTorsionForce.cpp.o.requires

platforms/reference/tests/CMakeFiles/TestReferenceCMAPTorsionForce.dir/TestReferenceCMAPTorsionForce.cpp.o.provides: platforms/reference/tests/CMakeFiles/TestReferenceCMAPTorsionForce.dir/TestReferenceCMAPTorsionForce.cpp.o.requires
	$(MAKE) -f platforms/reference/tests/CMakeFiles/TestReferenceCMAPTorsionForce.dir/build.make platforms/reference/tests/CMakeFiles/TestReferenceCMAPTorsionForce.dir/TestReferenceCMAPTorsionForce.cpp.o.provides.build
.PHONY : platforms/reference/tests/CMakeFiles/TestReferenceCMAPTorsionForce.dir/TestReferenceCMAPTorsionForce.cpp.o.provides

platforms/reference/tests/CMakeFiles/TestReferenceCMAPTorsionForce.dir/TestReferenceCMAPTorsionForce.cpp.o.provides.build: platforms/reference/tests/CMakeFiles/TestReferenceCMAPTorsionForce.dir/TestReferenceCMAPTorsionForce.cpp.o


# Object files for target TestReferenceCMAPTorsionForce
TestReferenceCMAPTorsionForce_OBJECTS = \
"CMakeFiles/TestReferenceCMAPTorsionForce.dir/TestReferenceCMAPTorsionForce.cpp.o"

# External object files for target TestReferenceCMAPTorsionForce
TestReferenceCMAPTorsionForce_EXTERNAL_OBJECTS =

TestReferenceCMAPTorsionForce: platforms/reference/tests/CMakeFiles/TestReferenceCMAPTorsionForce.dir/TestReferenceCMAPTorsionForce.cpp.o
TestReferenceCMAPTorsionForce: platforms/reference/tests/CMakeFiles/TestReferenceCMAPTorsionForce.dir/build.make
TestReferenceCMAPTorsionForce: libOpenMM.so
TestReferenceCMAPTorsionForce: /usr/lib/x86_64-linux-gnu/libdl.so
TestReferenceCMAPTorsionForce: platforms/reference/tests/CMakeFiles/TestReferenceCMAPTorsionForce.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/sam/github/openmm/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../../../TestReferenceCMAPTorsionForce"
	cd /home/sam/github/openmm/build/platforms/reference/tests && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/TestReferenceCMAPTorsionForce.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
platforms/reference/tests/CMakeFiles/TestReferenceCMAPTorsionForce.dir/build: TestReferenceCMAPTorsionForce

.PHONY : platforms/reference/tests/CMakeFiles/TestReferenceCMAPTorsionForce.dir/build

platforms/reference/tests/CMakeFiles/TestReferenceCMAPTorsionForce.dir/requires: platforms/reference/tests/CMakeFiles/TestReferenceCMAPTorsionForce.dir/TestReferenceCMAPTorsionForce.cpp.o.requires

.PHONY : platforms/reference/tests/CMakeFiles/TestReferenceCMAPTorsionForce.dir/requires

platforms/reference/tests/CMakeFiles/TestReferenceCMAPTorsionForce.dir/clean:
	cd /home/sam/github/openmm/build/platforms/reference/tests && $(CMAKE_COMMAND) -P CMakeFiles/TestReferenceCMAPTorsionForce.dir/cmake_clean.cmake
.PHONY : platforms/reference/tests/CMakeFiles/TestReferenceCMAPTorsionForce.dir/clean

platforms/reference/tests/CMakeFiles/TestReferenceCMAPTorsionForce.dir/depend:
	cd /home/sam/github/openmm/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/sam/github/openmm /home/sam/github/openmm/platforms/reference/tests /home/sam/github/openmm/build /home/sam/github/openmm/build/platforms/reference/tests /home/sam/github/openmm/build/platforms/reference/tests/CMakeFiles/TestReferenceCMAPTorsionForce.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : platforms/reference/tests/CMakeFiles/TestReferenceCMAPTorsionForce.dir/depend

