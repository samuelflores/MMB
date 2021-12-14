# MMB

## Build instructions

MMB uses [CMake](https://cmake.org/) build system to configure and build itself. MMB provides various build-time configuration options. Depending on how your particular system is set up, some of these options may have to be adjusted in order to build MMB.

MMB depends on a series of libraries that must be installed in order to build and run MMB. Current list of MMB dependencies is as follows:
- [Simbody](https://simtk.org/projects/simbody)
- [OpenMM](https://openmm.org/)
- [Gemmi](https://gemmi.readthedocs.io/en/latest/)
- [Molmodel](https://simtk.org/projects/molmodel)
- [SeqAn](https://github.com/seqan/seqan)

### List of configuration options
- `BUILD_SHARED` [`ON|OFF`] - Whether to build MMB library as dynamically linked (shared) library. Defaults to `ON`
- `BUILD_STATIC` [`ON|OFF`] - Whether to build MMB library as statically linked library. Defaults to `OFF`

Note that at least one of the options above shall be enabled

- `BUILD_PYTHON_WRAPPER` [`ON|OFF`] - Whether to build Python wrapper around MMB library API. Defaults to `OFF`
- `ASSUME_VERSIONED_SIMBODY` [`ON|OFF`] - Simbody libraries may be built with the Simbody version tag appended to the file names of the libraries. If your Simbody installation is built like this, enable this option to allow MMB to link properly against Simbody. Defaults to `OFF`
- `ENABLE_LEPTON` [`ON|OFF`] - Enable Lepton math parser. Defaults to `ON`
- `ENABLE_NTC` [`ON|OFF`] - Enable nucleotide conformers calculation code. Defaults to `ON`

- `GEMMI_DIR` - Path to [Gemmi](https://gemmi.readthedocs.io/en/latest/) library installation. This shall be specified only if Gemmi is not installed as a system-wide library.
- `SEQAN_DIR` - Path to [SeqAn](https://github.com/seqan/seqan) library installation. This shall be specified only if SeqAn is not installed as a system-wide library.
- `LEPTON_DIR` - Path to [Lepton](https://simtk.org/projects/lepton) math parser source code. MMB provides its own copy of Lepton source code which will be used if this variable is not set.

### Build and installation

Under most circumstances, MMB should be able to set itself up for build automatically. It may be necessary to pass some additional configuration parameters to CMake if some of MMB dependencies are not installed in standard system paths.

#### UNIX systems

To build Simbody on a UNIX system, `cd` into the directory with the MMB source code and issue the following commands

    mkdir build
    cd build
    cmake .. -DCMAKE_BUILD_TYPE=Release
    make
    make install

note that the `make install` command may require root privileges.

If either Simbody, OpenMM or Molmodel are not installed in standard system paths, it is necessary to tell CMake where to look for them. This can be done with the `CMAKE_PREFIX_PATH` variable. For example, if your Simbody is installed in `/opt/simbody`, use the following command to point CMake to the correct directory.

    cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=/opt/simbody

If you need to specify more than one path, use the `;` (semicolon) to input multiple paths, i.e.

    -DCMAKE_PREFIX_PATH=/path/to/simbody;/path/to/molmodel

#### Windows systems

Configuration and build process on Windows is similar to that on UNIX systems. Since there are no standard system paths on Windows, it is necessary to specify paths to all MMB dependencies. Note that CMake for Windows comes with a graphical user interface to input the configuration options.
