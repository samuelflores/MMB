

#############
# New Dockerfile for MMB 4.0.0.      
# 
#############

#Download base image 
FROM ubuntu:22.04 as OSSetup

#############
#Get some libraries and executables
#############
# cmake requires build-essential libssl-dev
RUN apt update
#RUN apt upgrade
RUN apt install -y git wget swig doxygen libblas-dev liblapack-dev  cmake-curses-gui zlib1g zlib1g-dev apt-utils snapd build-essential libssl-dev software-properties-common lsb-release vim
RUN mkdir /github
#############

#############
#for cmake:
#############
WORKDIR /github
RUN wget https://github.com/Kitware/CMake/archive/refs/tags/v3.25.0-rc3.tar.gz
RUN tar -zxvf v3.25.0-rc3.tar.gz
WORKDIR /github/CMake-3.25.0-rc3
# Build cmake
RUN ./bootstrap
RUN gmake -j8
RUN gmake install
#############

#############
# clone git repos
#############
WORKDIR /github
RUN git clone https://github.com/project-gemmi/gemmi.git /github/gemmi
#get a specific release of openmm:
RUN git clone https://github.com/pandegroup/openmm.git   /github/openmm --branch 7.7.0
RUN git clone  -b simbody-3.7 --single-branch https://github.com/simbody/simbody.git /github/simbody
RUN git clone https://github.com/seqan/seqan.git /github/seqan
RUN git clone  https://github.com/samuelflores/molmodel.git /github/molmodel
RUN git clone https://github.com/samuelflores/MMB.git /github/MMB
# make build directories
RUN mkdir /github/gemmi/build
RUN mkdir /github/simbody/build
RUN mkdir /github/molmodel/build
RUN mkdir /github/openmm/build
RUN mkdir /github/MMB/build
#############

#############
# build gemmi    
#############
WORKDIR /github/gemmi/build
#latest build is not compatible with MMB. checkout version 5.6, 21 sept 2022:
RUN git checkout cca2ac864f47c9690aec73c7702bb5247d665f5a
WORKDIR /github/gemmi/build
RUN cmake ..
RUN make -j8
RUN make install
#############

#############
# Build simbody 
#############
WORKDIR /github/simbody/build
RUN cmake -DUSE_GEMMI=TRUE -DGEMMI_PATH=/github/gemmi/include -DSimbody_DIR=/usr/local/lib/cmake/simbody/ -DCMAKE_INSTALL_PREFIX=/usr/local -DBUILD_EXAMPLES=OFF -DBUILD_TESTING=OFF -DBUILD_TESTING_SHARED=OFF -DBUILD_TESTING_STATIC=OFF ..
RUN make -j8
RUN make install
#############

#############
# build openmm   
#############
WORKDIR /github/openmm/build
RUN cmake ..
RUN make -j8
RUN make install
#############

#############
# Build molmodel
#############
WORKDIR  /github/molmodel/build
RUN cmake -DUSE_GEMMI=TRUE -DGEMMI_PATH=/github/gemmi/include -DSimbody_DIR=/usr/local/lib/cmake/simbody/ -DCMAKE_INSTALL_PREFIX=/usr/local -DBUILD_EXAMPLES=OFF -DBUILD_TESTING=OFF -DBUILD_TESTING_SHARED=OFF -DBUILD_TESTING_STATIC=OFF.  -DCMAKE_BUILD_TYPE=Release ..
RUN make -j8
RUN make install
#############

#############
# Build MMB     
#############
RUN mkdir /github/MMB/documentation
WORKDIR   /github/MMB/documentation
RUN wget http://pe1.scilifelab.se/MMB-annex//MMB.4.tutorial.pdf
RUN wget http://pe1.scilifelab.se/MMB-annex//MMB.4.Reference-Guide.pdf        
WORKDIR  /github/MMB/build
# checkout specifically version 4.0.0 of MMB:
RUN git checkout tags/v4.0.0
# Important to specify Release build. Otherwise will not link with SimTKCommon properly.
RUN cmake -DGEMMI_DIR=/github/gemmi -DSEQAN_DIR=/github/seqan -DCMAKE_BUILD_TYPE=Release  ..
RUN make -j8
RUN make install
#############

#############
# Set up environment
#############
ENV LD_LIBRARY_PATH /usr/local/lib
ENV PATH "/usr/local/bin:$PATH"
WORKDIR /work                  
RUN cp    /github/MMB/include/resources/parameters.csv .
#############
# Done!
#############
