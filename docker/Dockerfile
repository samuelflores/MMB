#Download base image ubuntu 18.04.3
FROM ubuntu:18.04 as OSSetup

run echo "Y" | apt-get update
run echo "Y" | apt-get install subversion
run echo "Y" | apt-get install cmake
run echo "Y" | apt-get install cmake-curses-gui
run echo "Y" | apt-get install g++
run echo "Y" | apt-get install vim
run echo "Y" | apt-get install doxygen
run echo "Y" | apt-get install swig

run echo "Y" | apt-get install python
run echo "Y" | apt-get install libblas-dev liblapack-dev
run echo "Y" | apt-get install git-core


run echo "Y" | apt-get install apache2
run echo "Y" | apt-get install subversion libapache2-mod-svn  libsvn-dev


run cd /
run mkdir github
run mkdir svn
run git clone https://github.com/pandegroup/openmm.git   /github/openmm
run git clone https://github.com/simbody/simbody.git /github/simbody
run git clone https://github.com/seqan/seqan.git /github/seqan
run svn checkout https://simtk.org/svn/molmodel/trunk /svn/molmodel
#run svn checkout https://simtk.org/svn/rnatoolbox/trunk /svn/RNAToolbox/



run mkdir /github/openmm/build
WORKDIR /github/openmm/build
run cmake ..
run make install 
# default install directory is /usr/usr/local/openmm

#run mkdir /github/simbody/build ; cd /github/simbody/build ; cmake .. ; make install 
run mkdir /github/simbody/build 
#; cd /github/simbody/build ; cmake .. ; make install 
WORKDIR  /github/simbody/build
run cmake ..
#install prefix is /usr/local
run make install 


run mkdir /svn/molmodel/build
WORKDIR /svn/molmodel/build
run cmake  -DSimbody_DIR=/usr/local/lib/cmake/simbody/ -DCMAKE_INSTALL_PREFIX=/usr/local -DBUILD_EXAMPLES=OFF -DBUILD_TESTING=OFF -DBUILD_TESTING_SHARED=OFF -DBUILD_TESTING_STATIC=OFF ..
run make install

# MMB part
run git clone https://github.com/samuelflores/MMB.git /github/MMB
run mkdir /github/MMB/build
WORKDIR /github/MMB/build
run git pull --all
run cmake   -DOpenMM_INSTALL_DIR=//usr/local/openmm  -DOpenMM_INCLUDE_DIR="/usr/local/openmm/include/openmm/reference/;/usr/local/openmm/include/openmm/;/usr/local/openmm/include"  -DCMAKE_BUILD_TYPE=Release -DSeqAn_INCLUDE_DIR=/github/seqan/include -DCMAKE_CXX_FLAGS="-std=c++14 -D BuildNtC -D USE_OPENMM" -DSimTK_INSTALL_DIR=/usr/local -DSimbody_DIR=/usr/local/lib/cmake/simbody/ -DCMAKE_PREFIX_PATH=/usr/local -DCMAKE_INSTALL_PREFIX=/usr/local ..
run touch /github/MMB/build/done-cmake.txt
run make install
run rm MMB libMMBlib.so
ENV LD_LIBRARY_PATH /usr/local/lib
ENV PATH "/usr/local/bin:$PATH"
WORKDIR /work                  
run cp    /github/MMB/include/resources/parameters.csv .

# Trying multi-stage build. first, get a clean ubuntu image:
#FROM ubuntu:18.04 as final
#run mkdir /github
#run mkdir /svn
#COPY --from=OSSetup /github/* /github/
#COPY --from=OSSetup /svn/* /svn/
#COPY --from=OSSetup /usr/local/* /usr/local/
#run echo "Y" | apt-get update
#run echo "Y" | apt-get install libblas-dev liblapack-dev
#run echo "Y" | apt-get install git-core
#ENV LD_LIBRARY_PATH /usr/local/lib
#ENV PATH "/usr/local/bin:$PATH"
