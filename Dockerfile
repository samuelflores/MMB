#Download base image ubuntu 18.04.3

FROM ubuntu:18.04 as OSSetup

#the way to override this is 
# docker build --build-arg GIT_COMMIT=`git rev-parse  HEAD`
#ARG GIT_COMMIT=83c4aeadeed0a9644b24f13808c068a4d51531a0
# this is an old commit
ARG GIT_COMMIT=
ARG GIT_COMMIT=b977ebfdd1928ce4f35d426215f68c09569bc1ce
LABEL git_commit=$GIT_COMMIT




RUN apt-get update && apt-get install -y  apt-utils     cmake    cmake-curses-gui    g++    vim    doxygen    swig python  libblas-dev liblapack-dev  git-core subversion mc wget
   
#RUN apt-get install -y mc

WORKDIR /
#RUN cd /
RUN mkdir /github /svn
# Removed openmm from here. This is because openmm is now forked into MMB, with an eye to debian packaging.
#RUN mkdir /svn
#//RUN git clone https://github.com/pandegroup/openmm.git   /github/openmm
RUN git clone https://github.com/simbody/simbody.git /github/simbody
RUN git clone https://github.com/seqan/seqan.git /github/seqan
#RUN svn checkout https://simtk.org/svn/molmodel/trunk /svn/molmodel
RUN git clone  https://github.com/samuelflores/molmodel.git /github/molmodel

RUN mkdir /Documentation
WORKDIR /Documentation
# Get the MMB reference guide from its repository on pe1:
RUN wget http://pe1.scilifelab.se/MMB-annex/Documentation/MMB.3_0.Reference-Guide.pdf
RUN wget http://pe1.scilifelab.se/MMB-annex/Documentation/MMB.3_2.tutorial.pdf       
RUN wget http://pe1.scilifelab.se/MMB-annex/Documentation/MMB.3_2.tutorial.docx        

#run mkdir /github/simbody/build ; cd /github/simbody/build ; cmake .. ; make install 
RUN mkdir /github/simbody/build 
#; cd /github/simbody/build ; cmake .. ; make install 
WORKDIR  /github/simbody/build
RUN cmake ..
#install prefix is /usr/local
RUN make install 

# MMB part
RUN git clone https://github.com/samuelflores/MMB.git /github/MMB

# Now openmm will install from MMB repo:
RUN mkdir /github/MMB/3rdparty/openmm/build
WORKDIR /github/MMB/3rdparty/openmm/build
RUN cmake ..
RUN make install 
# default install directory is /usr/local/openmm

#WORKDIR /github/MMB
# in the case of MMB 3.0:

##################
### CCP4 Maps Part
##################

# Install dependencies
RUN apt-get update && apt-get install -y python-dev gfortran m4

# Copy and install MMDB2
WORKDIR /github/MMB/3rdparty/mmdb2
RUN ./configure --enable-shared && make && make install

# Copy and install libccp4
WORKDIR /github/MMB/3rdparty/libccp4
RUN ./configure --enable-shared && make && make install

##################
### Molmodel
##################
#RUN mkdir /svn/molmodel/build
#WORKDIR /svn/molmodel/build
RUN mkdir /github/molmodel/build
WORKDIR /github/molmodel/build
RUN cmake -DADD_MMDB2_LIBRARY=TRUE -DSimbody_DIR=/usr/local/lib/cmake/simbody/ -DCMAKE_INSTALL_PREFIX=/usr/local -DBUILD_EXAMPLES=OFF -DBUILD_TESTING=OFF -DBUILD_TESTING_SHARED=OFF -DBUILD_TESTING_STATIC=OFF ..
RUN make install

#RUN git -C /github/MMB  checkout $GIT_COMMIT        

RUN mkdir /github/MMB/build
WORKDIR /github/MMB/build
RUN git pull --all
RUN cmake -DBuild_CCP4=TRUE -DLIBCCP4_INCLUDE_DIR=/github/MMB/3rdparty/include -DLIBCCP4_LIB_DIR=/github/MMB/3rdparty/lib  -DOpenMM_INSTALL_DIR=//usr/local/openmm  -DOpenMM_INCLUDE_DIR="/usr/local/openmm/include/openmm/reference/;/usr/local/openmm/include/openmm/;/usr/local/openmm/include"  -DCMAKE_BUILD_TYPE=Release -DSeqAn_INCLUDE_DIR=/github/seqan/include -DCMAKE_CXX_FLAGS="-std=c++14 -D BuildNtC -D USE_OPENMM -D  MMDB2_LIB_USAGE " -DSimTK_INSTALL_DIR=/usr/local -DSimbody_DIR=/usr/local/lib/cmake/simbody/ -DCMAKE_PREFIX_PATH=/usr/local -DCMAKE_INSTALL_PREFIX=/usr/local  ..
RUN touch /github/MMB/build/done-cmake.txt
RUN make install
RUN rm MMB libMMBlib.so
ENV LD_LIBRARY_PATH /usr/local/lib
ENV PATH "/usr/local/bin:$PATH"
WORKDIR /work                  
RUN cp    /github/MMB/include/resources/parameters.csv .
#RUN mv /github/MMB/docker/Dockerfile /
# Trying multi-stage build. first, get a clean ubuntu image:
#FROM ubuntu:18.04 as final
#run mkdir /github
#run mkdir /svn
#COPY --from=OSSetup /github/* /github/
#COPY --from=OSSetup /svn/* /svn/
#COPY --from=OSSetup /usr/local/* /usr/local/
#run echo "Y" | apt-get update
#run apt-get install -y  libblas-dev liblapack-dev
#run apt-get install -y  git-core
#ENV LD_LIBRARY_PATH /usr/local/lib
#ENV PATH "/usr/local/bin:$PATH"
