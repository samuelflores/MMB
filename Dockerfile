

#############
# New Dockerfile for MMB 4.0.0.      
# 
#############

#Download base image 

#FROM ubuntu:22.04 as OSSetup
FROM samuelflores/molmodel:v3.1.0 as MolModelImage
>>>>>>> 8e8d5967252cb34c982d83836b6b6efd2b183f34


#############
# Build MMB     
#############
# Clone specifically version 4.0.0 of MMB:
RUN git clone -b v4.0.0 https://github.com/samuelflores/MMB.git /github/MMB
RUN mkdir /github/MMB/build
RUN mkdir /github/MMB/documentation
WORKDIR   /github/MMB/documentation
RUN wget http://pe1.scilifelab.se/MMB-annex//MMB.4.tutorial.pdf
RUN wget http://pe1.scilifelab.se/MMB-annex//MMB.4.Reference-Guide.pdf        
WORKDIR  /github/MMB/build
#RUN git checkout tags/v4.0.0
# Important to specify Release build. Otherwise will not link with SimTKCommon properly.
RUN cmake -DGEMMI_DIR=/github/gemmi -DSEQAN_DIR=/github/seqan -DCMAKE_BUILD_TYPE=Release  ..
RUN make -j8
RUN make install
#############

#############
# Set up environment
#############
#inherit LD_LIBRARY_PATH directly from samuelflores/molmodel
#ENV LD_LIBRARY_PATH /usr/local/lib
ENV PATH "/usr/local/bin:$PATH"
RUN chmod a+x /github/MMB/scripts/05-MMB-splash
RUN cp /github/MMB/scripts/05-MMB-splash /etc/update-motd.d/
RUN cat "/etc/update-motd.d/00-header" >> /root/.bashrc
# Now on startup, /root/.bashrc should be read an a slightly informative splash message should appear:
RUN cat "/etc/update-motd.d/05-MMB-splash" >> /root/.bashrc
WORKDIR /work                  
RUN cp    /github/MMB/include/parameters.csv .
#############
# Done!
#############
