#
# Your previous .profile  (if any) is saved as .profile.dpsaved
# Setting the path for DarwinPorts.
export PATH=/Applications/TextEdit.app/Contents/MacOS:/Users/samuelflores/Desktop/NAMD_2.6_MacOSX-i686:/Applications/VMD1.8.6.app/Contents/MacOS:/usr/local/mysql/bin:/usr/local/bin:$PATH
#export PATH=/Applications/TextEdit.app/Contents/MacOS:/Users/samuelflores/Desktop/NAMD_2.6_MacOSX-i686:/Applications/VMD1.8.6.app/Contents/MacOS:/opt/local/bin:/opt/local/sbin:/usr/local/mysql/bin:/usr/local/bin:$PATH

alias rb='/Users/samuelflores/svn/RNAToolbox/src/RNABuilder.exec'
#alias rb='/Users/samuelflores/svn/RNAToolbox_v1_0/src/RNABuilder.exec'
alias mr='rm RNABuilder.o RNABuilder.exec; make RNABuilder.o; make RNABuilder.exec';
alias ct='cat  trajectory.1.pdb trajectory.2.pdb trajectory.3.pdb trajectory.4.pdb trajectory.5.pdb  trajectory.6.pdb trajectory.7.pdb  trajectory.8.pdb trajectory.9.pdb   trajectory.10.pdb trajectory.11.pdb trajectory.12.pdb trajectory.13.pdb trajectory.14.pdb trajectory.15.pdb  trajectory.16.pdb trajectory.17.pdb  trajectory.18.pdb trajectory.19.pdb   trajectory.20.pdb trajectory.21.pdb trajectory.22.pdb trajectory.23.pdb trajectory.24.pdb trajectory.25.pdb   trajectory.26.pdb trajectory.27.pdb trajectory.28.pdb trajectory.29.pdb trajectory.30.pdb trajectory.31.pdb trajectory.32.pdb trajectory.33.pdb trajectory.34.pdb trajectory.35.pdb   trajectory.36.pdb trajectory.37.pdb trajectory.38.pdb trajectory.39.pdb trajectory.40.pdb trajectory.41.pdb trajectory.42.pdb trajectory.43.pdb trajectory.44.pdb trajectory.45.pdb   trajectory.46.pdb trajectory.47.pdb trajectory.48.pdb trajectory.49.pdb trajectory.50.pdb trajectory.51.pdb trajectory.52.pdb trajectory.53.pdb trajectory.54.pdb trajectory.55.pdb   trajectory.56.pdb trajectory.57.pdb trajectory.58.pdb trajectory.59.pdb trajectory.60.pdb trajectory.61.pdb trajectory.62.pdb trajectory.63.pdb trajectory.64.pdb trajectory.65.pdb   trajectory.66.pdb trajectory.67.pdb trajectory.68.pdb trajectory.69.pdb trajectory.70.pdb trajectory.71.pdb trajectory.72.pdb trajectory.73.pdb trajectory.74.pdb trajectory.75.pdb   trajectory.76.pdb trajectory.77.pdb trajectory.78.pdb trajectory.79.pdb trajectory.80.pdb    > temp.pdb'

#export BIOX2_USER=scflores  
export RNATOOLBOX_PREFIX=/Users/samuelflores/svn/RNAToolbox
#export RNATOOLBOX_PREFIX=/Users/samuelflores/svn/RNAToolbox_v1_0
export SVN_EDITOR=/usr/bin/vi
# This is needed for SimTK
#export DYLD_LIBRARY_PATH=/Developer/SimTK/lib
export DYLD_LIBRARY_PATH=/usr/local/SimTK/lib:/usr/local/cuda/lib/:/Users/samuelflores/svn/OpenMMPreview2-Mac/lib/
#export SIMTK_INCLUDE_PATH=/Developer/SimTK 
export SIMTK_INCLUDE_PATH=/usr/local/SimTK/include
#export SIMTK_PREFIX=/Developer/SimTK
export SIMTK_PREFIX=/usr/local/SimTK
export PATH=/usr/local/cuda/bin:$PATH
