generate-base-pair-transform.o : ../src/generate-base-pair-transform.cpp
	g++ -o generate-base-pair-transform.o -c ../src/generate-base-pair-transform.cpp  -I/usr/local/SimTK/include/molmodel/internal  -I/usr/local/SimTK/include -I/usr/include/c++/4.2.1  -arch i386 


generate-base-pair-transform : generate-base-pair-transform.o
	g++ -o generate-base-pair-transform  generate-base-pair-transform.o   -I/usr/local/SimTK/include/SimTKcommon/internal  -I/usr/local/SimTK/include/molmodel/internal   -I/usr/local/SimTK/include -L/usr/local/SimTK/lib  -lSimTKcommon -lSimTKmolmodel -lSimTKsimbody -lSimTKmath -llapack -lblas -arch i386
#-02 -liostream -lfstream -lstdlib -lstring -I

