/* -------------------------------------------------------------------------- *
 *                           MMB (MacroMoleculeBuilder)                       *
 * -------------------------------------------------------------------------- *
 *                                                                            *
 * Copyright (c) 2011-12 by the Author.                                       *
 * Author: Samuel Flores                                                      *
 *                                                                            *
 * See RNABuilder.cpp for the copyright and usage agreement.                  *
 * -------------------------------------------------------------------------- */

#ifndef Spiral_H_
#define Spiral_H_
#include <regex>
//#include "ResidueStretchContainer.h"
//#include "ParameterReader.h"
#include "SimTKmolmodel.h"
#include "MMBLogger.h"    
#include "Utils.h"    
#include <vector>
#include <map>
#include <fstream>
using namespace SimTK;
using namespace std  ;
struct frequencyPhaseAmplitude {
    double frequency; // By convention, use rads. Please use fullest precision possible, because angle could be accumulative from beginning of DNA strand to the end.
    double phase    ; // also rads
    double amplitude; // also rads. Because this is angle in the theta-direction, that is to say, along the meridians.
    // We can also use the hard coded SimTK::Pi
};
class Spiral                                                            {
private:
    vector<frequencyPhaseAmplitude> frequencyPhaseAmplitudeVector; // This holds all the harmonic adjustments to be made to the spherical helix.
    double helixAdvancePerBasePair  ;  // in nm
    String spiralPdbFileName  ;
    String spiralCommandsFileName  ;
    Vec3   center	      ;
    double radius             ;
    double interStrandDistance; 
    double startTheta         ;
    double endTheta         ;
    double phiOffset          ;
    double harmonicThetaAdjustment(double inputPhi);
public:
    Spiral(); // Default constructor. Should set default values of parameters.
    void clear();
    void validate();
    void parseInput(String commandName);
    void parseInput(const map<const String,double> & userVariables, String parameterName, String parameterValue);
    void parseInput(const map<const String,double> & userVariables, String parameterName, String parameterValue1, String parameterValue2,String parameterValue3 );
    void parseInput(String parameterName, double X, double Y, double Z);
    void writeSyntax(); // To be written as a helpful aid to the user.
    void writeIonSpiralPdbFile();
    void writeDnaSpiralCommandfile();
};
#endif

