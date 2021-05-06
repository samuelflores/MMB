/* -------------------------------------------------------------------------- *
 *                           MMB (MacroMoleculeBuilder)                       *
 * -------------------------------------------------------------------------- *
 *                                                                            *
 * Copyright (c) 2011-12 by the Author.                                       *
 * Author: Samuel Flores                                                      *
 *                                                                            *
 * See RNABuilder.cpp for the copyright and usage agreement.                  *
 * -------------------------------------------------------------------------- */

#include "Spiral.h"
#include "MonoAtoms.h"

using namespace std;


std::string commonSpiralCommands = R"(
    numReportingIntervals 1
    reportingInterval .00001
    firstStage 2
    lastStage ZZZZ )";


double   phiFromXYZ( const Vec3 myPoint, const  Vec3 sphericalCenter) {
    double xDist = myPoint[0]-sphericalCenter[0];
    double yDist = myPoint[1]-sphericalCenter[1];
    double sliceRadius = sqrt(xDist*xDist+yDist*yDist);
    double myPhi = acos(xDist / sliceRadius);
    MMBLOG_FILE_FUNC_LINE(INFO, "myPhi = " << myPhi<<endl);
    return myPhi;
}
double thetaFromXYZ( const Vec3 myPoint, const  Vec3 sphericalCenter) {
    double xDist = myPoint[0]-sphericalCenter[0];
    double yDist = myPoint[1]-sphericalCenter[1];
    double zDist = myPoint[2]-sphericalCenter[2];
    double sliceRadius = sqrt(xDist*xDist+yDist*yDist);
    double computedSphericalRadius = sqrt(xDist*xDist+yDist*yDist+zDist*zDist);
    MMBLOG_FILE_FUNC_LINE(INFO, "slice radius = " << sliceRadius<<endl);
    double theta = asin(sliceRadius / computedSphericalRadius);
    MMBLOG_FILE_FUNC_LINE(INFO, "computedSphericalRadius = " << computedSphericalRadius<<endl);
    MMBLOG_FILE_FUNC_LINE(INFO, "theta = " << theta<<endl);
    return asin(sqrt(xDist*xDist+yDist*yDist) / sqrt(xDist*xDist+yDist*yDist+zDist*zDist));
}
double zSliceRadius ( const double mySphereRadius, const double myTheta){
    return mySphereRadius*sin(myTheta);
}
// computes arc length in a cylindrical helix, for a certain deltaTheta
// takes deltaPhi in radians
double helicalSpiralArcLength ( const double myCylindricalRadius, const double helicalPitch, const double deltaPhi){
    return deltaPhi * sqrt( myCylindricalRadius * myCylindricalRadius + helicalPitch * helicalPitch);
}
double deltaPhiFromCylindricalRadiusHelicalPitchAndHelicalArcLength( const double myCylindricalRadius, const double myHelicalPitch, const double myArcLength){
    return myArcLength / sqrt( myCylindricalRadius * myCylindricalRadius + myHelicalPitch * myHelicalPitch);
}
double deltaPhiFromThetaInterHelicalDistanceSphericalRadiusAndHelicalArcLength( const double myTheta, const double myInterHelicalDistance, const  double mySphericalRadius, const double myArcLength){
    double myHelicalPitch = myInterHelicalDistance*sin(myTheta);
    double myCylindricalRadius = mySphericalRadius * sin(myTheta);
    return  deltaPhiFromCylindricalRadiusHelicalPitchAndHelicalArcLength(myCylindricalRadius, myHelicalPitch, myArcLength);
}
double thetaFromPhi( const double phi, const double myInterHelicalDistance, const double mySphericalRadius, const double myPhiOffset){
    MMBLOG_FILE_FUNC_LINE(DEBUG, " locally, phiOffset = "<<myPhiOffset<<endl);
    return (phi - myPhiOffset) / (2*SimTK::Pi/myInterHelicalDistance*mySphericalRadius) ;
}
double phiFromTheta(const double theta,const double myInterHelicalDistance, const double mySphericalRadius, const double myPhiOffset){
    MMBLOG_FILE_FUNC_LINE(DEBUG, " locally, phiOffset = "<<myPhiOffset<<endl);
    return theta * (2*SimTK::Pi/myInterHelicalDistance*mySphericalRadius) + myPhiOffset ;
}


Spiral::Spiral(){clear();}
void   Spiral::clear(){
    frequencyPhaseAmplitudeVector.clear();    
    spiralPdbFileName = String("NOT-SET"); // was "spiral.pdb"
    MMBLOG_FILE_FUNC_LINE(INFO, " Setting spiralPdbFileName = " <<spiralPdbFileName                                   <<endl);
    spiralCommandsFileName = String("NOT-SET"); // was "commands.spiral.dat"
    center            = Vec3(-11111.,-11111.,-11111.);
    radius            = -11111.;
    interStrandDistance=-11111.;
    startTheta        = -11111.;
    endTheta          = -11111.;
    phiOffset         = -11111.;
    double helixAdvancePerBasePair = -11111. ;  // in nm // should be user specified.

}
void   Spiral::validate(){
    //frequencyPhaseAmplitudeVector.clear();    
    if (spiralPdbFileName == String("NOT-SET"))
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "Invalid value for spiralPdbFileName : "<<spiralPdbFileName<< endl);
    if (spiralCommandsFileName == String("NOT-SET"))
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "Invalid value for spiralCommandsFileName : "<<spiralCommandsFileName<< endl);
    if (radius         == -11111.)
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "Invalid value for radius : "<<radius<< endl);
    if (interStrandDistance         == -11111.)
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "Invalid value for interStrandDistance : "<<interStrandDistance<< endl);
    if (startTheta         == -11111.)
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "Invalid value for startTheta : "<<startTheta<< endl);
    if (endTheta         == -11111.)
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "Invalid value for endTheta : "<<endTheta<< endl);
    if (phiOffset         == -11111.)
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "Invalid value for phiOffset : "<<phiOffset<< endl);
    MMBLOG_FILE_FUNC_LINE(INFO, " Validation passed.           "<< endl);
}
double Spiral::harmonicThetaAdjustment(double inputPhi){
    double outputThetaAdjustment =0;
    for (int i =0; i < frequencyPhaseAmplitudeVector.size(); i++){
        outputThetaAdjustment +=  frequencyPhaseAmplitudeVector[i].amplitude * sin(inputPhi * frequencyPhaseAmplitudeVector[i].frequency + frequencyPhaseAmplitudeVector[i].phase);
        MMBLOG_FILE_FUNC_LINE(DEBUG, " frequencyPhaseAmplitudeVector i , frequency, amplitude, phase =" << i  <<" "<<  frequencyPhaseAmplitudeVector[i].frequency<<" "<<  frequencyPhaseAmplitudeVector[i].amplitude<<" " <<frequencyPhaseAmplitudeVector[i].phase<<" "    << endl);
        MMBLOG_FILE_FUNC_LINE(DEBUG, " outputThetaAdjustment = "<<outputThetaAdjustment<<endl);
    }
    MMBLOG_FILE_FUNC_LINE(DEBUG, " outputThetaAdjustment = "<<outputThetaAdjustment<<endl);
    return outputThetaAdjustment;
}                               	


void Spiral::writeDnaSpiralCommandfile(MonoAtomsContainer &monoAtomsContainer)
{
    validate();
    Vec3 priorXYZ(-9999.9, -9999.9, -9999.9);
    priorXYZ = Vec3(radius* sin(startTheta)*cos(phiOffset) , radius* sin(startTheta)*sin(phiOffset), radius* cos(startTheta));
    priorXYZ += center;
    MMBLOG_FILE_FUNC_LINE(INFO, "radius "<<radius<<endl);
    int n = 0; // counter   
    MMBLOG_FILE_FUNC_LINE(INFO, "priorXYZ "<<priorXYZ<<endl);
    FILE * spiralPdbFile;
    spiralPdbFile      = fopen (spiralPdbFileName,"w"); // file name should be specified by user. Was .."spiral.pdb","w"); 
    MMBLOG_FILE_FUNC_LINE(INFO, " Opened  spiralPdbFileName = " <<spiralPdbFileName                                   <<endl);
    fprintf (spiralPdbFile,"ATOM  %5d MG2+ MG  Z%4d    %8.3f%8.3f%8.3f \n",n,n,priorXYZ[0]*10, priorXYZ[1]*10,priorXYZ[2]*10  ); // Converting to Ångströms
    //double helixAdvancePerBasePair = 0.34 ;  // in nm // should be user specified.
    //double deltaTheta =  SimTK::Pi  / 2000;
    MMBLOG_FILE_FUNC_LINE(INFO, std::endl);
    MMBLOG_FILE_FUNC_LINE(INFO, " priorXYZ "<<priorXYZ<<endl);
    // might be good to run this just to confirm:
    MMBLOG_FILE_FUNC_LINE(INFO, "starting theta provided by user : "<<startTheta<<" compared to that gotten by trigonometry after converting spherical to cartesian and back to spherical: "<<thetaFromXYZ(priorXYZ, center)<<endl);
    MMBLOG_FILE_FUNC_LINE(INFO, "startTheta = "<< startTheta << ""<<endl);
    double currentTheta = -11111.1;
    // In the default spherical spiral, phi is just theta times a constant. However we need to be able to rotate the sphere and key the spiral wherever we want. So we need a phiOffset.
    double phiOffsetFromStartingTheta = phiFromTheta(startTheta, interStrandDistance, radius, 0.0); // use offset = 0 to retreive the original, non-offset phi
    // 21.05.04 SCF Here  "- phiOffsetFromStartingTheta" means that the first ion should be at phi=phiOffset.
    MMBLOG_FILE_FUNC_LINE(DEBUG, "phiOffset "<<phiOffset<<endl);
    MMBLOG_FILE_FUNC_LINE(DEBUG, "phiOffsetFromStartingTheta "<<phiOffsetFromStartingTheta<<endl);
    // Here we are DECLARING phiOffset, and immediately trying to use it!!! Why are we even being allowed to redeclare here???
    //double phiOffset = phiOffset - phiOffsetFromStartingTheta;
    MMBLOG_FILE_FUNC_LINE(DEBUG, "phiOffset "<<phiOffset<<endl);
    MMBLOG_FILE_FUNC_LINE(INFO, "priorXYZ "<<priorXYZ<<endl);
    double priorTheta = startTheta ; // thetaFromXYZ(priorXYZ, center) ;
    n = 1;
    std::string tetherCommands("");
    stringstream tetherCommandStream("");
    //                  chain ID, first residue #, number of ions (we start with an empty vector), name of ions.
    MonoAtoms monoAtoms(String("Z"),ResidueID(1),0,String("Mg+2"));  
    while (currentTheta < endTheta){
        MMBLOG_FILE_FUNC_LINE(INFO, "priorXYZ "<<priorXYZ<<endl);
        double priorPhi = phiFromXYZ(priorXYZ, center) ;
        MMBLOG_FILE_FUNC_LINE(INFO, "priorXYZ "<<priorXYZ<<endl);
        MMBLOG_FILE_FUNC_LINE(INFO, "priorPhi "<<priorPhi<<endl);
        double deltaPhi = deltaPhiFromThetaInterHelicalDistanceSphericalRadiusAndHelicalArcLength(priorTheta, interStrandDistance, radius, helixAdvancePerBasePair);
        MMBLOG_FILE_FUNC_LINE(INFO, "deltaPhi "<<deltaPhi<<endl);
        double deltaTheta = thetaFromPhi(deltaPhi, interStrandDistance, radius, 0.0); // SCF no need to include phiOffset here. This just gives us deltaTheta given deltaPhi
        MMBLOG_FILE_FUNC_LINE(INFO, "radius "<<radius<<endl);
        MMBLOG_FILE_FUNC_LINE(INFO, "deltaTheta = "<<deltaTheta<<" .. should be positive and tiny"<<endl);
        currentTheta = priorTheta + deltaTheta;
        MMBLOG_FILE_FUNC_LINE(INFO, "currentTheta = priorTheta + deltaTheta : "<<currentTheta<<" = "<<priorTheta<<" + "<<deltaTheta<<" "<<endl);
        priorTheta = currentTheta;
        MMBLOG_FILE_FUNC_LINE(INFO, "priorTheta "<<priorTheta<<endl);
        double currentPhi = priorPhi + deltaPhi ; // phiFromTheta(currentTheta, interStrandDistance, radius, phiOffset);
        MMBLOG_FILE_FUNC_LINE(INFO, "currentPhi "<<currentPhi<<endl);
        MMBLOG_FILE_FUNC_LINE(DEBUG, "phiOffset "<<phiOffset<<endl);
        Vec3 currentXYZ (
            radius * sin(currentTheta) * cos (phiFromTheta(currentTheta, interStrandDistance, radius, phiOffset)), 
            radius * sin(currentTheta) * sin (phiFromTheta(currentTheta, interStrandDistance, radius, phiOffset)),
            radius * cos(currentTheta)
            );
	double currentAdjustedTheta = currentTheta+harmonicThetaAdjustment(phiFromTheta(currentTheta, interStrandDistance, radius, phiOffset));
        Vec3 currentAdjustedXYZ (
            radius * sin(currentAdjustedTheta) * cos (phiFromTheta(currentTheta, interStrandDistance, radius, phiOffset)), 
            radius * sin(currentAdjustedTheta) * sin (phiFromTheta(currentTheta, interStrandDistance, radius, phiOffset)),
            radius * cos(currentAdjustedTheta)
            );
        
        MMBLOG_FILE_FUNC_LINE(INFO, "radius "<<radius<<endl);
        MMBLOG_FILE_FUNC_LINE(INFO, "interStrandDistance "<<interStrandDistance<<endl);
        currentXYZ += center;
        currentAdjustedXYZ += center;
        MMBLOG_FILE_FUNC_LINE(INFO, "currentXYZ "<<currentXYZ<<endl);
        MMBLOG_FILE_FUNC_LINE(INFO, "phi from the just-updated currentXYZ = "<< phiFromXYZ(currentXYZ, center)<<endl);
        MMBLOG_FILE_FUNC_LINE(INFO, " About to write to spiralPdbFileName = " <<spiralPdbFileName                                   <<endl);
        fprintf (spiralPdbFile,"ATOM  %5d MG2+ MG  Z%4d    %8.3f%8.3f%8.3f \n",n,n,currentAdjustedXYZ[0]*10, currentAdjustedXYZ[1]*10,currentAdjustedXYZ[2]*10  ); // Converting to Ångströms
	monoAtoms.addMonoAtom(Vec3(currentAdjustedXYZ[0], currentAdjustedXYZ[1],currentAdjustedXYZ[2])); // provide coords in nm.
	// I believe the above is enough to instantiate an ion in the MonoAtoms object. The question is, is this being passed on to the right systems in MMB?
        MMBLOG_FILE_FUNC_LINE(INFO, endl);
        MMBLOG_FILE_FUNC_LINE(INFO, "MMB-command: readAtStage "<<n<<endl);
        tetherCommandStream<<"readAtStage "<<n<<std::endl;
        // base-pair-at-origin.pdb contains a single base pair, with axis perpendicular to Z-axis.
        // We had a problem with a rotation command which was spitting out tiny values in scientific notation, which MMB later had a hard time parsing. Solution is to used fixed-format, with 6 places precision:
        tetherCommandStream.setf(std::ios_base::fixed, std::ios_base::floatfield);
        tetherCommandStream.precision(6);  
        tetherCommandStream<<"previousFrameFileName    base-pair."<<n<<".pdb   "<<std::endl;

        tetherCommandStream<<"DNA A "<<n<<" A "<<std::endl;
        tetherCommandStream<<"DNA B "<<9999-n<<" T "<<std::endl;
        tetherCommandStream<<"#renumberBiopolymerResidues A "<<n<<std::endl;
        tetherCommandStream<<"#renumberBiopolymerResidues B "<<9999-n<<std::endl;
        tetherCommandStream<<"mobilizer Rigid "<<std::endl;
        tetherCommandStream<<"initialDisplacement A "<<  currentAdjustedXYZ[0] <<" "<<  currentAdjustedXYZ[1]  <<" "<< currentAdjustedXYZ[2] <<std::endl;
        tetherCommandStream<<"initialDisplacement B "<<  currentAdjustedXYZ[0] <<" "<<  currentAdjustedXYZ[1]  <<" "<< currentAdjustedXYZ[2] <<std::endl;
        tetherCommandStream<<"rotation A  Z "<<  (SimTK::Pi *  2) / 10 * n    <<std::endl; // first, rotate the base pair by 360/10 degrees * number of base pairs.
        tetherCommandStream<<"rotation B  Z "<<  (SimTK::Pi *  2) / 10 * n    <<std::endl; // first, rotate the base pair by 360/10 degrees * number of base pairs.
        tetherCommandStream<<"rotation A X "<<  -atan(currentAdjustedTheta / phiFromTheta(currentTheta, interStrandDistance, radius, phiOffset*0.0)) - (SimTK::Pi / 2)    <<std::endl;  // Now, slope it so if follows the tangential slope of the helix.
        tetherCommandStream<<"rotation B X "<<  -atan(currentAdjustedTheta / phiFromTheta(currentTheta, interStrandDistance, radius, phiOffset*0.0)) - (SimTK::Pi / 2)   <<std::endl;  // Now, slope it so if follows the tangential slope of the helix.
        tetherCommandStream<<"rotation A Y "<<  -(SimTK::Pi / 2) +  currentAdjustedTheta  <<std::endl ; // tilt up 
        tetherCommandStream<<"rotation B Y "<<  -(SimTK::Pi / 2) +  currentAdjustedTheta  <<std::endl ; // tilt up 
        tetherCommandStream<<"rotation A Z "<<  phiFromTheta(currentTheta, interStrandDistance, radius, phiOffset)   <<std::endl;
        tetherCommandStream<<"rotation B Z "<<  phiFromTheta(currentTheta, interStrandDistance, radius, phiOffset)   <<std::endl;
        tetherCommandStream<<"readBlockEnd"<<std::endl;


        MMBLOG_FILE_FUNC_LINE(INFO, "MMB-command: tetherToGround A "<<n<<" N1 "<<currentXYZ[0]<<" "<<currentXYZ[1]<<" "<<currentXYZ[2] << " .5 30.0 "<<endl);
        MMBLOG_FILE_FUNC_LINE(INFO, "MMB-command: readBlockEnd "<<endl);
        MMBLOG_FILE_FUNC_LINE(INFO, "priorXYZ "<<priorXYZ<<endl);
        priorXYZ= currentXYZ;
        MMBLOG_FILE_FUNC_LINE(INFO, "priorXYZ "<<priorXYZ<<endl);
        n++;
    }
    fclose(spiralPdbFile); 
    monoAtoms.renumberPdbResidues(); // monoAtoms were above added one by one  with no residue numbers set. Now we fix that.
    MMBLOG_FILE_FUNC_LINE(INFO, " closed  spiralPdbFileName = " <<spiralPdbFileName                                   <<endl);
    MMBLOG_FILE_FUNC_LINE(INFO, endl);
    std::string commonSpiralCommandsAdjusted = std::regex_replace( commonSpiralCommands, std::regex(std::string("ZZZZ")), std::to_string(n-1 ) ); 
    commonSpiralCommandsAdjusted = std::regex_replace( commonSpiralCommandsAdjusted, std::regex(std::string("            ")), std::string("") );  // Now get rid of extra whitespace
    ofstream spiralCommandsFile2(spiralCommandsFileName); // was "commands.spiral.dat") ;
    spiralCommandsFile2<<commonSpiralCommandsAdjusted<<std::endl;
    spiralCommandsFile2<<tetherCommandStream.str()<<std::endl;
    spiralCommandsFile2.close();
    MMBLOG_FILE_FUNC_LINE(INFO, endl);
    monoAtomsContainer.addMonoAtoms(monoAtoms);
} // of procedure

void Spiral::writeSyntax()
{
    MMBLOG_FILE_FUNC_LINE(ALWAYS, "sphericalHelix creates a spherical spiral of MG2+ ions. Eventually it will be adaptive to the density. You need to provide the spherical center (3D), the spherical radius, the pitch (inter-duplex distance), all in nm. You need start and end theta in rads. There are other optional parameters. "<<endl);
    MMBLOG_FILE_FUNC_LINE(ALWAYS, "Syntax: "<<endl);
    MMBLOG_FILE_FUNC_LINE(ALWAYS, "To specify the center point of the sphere, in nm: "<<endl);
    MMBLOG_FILE_FUNC_LINE(ALWAYS, "sphericalHelix center <X>  <Y> <Z> "<<endl);
    MMBLOG_FILE_FUNC_LINE(ALWAYS, "To specify the radius of the sphere: "<<endl);
    MMBLOG_FILE_FUNC_LINE(ALWAYS, "sphericalHelix radius <radius, in nm>  "<<endl);
    MMBLOG_FILE_FUNC_LINE(ALWAYS, "To specify the separation between consecutive DNA duplexes : "<<endl);
    MMBLOG_FILE_FUNC_LINE(ALWAYS, "sphericalHelix interStrandDistance  <distance, in nm>  "<<endl);
    MMBLOG_FILE_FUNC_LINE(ALWAYS, "To specify the start theta (the angle from the 'north pole': "<<endl);
    MMBLOG_FILE_FUNC_LINE(ALWAYS, "sphericalHelix startTheta <angle, in rads>  "<<endl);
    MMBLOG_FILE_FUNC_LINE(ALWAYS, "To specify the end   theta (the angle from the 'north pole': "<<endl);
    MMBLOG_FILE_FUNC_LINE(ALWAYS, "sphericalHelix   endTheta <angle, in rads>  "<<endl);

    MMBLOG_FILE_FUNC_LINE(ALWAYS, "Next specify the offset in phi (the angle about the +Z axis).  Phi = 0 in the +X half of the XZ plane, and increases following the right-hand rule about the +Z-axis.  "<<endl);
    MMBLOG_FILE_FUNC_LINE(ALWAYS, "To specify the offset in phi : "<<endl);
    MMBLOG_FILE_FUNC_LINE(ALWAYS, "sphericalHelix phiOffset <angle, in rads>  "<<endl);
    MMBLOG_FILE_FUNC_LINE(ALWAYS, "This file will contain MG ions indicating the center of each base pair in the spiral. You will use this as a rough draft to check your spiral parameters: "<<endl);
    MMBLOG_FILE_FUNC_LINE(ALWAYS, "sphericalHelix spiralPdbFileName <string>  "<<endl);
    MMBLOG_FILE_FUNC_LINE(ALWAYS, "This file will contain the commands to place idealized base pairs at their positions in the ideal spiral you generated. You would rerun MMB with this as your input command file: "<<endl);
    MMBLOG_FILE_FUNC_LINE(ALWAYS, "sphericalHelix spiralCommandsFileName <string>  "<<endl);
    MMBLOG_FILE_FUNC_LINE(ALWAYS, "This will modulate the theta position as a functoin of phi. Use a frequency multiplier of 1 if you want this to oscillate by 2*pi over phi sweep of 2*pi. Or 2 if you want to oscillate by 4*pi, etc.  "<<endl);
    MMBLOG_FILE_FUNC_LINE(ALWAYS, "The effect is additive. Every time you call this, one element will be added to the frequencyPhaseAmplitudeVector, and all elements will be applied.            "<<endl);
    MMBLOG_FILE_FUNC_LINE(ALWAYS, "The output theta will have a correction added which looks like: amplitude * sin( multiplier * phi + phase) "<<endl);
    MMBLOG_FILE_FUNC_LINE(ALWAYS, "sphericalHelix frequencyPhaseAmplitude <frequency multiplier> <phase, rads> <amplitude, rads>  "<<endl);
    MMBLOG_FILE_FUNC_LINE(ALWAYS, "If you want to remove all elements from the frequencyPhaseAmplitudeVector, issue: "<<endl);
    MMBLOG_FILE_FUNC_LINE(ALWAYS, "sphericalHelix frequencyPhaseAmplitude clear                                                   "<<endl);
    MMBLOG_FILE_FUNC_LINE(ALWAYS, "And finally, to create the helix, issue:                        "<<endl);
    MMBLOG_FILE_FUNC_LINE(ALWAYS, "sphericalHelix writeCommands "<<endl);
    MMBLOG_FILE_FUNC_LINE(ALWAYS, "The above should be the last sphericalHelix command you issue. "<<endl);
}

    

// This function figures out which polymorphism of parseInput should be used, based on number of parameters provided, be it 0,1,or 3.
void Spiral::parseInput(const map<const String,double> & userVariables, String parameterName, String parameterValue1,String parameterValue2, String parameterValue3,MonoAtomsContainer &monoAtomsContainer){
    if (parameterValue3 != ""){
        if (parameterValue2 == "") MMBLOG_FILE_FUNC_LINE(CRITICAL, "Wrong number of parameters! Please check your syntax.       "<<endl); // if parameterValue3 is specified, then parameterValue2 must also be.
        parseInput(parameterName,myAtoF(userVariables,(parameterValue1).c_str()),myAtoF(userVariables,(parameterValue2).c_str()),myAtoF(userVariables,(parameterValue3).c_str()));
    } // Use is probably specifying a Vec3, namely a center.
    else if (parameterValue1 != ""){
        if (parameterName   == "") MMBLOG_FILE_FUNC_LINE(CRITICAL, " You are misusing the sphericalHelix family of commands! Please specify a parameter or command.       "<<endl);
        parseInput(userVariables,parameterName,parameterValue1);} // specifying a parameter with a certain single value
    else if (parameterName   != ""){parseInput(parameterName, monoAtomsContainer);} // specifying a command.                           
    else {
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "You are misusing the sphericalHelix family of commands! Please specify a parameter or command. "<<endl);
    
    }
}	

void Spiral::parseInput(String commandName, MonoAtomsContainer & monoAtomsContainer){
    if (commandName == "writeCommands"){
        writeDnaSpiralCommandfile(monoAtomsContainer);
    } else {
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "Non recognizable input: " << commandName << endl);
    }
}
void Spiral::parseInput(const  map<const String,double> & userVariables,String parameterName, String parameterValue){
    if ((parameterName).compare("radius") ==0)    {
        radius = myAtoF(userVariables,(parameterValue).c_str());
        return;
    } else if ((parameterName).compare("interStrandDistance") ==0)    {
        interStrandDistance = myAtoF(userVariables,(parameterValue).c_str());
        return;
    } else if ((parameterName).compare("startTheta") ==0)    {
        startTheta = myAtoF(userVariables,(parameterValue).c_str());
        return;
    } else if ((parameterName).compare("phiOffset") ==0)    {
        phiOffset = myAtoF(userVariables,(parameterValue).c_str());
        return;
    } else if ((parameterName).compare("endTheta") ==0)    {
        endTheta = myAtoF(userVariables,(parameterValue).c_str());
        return;
    } else if ((parameterName).compare("helixAdvancePerBasePair") ==0)    {
        helixAdvancePerBasePair = myAtoF(userVariables,(parameterValue).c_str());
        return;
    } else if ((parameterName).compare("spiralPdbFileName")==0){ // was "spiral.pdb"
        spiralPdbFileName = parameterValue;	    
        MMBLOG_FILE_FUNC_LINE(INFO, " Setting spiralPdbFileName = " <<spiralPdbFileName                                   <<endl);
        return;
    } else if ((parameterName).compare("spiralCommandsFileName")==0){ // was "commands.spiral.dat") ;
        spiralCommandsFileName = parameterValue;	    
        return;
    } else if ((parameterName).compare("frequencyPhaseAmplitude")==0){ 
        if ( parameterValue == "clear"){ // ATM there is only one supported command for frequencyPhaseAmplitude, and that is "clear". That is in addition to the setting of frequency multiplier, phase, and amplitude, which is treated separately.
            frequencyPhaseAmplitudeVector.clear();		
	} else {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "You have provided an unsupported parameter for "<<parameterName<<" : "<<parameterValue<<". Check your syntax.    "<<endl);
        }	    
        return;
    } else {
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "You are misusing the sphericalHelix family of commands! Check your syntax and try again. "<<endl);
    }

}
    vector<frequencyPhaseAmplitude> frequencyPhaseAmplitudeVector; // This holds all the harmonic adjustments to be made to the spherical helix.
void Spiral::parseInput(String parameterName, double X, double Y, double Z ){
    if ((parameterName).compare("center") ==0)    {
        center = Vec3(X,Y,Z);
    } else if ((parameterName).compare("frequencyPhaseAmplitude") ==0)    {
        frequencyPhaseAmplitude myFrequencyPhaseAmplitude;
        myFrequencyPhaseAmplitude.frequency = X;
        myFrequencyPhaseAmplitude.phase =     Y;
        myFrequencyPhaseAmplitude.amplitude = Z;	
	frequencyPhaseAmplitudeVector.push_back(myFrequencyPhaseAmplitude);
    } else {
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "Non recognizable input: " << parameterName << endl);
    }
}

