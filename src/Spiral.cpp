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
#include "Utils.h"  
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
    MMBLOG_FILE_FUNC_LINE(DEBUG, "myPhi = " << myPhi<<endl);
    return myPhi;
}
double thetaFromXYZ( const Vec3 myPoint, const  Vec3 sphericalCenter) {
    double xDist = myPoint[0]-sphericalCenter[0];
    double yDist = myPoint[1]-sphericalCenter[1];
    double zDist = myPoint[2]-sphericalCenter[2];
    double sliceRadius = sqrt(xDist*xDist+yDist*yDist);
    double computedSphericalRadius = sqrt(xDist*xDist+yDist*yDist+zDist*zDist);
    MMBLOG_FILE_FUNC_LINE(DEBUG, "slice radius = " << sliceRadius<<endl);
    double theta = asin(sliceRadius / computedSphericalRadius);
    MMBLOG_FILE_FUNC_LINE(DEBUG, "computedSphericalRadius = " << computedSphericalRadius<<endl);
    MMBLOG_FILE_FUNC_LINE(DEBUG, "theta = " << theta<<endl);
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

// Returns -1 for right handed spirals, +1 for left handed. This is for changing the sign of phi and related quantities.
double negativeWhenRightHanded(bool spiralIsRightHanded){
    return ( -2.0 * spiralIsRightHanded + 1);
}

double deltaPhiFromThetaInterHelicalDistanceSphericalRadiusAndHelicalArcLength( const double myTheta, const double myInterHelicalDistance, const  double mySphericalRadius, const double myArcLength , const bool spiralIsRightHanded){
    double myHelicalPitch = myInterHelicalDistance*sin(myTheta);
    double myCylindricalRadius = mySphericalRadius * sin(myTheta);
    //( -2.0 * spiralIsRightHanded + 1) returns +1 for left handed, -1 for right handed. So if it is right handed, we change sign of delta-phi.
    return  negativeWhenRightHanded(spiralIsRightHanded)*deltaPhiFromCylindricalRadiusHelicalPitchAndHelicalArcLength(myCylindricalRadius, myHelicalPitch, myArcLength);
}
double thetaFromPhi( const double phi, const double myInterHelicalDistance, const double mySphericalRadius, const double myPhiOffset, const bool spiralIsRightHanded){
    MMBLOG_FILE_FUNC_LINE(DEBUG, " locally, phiOffset = "<<myPhiOffset<<endl);
    return negativeWhenRightHanded(spiralIsRightHanded)*(phi - myPhiOffset) / (2*SimTK::Pi/myInterHelicalDistance*mySphericalRadius) ;
}
double phiFromTheta(const double theta,const double myInterHelicalDistance, const double mySphericalRadius, const double myPhiOffset, const bool spiralIsRightHanded){
    MMBLOG_FILE_FUNC_LINE(DEBUG, " locally, phiOffset = "<<myPhiOffset<<endl);
    return negativeWhenRightHanded(spiralIsRightHanded)*theta * (2*SimTK::Pi/myInterHelicalDistance*mySphericalRadius) + myPhiOffset ;
}


Spiral::Spiral(){clear();}
void   Spiral::clear(){
    frequencyPhaseAmplitudeVector.clear();    
    spiralPdbFileName = String("NOT-SET"); // was "spiral.pdb"
    MMBLOG_FILE_FUNC_LINE(INFO, " Setting spiralPdbFileName = " <<spiralPdbFileName                                   <<endl);
    spiralCommandsFileName = String("NOT-SET"); // was "commands.spiral.dat"
    center            = Vec3(-11111.,-11111.,-11111.);
    radius            = -11111.;
    pitch=-11111.;
    startTheta        = -11111.;
    endTheta          = -11111.;
    phiOffset         = -11111.;
    cylinderHeight    = -11111.;
    spiralIsRightHanded = 0; // defaults to left-handed
    double helixAdvancePerBasePair = -11111. ;  // in nm // should be user specified.
    String chainID    = "Z";
    MMBLOG_FILE_FUNC_LINE(DEBUG, " Setting chainID  = " << chainID                       <<endl);
    geometry = geometryEnum::sphere;

}
void   Spiral::validate(){
    //frequencyPhaseAmplitudeVector.clear();    
    if (spiralPdbFileName == String("NOT-SET"))
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "Invalid value for spiralPdbFileName : "<<spiralPdbFileName<< endl);
    if (spiralCommandsFileName == String("NOT-SET"))
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "Invalid value for spiralCommandsFileName : "<<spiralCommandsFileName<< endl);
    if (radius         == -11111.)
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "Invalid value for radius : "<<radius<< endl);
    if (pitch         == -11111.)
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "Invalid value for pitch : "<<pitch<< endl);
    if (geometry  == geometryEnum::sphere) {
        if (startTheta         == -11111.)
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "Invalid value for startTheta : "<<startTheta<< endl);
        if (endTheta         == -11111.)
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "Invalid value for endTheta : "<<endTheta<< endl);
    } else if (geometry  == geometryEnum::cylinder){
        if (cylinderHeight <=0){
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "Invalid value for cylinderHeight : "<<cylinderHeight<< endl);
	}
    }
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
double cylindricalSpiralSlopeAngle(const double pitch, const double radius){
    return atan(pitch/radius/2/SimTK::Pi);
}
void Spiral::writeCylindricalSpiralCommandFile(MonoAtomsContainer &monoAtomsContainer)
{
    MMBLOG_FILE_FUNC_LINE(DEBUG, "  "<<endl);
    validate();
    double deltaZScalar =  sin(cylindricalSpiralSlopeAngle(pitch,radius))*helixAdvancePerBasePair ;
    Vec3 deltaZ = Vec3(0,0,1) * sin(cylindricalSpiralSlopeAngle(pitch,radius))*helixAdvancePerBasePair ;
                        // arc length projected to "floor"                                     // divide by radius
    double deltaPhi = cos(cylindricalSpiralSlopeAngle(pitch,radius))*helixAdvancePerBasePair / radius;
    double currentZ ;
    Vec3   currentXYZ (0,0,0) ;
    double currentPhi = phiOffset;
    int n=1;
    stringstream tetherCommandStream("");
    MonoAtoms monoAtoms(chainID, ResidueID(1), 0, String("Mg+2"));
    for (currentZ = 0.; currentZ <= cylinderHeight; currentZ += deltaZScalar){
        currentXYZ[0] = center[0] + cos(currentPhi)*radius;
        currentXYZ[1] = center[1] + sin(currentPhi)*radius;
        currentXYZ[2] = center[2] + currentZ;                       
	monoAtoms.addMonoAtom(currentXYZ); // provide coords in nm.


        tetherCommandStream<<"readAtStage "<<n<<std::endl;
        // base-pair-at-origin.pdb contains a single base pair, with axis perpendicular to Z-axis.
        // We had a problem with a rotation command which was spitting out tiny values in scientific notation, which MMB later had a hard time parsing. Solution is to used fixed-format, with 6 places precision:
        tetherCommandStream.setf(std::ios_base::fixed, std::ios_base::floatfield);
        tetherCommandStream.precision(6);  
        tetherCommandStream<<"previousFrameFileName    base-pair."<<n<<".pdb   "<<std::endl;

        tetherCommandStream<<"DNA A "<<n<<" A "<<std::endl;
        tetherCommandStream<<"DNA B "<<9999-n<<" T "<<std::endl;
        tetherCommandStream<<"mobilizer Rigid "<<std::endl;
        tetherCommandStream<<"initialDisplacement A "<<  currentXYZ[0] <<" "<<  currentXYZ[1]  <<" "<< currentXYZ[2] <<std::endl;
        tetherCommandStream<<"initialDisplacement B "<<  currentXYZ[0] <<" "<<  currentXYZ[1]  <<" "<< currentXYZ[2] <<std::endl;
        tetherCommandStream<<"rotation A  Z "<<  (SimTK::Pi *  2) / 10 * n    <<std::endl; // first, rotate the base pair by 360/10 degrees * number of base pairs.
        tetherCommandStream<<"rotation B  Z "<<  (SimTK::Pi *  2) / 10 * n    <<std::endl; // first, rotate the base pair by 360/10 degrees * number of base pairs.
        tetherCommandStream<<"rotation A X "<<  atan( pitch / radius / SimTK::Pi / 2)*(-negativeWhenRightHanded(spiralIsRightHanded))    <<std::endl;  // Now, slope it so if follows the tangential slope of the helix. // here we are assuming right handed helix,  then  there is a sign change if it is actually left handed.
        tetherCommandStream<<"rotation A X "<<  1.0*spiralIsRightHanded*SimTK::Pi     <<std::endl;  // If the spiral is right-handed, we need an additional 180-degree rotation
        tetherCommandStream<<"rotation B X "<<  atan( pitch / radius / SimTK::Pi / 2)*(-negativeWhenRightHanded(spiralIsRightHanded))    <<std::endl;  // Now, slope it so if follows the tangential slope of the helix. // here we are assuming right handed helix,  then  there is a sign change if it is actually left handed.
        tetherCommandStream<<"rotation B X "<<  1.0*spiralIsRightHanded*SimTK::Pi     <<std::endl;  // If the spiral is right-handed, we need an additional 180-degree rotation
        tetherCommandStream<<"rotation A Z "<<  currentPhi   <<std::endl;
        tetherCommandStream<<"rotation B Z "<<  currentPhi   <<std::endl;
        tetherCommandStream<<"readBlockEnd"<<std::endl;


	currentPhi += deltaPhi;
	n++;
    }
    std::string commonSpiralCommandsAdjusted = std::regex_replace( commonSpiralCommands, std::regex(std::string("ZZZZ")), std::to_string(n-1 ) ); 
    commonSpiralCommandsAdjusted = std::regex_replace( commonSpiralCommandsAdjusted, std::regex(std::string("            ")), std::string("") );  // Now get rid of extra whitespace
    ofstream spiralCommandsFile2(spiralCommandsFileName); // was "commands.spiral.dat") ;
    spiralCommandsFile2<<commonSpiralCommandsAdjusted<<std::endl;
    spiralCommandsFile2<<tetherCommandStream.str()<<std::endl;
    spiralCommandsFile2.close();
    monoAtoms.renumberPdbResidues(); // monoAtoms were above added one by one  with no residue numbers set. Now we fix that.
    monoAtomsContainer.addMonoAtoms(monoAtoms);
} // of writeCylindricalSpiralCommandFile
/*
*/
void Spiral::writeSphericalSpiralCommandFile(MonoAtomsContainer &monoAtomsContainer)
{
    MMBLOG_FILE_FUNC_LINE(DEBUG, "  "<<endl);
    validate();
    Vec3 priorXYZ(-9999.9, -9999.9, -9999.9);
    priorXYZ = Vec3(radius* sin(startTheta)*cos(phiOffset) , radius* sin(startTheta)*sin(phiOffset), radius* cos(startTheta));
    priorXYZ += center;
    MMBLOG_FILE_FUNC_LINE(DEBUG, "radius "<<radius<<endl);
    int n = 0; // counter   
    MMBLOG_FILE_FUNC_LINE(DEBUG, "priorXYZ "<<priorXYZ<<endl);
    //FILE * spiralPdbFile;
    //spiralPdbFile      = fopen (spiralPdbFileName,"w"); // file name should be specified by user. Was .."spiral.pdb","w"); 
    double currentTheta = -11111.1;
    // In the default spherical spiral, phi is just theta times a constant. However we need to be able to rotate the sphere and key the spiral wherever we want. So we need a phiOffset.
    double phiOffsetFromStartingTheta = phiFromTheta(startTheta, pitch, radius, 0.0, spiralIsRightHanded); // use offset = 0 to retreive the original, non-offset phi
    // 21.05.04 SCF Here  "- phiOffsetFromStartingTheta" means that the first ion should be at phi=phiOffset.
    double priorTheta = startTheta ; // thetaFromXYZ(priorXYZ, center) ;
    n = 1;
    std::string tetherCommands("");
    stringstream tetherCommandStream("");
    //                  chain ID, first residue #, number of ions (we start with an empty vector), name of ions.
    MMBLOG_FILE_FUNC_LINE(DEBUG, "chainID "<<chainID<<endl);
    MonoAtoms monoAtoms(chainID,ResidueID(1), 0, String("Mg+2"));
    while (currentTheta < endTheta){
        double priorPhi = phiFromXYZ(priorXYZ, center) ;
        double deltaPhi = deltaPhiFromThetaInterHelicalDistanceSphericalRadiusAndHelicalArcLength(priorTheta, pitch, radius, helixAdvancePerBasePair,spiralIsRightHanded);
        double deltaTheta = thetaFromPhi(deltaPhi, pitch, radius, 0.0, spiralIsRightHanded); // SCF no need to include phiOffset here. This just gives us deltaTheta given deltaPhi
        MMBLOG_FILE_FUNC_LINE(DEBUG, "radius "<<radius<<endl);
        MMBLOG_FILE_FUNC_LINE(DEBUG, "deltaTheta = "<<deltaTheta<<" .. should be positive and tiny"<<endl);
        currentTheta = priorTheta + deltaTheta;
        priorTheta = currentTheta;
        double currentPhi = priorPhi + deltaPhi ; // phiFromTheta(currentTheta, pitch, radius, phiOffset);
        Vec3 currentXYZ (
            radius * sin(currentTheta) * cos (phiFromTheta(currentTheta, pitch, radius, phiOffset, spiralIsRightHanded)), 
            radius * sin(currentTheta) * sin (phiFromTheta(currentTheta, pitch, radius, phiOffset,spiralIsRightHanded)),
            radius * cos(currentTheta)
            );
	double currentAdjustedTheta = currentTheta+harmonicThetaAdjustment(phiFromTheta(currentTheta, pitch, radius, phiOffset,spiralIsRightHanded));
        Vec3 currentAdjustedXYZ (
            radius * sin(currentAdjustedTheta) * cos (phiFromTheta(currentTheta, pitch, radius, phiOffset,spiralIsRightHanded)), 
            radius * sin(currentAdjustedTheta) * sin (phiFromTheta(currentTheta, pitch, radius, phiOffset,spiralIsRightHanded)),
            radius * cos(currentAdjustedTheta)
            );
        
        MMBLOG_FILE_FUNC_LINE(DEBUG, "radius "<<radius<<endl);
        MMBLOG_FILE_FUNC_LINE(DEBUG, "pitch "<<pitch<<endl);
        currentXYZ += center;
        currentAdjustedXYZ += center;
        MMBLOG_FILE_FUNC_LINE(DEBUG, "currentXYZ "<<currentXYZ<<endl);
        MMBLOG_FILE_FUNC_LINE(DEBUG, "phi from the just-updated currentXYZ = "<< phiFromXYZ(currentXYZ, center)<<endl);
        MMBLOG_FILE_FUNC_LINE(DEBUG, " About to write to spiralPdbFileName = " <<spiralPdbFileName                                   <<endl);
        //fprintf (spiralPdbFile,"ATOM  %5d MG2+ MG  Z%4d    %8.3f%8.3f%8.3f \n",n,n,currentAdjustedXYZ[0]*10, currentAdjustedXYZ[1]*10,currentAdjustedXYZ[2]*10  ); // Converting to Ångströms
	monoAtoms.addMonoAtom(Vec3(currentAdjustedXYZ[0], currentAdjustedXYZ[1],currentAdjustedXYZ[2])); // provide coords in nm.
	// I believe the above is enough to instantiate an ion in the MonoAtoms object. The question is, is this being passed on to the right systems in MMB?
        MMBLOG_FILE_FUNC_LINE(DEBUG, endl);
        MMBLOG_FILE_FUNC_LINE(DEBUG, "MMB-command: readAtStage "<<n<<endl);
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
        tetherCommandStream<<"rotation A X "<<  -atan(currentAdjustedTheta / phiFromTheta(currentTheta, pitch, radius, phiOffset*0.0,spiralIsRightHanded)) - (SimTK::Pi / 2)    <<std::endl;  // Now, slope it so if follows the tangential slope of the helix.
        tetherCommandStream<<"rotation A X "<<  1.0*spiralIsRightHanded*SimTK::Pi     <<std::endl;  // If the spiral is right-handed, we need an additional 180-degree rotation
        tetherCommandStream<<"rotation B X "<<  -atan(currentAdjustedTheta / phiFromTheta(currentTheta, pitch, radius, phiOffset*0.0,spiralIsRightHanded)) - (SimTK::Pi / 2)   <<std::endl;  // Now, slope it so if follows the tangential slope of the helix.
        tetherCommandStream<<"rotation B X "<<  1.0*spiralIsRightHanded*SimTK::Pi     <<std::endl;  // If the spiral is right-handed, we need an additional 180-degree rotation
        tetherCommandStream<<"rotation A Y "<<  -(SimTK::Pi / 2) +  currentAdjustedTheta  <<std::endl ; // tilt up 
        tetherCommandStream<<"rotation B Y "<<  -(SimTK::Pi / 2) +  currentAdjustedTheta  <<std::endl ; // tilt up 
        tetherCommandStream<<"rotation A Z "<<  phiFromTheta(currentTheta, pitch, radius, phiOffset,spiralIsRightHanded)   <<std::endl;
        tetherCommandStream<<"rotation B Z "<<  phiFromTheta(currentTheta, pitch, radius, phiOffset,spiralIsRightHanded)   <<std::endl;
        tetherCommandStream<<"readBlockEnd"<<std::endl;


        //MMBLOG_FILE_FUNC_LINE(DEBUG, "MMB-command: tetherToGround A "<<n<<" N1 "<<currentXYZ[0]<<" "<<currentXYZ[1]<<" "<<currentXYZ[2] << " .5 30.0 "<<endl);
        //MMBLOG_FILE_FUNC_LINE(DEBUG, "MMB-command: readBlockEnd "<<endl);
        MMBLOG_FILE_FUNC_LINE(DEBUG, "priorXYZ "<<priorXYZ<<endl);
        priorXYZ= currentXYZ;
        MMBLOG_FILE_FUNC_LINE(DEBUG, "priorXYZ "<<priorXYZ<<endl);
        n++;
    }
    //fclose(spiralPdbFile); 
    monoAtoms.renumberPdbResidues(); // monoAtoms were above added one by one  with no residue numbers set. Now we fix that.
    MMBLOG_FILE_FUNC_LINE(INFO, " closed  spiralPdbFileName = " <<spiralPdbFileName                                   <<endl);
    MMBLOG_FILE_FUNC_LINE(DEBUG, endl);
    std::string commonSpiralCommandsAdjusted = std::regex_replace( commonSpiralCommands, std::regex(std::string("ZZZZ")), std::to_string(n-1 ) ); 
    commonSpiralCommandsAdjusted = std::regex_replace( commonSpiralCommandsAdjusted, std::regex(std::string("            ")), std::string("") );  // Now get rid of extra whitespace
    ofstream spiralCommandsFile2(spiralCommandsFileName); // was "commands.spiral.dat") ;
    spiralCommandsFile2<<commonSpiralCommandsAdjusted<<std::endl;
    spiralCommandsFile2<<tetherCommandStream.str()<<std::endl;
    spiralCommandsFile2.close();
    MMBLOG_FILE_FUNC_LINE(DEBUG, endl);
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
    MMBLOG_FILE_FUNC_LINE(ALWAYS, "sphericalHelix pitch  <distance, in nm>  "<<endl);
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
    MMBLOG_FILE_FUNC_LINE(ALWAYS, "You can create more than one geometry. Specify which one you want: "<<endl);
    MMBLOG_FILE_FUNC_LINE(ALWAYS, "geometry <sphere | cylinder>                                                   "<<endl);
    MMBLOG_FILE_FUNC_LINE(ALWAYS, "And finally, to create the helix, issue:                        "<<endl);
    MMBLOG_FILE_FUNC_LINE(ALWAYS, "sphericalHelix writeCommands "<<endl);
    MMBLOG_FILE_FUNC_LINE(ALWAYS, "The above should be the last sphericalHelix command you issue. "<<endl);
}

    

// This function figures out which polymorphism of parseInput should be used, based on number of parameters provided, be it 0,1,or 3.
void Spiral::parseInput(const map<const String,double> & userVariables, String parameterName, String parameterValue1,String parameterValue2, String parameterValue3,MonoAtomsContainer &monoAtomsContainer){
    MMBLOG_FILE_FUNC_LINE(DEBUG   , "  "<<endl);
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
    MMBLOG_FILE_FUNC_LINE(DEBUG   , "  "<<endl);
    if (commandName == "writeCommands"){
        MMBLOG_FILE_FUNC_LINE(DEBUG   , "  "<<endl);
	if (geometry == geometryEnum::sphere){
            MMBLOG_FILE_FUNC_LINE(DEBUG   , " About to write a sphere. "<<endl);
            writeSphericalSpiralCommandFile (monoAtomsContainer);}
        else if (geometry == geometryEnum::cylinder) {
            MMBLOG_FILE_FUNC_LINE(DEBUG   , " About to write a cylinder. "<<endl);
            writeCylindricalSpiralCommandFile (monoAtomsContainer);}
    }
    else {
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "Non recognizable input: " << commandName << endl);
    }
}
void Spiral::parseInput(const  map<const String,double> & userVariables,String parameterName, String parameterValue){
    MMBLOG_FILE_FUNC_LINE(DEBUG   , "  "<<endl);
    if ((parameterName).compare("radius") ==0)    {
        radius = myAtoF(userVariables,(parameterValue).c_str());
        return;
    } else if ((parameterName).compare("geometry") ==0)    {
        MMBLOG_FILE_FUNC_LINE(DEBUG   , " parameterValue "<<parameterValue<<endl);
	if (parameterValue.compare("sphere")==0){
            MMBLOG_FILE_FUNC_LINE(DEBUG   , "  "<<endl);
            geometry = geometryEnum::sphere;}
        else if (parameterValue.compare("cylinder") == 0) {
            MMBLOG_FILE_FUNC_LINE(DEBUG   , "  "<<endl);
            geometry = geometryEnum::cylinder;}
        return;
    } else if ((parameterName.compare("pitch") ==0) || (parameterName.compare("interStrandDistance") ==0) )    {
        pitch = myAtoF(userVariables,(parameterValue).c_str());
        return;
    } else if ((parameterName).compare("spiralIsRightHanded") ==0)    {
        spiralIsRightHanded = aToBool((parameterValue).c_str());
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
    } else if ((parameterName).compare("cylinderHeight") ==0)    {
        cylinderHeight          = myAtoF(userVariables,(parameterValue).c_str());
        return;
    } else if ((parameterName).compare("spiralPdbFileName")==0){ // was "spiral.pdb"
        spiralPdbFileName = parameterValue;	    
        MMBLOG_FILE_FUNC_LINE(INFO, " Setting spiralPdbFileName = " <<spiralPdbFileName                                   <<endl);
        return;
    } else if ((parameterName).compare("spiralCommandsFileName")==0){ // was "commands.spiral.dat") ;
        spiralCommandsFileName = parameterValue;	    
        return;
    } else if ((parameterName).compare("chainID")==0){ 
        chainID                = parameterValue;	    
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

