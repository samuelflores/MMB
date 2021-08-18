/* -------------------------------------------------------------------------- *
 *                           MMB (MacroMoleculeBuilder)                       *
 * -------------------------------------------------------------------------- *
 *                                                                            *
 * Copyright (c) 2011-12 by the Author.                                       *
 * Author: Samuel Flores                                                      *
 * Modifications by:  Alex Tek                                                           *  
 *                                                                            *
 * See RNABuilder.cpp for the copyright and usage agreement.                  *
 * -------------------------------------------------------------------------- */
//#include "MatlabEngine.hpp"
//#include "MatlabDataArray.hpp"

#define __STDCPP_MATH_SPEC_FUNCS__ 201003L        
#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1
#include <string>
#include "ParameterReader.h"
#include "SimTKmolmodel.h"
#include "BaseInteractionParameterReader.h"
#include "TetherForce.h"
#include "ProgressWriter.h"
#include <iostream>
#include <vector>
#include <cassert>
#include "BasePairContainer.h"
#include "MobilizerContainer.h"
#include "DisplacementContainer.h"
#include "ResidueStretchContainer.h"
#include "MoleculeContainer.h"
#include "NtC_Class_Container.h"
#include "NTC_FORCE_CLASS.h"
#include "NTC_PARAMETER_READER.h"
#include "Utils.h"
//#include "elliptic_integral.h"
// #define _DEBUG_FLAGS_ON_
#include <cerrno>
#include <cmath>
#include <algorithm>
#include <cctype>
//#include <tr1/cmath>
//#include <boost/math/special_functions/ellint_2.hpp>
#include <regex> // for search and replace
#ifdef Lepton_USAGE // This is included only if Lepton_USAGE is defined
#include "Lepton.h"
#endif

#ifdef _WINDOWS
#include <windows.h>

#define PATH_MAX MAX_PATH
#endif // _WINDOWS

using namespace SimTK;
using namespace std  ;


String get_and_set_working_path(String newPath = "RETRIEVE-ONLY" )
{
    MMBLOG_FILE_FUNC_LINE(DEBUG, "PATH_MAX = " << PATH_MAX <<endl);
    char temp [ PATH_MAX ];

    String currentDir = "RETRIEVE-ONLY";
    if (newPath.compare("RETRIEVE-ONLY") == 0) { // This means we do not wish to set the path, only retrieve it.
       if (GetCurrentDir(temp, PATH_MAX) != 0) {
	       currentDir = String (temp);
       }
	//if ( GetCurrentDir(temp, PATH_MAX) != 0) {
	//    currentDir = String (temp); 
        //    currentDir = "howdy2-";
        //    return currentDir; // If GetCurrentDir was successful, return immediately. For some reason Chimera crashes otherwise.
	//}
 
	if (currentDir.compare("/") == 0) {
	    MMBLOG_FILE_FUNC_LINE(INFO, "Current directory is : "<<currentDir<<" . This is unacceptable! We are changing to your home directory."<<endl);
	    char* myHome = getenv("HOME");
	    int chdirError = chdir (myHome);
	    if (chdirError) {
		    MMBLOG_FILE_FUNC_LINE(CRITICAL, " Unable to chdir to "<<myHome<<endl);
	    }
	    if (GetCurrentDir(temp, PATH_MAX) != 0) {
	        currentDir = String (temp); return currentDir;} // If GetCurrentDir was successful, return immediately. For some reason Chimera crashes otherwise.
            currentDir = "howdy1-";
	} else { // from spotlight, this is not called
            MMBLOG_FILE_FUNC_LINE(INFO, "Current directory is : "<<currentDir<<" ."<<endl);
            if (GetCurrentDir(temp, PATH_MAX) != 0) { // Run to check one more time that there is no error.  Otherwise we need to continue on to the error trapping below.
                return currentDir;}
            //currentDir = "howdy";
        }
    } else { // here we are setting the path, and changing into it
	    int chdirError = chdir (newPath.c_str());
	    if (chdirError) {
		    MMBLOG_FILE_FUNC_LINE(CRITICAL, "Unable to chdir to "<<newPath<<endl);
	    }
	    if (GetCurrentDir(temp, PATH_MAX) != 0) {
	    currentDir = String (temp); 
            currentDir = "howdy3-";
            return currentDir;} // If GetCurrentDir was successful, return immediately. For some reason Chimera crashes otherwise.
    } 
    
    int error = errno;
    
    switch ( error ) {
        // EINVAL can't happen - size argument > 0

        // PATH_MAX includes the terminating nul, 
        // so ERANGE should not be returned
        case 0: {return currentDir;} // Everything OK
        case EACCES:
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "Access denied to workingDirectory = "<<currentDir<<endl);
            break;
        case EINVAL:
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "The size argument is zero : "<<PATH_MAX<<endl);
            break;
        case ERANGE:
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "The size argument greater than zero, but smaller than the pathname+1 : "<<PATH_MAX<<endl);
            break;
        case ENOENT:
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "A component of the pathname no longer exists."<<endl);
            break;
        default:
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "Unrecognised error: " <<error<<endl);
            break;
    }
    //return currentDir;    
}

ChainResidueIndex::ChainResidueIndex(int myChainIndex,  int myResidueIndex) {
    residueIndex = myResidueIndex;
    chainIndex = myChainIndex;
};
ParameterReader::ParameterReader() {
    clearContainers();
    MMBLOG_FILE_FUNC_LINE(INFO, ""<<endl);
    myBiopolymerClassContainer.writePdbStructureMapDiagnostics(); //This is just a debug flag
};

void ParameterReader::addC1pSprings (LeontisWesthofClass myLeontisWesthofClass) {
    for (int p=0; p<(int)basePairContainer.numBasePairs(); p++) 
        if  ((myLeontisWesthofClass.myLeontisWesthofBondMatrix.myLeontisWesthofBondRow[basePairContainer.getBasePair(p).leontisWesthofBondRowIndex].distanceC1pC1p >= 0) &&
                applyC1pSprings)
        {

            AtomSpring myAtomSpring;
            atomSpringContainer.initializeAtomSpring(myAtomSpring);
            myAtomSpring.atom1Chain   = basePairContainer.getBasePair(p).FirstBPChain;
            myAtomSpring.atom1Residue = basePairContainer.getBasePair(p).FirstBPResidue;
            myAtomSpring.atom1Name    = "C1'";        
            myAtomSpring.atom2Chain   = basePairContainer.getBasePair(p).SecondBPChain;
            myAtomSpring.atom2Residue = basePairContainer.getBasePair(p).SecondBPResidue;
            myAtomSpring.atom2Name    = "C1'";         
            myAtomSpring.deadLength   = myLeontisWesthofClass.myLeontisWesthofBondMatrix.myLeontisWesthofBondRow[basePairContainer.getBasePair(p).leontisWesthofBondRowIndex].distanceC1pC1p;
            myAtomSpring.forceConstant= .30 * twoTransformForceMultiplier; 
            myAtomSpring.toGround = false  ; 
            myAtomSpring.tether   = false  ; 

            atomSpringContainer.add   (myAtomSpring);

        }
};


void ParameterReader::applyAtomSprings(SimbodyMatterSubsystem & matter, GeneralForceSubsystem & forces, State & state)

{
    MMBLOG_FILE_FUNC_LINE(INFO, "Convert alignmentForces instructions to springs."<<endl);
    atomSpringContainer.createSpringsFromThreading(myBiopolymerClassContainer);
    atomSpringContainer.createSpringsFromGappedThreading(myBiopolymerClassContainer);

    MMBLOG_FILE_FUNC_LINE(INFO, "Applying "<< atomSpringContainer.numAtomSprings()<<" atomSpring's."<<endl);
    for (int i = 0; i < atomSpringContainer.numAtomSprings() ; i++)
    {
        AtomSpring myAtomSpring = atomSpringContainer.getAtomSpring(i);

        MobilizedBody myMobilizedBody1 = matter.Ground();
        MobilizedBody myMobilizedBody2 = matter.Ground();
        Vec3          location1 = Vec3(1000.,1000.,1000.); // initialize to make it obvious when this is not set correctly
        Vec3          location2 = Vec3(1000.,1000.,1000.);
        Vec3          groundLocation1 = Vec3(1000.,1000.,1000.); // initialize to make it obvious when this is not set correctly
        Vec3          groundLocation2 = Vec3(1000.,1000.,1000.);
        // MMBLOG_FILE_FUNC_LINE(" : Applying atom spring : "<<myAtomSpring.atom1Chain<<" "<<myAtomSpring.atom1Residue.outString()<<" " <<myAtomSpring.atom1Name<<" "<<myAtomSpring.atom2Chain<<" "<<myAtomSpring.atom2Residue.outString()<<" " <<myAtomSpring.atom2Name<<endl;
        // get atom 1 location and mobilized body        
        if (myBiopolymerClassContainer.hasChainID(myAtomSpring.atom1Chain )) {
            //MMBLOG_FILE_FUNC_LINE(endl;
            myMobilizedBody1 =  myBiopolymerClassContainer.updAtomMobilizedBody (matter,myAtomSpring.atom1Chain, myAtomSpring.atom1Residue, myAtomSpring.atom1Name);
            location1 = myBiopolymerClassContainer.getAtomLocationInMobilizedBodyFrame (myAtomSpring.atom1Chain, myAtomSpring.atom1Residue, myAtomSpring.atom1Name);
            groundLocation1 = myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,myAtomSpring.atom1Chain,myAtomSpring.atom1Residue,myAtomSpring.atom1Name)  ;
        }
        else if (myMonoAtomsContainer.hasChainID(myAtomSpring.atom1Chain))
        {
            //MMBLOG_FILE_FUNC_LINE(endl;
            myMobilizedBody1 = myMonoAtomsContainer. getMonoAtoms(myAtomSpring.atom1Chain).updMobilizedBody(
                    myAtomSpring.atom1Residue,
                    matter
                    );
            location1 = myMonoAtomsContainer. getMonoAtoms(myAtomSpring.atom1Chain).
                getAtomLocationInMobilizedBodyFrame(myAtomSpring.atom1Residue);
        } else if (waterDropletContainer.hasChainID(myAtomSpring.atom1Chain ))
        {
            waterDropletContainer.validateWaterVectors(); 
            myMobilizedBody1 = waterDropletContainer.updWaterDroplet(myAtomSpring.atom1Chain).updOxygenMobilizedBody(
                    matter,
                    myAtomSpring.atom1Residue
                    );

            waterDropletContainer.validateWaterVectors(); 
            location1 = waterDropletContainer.updWaterDroplet(myAtomSpring.atom1Chain).
                getOxygenLocationInMobilizedBodyFrame(myAtomSpring.atom1Residue);
            if (myAtomSpring.atom1Name.compare("OW") != 0) {    
                MMBLOG_FILE_FUNC_LINE(CRITICAL, "Error! For water, atoms other than OW are not currently supported. You specified "<< myAtomSpring.atom1Name<<endl);
            }
        } else if ( moleculeClassContainer.hasChainID(myAtomSpring.atom1Chain )){
      
            Compound::AtomPathName myAtomPathName = /* myAtomSpring.atom1Residue.outString() + String ("/") + */ myAtomSpring.atom1Name;
            MMBLOG_FILE_FUNC_LINE(INFO, " myAtomPathName = >"<<myAtomPathName<<"< "<<endl);
            Compound::AtomIndex  myAtomIndex = moleculeClassContainer.updMoleculeClass(myAtomSpring.atom1Chain ).molecule. getAtomIndex(myAtomPathName );
            MMBLOG_FILE_FUNC_LINE(INFO, endl);
            MobilizedBodyIndex myAtomMobilizedBodyIndex = moleculeClassContainer.updMoleculeClass(myAtomSpring.atom1Chain ).molecule .getAtomMobilizedBodyIndex( myAtomIndex );
            MMBLOG_FILE_FUNC_LINE(INFO, endl);
            myMobilizedBody1 = matter.updMobilizedBody(myAtomMobilizedBodyIndex);
            MMBLOG_FILE_FUNC_LINE(INFO, endl);
            location1 = moleculeClassContainer.updMoleculeClass(myAtomSpring.atom1Chain ).molecule.getAtomLocationInMobilizedBodyFrame (myAtomIndex);
            MMBLOG_FILE_FUNC_LINE(INFO, "location1 = "<<location1<<endl);
        } else {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "You have tried to attach a spring or tether to atom "<<myAtomSpring.atom1Name<<" of chain "<<myAtomSpring.atom1Chain<<". It does not exist!"<<endl);
        }

        // MMBLOG_FILE_FUNC_LINE(endl;
        //get atom 2 location and mobilizedBody:
        if (myAtomSpring.toGround) {
            myMobilizedBody2 = matter.Ground();
	    
	    if (myAtomSpring.groundLocationIsRelativeToAtom1Location){
	        location2=myAtomSpring.groundLocation += groundLocation1 ; // If the ground location is RELATIVE to the atom, then add the atom's ground location.
	    }
	     
            location2 =  myAtomSpring.groundLocation;
        }  else if (myBiopolymerClassContainer.hasChainID(myAtomSpring.atom2Chain ))
        {
            myMobilizedBody2 =  myBiopolymerClassContainer.updAtomMobilizedBody (matter,myAtomSpring.atom2Chain, myAtomSpring.atom2Residue, myAtomSpring.atom2Name);
            location2 = myBiopolymerClassContainer.getAtomLocationInMobilizedBodyFrame (myAtomSpring.atom2Chain, myAtomSpring.atom2Residue, myAtomSpring.atom2Name);
            groundLocation2 = myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,myAtomSpring.atom2Chain,myAtomSpring.atom2Residue,myAtomSpring.atom2Name)  ;

        }  else if (myMonoAtomsContainer.hasChainID(myAtomSpring.atom2Chain))
        {
            myMobilizedBody2 = myMonoAtomsContainer.getMonoAtoms(myAtomSpring.atom2Chain).updMobilizedBody(
                    myAtomSpring.atom2Residue,
                    matter
                    );
            location2 = myMonoAtomsContainer. getMonoAtoms(myAtomSpring.atom2Chain).
                getAtomLocationInMobilizedBodyFrame(myAtomSpring.atom2Residue);
            //groundLocation2 = myMonoAtomsContainer. getMonoAtoms(myAtomSpring.atom2Chain).
            //    calcAtomLocationInGroundFrame(myAtomSpring.atom2Residue);

        } else if (moleculeClassContainer.hasChainID(myAtomSpring.atom2Chain )){
            Compound::AtomPathName myAtomPathName = /* myAtomSpring.atom2Residue.outString() + String ("/") + */ myAtomSpring.atom2Name;
            Compound::AtomIndex  myAtomIndex = moleculeClassContainer.updMoleculeClass(myAtomSpring.atom2Chain ).molecule. getAtomIndex(myAtomPathName );
            MobilizedBodyIndex myAtomMobilizedBodyIndex = moleculeClassContainer.updMoleculeClass(myAtomSpring.atom2Chain ).molecule .getAtomMobilizedBodyIndex( myAtomIndex );
            myMobilizedBody2 = matter.updMobilizedBody(myAtomMobilizedBodyIndex);
            location2 = moleculeClassContainer.updMoleculeClass(myAtomSpring.atom2Chain ).molecule.getAtomLocationInMobilizedBodyFrame (myAtomIndex);
            //groundLocation2 = moleculeClassContainer.updMoleculeClass(myAtomSpring.atom2Chain ).molecule.calcAtomLocationInMobilizedBodyFrame (myAtomIndex);
            MMBLOG_FILE_FUNC_LINE(INFO, "location2 = "<<location2<<endl);
        } else {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "Error! You have tried to attach a spring or tether to atom "<<myAtomSpring.atom2Name<<" of chain "<<myAtomSpring.atom2Chain<<". It does not exist!"<<endl);
        }
        // SCF: In the below, should determine initial spring length and determine whether to change the dead length to some fraction thereof:
        
        if (  myAtomSpring.tether  )             {
            Force::Custom(forces, new
                    TetherForce(
                        matter,
                        forces,
                        myMobilizedBody1.getMobilizedBodyIndex(),
                        location1,
                        myMobilizedBody2.getMobilizedBodyIndex(),
                        location2,
                        myAtomSpring.forceConstant,
                        myAtomSpring.deadLength
                        ));

        } else {
            double originalLength = (groundLocation2 - groundLocation1).norm(); // not sure .norm() can infer Vec3 type here, confirm.
            MMBLOG_FILE_FUNC_LINE(DEBUG, "groundLocation1 = "<<groundLocation1<<endl);
            MMBLOG_FILE_FUNC_LINE(DEBUG, "groundLocation2 = "<<groundLocation2<<endl);
            if (myAtomSpring.deadLengthIsFractionOfInitialLength) {myAtomSpring.deadLength = myAtomSpring.deadLengthFraction * originalLength;}
            MMBLOG_FILE_FUNC_LINE(DEBUG, "myAtomSpring.deadLengthIsFractionOfInitialLength = "<<myAtomSpring.deadLengthIsFractionOfInitialLength<<endl);
            MMBLOG_FILE_FUNC_LINE(DEBUG, "myAtomSpring.deadLengthFraction = "<<myAtomSpring.deadLengthFraction<<endl);
            MMBLOG_FILE_FUNC_LINE(DEBUG, "originalLength = "<<originalLength<<endl);
            MMBLOG_FILE_FUNC_LINE(DEBUG, "myAtomSpring.deadLength = "<<myAtomSpring.deadLength<<endl);
            Force::TwoPointLinearSpring myTwoPointLinearSpring(
                    forces,
                    myMobilizedBody1,
                    location1,
                    myMobilizedBody2,
                    location2,
                    myAtomSpring.forceConstant,
                    myAtomSpring.deadLength
                    );
            // MMBLOG_FILE_FUNC_LINE(": Adding a two point linear spring with location 1 = "<<location1<<" and location 2 = "<<location2<< endl;
        }
    }
    MMBLOG_FILE_FUNC_LINE(INFO, endl);
    atomSpringContainer.printAtomSprings();
    MMBLOG_FILE_FUNC_LINE(INFO, endl);

}

void ParameterReader::configureDumm(DuMMForceFieldSubsystem & dumm) {
    filebuf fb; 
    fb.open ((tinkerParameterFileName).c_str(),ios::in);
    istream is(&fb);
    if (loadTinkerParameterFile) {
        MMBLOG_FILE_FUNC_LINE(INFO, "You have specified tinkerParameterFileName = "<<tinkerParameterFileName<<".  Checking this file.."<<endl);
        SimTK_ERRCHK_ALWAYS(fb.is_open(),"[Repel.h]", "The Tinker parameter file you specified could not be opened.  Please check your tinkerParameterFileName parameter, or set \"loadTinkerParameterFile 0\" to use the hard-coded Tinker parameters instead.");
        dumm.populateFromTinkerParameterFile (is);
        Biotype::initializePopularBiotypes();
    }
    else 
        dumm.loadAmber99Parameters();
    fb.close();

    dumm.setUseOpenMMAcceleration(useOpenMMAcceleration);
    dumm.setTraceOpenMM(useOpenMMAcceleration);
    dumm.setCoulombGlobalScaleFactor(globalCoulombScaleFactor);
    dumm.setBondTorsionGlobalScaleFactor(globalBondTorsionScaleFactor);
    dumm.setGbsaGlobalScaleFactor(globalGbsaScaleFactor);
    dumm.setVdwGlobalScaleFactor(globalVdwScaleFactor);
    dumm.setBondStretchGlobalScaleFactor(globalBondStretchScaleFactor);
    dumm.setBondBendGlobalScaleFactor(globalBondBendScaleFactor);
    dumm.setAmberImproperTorsionGlobalScaleFactor(globalAmberImproperTorsionScaleFactor);
    dumm.setCustomBondStretchGlobalScaleFactor(0);
    dumm.setCustomBondBendGlobalScaleFactor(0);
    dumm.setUseMultithreadedComputation(useMultithreadedComputation);
    //vector<MagnesiumIon> myMagnesiumIonVec;



}


bool checkForDouble(String const& s) {
    std::istringstream ss(s);
    double d;
    return (ss >> d) && (ss >> std::ws).eof();
};




// a recursive algorithm for reading a double from a String.  This String may contain ints, user variables (begin with @), +, and -.  No whitespaces or additional characters should be in the String.
/*
double   myAtoF(  map<const String,double> myUserVariables,  const char* value){
    MMBLOG_FILE_FUNC_LINE(INFO, "inside myAtoF. converting string : >"<<value<<"<"<<endl);
    
#ifdef Lepton_USAGE
    map<string,double> leptonFormatUserVariables; // Wish this were not necessary. but userVariables uses the signature const SimTK::String,double . Lepton uses string, double.
    leptonFormatUserVariables.clear();
    for (auto  myUserVariablesIterator = myUserVariables.begin() ; myUserVariablesIterator !=myUserVariables.end(); myUserVariablesIterator++) {
        leptonFormatUserVariables[myUserVariablesIterator->first] = myUserVariablesIterator->second;	    
    }
    double leptonResult = Lepton::Parser::parse(std::string(value)).evaluate(leptonFormatUserVariables);
    MMBLOG_FILE_FUNC_LINE(INFO, " Lepton evaluation = >"<< leptonResult<<"< "<<std::endl);
    return leptonResult; // if Lepton_USAGE is defined, then we return here and the rest of the procedure is not used. 
#endif
    // If Lepton_USAGE is NOT defined, then we parse the formula the old dumb way, as follows.

    size_t plusPosition  = String(value).find_last_of('+'); // returns the position of the last '+' in value
    size_t minusPosition = String(value).find_last_of('-'); // ditto for '-'
    if ((plusPosition > minusPosition) && (plusPosition  != String::npos) )  minusPosition = String::npos; // If the plus sign is closer to the end of the string, pretend we didn't find any '-'
    if ((plusPosition < minusPosition) && (minusPosition != String::npos) )  plusPosition  = String::npos; // Conversely, if the '-' is closer to the end, pretend we didn't find any '+' .. actually we might not have found any '+' anyway.
    String baseDoubleString ;
    double          increment = -1111;
    double          decrement = -1111;
    MMBLOG_FILE_FUNC_LINE(INFO, endl);
    if (plusPosition != String::npos) { // We have a '+' to deal with
        MMBLOG_FILE_FUNC_LINE(INFO, endl);
        baseDoubleString = String(value).substr(0, (plusPosition + 0) );
        String incrementString = (String(value).substr(plusPosition+1,1000)); // the second parameter is ridiculously large, but will be truncated at the end of the input String.
        //MMBLOG_FILE_FUNC_LINE(" The increment String is : "<<incrementString<<endl;
        stringstream incrementStringStream(incrementString);
        increment = myAtoF(myUserVariables, incrementString.c_str() );
        decrement = 0;
        MMBLOG_FILE_FUNC_LINE(INFO, endl);
    } else if (minusPosition != String::npos ){ // we have a '-' to deal with
        MMBLOG_FILE_FUNC_LINE(INFO, endl);
        if (minusPosition== 0) {
            // If this is just a leading '-' sign, then put a zero to the left of that minus.
            MMBLOG_FILE_FUNC_LINE(INFO, "Detected a leading \'-\' sign. Will insert a zero to the left of the \'-\'."<<endl);
            baseDoubleString = "0.0";
        }
        else {
            // Otherwise, parse whatever is to the left of the minus sign:
            baseDoubleString = String(value).substr(0, (minusPosition + 0) ); }
        MMBLOG_FILE_FUNC_LINE(INFO, "baseDoubleString =  >"<<baseDoubleString  <<"< "<<endl);

        String decrementString = (String(value).substr(minusPosition+1,1000)); // the second parameter is ridiculously large, but will be truncated at the end of the input String.
        stringstream decrementStringStream(decrementString);
        MMBLOG_FILE_FUNC_LINE(INFO, "About to extract numerical decrement from the string >"<<decrementString<<"< "<<endl);
        decrement = myAtoF(myUserVariables, decrementString.c_str() );
        MMBLOG_FILE_FUNC_LINE(INFO, endl);
        increment = 0;
        //MMBLOG_FILE_FUNC_LINE(endl;
    } else { // no + or - found. This means we can return a result without further recursion.. i.e. we are at a leaf of the recursion tree.
        MMBLOG_FILE_FUNC_LINE(INFO, endl);
        if (!((increment == -1111 ) && (decrement == -1111 )  )) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "Unexplained error!"<<endl);
        }
        baseDoubleString = String(value);        
        increment = 0;
        decrement = 0;
        double baseDouble;
        {
            MMBLOG_FILE_FUNC_LINE(INFO, endl);
            if ((baseDoubleString.substr(0,1)).compare("@") ==0) {
                MMBLOG_FILE_FUNC_LINE(INFO, endl);
                if (myUserVariables.find(baseDoubleString.c_str()) == myUserVariables.end())
                {
                    MMBLOG_FILE_FUNC_LINE(CRITICAL, "Undefined user variable "<<value<<endl);
                }

                MMBLOG_FILE_FUNC_LINE(INFO, "Read user variable "<<baseDoubleString.c_str()<<"  which is set to : "<<myUserVariables[baseDoubleString.c_str()]<<endl);
                baseDouble = double(myUserVariables[baseDoubleString.c_str()]);
            }
            else {
                // This has no '+', '-', or leading '@'.  However it's still possible that the user gave scientific notation, so e.g. 1e-5 leads here to baseDoubleString = '1e'.  So we need to make sure there is nothing but [0-9],'+','-' in this string:
        if (isFixed(baseDoubleString)) {
                    baseDouble = (atof(baseDoubleString.c_str()));
                        MMBLOG_FILE_FUNC_LINE(INFO, "We appear to be in a leaf of the recursion tree for myAtoF. Parsed >"<<baseDoubleString<< "< as : "<<baseDouble<<endl);
        } else {
                        MMBLOG_FILE_FUNC_LINE(CRITICAL, "There was an error processing a putative floating point number."<<endl);
        }
            }
        }
        //MMBLOG_FILE_FUNC_LINE(endl;
        return baseDouble;
    }

    //MMBLOG_FILE_FUNC_LINE(endl;
    double baseDouble = myAtoF(myUserVariables,baseDoubleString.c_str() ) ;

    double finalDouble = baseDouble + increment - decrement;
    MMBLOG_FILE_FUNC_LINE(INFO, "Result of >"<< value  <<"< is : " << finalDouble <<endl);
    return finalDouble;
}*/
/*
bool ParameterReader::aToBool(  const char* value ) {

    String upperValue(value);
    for(int i=0;i<(int)upperValue.length();i++)  {
        upperValue[i] = toupper(value[i]);
    }    

    if (( upperValue ==  "TRUE" ) ||( upperValue ==  "1")) {
        MMBLOG_FILE_FUNC_LINE(DEBUG, "TRUE"<<endl);
        return true;
    }
    else if (( upperValue ==  "FALSE" ) ||( upperValue ==  "0")){
        MMBLOG_FILE_FUNC_LINE(DEBUG, "FALSE"<<endl);
        return false;
    }
    else {
        MMBLOG_FILE_FUNC_LINE(INFO, "Error -- you have specified"<<value<<endl);
        SimTK_ERRCHK_ALWAYS((upperValue == "TRUE" || upperValue == "FALSE" || upperValue == "1"  || upperValue == "0") ,"[ParameterReader.cpp]"," requires either True or False but was set to something else"); 
        return false;
    }    

}    

bool ParameterReader::aToBool( const String& name, const char* value ) {
    return aToBool(value);
}    
bool ParameterReader::compareUpper( const String& param, const char* symbol ) {

    String upperParam(param);
    String upperSym(symbol);

    if( upperParam.length() != upperSym.length() ) return false; 


    for(int i=0;i<(int)upperParam.length();i++)  {
        upperParam[i] = toupper(param[i]);
        upperSym[i] = toupper(symbol[i]);
    }    

    if( upperParam ==  upperSym )
        return true;
    else 
        return false;
}    

*/
// determines whether the provided residue is in any Rigid stretch in the base operation vector.
bool isInRigidStretch(const ResidueID & myResidueID, const String myChainID, const MobilizerContainer & mobilizerContainer    ){
    MMBLOG_FILE_FUNC_LINE(DEBUG, " Testing residueID "<<myResidueID.outString()<<endl);

    auto & residueStretchVector = mobilizerContainer.getResidueStretchVector();
    for (int i = 0 ; i< mobilizerContainer.getNumResidueStretches()           ; i++) {
	    if (residueStretchVector[i].getBondMobility() == SimTK::BondMobility::Rigid){
            MMBLOG_FILE_FUNC_LINE(DEBUG, "  mobilizer start residue = "<< residueStretchVector[i].getStartResidue().outString()<<endl);
            MMBLOG_FILE_FUNC_LINE(DEBUG, "  mobilizer end residue = "<< residueStretchVector[i].getEndResidue().outString()<<endl);
            if (residueStretchVector[i].getChain()        == myChainID){
            if (residueStretchVector[i].getStartResidue() <= myResidueID){
                if (residueStretchVector[i].getEndResidue() >= myResidueID){
                    MMBLOG_FILE_FUNC_LINE(DEBUG, " In a rigid stretch. return "<<1<<endl);
                    return true;}
	    }}
	}
    }
    MMBLOG_FILE_FUNC_LINE(DEBUG, " Not In a rigid stretch. return "<<0<<endl);
    return false; // If you got this far, you are not in any rigid stretch
}
// determines whether the provided residue is DEEP in any Rigid stretch in the base operation vector.
// The idea is that we do not want to remove base pairs that are in a rigid stretch , but bordering a flexible one.
// This way we can later add the helicalStacking and NtC forces across the rigid-flexible boundary. This should improve structural quality.
bool isDeepInRigidStretch(const ResidueID & myResidueID, const String myChainID, const MobilizerContainer & mobilizerContainer    ){
    MMBLOG_FILE_FUNC_LINE(DEBUG, " Testing residueID "<<myResidueID.outString()<<endl);

    auto & residueStretchVector = mobilizerContainer.getResidueStretchVector();
    for (int i = 0 ; i< mobilizerContainer.getNumResidueStretches()           ; i++) {
	    if (residueStretchVector[i].getBondMobility() == SimTK::BondMobility::Rigid){
            MMBLOG_FILE_FUNC_LINE(DEBUG, "  mobilizer start residue = "<< residueStretchVector[i].getStartResidue().outString()<<endl);
            MMBLOG_FILE_FUNC_LINE(DEBUG, "  mobilizer end residue = "<< residueStretchVector[i].getEndResidue().outString()<<endl);
            if (residueStretchVector[i].getChain()        == myChainID){
            if (residueStretchVector[i].getStartResidue() < myResidueID){  // would be <= in isInRigidStretch
                if (residueStretchVector[i].getEndResidue() > myResidueID){// would be >= in isInRigidStretch
                    MMBLOG_FILE_FUNC_LINE(DEBUG, " In a rigid stretch. return "<<1<<endl);
                    return true;}
	    }}
	}
    }
    MMBLOG_FILE_FUNC_LINE(DEBUG, " Not In a rigid stretch. return "<<0<<endl);
    return false; // If you got this far, you are not in any rigid stretch
}
// This function removes density forces from all residues in Rigid stretches.             
// replaced with selectivelyRemoveRigidMobilizerStretchesFromResidueStretchContainer , from BiopolymerClassContainer
/*
void ParameterReader::removeDensityForcesFromRigidStretches () { 
    for (int j = 0 ; j< (int)densityContainer.numDensityStretches() ; j++) {
        if ( isInRigidStretch(basePairContainer.getBasePair(j).FirstBPResidue, basePairContainer.getBasePair(j).FirstBPChain , mobilizerContainer) ) {
            if (isDeepInRigidStretch(basePairContainer.getBasePair(j).SecondBPResidue, basePairContainer.getBasePair(j).SecondBPChain ,mobilizerContainer)){
                MMBLOG_FILE_FUNC_LINE(DEBUG, " Deleting base pair."<< basePairContainer.getBasePair(j).FirstBPResidue.outString()<<" "<<basePairContainer.getBasePair(j).SecondBPResidue.outString() <<endl);
                basePairContainer.deleteBasePair(j);
                j--; 
	    } // of if Second	
	} // of if First 
    } // of for j
    MMBLOG_FILE_FUNC_LINE(INFO    , " Printing all base interactions   ,at the end of removeBasePairsAcrossRigidStretches."<<endl);
    basePairContainer.printBasePairs();
    MMBLOG_FILE_FUNC_LINE(INFO    , " Printing all mobilizer stretches, at the end of removeBasePairsAcrossRigidStretches."<<endl);
    mobilizerContainer.printMobilizerStretches();
}; // of function
*/


// This function removes base pairs when both residues of the base pair are in ANY rigid stretch. That means the two residues could be in DIFFERENT rigid stretches, or in the SAME rigid stretch.
void ParameterReader::removeBasePairsAcrossRigidStretches () { 
    //for (int j = 0 ; j< (int)mobilizerContainer.numMobilizerStretches()          ; j++) {
    for (int j = 0 ; j< (int)basePairContainer.numBasePairs() ; j++) {
        // baseOperationVector holds the endpoints of the rigid segment, while myBasePairVector holds the residues involved in the base pairing interaction.
        if ( isDeepInRigidStretch(basePairContainer.getBasePair(j).FirstBPResidue, basePairContainer.getBasePair(j).FirstBPChain , mobilizerContainer) ) {
            if (isDeepInRigidStretch(basePairContainer.getBasePair(j).SecondBPResidue, basePairContainer.getBasePair(j).SecondBPChain ,mobilizerContainer)){
                MMBLOG_FILE_FUNC_LINE(DEBUG, " Deleting base pair."<< basePairContainer.getBasePair(j).FirstBPResidue.outString()<<" "<<basePairContainer.getBasePair(j).SecondBPResidue.outString() <<endl);
                basePairContainer.deleteBasePair(j);
                j--; 
	    } // of if Second	
	} // of if First 
    } // of for j
    MMBLOG_FILE_FUNC_LINE(INFO    , " Printing all base interactions   ,at the end of removeBasePairsAcrossRigidStretches."<<endl);
    basePairContainer.printBasePairs();
    MMBLOG_FILE_FUNC_LINE(INFO    , " Printing all mobilizer stretches, at the end of removeBasePairsAcrossRigidStretches."<<endl);
    mobilizerContainer.printMobilizerStretches();
}; // of function

// This function only removes base pairs if both residues are in the SAME rigid stretch.
void ParameterReader::removeBasePairsInRigidStretch () { 
    MMBLOG_FILE_FUNC_LINE(CRITICAL, " This function probably does not work anymore."<<endl);
    for (int i = 0; i<(int)baseOperationVector.size(); i++) {
        //scf
        if (((baseOperationVector[i].BasePairIsTwoTransformForce).compare("mobilizer") == 0) && 
                ((baseOperationVector[i].FirstBPEdge).compare("Rigid"    ) == 0)) 
        {

            SimTK_ERRCHK_ALWAYS(
                    (baseOperationVector[i].FirstBPResidue <=baseOperationVector[i].SecondBPResidue), 
                    "[ParameterReader.cpp]","To use removeBasePairsInRigidStretch() you must have the start residue be lower-numbered than the end residue, for all mobilizer commands with the Rigid keyword. " );
            //MMBLOG_FILE_FUNC_LINE(DEBUG, " Will now delete all pairs where both members are in chain "<< baseOperationVector[i].FirstBPChain  <<" and "<<endl);

            for (int j = 0 ; j< (int)basePairContainer.numBasePairs() ; j++) {
                // baseOperationVector holds the endpoints of the rigid segment, while myBasePairVector holds the residues involved in the base pairing interaction.
                if ( 

                        (basePairContainer.getBasePair(j).FirstBPResidue  <= baseOperationVector[i].SecondBPResidue) &&
                        (basePairContainer.getBasePair(j).SecondBPResidue <= baseOperationVector[i].SecondBPResidue) &&
                        (basePairContainer.getBasePair(j).FirstBPResidue  >=  baseOperationVector[i].FirstBPResidue) &&
                        (basePairContainer.getBasePair(j).SecondBPResidue >= baseOperationVector[i].FirstBPResidue)  && 
                        (basePairContainer.getBasePair(j).FirstBPChain.compare(baseOperationVector[i].FirstBPChain) == 0)  && 
                        (basePairContainer.getBasePair(j).SecondBPChain.compare(baseOperationVector[i].SecondBPChain) == 0)  && 
                        (basePairContainer.getBasePair(j).FirstBPChain.compare(basePairContainer.getBasePair(j).SecondBPChain) == 0)   

                   ) {
                    //MMBLOG_FILE_FUNC_LINE(" about to delete BASE PAIR " <<j<<endl; 
                    basePairContainer.deleteBasePair(j);
                    //(basePairContainer.myBasePairVector).erase((basePairContainer.myBasePairVector).begin()+j);                          
                    j--; // need to still examine the next base pair, which is now lower numbered by 1.
                    //if (i > j) i--;
                }
            }
        }
    } 
};

void ParameterReader::printAllSettingsToMMCIF ( std::vector< std::pair < std::string, std::string > > &remarksVec ) {
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "addAllAtomSterics                      bool    " + std::to_string ( addAllAtomSterics ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "addAllHeavyAtomSterics                 bool    " + std::to_string ( addProteinBackboneSterics ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "addProteinBackboneSterics              bool    " + std::to_string ( addRNABackboneSterics ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "addRNABackboneSterics                  bool    " + std::to_string ( addRNABackboneSterics ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "addSelectedAtoms                       bool    " + std::to_string ( addSelectedAtoms ) + " : Add steric spheres to certain RNA atoms as specified in the RNABuilder parameter file. " ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "useCIFFileFormat                       bool    " + std::to_string ( useCIFFileFormat ) + " : Use mmCIF formatted files instead of PDB formatted files for internal and output files. " ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "alignmentForcesIsGapped                bool    " + std::to_string ( alignmentForcesIsGapped ) + " : Determines whether gaps are allowed in the alignment in alignmentForces command. Can vary through the course of the input commands file. This is only the final value." ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "alignmentForcesGapPenalty              double  " + std::to_string ( alignmentForcesGapPenalty ) + " : The penalty applied to gaps. The noGaps condition is enforced with a high value of this parameter. Can vary through the course of the input commands file. This is only the final value." ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "alignmentForcesDeadLength              double  " + std::to_string ( alignmentForcesDeadLengthFraction ) + " : The final length to which the alignmentSprings equilibrate. Defaults to 0.             " ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "alignmentForcesDeadLengthFraction      double  " + std::to_string ( alignmentForcesDeadLengthFraction ) + " : The fraction of the initial length to which the alignmentSprings equilibrate. Should be in the interval (0,1]. A nonzero value enables e.g. progressive morphing.  Can vary through the course of the input commands file. This is only the final value." ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "alignmentForcesForceConstant           double  " + std::to_string ( alignmentForcesForceConstant ) + " : Force constant for the  alignmentForces springs. Can vary through the course of the input commands file. This is only the final value." ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "applyC1pSprings                        bool    " + std::to_string ( applyC1pSprings ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "calcEnergy                             bool    " + std::to_string ( calcEnergy ) ) );
//    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "constrainRigidSegments                 bool    " + std::to_string ( constrainRigidSegments ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "constraintTolerance                    double  " + std::to_string ( constraintTolerance ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "cutoffRadius                           double  " + std::to_string ( cutoffRadius ) + " (nm)" ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "densityAtomFraction                    String  " + std::to_string ( densityAtomFraction ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "densityFileName                        String  " + std::string ( myDensityMap.getDensityFileName()  ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "densityFitPhosphates                   bool    " + std::to_string ( densityFitPhosphates ) + " : When set to False, this means phophate groups in DNA and RNA will feel zero density map fitting force. Be warned that this slows down your run A LOT -- proportional to the number of nucleic acid residues that have fitting forces turned on." ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "densityForceConstant                   double  " + std::to_string ( myDensityMap.getForceConstant( )  ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "densityNoiseComputeAutocorrelation     double  " + std::to_string ( densityNoiseComputeAutocorrelation ) + " : Compute the autocorrelation function for both the planck's law noise and input density. may only have effect if densityNoiseScale > 0." ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "densityReportAtEachAtomPosition        bool    " + std::to_string ( densityReportAtEachAtomPosition ) + " : Write out the local density observed at each atom position, and the corresponding atom name. Written to stdout. Only works when the density forces are active. " ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "densityNoiseTemperature                double  " + std::to_string ( densityNoiseTemperature ) + " : Temperature for the Planck's Law based noise generator for the density map." ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "densityNoiseScale                      double  " + std::to_string ( densityNoiseScale ) + " : Overall scale of the noise for the Planck's Law based noise generator for the density map. Note that this scales the noise amplitude, but the amplitude is squared prior to being added to the density map, and being used to compute signalToNoiseRatio. " ) );
//    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "densityMapActivate                     String  " + std::to_string ( densityMapActivate ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "addAllAtomSterics                      bool    " + std::to_string ( addAllAtomSterics ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "dutyCycle                              double  " + std::to_string ( dutyCycle ) + " : Must lie in (0,1)" ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "electroDensityFileName                 String  " + std::string ( myDensityMap.getDensityFileName() ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "electroDensityForceConstant            double  " + std::to_string ( myDensityMap.getForceConstant( )  ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "excludedVolumeRadius                   double  " + std::to_string ( excludedVolumeRadius ) + " : Radius (in nm) of contact spheres to be applied in AllHeavyAtomSterics and AllAtomSterics." ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "excludedVolumeStiffness                double  " + std::to_string ( excludedVolumeStiffness ) ) );
//    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "firstResidueMobilizerType              String  " + std::to_string ( firstResidueMobilizerType ) + " : use constrainToGround to set this to Weld." ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "firstStage                             int     " + std::to_string ( firstStage ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "fitDefaultTolerance                    double  " + std::to_string ( fitDefaultTolerance ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "baseInteractionScaleFactor             double  " + std::to_string ( twoTransformForceMultiplier ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "globalAmberImproperTorsionScaleFactor  double  " + std::to_string ( globalAmberImproperTorsionScaleFactor ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "globalBondBendScaleFactor              double  " + std::to_string ( globalBondBendScaleFactor ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "globalBondStretchScaleFactor           double  " + std::to_string ( globalBondStretchScaleFactor ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "globalBondTorsionScaleFactor           double  " + std::to_string ( globalBondTorsionScaleFactor ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "globalCoulombScaleFactor               double  " + std::to_string ( globalCoulombScaleFactor ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "globalGbsaScaleFactor                  double  " + std::to_string ( globalGbsaScaleFactor ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "globalVdwScaleFactor                   double  " + std::to_string ( globalVdwScaleFactor ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "guessCoordinates                       bool    " + std::to_string ( guessCoordinates ) + " : If true, invents coordinates for any atoms missing from the input PDB file." ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "inQVectorFileName                      String  " + std::string ( inQVectorFileName ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "initialSeparation                      double  " + std::to_string ( initialSeparation ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "integratorAccuracy                     double  " + std::to_string ( integratorAccuracy ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "integratorStepSize                     int     " + std::to_string ( integratorStepSize ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "integratorType                         String  " + std::string ( integratorType ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "lastStage                              int     " + std::to_string ( lastStage ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "leontisWesthofInFileName               String  " + std::string ( leontisWesthofInFileName ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "loadTinkerParameterFile                bool    " + std::to_string ( loadTinkerParameterFile ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "outQVectorFileName                     String  " + std::string ( outQVectorFileName ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "magnesiumIonChainId                    String  " + std::string ( magnesiumIonChainId ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "magnesiumIonRadius                     String  " + std::to_string ( magnesiumIonRadius ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "matchHydrogenAtomLocations             bool    " + std::to_string ( matchHydrogenAtomLocations ) + " : If false, do not read the hydrogen atom positions from the input pdb file.  Just guess new atom locations.  This is useful if the hydrogens are in bad (e.g. colinear) locations." ) );
//    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "matchProteinCarboxylOxygenLocations    bool    " + std::to_string ( matchProteinCarboxylOxygenLocations ) + " : If false, do not read the carboxyl oxygen atom positions of proteins from the input pdb file.  Just guess new atom locations.  This is useful if the Oxygens are in bad (e.g. non-planar) locations." ) );
//    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "matchingMinimizerTolerance             bool    " + std::to_string ( matchingMinimizerTolerance ) + " : This sets the tolerance for the minimizer used in optimizing the fitting of internal coordinates to the cartesian coordinates in the input structure file.   The default value typically leads to good accuracy, but you may wish to experiment with larger values to save compute time at the cost of accuracy." ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "matchExact                             bool    " + std::to_string ( matchExact ) + " : If True, this matches all bond lengths, angles, and dihedrals to the 2-, 3-, and 4- neighbor atom sets. Locally the match will be nearly perfect, but over a long biopolymer error can accumulate." ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "matchIdealized                         bool    " + std::to_string ( matchIdealized ) + " : If True, the bond lengths and angles will be set to default (idealized) values and the torsion angles will be iteratively adjusted to match the input structure.  Thus the global structure is likely to be good, but small-scale details will differ from those of the input structure.  This is much more expensive than matchExact." ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "matchOptimize                          bool    " + std::to_string ( matchOptimize ) + " : If True, this matches all bond lengths, angles, and dihedrals to the 2-, 3-, and 4- neighbor atom sets. Locally the match will be nearly perfect, but over a long biopolymer error can accumulate." ) );
//    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "matchPerfect                           bool    " + std::to_string ( matchPerfect ) + " : This is a macro.  If True, it sets BOTH matchExact and matchIdealized to True.  This means that all bond lengths, angles, and dihedrals will be matched locally, and then there will be a global refinement of torsion angles to correct for accumulated error.  This costs as much as matchIdealized, but generally gives better results." ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "minimize                               bool    " + std::to_string ( minimize ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "monteCarloTemperature                  int     " + std::to_string ( monteCarloTemperature ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "monteCarloTemperatureIncrement         int     " + std::to_string ( monteCarloTemperatureIncrement ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "nastGlobalBondTorsionScaleFactor       int     " + std::to_string ( nastGlobalBondTorsionScaleFactor ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "noseHooverTime                         double  " + std::to_string ( noseHooverTime ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "numReportingIntervals                  int     " + std::to_string ( numReportingIntervals ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "outMonteCarloFileName                  String  " + std::string ( outMonteCarloFileName ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "outTrajectoryFileName                  String  " + std::string ( outTrajectoryFileName ) ) );
//    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "physicsWhereYouWantIt                  bool    " + std::to_string ( physicsWhereYouWantIt ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "physicsRadius                          double  " + std::to_string ( physicsRadius ) + " : All residues within physicsRadius of \"flexible\" atoms are included in the physics zone. \"flexible\" is defined as belonging to a body of mass < 40." ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "planarityThreshold                     double  " + std::to_string ( planarityThreshold ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "potentialType                          int     " + std::string ( potentialType ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "prioritize                             int     " + std::to_string ( prioritize ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "proteinCapping                         bool    " + std::to_string ( proteinCapping ) + " : When true, adds terminal capping groups to protein chains. " ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "useNACappingHydroxyls                  bool    " + std::to_string ( useNACappingHydroxyls ) + " : When true (default) replaces the 5' phosphorus with an H5T." ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "randomizeInitialVelocities             bool    " + std::to_string ( setInitialVelocities ) + " : When true, adds a stochastic velocity to each body at the beginning of the stage. " ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "readPreviousFrameFile                  bool    " + std::to_string ( readPreviousFrameFile ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "readMagnesiumPositionsFromFile         bool    " + std::to_string ( readMagnesiumPositionsFromFile ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "removeMomentumPeriod                   double  " + std::to_string ( removeMomentumPeriod ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "removeRigidBodyMomentum                bool    " + std::to_string ( removeRigidBodyMomentum ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "reportingInterval                      double  " + std::to_string ( reportingInterval ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "restrainingForceConstant               double  " + std::to_string ( restrainingForceConstant ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "restrainingTorqueConstant              double  " + std::to_string ( restrainingTorqueConstant ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "rigidifyFormedHelices                  int     " + std::to_string ( rigidifyFormedHelices ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "scrubberPeriod                         double  " + std::to_string ( scrubberPeriod ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "safeParameters                         int     " + std::to_string ( safeParameters ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "setChiBondMobility                     int     " + std::to_string ( setChiBondMobility ) ) );
//    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "setDefaultMDParameters                 int     " + std::to_string ( setDefaultMDParameters ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "setForceAndStericScrubber              bool    " + std::to_string ( setForceAndStericScrubber ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "setForceScrubber                       bool    " + std::to_string ( setForceScrubber ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "setHelicalStacking                     bool    " + std::to_string ( setHelicalStacking ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "setRemoveBasePairsInRigidStretch       bool    " + std::to_string ( setRemoveBasePairsInRigidStretch ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "setTemperature                         bool    " + std::to_string ( setTemperature ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "smallGroupInertiaMultiplier            double  " + std::to_string ( smallGroupInertiaMultiplier ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "sphericalHelixCenter                   Vec3    " + std::to_string ( sphericalHelixCenter[0] ) + " ; " + std::to_string ( sphericalHelixCenter[1] ) + " ; " + std::to_string ( sphericalHelixCenter[2] ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "sphericalHelixRadius                   double  " + std::to_string ( sphericalHelixRadius ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "sphericalHelixStartTheta               double  " + std::to_string ( sphericalHelixStartTheta ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "sphericalHelixPhiOffset                double  " + std::to_string ( sphericalHelixPhiOffset ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "sphericalHelixInterStrandDistance      double  " + std::to_string ( sphericalHelixInterStrandDistance ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "stackAllHelicalResidues                bool    " + std::to_string ( stackAllHelicalResidues ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "temperature                            bool    " + std::to_string ( temperature ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "thermostatType                         String  " + std::string ( thermostatType ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "tinkerParameterFileName                String  " + std::string ( tinkerParameterFileName ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "useFixedStepSize                       bool    " + std::to_string ( useFixedStepSize ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "useMultithreadedComputation            bool    " + std::to_string ( useMultithreadedComputation ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "useOpenMMAcceleration                  bool    " + std::to_string ( useOpenMMAcceleration ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "vanderWallSphereRadius                 double  " + std::to_string ( vanderWallSphereRadius ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "velocityRescalingInterval              int     " + std::to_string ( velocityRescalingInterval ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "verbose                                int     " + std::to_string ( verbose ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "vmdOutput                              int     " + std::to_string ( vmdOutput ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "waterDropletMake                       bool    " + std::to_string ( waterDropletMake ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "waterInertiaMultiplier                 double  " + std::to_string ( waterInertiaMultiplier ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "writeCoordinates                       bool    " + std::to_string ( writeCoordinates ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "writeDoublePrecisionTrajectories       bool    " + std::to_string ( writeDoublePrecisionTrajectories ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "writeLastFrameFile                     bool    " + std::to_string ( writeLastFrameFile ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "workingDirectory                       String  " + std::string ( workingDirectory ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "helixBondMobility                      BondMobility::Mobility " + std::to_string ( helixBondMobility ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "loopBondMobility                       BondMobility::Mobility " + std::to_string ( loopBondMobility ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "overallBondMobility                    BondMobility::Mobility " + std::to_string ( overallBondMobility ) ) );
    remarksVec.push_back ( std::pair < std::string, std::string > ( "3", "chiBondMobility                        BondMobility::Mobility " + std::to_string ( chiBondMobility ) ) );
};

void ParameterReader::printAllSettings (ostream & myOstream, String remarkString ) { 
    myOstream << remarkString << "addAllAtomSterics                      bool    "<<addAllAtomSterics         <<endl;
    myOstream << remarkString << "addAllHeavyAtomSterics                 bool    "<<addAllHeavyAtomSterics         <<endl;
    myOstream << remarkString << "addProteinBackboneSterics              bool    "<<addProteinBackboneSterics         <<endl;
    myOstream << remarkString << "addRNABackboneSterics                  bool    "<<addRNABackboneSterics         <<endl;
    myOstream << remarkString << "addSelectedAtoms                       bool    "<<addSelectedAtoms<<" : Add steric spheres to certain RNA atoms as specified in the RNABuilder parameter file. " <<endl;
    myOstream << remarkString << "useCIFFileFormat                       bool    "<<useCIFFileFormat<<" : Use mmCIF formatted files instead of PDB formatted files for internal and output files. " <<endl;
    myOstream << remarkString << "alignmentForcesIsGapped                bool    "<<alignmentForcesIsGapped<<" : Determines whether gaps are allowed in the alignment in alignmentForces command. Can vary through the course of the input commands file. This is only the final value." <<endl;
    myOstream << remarkString << "alignmentForcesGapPenalty              double  "<<alignmentForcesGapPenalty<<" : The penalty applied to gaps. The noGaps condition is enforced with a high value of this parameter. Can vary through the course of the input commands file. This is only the final value." <<endl;
    myOstream << remarkString << "alignmentForcesDeadLength              double  "<<alignmentForcesDeadLength<<" : The equilibrium length to which the springs will equilibrate. In nm." <<endl;
    myOstream << remarkString << "alignmentForcesDeadLengthFraction      double  "<<alignmentForcesDeadLengthFraction<<" : The fraction of the initial length to which the alignmentSprings equilibrate. Should be in the interval (0,1]. A nonzero value enables e.g. progressive morphing.  Can vary through the course of the input commands file. This is only the final value." <<endl;
    myOstream << remarkString << "alignmentForcesForceConstant           double  "<<alignmentForcesForceConstant<<" : Force constant for the  alignmentForces springs. Can vary through the course of the input commands file. This is only the final value." <<endl;
    myOstream << remarkString << "applyC1pSprings                        bool    "<<applyC1pSprings              <<endl;
    myOstream << remarkString << "calcEnergy                             bool    "<<calcEnergy                        <<endl;
    //myOstream << remarkString << "constrainRigidSegments                 bool    "<<constrainRigidSegments            <<endl;
    myOstream << remarkString << "constraintTolerance                    double   "<<constraintTolerance               <<endl;
    myOstream << remarkString << "cutoffRadius                           double   "<<cutoffRadius<<" (nm)"    <<endl;
    myOstream << remarkString << "densityAtomFraction                    String  "<<densityAtomFraction<<endl;
    myOstream << remarkString << "densityFileName                        String  "<<myDensityMap.getDensityFileName()<<endl;
    myOstream << remarkString << "densityFitPhosphates                   bool    "<<densityFitPhosphates<<" : When set to False, this means phophate groups in DNA and RNA will feel zero density map fitting force. Be warned that this slows down your run A LOT -- proportional to the number of nucleic acid residues that have fitting forces turned on."<<endl;
    myOstream << remarkString << "densityForceConstant                   double   "<<myDensityMap.getForceConstant( ) <<endl;
    myOstream << remarkString << "densityNoiseComputeAutocorrelation                double   "<<densityNoiseComputeAutocorrelation<<" Compute the autocorrelation function for both the planck's law noise and input density. may only have effect if densityNoiseScale > 0."<<endl;
    myOstream << remarkString << "densityReportAtEachAtomPosition        bool    "<<densityReportAtEachAtomPosition<<" Write out the local density observed at each atom position, and the corresponding atom name. Written to stdout. Only works when the density forces are active. "<<endl;
    myOstream << remarkString << "densityNoiseTemperature                double   "<<densityNoiseTemperature<<" Temperature for the Planck's Law based noise generator for the density map."<<endl;
    myOstream << remarkString << "densityNoiseScale                      double   "<<densityNoiseScale<<" Overall scale of the noise for the Planck's Law based noise generator for the density map. Note that this scales the noise amplitude, but the amplitude is squared prior to being added to the density map, and being used to compute signalToNoiseRatio. "<< endl;
    //myOstream << remarkString << "densityMapActivate                     String  "<<densityMapActivate<<endl;
    myOstream << remarkString << "dutyCycle                              double   "<<dutyCycle<<" : Must lie in (0,1) "<<endl;
    myOstream << remarkString << "electroDensityFileName                        String  "<<myDensityMap.getDensityFileName()<<endl;
    myOstream << remarkString << "electroDensityForceConstant                   double   "<<myDensityMap.getForceConstant( ) <<endl;
    myOstream << remarkString << "excludedVolumeRadius                   double   "<<excludedVolumeRadius<<" : Radius (in nm) of contact spheres to be applied in AllHeavyAtomSterics and AllAtomSterics.  "        <<endl;
    myOstream << remarkString << "excludedVolumeStiffness                double   "<<excludedVolumeStiffness      <<endl;
    //myOstream << remarkString << "firstResidueMobilizerType              String  "<< firstResidueMobilizerType<<" : use constrainToGround to set this to Weld "<<endl;
    myOstream << remarkString << "firstStage                             int     "<<firstStage    <<endl;
    myOstream << remarkString << "fitDefaultTolerance                    double   "<<fitDefaultTolerance          <<endl;
    myOstream << remarkString << "baseInteractionScaleFactor             double   "<<twoTransformForceMultiplier <<endl;
    myOstream << remarkString << "globalAmberImproperTorsionScaleFactor  double   "<<globalAmberImproperTorsionScaleFactor    <<endl;
    myOstream << remarkString << "globalBondBendScaleFactor              double   "<<globalBondBendScaleFactor    <<endl;
    myOstream << remarkString << "globalBondStretchScaleFactor           double   "<<globalBondStretchScaleFactor <<endl;
    myOstream << remarkString << "globalBondTorsionScaleFactor           double   "<<globalBondTorsionScaleFactor <<endl;
    myOstream << remarkString << "globalCoulombScaleFactor               double   "<<globalCoulombScaleFactor     <<endl;
    myOstream << remarkString << "globalGbsaScaleFactor                  double   "<<globalGbsaScaleFactor        <<endl;
    myOstream << remarkString << "globalVdwScaleFactor                   double   "<<globalVdwScaleFactor         <<endl;
    myOstream << remarkString << "guessCoordinates                       bool     "<<guessCoordinates<<" : If true, invents coordinates for any atoms missing from the input PDB file."       <<endl;
    myOstream << remarkString << "inQVectorFileName                      String  "<<inQVectorFileName            <<endl;
    myOstream << remarkString << "initialSeparation                      double   "<<initialSeparation                         <<endl;
    myOstream << remarkString << "integratorAccuracy                     double   "<<integratorAccuracy                        <<endl;
    myOstream << remarkString << "integratorStepSize                     int     "<<integratorStepSize                        <<endl;
    myOstream << remarkString << "integratorType                         String  "<<integratorType               <<endl;
    myOstream << remarkString << "lastStage                              int     "<<lastStage     <<endl;
    myOstream << remarkString << "leontisWesthofInFileName               String  "<<leontisWesthofInFileName     <<endl;
    myOstream << remarkString << "loadTinkerParameterFile                bool    "<<loadTinkerParameterFile     <<endl;
    //myOstream << remarkString << "loggingSeverity                        string  "; 
	     //MMBLogger::instance().writeLoggingSeverity()
	     //myOstream<<" This sets how verbose your output will be. Options are DEBUG, INFO, WARNING, ALWAYS, CRITICAL. These are in order of decreasing verbosity. "<<endl;
    myOstream << remarkString << "outQVectorFileName                     String  "<<outQVectorFileName     <<endl;
    myOstream << remarkString << "magnesiumIonChainId                    String  "<<magnesiumIonChainId     <<endl;
    myOstream << remarkString << "magnesiumIonRadius                     String  "<<magnesiumIonRadius      <<endl;
    myOstream << remarkString << "matchHydrogenAtomLocations             bool  "<<matchHydrogenAtomLocations    <<" : If false, do not read the hydrogen atom positions from the input pdb file.  Just guess new atom locations.  This is useful if the hydrogens are in bad (e.g. colinear) locations."<<endl;
    //myOstream << remarkString << "matchProteinCarboxylOxygenLocations    bool  "<<matchProteinCarboxylOxygenLocations    <<" : If false, do not read the carboxyl oxygen atom positions of proteins from the input pdb file.  Just guess new atom locations.  This is useful if the Oxygens are in bad (e.g. non-planar) locations."<<endl;
    //myOstream << remarkString << "matchingMinimizerTolerance           bool  "<<matchingMinimizerTolerance<<" This sets the tolerance for the minimizer used in optimizing the fitting of internal coordinates to the cartesian coordinates in the input structure file.   The default value typically leads to good accuracy, but you may wish to experiment with larger values to save compute time at the cost of accuracy."<<endl;
    myOstream << remarkString << "matchExact                             bool  "<<matchExact     <<" If True, this matches all bond lengths, angles, and dihedrals to the 2-, 3-, and 4- neighbor atom sets. Locally the match will be nearly perfect, but over a long biopolymer error can accumulate."<<    endl;
    myOstream << remarkString << "matchIdealized                         bool  "<<matchIdealized <<" If True, the bond lengths and angles will be set to default (idealized) values and the torsion angles will be iteratively adjusted to match the input structure.  Thus the global structure is likely to be good, but small-scale details will differ from those of the input structure.  This is much more expensive than matchExact." <<endl;
    myOstream << remarkString << "matchOptimize                          bool  "<<matchOptimize     <<" If True, this matches all bond lengths, angles, and dihedrals to the 2-, 3-, and 4- neighbor atom sets. Locally the match will be nearly perfect, but over a long biopolymer error can accumulate."<<    endl;
    //myOstream << remarkString << "matchPerfect                         bool  "<<matchPerfect   <<" This is a macro.  If True, it sets BOTH matchExact and matchIdealized to True.  This means that all bond lengths, angles, and dihedrals will be matched locally, and then there will be a global refinement of torsion angles to correct for accumulated error.  This costs as much as matchIdealized, but generally gives better results."<</endl;
    myOstream << remarkString << "minimize                               bool  "<<minimize             <<endl;
    myOstream << remarkString << "monteCarloTemperature                  int   "<<monteCarloTemperature<<endl;
    myOstream << remarkString << "monteCarloTemperatureIncrement         int   "<<monteCarloTemperatureIncrement    <<endl;
    myOstream << remarkString << "nastGlobalBondTorsionScaleFactor       int   "<<nastGlobalBondTorsionScaleFactor     <<endl;
    myOstream << remarkString << "noseHooverTime                         double "<<noseHooverTime          <<endl;
    myOstream << remarkString << "numReportingIntervals                  int     "<<numReportingIntervals     <<endl;
    myOstream << remarkString << "outMonteCarloFileName                  String  "<<outMonteCarloFileName     <<endl;
    myOstream << remarkString << "outTrajectoryFileName                  String  "<<outTrajectoryFileName     <<endl;
    //myOstream << remarkString << "physicsWhereYouWantIt                  bool    "<<physicsWhereYouWantIt         <<endl;
    myOstream << remarkString << "physicsRadius                          double  "<<physicsRadius         <<" : All residues within physicsRadius of \"flexible\" atoms are included in the physics zone. \"flexible\" is defined as belonging to a body of mass < 40."<<endl;
    myOstream << remarkString << "planarityThreshold                     double  "<<planarityThreshold    <<endl;
    myOstream << remarkString << "potentialType                          int     "<<potentialType         <<endl;
    myOstream << remarkString << "prioritize                             int     "<<prioritize     <<endl;
    myOstream << remarkString << "proteinCapping                         bool    "<<proteinCapping     << " : When true, adds terminal capping groups to protein chains. "<<endl;
    myOstream << remarkString << "useNACappingHydroxyls                  bool    "<<useNACappingHydroxyls<<" : When true (default) replaces the 5' phosphorus with an H5T."<<endl;
    myOstream << remarkString << "randomizeInitialVelocities             bool    "<<setInitialVelocities<<" : When true, adds a stochastic velocity to each body at the beginning of the stage. " <<endl;
    myOstream << remarkString << "readPreviousFrameFile                  bool    "<<readPreviousFrameFile     <<endl;
    myOstream << remarkString << "readMagnesiumPositionsFromFile         bool    "<<readMagnesiumPositionsFromFile<<endl;
    myOstream << remarkString << "removeMomentumPeriod                   double   "<<removeMomentumPeriod<<endl;
    myOstream << remarkString << "removeRigidBodyMomentum                bool    "<<removeRigidBodyMomentum<<endl;
    myOstream << remarkString << "reportingInterval                      double   "<<reportingInterval     <<endl;
    myOstream << remarkString << "restrainingForceConstant               double   "<<restrainingForceConstant<<endl;
    myOstream << remarkString << "restrainingTorqueConstant              double   "<<restrainingTorqueConstant<<endl;
    myOstream << remarkString << "rigidifyFormedHelices                  int     "<<rigidifyFormedHelices     <<endl;
    myOstream << remarkString << "scrubberPeriod                         double   "<<scrubberPeriod     <<endl;
    myOstream << remarkString << "safeParameters                         int     "<<safeParameters     <<endl;
    myOstream << remarkString << "setChiBondMobility                     int     "<<setChiBondMobility     <<endl;
    //myOstream << remarkString << "setDefaultMDParameters               int     "<<setDefaultMDParameters     <<endl;
    myOstream << remarkString << "setForceAndStericScrubber              bool    "<<setForceAndStericScrubber     <<endl;
    myOstream << remarkString << "setForceScrubber                       bool    "<<setForceScrubber     <<endl;
    myOstream << remarkString << "setHelicalStacking                     bool    "<<setHelicalStacking     <<endl;
    myOstream << remarkString << "setRemoveBasePairsInRigidStretch       bool    "<<setRemoveBasePairsInRigidStretch     <<endl;
    myOstream << remarkString << "setRemoveBasePairsAcrossRigidStretches bool    "<<setRemoveBasePairsAcrossRigidStretches     <<endl;
    myOstream << remarkString << "setTemperature                         bool    "<<setTemperature     <<endl;
    myOstream << remarkString << "smallGroupInertiaMultiplier            double   "<<smallGroupInertiaMultiplier <<endl;
    myOstream << remarkString << "sphericalHelixCenter                   double   "<<sphericalHelixCenter        <<endl;
    myOstream << remarkString << "sphericalHelixRadius                   double   "<<sphericalHelixRadius        <<endl;
    myOstream << remarkString << "sphericalHelixStartTheta               double   "<<sphericalHelixStartTheta        <<endl;
    myOstream << remarkString << "sphericalHelixPhiOffset                double   "<<sphericalHelixPhiOffset         <<endl;
    myOstream << remarkString << "sphericalHelixInterStrandDistance      double   "<<sphericalHelixInterStrandDistance<<endl;
    myOstream << remarkString << "stackAllHelicalResidues                bool    "<<stackAllHelicalResidues     <<endl;
    myOstream << remarkString << "temperature                            double  "<<temperature                 <<endl;
    myOstream << remarkString << "thermostatType                         String  "<<thermostatType              <<endl;
    myOstream << remarkString << "tinkerParameterFileName                String  "<<tinkerParameterFileName     <<endl;
    myOstream << remarkString << "useFixedStepSize                       bool    "<<useFixedStepSize                <<endl;
    myOstream << remarkString << "useMultithreadedComputation            bool    "<<useMultithreadedComputation     <<endl;
    myOstream << remarkString << "useOpenMMAcceleration                  bool    "<<useOpenMMAcceleration           <<endl;
    myOstream << remarkString << "vanderWallSphereRadius                 double   "<<vanderWallSphereRadius     <<endl;
    myOstream << remarkString << "velocityRescalingInterval              int     "<<velocityRescalingInterval     <<endl;
    myOstream << remarkString << "verbose                                int     "<<verbose     <<endl;
    myOstream << remarkString << "vmdOutput                              int     "<<vmdOutput     <<endl;
    myOstream << remarkString << "waterDropletMake                       bool    "<<waterDropletMake     <<endl;
    myOstream << remarkString << "waterInertiaMultiplier                 double   "<<     waterInertiaMultiplier <<endl;
    myOstream << remarkString << "writeCoordinates                       bool    "<<writeCoordinates     <<endl;
    myOstream << remarkString << "writeDoublePrecisionTrajectories       bool    "<<writeDoublePrecisionTrajectories     <<endl;
    myOstream << remarkString << "writeLastFrameFile                     bool    "<<writeLastFrameFile       <<endl;
    myOstream << remarkString << "workingDirectory                       String  "<<workingDirectory         <<endl;
    myOstream << remarkString << "helixBondMobility                      BondMobility::Mobility"<<helixBondMobility     <<endl;
    myOstream << remarkString << "loopBondMobility                       BondMobility::Mobility"<<loopBondMobility     <<endl;
    myOstream << remarkString << "overallBondMobility                    BondMobility::Mobility"<<overallBondMobility     <<endl;
    myOstream << remarkString << "chiBondMobility                        BondMobility::Mobility"<<chiBondMobility     <<endl;

};


void ParameterReader::removeNonPriorityBasePairs (int priorityLevel) {  
    int oldDutyCyclePriority=0;             
    priority = priorityLevel; //not sure if this is such a great place to set this parameter.

    MMBLOG_FILE_FUNC_LINE(INFO, "At this stage, temperature = "<< temperature  <<endl);



};

void ParameterReader::updateBasePair(int index, 
                                     String ch1, int res1, String edge1, 
                                     String ch2, int res2, String edge2, 
                                     String orient)
{
    basePairContainer.updateBasePair(index,
                                     ch1, res1, edge1,
                                     ch2, res2, edge2,
                                     orient,
                                     myBiopolymerClassContainer,
                                     _leontisWesthofClass,
                                     setHelicalStacking);
}

void ParameterReader::updateMobilizerStretch(int index,
                                             String chainId,
                                             int startRes,
                                             int endRes,
                                             String bondMobility){
    mobilizerContainer.updateMobilizerStretch(index, chainId,
                                              ResidueID(startRes,' '),
                                              ResidueID(endRes,' '),
                                              bondMobility,
                                              myBiopolymerClassContainer);
}


void ParameterReader::addAllResiduesWithin(String chainID, int resID, double radius){
    BiopolymerClass & poly = myBiopolymerClassContainer.updBiopolymerClass(chainID);
    ResidueID res(resID, ' ');
    poly.validateResidueID(res);

    AllResiduesWithin physicsResWith   (chainID, res, radius);
    includeAllResiduesWithinVector.push_back(physicsResWith);
}

void ParameterReader::updateAllResiduesWithin(int index, String chainID, int resID, double radius){
    BiopolymerClass & poly = myBiopolymerClassContainer.updBiopolymerClass(chainID);
    ResidueID res(resID, ' ');
    poly.validateResidueID(res);
    if(index < 0 || index >= includeAllResiduesWithinVector.size()){
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "you tried to update a non existing includeAllResiduesWithin command." << endl);
    }
    includeAllResiduesWithinVector[index].setChain ( chainID);
    includeAllResiduesWithinVector[index].setResidue (res);
    includeAllResiduesWithinVector[index].setRadius ( radius);
}

void ParameterReader::deleteAllResiduesWithin(int index){
    if(index < 0 || index >= includeAllResiduesWithinVector.size()){
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "you tried to delete a non existing includeAllResiduesWithin command." << endl);
    }

    includeAllResiduesWithinVector.erase(includeAllResiduesWithinVector.begin()+index);
}


void ParameterReader::updateIncludeAllNonBondAtomsInResidue(int index,
                                                            const String & chainID,
                                                            int resID)
{
    BiopolymerClass & poly = myBiopolymerClassContainer.updBiopolymerClass(chainID);
    ResidueID res(resID, ' ');
    poly.validateResidueID(res);
    if(index < 0 || index >=  physicsContainer.getNumResidueStretches()){  // includeAllNonBondAtomsInResidueVector.size()){
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "you tried to update a non existing includeAllNonBondAtomsInResidue command." << endl);
    }

    auto & residueStretchVector = physicsContainer.updResidueStretchVector();
    residueStretchVector[index].setChain(chainID); //includeAllNonBondAtomsInResidueVector[index].setChain ( chainID);
    residueStretchVector[index].setResidue(res);   //includeAllNonBondAtomsInResidueVector[index].setStartResidue ( res);
}

void ParameterReader::deleteIncludeAllNonBondAtomsInResidue(int index)
{
    if(index < 0 || index >= physicsContainer.getNumResidueStretches()){ //includeAllNonBondAtomsInResidueVector.size()){
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "you tried to delete a non existing includeAllResiduesWithin command." << endl);
    }
    physicsContainer.deleteResidueStretch(index);
    //includeAllNonBondAtomsInResidueVector.erase(includeAllNonBondAtomsInResidueVector.begin()+index);
}

void ParameterReader::setLeontisWesthofBondRowIndex () 
{ 
    for (int i = 0; i<(int)baseOperationVector.size(); i++) // for restraint and restraintToGround and constraint
        if (!compareUpper("GROUND",(baseOperationVector[i].FirstBPChain).c_str()))
            if (!(baseOperationVector[i].BasePairIsTwoTransformForce.compare("mobilizer") ==0)) 
                if (!(baseOperationVector[i].BasePairIsTwoTransformForce.compare("constraint") ==0)) 
                    if (!(baseOperationVector[i].FirstBPEdge.compare("AllHeavyAtomSterics") ==0)) 
                        if (!(baseOperationVector[i].FirstBPEdge.compare("AllAtomSterics") ==0)) 
                        {
                            baseOperationVector[i].leontisWesthofBondRowIndex = _leontisWesthofClass.getLeontisWesthofBondRowIndex(
                                  myBiopolymerClassContainer.getPdbResidueName(baseOperationVector[i].FirstBPChain,baseOperationVector[i].FirstBPResidue),
                                  myBiopolymerClassContainer.getPdbResidueName(baseOperationVector[i].SecondBPChain,baseOperationVector[i].SecondBPResidue),
                                  //getResidueName(baseOperationVector[i].SecondBPChain,baseOperationVector[i].SecondBPResidue,myParameterReader,myMolecule),
                                  baseOperationVector[i].FirstBPEdge,
                                  baseOperationVector[i].SecondBPEdge,
                                  baseOperationVector[i].OrientationBP,
                                  baseOperationVector[i].BasePairIsTwoTransformForce
                            ) ;
                        }    
} ; 

void ParameterReader::parameterStringInterpreter(const String & paramstr)
{
    ParameterStringClass parameterStringClass( paramstr );
    parameterStringInterpreter(parameterStringClass);
}

void ParameterReader::parameterStringInterpreter(const ParameterStringClass & parameterStringClass,
                                                 const int readStage, 
                                                 const bool readAtOneStageOnly,
                                                 const bool readOnlyUntilStage,
                                                 const bool readExcept)
{
    parameterStringClass.print();	 
    if (   ((parameterStringClass.getString(0)).compare("-SQ") ==0 ) 
            || ((parameterStringClass.getString(0)).compare("sequence") ==0 )  
            || ((parameterStringClass.getString(0)).compare("rnaSequence") ==0 ) 
            || ((parameterStringClass.getString(0)).compare("RNA") ==0 ) 
            || ((parameterStringClass.getString(0)).compare("DNA") ==0 ) 
            || ((parameterStringClass.getString(0)).compare("proteinSequence") ==0 ) 
            || ((parameterStringClass.getString(0)).compare("protein") ==0 ) 
            || ((parameterStringClass.getString(0)).compare("sprnaSequence") ==0 ) 
            || ((parameterStringClass.getString(0)).compare("coarseNucleicAcidSequence") ==0 ) 
       )   
    {
        if (densityContainer.numDensityStretches() > 0) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "You have already declared "<<densityContainer.numDensityStretches()<<" biopolymer stretches to be fitted to the density map.  Please do this after you have created ALL biopolymer chains."<<endl);
        }

        if (myBiopolymerClassContainer.hasChainID(parameterStringClass.getString(1))) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "Can't use chain "<<parameterStringClass.getString(1) <<" to identify this "<<parameterStringClass.getString(0)<<", because it is already being used by a biopolymer in your system."<<endl);
        }
        //MMBLOG_FILE_FUNC_LINE(endl;

        for (int i = 0 ; i < (int)waterDropletAboutResidueVector.size(); i++) {
            if (waterDropletAboutResidueVector[i].waterDropletChainID.compare(parameterStringClass.getString(1) ) == 0) {    
                MMBLOG_FILE_FUNC_LINE(CRITICAL, "Can't use chain "<<parameterStringClass.getString(1)<<" to identify this "<<parameterStringClass.getString(0)<<", because it is already being used by a water droplet."<<endl);
            }
        }
        if (
                (addAllHeavyAtomSterics    ) ||
                (addProteinBackboneSterics ) ||
                (addSelectedAtoms          ) ||
                (addProteinBackboneSterics ) ||
                (addRNABackboneSterics     )
           )
        { 
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "You must specify any general sterics command (e.g. addAllHeavyAtomSterics, addProteinBackboneSterics, addSelectedAtoms, addProteinBackboneSterics, addRNABackboneSterics) AFTER the last biopolymer is specified.  "<<endl);

        }

        if (contactContainer.numContacts() > 0) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "You can only issue a contact command AFTER the last biopolymer is specified.  "<<endl);
        }


        //if this is the sequence parameter
        SimTK_ERRCHK_ALWAYS(
                (baseOperationVector.size()==0)&&(basePairContainer.numBasePairs()==0) , 
                "[ParameterReader.cpp]",": You must specify all chains prior to specifying any contact, mobilizer, or baseInteraction commands, including overall steric commands such as e.g. addSelectedAtoms.");

        SimTK_ERRCHK_ALWAYS(
                (atomSpringContainer.numAtomSprings()==0), 
                "[ParameterReader.cpp]",": You must specify all chains prior to specifying any atomSpring command.");


        if (parameterStringClass.getString(3).length()==0){
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "You have not provided enough parameters when specifying the sequence! You should use, for example :  \"RNA A 16 AUGC \" .. this creates an RNA with chain ID \"A\", first residue number = 16, and sequence AUGC."<<endl);
        }
        if (parameterStringClass.getString(4).length() >0){
            MMBLOG_FILE_FUNC_LINE(CRITICAL, " : You have provided too many parameters when specifying the sequence! You should use, for example :  \"RNA A 16 AUGC \" .. this creates an RNA with chain ID \"A\", first residue number = 16, and sequence AUGC."<<endl);
        }
        if (((parameterStringClass.getString(0)).compare("proteinSequence") ==0 ) || ((parameterStringClass.getString(0)).compare("protein") ==0 ) || ((parameterStringClass.getString(0)).compare("Protein") ==0 )  ) {
            MMBLOG_FILE_FUNC_LINE(INFO, endl);
            myBiopolymerClassContainer.addBiopolymerClass(parameterStringClass.getString(3),parameterStringClass.getString(1), 
                    ResidueID((parameterStringClass.getString(2))),
                    BiopolymerType::Protein, proteinCapping, previousFrameFileName, readPreviousFrameFile,useNACappingHydroxyls);
            //myBiopolymerClassContainer.updBiopolymerClass(parameterStringClass.getString(1)).modifyResidue();
            //MMBLOG_FILE_FUNC_LINE("This is temporary .. SCF"<<endl;
            //myBiopolymerClassContainer.updBiopolymerClass(String("A")).modifyResidue(_dumm); // can't do here, don't have dumm..
        } 
        else if (((parameterStringClass.getString(0)).compare("sprnaSequence") ==0 )) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "Unsupported sequence type."<<endl);
            //sequenceTypes[(parameterStringClass.getString(1)).c_str()] = "sprna";
        } 
        else if (((parameterStringClass.getString(0)).compare("coarseNucleicAcidSequence") ==0 )) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, ": Unsupported sequence type."<<endl);
            //sequenceTypes[(parameterStringClass.getString(1)).c_str()] = "CoarseNucleicAcid";
        } 
        else if ((((parameterStringClass.getString(0)).compare("RNA") ==0 )) ) 
        {
            /*ResidueID myResidue; // This makes a "plain" residue, not attached to an actual chain, and not validated in any way.
            if (parameterStringClass.getString(2).find( "@") != std::string::npos) {  // The user is trying to invoke a user variable.           
                //myAtoI(userVariables,(parameterStringClass.getString(2)).c_str());
                MMBLOG_FILE_FUNC_LINE(" Detected you are trying to invoke a user variable (starts with '@'). This will be interpreted as a string. No insertion code can be specified in this way."<<std::endl;
                myResidue = ResidueID(userVariables,(parameterStringClass.getString(2)).c_str()); 
///// s/ s/ d arguemnt is insertion code, will default to " " if left out.
            } else {
                MMBLOG_FILE_FUNC_LINE(" Detected you are NOT trying to invoke a user variable (starts with '@').  insertion code can be specified in this way."<<std::endl;
                myResidue = ResidueID(parameterStringClass.getString(2)); // I could have just done nothing, but this seems safer.
            }*/
    
            myBiopolymerClassContainer.addBiopolymerClass(parameterStringClass.getString(3),parameterStringClass.getString(1), 
                    //ResidueID( myAtoI(userVariables,(parameterStringClass.getString(2)).c_str())),
                    ResidueID(parameterStringClass.getString(2)) , 
                    //myResidue,
                    BiopolymerType::RNA, false, previousFrameFileName, readPreviousFrameFile, useNACappingHydroxyls);
        }
        else if ((((parameterStringClass.getString(0)).compare("DNA") ==0 )) ) 
        {

            /*ResidueID myResidue; // This makes a "plain" residue, not attached to an actual chain, and not validated in any way.
            if (parameterStringClass.getString(2).find( "@") != std::string::npos) {  // The user is trying to invoke a user variable.           
                //myAtoI(userVariables,(parameterStringClass.getString(2)).c_str());
                MMBLOG_FILE_FUNC_LINE(" Detected you are trying to invoke a user variable (starts with '@'). This will be interpreted as a string. No insertion code can be specified in this way."<<std::endl;
                myResidue = ResidueID(userVariables,(parameterStringClass.getString(2)).c_str()); 
                MMBLOG_FILE_FUNC_LINE(std::endl;
///// s/ s/ d arguemnt is insertion code, will default to " " if left out.
            } else {
                MMBLOG_FILE_FUNC_LINE(" Detected you are NOT trying to invoke a user variable (starts with '@').  insertion code can be specified in this way."<<std::endl;
                myResidue = ResidueID(parameterStringClass.getString(2)); // I could have just done nothing, but this seems safer.
                MMBLOG_FILE_FUNC_LINE(std::endl;
            }    
            MMBLOG_FILE_FUNC_LINE(std::endl; */


            //MMBLOG_FILE_FUNC_LINE(CRITICAL, ": Unsupported sequence type: DNA"<<endl;  
            myBiopolymerClassContainer.addBiopolymerClass(parameterStringClass.getString(3),parameterStringClass.getString(1), 
                    //ResidueID( userVariables,(parameterStringClass.getString(2)).c_str()),
                    ResidueID(parameterStringClass.getString(2)) ,
                    //myResidue,
                    BiopolymerType::DNA, false, previousFrameFileName, readPreviousFrameFile, useNACappingHydroxyls);
            //MMBLOG_FILE_FUNC_LINE(" myResidue = "<<myResidue.outString()<<std::endl;
        }
        else {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "Unsupported sequence type."<<endl);
        }
        myBiopolymerClassContainer.setOriginalSequence(parameterStringClass.getString(1), parameterStringClass.getString(3));
        myBiopolymerClassContainer.updBiopolymerClass(parameterStringClass.getString(1)).renumberPdbResidues(  ResidueID(parameterStringClass.getString(2))  );
        //myBiopolymerClassContainer.updBiopolymerClass(parameterStringClass.getString(1) ).printBiopolymerInfo();
        myBiopolymerClassContainer.updBiopolymerClass(parameterStringClass.getString(1)).loadResidueIDVectorAscending(ResidueID(parameterStringClass.getString(2)) ); // was at the end of postInitialize, now here and when reading loadSequencesFromPdb // this was previously in BiopolymerClass::initializeBiopolymer .. but I found addHelicalStacking was running too slowly in the absence of an initialized residueIDVector

        return;
    }

    if  ((parameterStringClass.getString(0)).compare(  "molecule") == 0)  {
        if ((parameterStringClass.getString(1)).compare(  "initialize") == 0) {
            MMBLOG_FILE_FUNC_LINE(INFO, "You are trying to create a custom molecule.  Start with the command : molecule initialize <chain ID> <residue name>" << endl);
	    parameterStringClass.validateNumFields(4);
	    String myChain = parameterStringClass.getString(2);
	    String myResidueName = parameterStringClass.getString(3);
	    MoleculeClass myMoleculeClass;
	    moleculeClassContainer.add(myChain,myResidueName,myMoleculeClass);
        } else {
            MMBLOG_FILE_FUNC_LINE(INFO, "Next, start issuing commands : molecule  <chain ID> <command> <parameter 1> <parameter 2> <parameter 3> ..." << endl);
            String myChain = parameterStringClass.getString(1);
            moleculeClassContainer.validateChainID(myChain);
            vector <String> commandVector;
            commandVector.clear();
            int i = 2;
            while (parameterStringClass.getString(i).length() > 0){
                //MMBLOG_FILE_FUNC_LINE(": reading in value i = "<<i<<" : " <<parameterStringClass.getString(i)<<endl;
                commandVector.push_back(parameterStringClass.getString(i)); 
                i++;
            }
            MMBLOG_FILE_FUNC_LINE(INFO, "commandVector has length "<<commandVector.size()<<endl);
            moleculeClassContainer.updMoleculeClass(myChain).addOneCommand(commandVector);
        }
        return;
    }
    if  ((parameterStringClass.getString(0)).compare(  "methane") == 0)  {
        parameterStringClass.validateNumFields(3);
        String myChain = parameterStringClass.getString(1);
        String myResidueName = parameterStringClass.getString(2);
        vector <String> myCommand; myCommand.clear();
        vector <vector <String> > myCommandVector; myCommandVector.clear();
        ///stringstream myBiotypeIndexString ;
        String tempArray[] = {"","","","",""};
        tempArray[0] = "setBaseCompound";
        tempArray[1] = "methyl";
    myCommand = vector <String> (tempArray, tempArray+2); 
        myCommandVector.push_back(myCommand);
        
        tempArray[0] = "convertInboardBondCenterToOutboard";
        tempArray[1] = "";
    myCommand = vector <String> (tempArray, tempArray+1); 
        myCommandVector.push_back(myCommand);

    tempArray[0] =        "bondAtom";
        tempArray[1] = "AliphaticHydrogen";
        tempArray[2] = "H4";
        tempArray[3] = String("methyl/bond")  ; 
        tempArray[4] = "0.1112";
        myCommand = vector <String>  (tempArray, tempArray+5);
        myCommandVector.push_back(myCommand);


        //myBiotypeIndexString << Biotype::MethaneC().getIndex();
    tempArray[0] =        "setBiotypeIndex";
        tempArray[1] = "C";
        tempArray[2] = String("MethaneC")  ; 
        myCommand = vector <String>  (tempArray, tempArray+3);
        myCommandVector.push_back(myCommand);
        

        //myBiotypeIndexString << Biotype::MethaneH().getIndex();
    tempArray[1] = "H1"; tempArray[2] = String("MethaneH")  ; 
        myCommand = vector <String>  (tempArray, tempArray+3);
        myCommandVector.push_back(myCommand);

    tempArray[1] = "H2";
        myCommand = vector <String>  (tempArray, tempArray+3);
        myCommandVector.push_back(myCommand);

    tempArray[1] = "H3";
        myCommand = vector <String>  (tempArray, tempArray+3);
        myCommandVector.push_back(myCommand);

    tempArray[1] = "H4";
        myCommand = vector <String>  (tempArray, tempArray+3);
        myCommandVector.push_back(myCommand);

    tempArray[0] = "defineAndSetChargedAtomType";
        tempArray[1] = "MethaneC";
    tempArray[2] = "1"  ; 
        tempArray[3] = "-0.18"  ;  // similar to Alanine CB
        myCommand = vector <String>  (tempArray, tempArray+4);
        myCommandVector.push_back(myCommand);

    tempArray[0] = "defineAndSetChargedAtomType";
        tempArray[1] = "MethaneH";
    tempArray[2] = "34"  ; 
        tempArray[3] = "0.06"  ; 
        myCommand = vector <String>  (tempArray, tempArray+4);
        myCommandVector.push_back(myCommand);


        MoleculeClass myMoleculeClass(myCommandVector);
        moleculeClassContainer.add(myChain,myResidueName,myMoleculeClass);
            //MMBLOG_FILE_FUNC_LINE(CRITICAL, ":  this command is not ready yet."<<endl;
        return;
         
    }
    if  (((parameterStringClass.getString(0)).compare(  "changeSequence") == 0))  
    {
        MMBLOG_FILE_FUNC_LINE(INFO, "This command changes the sequence of the specified chainID. The new sequence must be of same length than the current one. Syntax: changeSequence <chain ID> <new sequence> ." << endl);
        if (parameterStringClass.getString(3).size() != 0) 
        {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "You have specified too many parameters for this command."<<endl);
        }
        String myChain = parameterStringClass.getString(1);
        BiopolymerClass & myBpc = myBiopolymerClassContainer.updBiopolymerClass(myChain);
        myBpc.changeSequence(parameterStringClass.getString(2));
        return;
    }
    if  ((parameterStringClass.getString(0)) == (  "renumberBiopolymerResidues"))  {
        cout<<endl<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": This command renumbers biopolymers. "<<endl;//New residues numbers will be consecutively increasing integers starting with that provided. "<<endl; 
        MMBLOG_FILE_FUNC_LINE(INFO, "Syntax:  "<<endl);
        MMBLOG_FILE_FUNC_LINE(INFO, "renumberBiopolymerResidues                                     : This renumbers all biopolymers to start with residue number 1."<<endl);
        //MMBLOG_FILE_FUNC_LINE(": renumberBiopolymerResidues <first residue number> >          : This renumbers all biopolymers to start with the specified residue number."<<endl;
        //MMBLOG_FILE_FUNC_LINE(": renumberBiopolymerResidues <chain ID> <first residue number> : This renumbers the specified biopolymers to start with the specified residue number."<<endl<<endl;

        if (parameterStringClass.getString(1).length() != 0)  { 
            MMBLOG_FILE_FUNC_LINE(INFO, "  "<<endl);
        }
        if (parameterStringClass.getString(1).length() == 0) //  
        // This is the case that no parameters have been specified. We interpret this to mean that the user wants to renumber all biopolymer chains, to start with residue number 1.
        {
            MMBLOG_FILE_FUNC_LINE(INFO, "You have called renumberBiopolymerResidues with no parameters. We interpret this to mean that you want to renumber all biopolymers to start with 1."<<endl);
            myBiopolymerClassContainer.setRenumberPdbResidues(1);
        } else if (parameterStringClass.getString(2).length() == 0) 
        {
             MMBLOG_FILE_FUNC_LINE(INFO, "  "<<endl);
             MMBLOG_FILE_FUNC_LINE(CRITICAL, "Unexplained error!"<<endl
                                           <<"You have called renumberBiopolymerResidues with one parameter. We interpret this to mean that you want to renumber all biopolymers to start with the residue number : "
                                           << parameterStringClass.getString(1)  <<endl);
             myBiopolymerClassContainer.renumberPdbResidues(ResidueID((parameterStringClass.getString(1) )  )); 
        }      
        else if ((parameterStringClass.getString(3).length() == 0) ) // This is the case that the user has specified only a first
        {
             MMBLOG_FILE_FUNC_LINE(INFO, "  "<<endl);
             if (safeParameters) { 
                 MMBLOG_FILE_FUNC_LINE(CRITICAL, "Unexplained error!"<<endl);
             }
             MMBLOG_FILE_FUNC_LINE(INFO, "  "<<endl);
             myBiopolymerClassContainer.updBiopolymerClass(parameterStringClass.getString(1)).renumberPdbResidues(ResidueID(parameterStringClass.getString(2))); 
             MMBLOG_FILE_FUNC_LINE(INFO, "  "<<endl);
        } else {
            MMBLOG_FILE_FUNC_LINE(INFO, "  "<<endl);
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "Unexplained error!"<<endl);
        }
        return;
    }
        
    if  ((((parameterStringClass.getString(0)).compare(  "setPhiPsiAngles") == 0))  ||
     (((parameterStringClass.getString(0)).compare(  "setBackboneAngles") == 0)))  
    {
        MMBLOG_FILE_FUNC_LINE(ALWAYS,
                "This command sets the default backbone (for proteins: psi, phi, and peptide bond angles) for the specified chain and range of residues, to those of the specified secondary structure.  For secondary structure type \'Alpha\' (meaning alpha helical) that is -60 degrees for both. Note that these are the DEFAULT dihedrals, if there is any structure in the input file for any of the given residues, that structure prevails. "<<endl
                <<"At the moment only these backbone configurations are supported: Alpha, ParallelBeta, and AntiParallelBeta.  These work only on proteins, of course. "<<endl
                <<"Note also that you cannot apply this to the first or last residue on the chain, because the phi and psi are defined using atoms from preceding and succeeding residues. "<<endl
                <<"setBackboneAngles <chain ID> <start residue number> <end residue number> <Alpha | ParallelBeta | AntiParallelBeta>."<<endl);

        if ((parameterStringClass.getString(4).length() == 0) ) { MMBLOG_FILE_FUNC_LINE(CRITICAL, "Not enough parameters for this command."<<endl);}
        if ((parameterStringClass.getString(5).length() != 0) ) { MMBLOG_FILE_FUNC_LINE(CRITICAL, "Too many parameters for this command."<<endl);}

        String chainIDString = parameterStringClass.getString(1);
        ResidueID startResidue = myBiopolymerClassContainer.residueID(userVariables, parameterStringClass.getString(2), chainIDString); 
        ResidueID endResidue   = myBiopolymerClassContainer.residueID(userVariables, parameterStringClass.getString(3), chainIDString );
        MMBLOG_FILE_FUNC_LINE(INFO, "Start residue identified as : "<<startResidue.outString()<<endl);
        MMBLOG_FILE_FUNC_LINE(INFO, "End residue identified as : "<<endResidue.outString()<<endl);
        SecondaryStructureStretch mySecondaryStructureStretch;
        mySecondaryStructureStretch.setStartResidue(startResidue);
        mySecondaryStructureStretch.setEndResidue(endResidue);
        mySecondaryStructureStretch.setChain(chainIDString);
        mySecondaryStructureStretch.setSecondaryStructureType(parameterStringClass.getString(4));
        myBiopolymerClassContainer.secondaryStructureStretchVector.push_back( mySecondaryStructureStretch);
        //#myBiopolymerClassContainer.updBiopolymerClass(chainIDString).setAlphaHelicalDefaultBackboneAngles( startResidue, endResidue);
        return;
    }

    if  (((parameterStringClass.getString(0)).compare(  "substituteResidue") == 0))  {
        MMBLOG_FILE_FUNC_LINE(INFO, "This command generates a specified substitution in a specified chain, at a specified residue position. Syntax: substituteResidue <chain ID> <residue number> <new residue type> ."<<endl);
        if (parameterStringClass.getString(4).size() != 0) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "You have specified too many parameters for this command."<<endl);
        }
    
        String myChain = parameterStringClass.getString(1);
        ResidueID myResidue = myBiopolymerClassContainer.residueID(userVariables,parameterStringClass.getString(2),myChain);
        String mySubstitution = parameterStringClass.getString(3);
        Mutation myMutation(myChain,myResidue,mySubstitution);
        myBiopolymerClassContainer.substituteResidue(myMutation,safeParameters, matchPurineN1AtomLocations,proteinCapping);     
        return;
    }
    if  (((parameterStringClass.getString(0)).compare(  "insertResidue") == 0))  {
        MMBLOG_FILE_FUNC_LINE(INFO, "This command generates a specified insertion in a specified chain."<<endl);
        MMBLOG_FILE_FUNC_LINE(INFO, "Syntax: insertResidue <chain ID> <residue number> <new residue type> ."<<endl);
        parameterStringClass.validateNumFields(4);
        if (parameterStringClass.getString(4).size() != 0) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "You have specified too many parameters for this command."<<endl);
        }
    
        String myChain = parameterStringClass.getString(1);
        
        ResidueID myResidue; // This makes a "plain" residue, not attached to an actual chain, and not validated in any way.
        if (parameterStringClass.getString(2).find( "@") != std::string::npos) {  // The user is trying to invoke a user variable.           
            //myAtoI(userVariables,(parameterStringClass.getString(2)).c_str());
            MMBLOG_FILE_FUNC_LINE(INFO, "Detected you are trying to invoke a user variable (starts with '@'). This will be interpreted as a string. No insertion code can be specified in this way."<<endl);
            myResidue = ResidueID(userVariables,(parameterStringClass.getString(2)).c_str()); 
// second arguemnt is insertion code, will default to " " if left out.
        } else {
            MMBLOG_FILE_FUNC_LINE(INFO, "Detected you are NOT trying to invoke a user variable (starts with '@').  insertion code can be specified in this way."<<endl);
            myResidue = ResidueID(parameterStringClass.getString(2)); // I could have just done nothing, but this seems safer.
        }    
        String myResidueType = parameterStringClass.getString(3); 
        Mutation myMutation(myChain,myResidue,myResidueType);
        myBiopolymerClassContainer.insertResidue(myMutation,proteinCapping);        
        return;
    }

    if  (((parameterStringClass.getString(0)).compare(  "deleteResidue") == 0)  ||
             ((parameterStringClass.getString(0)).compare(  "deleteResidues") == 0))  {

        MMBLOG_FILE_FUNC_LINE(ALWAYS, "This command deletes the specified residue in the specified chain. The remaining residues retain the PDB residue numbers they had prior to the deletion (i.e. there may result a gap in the numbering)."<<endl);
        MMBLOG_FILE_FUNC_LINE(ALWAYS, "Syntax: deleteResidues <chain ID> <residue ID>  "<<endl);
        MMBLOG_FILE_FUNC_LINE(ALWAYS, "Syntax: deleteResidues <chain ID> <start residue ID> <end residue ID>  "<<endl);
        MMBLOG_FILE_FUNC_LINE(ALWAYS, "Or, you could just delete the entire chain, like this: "<<endl);
        MMBLOG_FILE_FUNC_LINE(ALWAYS, "Syntax: deleteResidues <chain ID>   "<<endl);
        MMBLOG_FILE_FUNC_LINE(ALWAYS, "Or, you could delete ALL chains, like this: "<<endl);
        MMBLOG_FILE_FUNC_LINE(ALWAYS, "Syntax: deleteResidues    "<<endl);
        MMBLOG_FILE_FUNC_LINE(ALWAYS, "But be careful with this command! Make sure you don't try to reference any residues that this command deletes, elsewhere in your command file."<<endl);
    
        String myChain = parameterStringClass.getString(1);
        if (safeParameters) if (!mobilizerContainer.isEmpty()) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "You have previously called some sort of mobilizer command. Please move all such commands, which reference residue numbers, below the deleteResidues command. Otherwise there  too much potential to reference a deleted residue. You can override this message by setting safeParameters False"<<endl);
        }
        if (parameterStringClass.getString(1).size() == 0) {
            MMBLOG_FILE_FUNC_LINE(INFO, "It appears you wish to delete ALL chains "<<endl);
            myBiopolymerClassContainer.deleteAllBiopolymerClasses(); 
            //MMBLOG_FILE_FUNC_LINE(CRITICAL, ": You have not specified enough parameters for this command."<<endl;
            return;
        }
        else if ((parameterStringClass.getString(1).size() != 0) && (parameterStringClass.getString(2).size() == 0)) { // only a chain ID was provided.  User wants to delete entire chain.
            MMBLOG_FILE_FUNC_LINE(INFO, "It would appears you wish to delete chain "<<myChain<<endl);
            myBiopolymerClassContainer.deleteBiopolymerClass(myChain);
            if (myBiopolymerClassContainer.hasChainID(myChain)) {
                MMBLOG_FILE_FUNC_LINE(CRITICAL, "For some reason chain "<<myChain<<" still exists!"<<endl);
            }
            else MMBLOG_FILE_FUNC_LINE(INFO, "Done! the chain "<<myChain<<" is gone!"<<endl);
            return;
        }
        else if (parameterStringClass.getString(4).size() != 0) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "You have specified too many parameters for this command."<<endl);
        }

        ResidueID startResidue = myBiopolymerClassContainer.residueID(userVariables,parameterStringClass.getString(2),myChain);
        ResidueID endResidue (-1111,' '); // Set to an absurd value to catch errors
        if (parameterStringClass.getString(3).size() != 0) {
            MMBLOG_FILE_FUNC_LINE(INFO, "You have specified a range of residues.."<<endl);
            endResidue = myBiopolymerClassContainer.residueID(userVariables,parameterStringClass.getString(3),myChain);

        } else {
            endResidue = myBiopolymerClassContainer.residueID(userVariables,parameterStringClass.getString(2),myChain); }
        ResidueID deletedResidue = startResidue;
        ResidueID nextResidue (-1111,' ');

        if ((startResidue == myBiopolymerClassContainer.updBiopolymerClass(myChain).getFirstResidueID()) &&
            (endResidue == myBiopolymerClassContainer.updBiopolymerClass(myChain).getLastResidueID())) {
                MMBLOG_FILE_FUNC_LINE(CRITICAL, "You can't delete ALL residues in the chain this way! Just delete the entire chain. See syntax above."<<endl);
        }

        while (deletedResidue<=endResidue) {
            if (myBiopolymerClassContainer.updBiopolymerClass(myChain).getLastResidueID() == myBiopolymerClassContainer.updBiopolymerClass(myChain).getFirstResidueID()) {
                MMBLOG_FILE_FUNC_LINE(CRITICAL, "You can't delete ALL residues in the chain!  not with this command, anyway."<<endl);
            }
            MMBLOG_FILE_FUNC_LINE(INFO, "About to delete "<<deletedResidue.outString()<<endl);
            nextResidue = myBiopolymerClassContainer.updBiopolymerClass(myChain).safeSum(deletedResidue,1); // Will only increment if it's deletedResidue is not the last in the chain.
        String mySubstitution = "?";
        Mutation myMutation(myChain,deletedResidue,mySubstitution);
        myBiopolymerClassContainer.deleteResidue(myMutation,proteinCapping);
                if (deletedResidue == nextResidue) {break;} // this occurs when we reach the end of the chain.
                deletedResidue = nextResidue;
                //if (deletedResidue == myBiopolymerClassContainer.updBiopolymerClass(myChain).getLastResidueID()) {break;}

            }
            return;
        }
        if  (((parameterStringClass.getString(0)).compare("renameChain") == 0))  {

            MMBLOG_FILE_FUNC_LINE(ALWAYS,
                    "Syntax:"<<endl
                    <<"renameChain <old chain ID> <new chain ID>"<<endl
                    <<"Works only with biopolymers for now.     "<<endl);
            parameterStringClass.validateNumFields(3);
            String oldChainID = parameterStringClass.getString(1);    
            String newChainID = parameterStringClass.getString(2);    
            if (!(myBiopolymerClassContainer.hasChainID(oldChainID))) {
                MMBLOG_FILE_FUNC_LINE(CRITICAL, "There is no biopolymer with chain >"<< oldChainID  <<"< !"<<endl);
            }    
            if (myBiopolymerClassContainer.hasChainID(newChainID)) {
                MMBLOG_FILE_FUNC_LINE(CRITICAL, "There is already a biopolymer with chain >"<< newChainID  <<"< !"<<endl);
            }    
            myBiopolymerClassContainer.renameChain(oldChainID,newChainID);
            //BiopolymerClass newBiopolymerClass = myBiopolymerClassContainer.updBiopolymerClass(oldChainID);
            //newBiopolymerClass .renameChain(newChainID);
            //myBiopolymerClassContainer.   addBiopolymerClass(newChainID,newBiopolymerClass);
            //myBiopolymerClassContainer.deleteBiopolymerClass(oldChainID);
             
            if ((myBiopolymerClassContainer.hasChainID(oldChainID))) {
                MMBLOG_FILE_FUNC_LINE(CRITICAL, "There is still a biopolymer with old chain >"<< oldChainID  <<"< !"<<endl);
            }    
            if (!(myBiopolymerClassContainer.hasChainID(newChainID))) {
                MMBLOG_FILE_FUNC_LINE(CRITICAL, "There is no biopolymer with chain >"<< newChainID  <<"< !"<<endl);
            }    
            return;
        }    

        if  (((parameterStringClass.getString(0)).compare(  "monoAtoms") == 0))  {

            MMBLOG_FILE_FUNC_LINE(ALWAYS,
                    "Syntax:"<<endl
                    <<"monoAtoms <chain ID> <first residue number>   <number of ions> <atom name>"<<endl);

            MonoAtoms myMonoAtoms(parameterStringClass.getString(1),
                    ResidueID(parameterStringClass.getString(2)),
                    myAtoI(userVariables,(parameterStringClass.getString(3)).c_str()),
                    parameterStringClass.getString(4));
            myMonoAtomsContainer.addMonoAtoms(myMonoAtoms);




            return;
        }


    if  (((parameterStringClass.getString(0)).compare(  "monoAtoms") == 0))  {

        MMBLOG_FILE_FUNC_LINE(ALWAYS,
                "Syntax:"<<endl
                <<"monoAtoms <chain ID> <first residue number>   <number of ions> <atom name>"<<endl);

        MonoAtoms myMonoAtoms(parameterStringClass.getString(1), 
                ResidueID(parameterStringClass.getString(2)),
                myAtoI(userVariables,(parameterStringClass.getString(3)).c_str()),
                parameterStringClass.getString(4)); 
        myMonoAtomsContainer.addMonoAtoms(myMonoAtoms); 
        return;
    }

    if(parameterStringClass.getString(0).compare("deleteChain") == 0)
    {
        MMBLOG_FILE_FUNC_LINE(ALWAYS,
                "This command allows to remove a chain from MMB. The chain will not be present in the output PDB files."<<endl
                <<"Syntax: deleteChain <chainID>"<<endl);
        parameterStringClass.validateNumFields(2);
        myBiopolymerClassContainer.deleteBiopolymerClass(parameterStringClass.getString(1));
        return;
    }
    if ((parameterStringClass.getString(0)).compare("atomSpring")      == 0) 
    {
        MMBLOG_FILE_FUNC_LINE(ALWAYS, "You have called the atomSpring command. This applies springs or tethers , either between two atoms, or between one atom and ground. There is a 'deadLength' parameter which specifies the zero-force distance for springs, or the radius of the zero-force sphere for tethers. there is a toGround parameter, which specifies that there is no atom 2, only atom 1 and a ground location.There is a 'groundLocationIsRelative' parameter, which specifies that -- in the case of a spring or tether to ground -- the coordinates given are relative to the first (and only) atom's location in the structure. "<<endl);
        MMBLOG_FILE_FUNC_LINE(ALWAYS, "Usage: ");
        MMBLOG_FILE_FUNC_LINE(ALWAYS, "The command file is read top to bottom. So first consider any parameters you wish to change, and set those. Then you will create the spring. Then change any paramters, create again, and so on.  : "<<endl);
        MMBLOG_FILE_FUNC_LINE(ALWAYS, " Parameters include: String parameters atom1Chain, atom1Residue, atom1Name, (ditto for atom2); boolean parameters tether, toGround forceConstant deadLengthIsFractionOfInitialLength, groundLocationIsRelativeToAtom1Location; double-precision parameters deadLengthFraction, deadLength, Vec3 parameter groundLocation. So one line might read e.g. atomSpring forceConstant 130000 .. The final line to create the spring would read atomSpring add  .. Issue the command or parameter of interest to get more context sensitive feedback. "<<endl);
        //MMBLOG_FILE_FUNC_LINE(ALWAYS, "atomSpring  [atom2Chain | atom2Residue | atom2Name ] : Just like the above, except for atom 2. These are ignored if toGround is true."<<endl);
        //MMBLOG_FILE_FUNC_LINE(ALWAYS, endl);
        //AtomSpring dummyAtomSpring is the default temporary AtomSpring, should have been declared in ParameterReader.h. No initialization necessary, it has a default constructor. Becomes and adult AtomSpring once it is added to atomSpringContainer.
        if ((parameterStringClass.getString(1)).compare("atom1Chain")==0){            
            MMBLOG_FILE_FUNC_LINE(ALWAYS, "atomSpring  atom1Chain [chain ID, string] : Chain ID of atom 1."<<endl);
            parameterStringClass.validateNumFields(3);   // expecting atomSpring, the parameter name, and parameter value.
            dummyAtomSpring.atom1Chain = (parameterStringClass.getString(2));
        }
	else if ((parameterStringClass.getString(1)).compare("atom2Chain")==0){            
            MMBLOG_FILE_FUNC_LINE(ALWAYS, "atomSpring  atom1Chain [chain ID, string] : Chain ID of atom 2. Ignored if toGround is true."<<endl);
            parameterStringClass.validateNumFields(3);   // expecting atomSpring, the parameter name, and parameter value.
            dummyAtomSpring.atom2Chain = (parameterStringClass.getString(2));
        }
        else if ((parameterStringClass.getString(1)).compare("atom1Name" )==0){            
            MMBLOG_FILE_FUNC_LINE(ALWAYS, "atomSpring  atom1Name  [name, string] : Name of atom 1."<<endl);
            parameterStringClass.validateNumFields(3);   // expecting atomSpring, the parameter name, and parameter value.
            dummyAtomSpring.atom1Name  = (parameterStringClass.getString(2));
        }
        else if ((parameterStringClass.getString(1)).compare("atom2Name" )==0){            
            MMBLOG_FILE_FUNC_LINE(ALWAYS, "atomSpring  atom1Name  [name, string] : Name of atom 2. Ignored if toGround is true."<<endl);
            parameterStringClass.validateNumFields(3);   // expecting atomSpring, the parameter name, and parameter value.
            dummyAtomSpring.atom2Name  = (parameterStringClass.getString(2));
        }
        else if ((parameterStringClass.getString(1)).compare("atom1Residue" )==0){            
            MMBLOG_FILE_FUNC_LINE(ALWAYS, "atomSpring  atom1Residue [residue ID, string] : Residue ID of atom 1. If there is an insertion code, put it right after the residue number, with no spaces, like this: 123C. You must first specify atom1Chain before specifying the atom1Residue "<<endl);
            MMBLOG_FILE_FUNC_LINE(DEBUG, " chain is :" <<dummyAtomSpring.atom1Chain<<endl);
            parameterStringClass.validateNumFields(3);   // expecting atomSpring, the parameter name, and parameter value.
            dummyAtomSpring.atom1Residue = myBiopolymerClassContainer.residueID(userVariables,parameterStringClass.getString(2), dummyAtomSpring.atom1Chain  );
            //myBiopolymerClassContainer.residueID(userVariables,parameterStringClass.getString(4), parameterStringClass.getString(2)  ),
            //dummyAtomSpring.atom1Residue = ResidueID (parameterStringClass.getString(2));
        }
        else if ((parameterStringClass.getString(1)).compare("atom2Residue" )==0){            
            MMBLOG_FILE_FUNC_LINE(ALWAYS, "atomSpring  atom2Residue [residue ID, string] : Residue ID of atom 2. If there is an insertion code, put it right after the residue number, with no spaces, like this: 123C. You must first specify atom1Chain before specifying the atom1Residue. Ignored if toGround is true. "<<endl);
            parameterStringClass.validateNumFields(3);   // expecting atomSpring, the parameter name, and parameter value.
            dummyAtomSpring.atom2Residue = myBiopolymerClassContainer.residueID(userVariables,parameterStringClass.getString(2), dummyAtomSpring.atom2Chain  );
            //dummyAtomSpring.atom2Residue = ResidueID (parameterStringClass.getString(2));
        }
        else if ((parameterStringClass.getString(1)).compare("tether")==0){            
            MMBLOG_FILE_FUNC_LINE(ALWAYS, "atomSpring tether [true | false] : Defaults to false. If false, behaves as a spring. If true, it is a tether, meaning the potential is flat bottomed, and becomes non-flat at deadLength."<<endl);
            parameterStringClass.validateNumFields(3);   // expecting atomSpring, the parameter name, and parameter value.
            dummyAtomSpring.tether = aToBool(parameterStringClass.getString(2));
        }
        else if ((parameterStringClass.getString(1)).compare("toGround")==0){         
            parameterStringClass.validateNumFields(3);   // expecting atomSpring, the parameter name, and parameter value.      
            dummyAtomSpring.toGround = aToBool(parameterStringClass.getString(2));
        }
        else if ((parameterStringClass.getString(1)).compare("groundLocationIsRelativeToAtom1Location")==0){         
            MMBLOG_FILE_FUNC_LINE(ALWAYS, "atomSpring groundLocationIsRelativeToAtom1Location [true | false] : If true, the groundLocation is relative to atom 1 locatoin. So if groundLocation is 0 0 1, one end of the spring will be fixed 1nm above atom 1's initial location. The other end will be fixed to atom1.  "<<endl);
            parameterStringClass.validateNumFields(3);   // expecting atomSpring, the parameter name, and parameter value.      
            dummyAtomSpring.groundLocationIsRelativeToAtom1Location = aToBool(parameterStringClass.getString(2));
        }
        else if  ((parameterStringClass.getString(1)).compare("deadLengthIsFractionOfInitialLength")==0){               
            MMBLOG_FILE_FUNC_LINE(ALWAYS, "atomSpring deadLengthIsFractionOfInitialLength [true | false] : If true, the dead length of the springs is <deadLengthFraction> * (initial length). Original applicatoin was progressive morphing.  "<<endl);
            parameterStringClass.validateNumFields(3);   // expecting atomSpring, the parameter name, and parameter value.
            dummyAtomSpring.deadLengthIsFractionOfInitialLength = aToBool(parameterStringClass.getString(2));
        }
        else if  ((parameterStringClass.getString(1)).compare("forceConstant")==0){        
            MMBLOG_FILE_FUNC_LINE(ALWAYS, "atomSpring forceConstant [spring constant] : The spring constant, in kJ/mol/nm/nm. Applies to both spring and tether. For reference, the spring constant of a carbon-carbon single bond is 129790.8 kJ/mol/nm/n"<<endl);
            parameterStringClass.validateNumFields(3);   // expecting atomSpring, the parameter name, and parameter value.       
            dummyAtomSpring.forceConstant = myAtoF(userVariables,parameterStringClass.getString(2));
        }
        else if  ((parameterStringClass.getString(1)).compare("deadLengthFraction")==0){               
            MMBLOG_FILE_FUNC_LINE(ALWAYS, "atomSpring deadLengthFraction <fraction> : This sets the dead length of the springs to <fraction> * (initial length). Original applicatoin was progressive morphing.  "<<endl);
            parameterStringClass.validateNumFields(3);   // expecting atomSpring, the parameter name, and parameter value.
            dummyAtomSpring.deadLengthFraction = myAtoF(userVariables,parameterStringClass.getString(2));
        }
        else if  ((parameterStringClass.getString(1)).compare("deadLength")==0){            
            MMBLOG_FILE_FUNC_LINE(ALWAYS, "atomSpring deadLength  [length] : Defaults to 0. In the case of a spring, it is the length at which force is zero. In the case of a tether, the MAXIMUM length at which force is zero. Ignored if deadLengthIsFractionOfInitialLength is true."<<endl);
            parameterStringClass.validateNumFields(3);   // expecting atomSpring, the parameter name, and parameter value.   
            dummyAtomSpring.deadLength = myAtoF(userVariables,parameterStringClass.getString(2));
        }

        else if  ((parameterStringClass.getString(1)).compare("groundLocation")==0){            
            MMBLOG_FILE_FUNC_LINE(ALWAYS, "atomSpring groundLocation  [X, Y, Z] : Location on ground to which one end of the spring will be fixed. The other end of the spring will be fixed to atom 1. Applies if toGround is true. Ignored if toGround is false."<<endl);
            parameterStringClass.validateNumFields(5);   // expecting atomSpring, the parameter name, and parameter value.   
            dummyAtomSpring.groundLocation = Vec3(
                  myAtoF(userVariables,parameterStringClass.getString(2)),
                  myAtoF(userVariables,parameterStringClass.getString(3)),
                  myAtoF(userVariables,parameterStringClass.getString(4))
            );
        }   
        else if  ((parameterStringClass.getString(1)).compare("add")==0){            
            MMBLOG_FILE_FUNC_LINE(ALWAYS, "atomSpring add  :  Once you have set the parameters, somewhere above this command in the file, issue this to create a spring with those parameters. Technically, you are adding an AtomSpring to atomSpringContainer.  You can then change one or more parameters, and issue the command again to create another spring with the new parameters. And so on. "<<endl);
            parameterStringClass.validateNumFields(2);   // expecting atomSpring, the command.   
	    dummyAtomSpring.print();
            atomSpringContainer.add(dummyAtomSpring);
        } else {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "Error! You have specified a parameter or command that is not related to atomSpring.  "<<endl);
        }  
        return;
    }  // of if atomSpring
       
    if ((parameterStringClass.getString(0)).compare("springToGround")      == 0) 
         {
                 MMBLOG_FILE_FUNC_LINE(CRITICAL, " This command is obsolete. Please use atomSpring, with toGround set to True.  "<<endl);
    }
     
    if ((parameterStringClass.getString(0)).compare("tetherToGround")      == 0) 
         {
                 MMBLOG_FILE_FUNC_LINE(CRITICAL, " This command is obsolete. Please use atomSpring, with toGround and tether set to True.  "<<endl);
    }
    if ((parameterStringClass.getString(0)).compare("atomTether")      == 0) 
         {
                 MMBLOG_FILE_FUNC_LINE(CRITICAL, " This command is obsolete. Please use atomSpring, with  tether set to True.  "<<endl);
    }
     
    // comment out the old atomTether,springToGround, etc case
     
     
     
    if ( ((parameterStringClass.getString(0)).compare("nucleicAcidDuplex")      == 0) ) 
    {

        MMBLOG_FILE_FUNC_LINE(ALWAYS, "format: nucleicAcidDuplex <Chain ID A> <start residue A> <end residue A> <Chain B ID> <start residue B> <end residue B> " <<  endl);
        //MMBLOG_FILE_FUNC_LINE(" : format:  nucleicAcidDuplex <Chain ID A> <start residue A> <end residue A> <Chain B ID> <start residue B> <end residue B> "<<std::cout;
        MMBLOG_FILE_FUNC_LINE(ALWAYS, "where start residue A will be Watson-Crick base paired with  end residue B, etc.  we require: (start residue A) <= (end residue A) and (start residue B) >= (end residue B), because duplexes are antiparallel."<<endl);
        String      chainA = parameterStringClass.getString(1);
        String      chainB = parameterStringClass.getString(4);
        ResidueID startResidueA = myBiopolymerClassContainer.residueID(userVariables,parameterStringClass.getString(2),chainA );
        ResidueID   endResidueA = myBiopolymerClassContainer.residueID(userVariables,parameterStringClass.getString(3),chainA );
        ResidueID startResidueB = myBiopolymerClassContainer.residueID(userVariables,parameterStringClass.getString(5),chainB );
        ResidueID   endResidueB = myBiopolymerClassContainer.residueID(userVariables,parameterStringClass.getString(6),chainB );
        int threadedLength = myBiopolymerClassContainer.updBiopolymerClass(chainA).difference(endResidueA , startResidueA) + 1;
        if (parameterStringClass.getString(7).length() >0) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL,
                    "You have specified too many parameters for this command.  Correct format is: "<<endl<<"nucleicAcidDuplex <Chain ID A> <start residue A> <end residue A> <Chain B ID> <start residue B> <end residue B> "<<endl
                    <<"Note that start Residue A <= end Residue A and start Residue B >= end Residue B, since duplexes are antiparallel."<<endl);
        }
        if (!(parameterStringClass.getString(6).length() >0)) {

            MMBLOG_FILE_FUNC_LINE(CRITICAL,
                    "You have not specified enough parameters for this command.  Correct format is: "<<endl<<"nucleicAcidDuplex <Chain ID A> <start residue A> <end residue A> <Chain B ID> <start residue B> <end residue B> "<<endl
                    <<"Note that start Residue A <= end Residue A and start Residue B >= end Residue B, since duplexes are antiparallel."<<endl);

        }
        if(startResidueA > endResidueA)
        {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "In the nucleicAcidDuplex command, the end residue must be greater than or equal to the start residue for each chain."<<endl);
        }
        if(startResidueB < endResidueB)
        {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "In the nucleicAcidDuplex command, the start residue must be greater than or equal to the end residue for each chain."<<endl);
        }
        if(myBiopolymerClassContainer.updBiopolymerClass(chainA).difference( endResidueA,startResidueA) != myBiopolymerClassContainer.updBiopolymerClass(chainB).difference (startResidueB, endResidueB))
        {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "In the nucleicAcidDuplex command, the two paired segments must be of the same length."<<endl);
        }

        // We use a temporary container to discard all the basepairs if one is not validated
        BasePairContainer temp;
        for (int i = 0; i < threadedLength; i++) {
            BaseInteraction myBasePair;
            myBasePair.rotationCorrection1 = Rotation(0.0,UnitVec3(0,0,1));
            myBasePair.rotationCorrection2 = Rotation(0.0,UnitVec3(0,0,1));

            myBasePair.FirstBPChain    = chainA;            
            myBasePair.SecondBPChain   = chainB;            
            myBasePair.FirstBPResidue  =myBiopolymerClassContainer.updBiopolymerClass(chainA).sum(  startResidueA ,i); 

            myBasePair.SecondBPResidue = myBiopolymerClassContainer.updBiopolymerClass(chainB).sum( startResidueB,( - i)); 

            myBasePair.FirstBPEdge     = "WatsonCrick";      
            myBasePair.SecondBPEdge    = "WatsonCrick";      

            myBasePair.OrientationBP = "Cis";                 
            //myBasePair.BasePairPriority = 1;
            //baseOperationVector.push_back(myBasePair); 
            temp.addBasePair(myBiopolymerClassContainer, _leontisWesthofClass, myBasePair, setHelicalStacking); 

        }
        // Copy from temp to the real container
        for(int i=0; i<temp.numBasePairs(); ++i)
        {
            basePairContainer.addBasePair(myBiopolymerClassContainer, _leontisWesthofClass, temp.getBasePair(i), setHelicalStacking);
        }
        return;
    } //End nucleicAcidDuplex


     if ((parameterStringClass.getString(0)).compare("alignmentForces")      == 0) 
     {
         MMBLOG_FILE_FUNC_LINE(ALWAYS, "You have called the alignmentForces command. This applies springs pulling together two user-defined stretches of biopolymer residues. The springs go between corresponding atoms in corresponding residues. You may provide aligned stretches of residues, or let MMB figure out the gapped alignment for you with Seqan. It uses a simple scoring function, which differentiates between \"match\" (residues the same), \"mismatch\" (residues not the same), gap opening, and gap extension."<<endl);
         MMBLOG_FILE_FUNC_LINE(ALWAYS, "Usage: ");
         MMBLOG_FILE_FUNC_LINE(ALWAYS, "First consider any parameters you wish to change : "<<endl);
         MMBLOG_FILE_FUNC_LINE(ALWAYS, "alignmentForces noGap : This sets the internal parameter alignmentForcesIsGapped to False, introduces a very high penalty for gaps in the alignment, and requires that the aligned fragments have the same number of residues. If you were looking for the old \"threading\" command, set this parameter."<<endl);
         MMBLOG_FILE_FUNC_LINE(ALWAYS, "alignmentForces gapped : This means that gappedAlignment is left at True, and Seqan is used to figure out the alignment. This is the default behavior. "<<endl);
         MMBLOG_FILE_FUNC_LINE(ALWAYS, "alignmentForces deadLengthFraction <fraction> : This sets the dead length of the springs to <fraction> * (initial length). Use this e.g. to do progressive morphing.  "<<endl);
         MMBLOG_FILE_FUNC_LINE(ALWAYS,"                                                For <fraction> in the interval (0, 1], this sets the dead length of the springs to <fraction> * (initial length). It also sets alignmentForcesDeadLengthIsFractionOfInitialLength = True. Set this to zero to set alignmentForcesDeadLengthIsFractionOfInitialLength = False and recover ordinary behavior. "<<endl);
         MMBLOG_FILE_FUNC_LINE(ALWAYS, "alignmentForces forceConstant <force constant (double)> : This sets the force constant to the specified value. Unless you are changing it now, it will be : "<< alignmentForcesForceConstant<<" ."<<endl);

         MMBLOG_FILE_FUNC_LINE(ALWAYS, endl);
         MMBLOG_FILE_FUNC_LINE(ALWAYS, "As you may know, MMB command files are read from top to bottom, so the last parameter setting prior to the execution command prevails. Any parameters set an the execution command do not apply to that execution command. Now issue the execution command: "<<endl);

         MMBLOG_FILE_FUNC_LINE(ALWAYS, "Syntax:  <alignmentForces> <Chain A>  <Chain-B> "<<endl);
         MMBLOG_FILE_FUNC_LINE(ALWAYS, "Or    :  <alignmentForces> <Chain A>  <start residue A> <end residue A>  <Chain-B>  <start residue B> <end residue B>   "<<endl);
         MMBLOG_FILE_FUNC_LINE(ALWAYS, endl);


         ThreadingStruct thread;
                   
         //#thread.gapPenalty = alignmentForcesGapPenalty;
         //#thread.isGapped   = alignmentForcesIsGapped;
         MMBLOG_FILE_FUNC_LINE(DEBUG, "Currently alignmentForcesIsGapped = "<<alignmentForcesIsGapped<<endl);
         MMBLOG_FILE_FUNC_LINE(DEBUG, "Currently alignmentForcesGapPenalty = "<<alignmentForcesGapPenalty<<endl);
         MMBLOG_FILE_FUNC_LINE(DEBUG, "Currently thread.isGapped = "<<thread.isGapped<<endl);

       
        
	 double myForceConstant = 1000000000.0; // set to super high value just to make sure it's being reset later.
         //double defaultForceConstant = 30.;

         if ((parameterStringClass.getString(2)).length()==0) // User intends to set a boolean  parameter
         {
             if ((parameterStringClass.getString(1)).compare("noGap")==0){
                    alignmentForcesIsGapped = false; 
                    alignmentForcesGapPenalty = -10000. ; // artificially high value
                    MMBLOG_FILE_FUNC_LINE(INFO, "You have set alignmentForcesIsGapped = "<<alignmentForcesIsGapped<<endl);
                    MMBLOG_FILE_FUNC_LINE(INFO, "You have set alignmentForcesGapPenalty = "<<alignmentForcesGapPenalty<<endl);

                    return; // done with this command, go on to next line in command file
              } else if ((parameterStringClass.getString(1)).compare("gapped")==0){
                  alignmentForcesIsGapped = true ; 
                  alignmentForcesGapPenalty = -1 ; // return to default
                  MMBLOG_FILE_FUNC_LINE(INFO, "You have set alignmentForcesIsGapped = "<<alignmentForcesIsGapped<<endl);
                  MMBLOG_FILE_FUNC_LINE(INFO, "You have set alignmentForcesGapPenalty = "<<alignmentForcesGapPenalty<<endl);
                  return; // done with this command, go on to next line in command file
              } else {
                  MMBLOG_FILE_FUNC_LINE(CRITICAL, "Syntax error!"<<endl);
             }

         // do deadLengthFraction
         } else if ((parameterStringClass.getString(3)).length()==0)  // User intends to set a parameter with a value
         {
              if ((parameterStringClass.getString(1)).compare("deadLengthFraction")==0){
                                       alignmentForcesDeadLengthFraction = myAtoF(userVariables,parameterStringClass.getString(2).c_str());
                         if ((alignmentForcesDeadLengthFraction < 0.0 ) || (alignmentForcesDeadLengthFraction > 1.0 )){
                           if (safeParameters) {
                  MMBLOG_FILE_FUNC_LINE(CRITICAL, "alignmentForcesDeadLengthFraction should lie in the interval (0, 1] if you want to apply progressive forces, or be set to 0 if you want to restore normal behavior. "<<endl);
                              }
                              alignmentForcesDeadLengthIsFractionOfInitialLength = true;
                              alignmentForcesDeadLengthFraction = myAtoF(userVariables,parameterStringClass.getString(2).c_str());
                              MMBLOG_FILE_FUNC_LINE(INFO, "You have set alignmentForcesDeadLengthIsFractionOfInitialLength = "<<alignmentForcesDeadLengthIsFractionOfInitialLength<<endl);
                              MMBLOG_FILE_FUNC_LINE(INFO, "You have set alignmentForcesDeadLengthFraction = " << alignmentForcesDeadLengthFraction <<endl);
                             return;
                                             
                         } else if (alignmentForcesDeadLengthFraction < 1E-14 ){
                              alignmentForcesDeadLengthIsFractionOfInitialLength = false;
                              MMBLOG_FILE_FUNC_LINE(INFO, "You have set alignmentForcesDeadLengthIsFractionOfInitialLength = "<<alignmentForcesDeadLengthIsFractionOfInitialLength<<endl);
                             return;
                         } else {
                              alignmentForcesDeadLengthIsFractionOfInitialLength = true;
                              alignmentForcesDeadLengthFraction = myAtoF(userVariables,parameterStringClass.getString(2).c_str());
                              MMBLOG_FILE_FUNC_LINE(INFO, "You have set alignmentForcesDeadLengthIsFractionOfInitialLength = "<<alignmentForcesDeadLengthIsFractionOfInitialLength<<endl);
             		 MMBLOG_FILE_FUNC_LINE(INFO, "You have set alignmentForcesDeadLengthFraction = " << alignmentForcesDeadLengthFraction <<endl);
 
                             return;
                          }
             } 
             else if ((parameterStringClass.getString(1)).compare("deadLength")==0){
                 alignmentForcesDeadLength         = myAtoF(userVariables,parameterStringClass.getString(2).c_str());
                 if (alignmentForcesDeadLength    <= 0.0) {
                     MMBLOG_FILE_FUNC_LINE(CRITICAL, "deadLength    must be greater than zero! You have specified : "<< alignmentForcesDeadLength  <<endl);
                 }
                 return;
             }
             else if ((parameterStringClass.getString(1)).compare("forceConstant")==0){
                 alignmentForcesForceConstant      = myAtoF(userVariables,parameterStringClass.getString(2).c_str());
                 if (alignmentForcesForceConstant <= 0.0) {
                     MMBLOG_FILE_FUNC_LINE(CRITICAL, "forceConstant must be greater than zero! You have specified : "<<alignmentForcesForceConstant<<endl);
                 }
                 return;
             }
             else {
		 MMBLOG_FILE_FUNC_LINE(CRITICAL, "alignmentForces : parameter "<< parameterStringClass.getString(1)<<" with value "<<parameterStringClass.getString(2)<<" not recognized."<<endl);
             }
	} // of parameter setting section
        else if ( myBiopolymerClassContainer.hasChainID( parameterStringClass.getString(1)) &&   myBiopolymerClassContainer.hasChainID(parameterStringClass.getString(2))) { 
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "This syntax is no longer supported. Please issue alignmentForces forceConstant [double]"<<endl);}
	 /*
	// Syntax:  <alignmentForces> <Chain A>  <Chain-B>
        else if ( myBiopolymerClassContainer.hasChainID( parameterStringClass.getString(1)) &&   myBiopolymerClassContainer.hasChainID(parameterStringClass.getString(2))) { 
            MMBLOG_FILE_FUNC_LINE(INFO, "Detected you wish to align chains "<< parameterStringClass.getString(1)<<" and "<<parameterStringClass.getString(2)<<endl);
	    String      chainA = parameterStringClass.getString(1);
	    String      chainB = parameterStringClass.getString(2);
            //MMBLOG_FILE_FUNC_LINE( " Detected you wish to align chains "<<chainA <<" and "<<chainB <<endl;

	    if ((parameterStringClass.getString(3)).length() >0) 
	    { 
		MMBLOG_FILE_FUNC_LINE(CRITICAL, "Too many parameters for this command! This command no longer takes an optional <force constant> parameter, if that is what you were trying to provide. Instead set the force constant using syntax:  alignmentForces forceConstant <force constant (double)> "<<endl);
         
		myForceConstant= myAtoF(userVariables,parameterStringClass.getString(3).c_str());
	    }
	    else
	    {
		myForceConstant= alignmentForcesForceConstant;
	    }
	    if ((parameterStringClass.getString(4)).length() >0) 
	    {
		MMBLOG_FILE_FUNC_LINE(CRITICAL, "Too many parameters for this command!"<<endl);
	    }
	     
	    MMBLOG_FILE_FUNC_LINE(INFO, "READER ThreadForceConstant "<<myForceConstant<<endl);
    
        MMBLOG_FILE_FUNC_LINE(INFO, "Creating a gapped threading for chain >"<<chainA<<"< versus chain >"<<chainB<<"< with force constant : "<<myForceConstant<<endl);
            thread =  atomSpringContainer.createGappedThreading(chainA, chainB , myForceConstant, false,  myBiopolymerClassContainer);
            //thread.gapPenalty = alignmentForcesGapPenalty;
            thread.isGapped   = alignmentForcesIsGapped;
            thread.deadLengthIsFractionOfInitialLength = alignmentForcesDeadLengthIsFractionOfInitialLength;
            thread.deadLengthFraction = alignmentForcesDeadLengthFraction;
            MMBLOG_FILE_FUNC_LINE(DEBUG, " Currently thread.isGapped = "<<thread.isGapped<<endl);
            if (safeParameters) if (thread.isGapped == false) {
                MMBLOG_FILE_FUNC_LINE(CRITICAL, "If you set \"noGap\" then this command requires residue numbers. Otherwise use a \"gapped\" alignment. "<<endl);
            }



        } // of alignmentFroces ChainA ChainB forceConstant
        */

        else if ((parameterStringClass.getString(4)).length()==0 ) { // Syntax:  <alignmentForces> <Chain A>  <Chain-B>  [forceConstant]

	    MMBLOG_FILE_FUNC_LINE(CRITICAL, "Wrong number of parameters for this command! This command no longer takes an optional <force constant> parameter, if that is what you were trying to provide. Instead set the force constant using syntax:  alignmentForces forceConstant <force constant (double)> "<<endl);


         } else if (((parameterStringClass.getString(7)).length()==0 ) &&  ((parameterStringClass.getString(6)).length()!=0 ))   { 
	 // In case the user u    sed the syntax  <alignmentForces> <Chain A>  <start residue A> <end residue A>  <Chain-B>  <start residue B> <end residue B> 
            MMBLOG_FILE_FUNC_LINE(ALWAYS, " You are invoking the syntax: "<<endl<<" alignmentForces <chain A> <start residue A> <end residue A> <chain B> <start residue B> <end residue B>       "<<endl);

	    if ((parameterStringClass.getString(7)).length() >0) // Should never happen. Just being ultra paranoid.
	    {
		MMBLOG_FILE_FUNC_LINE(CRITICAL, "Wrong number of parameters for this command! This command no longer takes an optional <force constant> parameter, if that is what you were trying to provide. Instead set the force constant using syntax:  alignmentForces forceConstant <force constant (double)> "<<endl);
		//myForceConstant= myAtoF(userVariables,parameterStringClass.getString(7).c_str());
                //MMBLOG_FILE_FUNC_LINE(" myForceConstant= "<<myForceConstant<<endl;
	    }
            MMBLOG_FILE_FUNC_LINE(INFO, "alignmentForcesForceConstant = "<<alignmentForcesForceConstant<<endl);
            thread.updThreadingPartner(0).biopolymerClass =  myBiopolymerClassContainer.updBiopolymerClass( parameterStringClass.getString(1));
            thread.updThreadingPartner(0).startResidue =  myBiopolymerClassContainer.residueID(userVariables, parameterStringClass.getString(2),thread.updThreadingPartner(0).biopolymerClass.getChainID());
//residueStart1   = myBiopolymerClassContainer.residueID(userVariables, parameterStringClass.getString(2),thread.chainID1);
            thread.updThreadingPartner(0).endResidue     = myBiopolymerClassContainer.residueID(userVariables, parameterStringClass.getString(3),thread.updThreadingPartner(0).biopolymerClass.getChainID());
            thread.forceConstant   = alignmentForcesForceConstant; // set to super high value just to make sure it's being reset later.
            thread.updThreadingPartner(1).biopolymerClass =  myBiopolymerClassContainer.updBiopolymerClass( parameterStringClass.getString(4));
            //thread.chainID2        = parameterStringClass.getString(4);
            
            thread.updThreadingPartner(1).startResidue   = myBiopolymerClassContainer.residueID(userVariables, parameterStringClass.getString(5),thread.updThreadingPartner(1).biopolymerClass.getChainID());
            thread.updThreadingPartner(1).  endResidue     = myBiopolymerClassContainer.residueID(userVariables, parameterStringClass.getString(6),thread.updThreadingPartner(1).biopolymerClass. getChainID());
            thread.isGapped   = alignmentForcesIsGapped;
            thread.deadLengthIsFractionOfInitialLength = alignmentForcesDeadLengthIsFractionOfInitialLength;
            thread.deadLengthFraction = alignmentForcesDeadLengthFraction;
            thread.deadLength         = alignmentForcesDeadLength        ;

         } else if ((parameterStringClass.getString(7)).length()>0 )   { // In case the user u    sed the syntax  <alignmentForces> <Chain A>  <start residue A> <end residue A>  <Chain-B>  <start residue B> <end residue B>   [forceConstant]
		MMBLOG_FILE_FUNC_LINE(CRITICAL, "Wrong number of parameters for this command! This command no longer takes an optional <force constant> parameter, if that is what you were trying to provide. Instead set the force constant using syntax:  alignmentForces forceConstant <force constant (double)> "<<endl);
        } else {
	    for (int i = 0 ; i < 10; i++){
	        MMBLOG_FILE_FUNC_LINE(INFO, " parameterStringClass.getString("<<i<<") = "<<parameterStringClass.getString(i)<<endl);
	    }	
	    MMBLOG_FILE_FUNC_LINE(CRITICAL, "Please check your syntax!"<<endl);
        }
        atomSpringContainer.addGappedThreading(thread, myBiopolymerClassContainer);
        return;
    } //End alignmentForces 

    if ( ((parameterStringClass.getString(0)).compare("proteinThreading")      == 0) || 
            ((parameterStringClass.getString(0)).compare("threading")      == 0) ||
            ((parameterStringClass.getString(0)).compare("gappedThreading")      == 0) ||
            ((parameterStringClass.getString(0)).compare("proteinBackboneThreading")      == 0) ) 
    {
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "Obsolete command! Please use alignmentForces. Issue with no parameters to get syntax information."<<endl);
        return;
    }   //End Threading

    if ( ((parameterStringClass.getString(0)).compare("RNAThreading")      == 0) || 
            ((parameterStringClass.getString(0)).compare("alignRNA")      == 0) ) {

        MMBLOG_FILE_FUNC_LINE(ALWAYS, "syntax:      RNAThreading ChainA start-residue-A end-residue-A Chain-B start-residue-B end-residue-B"<<endl);
        MMBLOG_FILE_FUNC_LINE(ALWAYS, "where start-residue-A will be aligned with start-residue-B, etc.  we require: end-residue > start-residue. "<<endl);
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "This command is being phased out! Try the \"threading\" command."<<endl);
        String      chainA = parameterStringClass.getString(1);
        String      chainB = parameterStringClass.getString(4);

        ResidueID startResidueA = myBiopolymerClassContainer.residueID(userVariables,parameterStringClass.getString(2),chainA );
        ResidueID   endResidueA = myBiopolymerClassContainer.residueID(userVariables,parameterStringClass.getString(3),chainA );

        ResidueID startResidueB = myBiopolymerClassContainer.residueID(userVariables,parameterStringClass.getString(5),chainB );
        ResidueID   endResidueB = myBiopolymerClassContainer.residueID(userVariables,parameterStringClass.getString(6),chainB );
        int threadedLength = myBiopolymerClassContainer.updBiopolymerClass(chainA).difference (endResidueA , startResidueA) + 1;
        SimTK_ERRCHK_ALWAYS(
                (( startResidueA <= endResidueA      )),
                "[ParameterReader.cpp]",  "In the proteinThreading command, the end residue must be greater than or equal to the start residue for each chain.");
        SimTK_ERRCHK_ALWAYS(
                (( startResidueB <= endResidueB      )),
                "[ParameterReader.cpp]",  "In the proteinThreading command, the end residue must be greater than or equal to the start residue for each chain.");
        SimTK_ERRCHK_ALWAYS(
                (myBiopolymerClassContainer.updBiopolymerClass(chainA).difference( endResidueA,startResidueA) == myBiopolymerClassContainer.updBiopolymerClass(chainB).difference ( endResidueB,startResidueB)),
                "[ParameterReader.cpp]",  "In the proteinThreading command, the two threaded segments must be of the same length.");


        for (int i = 0; i < threadedLength; i++) 
        {
            BaseInteraction myBasePair;
            myBasePair.rotationCorrection1 = Rotation(0.0,UnitVec3(0,0,1));
            myBasePair.rotationCorrection2 = Rotation(0.0,UnitVec3(0,0,1));

            myBasePair.FirstBPChain    = chainA;            
            myBasePair.SecondBPChain   = chainB;            
            myBasePair.FirstBPResidue  = myBiopolymerClassContainer.updBiopolymerClass(chainA).sum(startResidueA , i); 
            myBasePair.SecondBPResidue = myBiopolymerClassContainer.updBiopolymerClass(chainB).sum(startResidueB , i); 

            myBasePair.FirstBPEdge     = "Superimpose";      
            myBasePair.SecondBPEdge    = "Superimpose";      

            myBasePair.OrientationBP = "Cis";                 
            //myBasePair.BasePairPriority = 1;

            basePairContainer.addBasePair(myBiopolymerClassContainer, _leontisWesthofClass, myBasePair, setHelicalStacking); 

        }
        return;
    }

    if ( ((parameterStringClass.getString(0)).compare( "singleBondMobility") == 0)) {
        MMBLOG_FILE_FUNC_LINE(ALWAYS, "syntax:  singleBondMobility <chain1> <residue1> <atom1> <mobility> <chain2> <residue2> <atom2>"<<endl);

        SingleBondMobility mySingleBondMobility;
        mySingleBondMobility.chain1   = parameterStringClass.getString(1);
        mySingleBondMobility.residue1 = myBiopolymerClassContainer.residueID(userVariables,parameterStringClass.getString(2),mySingleBondMobility.chain1 );
        mySingleBondMobility.atom1    = parameterStringClass.getString(3);
        mySingleBondMobility.mobility = parameterStringClass.getString(4);
        mySingleBondMobility.chain2   = parameterStringClass.getString(5);
        mySingleBondMobility.residue2 = myBiopolymerClassContainer.residueID(userVariables,parameterStringClass.getString(6),mySingleBondMobility.chain2 );
        mySingleBondMobility.atom2    = parameterStringClass.getString(7);

        // atomPathString used to self-validate.  However it no longer does that, because occasionally we want to create paths that will be validated only later.  In any event, we need to validate explicitly here:
        Compound::AtomPathName myAtomPathName1 = myBiopolymerClassContainer.updBiopolymerClass(mySingleBondMobility.chain1).atomPathString(mySingleBondMobility.residue1, mySingleBondMobility.atom1);
        myBiopolymerClassContainer.updBiopolymerClass(mySingleBondMobility.chain1).validateAtomPathName(myAtomPathName1);
        Compound::AtomPathName myAtomPathName2 = myBiopolymerClassContainer.updBiopolymerClass(mySingleBondMobility.chain2).atomPathString(mySingleBondMobility.residue2, mySingleBondMobility.atom2);
        myBiopolymerClassContainer.updBiopolymerClass(mySingleBondMobility.chain2).validateAtomPathName(myAtomPathName2);

        mobilizerContainer.singleBondMobilityVector.push_back(mySingleBondMobility);

        return;
    } // End singleBondMobility

    if ( ((parameterStringClass.getString(0)).compare( "psiPhiMobility") == 0)) 
    {
        MMBLOG_FILE_FUNC_LINE(ALWAYS,
                "This command sets the bondMobility (Free, Torsion, or Rigid) of the Psi and Phi bonds for all residues in a given range."<<endl
                <<"syntax:  psiPhiMobility <chain> <start residue> <end residue> <mobility> "<<endl
                <<"You can also leave out the residue numbers, and the command will be applied to the whole chain:"<<endl
                <<"syntax:  psiPhiMobility <chain> <mobility> "<<endl
                <<"You can also leave out the chain ID, and the command will be applied to the whole chain, for every protein chain in the system:"<<endl
                <<"syntax:  psiPhiMobility <mobility> "<<endl);
        SingleBondMobility mySingleBondMobility;
        ResidueID startResidue;
        ResidueID endResidue;
        if (parameterStringClass.getString(5).length() != 0) {
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "You have specified too many parameters for this command. "<<endl);
        } // Catch error
        else if ((parameterStringClass.getString(5).length() == 0) &&
                 (parameterStringClass.getString(4).length() != 0)) {
        String chain   = parameterStringClass.getString(1);
        mobilizerContainer.addPhiPsiMobility(   chain, // chain ID
                        myBiopolymerClassContainer.residueID(userVariables,parameterStringClass.getString(2),chain ) , // start residue
                        myBiopolymerClassContainer.residueID(userVariables,parameterStringClass.getString(3),chain ), // end residue
                        parameterStringClass.getString(4), // bond mobility
                        myBiopolymerClassContainer);
        } //  User specified residue numbers explicitly.
        else if ((parameterStringClass.getString(4).length() == 0) &&
             (parameterStringClass.getString(3).length() != 0)) 
        {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "You have specified an incorrect number of parameters. "<<endl);
        }
        else if ((parameterStringClass.getString(3).length() == 0) &&
                 (parameterStringClass.getString(2).length() != 0)) {
            String chain = parameterStringClass.getString(1);
            mobilizerContainer.addPhiPsiMobility(   parameterStringClass.getString(1), // chain ID
                        myBiopolymerClassContainer.updBiopolymerClass(chain).getFirstResidueID(), // start residue
                        myBiopolymerClassContainer.updBiopolymerClass(chain).getLastResidueID(), // end residue
                        parameterStringClass.getString(2), // bond mobility
                        myBiopolymerClassContainer);
        }
        else if ((parameterStringClass.getString(2).length() == 0) &&
                 (parameterStringClass.getString(1).length() != 0)) {
            mobilizerContainer.addPhiPsiMobility(   
                        parameterStringClass.getString(1), // bond mobility
                        myBiopolymerClassContainer);
        } else {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "You have specified an incorrect number of parameters. "<<endl);
        }
        return;
    } //End psiPhiMobility

    if ( ((parameterStringClass.getString(0)).compare( "doubleHelix") == 0)) {
        MMBLOG_FILE_FUNC_LINE(ALWAYS, "Syntax: doubleHelix <chain A> <lowest numbered residue chain A> <highest numbered residue chain A> <chain B> <highest numbered residue chain B> <lowest numbered residue chain B> "<<endl);
        ResidueID lowerA = myBiopolymerClassContainer.residueID(userVariables,parameterStringClass.getString(2),parameterStringClass.getString(1)  );
        ResidueID higherA = myBiopolymerClassContainer.residueID(userVariables,parameterStringClass.getString(3),parameterStringClass.getString(1)  );
        ResidueID lowerB = myBiopolymerClassContainer.residueID(userVariables,parameterStringClass.getString(6),parameterStringClass.getString(4)  );
        ResidueID higherB = myBiopolymerClassContainer.residueID(userVariables,parameterStringClass.getString(5),parameterStringClass.getString(4)  );
        if (higherA <= myBiopolymerClassContainer.updBiopolymerClass(parameterStringClass.getString(1)).sum(lowerA,1)) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "for the first chain, specify first the lower numbered, then the higher numbered residue. Helix must be at least 3BP long."<<endl);
        }
        if (higherB <=  myBiopolymerClassContainer.updBiopolymerClass(parameterStringClass.getString(4)).sum (lowerB,1)) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "for the second chain, specify first the higher numbered, then the lower numbered residue. Helix must be at least 3BP long."<<endl);
        }
        for (int i = 0 ; i <= myBiopolymerClassContainer.updBiopolymerClass(parameterStringClass.getString(1)).difference(higherA,lowerA); i ++) {
            BaseInteraction myBasePair;
            myBasePair.FirstBPEdge = "WatsonCrick";
            myBasePair.SecondBPEdge = "WatsonCrick";
            myBasePair.FirstBPChain = parameterStringClass.getString(1);
            myBasePair.SecondBPChain = parameterStringClass.getString(4);
            myBasePair.FirstBPResidue =myBiopolymerClassContainer.updBiopolymerClass(myBasePair.FirstBPChain ).sum (lowerA , i);
            myBasePair.SecondBPResidue =myBiopolymerClassContainer.updBiopolymerClass(myBasePair.SecondBPChain ).sum (higherB,( - i));
            myBasePair.OrientationBP = "Cis";         
            myBasePair.basePairSatisfied = "False"; //initialize
            if (! myBiopolymerClassContainer.hasChainID(myBasePair.FirstBPChain))
            {
                MMBLOG_FILE_FUNC_LINE(CRITICAL, "No biopolymer found with chain ID "<<myBasePair.FirstBPChain<<".  Please declare all biopolymers before using them in this command.  "<<endl);
            }
            if (! myBiopolymerClassContainer.hasChainID(myBasePair.SecondBPChain))
            {
                MMBLOG_FILE_FUNC_LINE(CRITICAL, "No biopolymer found with chain ID "<<myBasePair.SecondBPChain<<".  Please declare all biopolymers before using them in this command.  "<<endl);
            }
            basePairContainer.addBasePair(myBiopolymerClassContainer, _leontisWesthofClass, myBasePair, setHelicalStacking); 
        }

        return;
    }
    if ( ((parameterStringClass.getString(0)).compare( "modifyBiopolymer") == 0)) {
        MMBLOG_FILE_FUNC_LINE(ALWAYS, "This command enables post transcriptional/translational modifications to canonical nucleic/amino acids.  You would first create a compound with the \"molecule\" command.  Then you would bond that compound to a specified biopolymer residue. "<<endl);
        MMBLOG_FILE_FUNC_LINE(ALWAYS, "Syntax: modifyBiopolymer <biopolymer chain to modify> <residue to modify> <atom to modify> <free bond (e.g. bond1, bond2 ..)> <chain ID of compound to add> <atom on compound to be bonded> <free bond on added compound> "<<endl);
        BiopolymerModification myBiopolymerModification;
        myBiopolymerModification.setChainToModify (parameterStringClass.getString(1));
        ResidueID residueToModify = myBiopolymerClassContainer.residueID(userVariables,parameterStringClass.getString(2),myBiopolymerModification.getChainToModify() );
        myBiopolymerModification.setResidueToModify(residueToModify);
        myBiopolymerModification.setAtomToModifyOnBiopolymer(parameterStringClass.getString(3));
        myBiopolymerModification.setFreeBondOnBiopolymer(parameterStringClass.getString(4));
        myBiopolymerModification.setChainToAdd(parameterStringClass.getString(5));
        myBiopolymerModification.setAtomOnAddedCompound(parameterStringClass.getString(6));
        myBiopolymerModification.setFreeBondOnAddedCompound(parameterStringClass.getString(7)); 
        biopolymerModificationVector.push_back(myBiopolymerModification);
        return;
    }

 //     cout << ((parameterStringClass.getString(0)) << endl;
  
    // PARSE NtC parameters from commands.dat file
    
    //cout << "here " << endl;
    #ifdef BuildNtC  
    if ( ((parameterStringClass.getString(0)).compare("NtC") == 0)) {
	//parameterStringClass.print();
        MMBLOG_FILE_FUNC_LINE(ALWAYS,
                "Syntax: NtC <chain> <start residue> <end residue> <NtC class> <force constant> [meta <secondary weight>] "<<endl
                <<"For example, if (DNA) chain A, residues 1 and 2 are in a B-form helix helix, and you want a force constant of 1.5, you can specify :  "<<endl
                <<"     NtC A 1 2 AA00 1.5  "<<endl);
      
        if (parameterStringClass.getString(5).length() == 0) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "Syntax error! You did not provide enough parameters for this command!"<<endl);
        }
        String myChain = parameterStringClass.getString(1);
        ResidueID firstNtCResidueInStretch = myBiopolymerClassContainer.residueID(userVariables,parameterStringClass.getString(2).c_str(), myChain);
        ResidueID lastNtCResidueInStretch  = myBiopolymerClassContainer.residueID(userVariables,parameterStringClass.getString(3).c_str(), myChain);
	String myNtCClassString = parameterStringClass.getString(4) ;
        double myNtCWeight = myAtoF(userVariables,parameterStringClass.getString(5).c_str());
	bool myMeta = 0;
	double myNtCWeight2=0.;
        if (parameterStringClass.getString(6).length() != 0){
            myMeta=1;
            myNtCWeight2= myAtoF(userVariables,parameterStringClass.getString(7).c_str());
	}
        ntc_class_container.add_NTC_Class(myBiopolymerClassContainer,ntc_par_class,myChain,firstNtCResidueInStretch,lastNtCResidueInStretch,myNtCClassString,myNtCWeight,myMeta,myNtCWeight2);
	/*
        if ( myBiopolymerClassContainer.updBiopolymerClass( myChain ).difference (firstNtCResidueInStretch, lastNtCResidueInStretch ) != -1) {
            MMBLOG_FILE_FUNC_LINE(DEBUG, "NtCs could previously only be applied between consecutive residues. "<<endl);
            //MMBLOG_FILE_FUNC_LINE(DEBUG, "Syntax error! NtCs can currently only be applied between consecutive residues. "<<endl);
        } 
        int firstNtCResidueIndexInStretch = myBiopolymerClassContainer.updBiopolymerClass( myChain ).getResidueIndex(firstNtCResidueInStretch);
        int lastNtCResidueIndexInStretch  = myBiopolymerClassContainer.updBiopolymerClass( myChain ).getResidueIndex(lastNtCResidueInStretch);
        for (int currentFirstResidueIndex = firstNtCResidueIndexInStretch; currentFirstResidueIndex <  lastNtCResidueIndexInStretch; currentFirstResidueIndex += 1 ) 
        {  
            //// The below was above "for":
             NTC_Classes NTC;
             NTC.NtC_FirstBPChain = (myChain);
             NTC.NtC_Class_String = (parameterStringClass.getString(4));
             //NTC.FirstBPResidue =  myBiopolymerClassContainer.residueID(userVariables,parameterStringClass.getString(2).c_str(), NTC.NtC_FirstBPChain);
             // Convert ResidueID to index:
             // These three can   be outside the loop:
             NTC.rotationCorrection1 = Rotation(0.0,UnitVec3(0,0,1));
             NTC.rotationCorrection2 = Rotation(0.0,UnitVec3(0,0,1));
             bool set_ntc_class = true;
             //
             if (! (firstNtCResidueIndexInStretch < lastNtCResidueIndexInStretch) ){
                 MMBLOG_FILE_FUNC_LINE(CRITICAL, "Syntax error! The residue numbers must be ascending! You specified "<<NTC.FirstBPResidue.outString()<<" followed by "<<NTC.SecondBPResidue.outString()<<" . These have residue indices from "<< firstNtCResidueIndexInStretch << " to " << lastNtCResidueIndexInStretch<<endl);
             }

             MMBLOG_FILE_FUNC_LINE(INFO, "About to start NtC loop. Overall residue stretch is from "<<firstNtCResidueInStretch.outString()<< " to "<< lastNtCResidueInStretch.outString()<<" . " <<endl);
             MMBLOG_FILE_FUNC_LINE(INFO, "Residue indices are from "<<firstNtCResidueIndexInStretch<< " to "<< lastNtCResidueIndexInStretch<<" . " <<endl);
             //int currentFirstResidueIndex = firstNtCResidueIndexInStretch;
             MMBLOG_FILE_FUNC_LINE(INFO, "currentFirstResidueIndex = "<<currentFirstResidueIndex<<endl);
                 ////
            MMBLOG_FILE_FUNC_LINE(INFO, "currentFirstResidueIndex = "<<currentFirstResidueIndex<<endl);
            MMBLOG_FILE_FUNC_LINE(INFO, "currentFirstResidueIndex + 1 = "<<currentFirstResidueIndex + 1<<endl);
            if ((currentFirstResidueIndex + 1 ) > lastNtCResidueIndexInStretch){ // Don't see how this could happen, but being ultra paranoid.
                MMBLOG_FILE_FUNC_LINE(INFO, "The index: "<<(currentFirstResidueIndex + 1 )<<" of the second residue in this NtC, is too high!"<<endl);
                MMBLOG_FILE_FUNC_LINE(INFO, "Compared to the index: "<<lastNtCResidueIndexInStretch<<" of the last residue in the range."<<endl);
                MMBLOG_FILE_FUNC_LINE(CRITICAL, "Unexplained error! The index: "<<(currentFirstResidueIndex + 1 )<<" of the second residue in this NtC, is too high! "<<endl);
            }
            NTC.FirstBPResidue  = myBiopolymerClassContainer.updBiopolymerClass(NTC.NtC_FirstBPChain).getResidueID(currentFirstResidueIndex     ); 
            NTC.SecondBPResidue = myBiopolymerClassContainer.updBiopolymerClass(NTC.NtC_FirstBPChain).getResidueID(currentFirstResidueIndex + 1 );
            NTC.NtC_step_ID = (NTC.FirstBPResidue.outString());
            MMBLOG_FILE_FUNC_LINE(INFO, "Starting NtC loop. Overall residue stretch is from "<<firstNtCResidueInStretch.outString()<< " to "<< lastNtCResidueInStretch.outString()<<" . In this round, NTC.FirstBPResidue = "<< NTC.FirstBPResidue.outString() << " , NTC.SecondBPResidue = "<< NTC.SecondBPResidue.outString() <<" . "<<endl);
            MMBLOG_FILE_FUNC_LINE(INFO, "NTC.FirstBPResidue = "<<NTC.FirstBPResidue.outString()<<endl);
            MMBLOG_FILE_FUNC_LINE(INFO, "NTC.SecondBPResidue = "<<NTC.SecondBPResidue.outString()<<endl);
            //NTC.weight          = stod(parameterStringClass.getString(5));
            NTC.weight = myAtoF(userVariables,parameterStringClass.getString(5).c_str());
            MMBLOG_FILE_FUNC_LINE(INFO, "NTC.weight = "<<NTC.weight<<endl);

            if (( NTC.weight > 20.0 ) && ( safeParameters )){ // Empirically found that a weight greater than 20 or so leads to strange jumpy nonconvergent behavior.
                 MMBLOG_FILE_FUNC_LINE(WARNING, "You have specified an NtC weight of " << NTC.weight <<" . The current NtC parameters are overconstrained, and unless these have been updated you may need to decrease the weight so as to avoid instability.  "<<endl);
             }
 
            NTC.meta            = 0;
            int metaPosition = 6;        
            if (parameterStringClass.getString(metaPosition).length() != 0){
                if (parameterStringClass.getString(metaPosition + 1).length() == 0){
                    MMBLOG_FILE_FUNC_LINE(CRITICAL, "Syntax error! Expected parameter at position "<<metaPosition + 1<<" since you have a parameter at position "<<metaPosition<<"."<<endl);
                }
                if((parameterStringClass.getString(metaPosition)).compare("meta") == 0){
                    NTC.meta = 1;
                    NTC.weight2     = stod(parameterStringClass.getString(metaPosition + 1));
                    // These three do need to be in this loop:
                    NTC.count = myBiopolymerClassContainer.count;
                    myBiopolymerClassContainer.count++;
                    MMBLOG_FILE_FUNC_LINE(INFO, NTC.count << " number of NTC meta input lines " << endl);
                    //
                } else {
                    MMBLOG_FILE_FUNC_LINE(CRITICAL, "Syntax error! Expected parameter \"meta\", or nothing at all, position "<<metaPosition<<"."<<endl);
                }
            } // of if
            
            ntc_class_container.add_NTC_Class(myBiopolymerClassContainer,ntc_par_class,NTC,set_ntc_class);
            MMBLOG_FILE_FUNC_LINE(INFO, "At end of loop. Just added NTC with NTC.FirstBPResidue  = "<<NTC.FirstBPResidue.outString() << " and NTC.SecondBPResidue = "<< NTC.SecondBPResidue.outString()<<endl);
            NTC.print();
            //currentFirstResidueIndex++;
        } 
        */
        return;
    };




    #endif 
    // END PARSE NtC parameters from commands.dat file   
    
    if ((  (parameterStringClass.getString(0)).compare("baseInteraction") == 0) || ((parameterStringClass.getString(0)).compare("aromatic") == 0) )    
    { // if this is a base pair
        BaseInteraction myBasePair;
        myBasePair.rotationCorrection1 = Rotation(0.0,UnitVec3(0,0,1));
        myBasePair.rotationCorrection2 = Rotation(0.0,UnitVec3(0,0,1));
        if (((parameterStringClass.getString(0)).compare("aromatic") == 0) ) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "The \"aromatic\" command is no longer supported. "<<endl);
        }
        if (((parameterStringClass.getString(0)).compare("baseInteraction") == 0) ) {
            if (parameterStringClass.getString(6).length() == 0) {
                MMBLOG_FILE_FUNC_LINE(CRITICAL, "You must specify the second interacting \"edge\" for the baseInteraction command. "<<endl);
            }
            if (parameterStringClass.getString(7).length() == 0) {
                MMBLOG_FILE_FUNC_LINE(CRITICAL, "You must specify the orientation (typically Cis or Trans) for the baseInteraction command. "<<endl);
            }
        }
        { 
            myBasePair.FirstBPChain = ((parameterStringClass.getString(1)));            
            myBasePair.FirstBPResidue =  myBiopolymerClassContainer.residueID(userVariables,parameterStringClass.getString(2).c_str(), myBasePair.FirstBPChain  );

            myBasePair.FirstBPEdge = parameterStringClass.getString(3);
            if ((myBasePair.FirstBPEdge.compare("ChiBondAnti") ==0) ){
                MMBLOG_FILE_FUNC_LINE(CRITICAL, "ChiBondAnti is not currently supported.  Making sure you don't use it.."<<endl);}
            myBasePair.SecondBPChain = ((parameterStringClass.getString(4)));            
            myBasePair.SecondBPResidue =  myBiopolymerClassContainer.residueID(userVariables,parameterStringClass.getString(5).c_str(), myBasePair.SecondBPChain  );
            myBasePair.SecondBPEdge = parameterStringClass.getString(6);            
            myBasePair.OrientationBP = parameterStringClass.getString(7);           
            if (parameterStringClass.getString(8).length() > 0) {  
                MMBLOG_FILE_FUNC_LINE(CRITICAL, "Too many parameters!  Priorities are no longer supported for baseInteraction\'s. "<<endl);
            }

            MMBLOG_FILE_FUNC_LINE(DEBUG, "Just read -BP:"<<myBasePair.FirstBPChain    <<" "<<myBasePair.FirstBPResidue.outString()<<endl);
            MMBLOG_FILE_FUNC_LINE(DEBUG, myBasePair.SecondBPChain <<" "<<myBasePair.SecondBPResidue.outString()<<endl);
            if ((parameterStringClass.getString(0).compare("constraint") == 0) || (parameterStringClass.getString(0).compare("restraint") == 0)) 
                if (myBasePair.FirstBPEdge.compare("Weld")==0){ // all is well
                } else {
                    MMBLOG_FILE_FUNC_LINE(CRITICAL, "The only type of "<<parameterStringClass.getString(0)<<" permitted is Weld.  You have specified : "<< myBasePair.FirstBPEdge<<endl);
                }
            if (safeParameters) {

                SimTK_ERRCHK_ALWAYS(
                        (matchDefaultSkipTopLevelTransform == 0),
                        "[ParameterReader.cpp]",
                        "We don't recommend matchDefaultSkipTopLevelTransform = 1   ");

                {

                    if (!
                            (
                             myMonoAtomsContainer.hasChainID(myBasePair.FirstBPChain) ||
                             myBiopolymerClassContainer.hasChainID(myBasePair.FirstBPChain)
                            ))
                    {MMBLOG_FILE_FUNC_LINE(CRITICAL, "You have not declared the chain "<<myBasePair.FirstBPChain<<" yet, or the chain is empty.  "<<endl);}

                    if (!
                            (
                             myMonoAtomsContainer.hasChainID(myBasePair.SecondBPChain) ||
                             myBiopolymerClassContainer.hasChainID(myBasePair.SecondBPChain)
                            ))
                    {MMBLOG_FILE_FUNC_LINE(CRITICAL, "You have not declared the chain "<<myBasePair.SecondBPChain<<" yet, or the chain is empty.  "<<endl);}
                }

                SimTK_ERRCHK_ALWAYS((myBasePair.FirstBPEdge.compare("BackboneOxygen") != 0 ),"[ParameterReader.cpp]","BackboneOxygen is no longer supported.");
            }
        }
        r++; 
        MMBLOG_FILE_FUNC_LINE(INFO, "r incremented to "<<r<<endl);
        myBasePair.basePairSatisfied = "False"; //initialize
        
        basePairContainer.addBasePair(myBiopolymerClassContainer, _leontisWesthofClass, myBasePair, setHelicalStacking); 
        return;
    } // End baseInteraction

    if ( ((parameterStringClass.getString(0)).compare("mobilizer") == 0)  )    
    { // if this is a mobilizer or constraint
        MMBLOG_FILE_FUNC_LINE(ALWAYS,
                "Syntax: mobilizer <Bond Mobility> <chain> <start residue> <end residue>"<<endl
                <<"Or:     mobilizer <Bond Mobility> <chain> .. to set all residues in <chain> to <Bond Mobility>"<<endl
                <<"Or:     mobilizer <Bond Mobility>  .. to set all residues in all chains to <Bond Mobility>"<<endl
                <<"Where Bond Mobility may be RNA, DNA, or Protein"<<endl);

        if (parameterStringClass.getString(1).length() == 0) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "You have not specified enough parameters for this command. "<<endl);

        } else if (parameterStringClass.getString(2).length() == 0) { // a command e.g. mobilizer Rigid will set all residues in all chains to BondMobility::Rigid.
            if (myBiopolymerClassContainer.getNumBiopolymers() == 0) {
		MMBLOG_FILE_FUNC_LINE(CRITICAL, "The number of biopolymers detected at this point is  "<<myBiopolymerClassContainer.getNumBiopolymers() <<". Please call the mobilizer command after the chain(s) in question have been instantitated."<<endl);
            
            }
            mobilizerContainer.setMobilizerTypeForAllChains(parameterStringClass.getString(1), myBiopolymerClassContainer); // That first argument is the mobilizer type string
            /*for (int i = 0 ; i < myBiopolymerClassContainer.getNumBiopolymers(); i++) {
                String myChainID = myBiopolymerClassContainer.updBiopolymerClass(i).getChainID();
                String myMobilizerString = parameterStringClass.getString(1);
                //MMBLOG_FILE_FUNC_LINE(" myMobilizerString = >"<<myMobilizerString<<"< "<<endl;
                MMBLOG_FILE_FUNC_LINE(" Adding mobilizer stretch to biopolymer index "<<i<<endl;
                mobilizerContainer.addMobilizerStretchToVector(
                        myChainID, 
                        myMobilizerString, //parameterStringClass.getString(1), 
                        myBiopolymerClassContainer
                        );
            }*/
        } else if (parameterStringClass.getString(3).length() == 0) {
            mobilizerContainer.addMobilizerStretchToVector(parameterStringClass.getString(2),parameterStringClass.getString(1),myBiopolymerClassContainer);
        } else if (parameterStringClass.getString(5).length() == 0) {
            MMBLOG_FILE_FUNC_LINE(INFO, "first residue ID: "<<    myBiopolymerClassContainer.residueID(userVariables,parameterStringClass.getString(3), parameterStringClass.getString(2)  ).outString()<<endl);
            MMBLOG_FILE_FUNC_LINE(INFO, "second residue ID: "<<    myBiopolymerClassContainer.residueID(userVariables,parameterStringClass.getString(4), parameterStringClass.getString(2)  ).outString()<<endl);
            String myChainID = parameterStringClass.getString(2);
            String myMobilizerString = parameterStringClass.getString(1);
            mobilizerContainer.addMobilizerStretchToVector( myChainID, //parameterStringClass.getString(2),

                    myBiopolymerClassContainer.residueID(userVariables,parameterStringClass.getString(3), parameterStringClass.getString(2)  ),
                    myBiopolymerClassContainer.residueID(userVariables,parameterStringClass.getString(4), parameterStringClass.getString(2)  ),
                    myMobilizerString, 
                    myBiopolymerClassContainer);
        } else {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "You have specified the wrong number of parameters for this command. "<<endl);
        }
        return;
    } // End mobilizer
    if ( ((parameterStringClass.getString(0)).compare("rotation") == 0)  )    { 
        MMBLOG_FILE_FUNC_LINE(ALWAYS,
                "This command applies a rotation to the specified chain, prior to start of the time integrator.> "<<endl
                <<"Syntax: rotation <chain> <axis about which to rotate, X|Y|Z> <angle, in radians.> "<<endl
                <<"Where <chain> is the chain to be displaced, the axis is that about which you will rotate, and the <angle> is the rotation angle, right handed around the named axis. "<<endl
                <<"Rotations will be applied in the order specified in the command file.        "<<endl<<endl);
        parameterStringClass.validateNumFields(4);
        String myChain = parameterStringClass.getString(1).c_str();
        String myAxis =  parameterStringClass.getString(2).c_str();
        CoordinateAxis myCoordinateAxis(0); // The only available constructors obligate me to choose an axis at contruct time. Make sure we later override this.
        if       ((myAxis == "X") || (myAxis == "x")){myCoordinateAxis = CoordinateAxis(0);}
        else if  ((myAxis == "Y") || (myAxis == "y")){myCoordinateAxis = CoordinateAxis(1);}
        else if  ((myAxis == "Z") || (myAxis == "z")){myCoordinateAxis = CoordinateAxis(2);}
        double myAngle = myAtoF(userVariables,parameterStringClass.getString(3).c_str());
        Rotation myRotation(myAngle, (myCoordinateAxis));
        MMBLOG_FILE_FUNC_LINE(INFO, " Current rotation : "<<displacementContainer.updDisplacement(myChain).rotation<<endl);
        displacementContainer.updDisplacement(myChain).rotation = myRotation*displacementContainer.updDisplacement(myChain).rotation; // updDisplacement will spit out an error if the displacement has not been created. Here, I am multiplying on the left by the user-supplied rotation. This means I can keep applying rotations and they will always be progressively multiplied from the left.
        MMBLOG_FILE_FUNC_LINE(INFO, " Just set rotation to : "<<displacementContainer.updDisplacement(myChain).rotation<<endl);
        return;
    
    }
    if ( ((parameterStringClass.getString(0)).compare("initialDisplacement") == 0)  )    { 
        MMBLOG_FILE_FUNC_LINE(ALWAYS,
                "Syntax: initialDisplacement <chain> <X> <Y> <Z> "<<endl
                <<"Where <chain> is the chain to be displaced, and the next 3 parameters are the displacement vector."<<endl
                <<"You can also leave out <chain>, and all chains will be spread out in a line: "<<endl
                <<"Syntax: initialDisplacement <X> <Y> <Z> "<<endl
                <<"Where the next 3 parameters are the displacement vector, and each Biopolymer chain will be displaced by i *  (X,Y,Z) , where i is the chain index. This is particularly useful for spreading out the chains to make rendering easier."<<endl
                <<"Note that in MMB 2.10 and earlier, we took the displacement in Å.  We are going back to nm for consistency, with apologies for the confusion."<<endl);
        //verified nm units are respected here

        if (parameterStringClass.getString(3).length() == 0) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "You have not specified enough parameters for this command. "<<endl);
        }
        else if (parameterStringClass.getString(4).length() == 0) {
            for (int i = 0 ; i < myBiopolymerClassContainer.getNumBiopolymers() ; i ++) {
                Displacement myDisplacement;
                // Warning! We are initializing the rotation matrix to the identity matrix. This means that if it has already been set, it will be overwritten.
                Rotation myRotation;
                myRotation.setRotationToIdentityMatrix ();
                myDisplacement.rotation = myRotation;
                myDisplacement.displacement = Vec3( 
                    myAtoF(userVariables,parameterStringClass.getString(1).c_str()),
                    myAtoF(userVariables,parameterStringClass.getString(2).c_str()),
                    myAtoF(userVariables,parameterStringClass.getString(3).c_str())
                    );
                myDisplacement.displacement *= i;
                myDisplacement.chain = myBiopolymerClassContainer.updBiopolymerClass(i).getChainID();
                displacementContainer.add(myDisplacement, myBiopolymerClassContainer);                
            } 
        }
        else if (parameterStringClass.getString(5).length() == 0) {
	    Displacement myDisplacement;
            Rotation myRotation;
            myRotation.setRotationToIdentityMatrix ();
            myDisplacement.rotation = myRotation;
	    myDisplacement.chain = parameterStringClass.getString(1);
	    myDisplacement.displacement = Vec3(
		    myAtoF(userVariables,parameterStringClass.getString(2).c_str()),
		    myAtoF(userVariables,parameterStringClass.getString(3).c_str()),
		    myAtoF(userVariables,parameterStringClass.getString(4).c_str())
		    );
	    displacementContainer.add(myDisplacement, myBiopolymerClassContainer); // This call validateDisplacement which among other things makes sure there is not already one displacement in the container.
        }
        else {//if (parameterStringClass.getString(6).length() != 0) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "You have specified too many parameters for this command. "<<endl);
        }
        return;
    } // End initialDisplacement


    if ((parameterStringClass.getString(0)).compare("constrainChainRigidSegments") == 0) {
        MMBLOG_FILE_FUNC_LINE(ALWAYS,
                "To constrain all rigid segments in chain <chain ID> to a specified residue in that chain, issue:  "<<endl
                <<"Syntax:  constrainChainRigidSegments <chain ID> <residue ID to constrain to>  "<<endl
                <<"To constrain all rigid segments in chain <chain ID> to ground, issue:  "<<endl
                <<"Syntax:  constrainChainRigidSegments <chain ID> Ground   "<<endl
                <<"To constrain all rigid segments in all chains to ground, issue:  "<<endl
                <<"Syntax:  constrainChainRigidSegments   "<<endl);

        if (parameterStringClass.getString(3).length()>0){
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "You have provided too many parameters for this command."<<endl);
        }
        MMBLOG_FILE_FUNC_LINE(INFO, endl);
        if (parameterStringClass.getString(2).length() > 0){
            ChainResidueToGround myChainResidue;
            myChainResidue.chainID = parameterStringClass.getString(1);
            myBiopolymerClassContainer.validateChainID(myChainResidue.chainID);
            if (parameterStringClass.getString(2).toUpper().compare("GROUND") != 0) { // if we're not constraining to Ground, but rather to some residue of the same chain:
                myChainResidue.residueID  = myBiopolymerClassContainer.residueID(userVariables, parameterStringClass.getString(2) , myChainResidue.chainID); //BiopolymerClassContainer.residueID supports e.g. FirstResidue, LastResidue.
                myChainResidue.toGround = false;
            } else {
                myChainResidue.toGround = true;
            }

            constraintToGroundContainer.queueConstrainChainRigidSegments (myChainResidue );
        } else if (parameterStringClass.getString(1).length() == 0) {  
            //MMBLOG_FILE_FUNC_LINE(endl;

	    //map<const String, BiopolymerClass>::iterator biopolymerClassMapIterator = myBiopolymerClassContainer.getBiopolymerClassMap().begin();
            MMBLOG_FILE_FUNC_LINE(INFO, "Detected that you wish to constrain rigid segments in the following :"<<endl);
            for (int i = 0; i < myBiopolymerClassContainer.getNumBiopolymers(); i++)
	    //for(biopolymerClassMapIterator = myBiopolymerClassContainer.getBiopolymerClassMap().begin(); biopolymerClassMapIterator != myBiopolymerClassContainer.getBiopolymerClassMap().end(); biopolymerClassMapIterator++) 
            {
                 MMBLOG_FILE_FUNC_LINE(INFO, "chain : "<<  myBiopolymerClassContainer.updBiopolymerClass(i).getChainID()  <<endl);
		 ChainResidueToGround myChainResidue;
                 myChainResidue.chainID = myBiopolymerClassContainer.updBiopolymerClass(i).getChainID();
                 myChainResidue.toGround = true;
		 constraintToGroundContainer.queueConstrainChainRigidSegments(myChainResidue );
                 //MMBLOG_FILE_FUNC_LINE(endl;
	    }
            //MMBLOG_FILE_FUNC_LINE(endl;
            //constrainRigidSegments = true; // discontinue this parameter
        } else if (parameterStringClass.getString(1).length() != 0) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "You have provided the wrong number of parameters for this command."<<endl);
        } else {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "Unexplained error! "<<endl);
        }
        
        MMBLOG_FILE_FUNC_LINE(INFO, endl);
        constraintToGroundContainer.printConstraintClasses(); // This might still be empty at this point..
        MMBLOG_FILE_FUNC_LINE(INFO, endl);
        return;
    }

    if ((parameterStringClass.getString(0)).compare("constraint") == 0) {
        MMBLOG_FILE_FUNC_LINE(ALWAYS,
                "Syntax:  constraint <chain> <atom name> Weld Ground "<<endl
                <<":      This one works for  custom molecules only. "<<endl
                <<": Syntax :  constraint <chain 1> <residue ID 1> Weld Ground "<<endl
                <<":      This one works for  biopolymers only. "<<endl
                <<": Or:  constraint <chain 1> <residue ID 1> <atom name 1> Weld Ground "<<endl
                <<":      This one works for  biopolymers only. "<<endl
                <<": Or:  constraint <chain 1> <residue ID 1> Weld <chain 2> <residue ID 2> "<<endl
                <<": Or:  constraint <chain 1> <residue ID 1> <atom name 1> Weld <chain 2> <residue ID 2> <atom name 2> "<<endl
                <<":      This one treats atoms from monoAtoms and biopolymers on an equal footing. It does not work for custom molecules though. "<<endl
                <<parameterStringClass.getString() <<  endl);

        if ( (parameterStringClass.getString(4).length()>0) &&(!(parameterStringClass.getString(5).length() > 0))){  // This means syntax:  constraint <chain 1> <residue ID 1> Weld Ground
            String myChain1 = parameterStringClass.getString(1);
            if (myBiopolymerClassContainer.hasChainID(myChain1)) {
                MMBLOG_FILE_FUNC_LINE(INFO, "detected syntax for biopolymers:  constraint <chain 1> <residue ID 1> Weld Ground "<<endl);
		if (!( (parameterStringClass.getString(3)).compare("Weld") == 0)) {
		    MMBLOG_FILE_FUNC_LINE(CRITICAL, "Syntax error!"<<endl);
		}
		if (!( (parameterStringClass.getString(4)).compare("Ground") == 0)) {
		    MMBLOG_FILE_FUNC_LINE(CRITICAL, "Syntax error!"<<endl);
		}
		if (removeRigidBodyMomentum) {
		    MMBLOG_FILE_FUNC_LINE(CRITICAL, "You must first set removeRigidBodyMomentum False before welding to ground ! "<<endl);
		}
		String myResidue1String = parameterStringClass.getString(2);
		myBiopolymerClassContainer.addConstraintToGround(userVariables, myResidue1String,myChain1,constraintToGroundContainer);
            }
            else if (moleculeClassContainer.hasChainID(myChain1)) {
                MMBLOG_FILE_FUNC_LINE(INFO, "Detected syntax for custom molecules:  constraint <chain> <atom name> Weld Ground "<<endl);
		if (!( (parameterStringClass.getString(3)).compare("Weld") == 0)) {
		    MMBLOG_FILE_FUNC_LINE(CRITICAL, "Syntax error!"<<endl);
		}
		if (!( (parameterStringClass.getString(4)).compare("Ground") == 0)) {
		    MMBLOG_FILE_FUNC_LINE(CRITICAL, "Syntax error!"<<endl);
		}
		if (removeRigidBodyMomentum) {
		    MMBLOG_FILE_FUNC_LINE(CRITICAL, "You must first set removeRigidBodyMomentum False before welding to ground ! "<<endl);
		}
		String myAtomName = parameterStringClass.getString(2);
		moleculeClassContainer.addConstraintToGround(userVariables,myChain1,myAtomName, constraintToGroundContainer);
	    } else {
		    MMBLOG_FILE_FUNC_LINE(CRITICAL, "Chain "<<myChain1<<" is neither a biopolymer nor a custom molecule!"<<endl);
		}
	    } // End Weld Ground syntax
            else if ( (parameterStringClass.getString(4) == "Weld") &&
                 (parameterStringClass.getString(5) == "Ground") &&
                 (!(parameterStringClass.getString(6).length() > 0))
               ){  // This means syntax:  constraint <chain ID> <residue ID> <atom name> Weld Ground 
                MMBLOG_FILE_FUNC_LINE(INFO, "detected syntax for biopolymers:  constraint <chain 1> <residue ID 1> Weld Ground "<<endl);
                String myChain1 = parameterStringClass.getString(1);
                String myResidue1String = parameterStringClass.getString(2);
                String myAtomNameString = parameterStringClass.getString(3);
                myBiopolymerClassContainer.validateChainID(myChain1);
                myBiopolymerClassContainer.updBiopolymerClass(myChain1).validateResidueID(ResidueID(myResidue1String)  );
                if (!(myBiopolymerClassContainer.updBiopolymerClass(myChain1).hasAtom(ResidueID(myResidue1String), myAtomNameString ))){
                    MMBLOG_FILE_FUNC_LINE(CRITICAL, " atom name "<< myAtomNameString << " not found for this chain "<< myChain1<<" and residues ID "<< myResidue1String<<endl);
                }
                myBiopolymerClassContainer.addConstraintToGround(userVariables, myResidue1String,myChain1,myAtomNameString, constraintToGroundContainer);
                MMBLOG_FILE_FUNC_LINE(INFO, "You now successfully added a constraint for atom name "<<myAtomNameString<<endl);

            }

	    else if (!(parameterStringClass.getString(5).length()>0)) { 
		MMBLOG_FILE_FUNC_LINE(CRITICAL, "Wrong number of parameters!"<<endl);
	    }
	    else if (parameterStringClass.getString(3).compare("Weld") == 0) 
	    { // this is teh case of "constraint <chain 1> <residue ID 1> Weld <chain 2> <residue ID 2>"
		    myBiopolymerClassContainer.addConstraint(
			    userVariables,
			    myBiopolymerClassContainer.updBiopolymerClass(parameterStringClass.getString(1)).getRepresentativeAtomName(), parameterStringClass.getString(2), parameterStringClass.getString(1), 
			    myBiopolymerClassContainer.updBiopolymerClass(parameterStringClass.getString(4)).getRepresentativeAtomName(), parameterStringClass.getString(5), parameterStringClass.getString(4), 
			    constraintToGroundContainer);
	    } // end syntax with two residue numbers and no atom names

	    else if (!(parameterStringClass.getString(7).length()>0)){
		MMBLOG_FILE_FUNC_LINE(CRITICAL, "Wrong number of parameters! This syntax is not supported for now."<<endl);
	    }    
	    else if (parameterStringClass.getString(7).length()>0){ // this is the case of "constraint <chain 1> <residue ID 1> <atom name 1> Weld <chain 2> <residue ID 2> <atom name 2> "
                MMBLOG_FILE_FUNC_LINE(INFO, endl);
		myBiopolymerClassContainer.addConstraint(
			userVariables,
			parameterStringClass.getString(3), parameterStringClass.getString(2), parameterStringClass.getString(1), 
			parameterStringClass.getString(7), parameterStringClass.getString(6), parameterStringClass.getString(5), 
			constraintToGroundContainer);
	    }
	    else if (parameterStringClass.getString(8).length()>0){
		MMBLOG_FILE_FUNC_LINE(CRITICAL, "Wrong number of parameters!"<<endl);
	    }
        return;
    }
    String constrainInterfacesCommand("constrainInterfaces");
    String detectInterChainClashesCommand("detectInterChainClashes");
    if ((((parameterStringClass.getString(0)).compare(constrainInterfacesCommand) == 0))  ||
        (((parameterStringClass.getString(0)).compare(detectInterChainClashesCommand) == 0)))  
        {
        MMBLOG_FILE_FUNC_LINE(ALWAYS,
                "Syntax variants and explanations follow :"<<endl
                <<"1. Syntax : "<<parameterStringClass.getString(0) <<"  <depth (nm)>   "<<endl
                <<"2. Syntax :  "<<parameterStringClass.getString(0) <<"  <depth (nm)>  <chain>  "<<endl
    
                <<"3. Syntax :  "<<parameterStringClass.getString(0) <<"  <depth (nm)>  <chain 1> [<chain 2> [<chain 3> [ ...etc]]]  "<<endl
                <<"4. Syntax :  "<<parameterStringClass.getString(0) <<"  <depth (nm)>  <chain 1> [<chain 2> [<chain 3> [ ...etc]]]  Versus <chain 1> [<chain 2> [<chain 3> [ ...etc]]]  "<<endl);
        if (((parameterStringClass.getString(0)).compare(constrainInterfacesCommand) == 0)) {
	    MMBLOG_FILE_FUNC_LINE(ALWAYS,
                "1. This constrains all  chains to all  chains, if they are within a distance of <depth (nm)>. The constraint (always of type Weld) is applied to the pair of atoms spanning the two chains at closest approach (in terms of internuclear distance).  "<<endl
	            <<"2. This constrains the specified <chain> to all other chains within a distance of <depth (nm)>. The constraint (always of type Weld) is applied to the pair of atoms spanning the two chains at closest approach (in terms of internuclear distance).  "<<endl
	            <<"3. You can also specify a multiple-chain set, and constrain all specified chains to all the remaining chains in the system:"<<endl
	            <<"4. You can similarly specify TWO multiple-chain sets, and apply constraints between all possible pairs spanning the two sets. This is useful e.g. in the ribosome, where 23S binds a different set of proteins than 16S. "<<endl);
        }
        if (((parameterStringClass.getString(0)).compare(detectInterChainClashesCommand) == 0)) {
	    MMBLOG_FILE_FUNC_LINE(ALWAYS,
                "1. This detects clashes between all vs. all chains, if they are within a distance of <depth (nm)>.   "<<endl
	            <<"2. This detects clahes between the specified <chain> and all other chains within a distance of <depth (nm)>.  "<<endl
	            <<"3. You can also specify a multiple-chain set, and detect clashes between  all specified chains to all the remaining chains in the system:"<<endl
	            <<"4. You can similarly specify TWO multiple-chain sets, and detect clashes between all possible pairs spanning the two sets.  "<<endl);
        }

        double depth = myAtoF(userVariables,parameterStringClass.getString(1).c_str());
        vector<String> chains;        chains.clear();
        vector<String> partnerChains; partnerChains.clear();

        if (parameterStringClass.getString(1).length() == 0) {
		MMBLOG_FILE_FUNC_LINE(CRITICAL, "Not enough parameters specified for this command!"<<endl);
        }
        if (parameterStringClass.getString(2).length() == 0) {
            MMBLOG_FILE_FUNC_LINE(INFO, "No chains specified! This must mean you want to constrain all vs. all."<<endl);
            // leave chains and partnerChains empty .. 
            /*for (int i = 0 ; i < myBiopolymerClassContainer.getNumBiopolymers(); i++) {
                chains.push_back       (myBiopolymerClassContainer.updBiopolymerClass(i).getChainID());
                MMBLOG_FILE_FUNC_LINE(" Added chain "<<myBiopolymerClassContainer.updBiopolymerClass(i).getChainID() <<" to first set "<<endl;
                partnerChains.push_back(myBiopolymerClassContainer.updBiopolymerClass(i).getChainID());
                MMBLOG_FILE_FUNC_LINE(" Added chain "<<myBiopolymerClassContainer.updBiopolymerClass(i).getChainID() <<" to second set "<<endl;
            }*/
        }
        int i = 2; // position of first possible chain ID.
        String versusString = String("Versus").toUpper();
        while ((parameterStringClass.getString(i).length() >0)  &&
               (parameterStringClass.getString(i).toUpper().compare(versusString ) != 0 )) { // if parameter i is "versus", stop and start reading partnerChains
            myBiopolymerClassContainer.validateChainID(parameterStringClass.getString(i)); // Make sure this is a valid chain
            chains.push_back( parameterStringClass.getString(i));
            MMBLOG_FILE_FUNC_LINE(INFO, "Added chain "<<parameterStringClass.getString(i)<<" to first set "<<endl);
            i++;
        }
        if (parameterStringClass.getString(i).toUpper().compare(versusString) == 0 ) { 
            i++;
            MMBLOG_FILE_FUNC_LINE(INFO, "Keyword "<<versusString<<" detected."<<endl);
            if (parameterStringClass.getString(i).length() == 0) {
		MMBLOG_FILE_FUNC_LINE(CRITICAL, "The \"versus\" keyword must be followed by at least one chain ID."<<endl);
            }
        } 
        while (parameterStringClass.getString(i).length() >0) {
            myBiopolymerClassContainer.validateChainID(parameterStringClass.getString(i)); // Make sure this is a valid chain
            partnerChains.push_back( parameterStringClass.getString(i));
            MMBLOG_FILE_FUNC_LINE(INFO, "Added chain "<<parameterStringClass.getString(i)<<" to second set "<<endl);
            i++;
        }
        
        if (partnerChains.size() == 0) {
            MMBLOG_FILE_FUNC_LINE(INFO, "No partner chains specified! This must mean you want to constrain the first set vs. all possible chains."<<endl);
        }
        if (((parameterStringClass.getString(0)).compare(constrainInterfacesCommand) == 0)) {
            constraintToGroundContainer.interfaceContainer.addInterface(chains, partnerChains, depth); // The remaining  parameter ( bondMobilityString) is  irrelevant and has  default values.   
        }
        if (((parameterStringClass.getString(0)).compare(detectInterChainClashesCommand) == 0)) {
            contactInterfaceContainer.addInterface(chains, partnerChains, depth); 
        }

        return;
    }

    if ((parameterStringClass.getString(0)).compare("coupleAtomMobilizers") == 0) {
        MMBLOG_FILE_FUNC_LINE(INFO, "This constrains the Torsion mobilizers on two given atoms to be equal.   "<<endl);
        MMBLOG_FILE_FUNC_LINE(INFO, "Syntax:  coupleAtomMobilizers <chain 1> <residue ID 1> <atom name 1> <chain 2> <residue ID 2> <atom name 2>  "<<endl);
        
        if ((parameterStringClass.getString(6).length() == 0) ||
            (parameterStringClass.getString(7).length() >  0))
        {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "Wrong number of parameters!"<<endl);
        }
        myBiopolymerClassContainer.addConstraint(userVariables, parameterStringClass.getString(3), parameterStringClass.getString(2), parameterStringClass.getString(1),
                                                                parameterStringClass.getString(6), parameterStringClass.getString(5), parameterStringClass.getString(4),
                                                                CoupledCoordinate,
                                constraintToGroundContainer);
                                

    // follows this prototype:
    /*void BiopolymerClassContainer::addConstraint(map<const String,double> myUserVariables,
               const String atomName1, const String inputResidueString1,const  String chain1, 
               const String atomName2, const String inputResidueString2,const  String chain2, 
               ConstraintType myConstraintType,
               ConstraintToGroundContainer & constraintToGroundContainer)*/




        return;
    }

    if (((parameterStringClass.getString(0)).compare("coupleMobilizers")      == 0) ||  
        ((parameterStringClass.getString(0)).compare("couplePsiPhiAngles")      == 0))  
    {   
        MMBLOG_FILE_FUNC_LINE(ALWAYS, "This constrains the Torsion mobilizers corresponding to Psi and Phi angles of <chain A>, residues  <start-residue-A> to <end-residue-A>, to be equal to the equivalent mobilizers in <Chain B>, residues <start-residue-B> to <end-residue-B>."<<endl);
        //MMBLOG_FILE_FUNC_LINE(" : This constrains the Torsion mobilizers of all atoms in <chain A>, residues  <start-residue-A> to <end-residue-A>, to be equal to the equivalent mobilizers in <Chain B>, residues <start-residue-B> to <end-residue-B>."<<endl;
        //MMBLOG_FILE_FUNC_LINE(" : Syntax:  coupleMobilizers <Chain A> <start-residue-A> <end-residue-A> <Chain-B> <start-residue-B> <end-residue-B> "<<endl;
        MMBLOG_FILE_FUNC_LINE(ALWAYS, "Syntax:  couplePsiPhiAngles <Chain A> <start-residue-A> <end-residue-A> <Chain-B> <start-residue-B> <end-residue-B> "<<endl);
        String      chainA = parameterStringClass.getString(1);
        String      chainB = parameterStringClass.getString(4);
        ResidueID startResidueA = myBiopolymerClassContainer.residueID(userVariables,parameterStringClass.getString(2),chainA );
        ResidueID   endResidueA = myBiopolymerClassContainer.residueID(userVariables,parameterStringClass.getString(3),chainA );
        ResidueID startResidueB = myBiopolymerClassContainer.residueID(userVariables,parameterStringClass.getString(5),chainB );
        ResidueID   endResidueB = myBiopolymerClassContainer.residueID(userVariables,parameterStringClass.getString(6),chainB );
        int threadedLength = myBiopolymerClassContainer.updBiopolymerClass(chainA).difference (endResidueA , startResidueA) + 1;
        SimTK_ERRCHK_ALWAYS(
                (( startResidueA <= endResidueA      )),
                "[ParameterReader.cpp]",  "In the threading command, the end residue must be greater than or equal to the start residue for each chain.");
        SimTK_ERRCHK_ALWAYS(
                (( startResidueB <= endResidueB      )),
                "[ParameterReader.cpp]",  "In the proteinThreading command, the end residue must be greater than or equal to the start residue for each chain.");
        SimTK_ERRCHK_ALWAYS(
                (myBiopolymerClassContainer.updBiopolymerClass(chainA).difference( endResidueA,startResidueA) ==myBiopolymerClassContainer.updBiopolymerClass(chainB).difference( endResidueB,startResidueB)),
                "[ParameterReader.cpp]",  "In the proteinThreading command, the two threaded segments must be of the same length.");

          
        if (((parameterStringClass.getString(0)).compare("coupleMobilizers")      == 0) ) 
        {
            ;
        }
        else if (((parameterStringClass.getString(0)).compare("couplePsiPhiAngles")      == 0)) 
        {
            for (int i = 0; i < threadedLength; i++) 
            {
                ResidueInfo myResidueInfoA = myBiopolymerClassContainer.updBiopolymerClass(chainA).updResidueInfo(myBiopolymerClassContainer.updBiopolymerClass(chainA).sum( startResidueA , i)) ;

                ResidueInfo myResidueInfoB = myBiopolymerClassContainer.updBiopolymerClass(chainB).updResidueInfo(   myBiopolymerClassContainer.updBiopolymerClass(chainB).sum(startResidueB , i)) ;
                for (int j = 0; j < 2; j++) 
                {
                    String atomNameA =  myResidueInfoA.getAtomName(ResidueInfo::AtomIndex (j));
                    if (j == 0) {
                        atomNameA = "CA";
                    } else if (j ==  1){ 
                        atomNameA = "C";
                    } 
                    else 
                    {
                        MMBLOG_FILE_FUNC_LINE(CRITICAL, "Something odd happened.   "<<endl);
                    }

                    if (myBiopolymerClassContainer.updBiopolymerClass(chainB).hasAtom( myBiopolymerClassContainer.updBiopolymerClass(chainB).sum (startResidueB , i),atomNameA)) 
                    {
                        myBiopolymerClassContainer.addConstraint(userVariables, 
                                                                 atomNameA,
                                                                 myBiopolymerClassContainer.updBiopolymerClass(chainA).sum(startResidueA , i).outString(), 
                                                                 chainA,
                                                                 atomNameA, 
                                                                 myBiopolymerClassContainer.updBiopolymerClass(chainB).sum(startResidueB,i).outString(),
                                                                 chainB,
                                                                 CoupledCoordinate,
                                                                 constraintToGroundContainer
                                                                );
                                
                        // follows this prototype:
                        //void BiopolymerClassContainer::addConstraint(map<const String,double> myUserVariables,
                        //     const String atomName1, const String inputResidueString1,const  String chain1, 
                        //     const String atomName2, const String inputResidueString2,const  String chain2, 
                        //     ConstraintType myConstraintType,
                        //     ConstraintToGroundContainer & constraintToGroundContainer)

                        }
                } // of for numatoms
            }        
        } else 
        {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "Something odd happened.  This is not a supported threading command. "<<endl);
        }

        return;
    };

    if ((parameterStringClass.getString(0)).compare("rootMobilizer") == 0) {
        MMBLOG_FILE_FUNC_LINE(ALWAYS,
                "Syntax: rootMobilizer <\"Free\" | \"Weld\"> "<<endl
                <<": ... if you want to apply this root mobilizer to all biopolymer chains in the system "<<endl
                <<":     or: rootMobilizer <Chain> <\"Free\" | \"Weld\"> "<<endl
                <<": ... if you want to apply this root mobilizer to the specified Chain.                "<<endl
                <<": \"Free\" is the default, and confers 6 degrees of freedom to the root atom.  \"Weld\" grants the root atom 0 DOFs -- this behaves much like welding the root atom to ground but is computationally MUCH cheaper and numerically more accurate (no constraint tolerance to worry about). "<<endl);
        if ((parameterStringClass.getString(1).length() == 0) ) { MMBLOG_FILE_FUNC_LINE(CRITICAL, "Not enough parameters for this command."<<endl);}
        if ((parameterStringClass.getString(2).length() == 0) ) {
            MMBLOG_FILE_FUNC_LINE(INFO, "You have specified that all chains will be joined to ground with a "<<parameterStringClass.getString(1)<<" mobilizer"<<endl);
            // This method will validate parameterStringClass.getString(1) :
            myBiopolymerClassContainer.setFirstResidueMobilizerType(parameterStringClass.getString(1));

        } else if ((parameterStringClass.getString(3).length() == 0) ) {
            MMBLOG_FILE_FUNC_LINE(INFO, "You have specified that chain "<<parameterStringClass.getString(1)<<" will be joined to ground with a Weld mobilizer"<<endl);
            myBiopolymerClassContainer.updBiopolymerClass(parameterStringClass.getString(1)).setFirstResidueMobilizerType(parameterStringClass.getString(2) );
        } else { MMBLOG_FILE_FUNC_LINE(CRITICAL, "Too many parameters for this command!"<<endl);}
        return;

    }

    if ((parameterStringClass.getString(0)).compare("constrainToGround") == 0) {
        MMBLOG_FILE_FUNC_LINE(ALWAYS, "Syntax: constrainToGround <chain> <ResidueID> "<<endl);
        //MMBLOG_FILE_FUNC_LINE(" : Or: constrainToGround  "<<endl;
        //MMBLOG_FILE_FUNC_LINE(" : In the latter case, all first residues will be have no degrees of freedom with respect to ground.  This is very efficient! It is not the same as constraining all first residues to ground.  "<<endl;

        if (removeRigidBodyMomentum) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "You must first set removeRigidBodyMomentum False before using this command! "<<endl);
        }
        if ((parameterStringClass.getString(1).length() == 0) ) {
            //firstResidueMobilizerType = "Weld";
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "You have tried to use an obsolete syntax.  Please use the rootMobilizer command."<<endl);

            //FIXME: Unreachable code ?!?!
            MMBLOG_FILE_FUNC_LINE(INFO, "You have specified that all chains will be joined to ground with a Weld mobilizer"<<endl);
            myBiopolymerClassContainer.setFirstResidueMobilizerType(String("Weld"));

        } else if ((parameterStringClass.getString(2).length() == 0) ) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "You have tried to use an obsolete syntax.  Please use the rootMobilizer command."<<endl);

            //FIXME: Unreachable code ?!?!
            MMBLOG_FILE_FUNC_LINE(INFO, "You have specified that chain "<<parameterStringClass.getString(2)<<" will be joined to ground with a Weld mobilizer"<<endl);
            myBiopolymerClassContainer.updBiopolymerClass(parameterStringClass.getString(1)).setFirstResidueMobilizerType(String("Weld"));
        }
        else {
            if ( myBiopolymerClassContainer.updBiopolymerClass(parameterStringClass.getString(1)).getFirstResidueID() == 
                    myBiopolymerClassContainer.residueID(userVariables, parameterStringClass.getString(2), parameterStringClass.getString(1)) )  
                        MMBLOG_FILE_FUNC_LINE(WARNING, "Looks like you're welding the first residue to ground.  Please consider using the rootMobilizer command -- it's more accurate and much cheaper."<<endl);
            // should figure out if parameterStringClass.getString(2) is the first residue number, and if so give it the same treatment as above.
            myBiopolymerClassContainer.addConstraintToGround(userVariables,  parameterStringClass.getString(2),parameterStringClass.getString(1),constraintToGroundContainer);
        }
        return;
    }

    if ((parameterStringClass.getString(0)).compare("addRingClosingBond") == 0) 
    {
        MMBLOG_FILE_FUNC_LINE(ALWAYS,
                "Syntax for biopolymers: addRingClosingBond  <chainID> <residueID1> <atomName1> <bondCenterName1>  <residueID2>  <atomName2> <bondCenterName2> "<<endl
                <<"Syntax for molecule objects: addRingClosingBond  <chainID> <atomName1> <bondCenterName1> <atomName2> <bondCenterName2> "<<endl);
	String chainID = parameterStringClass.getString(1); // I believe BiopolymerClassContainer will puke informatively if this chain doesn't exist.
	CovalentBondClass myBond;
        if (myBiopolymerClassContainer.hasChainID(chainID)){
                parameterStringClass.validateNumFields(8);
		ResidueID residueID1 = myBiopolymerClassContainer.residueID(userVariables,parameterStringClass.getString(2), chainID);
		ResidueID residueID2 = myBiopolymerClassContainer.residueID(userVariables,parameterStringClass.getString(5), chainID);
		String atomName1 = parameterStringClass.getString(3); // should we validate now?
		String atomName2 = parameterStringClass.getString(6); // should we validate now?
		myBiopolymerClassContainer.updBiopolymerClass(chainID).validateAtomPathName(myBiopolymerClassContainer.updBiopolymerClass(chainID).atomPathString(residueID1,atomName1));
		myBiopolymerClassContainer.updBiopolymerClass(chainID).validateAtomPathName(myBiopolymerClassContainer.updBiopolymerClass(chainID).atomPathString(residueID2,atomName2));
		String bondCenterName1 = parameterStringClass.getString(4);
		String bondCenterName2 = parameterStringClass.getString(7);
	     
		myBond.setResidueID1(residueID1 ) ;
		myBond.setResidueID2(residueID2 ) ;
		myBond.setAtomName1(atomName1); 
		myBond.setAtomName2(atomName2); 
		myBond.setChain1(chainID);
		myBond.setChain2(chainID);
		myBond.setBondCenterName1(bondCenterName1);
		myBond.setBondCenterName2(bondCenterName2);
		additionalCovalentBondVector.push_back(myBond);
        } else if (moleculeClassContainer.hasChainID(chainID)) {
		ResidueID residueID1 = ResidueID(-11111,'X');
		ResidueID residueID2 = ResidueID(-11111,'X');
                parameterStringClass.validateNumFields(6);
		String atomName1 = parameterStringClass.getString(2); // should we validate now?
		String atomName2 = parameterStringClass.getString(4); // should we validate now?
		String bondCenterName1 = parameterStringClass.getString(3);
		String bondCenterName2 = parameterStringClass.getString(5);
		myBond.setAtomName1(atomName1); 
		myBond.setAtomName2(atomName2); 
		myBond.setChain1(chainID);
		myBond.setChain2(chainID);
		myBond.setBondCenterName1(bondCenterName1);
		myBond.setBondCenterName2(bondCenterName2);
		additionalCovalentBondVector.push_back(myBond);
                        
        } else {
                MMBLOG_FILE_FUNC_LINE(CRITICAL, "Error! chain >"<<chainID<<"< is not a biopolymer or molecule! Or possibly not enough parameters specified."<<endl);
        }
        //addRingClosingBond (chainID, residueID1 , atomName1, bondCenterName1, residueID2, atomName2, bondCenterName2);
        return;
    }

    // command is problematic after createing ParameterString class, leaving out for now.
    if ( ((parameterStringClass.getString(0)).compare( "restraint") == 0) || ((parameterStringClass.getString(0)).compare("restrainToGround") == 0)  )    
    { // if this is a restraint or constraint
        BasePair        myBasePair;
        myBasePair.rotationCorrection1 = Rotation(0.0,UnitVec3(0,0,1));
        myBasePair.rotationCorrection2 = Rotation(0.0,UnitVec3(0,0,1));
        myBasePair.BasePairIsTwoTransformForce = parameterStringClass.getString(0);
        if (((parameterStringClass.getString(0)).compare("restrainToGround") == 0)  ) 
        {
            if (removeRigidBodyMomentum)
            {
                MMBLOG_FILE_FUNC_LINE(CRITICAL, "Please set \"removeRigidBodyMomentum FALSE\" before issuing the constrainToGround or restrainToGround commands!"<<endl);
            }

            SimTK_ERRCHK_ALWAYS(
                (!(compareUpper("TRUE",(parameterStringClass.getString(1)).c_str()) || compareUpper("FALSE",(parameterStringClass.getString(1)).c_str()))),
                "[ParameterReader.cpp]", "The weldToGround command has changed.  The new syntax is \"restrainToGround chain-ID residue-number stage\" where the middle two parameters are the chain ID and residue number which is to be welded to the ground.  The stage is that stage at which the constraint is first applied.");
            myBasePair.FirstBPChain = "GROUND"       ;            
            myBasePair.FirstBPResidue = ResidueID();                           
            myBasePair.FirstBPEdge = "Weld";             
            myBasePair.SecondBPChain = ((parameterStringClass.getString(1)));            
            myBasePair.SecondBPResidue = myBiopolymerClassContainer.residueID(userVariables,parameterStringClass.getString(2),myBasePair.SecondBPChain);
            if (!
                (
                 myMonoAtomsContainer.hasChainID(myBasePair.SecondBPChain) ||
                 myBiopolymerClassContainer.hasChainID(myBasePair.SecondBPChain)
                ))
            {
                MMBLOG_FILE_FUNC_LINE(CRITICAL, "You have not declared the chain "<<myBasePair.SecondBPChain<<" yet, or the chain is empty.  "<<endl);
            }
            myBasePair.SecondBPEdge = "Weld";                 
            myBasePair.OrientationBP = "Cis";                 
            myBasePair.BasePairPriority = myAtoI(userVariables,(parameterStringClass.getString(3)).c_str());
        } else 
        { 
            myBasePair.FirstBPChain = ((parameterStringClass.getString(1)));            
            myBasePair.FirstBPResidue = myBiopolymerClassContainer.residueID(userVariables,(parameterStringClass.getString(2)),myBasePair.FirstBPChain ); 
            myBasePair.FirstBPEdge = parameterStringClass.getString(3);
            myBasePair.SecondBPChain = ((parameterStringClass.getString(4)));            
            myBasePair.SecondBPResidue =  myBiopolymerClassContainer.residueID(userVariables,(parameterStringClass.getString(5)),myBasePair.SecondBPChain );            
            myBasePair.SecondBPEdge = parameterStringClass.getString(6);            
            myBasePair.OrientationBP = parameterStringClass.getString(7);           
            if (parameterStringClass.getString(8).length() > 0) {  
                MMBLOG_FILE_FUNC_LINE(CRITICAL, "Too many parameters!  Check your syntax.  This command no longer accepts a Stage parameter. "<<endl);
            }
            else { 
                MMBLOG_FILE_FUNC_LINE(INFO, "Stage of command automatically set to stage 1."<<endl);
                myBasePair.BasePairPriority = 1; 
            }
            if ((parameterStringClass.getString(9)).length() >0) 
            {
                MMBLOG_FILE_FUNC_LINE(CRITICAL, "Too many parameters!  Check your syntax.  This command no longer accepts a Stage parameter. "<<endl);
            } 
            else 
                myBasePair.BasePairTemporary = 0;
            if ((parameterStringClass.getString(0).compare("constraint") == 0) || (parameterStringClass.getString(0).compare("restraint") == 0)) 
                if (myBasePair.FirstBPEdge.compare("Weld")==0){ // all is well
                } else {
                    MMBLOG_FILE_FUNC_LINE(CRITICAL, "The only type of "<<parameterStringClass.getString(0)<<" permitted is Weld.  You have specified : "<< myBasePair.FirstBPEdge<<endl);
                }
        }
        if (safeParameters) 
        {

            if ((parameterStringClass.getString(0).compare("mobilizer") == 0)) 
            {
                SimTK_ERRCHK_ALWAYS(
                        (myBasePair.FirstBPResidue <= myBasePair.SecondBPResidue  ), 
                        "[ParameterReader.cpp]",                    "For contact and mobilizer commands, please specify the lower numbered residue first.  ");
            }
            if ((parameterStringClass.getString(0).compare("constraint") == 0) || (parameterStringClass.getString(0).compare("mobilizer") == 0)  ) {
                if (parameterStringClass.getString(5).length() == 0) {
                    MMBLOG_FILE_FUNC_LINE(CRITICAL, "You have not specified enough parameters for this command. "<<endl);
                };
                SimTK_ERRCHK1_ALWAYS(
                        (parameterStringClass.getString(6).size() == 0) ,
                        "[ParameterReader.cpp]","For constraint, contact, and mobilizer commands, it is redundant to specify a second type, you have specified : %s.    ", parameterStringClass.getString(6).c_str());
            };                   
        }

        if (compareUpper("GROUND",((myBasePair).FirstBPChain).c_str())) 
        {
        }
        else {

            if (!
                    (
                     myMonoAtomsContainer.hasChainID(myBasePair.FirstBPChain) ||
                     myBiopolymerClassContainer.hasChainID(myBasePair.FirstBPChain)
                    ))
            {MMBLOG_FILE_FUNC_LINE(CRITICAL, "You have not declared the chain "<<myBasePair.FirstBPChain<<" yet, or the chain is empty.  "<<endl);}

            if (!
                    (
                     myMonoAtomsContainer.hasChainID(myBasePair.SecondBPChain) ||
                     myBiopolymerClassContainer.hasChainID(myBasePair.SecondBPChain)
                    ))
            {MMBLOG_FILE_FUNC_LINE(CRITICAL, "You have not declared the chain "<<myBasePair.SecondBPChain<<" yet, or the chain is empty.  "<<endl);}
        }

        r++; 
        MMBLOG_FILE_FUNC_LINE(DEBUG, " r incremented to "<<r<<endl);
        myBasePair.basePairSatisfied = "False"; //initialize
        baseOperationVector.push_back(myBasePair); 
        return;
    } //End restraint

    if  (((parameterStringClass.getString(0)).compare("contact") == 0))  {
        MMBLOG_FILE_FUNC_LINE(ALWAYS,
                "Current syntax of the \"contact\" command is:  contact <contact scheme> <chain> <start residue> <end residue> "<<endl
                <<"Or:  contact <contact scheme> <chain> <residue>  "<<endl
                <<"Or:  contact <contact scheme> <chain>  "<<endl);

        if (parameterStringClass.getString(5).length() >0) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "You have specified too many parameters for this command. "<<endl);
        } else if (parameterStringClass.getString(4).length() >0 ) {
            ContactStretch myContactStretch;
            myContactStretch.ContactScheme = parameterStringClass.getString(1);
            myContactStretch.setChain ( parameterStringClass.getString(2));
            myContactStretch.setStartResidue (myBiopolymerClassContainer.residueID(userVariables,parameterStringClass.getString(3),myContactStretch.getChain()));
            myContactStretch.setEndResidue ( myBiopolymerClassContainer.residueID(userVariables,parameterStringClass.getString(4),myContactStretch.getChain() ));
            //myBiopolymerClassContainer.updBiopolymerClass(myContactStretch.getChain()).validateResidueID(myContactStretch.StartResidue);
            //myBiopolymerClassContainer.updBiopolymerClass(myContactStretch.getChain()).validateResidueID(myContactStretch.getEndResidue());
            contactContainer.addContactToVector(myContactStretch,myBiopolymerClassContainer);
        } else if (parameterStringClass.getString(3).length() >0 ) { // User has specified a range of residue numbers
            ContactStretch myContactStretch;
            myContactStretch.ContactScheme = parameterStringClass.getString(1);
            myContactStretch.setChain ( parameterStringClass.getString(2));
            myContactStretch.setStartResidue (myBiopolymerClassContainer.residueID(userVariables,parameterStringClass.getString(3),myContactStretch.getChain()));
            myContactStretch.setEndResidue ( myBiopolymerClassContainer.residueID(userVariables,parameterStringClass.getString(3),myContactStretch.getChain() ));
        } else if ((parameterStringClass.getString(2).length() >0 ))   
        {
            ContactStretch myContactStretch;
            myContactStretch.ContactScheme = parameterStringClass.getString(1);
            myContactStretch.setChain ( parameterStringClass.getString(2));
            myContactStretch.setStartResidue (myBiopolymerClassContainer.residueID(userVariables, String("FirstResidue"), myContactStretch.getChain()));
            myContactStretch.setEndResidue (myBiopolymerClassContainer.residueID(userVariables, String("LastResidue"), myContactStretch.getChain()));
            //myBiopolymerClassContainer.updBiopolymerClass(myContactStretch.getChain()).validateResidueID(myContactStretch.StartResidue);
            //myBiopolymerClassContainer.updBiopolymerClass(myContactStretch.getChain()).validateResidueID(myContactStretch.getEndResidue());
            contactContainer.addContactToVector(myContactStretch,myBiopolymerClassContainer);

        } else {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "You have specified the wrong number of parameters for this command. "<<endl);
        }
        return;
    } //End  contact

    if  (((parameterStringClass.getString(0)).compare("fitToDensity") == 0))  {

        MMBLOG_FILE_FUNC_LINE(ALWAYS,
                "Current syntax of  command is:  fitToDensity  <chain> "<<endl
                <<"Or :  fitToDensity  <chain> <first residue> <last residue> "<<endl);

        if ((dutyCycle<1.0) ) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "you can't set dutyCycle<1.0 if you also issue the fitToDensity command"<<endl);
        }
        if (myDensityMap.getDensityFileName() == "densityFileName-NOT-SET") {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "you can't issue fitToDensity before setting the densityFileName parameter, Currently densityFileName is set to: "<<myDensityMap.getDensityFileName()<<". Perhaps you issued 'density clear' AFTER 'density densityFileName' ? If so flip the order. "<<std::endl);
        }
        DensityStretch myDensityStretch;
        myDensityStretch.setChain ( parameterStringClass.getString(1));

        if (parameterStringClass.getString(1).length() == 0) {
            MMBLOG_FILE_FUNC_LINE(INFO, "You have requested that all biopolymer chains be fitted to the density map. "<<endl);
            if (myBiopolymerClassContainer.getNumBiopolymers() == 0) {
                MMBLOG_FILE_FUNC_LINE(CRITICAL, "This command can only be issued AFTER you have declared or instantiated your biopolymer chains."<<endl);
            }

            densityContainer.stuffDensityStretchVector(myBiopolymerClassContainer);
        }
        else if (parameterStringClass.getString(2).length() == 0) {
            if (parameterStringClass.getString(1).length() == 0) {MMBLOG_FILE_FUNC_LINE(CRITICAL, "You have not specified enough parameters for this command. "<<endl);}
	    if (myBiopolymerClassContainer.hasChainID(parameterStringClass.getString(1))) {
                myDensityStretch.setStartResidue ( myBiopolymerClassContainer.updBiopolymerClass(myDensityStretch.getChain()).getFirstResidueID());
                myDensityStretch.setEndResidue ( myBiopolymerClassContainer.updBiopolymerClass(myDensityStretch.getChain()).getLastResidueID());
                myBiopolymerClassContainer.updBiopolymerClass(myDensityStretch.getChain()).validateResidueID(myDensityStretch.getStartResidue());
                myBiopolymerClassContainer.updBiopolymerClass(myDensityStretch.getChain()).validateResidueID(myDensityStretch.getEndResidue());
                densityContainer.add(myDensityStretch,myBiopolymerClassContainer );}
	    else if (myMonoAtomsContainer.hasChainID(parameterStringClass.getString(1))) {
                myDensityStretch.setStartResidue ( myMonoAtomsContainer.getMonoAtoms(myDensityStretch.getChain()).getFirstResidueID());
                myDensityStretch.setEndResidue   ( myMonoAtomsContainer.getMonoAtoms(myDensityStretch.getChain()).getLastResidueID());
                densityContainer.addStretch(myDensityStretch );//}
	    }

        } else if (parameterStringClass.getString(4).length() == 0 ) {
            if (parameterStringClass.getString(3).length() == 0) {MMBLOG_FILE_FUNC_LINE(CRITICAL, "You have not specified enough parameters for this command. "<<endl);}
            myDensityStretch.setStartResidue ( myBiopolymerClassContainer.residueID(userVariables, parameterStringClass.getString(2), myDensityStretch.getChain()));
            myDensityStretch.setEndResidue   ( myBiopolymerClassContainer.residueID(userVariables, parameterStringClass.getString(3), myDensityStretch.getChain()));
            myBiopolymerClassContainer.updBiopolymerClass(myDensityStretch.getChain()).validateResidueID(myDensityStretch.getStartResidue() );
            myBiopolymerClassContainer.updBiopolymerClass(myDensityStretch.getChain()).validateResidueID(myDensityStretch.getEndResidue());
            densityContainer.add(myDensityStretch,myBiopolymerClassContainer );

        } else {MMBLOG_FILE_FUNC_LINE(CRITICAL, "You have specified too many parameters for this command."<<endl); } // field 5 is not empty
        MMBLOG_FILE_FUNC_LINE(INFO, "During reading of input file, printing updated list of density stretches: "<<endl);
        densityContainer.printDensityStretches    (); 
        return;
    } //End fitToDensity

    if  (((parameterStringClass.getString(0)).compare("fitElectrostaticDensity") == 0))  {

        MMBLOG_FILE_FUNC_LINE(ALWAYS,
                "Current syntax of  command is:  fitElectrostaticDensity  <chain> "<<endl
                <<"Or :  fitElectrostaticDensity  <chain> <first residue> <last residue> "<<endl);

        if ((dutyCycle<1.0) ) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "you can't set dutyCycle<1.0 if you also issue the fitElectrostaticDensity command"<<endl);
        }

        DensityStretch myDensityStretch;
        myDensityStretch.setChain ( parameterStringClass.getString(1));

        if (parameterStringClass.getString(1).length() == 0) {
            MMBLOG_FILE_FUNC_LINE(INFO, "You have requested that all biopolymer chains be fitted to the density map. "<<endl);
            if (myBiopolymerClassContainer.getNumBiopolymers() == 0) {
                MMBLOG_FILE_FUNC_LINE(CRITICAL, "This command can only be issued AFTER you have declared or instantiated your biopolymer chains."<<endl);
            }

            electroDensityContainer.stuffDensityStretchVector(myBiopolymerClassContainer);
        }
        else if (parameterStringClass.getString(2).length() == 0) {
            if (parameterStringClass.getString(1).length() == 0) {MMBLOG_FILE_FUNC_LINE(CRITICAL, "You have not specified enough parameters for this command. "<<endl);}
            myDensityStretch.setStartResidue ( myBiopolymerClassContainer.updBiopolymerClass(myDensityStretch.getChain()).getFirstResidueID());
            myDensityStretch.setEndResidue ( myBiopolymerClassContainer.updBiopolymerClass(myDensityStretch.getChain()).getLastResidueID());
            myBiopolymerClassContainer.updBiopolymerClass(myDensityStretch.getChain()).validateResidueID(myDensityStretch.getStartResidue());
            myBiopolymerClassContainer.updBiopolymerClass(myDensityStretch.getChain()).validateResidueID(myDensityStretch.getEndResidue());
            electroDensityContainer.add(myDensityStretch,myBiopolymerClassContainer );
        } else if (parameterStringClass.getString(4).length() == 0 ) {
            if (parameterStringClass.getString(3).length() == 0) {MMBLOG_FILE_FUNC_LINE(CRITICAL, "You have not specified enough parameters for this command. "<<endl);}
            myDensityStretch.setStartResidue ( myBiopolymerClassContainer.residueID(userVariables, parameterStringClass.getString(2), myDensityStretch.getChain()));
            myDensityStretch.setEndResidue   ( myBiopolymerClassContainer.residueID(userVariables, parameterStringClass.getString(3), myDensityStretch.getChain()));
            myBiopolymerClassContainer.updBiopolymerClass(myDensityStretch.getChain()).validateResidueID(myDensityStretch.getStartResidue() );
            myBiopolymerClassContainer.updBiopolymerClass(myDensityStretch.getChain()).validateResidueID(myDensityStretch.getEndResidue());
            electroDensityContainer.add(myDensityStretch,myBiopolymerClassContainer );

        } else {MMBLOG_FILE_FUNC_LINE(CRITICAL, "You have specified too many parameters for this command."<<endl); } // field 5 is not empty
        MMBLOG_FILE_FUNC_LINE(INFO, "During reading of input file, printing updated list of density stretches: "<<endl);
        electroDensityContainer.printDensityStretches    (); 
        return;
    } //End fitElectrostaticDensity

    if  (((parameterStringClass.getString(0)).compare("applyContactsWithin") == 0))  {
        #ifndef USE_OPENMM
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "You have compiled without USE_OPENMM ."<<endl);
        #endif
        MMBLOG_FILE_FUNC_LINE(ALWAYS,
                "Current syntax of  command is:  applyContactsWithin <radius (nm)> <contact scheme> <chain> <residue>  "<<endl
                <<"Note that in MMB 2.10 and earlier, we took the radius in Å.  We are going back to nm for consistency, with apologies for the confusion."<<endl);
        ContactWithin contactWithin;
        contactWithin.Chain = parameterStringClass.getString(3);
        contactWithin.Residue = myBiopolymerClassContainer.residueID(userVariables,parameterStringClass.getString(4), contactWithin.Chain);
        contactWithin.Radius =  myAtoF(userVariables,parameterStringClass.getString(1).c_str());
        contactWithin.ContactScheme = parameterStringClass.getString(2);
        contactContainer.pushContactWithin(contactWithin,myBiopolymerClassContainer);

        return;
    }
    if  (((parameterStringClass.getString(0)).compare("applyMobilizersWithin") == 0))  {
        MMBLOG_FILE_FUNC_LINE(ALWAYS,
                "Current syntax of  command is:  applyMobilizersWithin <Bond Mobility> <radius (nm)>  <chain> <residue>  "<<endl
                <<"Note that in MMB 2.10 and earlier, we took the radius in Å.  We are going back to nm for consistency, with apologies for the confusion."<<endl);
        MobilizerWithin mobilizerWithin;
        mobilizerWithin.setChain ( parameterStringClass.getString(3));
        mobilizerWithin.setResidue ( myBiopolymerClassContainer.residueID(userVariables,parameterStringClass.getString(4), mobilizerWithin.getChain()));
        mobilizerWithin.setRadius (  myAtoF(userVariables,parameterStringClass.getString(2).c_str()));
        mobilizerWithin.setBondMobilityString ( parameterStringClass.getString(1));
        mobilizerContainer.pushMobilizerWithin(mobilizerWithin,myBiopolymerClassContainer);
        return;
    }

    if  (((parameterStringClass.getString(0)).compare("mobilizeInterfaces") == 0))  {
        MMBLOG_FILE_FUNC_LINE(ALWAYS,
                "Syntax is:  mobilizeInterfaces <depth (nm)> <bond mobility> <chain>  "<<endl
                <<"This sets bond mobility to <bond mobility>  for all residues within <depth> of <chain> in all other chains.  It also sets the bond mobility for residues in <chain> that are within <depth> of that other chain -- so the interface is treated symmetrically."<<endl
                <<"You can also specify a multiple-chain set, and find all interfaces between that set and the remainder of chains in the system:"<<endl);

        MMBLOG_FILE_FUNC_LINE(ALWAYS,
                "Syntax is:  mobilizeInterfaces <depth (nm)> <bond mobility> <chain 1> [<chain 2> [<chain 3> [ ...etc]]]  "<<endl
                <<"Lastly, you can specify TWO multiple-chain sets, and find only the all interface between those two sets.  This is particularly useful if you have chains (e.g. threading templates) in your system which need to be ignored:"<<endl);
    
        MMBLOG_FILE_FUNC_LINE(ALWAYS,
                "Syntax is:  mobilizeInterfaces <depth (nm)> <bond mobility> <chain 1> [<chain 2> [<chain 3> [ ...etc]]]  Versus <chain 1> [<chain 2> [<chain 3> [ ...etc]]]  "<<endl);

        if (parameterStringClass.getString(3).length() == 0) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "Not enough parameters for this command! You must specify at least one chain."<<endl);
        }
        double depth = myAtoF(userVariables,parameterStringClass.getString(1).c_str());
        String bondMobilityString = parameterStringClass.getString(2);
        vector<String> chains;        chains.clear();
        vector<String> partnerChains; partnerChains.clear();
        int i = 3;
        String versusString = String("Versus").toUpper();
        while ((parameterStringClass.getString(i).length() >0)  &&
               (parameterStringClass.getString(i).toUpper().compare(versusString ) != 0 )) { // if parameter i is "versus", stop and start reading partnerChains
            myBiopolymerClassContainer.validateChainID(parameterStringClass.getString(i)); // Make sure this is a valid chain
            chains.push_back( parameterStringClass.getString(i));
            MMBLOG_FILE_FUNC_LINE(INFO, "Added chain "<<parameterStringClass.getString(i)<<" to first set "<<endl);
            i++;
        }
        if (parameterStringClass.getString(i).toUpper().compare(versusString) == 0 ) { 
            i++;
            MMBLOG_FILE_FUNC_LINE(INFO, "Keyword "<<versusString<<" detected."<<endl);
            if (parameterStringClass.getString(i).length() == 0) {
		MMBLOG_FILE_FUNC_LINE(CRITICAL, "The \"versus\" keyword must be followed by at least one chain ID."<<endl);
            }
        }
        while (parameterStringClass.getString(i).length() >0) {
            myBiopolymerClassContainer.validateChainID(parameterStringClass.getString(i)); // Make sure this is a valid chain
            partnerChains.push_back( parameterStringClass.getString(i));
            MMBLOG_FILE_FUNC_LINE(INFO, "Added chain "<<parameterStringClass.getString(i)<<" to second set "<<endl);
            i++;
        }
        
        for (size_t n = 0; n < chains.size(); n++) for (size_t m = 0; m < partnerChains.size(); m++) if (chains[n].compare(partnerChains[m])==0) {
	    MMBLOG_FILE_FUNC_LINE(CRITICAL, "The reference chain "<<chains[n]<<" is the same as partner chain "<<partnerChains[m]<<". This is not kosher!"<<endl);
        }
        
        mobilizerContainer.interfaceContainer.addInterface(chains, partnerChains, depth, bondMobilityString);

        return;
    }

    if(((parameterStringClass.getString(0)).compare("mobilizeDomainsInterface") == 0))
    {
        MMBLOG_FILE_FUNC_LINE(ALWAYS,
                "Syntax is:  mobilizeDomainsInterface <depth (nm)> <bond mobility> <chain domain 1> <resStart1> <resEnd1> <chain domain 2> <resStart2> <resEnd2> [<'backbone-rigid'>] "<<endl);
        if (parameterStringClass.getString(8).length() == 0) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "Not enough parameters for this command!"<<endl);
        }
        double range = myAtoF(userVariables,parameterStringClass.getString(1).c_str());
        String bondMobilityString = parameterStringClass.getString(2);
        String chain1 = parameterStringClass.getString(3);
        ResidueID resStart1 = myBiopolymerClassContainer.residueID(userVariables,parameterStringClass.getString(4), chain1);
        ResidueID resEnd1   = myBiopolymerClassContainer.residueID(userVariables,parameterStringClass.getString(5), chain1);
        String chain2 = parameterStringClass.getString(6);
        ResidueID resStart2 = myBiopolymerClassContainer.residueID(userVariables,parameterStringClass.getString(7), chain2);
        ResidueID resEnd2   = myBiopolymerClassContainer.residueID(userVariables,parameterStringClass.getString(8), chain2);
        bool rigidBackbone = !parameterStringClass.getString(9).compare("backbone-rigid");

        MobilizerDomainsInterface mDI;
        mDI.domain1 = ResidueStretch(chain1, resStart1, resEnd1);
        mDI.domain2 = ResidueStretch(chain2, resStart2, resEnd2);
        mDI.range = range;
        mDI.MobilizerString = bondMobilityString;
        mDI.rigidBackbone = rigidBackbone;
        mobilizerDomainsInterfaceVector.push_back(mDI);
        return;       
    }
    if  (((parameterStringClass.getString(0)).compare("readInQVector") == 0))  {
        parameterStringClass.validateNumFields(2);
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "The readInQVector is obsolete."<<endl);
        readInQVector = myAtoI(userVariables,(parameterStringClass.getString(1)).c_str());     
        readPreviousFrameFile = (!(readInQVector));
        MMBLOG_FILE_FUNC_LINE(DEBUG, "readInQVector ="<<readInQVector<<endl);
        return;
    }
    if  (((parameterStringClass.getString(0)).compare("readMagnesiumPositionsFromFile") == 0))  {
        parameterStringClass.validateNumFields(2);
        readMagnesiumPositionsFromFile = myAtoI(userVariables,(parameterStringClass.getString(1)).c_str());     
        return;
    }
    if  (((parameterStringClass.getString(0)).compare("-FS") == 0) || ((parameterStringClass.getString(0)).compare("firstStage") == 0))  {
        parameterStringClass.validateNumFields(2);
        firstStage = myAtoI(userVariables,(parameterStringClass.getString(1)).c_str());     

        if (readStage || readAtOneStageOnly ||readOnlyUntilStage || readExcept ) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "You are not allowed to set the firstStage or lastStage parameter in a conditional block."<<endl);
        }
        MMBLOG_FILE_FUNC_LINE(DEBUG, " firstStage ="<<firstStage<<endl);
        return;
    }
    if  (((parameterStringClass.getString(0)).compare( "periodicallyUpdateParameters") == 0))  {
        parameterStringClass.validateNumFields(2);
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "The periodicallyUpdateParameters is obsolete."<<endl);
        periodicallyUpdateParameters = myAtoI(userVariables,(parameterStringClass.getString(1)).c_str());     
        MMBLOG_FILE_FUNC_LINE(DEBUG, "periodicallyUpdateParameters ="<<periodicallyUpdateParameters<<endl);
        return;
    }
    if  (((parameterStringClass.getString(0)).compare( "lastStage") == 0))  {
        parameterStringClass.validateNumFields(2);
        lastStage = myAtoI(userVariables,(parameterStringClass.getString(1)).c_str());     

        if (readStage || readAtOneStageOnly ||readOnlyUntilStage || readExcept ) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "You are not allowed to set the firstStage or lastStage parameter in a conditional block."<<endl);
        }
        MMBLOG_FILE_FUNC_LINE(DEBUG, "lastStage ="<< lastStage<<endl);
        return;
    }
    if  (((parameterStringClass.getString(0)).compare( "physicsRadius"  ) == 0))  {
        MMBLOG_FILE_FUNC_LINE(ALWAYS,
                "If >0, MMB will include all residues within the provided radius of 'flexible' atoms, in the physics zone. 'flexible' is defined as atoms belonging to low-mass (<40 au) groups."<<endl
                <<"Distance is determined using OpenMM's neighborlisting function, so is on basis of nearest atoms, and is quite efficient"<<endl);
        MMBLOG_FILE_FUNC_LINE(ALWAYS, "Syntax:  physicsRadius <radius>"<<endl);
        parameterStringClass.validateNumFields(2);
        physicsRadius = myAtoF(userVariables,(parameterStringClass.getString(1)));     
        if (physicsRadius < 0.0) {MMBLOG_FILE_FUNC_LINE(CRITICAL, "physicsRadius must be >0 !"<<endl);}
        MMBLOG_FILE_FUNC_LINE(DEBUG, "physicsRadius ="<< physicsRadius<<endl);
        return;
    }

    if(parameterStringClass.getString(0).compare("deactivatePhysics") == 0) 
    {
        parameterStringClass.validateNumFields(2);
        String chain = parameterStringClass.getString(1);
        BiopolymerClass & bpc = myBiopolymerClassContainer.updBiopolymerClass(chain);
        MMBLOG_FILE_FUNC_LINE(INFO, "will not be included in any physics or steric interactions." << endl);
        bpc.setActivePhysics(false);
        return;
    }
    if(parameterStringClass.getString(0).compare("activatePhysics") == 0) 
    {
        parameterStringClass.validateNumFields(2);
        String chain = parameterStringClass.getString(1);
        BiopolymerClass & bpc = myBiopolymerClassContainer.updBiopolymerClass(chain);
        MMBLOG_FILE_FUNC_LINE(INFO, " Chain " << chain << " will now be included in requested physics or steric interactions." << endl);
        bpc.setActivePhysics(true);
        return;
    }

    if  (((parameterStringClass.getString(0)).compare( "piecewiseRigidify") == 0))  {
        parameterStringClass.validateNumFields(2);
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "This is no longer a supported parameter.   "<<endl);
        piecewiseRigidify =  aToBool(parameterStringClass.getString(0),(parameterStringClass.getString(1)).c_str());     
        MMBLOG_FILE_FUNC_LINE(DEBUG, "piecewiseRigidify ="<< piecewiseRigidify<<endl);
        return;
    }

    if  (((parameterStringClass.getString(0)).compare( "planarityThreshold") == 0))  {
        MMBLOG_FILE_FUNC_LINE(ALWAYS, "This parameter sets the threshold for considering the chirality to be of type Planar.  Units: radians."<<endl);
        parameterStringClass.validateNumFields(2);
        planarityThreshold = myAtoF(userVariables,(parameterStringClass.getString(1)));     
        if (planarityThreshold <= 0) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "This angle should be positive."<<endl);
        }

        if (planarityThreshold >45 ) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "Doesn't make sense for this angle to be so large."<<endl);
        }
        MMBLOG_FILE_FUNC_LINE(DEBUG, "planarityThreshold="<<planarityThreshold<<endl);
        return;
    }
    if  (((parameterStringClass.getString(0)).compare( "potentialType"    ) == 0))  {
        parameterStringClass.validateNumFields(2);
        potentialType     = ((parameterStringClass.getString(1)));     
        MMBLOG_FILE_FUNC_LINE(DEBUG, "  potentialType     ="<< potentialType    <<endl);
        return;
    }
    if  (((parameterStringClass.getString(0)).compare( "hardSphereStiffnessMultiplier") == 0))  {
        parameterStringClass.validateNumFields(2);
        hardSphereStiffnessMultiplier = myAtoF(userVariables,(parameterStringClass.getString(1)).c_str());     
        return;
    }
    if  (((parameterStringClass.getString(0)).compare( "magnesiumIonRadius") == 0))  {
        parameterStringClass.validateNumFields(2);
        magnesiumIonRadius = myAtoF(userVariables,(parameterStringClass.getString(1)).c_str());     
        MMBLOG_FILE_FUNC_LINE(DEBUG, "  magnesiumIonRadius ="<< magnesiumIonRadius<<endl);
        return;
    }
    if  (((parameterStringClass.getString(0)).compare( "matchDefaultSkipTopLevelTransform") == 0))  {
        parameterStringClass.validateNumFields(2);
        matchDefaultSkipTopLevelTransform = myAtoI(userVariables,(parameterStringClass.getString(1)).c_str());     
        MMBLOG_FILE_FUNC_LINE(DEBUG, "  matchDefaultSkipTopLevelTransform ="<< matchDefaultSkipTopLevelTransform<<endl);
        return;
    }
    if  (((parameterStringClass.getString(0)).compare( "numMagnesiumIons") == 0))  {
        parameterStringClass.validateNumFields(2);
        MMBLOG_FILE_FUNC_LINE(CRITICAL, " The numMagnesiumIons parameter is obsolete.  Please use the new monoAtoms command."<<endl);
        //FIXME: Unreachable code ?!?!?
        numMagnesiumIons = myAtoI(userVariables,(parameterStringClass.getString(1)).c_str());     
        MMBLOG_FILE_FUNC_LINE(DEBUG, "  numMagnesiumIons ="<< numMagnesiumIons<<endl);
        return;
    }
    if  (((parameterStringClass.getString(0)).compare( "safeParameters") == 0))  {
        parameterStringClass.validateNumFields(2);
        safeParameters = aToBool(parameterStringClass.getString(0),(parameterStringClass.getString(1)).c_str());     
        MMBLOG_FILE_FUNC_LINE(DEBUG, "  safeParameters ="<<safeParameters <<endl);
        return;
    }
    if  (((parameterStringClass.getString(0)).compare( "setDefaultMDParameters") == 0))  {
        parameterStringClass.validateNumFields(1);
        addBackboneOxygenForces =0;
        setChiBondMobility =0;
        globalBondTorsionScaleFactor   = 1.0 ;
        globalAmberImproperTorsionScaleFactor = 1.0;
        globalBondBendScaleFactor      = 1.0 ;
        globalBondStretchScaleFactor   = 1.0 ;
        globalBondTorsionScaleFactor   = 1.0 ;
        globalCoulombScaleFactor       = 1.0 ;
        globalVdwScaleFactor           = 1.0 ;
        globalAmberImproperTorsionScaleFactor = 1.0;
        weldToGround = false;
        wkdpGlobalBondTorsionScaleFactor = 0;
        setOverallBondMobility = 0;
        gcsfi++;
        gvsfi++;
        //}   
        return;
    }
    if  (((parameterStringClass.getString(0)).compare( "setOverallBondMobility") == 0))  {
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "This is no longer a supported parameter.   "<<endl);
        // FIXME: Unreachable code ?!?!?!
        setOverallBondMobility = aToBool(parameterStringClass.getString(0),(parameterStringClass.getString(1)).c_str());     
        MMBLOG_FILE_FUNC_LINE(INFO, "  setOverallBondMobility ="<< setOverallBondMobility<<endl);
        return;
    }
    if  (((parameterStringClass.getString(0)).compare( "integratorAccuracy") == 0) ) {
        parameterStringClass.validateNumFields(2);
        integratorAccuracy = myAtoF(userVariables,(parameterStringClass.getString(1)).c_str());     
        return;
    }
    if  (((parameterStringClass.getString(0)).compare( "initialSeparation" ) == 0) ) {
        parameterStringClass.validateNumFields(2);
        initialSeparation  = myAtoF(userVariables,(parameterStringClass.getString(1)).c_str());     
        return;
    }
    if  (((parameterStringClass.getString(0)).compare( "integratorStepSize") == 0) || ((parameterStringClass.getString(0)).compare( "stepSize") == 0))  {
        parameterStringClass.validateNumFields(2);
        integratorStepSize = myAtoF(userVariables,(parameterStringClass.getString(1)).c_str());     
        return;
    }
    if  (((parameterStringClass.getString(0)).compare( "kbBackboneTorsionGlobalScaleFactor") == 0))  {
        parameterStringClass.validateNumFields(2);
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "The kbBackboneTorsionGlobalScaleFactor obsolete."<<endl);
        kbBackboneTorsionGlobalScaleFactor = myAtoF(userVariables,(parameterStringClass.getString(1)).c_str());     
        if (safeParameters) SimTK_ERRCHK1_ALWAYS( ((kbBackboneTorsionGlobalScaleFactor ) >= .0  )    ,"[ParameterReader.cpp]","You have selected kbBackboneTorsionGlobalScaleFactor =  %f, whereas this value should be >=0.",kbBackboneTorsionGlobalScaleFactor);

        return;
    }
    if  (((parameterStringClass.getString(0)).compare( "twoTransformForceMultiplier") == 0) ||((parameterStringClass.getString(0)).compare( "forceMultiplier") == 0) || ((parameterStringClass.getString(0)).compare( "baseInteractionScaleFactor") == 0) || ((parameterStringClass.getString(0)).compare( "baseInteractionForceMultiplier") == 0) )  {
        parameterStringClass.validateNumFields(2);
        twoTransformForceMultiplier = myAtoF(userVariables,(parameterStringClass.getString(1)).c_str()); 
        return;
    }
    if  (((parameterStringClass.getString(0)).compare( "tinkerParameterFileName") == 0))  {
        parameterStringClass.validateNumFields(2);
        if (myBiopolymerClassContainer.getNumBiopolymers() > 0){
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "The tinkerParameterFileName parameter can only be set prior to instantiating the first biopolymer. "<<endl);
        } 
        tinkerParameterFileName =((parameterStringClass.getString(1)));     
        MMBLOG_FILE_FUNC_LINE(DEBUG, " tinkerParameterFileName ="<<tinkerParameterFileName<<endl);
        loadTinkerParameterFile = 1   ;
        MMBLOG_FILE_FUNC_LINE(INFO, "Since you specified tinkerParameterFileName, we are setting loadTinkerParameterFile = "<<loadTinkerParameterFile<<endl);

        return;
    }
    if  (((parameterStringClass.getString(0)).compare( "integratorType"         ) == 0))  {
        parameterStringClass.validateNumFields(2);
        integratorType          =((parameterStringClass.getString(1)));     
        MMBLOG_FILE_FUNC_LINE(DEBUG, " integratorType          ="<<integratorType         <<endl);

        return;
    }
    if  (((parameterStringClass.getString(0)).compare("readPreviousFrameFile") == 0))  {
        parameterStringClass.validateNumFields(2);
        //readPreviousFrameFile = myAtoF(userVariables,(parameterStringClass.getString(1)).c_str());
        readPreviousFrameFile = aToBool(parameterStringClass.getString(0), (parameterStringClass.getString(1)).c_str());      
        MMBLOG_FILE_FUNC_LINE(DEBUG, "  readPreviousFrameFile ="<< readPreviousFrameFile<<endl);
        return;
    }
    if  (((parameterStringClass.getString(0)).compare("-T")  == 0) || ((parameterStringClass.getString(0)).compare("temperature")  == 0))   {
        MMBLOG_FILE_FUNC_LINE(ALWAYS, "syntax: temperature <temperature>"<<endl);
        parameterStringClass.validateNumFields(2);
        if (parameterStringClass.getString(2).length() >  0) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "You have provided too many parameters for this command. temperature is no longer a staged parameter."<<endl);
        }
        temperature = myAtoF(userVariables,(parameterStringClass.getString(1)).c_str());
        return;
    }
    if  (((parameterStringClass.getString(0)).compare("-DC") == 0) || ((parameterStringClass.getString(0)).compare("dutyCycle") == 0))  {
        parameterStringClass.validateNumFields(2);
        dutyCycle = (myAtoF(userVariables,(parameterStringClass.getString(1)).c_str()));
        //dutyCycleArray.push_back(myAtoF(userVariables,(parameterStringClass.getString(1)).c_str()));
        //dutyCyclePriority.push_back(     myAtoI(userVariables,(parameterStringClass.getString(2)).c_str()));
        //dutyCycleArray[d]    = myAtoF(userVariables,(parameterStringClass.getString(1)).c_str());
        //dutyCyclePriority[d] = myAtoF(userVariables,(parameterStringClass.getString(2)).c_str());
        if (dutyCycle < 1.00) {
            setForceAndStericScrubber = true; //setForceAndStericScrubber is set automatically now, based on the value of dutyCycle;
            MMBLOG_FILE_FUNC_LINE(INFO, "You have set dutyCycle to "<<dutyCycle<<" .  Therefore, setForceAndStericScrubber has been set to "<<setForceAndStericScrubber <<endl);
        }
        if (dutyCycle == 1.00) {
            setForceAndStericScrubber = false; //setForceAndStericScrubber is set automatically now, based on the value of dutyCycle;
            MMBLOG_FILE_FUNC_LINE(INFO, "You have set dutyCycle to "<<dutyCycle<<" .  Therefore, setForceAndStericScrubber has been set to "<<setForceAndStericScrubber <<endl);
        }
        if ((dutyCycle > 1.00) ||
                (fabs(dutyCycle) < 1E-14) )
        {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "dutyCycle must lie in the open interval (0,1) ! You have tried to set it to : "<<dutyCycle <<endl);
        }
        if ((dutyCycle<1.0) &&   (densityContainer.numDensityStretches()>0)) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "you can't set dutyCycle<1.0 if you also issue the fitToDensity command"<<endl);
        }
        if ((dutyCycle<1.0) &&   (electroDensityContainer.numDensityStretches()>0)) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "you can't set dutyCycle<1.0 if you also issue the fitElectrostaticDensity command"<<endl);
        }
        //if (verbose) MMBLOG_FILE_FUNC_LINE(" dutyCycleArray["<<d <<"] ="<<dutyCycleArray[d]<< "from String :"<<parameterStringClass.getString(1)<<endl;  
        //if (verbose) MMBLOG_FILE_FUNC_LINE(" dutyCycleArray["<<d <<"] ="<<dutyCycleArray[d]<< "from String :"<<parameterStringClass.getString(1)<<endl;  
        //d++;
        return;
    }
    if ( ((parameterStringClass.getString(0)).compare("-GCSF") == 0) || ((parameterStringClass.getString(0)).compare("globalCoulombScaleFactor") == 0) )  {  
        parameterStringClass.validateNumFields(2);
        globalCoulombScaleFactor = (myAtoF(userVariables,(parameterStringClass.getString(1)).c_str()));
        SimTK_ERRCHK_ALWAYS((parameterStringClass.getString(2).size() == 0), "[ParameterReader.cpp] ","globalCoulombScaleFactor takes only a single parameter, please remove excess input from this line.  This is now a global rather than staged parameter.  This is a new syntax introduced in revision 201, sorry for the inconvenience.");
        return;
    }
    if  (((parameterStringClass.getString(0)).compare("-GVSF") == 0) ||((parameterStringClass.getString(0)).compare("globalVdwScaleFactor") == 0) )  {
        parameterStringClass.validateNumFields(2);
        globalVdwScaleFactor = ( myAtoF(userVariables,(parameterStringClass.getString(1)).c_str()));
        SimTK_ERRCHK_ALWAYS((parameterStringClass.getString(2).size() == 0), "[ParameterReader.cpp] ","globalVdwScaleFactor takes only a single parameter, please remove excess input from this line. This is now a global rather than staged parameter.  This is a new syntax introduced in revision 201, sorry for the inconvenience.");
        return;
    }
    if (((parameterStringClass.getString(0)).compare("-PRIOR") ==0) || ((parameterStringClass.getString(0)).compare("prioritize") ==0))  {
        parameterStringClass.validateNumFields(2);
        prioritize = aToBool(parameterStringClass.getString(0), (parameterStringClass.getString(1)).c_str());     
        return;
    }
    if (((parameterStringClass.getString(0)).compare("useNACappingHydroxyls") ==0))  {
        parameterStringClass.validateNumFields(2);
        if (myBiopolymerClassContainer.getNumBiopolymers() >0) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "You must specify the useNACappingHydroxyls parameter BEFORE you specify the first biopolymer (DNA, RNA, protein, etc.)"<<endl);
        }
        useNACappingHydroxyls = aToBool(parameterStringClass.getString(0), (parameterStringClass.getString(1)).c_str());     
        return;
    }
    if (((parameterStringClass.getString(0)).compare("proteinCapping") ==0))  {
        parameterStringClass.validateNumFields(2);
        if (myBiopolymerClassContainer.getNumBiopolymers() >0) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "You must specify the proteinCapping parameter BEFORE you specify the first biopolymer (DNA, RNA, protein, etc.)"<<endl);
        }
        proteinCapping = aToBool(parameterStringClass.getString(0), (parameterStringClass.getString(1)).c_str());     
        return;
    }
    if (((parameterStringClass.getString(0)).compare("excludedVolumeRadius") ==0))  {
        parameterStringClass.validateNumFields(2);
        excludedVolumeRadius = myAtoF(userVariables,(parameterStringClass.getString(1)).c_str());     
        return;
    }
    if (((parameterStringClass.getString(0)).compare("excludedVolumeStiffness") ==0))  {
        parameterStringClass.validateNumFields(2);
        excludedVolumeStiffness = myAtoF(userVariables,(parameterStringClass.getString(1)).c_str());      

        return;
    }
    if (((parameterStringClass.getString(0)).compare("cutoffAngle") ==0)  )  {
        parameterStringClass.validateNumFields(2);
        cutoffAngle = myAtoF(userVariables,(parameterStringClass.getString(1)).c_str());      
        MMBLOG_FILE_FUNC_LINE(DEBUG, "  cutoffAngle ="<<cutoffAngle<<endl);
        return;
    }
    if (((parameterStringClass.getString(0)).compare("cutoffRadius") ==0)  )  {
        parameterStringClass.validateNumFields(2);
        cutoffRadius = myAtoF(userVariables,(parameterStringClass.getString(1)).c_str()); 
        MMBLOG_FILE_FUNC_LINE(DEBUG, "  cutoffRadius ="<<cutoffRadius<<" nm "<<endl);
        return;
    }
    if (((parameterStringClass.getString(0)).compare("densityAtomFraction") ==0)  )  {
        parameterStringClass.validateNumFields(2);
        densityAtomFraction = myAtoF(userVariables,(parameterStringClass.getString(1)).c_str());
        if ((densityAtomFraction < 0.0) || (densityAtomFraction > 1.0)) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, " densityAtomFraction must be in the range of 0 to 1."<<endl);
        }
        return;
    }
    if (  ((parameterStringClass.getString(0)).compare("densityFileName") ==0)  )  {
        MMBLOG_FILE_FUNC_LINE(CRITICAL, " This syntax is no longer supported. Please instead use: density densityFileName <name of file>"<<endl);
    }
    if (parameterStringClass.getString(0).compare("density") == 0) {
        MMBLOG_FILE_FUNC_LINE(ALWAYS, "The density family of commands and parameters supports the following >"<<endl);
        MMBLOG_FILE_FUNC_LINE(ALWAYS, "Syntax:  "<<endl);
        MMBLOG_FILE_FUNC_LINE(ALWAYS, "    density densityFileName <name of file, in any of the supported formats, use right suffix>"<<endl);
        MMBLOG_FILE_FUNC_LINE(ALWAYS, "        sets the name of the input density file."<<endl);
        MMBLOG_FILE_FUNC_LINE(ALWAYS, "    density forceConstant <double>"<<endl);
        MMBLOG_FILE_FUNC_LINE(ALWAYS, "        sets the scaling factor for the density fitting forces. Defaults to 1.0. "<<endl);
        MMBLOG_FILE_FUNC_LINE(ALWAYS, "    density clear"<<endl);
        MMBLOG_FILE_FUNC_LINE(ALWAYS, "        clears the existing, loaded density map. This forces a reload at the next MMB stage. If you do not call this, the density map loaded in a previous stage will be reused for time savings."<<endl);
        if (parameterStringClass.getString(1).compare("densityFileName") ==0) {
            parameterStringClass.validateNumFields(3);
            myDensityMap.setDensityFileName(parameterStringClass.getString(2));
            return;
        }
        if (parameterStringClass.getString(1).compare("forceConstant") ==0) {
            parameterStringClass.validateNumFields(3);
            myDensityMap.setForceConstant(myAtoF(userVariables,parameterStringClass.getString(2).c_str()));
            return;
        }
        if (parameterStringClass.getString(1).compare("clear") ==0) {
            parameterStringClass.validateNumFields(2);
            myDensityMap.initializeMap();
            return;
	}
    }
    if (((parameterStringClass.getString(0)).compare("densityForceConstant") ==0)  )  {
        MMBLOG_FILE_FUNC_LINE(CRITICAL, " This syntax is no longer supported. Please instead use: density forceConstant                 "<<endl);
        //parameterStringClass.validateNumFields(2);
        //densityForceConstant = myAtoF(userVariables,(parameterStringClass.getString(1)).c_str());
        return;
    }
    if (((parameterStringClass.getString(0)).compare("densityNoiseTemperature") ==0)  )  {
        MMBLOG_FILE_FUNC_LINE(INFO, "densityNoiseTemperature provides the temperature for the Planck Law simulated density noise. Syntax: densityNoiseTemperature <temperature, double>."<<endl);
        parameterStringClass.validateNumFields(2);
        densityNoiseTemperature  = myAtoF(userVariables,(parameterStringClass.getString(1)).c_str());
        if ((densityNoiseTemperature < 0.0) ) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "densityNoiseTemperature should be > 0.0. You specified : "<<densityNoiseTemperature <<endl);
        }
        return;
    }
    if (((parameterStringClass.getString(0)).compare("densityNoiseScale") ==0)  )  {
        MMBLOG_FILE_FUNC_LINE(INFO, "densityNoiseScale provides an overall scale for the  Planck Law simulated density noise. Syntax: densityNoiseScale <double>."<<endl);
        parameterStringClass.validateNumFields(2);
        densityNoiseScale  = myAtoF(userVariables,(parameterStringClass.getString(1)).c_str());
        if ((densityNoiseScale < 0.0) ) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "densityNoiseScale should be > 0.0. You specified : "<<densityNoiseScale <<endl);
        }
        return;
    }
    if ((parameterStringClass.getString(0)).compare("densityReportAtEachAtomPosition" ) ==0) {
        MMBLOG_FILE_FUNC_LINE(INFO, "densityReportAtEachAtomPosition provides the local density at each atom position. This could get  somewhat verbose so is set to 0 by default. It is intended to give us an idea of how to override the atomicNumber property for better fitting performance. Syntax: densityReportAtEachAtomPosition <bool>."<<endl);
        parameterStringClass.validateNumFields(2);
        densityReportAtEachAtomPosition = aToBool(parameterStringClass.getString(0),(parameterStringClass.getString(1)).c_str());    
        return;
    }
    if ((parameterStringClass.getString(0)).compare("densityNoiseComputeAutocorrelation" ) ==0) {
        parameterStringClass.validateNumFields(2);
        if (parameterStringClass.getString(1).length() == 0) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "You have not provided enough parameters for this command."<<endl);
        }
        if (parameterStringClass.getString(2).length() >  0) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "You have provided too many parameters for this command."<<endl);
        }
        densityNoiseComputeAutocorrelation = aToBool(parameterStringClass.getString(0),(parameterStringClass.getString(1)).c_str());    
        return;
    }
    if ( ((parameterStringClass.getString(0)).compare("overrideAtomicProperty") == 0)  )    { 
        MMBLOG_FILE_FUNC_LINE(ALWAYS,
                "Syntax: overrideAtomicProperty <atom name>  <property> <value>  "<<endl
                <<"Where <atom name> is the name of the atom whose property you want to override. All atoms in the system will have their <property> set to <value>." <<endl
                <<"currently the only <property> supported is atomicNumber. "<<endl
                <<"<value> is a double. In the case that <property> is an int, the <value> will be cast to int."<<endl
                <<"The expected use of this is to override atomicNumber for density map fitting. If phosphate are not visible, you may wish to set atoms such as P, OP1, etc, to zero or even a negative value. This means they will feel zero or negative fitting force."<<endl);
        parameterStringClass.validateNumFields(4);
        AtomicPropertyOverrideStruct atomicPropertyOverrideStruct ;
        atomicPropertyOverrideStruct.atomName = parameterStringClass.getString(1);
        atomicPropertyOverrideStruct.property = parameterStringClass.getString(2);
        atomicPropertyOverrideStruct.value    = myAtoF(userVariables,parameterStringClass.getString(3).c_str());
         myBiopolymerClassContainer.atomicPropertyOverrideVector.push_back(atomicPropertyOverrideStruct);       
        return;
    } // overrideAtomicProperty  

    if ((parameterStringClass.getString(0)).compare("densityFitPhosphates" ) ==0) {
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "This parameter is no longer supported.  You should now just use the overrideAtomicProperty command. Just issue overrideAtomicProperty to get the syntax."<<endl);
        parameterStringClass.validateNumFields(2);
        if (parameterStringClass.getString(1).length() == 0) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "You have not provided enough parameters for this command."<<endl);
        }
        if (parameterStringClass.getString(2).length() >  0) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "You have provided too many parameters for this command."<<endl);
        }
        densityFitPhosphates = aToBool(parameterStringClass.getString(0),(parameterStringClass.getString(1)).c_str());    
        if (!(densityFitPhosphates)){
            AtomicPropertyOverrideStruct atomicPropertyOverrideStruct ;
            atomicPropertyOverrideStruct.property = "atomicNumber";
            atomicPropertyOverrideStruct.value = 0.0;   
            atomicPropertyOverrideStruct.atomName = "P";   myBiopolymerClassContainer.atomicPropertyOverrideVector.push_back(atomicPropertyOverrideStruct);
            atomicPropertyOverrideStruct.atomName = "OP1"; myBiopolymerClassContainer.atomicPropertyOverrideVector.push_back(atomicPropertyOverrideStruct);
            atomicPropertyOverrideStruct.atomName = "OP2"; myBiopolymerClassContainer.atomicPropertyOverrideVector.push_back(atomicPropertyOverrideStruct);
            atomicPropertyOverrideStruct.atomName = "O5'"; myBiopolymerClassContainer.atomicPropertyOverrideVector.push_back(atomicPropertyOverrideStruct);
            //atomicPropertyOverrideStruct.atomName = "O5'"; myBiopolymerClassContainer.atomicPropertyOverrideVector.push_back(atomicPropertyOverrideStruct);
            atomicPropertyOverrideStruct.atomName = "O3'"; myBiopolymerClassContainer.atomicPropertyOverrideVector.push_back(atomicPropertyOverrideStruct);
            //atomicPropertyOverrideStruct.atomName = "O3'"; myBiopolymerClassContainer.atomicPropertyOverrideVector.push_back(atomicPropertyOverrideStruct);
        } else {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "You have tried to set densityFitPhosphates to TRUE. please don't do this! Just leave out the command altogether. Actually this whole command is obsolete."<<endl);
        }
        return;
    }
    if (((parameterStringClass.getString(0)).compare("densityMapActivate") ==0)  )  {
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "This parameter is no longer supported.  You should now just use the fitToDensity command.. see the reference guide."<<endl);
        return;
    }
    if (((parameterStringClass.getString(0)).compare("electroDensityFileName") ==0)  )  {
        parameterStringClass.validateNumFields(2);
        electroDensityFileName = parameterStringClass.getString(1);    
        return;
    }
    if (((parameterStringClass.getString(0)).compare("electroDensityForceConstant") ==0)  )  {
        parameterStringClass.validateNumFields(2);
        electroDensityForceConstant = myAtoF(userVariables,(parameterStringClass.getString(1)).c_str());
        return;
    }







    if (((parameterStringClass.getString(0)).compare("sphericalHelix") ==0)  )  {
	//MonoAtoms myMonoAtoms;
	//myMonoAtoms.clear();
	spiral.writeSyntax();    
	spiral.parseInput(userVariables,parameterStringClass.getString(1),parameterStringClass.getString(2), parameterStringClass.getString(3), parameterStringClass.getString(4),myMonoAtomsContainer);
	//spiral.monoAtoms;
        return;
            //myMonoAtomsContainer.addMonoAtoms(myMonoAtoms);
    }

    if (((parameterStringClass.getString(0)).compare("vanderWallSphereRadius") ==0)  )  {
        parameterStringClass.validateNumFields(2);
        vanderWallSphereRadius = myAtoF(userVariables,(parameterStringClass.getString(1)).c_str());   
        MMBLOG_FILE_FUNC_LINE(DEBUG, "vanderWallSphereRadius ="<<vanderWallSphereRadius<<endl);
        return;
    }
    if (((parameterStringClass.getString(0)).compare("-REMARK") ==0) || (((parameterStringClass.getString(0)).substr(0,1)).compare("#"      ) ==0) || ((parameterStringClass.getString(0)).compare("") ==0) )  { 
        // do nothing; this is a comment or blank line
        return;
    }
    if (((parameterStringClass.getString(0)).substr(0,1)).compare("@") ==0)  { 
        parameterStringClass.validateNumFields(2);
        // this is a user specified variable.  For now it works only for doubles and ints.          
        userVariables[parameterStringClass.getString(0)] =  myAtoF(userVariables,(parameterStringClass.getString(1)).c_str());//  (double)atof((parameterStringClass.getString(1)).c_str());  // this maybe should be the only atof in the whole program.   
        if (((parameterStringClass.getString(0)).compare("@CURRENTSTAGE") == 0)) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "@CURRENTSTAGE is a reserved variable, you cannot change it!"<<endl);
        }
        MMBLOG_FILE_FUNC_LINE(INFO, "read variable "<<parameterStringClass.getString(0) <<", successfully set to "<<userVariables[parameterStringClass.getString(0)]<<endl);
        //#userVariables[parameterStringClass.getString(0)] = (double)atof((parameterStringClass.getString(1)).c_str());  // this maybe should be the only atof in the whole program.   

        return;
    }
    if ((parameterStringClass.getString(0)).compare("setChiBondMobility") ==0) {
        parameterStringClass.validateNumFields(2);
        setChiBondMobility  = myAtoI(userVariables,(parameterStringClass.getString(1)).c_str());        
        return;
    }
    if ((parameterStringClass.getString(0)).compare("fitDefaultTolerance"         ) ==0) {

        MMBLOG_FILE_FUNC_LINE(CRITICAL, "read variable "<<parameterStringClass.getString(0) <<".  This is an obsolete parameter!"<<endl);
        parameterStringClass.validateNumFields(2);
        fitDefaultTolerance = myAtoF(userVariables,(parameterStringClass.getString(1)).c_str());        
        return;
    }
    if ((parameterStringClass.getString(0)).compare("globalBondStretchScaleFactor") ==0) {
        parameterStringClass.validateNumFields(2);
        globalBondStretchScaleFactor = myAtoF(userVariables,(parameterStringClass.getString(1)).c_str());        
        return;
    }
    if ((parameterStringClass.getString(0)).compare("globalBondTorsionScaleFactor") ==0) {
        parameterStringClass.validateNumFields(2);
        globalBondTorsionScaleFactor = myAtoF(userVariables,(parameterStringClass.getString(1)).c_str());        
        return;
    }
    if ((parameterStringClass.getString(0)).compare("globalAmberImproperTorsionScaleFactor") ==0) {
        parameterStringClass.validateNumFields(2);
        globalAmberImproperTorsionScaleFactor = myAtoF(userVariables,(parameterStringClass.getString(1)).c_str());        
        return;
    }
    if ((parameterStringClass.getString(0)).compare("globalBondBendScaleFactor") ==0) {
        parameterStringClass.validateNumFields(2);
        globalBondBendScaleFactor = myAtoF(userVariables,(parameterStringClass.getString(1)).c_str());        
        return;
    }
    if ((parameterStringClass.getString(0)).compare("globalGbsaScaleFactor") ==0) {
        parameterStringClass.validateNumFields(2);
        globalGbsaScaleFactor = myAtoF(userVariables,(parameterStringClass.getString(1)).c_str());        
        if ((!(fabs(globalGbsaScaleFactor) <= 0.0000001)) && safeParameters ) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "Gbsa isn't supported in this release.  Please set globalGbsaScaleFactor 0.  You can set safeParameters FALSE to override this message. "<<endl);
        }   
        return;
    }
    if ((parameterStringClass.getString(0)).compare("vmdOutput")    ==0) {
        parameterStringClass.validateNumFields(2);
        vmdOutput    = myAtoI(userVariables,(parameterStringClass.getString(1)).c_str());    
        return;
    }
    if ((parameterStringClass.getString(0)).compare("waterDroplet") ==0) {
        MMBLOG_FILE_FUNC_LINE(ALWAYS, "syntax:  waterDroplet <droplet chainID>  X Y Z radius [tetherStrength]"<<endl);
        //MMBLOG_FILE_FUNC_LINE("Alternatively, one can specify an atom as the center. The tethers will not be applied to the atom, but rather to the Ground at the atom's initial position."<<endl;
        //MMBLOG_FILE_FUNC_LINE("syntax:  waterDroplet <droplet chainID> AtomLocationInGround <atom chain ID> <atom residue number> <atom name> radius [tetherStrength [tetherType]]"<<endl;
        MMBLOG_FILE_FUNC_LINE(ALWAYS, "Note that in MMB 2.10 and earlier, we took the radius in Å.  We are going back to nm for consistency, with apologies for the confusion."<<endl);

        WaterDroplet myWaterDroplet;
        myWaterDroplet.chainID   = parameterStringClass.getString(1);  
        //if (parameterStringClass.getString(2).compare ("Coordinates") == 0) {        
	myWaterDroplet.center[0] = myAtoF(userVariables,(parameterStringClass.getString(2)).c_str());
	myWaterDroplet.center[1] = myAtoF(userVariables,(parameterStringClass.getString(3)).c_str());
	myWaterDroplet.center[2] = myAtoF(userVariables,(parameterStringClass.getString(4)).c_str());
        //}
        //else if (parameterStringClass.getString(1).compare ("AtomLocationInGround") == 0) {
        //}
        myWaterDroplet.setRadius ( myAtoF(userVariables,(parameterStringClass.getString(5)).c_str()));
        MMBLOG_FILE_FUNC_LINE(INFO, "Setting water droplet radius to : "<<myWaterDroplet.getRadius()<<endl);
        if (myWaterDroplet.getRadius() <=0) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "Droplet radius = "<<myWaterDroplet.getRadius()<<" .. less than zero!"<<endl);
        }
        if (parameterStringClass.getString(6).length() > 0) myWaterDroplet.tetherStrength = myAtoF(userVariables,(parameterStringClass.getString(6)).c_str());
        MMBLOG_FILE_FUNC_LINE(INFO, "Setting tether strength to "<<myWaterDroplet.tetherStrength<<endl);

        if (!(parameterStringClass.getString(5).length() > 0)) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "Invalid syntax!  You don't have enough parameters."<<endl);
        }

        if (myBiopolymerClassContainer.hasChainID(myWaterDroplet.chainID)) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "Can't use chain "<<myWaterDroplet.chainID<<" to identify this water droplet, because it is already being used by a biopolymer in your system."<<endl);
        }
        for (int i = 0 ; i < (int)waterDropletAboutResidueVector.size(); i++) {
            if (waterDropletAboutResidueVector[i].waterDropletChainID.compare(myWaterDroplet.chainID) == 0) {    
                MMBLOG_FILE_FUNC_LINE(CRITICAL, "Can't use chain "<<myWaterDroplet.chainID<<" to identify this water droplet, because it is already being used by another water droplet."<<endl);
            }
        } 
        if (waterDropletContainer.hasChainID(myWaterDroplet.chainID)) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "Can't use chain "<<myWaterDroplet.chainID<<" to identify this water droplet, because it is already being used by another water droplet."<<endl);
        }
        waterDropletContainer.add(myWaterDroplet);
        return;
    }
    if ((parameterStringClass.getString(0)).compare("waterDropletAboutResidue") ==0) {
        cout<<"syntax:  waterDropletAboutResidue <biopolymer Chain ID> <biopolymer Residue Number> <radius> <tether Strength> <water Droplet Chain ID>"<<endl;
        WaterDropletAboutResidueStruct myWaterDropletAboutResidueStruct;
        myWaterDropletAboutResidueStruct.biopolymerChainID   = parameterStringClass.getString(1);
        myWaterDropletAboutResidueStruct.residue = myBiopolymerClassContainer.residueID(userVariables,(parameterStringClass.getString(2)),myWaterDropletAboutResidueStruct.biopolymerChainID  );
        myWaterDropletAboutResidueStruct.radius = myAtoF(userVariables,(parameterStringClass.getString(3)).c_str());
        myWaterDropletAboutResidueStruct.tetherStrength = myAtoF(userVariables,(parameterStringClass.getString(4)).c_str());
        if (! (parameterStringClass.getString(5)).length()) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "Please specify a valid chain ID to designate this water droplet."<<endl);
        }
        myWaterDropletAboutResidueStruct.waterDropletChainID  = parameterStringClass.getString(5);

        if ( (parameterStringClass.getString(6)).length()) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, " : You have specified too many parameters for this water droplet."<<endl);
        }
        if (myBiopolymerClassContainer.hasChainID(myWaterDropletAboutResidueStruct.waterDropletChainID)) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, " : Can't use chain "<<myWaterDropletAboutResidueStruct.waterDropletChainID<<" to identify this water droplet, because it is already being used by a biopolymer in your system."<<endl);
        }
        for (int i = 0 ; i < (int)waterDropletAboutResidueVector.size(); i++) {
            if (waterDropletAboutResidueVector[i].waterDropletChainID.compare(myWaterDropletAboutResidueStruct.waterDropletChainID) == 0) {    
                MMBLOG_FILE_FUNC_LINE(CRITICAL, " : Can't use chain "<<myWaterDropletAboutResidueStruct.waterDropletChainID<<" to identify this water droplet, because it is already being used by another water droplet."<<endl);
            } 
        } 
        waterDropletAboutResidueVector.push_back(myWaterDropletAboutResidueStruct);
        return;
    }
    if ((parameterStringClass.getString(0)).compare("waterDropletMake") ==0) {
        parameterStringClass.validateNumFields(2);
        waterDropletMake = aToBool(parameterStringClass.getString(0),(parameterStringClass.getString(1)).c_str());    
        return;
    }
    if ((parameterStringClass.getString(0)).compare("waterInertiaMultiplier") ==0) {
        parameterStringClass.validateNumFields(2);
        waterInertiaMultiplier = myAtoF   (userVariables,(parameterStringClass.getString(1)).c_str());    
        return;
    }
    if ((parameterStringClass.getString(0)).compare("weldToGround") ==0) {
        parameterStringClass.validateNumFields(2);
        weldToGround = aToBool(parameterStringClass.getString(0),(parameterStringClass.getString(1)).c_str());    
        return;
    }
    if ((parameterStringClass.getString(0)).compare("wkdpGlobalBondTorsionScaleFactor") ==0) {
        parameterStringClass.validateNumFields(2);
        MMBLOG_FILE_FUNC_LINE(CRITICAL, " The wkdpGlobalBondTorsionScaleFactor is obsolete."<<endl);
        wkdpGlobalBondTorsionScaleFactor = myAtoF(userVariables,(parameterStringClass.getString(1)).c_str());    
        return;
    }
    if ((parameterStringClass.getString(0)).compare("writeCoordinates") ==0) {
        parameterStringClass.validateNumFields(2);
        writeCoordinates = aToBool(parameterStringClass.getString(0),(parameterStringClass.getString(1)).c_str());    
        return;
    }
    if ((parameterStringClass.getString(0)).compare("writeDoublePrecisionTrajectories") ==0) {
        parameterStringClass.validateNumFields(2);
        writeDoublePrecisionTrajectories = aToBool(parameterStringClass.getString(0),(parameterStringClass.getString(1)).c_str());    
        return;
    }
    if ((parameterStringClass.getString(0)).compare("writeFrameFile"  ) ==0) {
        parameterStringClass.validateNumFields(2);
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "writeFrameFile is obsolete."<<endl);
        writeFrameFile   = aToBool(parameterStringClass.getString(0),(parameterStringClass.getString(1)).c_str());    

        return;
    }
    if ((parameterStringClass.getString(0)).compare("writeLastFrameFile"  ) ==0) {
        parameterStringClass.validateNumFields(2);
        writeLastFrameFile   = aToBool(parameterStringClass.getString(0),(parameterStringClass.getString(1)).c_str());    
        return;
    }
    if ((parameterStringClass.getString(0)).compare("workingDirectory"    ) ==0) {
        parameterStringClass.validateNumFields(2);
        workingDirectory     = get_and_set_working_path (parameterStringClass.getString(1)); // So now we have set workingDirectory, and cd'd into it.   
        return;
    }
    if ((parameterStringClass.getString(0)).compare("removeDensityForcesFromRigidStretches"  ) ==0) {
        MMBLOG_FILE_FUNC_LINE(INFO, " when set to True, this tells MMB to remove all residues in mobilizer stretches with BondMobility::Rigid, from the densityContainer. That is to say, any residues that have been rigidified will not have density gradient forces applied  .  "<<endl);
        MMBLOG_FILE_FUNC_LINE(INFO, " Syntax : "<<endl); 
        MMBLOG_FILE_FUNC_LINE(INFO, " removeDensityForcesFromRigidStretches <bool>  "<<endl); 
        parameterStringClass.validateNumFields(2);
        removeDensityForcesFromRigidStretches   = aToBool(parameterStringClass.getString(0),(parameterStringClass.getString(1)).c_str());    
        return;
    }
    if ((parameterStringClass.getString(0)).compare("removeRigidBodyMomentum"  ) ==0) {
        parameterStringClass.validateNumFields(2);
        removeRigidBodyMomentum   = aToBool(parameterStringClass.getString(0),(parameterStringClass.getString(1)).c_str());    
        return;
    }
    if ((parameterStringClass.getString(0)).compare("removeMomentumPeriod"  ) ==0) {
        parameterStringClass.validateNumFields(2);
        removeMomentumPeriod   = myAtoF(userVariables,(parameterStringClass.getString(1)).c_str());
        SimTK_ERRCHK_ALWAYS(
                (removeMomentumPeriod > 0),
                "[ParameterReader.cpp]","Parameter removeMomentumPeriod must be greater than zero.");
        return;
    }
    if ((parameterStringClass.getString(0)).compare("rigidifyFormedHelices") ==0) {
        parameterStringClass.validateNumFields(2);
        rigidifyFormedHelices = aToBool(parameterStringClass.getString(0),(parameterStringClass.getString(1)).c_str());    
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "This is no longer a user configurable parameter.  .  "<<endl);
        return;
    }
    if ((parameterStringClass.getString(0)).compare("velocityRescalingInterval"              ) ==0) {
        parameterStringClass.validateNumFields(2);
        velocityRescalingInterval = myAtoF(userVariables,(parameterStringClass.getString(1)).c_str());    
        return;
    }
    if ((parameterStringClass.getString(0)).compare("verbose"              ) ==0) {
        parameterStringClass.validateNumFields(2);
        //verbose = myAtoI(userVariables,(parameterStringClass.getString(1)).c_str());    
        verbose = aToBool(parameterStringClass.getString(0),(parameterStringClass.getString(1)).c_str());    
        return;
    }

    if 
        (
         ((parameterStringClass.getString(0)).compare("addSelectedAtoms")          ==0) || 
         ((parameterStringClass.getString(0)).compare("addProteinBackboneSterics") ==0) ||
         ((parameterStringClass.getString(0)).compare("addRNABackboneSterics")     ==0)  
        ) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "This command is no longer supported!  Instead, use the contact command, for each chain and stretch of residues you want to apply it to."<<endl);
        }
    if 
        (
         ((parameterStringClass.getString(0)).compare("addAllHeavyAtomSterics")    ==0) || 
         ((parameterStringClass.getString(0)).compare("addAllAtomSterics")         ==0) //|| 
         //((parameterStringClass.getString(0)).compare("addSelectedAtoms")          ==0) || 
         //((parameterStringClass.getString(0)).compare("addProteinBackboneSterics") ==0) ||
         //((parameterStringClass.getString(0)).compare("addRNABackboneSterics")     ==0)  
        ) {

            parameterStringClass.validateNumFields(2);
            if (addAllHeavyAtomSterics && ((parameterStringClass.getString(0)).compare("addAllHeavyAtomSterics") ==0))      
            {MMBLOG_FILE_FUNC_LINE(CRITICAL, "You may only set the addAllHeavyAtomSterics parameter to TRUE once, and afterwards cannot set it to FALSE."<<endl); }
            if (addAllAtomSterics &&  ((parameterStringClass.getString(0)).compare("addAllAtomSterics") ==0))                
            {MMBLOG_FILE_FUNC_LINE(CRITICAL, "You may only set the addAllAtomSterics         parameter to TRUE once, and afterwards cannot set it to FALSE."<<endl); }
            if (addSelectedAtoms  &&  ((parameterStringClass.getString(0)).compare("addSelectedAtoms") ==0))             
            {MMBLOG_FILE_FUNC_LINE(CRITICAL, "You may only set the addSelectedAtoms          parameter to TRUE once, and afterwards cannot set it to FALSE."<<endl); }
            if (addProteinBackboneSterics &&   ((parameterStringClass.getString(0)).compare("addProteinBackboneSterics") ==0)) 
            {MMBLOG_FILE_FUNC_LINE(CRITICAL, "You may only set the addProteinBackboneSterics parameter to TRUE once, and afterwards cannot set it to FALSE."<<endl); }
            if (addRNABackboneSterics &&   ((parameterStringClass.getString(0)).compare("addRNABackboneSterics") ==0)) 
            {MMBLOG_FILE_FUNC_LINE(CRITICAL, "You may only set the addRNABackboneSterics parameter to TRUE once, and afterwards cannot set it to FALSE."<<endl); }

            String myEdge = "";
            if  ((parameterStringClass.getString(0)).compare("addAllHeavyAtomSterics") ==0)    {
                addAllHeavyAtomSterics =aToBool(parameterStringClass.getString(0),(parameterStringClass.getString(1)).c_str());    
                myEdge = "AllHeavyAtomSterics";  }  
                if  ((parameterStringClass.getString(0)).compare("addAllAtomSterics") ==0)         {
                    addAllAtomSterics =aToBool(parameterStringClass.getString(0),(parameterStringClass.getString(1)).c_str());         
                    myEdge = "AllAtomSterics";  }  
                    // these three are confusing, so I'm retiring them:
                    if  ((parameterStringClass.getString(0)).compare("addSelectedAtoms") ==0)          {
                        addSelectedAtoms =aToBool(parameterStringClass.getString(0),(parameterStringClass.getString(1)).c_str());          
                        myEdge = "SelectedAtoms";  }  
                        if  ((parameterStringClass.getString(0)).compare("addProteinBackboneSterics") ==0) {
                            addProteinBackboneSterics =aToBool(parameterStringClass.getString(0),(parameterStringClass.getString(1)).c_str()); 
                            myEdge = "ProteinBackboneSterics";  }  
                            if  ((parameterStringClass.getString(0)).compare("addRNABackboneSterics") ==0)     {
                                addRNABackboneSterics =aToBool(parameterStringClass.getString(0),(parameterStringClass.getString(1)).c_str());     
                                myEdge = "RNABackboneSterics";  }  


                                if (
                                        (addAllHeavyAtomSterics    &&  ((parameterStringClass.getString(0)).compare("addAllHeavyAtomSterics") ==0)   ) ||
                                        (addAllAtomSterics         &&  ((parameterStringClass.getString(0)).compare("addAllAtomSterics") ==0)   )      ||
                                        (addProteinBackboneSterics &&  ((parameterStringClass.getString(0)).compare("addProteinBackboneSterics") ==0)) ||
                                        (addSelectedAtoms          &&  ((parameterStringClass.getString(0)).compare("addSelectedAtoms") ==0)         ) ||
                                        (addProteinBackboneSterics &&  ((parameterStringClass.getString(0)).compare("addProteinBackboneSterics") ==0)) ||
                                        (addRNABackboneSterics     &&  ((parameterStringClass.getString(0)).compare("addRNABackboneSterics") ==0))
                                   )
                                {
                                    parameterStringClass.validateNumFields(2);
                                    int t = 0;
                                    for (t = 0; t<myBiopolymerClassContainer.getNumBiopolymers(); t++) {

                                        BiopolymerClass tempBiopolymerClass = (myBiopolymerClassContainer.updBiopolymerClass(t));
                                        ContactStretch        myContactStretch;
                                        //ContactStretch.ContactStretchIsTwoTransformForce = "contact";   
                                        myContactStretch.ContactScheme= myEdge;                           
                                        myContactStretch.setChain ( tempBiopolymerClass.getChainID());
                                        myContactStretch.setStartResidue ( tempBiopolymerClass.getFirstResidueID());
                                        //myContactStretch.SecondBPChain = tempBiopolymerClass.getChainID();
                                        myContactStretch.setEndResidue  (tempBiopolymerClass.getLastResidueID());
                                        myBiopolymerClassContainer.updBiopolymerClass(myContactStretch.getChain()).validateResidueID(myContactStretch.getStartResidue());
                                        myBiopolymerClassContainer.updBiopolymerClass(myContactStretch.getChain()).validateResidueID(myContactStretch.getEndResidue());
                                        //myContactStretch.SecondBPEdge = "";
                                        //myContactStretch.ContactStretchPriority = 1;
                                        contactContainer.addContactToVector(myContactStretch,myBiopolymerClassContainer);

                                    }
                                }

                                return;
        }
    if ((parameterStringClass.getString(0)).compare("addNASTForces") ==0) {
        parameterStringClass.validateNumFields(2);
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "The addNASTForces command has been discontinued."<<endl);
        //addNASTForces = aToBool(parameterStringClass.getString(0),  (parameterStringClass.getString(1)).c_str());    
        return;
    }
    if ((parameterStringClass.getString(0)).compare("addBackboneOxygenForces") ==0) {
        parameterStringClass.validateNumFields(2);
        addBackboneOxygenForces = aToBool(parameterStringClass.getString(0),(parameterStringClass.getString(1)).c_str());
        if (safeParameters) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "addBackboneOxygenForces is no longer a supported parameter."<<endl);
        }
        return;
    }
    if ((parameterStringClass.getString(0)).compare("addExcludedVolume"      ) ==0) {
        parameterStringClass.validateNumFields(2);
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "The addExcludedVolume command is no longer supported."<<endl);
        return;
    }
    if ((parameterStringClass.getString(0)).compare("loadTinkerParameterFile") ==0) {
	MMBLOG_FILE_FUNC_LINE(CRITICAL, "The loadTinkerParameterFile parameter is false by default, but is set to TRUE upon setting tinkerParameterFileName . You are no longer permitted to set loadTinkerParameterFile directly!. "<<endl);
        
        return;
    }

    if ( ( parameterStringClass.getString(0)).compare("useCIFFiles") == 0 )
    {
        //============================================ Set the CIF files usage
        parameterStringClass.validateNumFields        ( 2 );
        this->useCIFFileFormat                        = aToBool ( parameterStringClass.getString ( 0 ), ( parameterStringClass.getString ( 1 ) ).c_str() );
        return;
    }

    if (((parameterStringClass.getString(0)).compare("previousFrameFileName") ==0) ) {
        MMBLOG_FILE_FUNC_LINE(ALWAYS, "This command overrides the previous structure file name -- which defaults to last.[n-1].pdb"<<endl);
        MMBLOG_FILE_FUNC_LINE(ALWAYS, "Syntax is:  previousFrameFileName [input structure file name] "<<endl);
        parameterStringClass.validateNumFields(2);
        previousFrameFileName = parameterStringClass.getString(1);
        return;       
    }
    if (((parameterStringClass.getString(0)).compare("lastFrameFileName") ==0) ) {
        MMBLOG_FILE_FUNC_LINE(ALWAYS, "Syntax is:  lastFrameFileName [input structure file name] "<<endl);
        if (safeParameters) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "This command is not currently allowed!"<<endl);
        }
        parameterStringClass.validateNumFields(2);
        lastFrameFileName = parameterStringClass.getString(1);
        return;       
    }
    if (((parameterStringClass.getString(0)).compare("loadSequencesFromPdb") ==0) ||
        ((parameterStringClass.getString(0)).compare("loadSequencesFromPdbAndRenumber") ==0)) {
        MMBLOG_FILE_FUNC_LINE(ALWAYS, "Syntax is:  loadSequencesFromPdb [input structure file name] "<<endl);
        MMBLOG_FILE_FUNC_LINE(ALWAYS, "If the optional structure file is omitted, sequences will be read from last.(n-1).pdb (n being the current stage). "<<endl);
        
        if (currentStage == 1) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "You can't call loadSequencesFromPdb at stage 1! This stage does not use any input structure file."<<endl);
        }
        if(parameterStringClass.getString(3).length() > 0) 
        {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "loadSequencesFromPdb takes only two  parameters at most,  a pdb file name and (optionally) a chains prefix."<<endl);
        }
        if ((currentStage > firstStage) && (myBiopolymerClassContainer.getNumBiopolymers() >0) )
        {
            //MMBLOG_FILE_FUNC_LINE(" :We are not on first stage and some polymers are already been instantied, "<<previousFrameFileName<<" from previous stage has already been read. No need to do it again." << endl; 
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "We are not on first stage and some polymers are already been instantied, "<<previousFrameFileName<<" from previous stage has already been read. No need to do it again." << endl);// at stage 1! This stage does not use any input structure file."<<endl; 
            //return;
        }
        bool tempRenumberPdbResidues = 0; 
        if (parameterStringClass.getString(0) == "loadSequencesFromPdbAndRenumber") { tempRenumberPdbResidues = 1;}
        else if (parameterStringClass.getString(0) == "loadSequencesFromPdb")       { tempRenumberPdbResidues = 0;}
        else {MMBLOG_FILE_FUNC_LINE(CRITICAL, "Unexplained error!"<<endl); }

        // If not on first stage and some polymers are already been instantied, last.?.pd from previous stage has already been read. No need to do it again.

        // If no pdb file name provided or not on first stage, load last.?.pdb from previous stage.
        //if(currentStage > firstStage || parameterStringClass.getString(1).length() == 0)
        String chainsPrefix = parameterStringClass.getString(2);
        String pdbFileName = parameterStringClass.getString(1);
        if( pdbFileName.length() == 0)
        {
            MMBLOG_FILE_FUNC_LINE(INFO, "We load " << previousFrameFileName << " from previous stage." << endl);
            pdbFileName = previousFrameFileName; 
            //MMBLOG_FILE_FUNC_LINE(" : We are not on first stage so we load " << previousFrameFileName << " from previous stage." << endl;
            //loadSequencesFromPdb(previousFrameFileName,chainsPrefix, tempRenumberPdbResidues);
            //return;
        }
        else if (parameterStringClass.getString(1).length() > 0) 
        { 
            if ( chainsPrefix.length() > 1 ) {
                MMBLOG_FILE_FUNC_LINE(CRITICAL, "loadSequencesFromPdb: chains prefix must be empty or exactly one character."<<endl);
            }
            MMBLOG_FILE_FUNC_LINE(INFO, endl);
            readPreviousFrameFile = false;
            //return;
        }
        loadSequencesFromPdb(pdbFileName, chainsPrefix, tempRenumberPdbResidues);


        return;
    }

    if ((parameterStringClass.getString(0)).compare("calcBaseBodyFramesAtEveryTimeStep" ) ==0) {
        parameterStringClass.validateNumFields(2);
        if (parameterStringClass.getString(1).length() == 0) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "You have not provided enough parameters for this command."<<endl);
        }
        if (parameterStringClass.getString(2).length() >  0) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "You have provided too many parameters for this command."<<endl);
        }
        calcBaseBodyFramesAtEveryTimeStep  = myAtoI(userVariables,(parameterStringClass.getString(1)).c_str());    
        return;
    }
    if ((parameterStringClass.getString(0)).compare("calcEnergy" ) ==0) {
        parameterStringClass.validateNumFields(2);
        if (parameterStringClass.getString(1).length() == 0) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "You have not provided enough parameters for this command."<<endl);
        }
        if (parameterStringClass.getString(2).length() >  0) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "You have provided too many parameters for this command."<<endl);
        }
        calcEnergy = aToBool(parameterStringClass.getString(0),(parameterStringClass.getString(1)).c_str());    
        return;
    }
    if ((parameterStringClass.getString(0)).compare("checkSatisfied" ) ==0) {
        parameterStringClass.validateNumFields(2);
        if (parameterStringClass.getString(1).length() == 0) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "You have not provided enough parameters for this command."<<endl);
        }
        if (parameterStringClass.getString(2).length() >  0) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "You have provided too many parameters for this command."<<endl);
        }
        checkSatisfied  = aToBool(parameterStringClass.getString(0),(parameterStringClass.getString(1)).c_str());    
        return;
    } 

    if ((parameterStringClass.getString(0)).compare("constrainRigidSegments" ) ==0) {
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "This parameter is discontinued.  Use constrainChainRigidSegments instead."<<endl);

        /*MMBLOG_FILE_FUNC_LINE(" : Syntax: constrainRigidSegments <True | False>"<<endl; 

        parameterStringClass.validateNumFields(2);

        constrainRigidSegments  = aToBool(parameterStringClass.getString(0),(parameterStringClass.getString(1)).c_str());    */
        return;
    }
    if ((parameterStringClass.getString(0)).compare("constraintTolerance") ==0) {
        parameterStringClass.validateNumFields(2);
        constraintTolerance                = myAtoF(userVariables,(parameterStringClass.getString(1)).c_str());    
        return;
    }
    if ( ((parameterStringClass.getString(0)).compare( "includeAllNonBondAtomsInResidues") == 0) ||
        ((parameterStringClass.getString(0)).compare( "includeResidues") == 0) )  
    {
        cout<<"syntax:  includeResidues chain1 residue1    residue2 "<<endl;
        if (parameterStringClass.getString(3).length() == 0) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "You have not provided enough parameters for this command."<<endl);
        }
        if (parameterStringClass.getString(4).length() >  0) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "You have provided too many parameters for this command."<<endl);
        }
        for (ResidueID i = myBiopolymerClassContainer.residueID(userVariables,(parameterStringClass.getString(2)),parameterStringClass.getString(1) )  ; i <=myBiopolymerClassContainer.residueID (userVariables, (parameterStringClass.getString(3)),parameterStringClass.getString(1)); myBiopolymerClassContainer.updBiopolymerClass(parameterStringClass.getString(1)).incrementResidueID( i)  ) {
            IncludeAllNonBondAtomsInResidue myIncludeAllNonBondAtomsInResidue;
            myIncludeAllNonBondAtomsInResidue.setChain( parameterStringClass.getString(1));
            myIncludeAllNonBondAtomsInResidue.setResidue(i);  
            myBiopolymerClassContainer.updBiopolymerClass(myIncludeAllNonBondAtomsInResidue.getChain() ).validateResidueID( myIncludeAllNonBondAtomsInResidue.getResidue() );
            //includeAllNonBondAtomsInResidueVector.push_back(myIncludeAllNonBondAtomsInResidue);     
            physicsContainer.addStretch( myIncludeAllNonBondAtomsInResidue );
            if (i == myBiopolymerClassContainer.residueID (userVariables, (parameterStringClass.getString(3)),parameterStringClass.getString(1)))
                {break;}
        }
        return;
    }

    if ((parameterStringClass.getString(0)).compare( "includeIntraChainInterfaceResidues") == 0) {       
        MMBLOG_FILE_FUNC_LINE(ALWAYS,
                "syntax:  includeIntraChainInterfaceResidues <chain> <depth>    "<<endl
                <<"This will add all inter-BODY interfaces, to a <depth>, in the specified <chain>, to the physics zone.    "<<endl);
        //MMBLOG_FILE_FUNC_LINE(" You can omit <depth>, and it will be assumed that you want to set the depthain>, to the physics zone.    "<<endl;
        parameterStringClass.validateNumFields(3);       
        IncludeIntraChainInterface myIncludeIntraChainInterface;         
        myIncludeIntraChainInterface.Chain = parameterStringClass.getString(1);          
        myIncludeIntraChainInterface.Depth = myAtoF(userVariables,(parameterStringClass.getString(2)).c_str());          
        if (myIncludeIntraChainInterface.Depth < 0.0000001) {MMBLOG_FILE_FUNC_LINE(CRITICAL, "depth must be >0 !"<<endl); }       
        myBiopolymerClassContainer.validateChainID(myIncludeIntraChainInterface.Chain );         
        for (size_t i = 0 ; i < includeIntraChainInterfaceVector.size() ; i++) {
            if (includeIntraChainInterfaceVector[i].Chain.compare(myIncludeIntraChainInterface.Chain ) == 0) {       
                MMBLOG_FILE_FUNC_LINE(CRITICAL, "The chain "<< myIncludeIntraChainInterface.Chain <<" is already in the includeIntraChainInterfaceVector! Please do not call this twice!"<<endl);
            }        
        }        
        includeIntraChainInterfaceVector.push_back(myIncludeIntraChainInterface); // If myChain survived the tests, add to vector.  This will be dealt with in Repel.cpp:  ConstrainedDynamics::initializeCustomForcesConstraints()          
        return;          
    }
    // SCF: For some reason lines 2290 to 3152 had been duplicated here. Deleted duplicate 4 MAR 2015.

    if  ((parameterStringClass.getString(0)).compare( "physicsInterfaces") == 0) {
        //exit(1);
        MMBLOG_FILE_FUNC_LINE(ALWAYS, "Syntax is:  physicsInterfaces <depth (nm)> <chain>  "<<endl);
        MMBLOG_FILE_FUNC_LINE(ALWAYS, "This adds all residues to the physics zone, that are within <depth> of <chain> in all other chains.  It also sets the bond mobility for residues in <chain> that are within <depth> of that other chain -- so the interface is treated symmetrically."<<endl);
        MMBLOG_FILE_FUNC_LINE(ALWAYS, "You can also specify a multiple-chain set, and find all interfaces between that set and the remainder of chains in the system:"<<endl);
    
        MMBLOG_FILE_FUNC_LINE(ALWAYS, "Syntax is:  physicsInterfaces <depth (nm)> <chain 1> [<chain 2> [<chain 3> [ ...etc]]]  "<<endl);
        MMBLOG_FILE_FUNC_LINE(ALWAYS, "Lastly, you can specify TWO multiple-chain sets, and find only the interfaces between those two sets.  This is particularly useful if you have chains (e.g. threading templates) in your system which need to be ignored:"<<endl);
    
        MMBLOG_FILE_FUNC_LINE(ALWAYS, "Syntax is:  physicsInterfaces <depth (nm)> <chain 1> [<chain 2> [<chain 3> [ ...etc]]]  Versus <chain 1> [<chain 2> [<chain 3> [ ...etc]]]  "<<endl);
        if (parameterStringClass.getString(2).length() == 0) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "Not enough parameters for this command! You must specify at least one chain."<<endl);
        }
        double depth = myAtoF(userVariables,parameterStringClass.getString(1).c_str());
        vector<String> chains;        chains.clear();
        vector<String> partnerChains; partnerChains.clear();
        int i = 2;
        String versusString = String("Versus").toUpper();
        while ((parameterStringClass.getString(i).length() >0)  &&
               (parameterStringClass.getString(i).toUpper().compare(versusString ) != 0 )) { // if parameter i is "versus", stop and start reading partnerChains
            myBiopolymerClassContainer.validateChainID(parameterStringClass.getString(i)); // Make sure this is a valid chain
            chains.push_back( parameterStringClass.getString(i));
            MMBLOG_FILE_FUNC_LINE(INFO, "Added chain "<<parameterStringClass.getString(i)<<" to first set "<<endl);
            i++;
        }
        if (parameterStringClass.getString(i).toUpper().compare(versusString) == 0 ) { 
            i++;
            MMBLOG_FILE_FUNC_LINE(INFO, "Keyword "<<versusString<<" detected."<<endl);
            if (parameterStringClass.getString(i).length() == 0) {
                MMBLOG_FILE_FUNC_LINE(CRITICAL, "The \"versus\" keyword must be followed by at least one chain ID."<<endl);
            }
        }
        while (parameterStringClass.getString(i).length() >0) {
            myBiopolymerClassContainer.validateChainID(parameterStringClass.getString(i)); // Make sure this is a valid chain
            partnerChains.push_back( parameterStringClass.getString(i));
            MMBLOG_FILE_FUNC_LINE(INFO, "Added chain "<<parameterStringClass.getString(i)<<" to second set "<<endl);
            i++;
        }
        
        for (size_t n = 0; n < chains.size(); n++) for (size_t m = 0; m < partnerChains.size(); m++) if (chains[n].compare(partnerChains[m])==0) {
	    MMBLOG_FILE_FUNC_LINE(CRITICAL, "The reference chain "<<chains[n]<<" is the same as partner chain "<<partnerChains[m]<<". This is not kosher!"<<endl);
        }
        physicsContainer.interfaceContainer.addInterface(chains, partnerChains, depth);
        return;
    }


    if (( ((parameterStringClass.getString(0)).compare( "includeAllResiduesWithin") == 0)) ||
        (  ((parameterStringClass.getString(0)).compare( "includeResiduesWithin") == 0))) {

        if ((globalCoulombScaleFactor == 0) && (globalVdwScaleFactor == 0)) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "Please set globalCoulombScaleFactor and/or globalVdwScaleFactor to something other than zero, before issuing this command.  Otherwise, there is no point. You can set safeParameters FALSE to override this."<<endl);
        }   
        if ((!(fabs(globalGbsaScaleFactor) <= 0.0000001) ) && (safeParameters)) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "Please set globalGbsaScaleFactor 0 if you are going to run this command."<<endl);
        }

        MMBLOG_FILE_FUNC_LINE(ALWAYS, " syntax:  includeAllResiduesWithin <radius (nm)> <chain> <residue>    "<<endl);
        MMBLOG_FILE_FUNC_LINE(ALWAYS, "In prior releases of MMB, including 2.10, length was taken in Å -- however we are going back to nm length units, with apologies for the confusion.  "<<endl);
        if (parameterStringClass.getString(2).length() == 0) 
        {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "You have not provided enough parameters for this command."<<endl);
        }
        if (parameterStringClass.getString(4).length() >  0) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "You have provided too many parameters for this command."<<endl);
        }
        {
            //MMBLOG_FILE_FUNC_LINE(" : About to print sequence for chain "<< parameterStringClass.getString(2) <<endl;
            AllResiduesWithin myAllResiduesWithin ; //(parameterStringClass.getString(2), myBiopolymerClassContainer.residueID(userVariables,parameterStringClass.getString(3),myAllResiduesWithin.chain ), myAtoF(userVariables,parameterStringClass.getString(1).c_str()));
            myAllResiduesWithin.setRadius   ( myAtoF(userVariables,parameterStringClass.getString(1).c_str()));
            myAllResiduesWithin.setChain   ( parameterStringClass.getString(2));
            if (myBiopolymerClassContainer.hasChainID(myAllResiduesWithin.getChain())) {
                myAllResiduesWithin.setResidue (  myBiopolymerClassContainer.residueID(userVariables,parameterStringClass.getString(3),myAllResiduesWithin.getChain() ));
		if (parameterStringClass.getString(3).length() == 0) 
		{
		    MMBLOG_FILE_FUNC_LINE(CRITICAL, "You have not provided enough parameters for this command."<<endl);
		}
            }
            else if (moleculeClassContainer.hasChainID(myAllResiduesWithin.getChain())) {
                myAllResiduesWithin.setResidue   ( ResidueID(1,' '));//(userVariables,parameterStringClass.getString(3) ));
            }
            includeAllResiduesWithinVector.push_back(myAllResiduesWithin);             }
            return;
    }

    if ( ((parameterStringClass.getString(0)).compare( "includeNonBondAtomInBiopolymer") == 0) )
    {
        cout<<"syntax:  includeNonBondAtomInBiopolymer chain residue atomName   "<<endl;
        if ((globalCoulombScaleFactor == 0) && (globalVdwScaleFactor == 0)) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "Please set globalCoulombScaleFactor and/or globalVdwScaleFactor to something other than zero, before issuing this command.  Otherwise, there is no point. You can set safeParameters FALSE to override this."<<endl);
        }   
        if ((!(fabs(globalGbsaScaleFactor) <= 0.0000001) ) && (safeParameters)) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "Please set globalGbsaScaleFactor 0 if you are going to run this command."<<endl);
        }   
        if (parameterStringClass.getString(3).length() == 0) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "You have not provided enough parameters for this command."<<endl);
        }
        if (parameterStringClass.getString(4).length() >  0) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "You have provided too many parameters for this command."<<endl);
        }
        {
            IncludeNonBondAtomInBiopolymerStruct myIncludeNonBondAtomInBiopolymerStruct;
            myIncludeNonBondAtomInBiopolymerStruct.chain   = parameterStringClass.getString(1);
            myIncludeNonBondAtomInBiopolymerStruct.residue   =myBiopolymerClassContainer.residueID(userVariables,(parameterStringClass.getString(2)),myIncludeNonBondAtomInBiopolymerStruct.chain  );
            myIncludeNonBondAtomInBiopolymerStruct.atomName   = parameterStringClass.getString(3);
            String myAtomPath = myBiopolymerClassContainer.updBiopolymerClass( myIncludeNonBondAtomInBiopolymerStruct.chain).atomPathString(myIncludeNonBondAtomInBiopolymerStruct.residue, myIncludeNonBondAtomInBiopolymerStruct.atomName); // checking that atom exists
            myBiopolymerClassContainer.updBiopolymerClass( myIncludeNonBondAtomInBiopolymerStruct.chain).validateAtomPathName(myAtomPath);
            includeNonBondAtomInBiopolymerVector.push_back(myIncludeNonBondAtomInBiopolymerStruct);             
        }

            return;
    }
    if ((parameterStringClass.getString(0)).compare("nastGlobalBondTorsionScaleFactor") ==0) {
        parameterStringClass.validateNumFields(2);
        nastGlobalBondTorsionScaleFactor = myAtoF(userVariables,(parameterStringClass.getString(1)).c_str());    
        return;
    }
    if ((parameterStringClass.getString(0)).compare("noseHooverTime") ==0) {
        parameterStringClass.validateNumFields(2);
        noseHooverTime = myAtoF(userVariables,(parameterStringClass.getString(1)).c_str());    

        return;
    }
    if ( ((parameterStringClass.getString(0)).compare( "physicsWhereYouWantIt") == 0)) {
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "This parameter is obsolete!  You do not need to set it explicitly anymore.  It is set implicitly when you issue any command to add atoms to the force field.  For example, includeAllResiduesWithin or includeAllNonBondAtomsInResidues will do this."<<endl);
        parameterStringClass.validateNumFields(2);
        //physicsWhereYouWantIt = aToBool(parameterStringClass.getString(0),(parameterStringClass.getString(1)).c_str());
        if (safeParameters ) 
            if ((globalCoulombScaleFactor == 0) && (globalVdwScaleFactor == 0)) {
                MMBLOG_FILE_FUNC_LINE(CRITICAL, "Please set globalCoulombScaleFactor and/or globalVdwScaleFactor to something other than zero, before issuing this command.  Otherwise, there is no point. You can set safeParameters FALSE to override this."<<endl);
            }   
        if ((!(fabs(globalGbsaScaleFactor) <= 0.0000001) ) && (safeParameters)) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "Please set globalGbsaScaleFactor 0 if you are going to run this command."<<endl);
        }   
        return;
    }
    if (((parameterStringClass.getString(0)).compare("addHuntCrossleySpheres") ==0)||((parameterStringClass.getString(0)).compare("addLargeHuntCrossleySpheres") ==0) ||((parameterStringClass.getString(0)).compare("addHardSphere") ==0) ) {
        parameterStringClass.validateNumFields(2);
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "The commands addHuntCrossleySpheres, addLargeHuntCrossleySpheres, and addHardSphere are no longer supported."<<endl);
        return;
    }
    if ((parameterStringClass.getString(0)).compare("addTestSpring") ==0) {
        parameterStringClass.validateNumFields(2);
        addTestSpring               = aToBool(parameterStringClass.getString(0),(parameterStringClass.getString(1)).c_str());    
        return;
    }
    if ((parameterStringClass.getString(0)).compare("applyC1pSprings") ==0) {
        parameterStringClass.validateNumFields(2);
        applyC1pSprings             = aToBool(parameterStringClass.getString(0),(parameterStringClass.getString(1)).c_str());    
        return;
    }
    if ((parameterStringClass.getString(0)).compare("stackAllHelicalResidues") ==0) {
        parameterStringClass.validateNumFields(2);
        stackAllHelicalResidues = aToBool(parameterStringClass.getString(0),  (parameterStringClass.getString(1)).c_str());    
        return;
    }
    if ((parameterStringClass.getString(0)).compare("thermostatType"         ) ==0) {
        parameterStringClass.validateNumFields(2);
        thermostatType          = ((parameterStringClass.getString(1)));    
        return;
    }
    if ((parameterStringClass.getString(0)).compare("setForceAndStericScrubber") ==0) {

        MMBLOG_FILE_FUNC_LINE(CRITICAL, "This is no longer a user configurable parameter.  setForceAndStericScrubber is set to TRUE automatically by setting dutyCycle to anythin other than 1.  "<<endl);
        parameterStringClass.validateNumFields(2);
        setForceAndStericScrubber = aToBool(parameterStringClass.getString(0),(parameterStringClass.getString(1)).c_str());    
        if ( setForceAndStericScrubber  &&
                (densityContainer.numDensityStretches()>0 || electroDensityContainer.numDensityStretches()>0)) {
            // MMBLOG_FILE_FUNC_LINE(": You can't use fitToDensity if you are also going to issue setForceAndStericScrubber TRUE "<<parameterFileName<<endl;
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "You can't use fitToDensity if you are also going to issue setForceAndStericScrubber TRUE "<<endl);
        }
        return;
    }
    if ((parameterStringClass.getString(0)).compare("setForceScrubber") ==0) {
        parameterStringClass.validateNumFields(2);
        //MMBLOG_FILE_FUNC_LINE(": setForceScrubber is no longer a supported parameter "<<parameterFileName<<endl;
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "setForceScrubber is no longer a supported parameter "<<endl);
        setForceScrubber = aToBool(parameterStringClass.getString(0),(parameterStringClass.getString(1)).c_str());    
        return;
    }
    if ((parameterStringClass.getString(0)).compare("setHelicalStacking") ==0) {
        parameterStringClass.validateNumFields(2);
        setHelicalStacking = aToBool(parameterStringClass.getString(0),(parameterStringClass.getString(1)).c_str());    
        return;
    }
    if (((parameterStringClass.getString(0)).compare("setInitialVelocities") ==0) || ((parameterStringClass.getString(0)).compare("randomizeInitialVelocities") ==0)) {
        parameterStringClass.validateNumFields(2);
        setInitialVelocities = aToBool(parameterStringClass.getString(0),(parameterStringClass.getString(1)).c_str());    
        return;
    }
    if ((parameterStringClass.getString(0)).compare("setRemoveBasePairsInRigidStretch") ==0) {
        parameterStringClass.validateNumFields(2);
        setRemoveBasePairsInRigidStretch = aToBool(parameterStringClass.getString(0),(parameterStringClass.getString(1)).c_str());    
        return;
    }
    if ((parameterStringClass.getString(0)).compare("setRemoveBasePairsAcrossRigidStretches") ==0) {
        parameterStringClass.validateNumFields(2);
        setRemoveBasePairsAcrossRigidStretches = aToBool(parameterStringClass.getString(0),(parameterStringClass.getString(1)).c_str());    
        return;
    }
    if ((parameterStringClass.getString(0)).compare("setTemperature") ==0) {
        parameterStringClass.validateNumFields(2);
        setTemperature = aToBool(parameterStringClass.getString(0),(parameterStringClass.getString(1)).c_str());    
        return;

    } //else  SCF
    // KLUDGE -- should use the if-continue method for all above
    if ((parameterStringClass.getString(0)).compare("smallGroupInertiaMultiplier") ==0) {
        parameterStringClass.validateNumFields(2);
        smallGroupInertiaMultiplier = myAtoF   (userVariables,(parameterStringClass.getString(1)).c_str());    
        return;
    }
    if ((parameterStringClass.getString(0)).compare("reportingInterval") ==0) {
        MMBLOG_FILE_FUNC_LINE(ALWAYS, "This command sets the simulation time per reporting interval, in ps. Total simulation time = reportingInterval * numReportingIntervals.." << endl);
        MMBLOG_FILE_FUNC_LINE(ALWAYS, "Syntax: reportingInterval <time (ps)> " << endl);
        MMBLOG_FILE_FUNC_LINE(ALWAYS, endl);
        parameterStringClass.print();

        MMBLOG_FILE_FUNC_LINE(ALWAYS, endl);
        parameterStringClass.validateNumFields(2);
        reportingInterval= myAtoF(userVariables,(parameterStringClass.getString(1)).c_str());    
        return;
    }
    if ((parameterStringClass.getString(0)).compare("restrainingForceConstant") ==0) {
        parameterStringClass.validateNumFields(2);
        restrainingForceConstant= myAtoF(userVariables,(parameterStringClass.getString(1)).c_str());    
        return;
    }
    if ((parameterStringClass.getString(0)).compare("restrainingTorqueConstant") ==0) {
        parameterStringClass.validateNumFields(2);
        restrainingTorqueConstant= myAtoF(userVariables,(parameterStringClass.getString(1)).c_str());    
        //} else if ((parameterStringClass.getString(0)).compare("resetBases"       ) ==0) {
        //    resetBases       = myAtoI(userVariables,(parameterStringClass.getString(1)).c_str());    
        return;
    }
    if ((parameterStringClass.getString(0)).compare("scrubberPeriod") ==0) {
        scrubberPeriod = myAtoF(userVariables,(parameterStringClass.getString(1)).c_str());    
        return;
    }
    if (((parameterStringClass.getString(0)).compare("maxReportingIntervals") ==0) || ((parameterStringClass.getString(0)).compare("numReportingIntervals") ==0)) {
        MMBLOG_FILE_FUNC_LINE(ALWAYS, "This command sets the number of reporting intervals, an integer. Total simulation time = reportingInterval * numReportingIntervals.." << endl);
        MMBLOG_FILE_FUNC_LINE(ALWAYS, "Syntax: numReportingIntervals <number of intervals> " << endl);
        parameterStringClass.validateNumFields(2);
        numReportingIntervals = myAtoI(userVariables,(parameterStringClass.getString(1)).c_str());    

        GlobalProgressWriter::get().setTotalSteps(numReportingIntervals);
        return;
    }
    if ((parameterStringClass.getString(0)).compare("leontisWesthofInFileName") ==0) {
        if (safeParameters) {
	    MMBLOG_FILE_FUNC_LINE(CRITICAL, "leontisWesthofInFileName is not a supported parameter at the moment. Please just have a Leontis-Westhof parameter file called parameters.csv in your working directory. "<<endl);
        }
        parameterStringClass.validateNumFields(2);
        leontisWesthofInFileName = parameterStringClass.getString(1);
        if (verbose) MMBLOG_FILE_FUNC_LINE(INFO, "leontisWesthofInFileName ="<<leontisWesthofInFileName<<endl);
        return;
    }
    if ((parameterStringClass.getString(0)).compare("rigidifyTermini") ==0) {
        parameterStringClass.validateNumFields(2);
        // MMBLOG_FILE_FUNC_LINE(": rigidifyTermini is no longer a supported parameter "<<parameterFileName<<endl;
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "rigidifyTermini is no longer a supported parameter "<<endl);
        //rigidifyTermini = aToBool(parameterStringClass.getString(0),(parameterStringClass.getString(1)).c_str());    
        return;
    }
    if ((parameterStringClass.getString(0)).compare("integratorStepSize") ==0) {
        parameterStringClass.validateNumFields(2);
        integratorStepSize                    = myAtoF(userVariables,(parameterStringClass.getString(1)).c_str());    
        return;
    }
    if ((parameterStringClass.getString(0)).compare("useFixedStepSize") ==0) {
        parameterStringClass.validateNumFields(2);
        useFixedStepSize            = aToBool(parameterStringClass.getString(0),(parameterStringClass.getString(1)).c_str());    
        return;
    }
    if ((parameterStringClass.getString(0)).compare("useMultithreadedComputation") ==0) {
        parameterStringClass.validateNumFields(2);
        useMultithreadedComputation =  aToBool(parameterStringClass.getString(0),(parameterStringClass.getString(1)).c_str());    
        return;
    }
    if ((parameterStringClass.getString(0)).compare("useOpenMMAcceleration") ==0) {
        parameterStringClass.validateNumFields(2);
        useOpenMMAcceleration =  aToBool(parameterStringClass.getString(0),(parameterStringClass.getString(1)).c_str());    
        return;
    } 
    if ((parameterStringClass.getString(0)).compare("matchHydrogenAtomLocations") ==0) {
        parameterStringClass.validateNumFields(2);
        matchHydrogenAtomLocations = aToBool(parameterStringClass.getString(0),(parameterStringClass.getString(1)).c_str());
        return;
    }
 
    if ((parameterStringClass.getString(0)).compare("matchPurineN1AtomLocations") ==0) {
        parameterStringClass.validateNumFields(2);
        matchPurineN1AtomLocations = aToBool(parameterStringClass.getString(0),(parameterStringClass.getString(1)).c_str());
        return;
    }

    /*if ((parameterStringClass.getString(0)).compare("matchProteinCarboxylOxygenLocations") ==0) {
      parameterStringClass.validateNumFields(2);
      matchProteinCarboxylOxygenLocations = aToBool(parameterStringClass.getString(0),(parameterStringClass.getString(1)).c_str());
      return;
      }*/
    if ((parameterStringClass.getString(0)).compare("ignoreAtomLocation") ==0) {
        MMBLOG_FILE_FUNC_LINE(ALWAYS, "This command tells MMB to not match the atom coordinates for a specified atom. Use this, for example, if there is a certain atom engaged in some wonky planarity/nonplanarity."<<endl);
        MMBLOG_FILE_FUNC_LINE(ALWAYS, "At the moment, this works only for biopolymers."<<endl);
        MMBLOG_FILE_FUNC_LINE(ALWAYS, "Syntax: ignoreAtomLocation <chain> <residue ID> <atom name> ." << endl);
        parameterStringClass.validateNumFields(4);
        String myChain = (parameterStringClass.getString(1)); 
        myBiopolymerClassContainer.validateChainID(myChain);
        ResidueID myResidueID = myBiopolymerClassContainer.updBiopolymerClass(myChain).residueID(parameterStringClass.getString(2));
        String myAtom  = (parameterStringClass.getString(3)); 
        myBiopolymerClassContainer.updBiopolymerClass(myChain).validateAtomPathName(  myBiopolymerClassContainer.updBiopolymerClass(myChain). atomPathString(myResidueID,myAtom)); 
        ResidueInfo::Index myResidueIndex = myBiopolymerClassContainer.updBiopolymerClass(myChain).getResidueIndex( myResidueID); 
        MMBAtomInfo atomPositionToIgnore(myChain,myResidueID,myResidueIndex, myAtom);
        myBiopolymerClassContainer.updBiopolymerClass(myChain).addAtomPositionToIgnore(atomPositionToIgnore );
        return;
    }
    if ((parameterStringClass.getString(0)).compare("matchingMinimizerTolerance") ==0) {
        parameterStringClass.validateNumFields(2);
        matchingMinimizerTolerance = myAtoF(userVariables,(parameterStringClass.getString(1)).c_str());
        return;
    }



    if ((parameterStringClass.getString(0)).compare("guessCoordinates") ==0) {
        if (safeParameters) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "This parameter is no longer user configurable.  Use the matchFast, matchGapped, or matchGappedNoHeal macros instead. "<<endl);
        }
        parameterStringClass.validateNumFields(2);
        guessCoordinates = aToBool(parameterStringClass.getString(0),(parameterStringClass.getString(1)).c_str());
        return;
    }
    if ((parameterStringClass.getString(0)).compare("matchOptimize") ==0) {

        if (safeParameters) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "This parameter is no longer user configurable.  Use the matchFast, matchGapped, or matchGappedNoHeal macros instead. "<<endl);
        }
        parameterStringClass.validateNumFields(2);
        matchOptimize      =  aToBool(parameterStringClass.getString(0),(parameterStringClass.getString(1)).c_str());    
        return;
    }
    if ((parameterStringClass.getString(0)).compare("matchExact") ==0) {

        if (safeParameters) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "This parameter is no longer user configurable.  Use the matchFast, matchGapped, or matchGappedNoHeal macros instead. "<<endl);
        }
        parameterStringClass.validateNumFields(2);
        matchExact      =  aToBool(parameterStringClass.getString(0),(parameterStringClass.getString(1)).c_str());    

        if (! (matchExact || matchIdealized)) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "You must set matchExact and/or matchIdealized to true. Right now both are false. "<<endl);
        }

        return;
    }
    if ((parameterStringClass.getString(0)).compare("matchIdealized") ==0) 
    {

        if (safeParameters) 
        {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "This parameter is no longer user configurable.  Use the matchFast, matchGapped, or matchGappedNoHeal macros instead. "<<endl);
        }
        parameterStringClass.validateNumFields(2);
        matchIdealized      =  aToBool(parameterStringClass.getString(0),(parameterStringClass.getString(1)).c_str());    

        if (! (matchExact || matchIdealized)) 
        {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "You must set matchExact and/or matchIdealized to true. Right now both are false. "<<endl);
        }
        return;
    }

    if ((parameterStringClass.getString(0)).compare("matchFast") ==0) 
    {
        parameterStringClass.validateNumFields(1);
        MMBLOG_FILE_FUNC_LINE(ALWAYS, "This is a macro which sets matchExact to true,  matchIdealized to false.  This will fit to all atoms provided in the input structure file, but will not be able to handle any missing backbone atoms between regions of known structure.  Use this when there are no gaps in your known structure.  missing residues or atoms at the beginning and end of the known region are OK; and excess length will be instantiated as an extended chain.  This is the fastest fitting protocol -- but BEWARE:  you should only do this if you have an MMB-generated double-precision input PDB file (which contains REMARK-SIMTK-COORDS records).  Otherwise the error is likely to accumulate over any chain of nontrivial length."<<endl);

        matchExact         =  true ; 
        matchIdealized     =  false; 
        matchOptimize      =  false; 
    guessCoordinates   =  false;
        return;
    }

    if ((parameterStringClass.getString(0)).compare("matchNoGaps") ==0) 
    {
        parameterStringClass.validateNumFields(1);
        MMBLOG_FILE_FUNC_LINE(ALWAYS, "This is a macro which sets matchExact, matchIdealized to true.  This will fit to all atoms provided in the input structure file, and will guess atom positions for all fragments missing from that file, consistent with default bond lengths and angles.  You should generally use this macro when you have atoms missing along the backbone.   "<<endl);

        matchExact         =  true ; 
        matchIdealized     =  true ; 
        matchOptimize      =  false; 
        guessCoordinates   =  false;  
        return;
    }
    if ((parameterStringClass.getString(0)).compare("matchGapped") ==0) 
    {
        parameterStringClass.validateNumFields(1);
        MMBLOG_FILE_FUNC_LINE(ALWAYS, "This is a macro which sets matchExact, matchIdealized, matchOptimize, and guessCoordinates to true.  This will fit to all atoms provided in the input structure file, and will guess atom positions for all fragments missing from that file, consistent with default bond lengths and angles.  You should generally use this macro when you have atoms missing along the backbone. The missing atoms can be in the middle of a chain, although unnatural bond geometries may occur as MMB does its darnedest to achieve loop closure and match all it can.   "<<endl);

        matchExact         =  true ; 
        matchIdealized     =  true ; 
        matchOptimize      =  true ; 
        guessCoordinates   =  true;
        return;
    }

    if ((parameterStringClass.getString(0)).compare("matchGappedNoHeal") ==0) 
    {
        if (safeParameters) {MMBLOG_FILE_FUNC_LINE(CRITICAL, "This macro is not currently supported!"<<endl);}
        MMBLOG_FILE_FUNC_LINE(ALWAYS, "This is a macro which sets matchExact to true, and matchIdealized to false.  This will fit to all atoms provided in the input structure file, and will guess atom positions for all fragments missing from that file.  The connection between known and unknown fragments will in general have an unphysical bond geometry.  Use this in unusual cases where physically reasonable connection between fragments of known structure is impossible. Use matchGapped for cases when such a connection IS possible.    "<<endl);

        parameterStringClass.validateNumFields(1);
        matchExact         =  true ; 
        matchIdealized     =  false; 
        matchOptimize      =  false;       
        guessCoordinates   =  true ;
        return;
    }
    if ((parameterStringClass.getString(0)).compare("matchPerfect") ==0) 
    {
      
        if (safeParameters) {
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "This macro is no longer user configurable.  Use the matchFast, matchGapped, or matchGappedNoHeal macros instead. "<<endl);}
        MMBLOG_FILE_FUNC_LINE(ALWAYS, "This is a macro.  If True, it sets BOTH matchExact and matchIdealized to True.  This means that all bond lengths, angles, and dihedrals will be matched locally, and then there will be a global refinement of torsion angles to correct for accumulated error.  This costs as much as matchIdealized, but generally gives better results."<<endl);

        matchExact          =  bool(true);    
        matchIdealized      =  bool(true);    
        matchOptimize      =   true;  
        guessCoordinates    =  bool(true);
        if (parameterStringClass.getString(1).length() >0){
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "You have provided too many parameters! matchPerfect is a macro which takes no parameters."<<endl);
        }
        return;
    }
    if ((parameterStringClass.getString(0)).compare("minimize") ==0) {
        parameterStringClass.validateNumFields(2);
        //MMBLOG_FILE_FUNC_LINE(CRITICAL, " minimize is obsolete."<<endl;
        minimize = aToBool(parameterStringClass.getString(0),(parameterStringClass.getString(1)).c_str()); //aToBool(userVariables,(parameterStringClass.getString(1)).c_str());    
        return;
    }
    if ((parameterStringClass.getString(0)).compare("monteCarloRun") ==0) {
        parameterStringClass.validateNumFields(2);
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "The monteCarloRun is obsolete."<<endl);
        monteCarloRun = myAtoI(userVariables,(parameterStringClass.getString(1)).c_str());    
        return;
    }
    if ((parameterStringClass.getString(0)).compare("monteCarloTemperature") ==0) {
        parameterStringClass.validateNumFields(2);
        monteCarloTemperature = myAtoF(userVariables,(parameterStringClass.getString(1)).c_str()); 
        return;
    }
    if ((parameterStringClass.getString(0)).compare("monteCarloTemperatureIncrement") ==0) {
        parameterStringClass.validateNumFields(2);
        monteCarloTemperatureIncrement = myAtoF(userVariables,(parameterStringClass.getString(1)).c_str());    
        return;
    }
    if ((parameterStringClass.getString(0)).compare("setLoopBondMobility") ==0) {
        parameterStringClass.validateNumFields(2);
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "This is no longer a supported parameter.   "<<endl);
        setLoopBondMobility =  aToBool(parameterStringClass.getString(0),(parameterStringClass.getString(1)).c_str());    
        return;
    }
    if ((parameterStringClass.getString(0)).compare("loopBondMobility") ==0) {
        parameterStringClass.validateNumFields(2);
        loopBondMobility = BondMobility::Mobility(myAtoI(userVariables,(parameterStringClass.getString(1)).c_str())); 
        return;
    }
    if ((parameterStringClass.getString(0)).compare("helixBondMobility") ==0) {
        parameterStringClass.validateNumFields(2);
        helixBondMobility = BondMobility::Mobility(myAtoI(userVariables,(parameterStringClass.getString(1)).c_str()));  
        return;
    }
    if ((parameterStringClass.getString(0)).compare("overallBondMobility") ==0) {
        parameterStringClass.validateNumFields(2);
        overallBondMobility = BondMobility::Mobility(myAtoI(userVariables,(parameterStringClass.getString(1)).c_str()));
        return;
    }
    if ((parameterStringClass.getString(0)).compare(    "chiBondMobility") ==0) {
        parameterStringClass.validateNumFields(2);
        chiBondMobility = BondMobility::Mobility(myAtoI(userVariables,(parameterStringClass.getString(1)).c_str())); 
        return;
    }
    if ((parameterStringClass.getString(0)).compare("") ==0) {
        parameterStringClass.validateNumFields(2);
        if (verbose) cout << __FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" blank line in input file."<<parameterStringClass.getString(0)<<endl;
        return;
    }

    if(parameterStringClass.getString(0).compare("detectConvergence")==0)
    {
        MMBLOG_FILE_FUNC_LINE(ALWAYS, "detectConvergence: stop the stage if the simulation converged, i.e. if the difference in total energy for the few last frames (default 5) is inferior to an epsilon value (default 0.5 kJ/mol)"<< endl);
        parameterStringClass.validateNumFields(2);
        detectConvergence = aToBool(parameterStringClass.getString(0),(parameterStringClass.getString(1)).c_str());    
        return;
    }

    if(parameterStringClass.getString(0).compare("convergenceTimeout")==0)
    {
        parameterStringClass.validateNumFields(2);
        convergenceTimeout = myAtoI(userVariables,(parameterStringClass.getString(1)).c_str());    
        return;
    }
    if(parameterStringClass.getString(0).compare("convergenceEpsilon")==0)
    {
        parameterStringClass.validateNumFields(2);
        convergenceEpsilon = myAtoF(userVariables,(parameterStringClass.getString(1)).c_str());    
        return;
    }

    if ((parameterStringClass.getString(0)).compare("loggingSeverity") == 0) {
        parameterStringClass.validateNumFields(2);

        std::string param = parameterStringClass.getString(1);
        std::transform(param.begin(), param.end(), param.begin(), [](auto &&ch){ return std::tolower(ch); });

        if (param == "debug")
            MMBLogger::instance().setLoggingSeverity(MMBLogger::Severity::DEBUG);
        else if (param == "info")
            MMBLogger::instance().setLoggingSeverity(MMBLogger::Severity::INFO);
        else if (param == "warning")
            MMBLogger::instance().setLoggingSeverity(MMBLogger::Severity::WARNING);
        else {
            MMBLOG_FILE_FUNC_LINE(
                CRITICAL,
                "Unknown logging severity "<<param<<endl
                <<"Allowed values are: debug, info, warning"<<endl
            );
        }

        return;
    }

    // if non recognizable string
    MMBLOG_FILE_FUNC_LINE(CRITICAL, "Non recognizable input: " << parameterStringClass.getString() << endl);
}



void ParameterReader::initializeFromFileOnly(const char * parameterFileName ) {
    basePairContainer.printBasePairs();
    MMBLOG_FILE_FUNC_LINE(INFO, endl);
    ifstream inFile(parameterFileName,ifstream::in);
    if (!(inFile.good())){
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "Could not open file "<<parameterFileName<<endl);
    }
    int readStage = 0;
    // added to initializeDefaults() / clearContainers().  removing from here.
    bool readAtOneStageOnly = false;                  
    bool readOnlyUntilStage = false;                  
    bool readExcept         = false;                  

    userVariables[String("@CURRENTSTAGE")] = currentStage;
    MMBLOG_FILE_FUNC_LINE(INFO, "Just set variable @CURRENTSTAGE to "<<userVariables[String("@CURRENTSTAGE")]<<endl);

    while (inFile.good()) {
        String tempString;
        getline(inFile, tempString);
        ParameterStringClass parameterStringClass( tempString );
        //parameterStringClass.print();
        // now start detecting and dealing with parameter flags
        if (((parameterStringClass.getString(0)).compare("readFromStage") == 0) || ((parameterStringClass.getString(0)).compare("readAtStage") == 0) || ((parameterStringClass.getString(0)).compare("readToStage") == 0) || ((parameterStringClass.getString(0)).compare("readExceptAtStage") == 0)  ) 
        {
            SimTK_ERRCHK_ALWAYS(
                    ((readStage == 0)), 
                    "[ParameterReader.cpp]", "You are not allowed to nest readFromStage -- readBlockEnd or readAtStage -- readBlockEnd blocks.");
            readStage =myAtoI(userVariables,(parameterStringClass.getString(1)).c_str()); 
            MMBLOG_FILE_FUNC_LINE(INFO, "This statement applies to stage readStage = "<<readStage<<endl);

            //MMBLOG_FILE_FUNC_LINE(": "<<(parameterStringClass.getString(0)).compare("readAtStage")<<endl;
            if ((parameterStringClass.getString(0)).compare("readAtStage") == 0)
            {
                readAtOneStageOnly = true; 
                MMBLOG_FILE_FUNC_LINE(INFO, "readAtOneStageOnly: "<<readAtOneStageOnly<<endl);
            }
            else if ((parameterStringClass.getString(0)).compare("readToStage") == 0)
            {readOnlyUntilStage = true;}
            else if ((parameterStringClass.getString(0)).compare("readExceptAtStage") == 0)
            {
                readExcept = true;
            }
        }
        else if (((parameterStringClass.getString(0)).compare("readBlockEnd") == 0)) {
            readStage = 0; 
            readAtOneStageOnly = false; 
            readOnlyUntilStage = false; 
            readExcept = false;
        }
        else if (((readStage < currentStage ) &&  !(readAtOneStageOnly ) && !(readOnlyUntilStage)) 
                || ((readStage > currentStage ) && (readOnlyUntilStage || readExcept) )  
                || ((readStage == currentStage) && !(readExcept)  ))  // if current stage is higher than or equal to readStage in the case of a readFromStage block, or equal in the case of a readAtStage block.
        {
            MMBLOG_FILE_FUNC_LINE(INFO, "read line: "<<tempString.c_str()<<endl);
            parameterStringInterpreter(parameterStringClass, readStage, readAtOneStageOnly, readOnlyUntilStage, readExcept);
        }

    } // end of while loop

    if (safeParameters && removeRigidBodyMomentum) {
        for (int i = 0; i < atomSpringContainer.numAtomSprings(); i++){
            SimTK_ERRCHK_ALWAYS(
                    (!atomSpringContainer.getAtomSpring(i).toGround  ), 
                    "[ParameterReader.cpp]",                    "springToGround and tetherToGround probably won't work properly unless you set removeRigidBodyMomentum False .     ");
        }
    }

    MMBLOG_FILE_FUNC_LINE(DEBUG, "base pairs just after reading from commands.dat:"<<endl);
    //if (verbose) 
    basePairContainer.printBasePairs();
    MMBLOG_FILE_FUNC_LINE(INFO, "After reading command file, here is the list of atom springs: "<<endl);
    atomSpringContainer.printAtomSprings();
    MMBLOG_FILE_FUNC_LINE(INFO, "Done printing atom springs. "<<endl);
    //printBiopolymerSequenceInfo(myBiopolymerClassContainer.updBiopolymerClass("g").myBiopolymer);

} ;


void ParameterReader::loadSequencesFromPdb(const char * pdbFileName, const string & chainsPrefix, const bool tempRenumberPdbResidues ) { //, vector<std::string> deletedResidueVector ){
    MMBLOG_FILE_FUNC_LINE(INFO, "This command will now extract all DNA, RNA, and protein sequences found in the input structure file, "<<pdbFileName<<" and instantiate the corresponding biopolymers. "<<endl);

    if (densityContainer.numDensityStretches() > 0) {
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "You have already declared "<<densityContainer.numDensityStretches()<<" biopolymer stretches to be fitted to the density map.  Please do this after you have created ALL biopolymer chains."<<endl);
    }
    myBiopolymerClassContainer.loadSequencesFromPdb( pdbFileName, proteinCapping, chainsPrefix, tempRenumberPdbResidues, useNACappingHydroxyls  );
    //myBiopolymerClassContainer.loadResidueIDVector();  // is now being done by setResidueIDsAndInsertionCodesFromBiopolymer
    return;
}


void ParameterReader::setFirstAndLastStageAndUseCifFiles(const char * parameterFileName ) {
//void ParameterReader::setFirstAndLastStage(const char * parameterFileName ) {
    MMBLOG_FILE_FUNC_LINE(INFO, endl);
    ifstream inFile(parameterFileName,ifstream::in);
    if (! inFile.good()){
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "Unable to open command file: "<<parameterFileName<<endl);
    }
    const int numFields = 10 ;
    String mystring[numFields];//1;   
    //int readStage = 0;
    baseOperationVector.clear();
    basePairContainer.clear();
    contactContainer.clear();
    densityContainer.clear();
    electroDensityContainer.clear();
    myMonoAtomsContainer.clear();
    myBiopolymerClassContainer.clear();
    waterDropletContainer.clear();
    //bool readAtOneStageOnly = false;                  
    //bool readOnlyUntilStage = false;                  
    //bool readExcept         = false;                  
    // make sure firstStage and lastStage are initialized
    firstStage = 0;
    lastStage = 0;
    useCIFFileFormat = 0;
    while (inFile.good()) {
        for (int w = 0;w<numFields;w++) mystring[w]="";
        String tempString;
        getline(inFile, tempString);
        stringstream u;
        u.str(tempString); 
        MMBLOG_FILE_FUNC_LINE(DEBUG, "read line: "<<u.str()<<endl);
        u>>mystring[0]>>mystring[1]>>mystring[2]>>mystring[3]>>mystring[4]>>mystring[5]>>mystring[6]>>mystring[7]>>mystring[8]>>mystring[9]; 
        // now start detecting and dealing with parameter flags

        if  (((mystring[0]).compare("-FS") == 0) || ((mystring[0]).compare("firstStage") == 0))  {


            size_t position  = String(mystring[1] ).find_last_of('@');
            if (    position != String::npos) { 
                MMBLOG_FILE_FUNC_LINE(WARNING, "User variables are not permitted for setting firstStage or lastStage. "<<endl);
            }


            firstStage = myAtoI(userVariables,(mystring[1]).c_str());     
        } // of if
        else if (((mystring[0]).compare("-LS") == 0) || ((mystring[0]).compare("lastStage") == 0))  {
            size_t position  = String(mystring[1] ).find_last_of('@');
            if (    position != String::npos) { 
                MMBLOG_FILE_FUNC_LINE(WARNING, "User variables are not permitted for setting firstStage or lastStage. "<<endl);
            }
            lastStage = myAtoI(userVariables,(mystring[1]).c_str());     
        }
	else if ( ( mystring[0]).compare("useCIFFiles") == 0 )
	{
	    //============================================ Set the CIF files usage
	    //parameterStringClass.validateNumFields        ( 2 );
	    useCIFFileFormat                        = aToBool ( mystring[0], (mystring[1]).c_str() );
            MMBLOG_FILE_FUNC_LINE(INFO, "useCIFFileFormat now set to "<<useCIFFileFormat<<endl);
	}
    } // of while inFile 
} ;



void ParameterReader::postInitialize(){
    MMBLOG_FILE_FUNC_LINE(INFO, "Printing original and renumbered residue numbers: "<<endl);
    myBiopolymerClassContainer.printOriginalAndRenumberedResidueIDs();
    MMBLOG_FILE_FUNC_LINE(INFO, endl);
    MMBLOG_FILE_FUNC_LINE(INFO, "In ParameterReader::postInitialize, printing out density stretches: "<<endl);
    //densityContainer.printDensityStretches    (); 
    //densityContainer.stuffDensityStretchVector(myBiopolymerClassContainer); 
    densityContainer.printDensityStretches    (); 
    if (safeParameters) {
        MMBLOG_FILE_FUNC_LINE(ALWAYS, ".1 - integratorAccuracy "<<.1f - integratorAccuracy<<endl);
        SimTK_ERRCHK_ALWAYS( ( (integratorType.compare("RungeKuttaMerson") == 0  ))    ,"[ParameterReader.cpp]","RungeKuttaMerson is the only integratorType that conserves angular momentum.  You have selected something different.  To override this message, set safeParameters False"      );
        //if you don't think this is necessary try running P5abc-softspheres, see if it works.  Maybe this was not a condition we needed to check.  I suspect root cause had to do with Torsions.
        SimTK_ERRCHK1_ALWAYS( !(((integratorAccuracy - 0.001) > .000000001  ) && (integratorType.compare("Verlet") == 0  ))    ,"[ParameterReader.cpp]","If you choose integratorType Verlet you should set integratorAccuracy <= .001. You are using: %f",integratorAccuracy);
        SimTK_ERRCHK_ALWAYS(((potentialType.compare("HarmonicInverse") == 0 )||(potentialType.compare("HarmonicInverseLorentzian") == 0 ) )    ,"[ParameterReader.cpp]","Make sure you use a supported potentialType, e.g. HarmonicInverse or HarmonicInverseLorentzian."       );
        SimTK_ERRCHK1_ALWAYS(fabs(cutoffRadius      -  0.10f) < .000001 ,"[ParameterReader.cpp]","Making sure you are using cutoffRadius = .1 nm.  You are using: %f nm.  If you want to override this, set safeParameters False.",cutoffRadius      );
        SimTK_ERRCHK1_ALWAYS((integratorAccuracy-  0.10f) < .000001 ,"[ParameterReader.cpp]","Making sure your integratorAccuracy of %f is less than or equal to .1",integratorAccuracy);
        SimTK_ERRCHK1_ALWAYS(!(piecewiseRigidify) ,"[ParameterReader.cpp]","Making sure you set piecewiseRigidify to 0 .. it's :%d",piecewiseRigidify);

        MMBLOG_FILE_FUNC_LINE(ALWAYS, "RNABuilder will not reorder base pairs (prioritize 1) for multi-chain jobs.  If you have multiple chains make sure to set prioritize 0."<<endl);
        assert(!((prioritize==1) && (myBiopolymerClassContainer.getNumBiopolymers() >1)));
        assert(!((monteCarloRun==1) && (setChiBondMobility == 1)));
    } else MMBLOG_FILE_FUNC_LINE(WARNING, "You have set safeParameters = 0 and are now syntactically and semantically on your own.  May the RNA gods have mercy on you."<<endl);
    //numDutyCycles = dutyCycleArray.size();//d;


    MMBLOG_FILE_FUNC_LINE(DEBUG, "check 14.5"<<endl);
    MMBLOG_FILE_FUNC_LINE(DEBUG, "check 15"<<endl);
    if (lastStage == 0) {
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "For some reason lastStage is set to 0."<<endl);
    }//lastStage = calcHighestPriority();
    //MMBLOG_FILE_FUNC_LINE(DEBUG, "All base pairs, before removeBasePairsAcrossRigidStretches"<<endl);
    //if (setRemoveBasePairsAcrossRigidStretches) {removeBasePairsAcrossRigidStretches();}
    MMBLOG_FILE_FUNC_LINE(DEBUG, "All base pairs, before addHelicalStacking"<<endl);
    basePairContainer.printBasePairs();
    //if (setHelicalStacking){ 
    //    MMBLOG_FILE_FUNC_LINE(DEBUG, "check 17"<<endl);
    //    basePairContainer.addHelicalStacking(myBiopolymerClassContainer, _leontisWesthofClass, ntc_par_class,ntc_class_container);
    //}
    MMBLOG_FILE_FUNC_LINE(DEBUG, "All base pairs, before addHelicalStacking"<<endl);
    basePairContainer.printBasePairs();
    MMBLOG_FILE_FUNC_LINE(DEBUG, endl);
}

void ParameterReader::clearContainers(){
    contactInterfaceContainer.clear();
    clearConstraints();
    clearForces();
    clearBiopolymers();
    moleculeClassContainer.clear();
    // secondaryStructureStretchVector.clear() is now included in :
    myBiopolymerClassContainer.clear(); //secondaryStructureStretchVector.clear();
    MMBLOG_FILE_FUNC_LINE(INFO, ""<<endl);
    myBiopolymerClassContainer.writePdbStructureMapDiagnostics(); //This is just a debug flag
    displacementContainer.clear();
    atomSpringContainer.clear();
}

void ParameterReader::clearBiopolymers(){
    myBiopolymerClassContainer.clear();
}

void ParameterReader::clearForces(){
    atomSpringContainer.clear();
    myMonoAtomsContainer.clear();

    baseOperationVector.clear();
    basePairContainer.clear();
    basePairPartners.clear();  
    
    contactContainer.clear();
    
    densityContainer.clear();
    electroDensityContainer.clear();
    
    waterDropletAboutResidueVector.clear();
    waterDropletContainer.clear();
    physicsContainer.clear();
    //includeAllNonBondAtomsInResidueVector.clear(); // obsolete, replaced with physicsContainer
    includeAllResiduesWithinVector.clear();
    additionalCovalentBondVector.clear();
    includeIntraChainInterfaceVector.clear();
    
}

void ParameterReader::clearConstraints(){
    mobilizerContainer.clear();
    mobilizerContainer.interfaceContainer.clear();

    constraintToGroundContainer.clear();
}


void ParameterReader::initializeDefaults(const char * leontisWesthofInFileName){
    MMBLOG_FILE_FUNC_LINE(INFO, endl);
    ifstream inFile(leontisWesthofInFileName,ifstream::in);
    if (!inFile.good()) {
        inFile.close();
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "Failed to open leontisWesthofInFileName "<<leontisWesthofInFileName<<endl);
    }
    inFile.close();

    this->leontisWesthofInFileName = leontisWesthofInFileName;
    
    // SCF start rationalized initialization scheme
    clearContainers();

    //variables previously declared and initialized in Repel.h:
    ///useCIFFileFormat         = false;
    addAllAtomSterics        = false;
    addAllHeavyAtomSterics   = false;
    addBackboneOxygenForces  = false;
    addProteinBackboneSterics= false;
    addRNABackboneSterics= false;
    addSelectedAtoms         = false;
    addTestSpring = false;
    alignmentForcesIsGapped = true;
    alignmentForcesGapPenalty = -1;
    alignmentForcesDeadLengthFraction = 0;
    alignmentForcesDeadLength = 0;
    alignmentForcesDeadLengthIsFractionOfInitialLength = false;
    alignmentForcesForceConstant = 30.0;
    applyC1pSprings = true;
    biopolymerModificationVector.clear();
    calcBaseBodyFramesAtEveryTimeStep = true ; //scf maybe should set this to false?
    calcEnergy = true      ;
    totalEnergy = 0.0;
    potentialEnergy = 0.0;
    kineticEnergy = 0.0;
    checkSatisfied  = false; //scf should be true?
    //constrainRigidSegments = false;
    constraintTolerance = .00005;
    cutoffRadius  = .1; // in nm, of course
    cutoffAngle = 0    ;
    densityAtomFraction = 1.0;
    //densityFileName = "densityFileName-NOT-SET";
    //densityForceConstant = 1.;
    densityNoiseTemperature = 0.;
    densityNoiseScale    = 0.;
    densityNoiseComputeAutocorrelation = 0;
    densityReportAtEachAtomPosition = 0;
    densityFitPhosphates  = 1;
    //densityMapActivate = false;
    electroDensityFileName = "electroDensityFileName-NOT-SET";
    electroDensityForceConstant = 1.;
    excludedVolumeRadius = .06 ; //24; // in nanometers
    excludedVolumeStiffness = 10000000;
    //firstResidueMobilizerType = "Free";
    // firstStage = 0;
    fitDefaultTolerance =  0.3;
    globalAmberImproperTorsionScaleFactor = 0;
    globalBondBendScaleFactor =1 ;
    globalBondStretchScaleFactor = 1 ;
    globalBondTorsionScaleFactor =1;
    globalCoulombScaleFactor =0;
    globalGbsaScaleFactor =0    ;
    globalVdwScaleFactor=0;
    guessCoordinates = false;
    //halfSpaceMaxX=100;
    hardSphereStiffnessMultiplier= 1. ;
    inQVectorFileName = "NOT-SET";
    initialSeparation = 2; // nm
    integratorAccuracy =.0001;      
    integratorStepSize= .001;      
    integratorType= "RungeKuttaMerson";
    kbBackboneTorsionGlobalScaleFactor= 0.0;
    // lastStage=0;
    // leontisWesthofInFileName   = "./parameters.csv";
    loadTinkerParameterFile = false;
    outQVectorFileName  = "NOT-SET";
    magnesiumIonChainId= "X";
    magnesiumIonRadius  = 25;
    matchDefaultSkipTopLevelTransform = false;
    matchHydrogenAtomLocations = false;
    matchPurineN1AtomLocations = true ;
    matchingMinimizerTolerance = .1500;
    matchOptimize =  false;
    matchExact =true ;
    matchIdealized =false; // was false;
    minimize = false;
    monteCarloRun = false;
    monteCarloTemperature  = 400;
    monteCarloTemperatureIncrement = .1;

    nastGlobalBondTorsionScaleFactor = 10;
    noseHooverTime               = 10    ;
    numMagnesiumIons =0;
    numReportingIntervals =   90 ;
    outMonteCarloFileName="./MonteCarlo.pdb" ;
    outTrajectoryFileName = "NOT-SET";  //"/Users/samuelflores/svn/tar-dynamics/mymovie.pdb";
    //physicsWhereYouWantIt = false;
    physicsRadius = 0.0;
    piecewiseRigidify = false;
    planarityThreshold = .07;
    potentialType  = "HarmonicInverse";
    prioritize = false ;
    proteinCapping = false;
    useNACappingHydroxyls = true ;
    readInQVector = false;
    readPreviousFrameFile = true;
    readMagnesiumPositionsFromFile = true;
    removeRigidBodyMomentum = false; // caused us integration errors too many times ..
    removeDensityForcesFromRigidStretches = false; // Default is to allow density forces on rigid segments.
    removeMomentumPeriod =1;
    reportingInterval = 1.;// ps
    restrainingForceConstant  = 1000000;
    restrainingTorqueConstant= 10000;
    //resetBases;
    rigidifyFormedHelices = false;
    rigidifyTermini= false;
    satisfiedBasePairs =0;
    unSatisfiedBasePairs =0;
    scrubberPeriod  = reportingInterval * 40 ;
    safeParameters = true;
    //setChiBondAnti = false;
    setChiBondMobility = false;
    //setDefaultMDParameters = false;
    setDefaultStructurePredictionParameters =false;
    setDefaultThreadingParameters=false;

    setForceAndStericScrubber   =false ;
    setForceScrubber   =false ;
    setHelicalStacking =true;
    setInitialVelocities=false;
    setLoopBondMobility = false;//minimize;
    setOverallBondMobility = false;
    setRemoveBasePairsInRigidStretch = false; // This should be set to default to true, as soon as we can confirm this feature is working again.
    setRemoveBasePairsAcrossRigidStretches = false;
    setRepulsiveForce =false;
    setTemperature=true;
    smallGroupInertiaMultiplier = 1.0;
    waterInertiaMultiplier = 1.0;

    sphericalHelixCenter     = Vec3(0.,0.,0.);                                               
    sphericalHelixRadius     = 100.0;           
    sphericalHelixStartTheta = SimTK::Pi / 4.0; 
    sphericalHelixPhiOffset  = 0.0;             

    stackAllHelicalResidues = true ;
    thermostatType ="NoseHoover";
    tinkerParameterFileName = "NOT-SET";
    twoTransformForceMultiplier = 10;
    useFixedStepSize = false;
    useMultithreadedComputation = false;
    useOpenMMAcceleration = false;
    vanderWallSphereRadius = 300;
    velocityRescalingInterval= reportingInterval;
    verbose = false;
    vmdOutput = false;
    waterDropletMake = false;
    weldToGround = false;
    wkdpGlobalBondTorsionScaleFactor=0;
    writeCoordinates=true;
    writeDoublePrecisionTrajectories = false;
    writeFrameFile  =false;
    writeLastFrameFile = true;
    workingDirectory =  get_and_set_working_path();
    helixBondMobility  = BondMobility::Rigid ;
    loopBondMobility = BondMobility::Default;
    overallBondMobility = BondMobility::Free;
    chiBondMobility  = BondMobility::Free;
    qVector.clear();
    lastFrameFileName  = "NOT-SET"; //  "/Users/samuelflores/svn/tar-dynamics/last.pdb";
    previousFrameFileName= "NOT-SET";///Users/samuelflores/svn/tar-dynamics/last.pdb" ;
    enforceParallelness  = false;
    // end of variables improted from Repel.h

    sequence = "";
    proteinSequence = "";
    coarseNucleicAcidSequence = "";
    //numChains =0;
    numFirstResidues = 0;
    numResetBases    = 0;
    numProteinFirstResidues = 0;
    numProteinChains =0;
    numTemperatures =0;      
    numGlobalCoulombScaleFactors =0;
    numGlobalVdwScaleFactors =0;
    //numDutyCycles =0;    
    temperature = 10; 
    dutyCycle  = 1 ; 
    periodicallyUpdateParameters = false;
    // currentStage = 0;
    priority =0;
    //dutyCycleArray.clear();
    //dutyCyclePriority.clear();

    userVariables.clear();
    numRigidSegments.clear()   ;

    baseOperationVector.clear();
    basePairContainer.clear();
    #ifdef BuildNtC
    ntc_class_container.clear();
    #endif
    contactContainer.clear();
    mobilizerContainer.singleBondMobilityVector.clear();
    basePairPartnerVector.clear();
    densityContainer.clear();
    // SCF end `new rationalized initialization scheme

    detectConvergence = false;
    converged = false;
    convergenceTimeout = 5;
    convergenceEpsilon = 0.5;

    if (this->firstStage == this->currentStage) {
        _leontisWesthofClass.initialize(leontisWesthofInFileName);

        #ifdef BuildNtC
        ifstream ifile("parameters_user.csv");
        if(ifile){
             ntc_par_class.initialize("parameters_user.csv");
             ntc_par_class.initialize(leontisWesthofInFileName);
        } else {
             ntc_par_class.initialize(leontisWesthofInFileName);
        }
        #endif
    }
}


void ParameterReader::initialize(const char * parameterFileName ) {
    ifstream inFile(parameterFileName,ifstream::in);
    if (!inFile.good()) {
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "Failed to open user input file "<<parameterFileName<<endl);
    }

    stringstream u;
    //MMBLOG_FILE_FUNC_LINE(" About to initializeDefaults with parameterFileName = >"<<parameterFileName<<"< "<<endl; 
    initializeDefaults();
    initializeFromFileOnly(parameterFileName);

};


