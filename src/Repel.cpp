/* -------------------------------------------------------------------------- *
 *                           MMB (MacroMoleculeBuilder)                       *
 * -------------------------------------------------------------------------- *
 *                                                                            *
 * Copyright (c) 2011-12 by the Author.                                       *
 * Author: Samuel Flores                                                      *
 *                                                                            *
 * See RNABuilder.cpp for the copyright and usage agreement.                  *
 * -------------------------------------------------------------------------- */


#include "Repel.h"
#include "ProgressWriter.h"
#include "Utils.h"
#include <stdio.h>
#include <string.h>
#include   "MonoAtoms.h"
#include "BiopolymerClass.h"
#include "ConstraintContainer.h"
#include "ElectrostaticPotentialGridForce.h"
#include "PeriodicPdbAndEnergyWriter.h"
#include "CifOutput.h"

using namespace SimTK;
using namespace std;
/**
 * /brief This method sets the BondMobility for all bonds within each residue in a certain stretch of residues of any biopolymer   chain
 *
 * added by scf
 */
    void setBiopolymerResidueBondMobility (Biopolymer & myChain , BondMobility::Mobility  mobility, ResidueInfo::Index startResidue, ResidueInfo::Index endResidue){
        for (ResidueInfo::Index q =startResidue; q<=endResidue; q++) {
            myChain.setResidueBondMobility(q, mobility);
        }
    }

/**
 * /brief This function retrieves the MobilizedBody corresponding to a given chain ID, residue ID, and atom name. It figures out whether the chain ID belongs to a BiopolymerClass or MonoAtoms. It does not yet handle MoleculeClass.
 *
 * added by scf
 */
    MobilizedBody getMobilizedBody(SimbodyMatterSubsystem & myMatter, ParameterReader & parameterReader, BiopolymerClassContainer & biopolymerClassContainer, String myChainID, ResidueID myResidue, String myAtom) {
	    MobilizedBody mobilizedBody;
	if (biopolymerClassContainer.hasChainID(myChainID)) {
	    mobilizedBody = biopolymerClassContainer.updAtomMobilizedBody(myMatter,myChainID,myResidue,myAtom);
                //myMobilizedBody1 = myBiopolymerClassContainer.updAtomMobilizedBody(matter,myChain, myResidueID , myAtomName );
	} else if (parameterReader.moleculeClassContainer.hasChainID(myChainID)) {
	    // This is pending. I don't know yet what to do with this residue number..
	    MMBLOG_FILE_FUNC_LINE(CRITICAL, "You are trying to weld an atom of a custom molecule (of chain ID: >"<<myChainID <<"<) to another atom. Unfortunately this feature is not yet supported. You should, however, be able to weld a custom molecule to Ground."<<endl);
	} else if (parameterReader.myMonoAtomsContainer.hasChainID(myChainID)) {
	    mobilizedBody = parameterReader.myMonoAtomsContainer.getMonoAtoms(myChainID).updMobilizedBody(myResidue, myMatter); 
	}
	return mobilizedBody;
    };

/**
 * /brief This method sets the BondMobility for all bonds in a certain stretch of residues of a Protein        chain
 *
 * added by scf
 */
    void setProteinBondMobility (Biopolymer & myProteinChain , BondMobility::Mobility  mobility, ResidueInfo::Index startResidue, ResidueInfo::Index endResidue){
        for (ResidueInfo::Index q =startResidue; q<=endResidue; q++) {
            myProteinChain.setResidueBondMobility(q, mobility);
            if (q>startResidue) {
                std::stringstream ss1;
                ss1<<q-1<<"/C"  ;
                std::stringstream ss2;
                ss2<<q<<"/N";
                myProteinChain.setBondMobility(mobility ,ss1.str() ,ss2.str()  );
            }
        }
    }



double PointToPlaneDistance (Vec3 Point1, Vec3 Normal1, Vec3 Point2) {
   double d = (double)(-dot(Point1,Normal1));
   double D = (double)((d + dot(Point2,Normal1))/Normal1.norm());
   return D;
};

/**
 * 
 * 
 * /param 
 * myPdbResidueName1,2 must be one of "A","C","G","U".
 * bondingEdge1,2 must be one of "WatsonCrick","Hoogsteen","Sugar","Bifurcated".
 * glycosidicBondOrientation must be either "Cis" or "Trans".
 *
 */

/**
 * /brief This utility is used for the Concentricity constraint.  In this constraint, a certain atom on residue 1 is connected by a LinearSpring to a certain point OUTSIDE of an atom on residue 2.  location2 is that point, in the body frame of the atom on residue 2.
 * 
 */

    static bool compareUpper( const string& param, const char* symbol ) {

         std::string upperParam(param);
         std::string upperSym(symbol);

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





    void ConstrainedDynamics::constraintsAndRestraints  (ParameterReader & myParameterReader,BiopolymerClassContainer & myBiopolymerClassContainer,GeneralForceSubsystem & forces,SimbodyMatterSubsystem & matter, State & state, CompoundSystem & system) {
        MMBLOG_FILE_FUNC_LINE(INFO, endl);
        //#ifdef BUILD_MMB_SHARED_LIB
        //#pragma message ("BUILD_MMB_SHARED_LIB is currently defined")
        //#else
        //#pragma message ("BUILD_MMB_SHARED_LIB is NOT currently defined")
        //#endif
        MMBLOG_FILE_FUNC_LINE(INFO, endl);
 
        #if defined(USE_OPENMM) && defined(WARN_USE_OPENMM)
        #pragma message ("USE_OPENMM is defined")
        #elif !defined(USE_OPENMM) && defined(WARN_USE_OPENMM)
        #pragma message ("USE_OPENMM is NOT defined")
        #endif 
        MMBLOG_FILE_FUNC_LINE(INFO, endl);
        // here we turn the list of interface constraints from myParameterReader.constraintToGroundContainer.interfaceContainer into pairs of constraints, at most one pair for each pair of constrained chains.
        #ifdef USE_OPENMM    
        myParameterReader.constraintToGroundContainer.addSingleWeldConstraintPerInterfaceChainPair ( myBiopolymerClassContainer);
        #endif
        MMBLOG_FILE_FUNC_LINE(INFO, endl);
        myParameterReader.constraintToGroundContainer.applyConstrainChainRigidSegments ( myBiopolymerClassContainer,  system,  matter, state);
        
        MMBLOG_FILE_FUNC_LINE(INFO, "Now printing all constraints :"<<endl);
        myParameterReader.constraintToGroundContainer.printConstraintClasses();
        MMBLOG_FILE_FUNC_LINE(INFO, endl);
        vector<Constraint> myWeld2;        
        vector<Force> myLinearBushing;
       
        MobilizedBody myMobilizedBody1, myMobilizedBody2;
        Transform myTransform1, myTransform2;
        
        for (int i = 0; i < myParameterReader.constraintToGroundContainer.numConstraintClasses(); i++){
            ConstraintClass  myConstraintClass = myParameterReader.constraintToGroundContainer.getConstraintClass(i);
            MMBLOG_FILE_FUNC_LINE(INFO, "Applying constraint index : "<<i<<endl);
            myConstraintClass.print();
            if (myConstraintClass.getConstraintType() == WeldToGround) { // constrain an atom to ground
                String    myChain     = myConstraintClass.getChain1();
                if (myBiopolymerClassContainer.hasChainID(myChain)) { // We are constraining an atom on a biopolymer
			myTransform1 = Transform(Vec3(0));          
			myMobilizedBody1 =  matter.Ground();

			String myAtomName = myConstraintClass.getAtomName1();
			ResidueID myResidueID = myConstraintClass.getResidueID1();
			String    myChain     = myConstraintClass.getChain1();
			myMobilizedBody2 = myBiopolymerClassContainer.updAtomMobilizedBody(matter,myChain, myResidueID , myAtomName );
			myTransform2 = myMobilizedBody1.findBodyTransformInAnotherBody(state,myMobilizedBody2);            
			myWeld2.push_back( 
			   Constraint::Weld  (
			       myMobilizedBody1, // Ground body
			       myTransform1,     // Ground transform = Vec3(0)
			       myMobilizedBody2, // body of atom to be constrained to ground
                               myTransform2      // ground in the reference frame of the atom to be constrained
                           )
                        );
                } else if (myParameterReader.moleculeClassContainer.hasChainID(myChain)) {

                    
		    myTransform1 = Transform(Vec3(0));          
		    myMobilizedBody1 =  matter.Ground();
		    String myAtomName = myConstraintClass.getAtomName1();
		    ResidueID myResidueID = myConstraintClass.getResidueID1();
		    Compound::AtomIndex myAtomIndex = myParameterReader.moleculeClassContainer.updMoleculeClass(myChain).molecule.getAtomIndex(myAtomName);
		    MobilizedBodyIndex myAtomMobilizedBodyIndex = myParameterReader.moleculeClassContainer.updMoleculeClass(myChain).molecule.getAtomMobilizedBodyIndex(myAtomIndex);
		    myMobilizedBody2 = matter.updMobilizedBody(myAtomMobilizedBodyIndex); 
		    myTransform2 = myMobilizedBody1.findBodyTransformInAnotherBody(state,myMobilizedBody2);            
		    myWeld2.push_back( 
		       Constraint::Weld  (
		       myMobilizedBody1, // Ground body
		       myTransform1,     // Ground transform = Vec3(0)
		       myMobilizedBody2, // body of atom to be constrained to ground
		       myTransform2      // ground in the reference frame of the atom to be constrained
		       )
                    );
                } else {
                    MMBLOG_FILE_FUNC_LINE(CRITICAL, "You are trying to constriain chain "<<myChain<<" which is neither a biopolymer nor a custom molecule"<<endl);
                }
            } else if (myConstraintClass.getConstraintType() == WeldToAtom) { // constrain two atoms to each other
                String myAtomName = myConstraintClass.getAtomName1();
                ResidueID myResidueID = myConstraintClass.getResidueID1();
                String    myChain     = myConstraintClass.getChain1();
                // SCF : check whether we are dealing with biopolymer, monoAtom, or custom molecule, then generalize to handle all three. No further change should be necessary:
                
                //MobilizedBody getMobilizedBody(SimbodyMatterSubsystem & myMatter, ParameterReader & parameterReader, BiopolymerClassContainer & biopolymerClassContainer, String myChainID, ResidueID myResidue, String myAtom) {
                MMBLOG_FILE_FUNC_LINE(INFO, endl);
                myMobilizedBody1 = getMobilizedBody(matter,myParameterReader,myBiopolymerClassContainer, myChain, myResidueID, myAtomName); 
                // was:
                //myMobilizedBody1 = myBiopolymerClassContainer.updAtomMobilizedBody(matter,myChain, myResidueID , myAtomName );
                MobilizedBodyIndex myMobilizedBodyIndex1 = MobilizedBodyIndex(myMobilizedBody1);
                // was:
                //MobilizedBodyIndex myMobilizedBodyIndex1 = myBiopolymerClassContainer.updBiopolymerClass(myChain).getAtomMobilizedBodyIndex(matter, myResidueID , myAtomName );
                myTransform1 = Transform(Vec3(0));

                String myAtomName2 = myConstraintClass.getAtomName2();
                ResidueID myResidueID2 = myConstraintClass.getResidueID2();
                String    myChain2     = myConstraintClass.getChain2();
 		myMobilizedBody2 =  getMobilizedBody(matter,myParameterReader,myBiopolymerClassContainer, myChain2,myResidueID2,myAtomName2);
		// was:
                //myMobilizedBody2 = myBiopolymerClassContainer.updAtomMobilizedBody(matter,myChain2, myResidueID2 , myAtomName2 );
   		MobilizedBodyIndex myMobilizedBodyIndex2 = MobilizedBodyIndex(myMobilizedBody2);
		// was:
                //MobilizedBodyIndex myMobilizedBodyIndex2 = myBiopolymerClassContainer.updBiopolymerClass(myChain2).getAtomMobilizedBodyIndex(matter, myResidueID2 , myAtomName2 );
                if (myMobilizedBodyIndex1 == myMobilizedBodyIndex2) {
                    MMBLOG_FILE_FUNC_LINE(CRITICAL, "You are trying to weld a mobilized body to itself.  This could happen, for instance, if you are welding two atoms from the same rigid stretch to each other."<<endl);
                }
                myTransform2 = myMobilizedBody1.findBodyTransformInAnotherBody(state,myMobilizedBody2);
                myWeld2.push_back(
                Constraint::Weld  (
                   myMobilizedBody1, // first atom body
                   myTransform1,     // first atom transform = Vec3(0)
                   myMobilizedBody2, 
                   myTransform2  
                   )
                );

            } // of else if WeldToAtom

                     
            else if (myConstraintClass.getConstraintType() == CoupledCoordinate)  {
		Array_<MobilizedBodyIndex> myMobilizedBodyIndexArray_(2);
		myMobilizedBodyIndexArray_[0] = myBiopolymerClassContainer.updBiopolymerClass(myConstraintClass.getChain1()).getAtomMobilizedBodyIndex(matter,myConstraintClass.getResidueID1()    ,myConstraintClass.getAtomName1() );
		myMobilizedBodyIndexArray_[1] = myBiopolymerClassContainer.updBiopolymerClass(myConstraintClass.getChain2()).getAtomMobilizedBodyIndex(matter,myConstraintClass.getResidueID2()    ,myConstraintClass.getAtomName2() );
		Array_<MobilizerQIndex> myMobilizerQIndexArray_(2);
			myMobilizerQIndexArray_[0] = MobilizerQIndex(0);
			myMobilizerQIndexArray_[1] = MobilizerQIndex(0);
		    ///// adding SimTK::Constraint::CoordinateCoupler 
			Constraint::CoordinateCoupler (
			    matter,
				new ConstraintFunction(), // Defined in Util.h
			    //constraintEquation(coordMobod, coordQIndex, state),
			    myMobilizedBodyIndexArray_, //myConstraintClass.fetchMobilizedBodyIndexArray_ (myBiopolymerClassContainer, matter),
			    myMobilizerQIndexArray_ //myConstraintClass.fetchMobilizerQIndexArray_    (myBiopolymerClassContainer, matter, state)
                        );

		}  else {
			MMBLOG_FILE_FUNC_LINE(CRITICAL, "You have specified an unsupported ConstraintType : "<<myConstraintClass.getConstraintType() <<endl);
		    }


        } // of for

        for (int q=0;q<(int)myParameterReader.baseOperationVector.size();q++) {
            if (((myParameterReader.baseOperationVector[q]).BasePairIsTwoTransformForce).compare("constraint" ) == 0) {
                MMBLOG_FILE_FUNC_LINE(CRITICAL, "An error has occurred!  Constraints are no longer specified with this data structure. "<<endl);
            }
  
            if ((   //(((myParameterReader.baseOperationVector[q]).BasePairIsTwoTransformForce).compare("constraint" ) == 0)||  // this no longer how we do constraints
                (((myParameterReader.baseOperationVector[q]).BasePairIsTwoTransformForce).compare("restraint" ) == 0) ||
                (((myParameterReader.baseOperationVector[q]).BasePairIsTwoTransformForce).compare("restrainToGround" ) == 0) //||
                //(((myParameterReader.baseOperationVector[q]).BasePairIsTwoTransformForce).compare("constrainToGround" ) == 0) 
                ) 
                &&   (((myParameterReader.baseOperationVector[q]).FirstBPEdge).compare("Weld" ) == 0))
           {
                MobilizedBody myMobilizedBody1, myMobilizedBody2;
                Transform myTransform1, myTransform2;
                string myAtom1Name, myAtom2Name;

                if (compareUpper( myParameterReader.baseOperationVector[q].FirstBPChain,"GROUND")) {
                    myTransform1 = Transform(Vec3(0));          
                    myMobilizedBody1 =  matter.Ground();
                }
                else if (myBiopolymerClassContainer.hasChainID(myParameterReader.baseOperationVector[q].FirstBPChain)) {
                    myAtom1Name = myBiopolymerClassContainer.updBiopolymerClass(myParameterReader.baseOperationVector[q].FirstBPChain).getRepresentativeAtomName();
                    myMobilizedBody1 = myBiopolymerClassContainer.updAtomMobilizedBody(matter,myParameterReader.baseOperationVector[q].FirstBPChain, myParameterReader.baseOperationVector[q].FirstBPResidue, myAtom1Name );
                    myTransform1 = Transform(Vec3(0));          
                } else if (myParameterReader.myMonoAtomsContainer.hasChainID(myParameterReader.baseOperationVector[q].FirstBPChain ))         {
                    myMobilizedBody1 = myParameterReader.myMonoAtomsContainer.getMonoAtoms(myParameterReader.baseOperationVector[q].FirstBPChain).updMobilizedBody(myParameterReader.baseOperationVector[q].FirstBPResidue,
                        matter
                        );
                    myTransform1 = Transform(Vec3(0));          
                } else 
                {MMBLOG_FILE_FUNC_LINE(CRITICAL, "The first chain in the constraint or restraint was neither a biopolymer nor a monoAtoms object, and you are not constraining to ground."<<endl); }

                if (myBiopolymerClassContainer.hasChainID(myParameterReader.baseOperationVector[q].SecondBPChain)) {
                    myAtom2Name = myBiopolymerClassContainer.updBiopolymerClass(myParameterReader.baseOperationVector[q].SecondBPChain).getRepresentativeAtomName();
                    myMobilizedBody2 = myBiopolymerClassContainer.updAtomMobilizedBody(matter,myParameterReader.baseOperationVector[q].SecondBPChain, myParameterReader.baseOperationVector[q].SecondBPResidue, myAtom2Name );
                    myTransform2 = myMobilizedBody1.findBodyTransformInAnotherBody(state,myMobilizedBody2);            
                } else if (myParameterReader.myMonoAtomsContainer.hasChainID((myParameterReader.baseOperationVector[q]).SecondBPChain ))         {
                    myMobilizedBody2 = myParameterReader.myMonoAtomsContainer. getMonoAtoms(myParameterReader.baseOperationVector[q].SecondBPChain).updMobilizedBody(myParameterReader.baseOperationVector[q].SecondBPResidue,
                        matter
                        );
                    myTransform2 = myMobilizedBody2.findBodyTransformInAnotherBody(state,myMobilizedBody1);            
                } else 
                {MMBLOG_FILE_FUNC_LINE(CRITICAL, "Error! The second chain in the constraint or restraint was neither a biopolymer nor a monoAtoms object."<<endl); }

               MMBLOG_FILE_FUNC_LINE(INFO, "system.realizeTopology() "<<endl);
               state = system.realizeTopology();                                        
               system.realize(state,Stage::Position);

               MMBLOG_FILE_FUNC_LINE(INFO, "Enforcing constraint or restraint from : "<<string(myParameterReader.baseOperationVector[q].FirstBPChain) <<", " << myParameterReader.baseOperationVector[q].FirstBPResidue.outString() << " to : " << string( myParameterReader.baseOperationVector[q].SecondBPChain)<<", "<<myParameterReader.baseOperationVector[q].SecondBPResidue.outString()<<  std::endl);
               if ((myParameterReader.baseOperationVector[q].BasePairIsTwoTransformForce.compare("restrainToGround" ) == 0) ||
                    (myParameterReader.baseOperationVector[q].BasePairIsTwoTransformForce.compare("restraint" ) == 0) )   {
                    myLinearBushing.push_back (
                        Force::LinearBushing (
                           forces,                                            
                           myMobilizedBody1,
                           myTransform1,
                           myMobilizedBody2,
                           myTransform2,
                           Vec6(
                               myParameterReader.restrainingTorqueConstant,myParameterReader.restrainingTorqueConstant,myParameterReader.restrainingTorqueConstant,
                               myParameterReader.restrainingForceConstant,myParameterReader.restrainingForceConstant,myParameterReader.restrainingForceConstant

                               ),
                           Vec6(.0)
                           )
                           );
                }
                else if ((myParameterReader.baseOperationVector[q].BasePairIsTwoTransformForce.compare("constrainToGround" ) == 0) ||
                    (myParameterReader.baseOperationVector[q].BasePairIsTwoTransformForce.compare("constraint" ) == 0) ) {

                   myWeld2.push_back( 
                       Constraint::Weld  (
                           myMobilizedBody1,
                           myTransform1,
                           myMobilizedBody2,
                           myTransform2
                           )
                       );
                }
                else           
                {MMBLOG_FILE_FUNC_LINE(CRITICAL, "Error! You have specified an unsupported constraint or restraint."<<endl);}

                MMBLOG_FILE_FUNC_LINE(INFO, "system.realizeTopology() "<<endl);
                state = system.realizeTopology();
                system.realize(state,Stage::Position);
           } // of if constraint,restraint, etc.
        } // of for q

      };  // of class

ConstrainedDynamics::ConstrainedDynamics(ParameterReader * myParameterReader) : _parameterReader(myParameterReader),
                                                                                _matter(    _system),
                                                                                _forces(    _system),
                                                                                _contacts(  _system),
                                                                                _dumm(_system) {
    _study = NULL;
    _ts = NULL;
    _nextFrame = 1;
    _previousTimeStep = 0.0;
}

ConstrainedDynamics::~ConstrainedDynamics()
{
    if(_study)
        delete _study;

    if(_ts)
        delete _ts;
}

void ConstrainedDynamics::initializeDumm(){
    _parameterReader->configureDumm(_dumm);
}

int  ConstrainedDynamics::initializeBiopolymersAndCustomMolecules()
{
    return initializeBiopolymersAndCustomMolecules(_system);
}

int  ConstrainedDynamics::initializeBiopolymersAndCustomMolecules(CompoundSystem & system ){
    _parameterReader->moleculeClassContainer.initializeCompounds (_dumm);
    MMBLOG_FILE_FUNC_LINE(INFO, endl);
    _parameterReader->moleculeClassContainer.matchDefaultConfiguration (_parameterReader->readPreviousFrameFile, _parameterReader->previousFrameFileName, _parameterReader->matchExact,  _parameterReader->matchIdealized  );
    MMBLOG_FILE_FUNC_LINE(INFO, endl);
    _parameterReader->moleculeClassContainer.adoptCompounds (system);
    MMBLOG_FILE_FUNC_LINE(INFO, endl);

    //MMBLOG_FILE_FUNC_LINE("This is temporary .. SCF"<<endl;
    int returnValue = 0;
    returnValue = _parameterReader->myBiopolymerClassContainer.initializeBiopolymers(system,
                                                                       _parameterReader->proteinCapping, 
                                                                       _parameterReader->matchExact,  
                                                                       _parameterReader->matchIdealized, 
                                                                       _parameterReader->matchOptimize , 
                                                                       _parameterReader->matchHydrogenAtomLocations, 
                                                                       _parameterReader->matchPurineN1AtomLocations, 
                                                                       _parameterReader->guessCoordinates, 
                                                                       _parameterReader->initialSeparation, 
                                                                       _parameterReader->displacementContainer.getDisplacementVector(), 
                                                                       _parameterReader->matchingMinimizerTolerance, 
                                                                       _parameterReader->planarityThreshold );
    MMBLOG_FILE_FUNC_LINE(INFO, "returnValue = "<<returnValue<<endl);
    if (returnValue)
    {
        MMBLOG_FILE_FUNC_LINE(WARNING, "There was an error in initializeBiopolymers"<<endl);
        //return 1; 
    };
     
    #ifdef USE_OPENMM 
    _parameterReader->myBiopolymerClassContainer.initializeAtomInfoVectors(_matter); // May not work to use the DuMM version of this, due to possible topology not being realized.
    #endif
    return returnValue;
}

void ConstrainedDynamics::initializeBiopolymer(String chainID, String inputPDBFile){
    bool readFromFile = false;
    if(inputPDBFile != "")
        readFromFile = true;
    _parameterReader->myBiopolymerClassContainer.initializeBiopolymer(chainID, 
                                                                      _system,
                                                                      _parameterReader->proteinCapping, 
                                                                      _parameterReader->matchExact,  
                                                                      _parameterReader->matchIdealized, 
                                                                      _parameterReader->matchOptimize , 
                                                                      _parameterReader->matchHydrogenAtomLocations, 
                                                                      _parameterReader->matchPurineN1AtomLocations,
                                                                      _parameterReader->guessCoordinates,
                                                                      _parameterReader->initialSeparation, 
                                                                      _parameterReader->displacementContainer.getDisplacementVector(), 
                                                                      _parameterReader->matchingMinimizerTolerance, 
                                                                      _parameterReader->planarityThreshold , 
                                                                      _parameterReader->myBiopolymerClassContainer.secondaryStructureStretchVector);
}

void ConstrainedDynamics::initializeMoleculesAndBonds()
{
    initializeMoleculesAndBonds(_system, _dumm, _matter);
}

void ConstrainedDynamics::initializeMoleculesAndBonds(CompoundSystem & system, DuMMForceFieldSubsystem & dumm, SimbodyMatterSubsystem & matter){
    _parameterReader->myMonoAtomsContainer.initialize(system ,  _parameterReader->readPreviousFrameFile, _parameterReader->previousFrameFileName, _parameterReader->matchExact,  _parameterReader->matchIdealized  );
    // This is actually broken now. Fix in void BiopolymerClassContainer::waterDropletAboutResidues 
    if (_parameterReader->waterDropletAboutResidueVector.size() > 0) _parameterReader->myBiopolymerClassContainer.waterDropletAboutResidues( _parameterReader->waterDropletAboutResidueVector,    _parameterReader->waterDropletContainer) ;

    //MMBLOG_FILE_FUNC_LINE(endl;
    _parameterReader->waterDropletContainer.addWaterMolecules(system, dumm, matter, _parameterReader->myBiopolymerClassContainer);
    //MMBLOG_FILE_FUNC_LINE(endl;
    _parameterReader->waterDropletContainer.validateWaterVectors();
    //MMBLOG_FILE_FUNC_LINE(endl;
    _parameterReader->waterDropletContainer.matchDefaultConfiguration (_parameterReader->readPreviousFrameFile, _parameterReader->previousFrameFileName, _parameterReader->matchExact,  _parameterReader->matchIdealized  );
    //MMBLOG_FILE_FUNC_LINE(endl;
    _parameterReader->waterDropletContainer.adopt (system,_parameterReader->readPreviousFrameFile);
    //MMBLOG_FILE_FUNC_LINE(endl;
    // Moved this from ParameterReader.cpp:
    for (size_t i = 0; i < _parameterReader->additionalCovalentBondVector.size(); i++) {
        CovalentBondClass myBond = _parameterReader->additionalCovalentBondVector[i];
        if (myBond.getChain1().compare(myBond.getChain2()) != 0 ){
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "Chain IDs of both atoms must be the same."<<endl);
        }
	if (_parameterReader->myBiopolymerClassContainer.hasChainID(myBond.getChain1())) {
	    _parameterReader->myBiopolymerClassContainer.updBiopolymerClass(myBond.getChain1() ).addRingClosingBond(myBond);
	} else if (_parameterReader->moleculeClassContainer.hasChainID(myBond.getChain1())) {
            _parameterReader->moleculeClassContainer.updMoleculeClass(myBond.getChain1()).addRingClosingBond(myBond);
        } else {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "Unexplained error!"<<endl);
        }

    }
    MMBLOG_FILE_FUNC_LINE(INFO, "If you want to close the cysteine bridges, you can issue:"<< endl);
    //_parameterReader->myBiopolymerClassContainer.createDisulphideBridges(std::cout); // Should make this user controllable. // there is an issue with immediately mutating.. consider adding to mutation vector..
}

#ifdef USE_OPENMM
void ConstrainedDynamics::setInterfaceMobilizers(){
    MMBLOG_FILE_FUNC_LINE(INFO, endl);
    setInterfaceMobilizers(_system, _matter, _state);
}

void ConstrainedDynamics::setInterfaceMobilizers(CompoundSystem & system, SimbodyMatterSubsystem & matter, State & state){
    MMBLOG_FILE_FUNC_LINE(INFO, "Checking interface mobilizers .."<<endl);
    _parameterReader->myBiopolymerClassContainer.validateAtomInfoVectors();
    _parameterReader->mobilizerContainer.addStretchesToVectorFromInterfaceContainer(_parameterReader->myBiopolymerClassContainer);
    MMBLOG_FILE_FUNC_LINE(INFO, "Checking domains interface mobilizers .."<<endl);
    _parameterReader->mobilizerContainer.addMobilizerDomainsInterfacesToVector(_parameterReader->mobilizerDomainsInterfaceVector, _parameterReader->myBiopolymerClassContainer);
    MMBLOG_FILE_FUNC_LINE(INFO, "Domains interface mobilizers checked"<<endl);
    _parameterReader->mobilizerContainer.printMobilizerStretches();
    MMBLOG_FILE_FUNC_LINE(INFO, "system.getNumRealizeCalls () = "<< system.getNumRealizeCalls () <<endl);
    MMBLOG_FILE_FUNC_LINE(INFO, "system.realizeTopology() "<<endl);
    state = system.realizeTopology();
    MMBLOG_FILE_FUNC_LINE(INFO, endl);
    system.realize(state,Stage::Position);
}
#endif
void ConstrainedDynamics::setMobilizers()
{
    MMBLOG_FILE_FUNC_LINE(DEBUG, endl);
    #ifdef USE_OPENMM
    _parameterReader->mobilizerContainer.createMobilizersWithin(_parameterReader->myBiopolymerClassContainer,_state);
    #endif
    MMBLOG_FILE_FUNC_LINE(INFO, "Done adding MobilizerStretch's .  However we haven't yet dealt with 'Default' MobilizerStretch's. "<<endl);
    MMBLOG_FILE_FUNC_LINE(INFO, "Printing all mobilizer stretches: "<<endl);
    _parameterReader->mobilizerContainer.printMobilizerStretches();

    MMBLOG_FILE_FUNC_LINE(INFO, "About to clear MobilizerStretch's that have been overridden, or specified to remain at BondMobility 'Default' "<<endl);

    // removing all MobilizerStretch's that modify residues specified to be at BondMobility = "Default"
    for (int i = 0; i < _parameterReader->mobilizerContainer.getNumResidueStretches(); i++) {
        MMBLOG_FILE_FUNC_LINE(INFO, "Examining mobilizer stretch "<<i<<" to see if it has bondMobility Default :"<<endl);

        const auto & residueStretchVector = _parameterReader->mobilizerContainer.getResidueStretchVector();
        residueStretchVector[i].printStretch();
        if (residueStretchVector[i].getBondMobilityString().compare("Default") == 0){
            ResidueStretch tempResidueStretch = _parameterReader->mobilizerContainer.getResidueStretch(i);
            MMBLOG_FILE_FUNC_LINE(INFO, "Detected that BondMobility : >"<<residueStretchVector[i].getBondMobilityString()<<"< is Default."<<endl);
            _parameterReader->myBiopolymerClassContainer.updBiopolymerClass(tempResidueStretch.getChain()).selectivelyRemoveResidueStretchFromContainer(tempResidueStretch,_parameterReader->mobilizerContainer);
            i--; // Compensating for the fact that the Default stretch i has been deleted.
            
        } else {
            MMBLOG_FILE_FUNC_LINE(WARNING, "Detected that BondMobility : >"<<residueStretchVector[i].getBondMobilityString()<<"< is NOT Default."<<endl);
        }    
    }  
    // Done dealing with "Default" mobilizer stretches
    MMBLOG_FILE_FUNC_LINE(INFO, "Done dealing with MobilizerStretch's specified to remain at BondMobility 'Default' "<<endl);
    MMBLOG_FILE_FUNC_LINE(INFO, "Printing all mobilizer stretches: "<<endl);
    _parameterReader->mobilizerContainer.printMobilizerStretches();
    
    _parameterReader->mobilizerContainer.setBiopolymerBondMobility(_parameterReader->myBiopolymerClassContainer);
    _parameterReader->myBiopolymerClassContainer.setSingleBondMobility(_parameterReader->mobilizerContainer.singleBondMobilityVector);
}

void ConstrainedDynamics::createMultibodyTree()
{
    MMBLOG_FILE_FUNC_LINE(INFO, endl);
    createMultibodyTree(_system, _state);
}

void ConstrainedDynamics::createMultibodyTree(CompoundSystem & system, State & state){
    MMBLOG_FILE_FUNC_LINE(INFO, "Calling modelOneCompound "<< system.getNumCompounds()<<" times."<<endl);
    for (CompoundSystem::CompoundIndex c(0); c < system.getNumCompounds(); ++c) {
        String currentChain = String(system.getCompound(c).getPdbChainId());
        if (_parameterReader->myBiopolymerClassContainer.hasChainID(currentChain)){
            String tempMobilizerType =  _parameterReader->myBiopolymerClassContainer.updBiopolymerClass(currentChain).getFirstResidueMobilizerType();
            MMBLOG_FILE_FUNC_LINE(INFO, "Attaching biopolymer chain "<<currentChain<<" using mobilizer type = "<<tempMobilizerType<<endl);
            system.modelOneCompound(c,tempMobilizerType  );
        } else {
            String tempMobilizerType ("Free");
            MMBLOG_FILE_FUNC_LINE(INFO, "Attaching non-biopolymer chain "<<currentChain<<" using mobilizer type = "<<tempMobilizerType<<endl);
            system.modelOneCompound(c, tempMobilizerType );  // for non-biopolymers, use the default behavior, which is for root atom to be connected with a Free mobilizer.
        }
    }
    MMBLOG_FILE_FUNC_LINE(INFO, "system.realizeTopology() "<<endl);
    state = system.realizeTopology();
    system.realize(state,Stage::Position);
}

void ConstrainedDynamics::initializeCustomForcesConstraints(){
    _output.open(_parameterReader->outTrajectoryFileName.c_str(),  ios_base::app);

    _parameterReader->waterDropletContainer.addTethers(_parameterReader->atomSpringContainer);

    _parameterReader->waterDropletContainer.validateWaterVectors();
    MMBLOG_FILE_FUNC_LINE(INFO, endl);
    _parameterReader->myBiopolymerClassContainer.multiplySmallGroupInertia( _parameterReader->smallGroupInertiaMultiplier, _system,_matter, _state      );
    _parameterReader->waterDropletContainer.multiplySmallGroupInertia( _parameterReader->waterInertiaMultiplier, _system,_matter, _state      );
    //MMBLOG_FILE_FUNC_LINE(endl;
    #ifdef USE_OPENMM  
    _parameterReader->myBiopolymerClassContainer.initializeAtomInfoVectors(_matter,_dumm); // investigate moving this outside the conditional, so it's available for DensityForce below.
    #endif
    MMBLOG_FILE_FUNC_LINE(INFO, "You have specified _parameterReader->physicsRadius = "<<_parameterReader->physicsRadius<<endl);
    //_parameterReader->myBiopolymerClassContainer.printAtomInfoVector(); // Looks fine at this point ..
    if ( _parameterReader->physicsRadius > 0.0000001) 
    {
        MMBLOG_FILE_FUNC_LINE(INFO, endl);
        _parameterReader->myBiopolymerClassContainer.physicsZone (_parameterReader->includeAllResiduesWithinVector, _parameterReader->physicsRadius,_matter,_state); 
    }
    MMBLOG_FILE_FUNC_LINE(INFO, "_parameterReader->includeAllResiduesWithinVector.size() = "<<_parameterReader->includeAllResiduesWithinVector.size() << endl);
    #ifdef USE_OPENMM
    // Act on physicsInterfaces command. Turn the interfaces into specific lists of residues:
    _parameterReader->  physicsContainer.addStretchesToVectorFromInterfaceContainer(_parameterReader->myBiopolymerClassContainer);
    #endif
    MMBLOG_FILE_FUNC_LINE(INFO, "_parameterReader->includeAllResiduesWithinVector.size() = "<<_parameterReader->includeAllResiduesWithinVector.size() << endl);
    // SCF this is a good place to insert includeInterChainInterface processing. 
    if (_parameterReader->includeIntraChainInterfaceVector.size() >0) 
    {
        //_parameterReader->myBiopolymerClassContainer.initializeAtomInfoVectors(_matter,_dumm); // investigate moving this outside the conditional, so it's available for DensityForce below.
       for (size_t i = 0; i < _parameterReader->includeIntraChainInterfaceVector.size() ; i++) 
       {
           MMBLOG_FILE_FUNC_LINE(INFO, endl);
           _parameterReader->myBiopolymerClassContainer.addIntraChainInterfaceResidues(_parameterReader->includeIntraChainInterfaceVector[i].Chain,  
               _parameterReader->physicsContainer.updResidueStretchVector() //includeAllNonBondAtomsInResidueVector
               , _parameterReader->includeIntraChainInterfaceVector[i].Depth ,  
	       _matter,_state) ;
           MMBLOG_FILE_FUNC_LINE(INFO, endl);
       }
    }
    if (_parameterReader->includeIntraChainInterfaceVector.size() >0) {
        //_parameterReader->myBiopolymerClassContainer.initializeAtomInfoVectors(_matter,_dumm); // investigate moving this outside the conditional, so it's available for DensityForce below.
    for (size_t i = 0; i < _parameterReader->includeIntraChainInterfaceVector.size() ; i++) {
        _parameterReader->myBiopolymerClassContainer.addIntraChainInterfaceResidues(_parameterReader->includeIntraChainInterfaceVector[i].Chain,  
	    _parameterReader->physicsContainer.updResidueStretchVector() , //includeAllNonBondAtomsInResidueVector ,
	    _parameterReader->includeIntraChainInterfaceVector[i].Depth ,  
            _matter,_state) ;
    }
    }
    constraintsAndRestraints(*_parameterReader, _parameterReader->myBiopolymerClassContainer, _forces, _matter, _state,_system);
    _parameterReader->myBiopolymerClassContainer.computeCorrection(_parameterReader->_leontisWesthofClass, _parameterReader->basePairContainer.myBasePairVector, _state, _matter);
    _parameterReader->setLeontisWesthofBondRowIndex(); 
    
    _parameterReader->addC1pSprings(_parameterReader->_leontisWesthofClass);         
    _parameterReader->waterDropletContainer.validateWaterVectors();

    _parameterReader->applyAtomSprings(_matter,_forces, _state);
    #ifdef USE_OPENMM
    _parameterReader->contactContainer.createContactsWithin(_parameterReader->myBiopolymerClassContainer,_state);
    #endif
    _parameterReader->contactContainer.printContacts();
    _parameterReader->contactContainer.applyContactsToBiopolymers (_parameterReader->myBiopolymerClassContainer,   _contacts,  _forces,_matter, _parameterReader->_leontisWesthofClass, _parameterReader->excludedVolumeRadius, _parameterReader->excludedVolumeStiffness);
    AllTwoTransformLinearSprings * myAllTwoTransformLinearSpringsPointer = 
        new AllTwoTransformLinearSprings( _matter,  *_parameterReader,  _parameterReader->_leontisWesthofClass, _parameterReader->myBiopolymerClassContainer, _output);
    Force::Custom(_forces, myAllTwoTransformLinearSpringsPointer);
    #ifdef BuildNtC    
    NTC_Torque * myNTC_Torque = new NTC_Torque( _matter,  *_parameterReader,  _parameterReader->ntc_par_class, _parameterReader->myBiopolymerClassContainer, _output);
    Force::Custom(_forces, myNTC_Torque);
    #endif
    MMBLOG_FILE_FUNC_LINE(INFO, endl);
    if (_parameterReader->densityContainer.numDensityStretches() > 0) 
    {
        MMBLOG_FILE_FUNC_LINE(INFO, endl);
        _parameterReader->myDensityMap.setNoiseTemperature(_parameterReader->densityNoiseTemperature);
        MMBLOG_FILE_FUNC_LINE(INFO, "setting noiseScale to :" <<_parameterReader->densityNoiseScale<<endl);
        _parameterReader->myDensityMap.setNoiseScale(_parameterReader->densityNoiseScale);
        MMBLOG_FILE_FUNC_LINE(INFO, "set noiseScale to :" <<_parameterReader->myDensityMap.getNoiseScale() <<endl);
        _parameterReader->myDensityMap.loadParametersAndDensity(_parameterReader->densityFileName);
        MMBLOG_FILE_FUNC_LINE(INFO, endl);
        if (_parameterReader->myDensityMap.getNoiseScale() > 0.0000001) {
            MMBLOG_FILE_FUNC_LINE(INFO, endl);
            _parameterReader->myDensityMap.populateNoiseMap();
            _parameterReader->myDensityMap.writeDensityMapXplor("./noise.xplor", 0, 1);   // write noise only
            _parameterReader->myDensityMap.writeDensityMapXplor("./density.xplor", 1,0);  // write noiseless density only
            _parameterReader->myDensityMap.writeDensityMapXplor("./noisyMap.xplor", 1,1); // write the sum of density + noise
            if (_parameterReader->densityNoiseComputeAutocorrelation){
                MMBLOG_FILE_FUNC_LINE(INFO, "First, we write out the noise autocorrelation: "<<endl);
                _parameterReader->myDensityMap.densityAutocorrelation(1,0); // Arguments are : calculate noise correlation = 1, calculate density corrleation = 0
                MMBLOG_FILE_FUNC_LINE(INFO, "Second, we write out the autocorrelation of the input density map: "<<endl);
                _parameterReader->myDensityMap.densityAutocorrelation(0,1);
                MMBLOG_FILE_FUNC_LINE(INFO, endl);
            }
        }
        MMBLOG_FILE_FUNC_LINE(INFO, endl);
        _parameterReader->myDensityMap.precomputeGradient();
        MMBLOG_FILE_FUNC_LINE(INFO, endl);
        _parameterReader->myDensityMap.precomputeGradientDerivatives();
        MMBLOG_FILE_FUNC_LINE(INFO, endl);
        MMBLOG_FILE_FUNC_LINE(INFO, endl);
        // This needs to be re-done here, with the "dumm" version of initializeAtomInfoVectors, which sets atomic numbers, masses, etc. Consider making a cheaper version of initializeAtomInfoVectors that only updates the existing MMBAtomInfo's, though I'm not sure this would save much.
        //_parameterReader->myBiopolymerClassContainer.initializeAtomInfoVectors(_matter,_dumm); // SCF I'm not sure this is ready to be moved up .. trying dumm version earlier causes a crash
        _parameterReader->myBiopolymerClassContainer.validateAtomInfoVectors();
        DensityForce * myDensityForce = new DensityForce( _matter, 
                                                          *_parameterReader,  
                                                          _parameterReader->myDensityMap, 
                                                          _dumm, 
                                                          _parameterReader->myBiopolymerClassContainer, 
                                                          _output
                                                         );
        Force::Custom(_forces, myDensityForce);
    }
    if (_parameterReader->electroDensityContainer.numDensityStretches() > 0) 
    {
        _parameterReader->myElectroDensityMap.loadParametersAndDensity(_parameterReader->electroDensityFileName);
        _parameterReader->myElectroDensityMap.precomputeGradient();
        _parameterReader->myElectroDensityMap.precomputeGradientDerivatives();
        // This needs to be re-done here, with the "dumm" version of initializeAtomInfoVectors, which sets atomic numbers, masses, etc. Consider making a cheaper version of initializeAtomInfoVectors that only updates the existing MMBAtomInfo's, though I'm not sure this would save much.
        //_parameterReader->myBiopolymerClassContainer.initializeAtomInfoVectors(_matter,_dumm); 
        _parameterReader->myBiopolymerClassContainer.validateAtomInfoVectors();
        ElectrostaticPotentialGridForce * myElectroForce = new ElectrostaticPotentialGridForce( _matter, 
                                                          *_parameterReader,  
                                                          _parameterReader->myElectroDensityMap, 
                                                          _dumm, 
                                                          _parameterReader->myBiopolymerClassContainer, 
                                                          _output
                                                         );
        Force::Custom(_forces, myElectroForce);
    }
    MMBLOG_FILE_FUNC_LINE(INFO, endl);
    //_parameterReader->myBiopolymerClassContainer.printAtomInfoVector(); //  at this point.. 
    MMBLOG_FILE_FUNC_LINE(INFO, "_parameterReader->includeAllResiduesWithinVector.size() = "<<_parameterReader->includeAllResiduesWithinVector.size() << endl);
    MMBLOG_FILE_FUNC_LINE(INFO, "system.realizeTopology() "<<endl);
    _state = _system.realizeTopology();
    _system.realize(_state,Stage::Position);

    MMBLOG_FILE_FUNC_LINE(INFO, "_parameterReader->includeAllResiduesWithinVector.size() = "<<_parameterReader->includeAllResiduesWithinVector.size() << endl);


    MMBLOG_FILE_FUNC_LINE(INFO, "Did you intend to use PhysicsWhereYouWantIt? Checking size of _parameterReader->physicsContainer.getNumResidueStretches() " //includeAllNonBondAtomsInResidueVector.size() = "
        <<_parameterReader->physicsContainer.getNumResidueStretches()<<endl);
    MMBLOG_FILE_FUNC_LINE(INFO, "About to work on physicsWhereYouWantIt. DuMM now has included : "<< _dumm.getNumIncludedAtoms () <<" atoms. "<<endl);


    MMBLOG_FILE_FUNC_LINE(INFO, endl);
    //_parameterReader->myBiopolymerClassContainer.printAtomInfoVector(); //  at this point.. 


    bool myPhysicsWhereYouWantItActive =  
       ((_parameterReader->includeNonBondAtomInBiopolymerVector.size() > 0) ||
    (_parameterReader->physicsContainer.getNumResidueStretches() /*includeAllNonBondAtomsInResidueVector.size()*/ > 0) ||
        (_parameterReader->includeAllResiduesWithinVector.size() > 0)) ;
 
    if (myPhysicsWhereYouWantItActive) {
        MMBLOG_FILE_FUNC_LINE(INFO, "physicsWhereYouWantIt now has included : "<< _dumm.getNumIncludedAtoms () <<" atoms. Clearing list .."<<endl);
        _dumm.clearIncludedNonbondAtomList();
        MMBLOG_FILE_FUNC_LINE(INFO, "system.realizeTopology() "<<endl);
        _state = _system.realizeTopology();
        _system.realize(_state,Stage::Position);
        MMBLOG_FILE_FUNC_LINE(INFO, "Cleared non-bonded atom list. Realized _state to the position stage. physicsWhereYouWantIt now has included : "<< _dumm.getNumIncludedAtoms () <<" atoms. Note that this number might not be zero -- we do not clear the bonded atom list."<<endl);
        // Could end this loop here. The only point is to clear the included atom list, in cases that any physicsWhereYouWantIt is active.
    }  
    else {
        // Do nothing.  Don't clear included atom list. Non-bonded forces are off, so the simulation should be cheap.
    // Note that bonded forces may still be more than needed depending on mobilizers.  Worth checking to see if futher optimization is needed.
    }
   
    MMBLOG_FILE_FUNC_LINE(INFO, endl);
    //_parameterReader->myBiopolymerClassContainer.printAtomInfoVector(); //  at this point.. 


    if (_parameterReader->includeNonBondAtomInBiopolymerVector.size() > 0) 
    {
        MMBLOG_FILE_FUNC_LINE(INFO, "About to read a list of included atoms in biopolymer chains. physicsWhereYouWantIt now has included : "<< _dumm.getNumIncludedAtoms () <<" atoms. "<<endl);
        _parameterReader->myBiopolymerClassContainer.includeNonBondAtoms(  _parameterReader->includeNonBondAtomInBiopolymerVector,  _state,  _dumm);
            MMBLOG_FILE_FUNC_LINE(INFO, "system.realizeTopology() "<<endl);
            _state = _system.realizeTopology();
        MMBLOG_FILE_FUNC_LINE(INFO, "Read a list of included atoms in biopolymer chains. physicsWhereYouWantIt now has included : "<< _dumm.getNumIncludedAtoms () <<" atoms. "<<endl);
    }
    MMBLOG_FILE_FUNC_LINE(INFO, "About to add all residues within the specified radius of of the residues specified in _parameterReader->includeAllResiduesWithinVector. Right now DuMM has "<< _dumm.getNumIncludedAtoms () <<" included atoms. "<<endl);
    // adds to includeAllNonBondAtomsInResidueVector, all residues within the specified radius of the residues specified in _parameterReader->includeAllResiduesWithinVector.
    // this doesn't need a check to make sure it contains something. Because if  _parameterReader->includeAllResiduesWithinVector is empty, this call does nothing:
    MMBLOG_FILE_FUNC_LINE(INFO, "_parameterReader->includeAllResiduesWithinVector.size() = "<<_parameterReader->includeAllResiduesWithinVector.size() << endl);
    #ifdef USE_OPENMM
    for (size_t i = 0 ;  i < _parameterReader->includeAllResiduesWithinVector.size() ; i++){_parameterReader->includeAllResiduesWithinVector[i].print(); }
    // this is where we add residues within a certain radius of a residue of interest:
    _parameterReader->myBiopolymerClassContainer.includeAllResiduesWithin(_parameterReader->includeAllResiduesWithinVector, _parameterReader->physicsContainer.updResidueStretchVector() /*includeAllNonBondAtomsInResidueVector */ , _state);
    MMBLOG_FILE_FUNC_LINE(INFO, "Just finished adding all residues within the specified radius of of the residues specified in _parameterReader->includeAllResiduesWithinVector. Have not actually added _parameterReader->includeAllResiduesWithinVector atoms to DuMM. Right now DuMM has "<< _dumm.getNumIncludedAtoms () <<" included atoms. "<<endl);
    MMBLOG_FILE_FUNC_LINE(INFO, "_parameterReader->includeAllResiduesWithinVector.size() = "<<_parameterReader->includeAllResiduesWithinVector.size() << endl);
    #endif


    if ((_parameterReader->physicsContainer.getNumResidueStretches() /* includeAllNonBondAtomsInResidueVector.size()*/  >  0)) 
    {
        _parameterReader->myBiopolymerClassContainer.includeAllNonBondAtomsInResidues(_parameterReader->physicsContainer.updResidueStretchVector() /*includeAllNonBondAtomsInResidueVector*/ , _state, _dumm);
        // Made a change that affected topology. DuMM can't count atoms unless I update the topology:
        MMBLOG_FILE_FUNC_LINE(INFO, "system.realizeTopology() "<<endl);
        _state = _system.realizeTopology();
        //_system.realize(_state,Stage::Position);
        MMBLOG_FILE_FUNC_LINE(INFO, "Turned on Physics Where You Want It, for the following residues: "<<endl);
        _parameterReader->myBiopolymerClassContainer.printAllIncludedResidues ( _parameterReader->physicsContainer.getResidueStretchVector());   //includeAllNonBondAtomsInResidueVector);
    }
    MMBLOG_FILE_FUNC_LINE(INFO, "_parameterReader->includeAllResiduesWithinVector.size() = "<<_parameterReader->includeAllResiduesWithinVector.size() << endl);
    bool myNonBondedOn = 
    ((_parameterReader->globalCoulombScaleFactor >0) ||
     (_parameterReader->globalVdwScaleFactor > 0));
    MMBLOG_FILE_FUNC_LINE(INFO, "You have Coulomb and/or VdW forces : "<<myNonBondedOn<<" and some sort of PhysicsWhereYouWantIt : "<<myPhysicsWhereYouWantItActive<<endl);
    if (    myPhysicsWhereYouWantItActive && myNonBondedOn )
    {
        MMBLOG_FILE_FUNC_LINE(INFO, "This is fine!"<<endl);
    } else if (myPhysicsWhereYouWantItActive && (!( myNonBondedOn)) ) 
    {
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "You turned on Physics Where You Want It but don't have VdW or electrostatic interactions turned on!"<<endl);
    } else if ((!(myPhysicsWhereYouWantItActive)) && myNonBondedOn) 
    {
        MMBLOG_FILE_FUNC_LINE(WARNING, " you are using non-bonded force field terms, but have not specified a physics zone.  Therefore all atoms in your _system are included in the physics zone.  This could be expensive!"<<endl);
    } else if ((!(myPhysicsWhereYouWantItActive)) && (!( myNonBondedOn)) ) 
    {
        MMBLOG_FILE_FUNC_LINE(INFO, "This is fine!"<<endl);
    } else 
    {
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "Unexplained error!"<<endl); // Just being paranoid..
    }

    //_parameterReader->myBiopolymerClassContainer.printAllIncludedResidues ( _parameterReader->physicsContainer.residueStretchVector ); //includeAllNonBondAtomsInResidueVector);

    ////////////// Reading included water droplet atoms //////////////
    MMBLOG_FILE_FUNC_LINE(INFO, "About to read water droplet atoms, physicsWhereYouWantIt up to now has included : "<< _dumm.getNumIncludedAtoms () <<" atoms. "<<endl);
    _parameterReader->waterDropletContainer.includeAllAtoms(_dumm);
    // this may be redundant .. should just do it once below:
    MMBLOG_FILE_FUNC_LINE(INFO, "system.realizeTopology() "<<endl);
    _state = _system.realizeTopology();

    ////////////// Reading included monoAtoms (e.g. ions) //////////////
    MMBLOG_FILE_FUNC_LINE(INFO, "About to read monoAtoms, physicsWhereYouWantIt up to now has included : "<< _dumm.getNumIncludedAtoms () <<" atoms. "<<endl);
    _parameterReader->myMonoAtomsContainer.includeAllAtoms(_dumm);
    MMBLOG_FILE_FUNC_LINE(INFO, "system.realizeTopology() "<<endl);
    _state = _system.realizeTopology();
    ////////////// Reading moleculeClassContainer //////////////
    MMBLOG_FILE_FUNC_LINE(INFO, "About to read moleculeClassContainer, physicsWhereYouWantIt up to now has included : "<< _dumm.getNumIncludedAtoms () <<" atoms. "<<endl);
    _parameterReader->moleculeClassContainer.includeAllAtoms(_dumm);
    MMBLOG_FILE_FUNC_LINE(INFO, "system.realizeTopology() "<<endl);
    _state = _system.realizeTopology();

    MMBLOG_FILE_FUNC_LINE(INFO, "read directly specified residues + water droplet atoms, physicsWhereYouWantIt now has included : "<< _dumm.getNumIncludedAtoms () <<" atoms. "<<endl);
    //////////////////////////////////////////////////////////////////
    //DuMM::IncludedAtom 
    if (_parameterReader->verbose){
    MMBLOG_FILE_FUNC_LINE(INFO, "Dumping "<<_dumm.getNumIncludedAtoms()<<" charged atom types: "<<endl);
    //for (DuMM::NonbondAtomIndex  nbax (0); nbax<_dumm.getNumNonbondAtoms (); nbax++)
    for (DuMM::IncludedAtomIndex iax1 = DuMM::IncludedAtomIndex(0); iax1<_dumm.getNumIncludedAtoms (); iax1++)
    {
                   
        //DuMM::IncludedAtomIndex iax1 = _dumm.getIncludedAtomIndexOfNonbondAtom(nbax);
        MMBLOG_FILE_FUNC_LINE(INFO, "partial charge:  "<<_dumm.getPartialCharge(iax1)<<" ");//<<endl;
        MMBLOG_FILE_FUNC_LINE(INFO, "atom type name:  "<<_dumm.getChargedAtomName(iax1) <<endl);
    } 
    }
}


// !!!!!! WARNING The destruction of the event's pointers must be checked
void ConstrainedDynamics::createEventHandlers(){
    if (_parameterReader->setTemperature){

        // last param is relax time in ps.  if this number is small, integrator will be slow.  this is the time it takes to enforce temperature.
        // new nose-hoover.  uncomment this when Simbody has been updated. then get rid of velocityrescaling thermostat
        //
        if ((_parameterReader->thermostatType).compare("NoseHoover") == 0 ) { 
            Force::Thermostat  (_forces, _matter,SimTK_BOLTZMANN_CONSTANT_MD,_parameterReader->temperature,_parameterReader->noseHooverTime); //don't forget to set time constant back to 1.0.
            MMBLOG_FILE_FUNC_LINE(INFO, "Created Nose-Hoover thermostat."<<endl);
        }
        else if ((_parameterReader->thermostatType).compare("VelocityRescaling") == 0 ) { 
            VelocityRescalingThermostat * myVelocityRescalingThermostat = new VelocityRescalingThermostat(_system,_parameterReader->temperature, _parameterReader->velocityRescalingInterval);
            _system.updDefaultSubsystem().addEventHandler(myVelocityRescalingThermostat);
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "Created velocity rescaling thermostat."<<endl); // FIXME: Is the CRITICAL here actually correct?
        } else {
          MMBLOG_FILE_FUNC_LINE(CRITICAL, "Unrecognized thermostatType = "<<_parameterReader->thermostatType<<endl);
        }
    }
    //Force::GlobalDamper globalDamper(_forces,_matter,.001); // didn't seem to have an effect.  At .1 damping, it actually was sl wer.

    if (_parameterReader->removeRigidBodyMomentum  ) {
        MassCenterMotionRemover  * myMassCenterMotionRemover = new MassCenterMotionRemover(_system, _parameterReader->removeMomentumPeriod);
        _system.updDefaultSubsystem().addEventHandler (myMassCenterMotionRemover );
        (*myMassCenterMotionRemover).enableMassCenterCorrection(false); // this should turn off the displacement of _system center of mass to origin
        // was:
        //(*myMassCenterMotionRemover).enableMassCenterCorrection(true);
    }
    // this is actually set in PeriodicPdbAndEnergyWriter
    //if (_parameterReader->writeDoublePrecisionTrajectories) {
    //        PdbAtom::setWriteFullPrecisionLocation(true);// write trajectory file with full precision (REMARK-SIMTK-COORDS lines added).
    //}

    PeriodicPdbAndEnergyWriter * myPeriodicPdbWriter = new PeriodicPdbAndEnergyWriter(_system,_dumm,_output,  _parameterReader->reportingInterval, (*_parameterReader),_parameterReader->myBiopolymerClassContainer );
    _system.updDefaultSubsystem().addEventHandler(myPeriodicPdbWriter);

    if (_parameterReader->setForceAndStericScrubber) {
        MMBLOG_FILE_FUNC_LINE(INFO, "_parameterReader->setForceAndStericScrubber = "<<_parameterReader->setForceAndStericScrubber<<", now activating scrubber." <<endl);
       
        PeriodicScrubber  * myPeriodicScrubber  = new PeriodicScrubber (_system,_forces,_dumm,(*_parameterReader),_parameterReader->scrubberPeriod/2);//,contacts,hc );//,  _parameterReader->scrubberPeriod/2,_parameterReader->dutyCycle);
        _system.updDefaultSubsystem().addEventHandler (myPeriodicScrubber );
    }
        
    // Send coordinates to VMD 
    if (_parameterReader->vmdOutput == 1) 
        _system.updDefaultSubsystem().addEventReporter(new PeriodicVmdReporter(
            _system, 0.015, 3000, true));
}

void ConstrainedDynamics::initializeIntegrator(){
    ////////////////////////////////////////////
    /// Create and configure time integrator ///
    ////////////////////////////////////////////

    

    // Integrator *study;
    if(_study) delete _study;
    if ((_parameterReader->integratorType).compare("RungeKuttaMerson")==0) 
        _study = new RungeKuttaMersonIntegrator(_system);
    else 
        if ((_parameterReader->integratorType).compare("Verlet"          )==0)
            _study = new VerletIntegrator(_system);
        else {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "You must specify an integratorType of one of the supported types!!!" << endl);
        }

    if (_parameterReader->useFixedStepSize)
        _study->setFixedStepSize(_parameterReader->integratorStepSize);

    _study->setAccuracy(_parameterReader->integratorAccuracy);
    _study->setConstraintTolerance(_parameterReader->constraintTolerance);//_parameterReader->integratorAccuracy);

    ////////////////////////////////////////////

    // have to realize topology since we've added periodic event handlers and reporters and other system objects.
    MMBLOG_FILE_FUNC_LINE(INFO, "system.realizeTopology() "<<endl);
    _state = _system.realizeTopology();

    ////////////////////////////////////////////
    /// Create and call time stepper         ///
    ////////////////////////////////////////////

    if(_ts) delete _ts;
    _ts = new TimeStepper(_system,*_study);
    // added an initial momentum
    srand((unsigned)time(0));
    if (_parameterReader->setInitialVelocities) 
    for (int i =0; i< _state.getNU(); i++)
            _state.updU()[i] =  rand() * .0001 / RAND_MAX ;
    MMBLOG_FILE_FUNC_LINE(INFO, "Starting dynamics now."<<endl);
    //cout<<"[Repel.h:ConstrainedDynamics] _study.getPredictedNextStepSize() " << (*_study).getPredictedNextStepSize()<<endl;
    MMBLOG_FILE_FUNC_LINE(INFO, "size of _matter subsystem Q vector or number of degrees of freedom: "<<(_matter.getQ(_state)).size()<<endl);
    if (_parameterReader->minimize) {
        if (_parameterReader->basePairContainer.numBasePairs() >0) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "If you want to minimize, you can't have any baseInteraction\'s!"<<endl);
        } 
        MMBLOG_FILE_FUNC_LINE(INFO, "Starting minimizer ..."<<endl);
        LocalEnergyMinimizer::minimizeEnergy(_system, _state, .001);
    } 
    _ts->initialize(_state);
    MMBLOG_FILE_FUNC_LINE(INFO, "Force subsystem topology has been realized: "<<_forces.subsystemTopologyHasBeenRealized()    <<endl);

    MMBLOG_FILE_FUNC_LINE(INFO, "Matter subsystem mass = "<<_matter.calcSystemMass(_state)<<std::endl);
    
    _nextFrame = 1;
    _previousTimeStep = 0.0;

}

void ConstrainedDynamics::postDynamics(){
    ////////////////////////////////////////////
    /// Post-dynamics tasks                  ///
    ////////////////////////////////////////////

    _output.close();

    MMBLOG_FILE_FUNC_LINE(INFO, "_study.getPreviousStepSizeTaken() " << (*_study).getPreviousStepSizeTaken()<<endl);
    MMBLOG_FILE_FUNC_LINE(INFO, "(*_study).getPredictedNextStepSize() " <<(*_study).getPredictedNextStepSize()<<endl);
    MMBLOG_FILE_FUNC_LINE(INFO, "(*_study).getActualInitialStepSizeTaken() " <<(*_study).getActualInitialStepSizeTaken()<<endl);
    MMBLOG_FILE_FUNC_LINE(INFO, "(*_study).getNumStepsTaken() " << (*_study).getNumStepsTaken()<<endl);
    MMBLOG_FILE_FUNC_LINE(INFO, "(*_study).getNumStepsAttempted() " <<(*_study).getNumStepsAttempted()<<endl);
    _state = _ts->getState();
    if (_parameterReader->writeLastFrameFile)
    {
        if ( _parameterReader->useCIFFileFormat )
        {
            //======================================== Initialise internal variables
            gemmi::Structure outStruct;
            const auto& biopolymers                   = _parameterReader->myBiopolymerClassContainer.getBiopolymerClassMap();

            //======================================== Determine structure name and save it
            std::string strName                       = _parameterReader->lastFrameFileName;
            istringstream iss                         ( strName );
            std::vector<std::string> tokens; std::string token;
            while ( std::getline ( iss, token, '.') ) { if ( !token.empty() ) { tokens.push_back ( token ); } }
            if ( tokens.size() > 2 )                  { strName = std::string ( tokens.at(tokens.size()-3) + tokens.at(tokens.size()-2) ); }
            else                                      { strName = "LASTX"; }
            if ( strName[0] == '/' ) { strName.erase(0, 1); }
            outStruct.name                            = strName;

            //======================================== Set output structure cell
            outStruct.spacegroup_hm                   = "P1";

            //==================================== Create new gemmi model for this compound
            gemmi::Model gModel                       ( "1" );

            SimTK::CIFOut::buildModel                 ( _state, gModel, biopolymers, _system, 17 );
            outStruct.models.emplace_back             ( std::move( gModel ) );

            gemmi::setup_entities                     ( outStruct );
            gemmi::assign_label_seq_id                ( outStruct, true );
            gemmi::assign_subchains                   ( outStruct, true );

            SimTK::CIFOut::assignEntities             ( outStruct, biopolymers );

            //======================================== Write out CIF
            SimTK::CIFOut::writeOutCif                ( outStruct, _parameterReader->lastFrameFileName, _parameterReader->lastFileRemarks );
        }
        else
        {
            PdbAtom::setWriteFullPrecisionLocation(true);// write last frame file with full precision (REMARK-SIMTK-COORDS lines added).
            ofstream lastFrameFile(_parameterReader->lastFrameFileName.c_str());
            for (SimTK::CompoundSystem::CompoundIndex c(0); c < _system.getNumCompounds(); ++c)
                (_system.getCompound(c)).writePdb(_state, lastFrameFile,Transform(Vec3(0)));
            lastFrameFile.close();
        }
    }
    ////////////////////////////////////////////

    if(_ts) { 
        delete _ts;
        _ts = NULL;
    }
    if(_study) { 
        delete _study;
        _study = NULL;
    }

    MMBLOG_FILE_FUNC_LINE(INFO, "Stage "<<_parameterReader->currentStage<<" completed successfully! "<<endl);
    GlobalProgressWriter::get().update(ProgressWriter::State::FINISHED);
}

void ConstrainedDynamics::runAllSteps() {
    if(getRemainingFramesNum() <= 0){
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "Cannot run more steps" << endl);
    }

    // Normal way to compute dynamics
    //_parameterReader->myBiopolymerClassContainer.printBiopolymerInfo();
    // _ts->stepTo(_parameterReader->numReportingIntervals*_parameterReader->reportingInterval);
    MMBLOG_FILE_FUNC_LINE(INFO, "About to run time integrator for "<<_parameterReader->numReportingIntervals<< " reporting intervals, each of duration "<<_parameterReader->reportingInterval<<" ps "<<endl);

    // New way, to allow a step by step control later
        
    for(; _nextFrame <= _parameterReader->numReportingIntervals; _nextFrame++){
        //MMBLOG_FILE_FUNC_LINE(std::endl; 
        _ts->stepTo(_nextFrame*_parameterReader->reportingInterval);
        _state = _ts->getState();
        //MMBLOG_FILE_FUNC_LINE(std::endl; 
        if(_parameterReader->detectConvergence && _parameterReader->converged){
            //MMBLOG_FILE_FUNC_LINE(std::endl; 
            _nextFrame = _parameterReader->numReportingIntervals;
        }
    }

}

unsigned int ConstrainedDynamics::runOneStep(){
    if(getRemainingFramesNum() <= 0){
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "Cannot run more steps" << endl);
    }

    double nextStep = _previousTimeStep + _parameterReader->reportingInterval;
    _ts->stepTo(nextStep);
    _state = _ts->getState();

    _previousTimeStep = nextStep;
    _nextFrame ++;

    return getRemainingFramesNum();
}

void ConstrainedDynamics::initializeBodies(){
    //MMBLOG_FILE_FUNC_LINE(endl;
    initializeMoleculesAndBonds();
    MMBLOG_FILE_FUNC_LINE(INFO, endl);
    #ifdef USE_OPENMM
    setInterfaceMobilizers();
    #endif
    setMobilizers();
    //_parameterReader->removeBasePairsAcrossRigidStretches(); //SCF
    createMultibodyTree();
}

void ConstrainedDynamics::initializeDynamics(){
    // This is where we call mobilizerContainer.addMobilizerStretchesToVector
    // This is where initializeAtomInfoVectors is called (at least when density map fitting is active): 
    // but note that initializeMoleculesAndBonds() precedes createState()
    MMBLOG_FILE_FUNC_LINE(INFO, endl);
    //_parameterReader->myBiopolymerClassContainer.printAtomInfoVector(); //   at this point ..
    initializeCustomForcesConstraints();
    MMBLOG_FILE_FUNC_LINE(INFO, endl);
    createEventHandlers();    
    MMBLOG_FILE_FUNC_LINE(INFO, endl);
    initializeIntegrator();    
    MMBLOG_FILE_FUNC_LINE(INFO, "Net charge over "<<_dumm.getNumIncludedAtoms ()<<" atoms included in physics zone = "<< _dumm.getTotalIncludedCharge()<<endl);
}

void ConstrainedDynamics::writeMMBPDB(std::ofstream & filestream){
    CompoundSystem system;
    State state;
    SimbodyMatterSubsystem matter(system);
    DuMMForceFieldSubsystem dumm(system);
    _parameterReader->configureDumm(dumm);
    MMBLOG_FILE_FUNC_LINE(INFO, endl);
    initializeBiopolymersAndCustomMolecules(system);
    initializeMoleculesAndBonds(system, dumm, matter);
    //MMBLOG_FILE_FUNC_LINE(endl;
    //MMBLOG_FILE_FUNC_LINE(endl;
    // setInterfaceMobilizers(system, matter, state);
    // setMobilizers();
    createMultibodyTree(system, state);
    _parameterReader->myBiopolymerClassContainer.writePdb(state, system, filestream);
}

void ConstrainedDynamics::forceAdjustmentsWithFinalMobilizers(){
    MMBLOG_FILE_FUNC_LINE(INFO, " About to start myBiopolymerClassContainer.selectivelyRemoveRigidMobilizerStretchesFromResidueStretchContainer"<<endl);
    if (_parameterReader->removeDensityForcesFromRigidStretches) {_parameterReader->myBiopolymerClassContainer.selectivelyRemoveRigidMobilizerStretchesFromResidueStretchContainer(_parameterReader->mobilizerContainer, _parameterReader->densityContainer);}    
    MMBLOG_FILE_FUNC_LINE(INFO, " About to start removeBasePairsAcrossRigidStretches"<<endl);
    if (_parameterReader->setRemoveBasePairsAcrossRigidStretches) {_parameterReader->removeBasePairsAcrossRigidStretches();}    
    MMBLOG_FILE_FUNC_LINE(INFO, " About to start basePairContainer.addHelicalStacking"<<endl);
    if (_parameterReader->setHelicalStacking){_parameterReader->basePairContainer.addHelicalStacking(_parameterReader->myBiopolymerClassContainer, _parameterReader->_leontisWesthofClass, _parameterReader->ntc_par_class,_parameterReader->ntc_class_container);}
}

void ConstrainedDynamics::runDynamics() {
    MMBLOG_FILE_FUNC_LINE(INFO, endl);
    if (initializeBiopolymersAndCustomMolecules()){
        MMBLOG_FILE_FUNC_LINE(WARNING, "Returned an error from initializeBiopolymersAndCustomMolecules"<<endl);
        //exit(1) ;
    }
    MMBLOG_FILE_FUNC_LINE(INFO, endl);
    initializeBodies();
    //MMBLOG_FILE_FUNC_LINE(endl;
     
    //MMBLOG_FILE_FUNC_LINE(endl;
    //_parameterReader->myBiopolymerClassContainer.printAtomInfoVector(); //  Looks fine at this point ..
    // This should be done after initializeBodies() because that is when we are reverting residues back to Default BondMobility
    forceAdjustmentsWithFinalMobilizers();
    initializeDynamics();

    runAllSteps();

    postDynamics();
};

