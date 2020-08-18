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
#include "Utils.h"
#include <stdio.h>
#include <string.h>
#include   "MonoAtoms.h"
#include "BiopolymerClass.h"
#include "ConstraintContainer.h"
#include "ElectrostaticPotentialGridForce.h"
#include "PeriodicPdbAndEnergyWriter.h"

#ifdef CPP4_MAPS_USAGE
  #include <mmdb2/mmdb_manager.h>
  #include <mmdb2/mmdb_cifdefs.h>
#endif

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
	    ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" You are trying to weld an atom of a custom molecule (of chain ID: >"<<myChainID <<"<) to another atom. Unfortunately this feature is not yet supported. You should, however, be able to weld a custom molecule to Ground."<<endl;
	    ErrorManager::instance.treatError();
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
	cout<<__FILE__<<":"<<__LINE__<< endl;
        //#ifdef BUILD_MMB_SHARED_LIB
        //#pragma message ("BUILD_MMB_SHARED_LIB is currently defined")
        //#else
        //#pragma message ("BUILD_MMB_SHARED_LIB is NOT currently defined")
        //#endif
	cout<<__FILE__<<":"<<__LINE__<< endl;
 
        #ifdef USE_OPENMM
        #pragma message ("USE_OPENMM is defined")
        #else 
        #pragma message ("USE_OPENMM is NOT defined")
        #endif 
        cout<<__FILE__<<":"<<__LINE__<< endl;
        // here we turn the list of interface constraints from myParameterReader.constraintToGroundContainer.interfaceContainer into pairs of constraints, at most one pair for each pair of constrained chains.
        #ifdef USE_OPENMM    
        myParameterReader.constraintToGroundContainer.addSingleWeldConstraintPerInterfaceChainPair ( myBiopolymerClassContainer);
        #endif
        cout<<__FILE__<<":"<<__LINE__<<endl;
        myParameterReader.constraintToGroundContainer.applyConstrainChainRigidSegments ( myBiopolymerClassContainer,  system,  matter, state);
        
        cout<<__FILE__<<":"<<__LINE__<<" Now printing all constraints :"<<endl;
        myParameterReader.constraintToGroundContainer.printConstraintClasses();
        cout<<__FILE__<<":"<<__LINE__<<endl;
        vector<Constraint> myWeld2;        
        vector<Force> myLinearBushing;
       
        MobilizedBody myMobilizedBody1, myMobilizedBody2;
        Transform myTransform1, myTransform2;
        
        for (int i = 0; i < myParameterReader.constraintToGroundContainer.numConstraintClasses(); i++){
            ConstraintClass  myConstraintClass = myParameterReader.constraintToGroundContainer.getConstraintClass(i);
            cout<<__FILE__<<":"<<__LINE__<<" Applying constraint index : "<<i<<endl;
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
                    ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" You are trying to constriain chain "<<myChain<<" which is neither a biopolymer nor a custom molecule"<<endl;
                    ErrorManager::instance.treatError();
                }
            } else if (myConstraintClass.getConstraintType() == WeldToAtom) { // constrain two atoms to each other
                String myAtomName = myConstraintClass.getAtomName1();
                ResidueID myResidueID = myConstraintClass.getResidueID1();
                String    myChain     = myConstraintClass.getChain1();
                // SCF : check whether we are dealing with biopolymer, monoAtom, or custom molecule, then generalize to handle all three. No further change should be necessary:
                
                //MobilizedBody getMobilizedBody(SimbodyMatterSubsystem & myMatter, ParameterReader & parameterReader, BiopolymerClassContainer & biopolymerClassContainer, String myChainID, ResidueID myResidue, String myAtom) {
                cout<<__FILE__<<":"<<__LINE__<<endl;
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
                    ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" You are trying to weld a mobilized body to itself.  This could happen, for instance, if you are welding two atoms from the same rigid stretch to each other."<<endl;
                    ErrorManager::instance.treatError();
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
			ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" : You have specified an unsupported ConstraintType : "<<myConstraintClass.getConstraintType() <<endl;
			ErrorManager::instance.treatError();
		    }


        } // of for

        for (int q=0;q<(int)myParameterReader.baseOperationVector.size();q++) {
            if (((myParameterReader.baseOperationVector[q]).BasePairIsTwoTransformForce).compare("constraint" ) == 0) {
                ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" : An error has occurred!  Constraints are no longer specified with this data structure. "<<endl;
                ErrorManager::instance.treatError();
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
                {ErrorManager::instance <<__FILE__<<":"<<__LINE__<<"Error! The first chain in the constraint or restraint was neither a biopolymer nor a monoAtoms object, and you are not constraining to ground."<<endl; ErrorManager::instance.treatError();}

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
                {ErrorManager::instance <<__FILE__<<":"<<__LINE__<<"Error! The second chain in the constraint or restraint was neither a biopolymer nor a monoAtoms object."<<endl; ErrorManager::instance.treatError();}

               cout<<__FILE__<<":"<<__LINE__<<" system.realizeTopology() "<<endl;
               state = system.realizeTopology();                                        
               system.realize(state,Stage::Position);

               std::cout<<__FILE__<<":"<<__LINE__<<" : Enforcing constraint or restraint from : "<<string(myParameterReader.baseOperationVector[q].FirstBPChain) <<", " << myParameterReader.baseOperationVector[q].FirstBPResidue.outString() << " to : " << string( myParameterReader.baseOperationVector[q].SecondBPChain)<<", "<<myParameterReader.baseOperationVector[q].SecondBPResidue.outString()<<  std::endl;
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
                {ErrorManager::instance <<__FILE__<<":"<<__LINE__<<"Error! You have specified an unsupported constraint or restraint."<<endl; ErrorManager::instance.treatError();}

                cout<<__FILE__<<":"<<__LINE__<<" system.realizeTopology() "<<endl;
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

int  ConstrainedDynamics::initializeBiopolymersAndCustomMolecules(CompoundSystem & system){
    //================================================ Depending on settings, set the PDB file name to CIF instead of PDB
    _parameterReader->myBiopolymerClassContainer.resetAllPdbFileNames ( _parameterReader->previousFrameFileName );
    
    _parameterReader->moleculeClassContainer.initializeCompounds (_dumm);
    cout<<__FILE__<<":"<<__LINE__<<endl;
    _parameterReader->moleculeClassContainer.matchDefaultConfiguration (_parameterReader->readPreviousFrameFile, _parameterReader->previousFrameFileName, _parameterReader->matchExact,  _parameterReader->matchIdealized  );
    cout<<__FILE__<<":"<<__LINE__<<endl;
    _parameterReader->moleculeClassContainer.adoptCompounds (system);
    cout<<__FILE__<<":"<<__LINE__<<endl;

    //cout<<__FILE__<<":"<<__LINE__<<"This is temporary .. SCF"<<endl;
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
    cout<<__FILE__<<":"<<__LINE__<<" returnValue = "<<returnValue<<std::endl;
    if (returnValue)
    {
        cout<<__FILE__<<":"<<__LINE__<<" Warning : There was an error in initializeBiopolymers"<<endl;
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

    //cout<<__FILE__<<":"<<__LINE__<<endl;
    _parameterReader->waterDropletContainer.addWaterMolecules(system, dumm, matter, _parameterReader->myBiopolymerClassContainer);
    //cout<<__FILE__<<":"<<__LINE__<<endl;
    _parameterReader->waterDropletContainer.validateWaterVectors();
    //cout<<__FILE__<<":"<<__LINE__<<endl;
    _parameterReader->waterDropletContainer.matchDefaultConfiguration (_parameterReader->readPreviousFrameFile, _parameterReader->previousFrameFileName, _parameterReader->matchExact,  _parameterReader->matchIdealized  );
    //cout<<__FILE__<<":"<<__LINE__<<endl;
    _parameterReader->waterDropletContainer.adopt (system,_parameterReader->readPreviousFrameFile);
    //cout<<__FILE__<<":"<<__LINE__<<endl;
    // Moved this from ParameterReader.cpp:
    for (int i = 0; i < _parameterReader->additionalCovalentBondVector.size(); i++) {
        CovalentBondClass myBond = _parameterReader->additionalCovalentBondVector[i];
        if (myBond.getChain1().compare(myBond.getChain2()) != 0 ){
         ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" : Chain IDs of both atoms must be the same."<<endl;
             ErrorManager::instance.treatError();
        }
	if (_parameterReader->myBiopolymerClassContainer.hasChainID(myBond.getChain1())) {
	    _parameterReader->myBiopolymerClassContainer.updBiopolymerClass(myBond.getChain1() ).addRingClosingBond(myBond);
	} else if (_parameterReader->moleculeClassContainer.hasChainID(myBond.getChain1())) {
            _parameterReader->moleculeClassContainer.updMoleculeClass(myBond.getChain1()).addRingClosingBond(myBond);
        } else {
            ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" Unexplained error!"<<endl;
            ErrorManager::instance.treatError();
        }

    }
    cout<<__FILE__<<":"<<__LINE__<<" If you want to close the cysteine bridges, you can issue:"<<std::endl;
    //_parameterReader->myBiopolymerClassContainer.createDisulphideBridges(std::cout); // Should make this user controllable. // there is an issue with immediately mutating.. consider adding to mutation vector..
}

#ifdef USE_OPENMM
void ConstrainedDynamics::setInterfaceMobilizers(){
    cout<<__FILE__<<":"<<__LINE__<<endl;
    setInterfaceMobilizers(_system, _matter, _state);
}

void ConstrainedDynamics::setInterfaceMobilizers(CompoundSystem & system, SimbodyMatterSubsystem & matter, State & state){
    cout<<__FILE__<<":"<<__LINE__<<" Checking interface mobilizers .."<<endl;
    _parameterReader->myBiopolymerClassContainer.validateAtomInfoVectors();
    _parameterReader->mobilizerContainer.addStretchesToVectorFromInterfaceContainer(_parameterReader->myBiopolymerClassContainer);
    cout<<__FILE__<<":"<<__LINE__<<" Checking domains interface mobilizers .."<<endl;
    _parameterReader->mobilizerContainer.addMobilizerDomainsInterfacesToVector(_parameterReader->mobilizerDomainsInterfaceVector, _parameterReader->myBiopolymerClassContainer);
    cout<<__FILE__<<":"<<__LINE__<<" Domains interface mobilizers checked"<<endl;
    _parameterReader->mobilizerContainer.printMobilizerStretches();
    cout<<__FILE__<<":"<<__LINE__<<" system.getNumRealizeCalls () = "<< system.getNumRealizeCalls () <<endl;
    cout<<__FILE__<<":"<<__LINE__<<" system.realizeTopology() "<<endl;
    state = system.realizeTopology();
    cout<<__FILE__<<":"<<__LINE__<<endl;
    system.realize(state,Stage::Position);
}
#endif
void ConstrainedDynamics::setMobilizers()
{
    cout<<__FILE__<<":"<<__LINE__<<endl; 
    #ifdef USE_OPENMM
    _parameterReader->mobilizerContainer.createMobilizersWithin(_parameterReader->myBiopolymerClassContainer,_state);
    #endif
    cout<<__FILE__<<":"<<__LINE__<<" Done adding MobilizerStretch's .  However we haven't yet dealt with 'Default' MobilizerStretch's. "<<endl;
    cout<<__FILE__<<":"<<__LINE__<<" Printing all mobilizer stretches: "<<endl;
    _parameterReader->mobilizerContainer.printMobilizerStretches();

    cout<<__FILE__<<":"<<__LINE__<<" About to clear MobilizerStretch's that have been overridden, or specified to remain at BondMobility 'Default' "<<endl;

    // removing all MobilizerStretch's that modify residues specified to be at BondMobility = "Default"
    for (int i = 0; i < _parameterReader->mobilizerContainer.getNumResidueStretches(); i++) {
        cout<<__FILE__<<":"<<__LINE__<<" Examining mobilizer stretch "<<i<<" to see if it has bondMobility Default :"<<endl;
        _parameterReader->mobilizerContainer.residueStretchVector[i].printStretch();
        if (_parameterReader->mobilizerContainer.residueStretchVector[i].getBondMobilityString().compare("Default") == 0){
            ResidueStretch tempResidueStretch = _parameterReader->mobilizerContainer.getResidueStretch(i);
            cout<<__FILE__<<":"<<__LINE__<<" Detected that BondMobility : >"<<_parameterReader->mobilizerContainer.residueStretchVector[i].getBondMobilityString()<<"< is Default."<<endl;
            _parameterReader->myBiopolymerClassContainer.updBiopolymerClass(tempResidueStretch.getChain() ).selectivelyRemoveResidueStretchFromContainer(tempResidueStretch,_parameterReader->mobilizerContainer);
            i--; // Compensating for the fact that the Default stretch i has been deleted.
            
        } else {
            cout<<__FILE__<<":"<<__LINE__<<" Detected that BondMobility : >"<<_parameterReader->mobilizerContainer.residueStretchVector[i].getBondMobilityString()<<"< is NOT Default."<<endl;
        }    
    }  
    // Done dealing with "Default" mobilizer stretches
    cout<<__FILE__<<":"<<__LINE__<<" Done dealing with MobilizerStretch's specified to remain at BondMobility 'Default' "<<endl;
    cout<<__FILE__<<":"<<__LINE__<<" Printing all mobilizer stretches: "<<endl;
    _parameterReader->mobilizerContainer.printMobilizerStretches();
    
    _parameterReader->mobilizerContainer.setBiopolymerBondMobility(_parameterReader->myBiopolymerClassContainer);
    _parameterReader->myBiopolymerClassContainer.setSingleBondMobility(_parameterReader->mobilizerContainer.singleBondMobilityVector);
}

void ConstrainedDynamics::createMultibodyTree()
{
    cout<<__FILE__<<":"<<__LINE__<<endl;
    createMultibodyTree(_system, _state);
}

void ConstrainedDynamics::createMultibodyTree(CompoundSystem & system, State & state){
    cout<<__FILE__<<":"<<__LINE__<<" Calling modelOneCompound "<< system.getNumCompounds()<<" times."<<endl;
    for (CompoundSystem::CompoundIndex c(0); c < system.getNumCompounds(); ++c) {
        String currentChain = String(system.getCompound(c).getPdbChainId());
        if (_parameterReader->myBiopolymerClassContainer.hasChainID(currentChain)){
            String tempMobilizerType =  _parameterReader->myBiopolymerClassContainer.updBiopolymerClass(currentChain).getFirstResidueMobilizerType();
        cout<<__FILE__<<":"<<__LINE__<<" Attaching biopolymer chain "<<currentChain<<" using mobilizer type = "<<tempMobilizerType<<endl;
            system.modelOneCompound(c,tempMobilizerType  );
        } else {
            String tempMobilizerType ("Free");
        cout<<__FILE__<<":"<<__LINE__<<" Attaching non-biopolymer chain "<<currentChain<<" using mobilizer type = "<<tempMobilizerType<<endl;
            system.modelOneCompound(c, tempMobilizerType );  // for non-biopolymers, use the default behavior, which is for root atom to be connected with a Free mobilizer.
        }
    }
    cout<<__FILE__<<":"<<__LINE__<<" system.realizeTopology() "<<endl;
    state = system.realizeTopology();
    system.realize(state,Stage::Position);
}

void ConstrainedDynamics::initializeCustomForcesConstraints(){
    _output.open(_parameterReader->outTrajectoryFileName.c_str(),  ios_base::app);

    _parameterReader->waterDropletContainer.addTethers(_parameterReader->atomSpringContainer);

    _parameterReader->waterDropletContainer.validateWaterVectors();
    cout<<__FILE__<<":"<<__LINE__<<endl; 
    _parameterReader->myBiopolymerClassContainer.multiplySmallGroupInertia( _parameterReader->smallGroupInertiaMultiplier, _system,_matter, _state      );
    _parameterReader->waterDropletContainer.multiplySmallGroupInertia( _parameterReader->waterInertiaMultiplier, _system,_matter, _state      );
    //cout<<__FILE__<<":"<<__LINE__<<endl;
    #ifdef USE_OPENMM  
    _parameterReader->myBiopolymerClassContainer.initializeAtomInfoVectors(_matter,_dumm); // investigate moving this outside the conditional, so it's available for DensityForce below.
    #endif
    cout<<__FILE__<<":"<<__LINE__<<" You have specified _parameterReader->physicsRadius = "<<_parameterReader->physicsRadius<<endl;
    //_parameterReader->myBiopolymerClassContainer.printAtomInfoVector(); // Looks fine at this point ..
    if ( _parameterReader->physicsRadius > 0.0000001) 
    {
        cout<<__FILE__<<":"<<__LINE__<<endl;
        _parameterReader->myBiopolymerClassContainer.physicsZone (_parameterReader->includeAllResiduesWithinVector, _parameterReader->physicsRadius,_matter,_state); 
    }
    cout<<__FILE__<<":"<<__LINE__<< "_parameterReader->includeAllResiduesWithinVector.size() = "<<_parameterReader->includeAllResiduesWithinVector.size() << endl;
    #ifdef USE_OPENMM
    // Act on physicsInterfaces command. Turn the interfaces into specific lists of residues:
    _parameterReader->  physicsContainer.addStretchesToVectorFromInterfaceContainer(_parameterReader->myBiopolymerClassContainer);
    #endif
    cout<<__FILE__<<":"<<__LINE__<< "_parameterReader->includeAllResiduesWithinVector.size() = "<<_parameterReader->includeAllResiduesWithinVector.size() << endl;
    // SCF this is a good place to insert includeInterChainInterface processing. 
    if (_parameterReader->includeIntraChainInterfaceVector.size() >0) 
    {
        //_parameterReader->myBiopolymerClassContainer.initializeAtomInfoVectors(_matter,_dumm); // investigate moving this outside the conditional, so it's available for DensityForce below.
       for (int i = 0; i < _parameterReader->includeIntraChainInterfaceVector.size() ; i++) 
       {
           cout<<__FILE__<<":"<<__LINE__<<endl;
           _parameterReader->myBiopolymerClassContainer.addIntraChainInterfaceResidues(_parameterReader->includeIntraChainInterfaceVector[i].Chain,  
               _parameterReader->physicsContainer.residueStretchVector //includeAllNonBondAtomsInResidueVector 
               , _parameterReader->includeIntraChainInterfaceVector[i].Depth ,  
	       _matter,_state) ;
           cout<<__FILE__<<":"<<__LINE__<<endl;
       }
    }
    if (_parameterReader->includeIntraChainInterfaceVector.size() >0) {
        //_parameterReader->myBiopolymerClassContainer.initializeAtomInfoVectors(_matter,_dumm); // investigate moving this outside the conditional, so it's available for DensityForce below.
    for (int i = 0; i < _parameterReader->includeIntraChainInterfaceVector.size() ; i++) {
        _parameterReader->myBiopolymerClassContainer.addIntraChainInterfaceResidues(_parameterReader->includeIntraChainInterfaceVector[i].Chain,  
	    _parameterReader->physicsContainer.residueStretchVector , //includeAllNonBondAtomsInResidueVector , 
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
    cout<<__FILE__<<":"<<__LINE__<<endl;
    if (_parameterReader->densityContainer.numDensityStretches() > 0) 
    {
        cout<<__FILE__<<":"<<__LINE__<<endl;
        _parameterReader->myDensityMap.setNoiseTemperature(_parameterReader->densityNoiseTemperature);
        cout<<__FILE__<<":"<<__LINE__<<" setting noiseScale to :" <<_parameterReader->densityNoiseScale<<endl;
        _parameterReader->myDensityMap.setNoiseScale(_parameterReader->densityNoiseScale);
        cout<<__FILE__<<":"<<__LINE__<<" set noiseScale to :" <<_parameterReader->myDensityMap.getNoiseScale() <<endl;
        _parameterReader->myDensityMap.loadParametersAndDensity(_parameterReader->densityFileName);
        cout<<__FILE__<<":"<<__LINE__<<endl;
        if (_parameterReader->myDensityMap.getNoiseScale() > 0.0000001) {
            cout<<__FILE__<<":"<<__LINE__<<endl;
            _parameterReader->myDensityMap.populateNoiseMap();
            _parameterReader->myDensityMap.writeDensityMapXplor("./noise.xplor", 0, 1);   // write noise only
            _parameterReader->myDensityMap.writeDensityMapXplor("./density.xplor", 1,0);  // write noiseless density only
            _parameterReader->myDensityMap.writeDensityMapXplor("./noisyMap.xplor", 1,1); // write the sum of density + noise
            if (_parameterReader->densityNoiseComputeAutocorrelation){
                std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" First, we write out the noise autocorrelation: "<< std::endl;
                _parameterReader->myDensityMap.densityAutocorrelation(1,0); // Arguments are : calculate noise correlation = 1, calculate density corrleation = 0
                std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Second, we write out the autocorrelation of the input density map: "<< std::endl;
                _parameterReader->myDensityMap.densityAutocorrelation(0,1);
                std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<< std::endl;
            }
        }
        cout<<__FILE__<<":"<<__LINE__<<endl;
        _parameterReader->myDensityMap.precomputeGradient();
        cout<<__FILE__<<":"<<__LINE__<<endl;
        _parameterReader->myDensityMap.precomputeGradientDerivatives();
        cout<<__FILE__<<":"<<__LINE__<<endl;
        cout<<__FILE__<<":"<<__LINE__<<endl;
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
    cout<<__FILE__<<":"<<__LINE__<<endl;
    //_parameterReader->myBiopolymerClassContainer.printAtomInfoVector(); //  at this point.. 
    cout<<__FILE__<<":"<<__LINE__<< "_parameterReader->includeAllResiduesWithinVector.size() = "<<_parameterReader->includeAllResiduesWithinVector.size() << endl;
    cout<<__FILE__<<":"<<__LINE__<<" system.realizeTopology() "<<endl;
    _state = _system.realizeTopology();
    _system.realize(_state,Stage::Position);

    cout<<__FILE__<<":"<<__LINE__<< "_parameterReader->includeAllResiduesWithinVector.size() = "<<_parameterReader->includeAllResiduesWithinVector.size() << endl;


    cout<<__FILE__<<":"<<__LINE__<<" Did you intend to use PhysicsWhereYouWantIt? Checking size of _parameterReader->physicsContainer.getNumResidueStretches() " //includeAllNonBondAtomsInResidueVector.size() = "
        <<_parameterReader->physicsContainer.getNumResidueStretches()<<endl;
    cout<<__FILE__<<":"<<__LINE__<<" About to work on physicsWhereYouWantIt. DuMM now has included : "<< _dumm.getNumIncludedAtoms () <<" atoms. "<<endl;


    cout<<__FILE__<<":"<<__LINE__<<endl;
    //_parameterReader->myBiopolymerClassContainer.printAtomInfoVector(); //  at this point.. 


    bool myPhysicsWhereYouWantItActive =  
       ((_parameterReader->includeNonBondAtomInBiopolymerVector.size() > 0) ||
    (_parameterReader->physicsContainer.getNumResidueStretches() /*includeAllNonBondAtomsInResidueVector.size()*/ > 0) ||
        (_parameterReader->includeAllResiduesWithinVector.size() > 0)) ;
 
    if (myPhysicsWhereYouWantItActive) {
        cout<<__FILE__<<":"<<__LINE__<<" physicsWhereYouWantIt now has included : "<< _dumm.getNumIncludedAtoms () <<" atoms. Clearing list .."<<endl;
        _dumm.clearIncludedNonbondAtomList();
        cout<<__FILE__<<":"<<__LINE__<<" system.realizeTopology() "<<endl;
        _state = _system.realizeTopology();
        _system.realize(_state,Stage::Position);
        cout<<__FILE__<<":"<<__LINE__<<" Cleared non-bonded atom list. Realized _state to the position stage. physicsWhereYouWantIt now has included : "<< _dumm.getNumIncludedAtoms () <<" atoms. Note that this number might not be zero -- we do not clear the bonded atom list."<<endl;
        // Could end this loop here. The only point is to clear the included atom list, in cases that any physicsWhereYouWantIt is active.
    }  
    else {
        // Do nothing.  Don't clear included atom list. Non-bonded forces are off, so the simulation should be cheap.
    // Note that bonded forces may still be more than needed depending on mobilizers.  Worth checking to see if futher optimization is needed.
    }
   
    cout<<__FILE__<<":"<<__LINE__<<endl;
    //_parameterReader->myBiopolymerClassContainer.printAtomInfoVector(); //  at this point.. 


    if (_parameterReader->includeNonBondAtomInBiopolymerVector.size() > 0) 
    {
        cout<<__FILE__<<":"<<__LINE__<<" About to read a list of included atoms in biopolymer chains. physicsWhereYouWantIt now has included : "<< _dumm.getNumIncludedAtoms () <<" atoms. "<<endl;
        _parameterReader->myBiopolymerClassContainer.includeNonBondAtoms(  _parameterReader->includeNonBondAtomInBiopolymerVector,  _state,  _dumm);
            cout<<__FILE__<<":"<<__LINE__<<" system.realizeTopology() "<<endl;
            _state = _system.realizeTopology();
        cout<<__FILE__<<":"<<__LINE__<<" Read a list of included atoms in biopolymer chains. physicsWhereYouWantIt now has included : "<< _dumm.getNumIncludedAtoms () <<" atoms. "<<endl;
    }
    cout<<__FILE__<<":"<<__LINE__<<" About to add all residues within the specified radius of of the residues specified in _parameterReader->includeAllResiduesWithinVector. Right now DuMM has "<< _dumm.getNumIncludedAtoms () <<" included atoms. "<<endl;
    // adds to includeAllNonBondAtomsInResidueVector, all residues within the specified radius of the residues specified in _parameterReader->includeAllResiduesWithinVector.
    // this doesn't need a check to make sure it contains something. Because if  _parameterReader->includeAllResiduesWithinVector is empty, this call does nothing:
    cout<<__FILE__<<":"<<__LINE__<< "_parameterReader->includeAllResiduesWithinVector.size() = "<<_parameterReader->includeAllResiduesWithinVector.size() << endl;
    #ifdef USE_OPENMM
    for (int i = 0 ;  i < _parameterReader->includeAllResiduesWithinVector.size() ; i++){_parameterReader->includeAllResiduesWithinVector[i].print(); }
    // this is where we add residues within a certain radius of a residue of interest:
    _parameterReader->myBiopolymerClassContainer.includeAllResiduesWithin(_parameterReader->includeAllResiduesWithinVector, _parameterReader->physicsContainer.residueStretchVector /*includeAllNonBondAtomsInResidueVector */ , _state); 
    cout<<__FILE__<<":"<<__LINE__<<" Just finished adding all residues within the specified radius of of the residues specified in _parameterReader->includeAllResiduesWithinVector. Have not actually added _parameterReader->includeAllResiduesWithinVector atoms to DuMM. Right now DuMM has "<< _dumm.getNumIncludedAtoms () <<" included atoms. "<<endl;
    cout<<__FILE__<<":"<<__LINE__<< "_parameterReader->includeAllResiduesWithinVector.size() = "<<_parameterReader->includeAllResiduesWithinVector.size() << endl;
    #endif


    if ((_parameterReader->physicsContainer.getNumResidueStretches() /* includeAllNonBondAtomsInResidueVector.size()*/  >  0)) 
    {
        _parameterReader->myBiopolymerClassContainer.includeAllNonBondAtomsInResidues(_parameterReader->physicsContainer.residueStretchVector /*includeAllNonBondAtomsInResidueVector*/ , _state, _dumm);
        // Made a change that affected topology. DuMM can't count atoms unless I update the topology:
        cout<<__FILE__<<":"<<__LINE__<<" system.realizeTopology() "<<endl;
        _state = _system.realizeTopology();
        //_system.realize(_state,Stage::Position);
        cout<<__FILE__<<":"<<__LINE__<<" Turned on Physics Where You Want It, for the following residues: "<<endl;
        _parameterReader->myBiopolymerClassContainer.printAllIncludedResidues ( _parameterReader->physicsContainer.residueStretchVector);   //includeAllNonBondAtomsInResidueVector);
    }
    cout<<__FILE__<<":"<<__LINE__<< "_parameterReader->includeAllResiduesWithinVector.size() = "<<_parameterReader->includeAllResiduesWithinVector.size() << endl;
    bool myNonBondedOn = 
    ((_parameterReader->globalCoulombScaleFactor >0) ||
     (_parameterReader->globalVdwScaleFactor > 0));
    cout<<__FILE__<<":"<<__LINE__<<" You have Coulomb and/or VdW forces : "<<myNonBondedOn<<" and some sort of PhysicsWhereYouWantIt : "<<myPhysicsWhereYouWantItActive<<endl;
    if (    myPhysicsWhereYouWantItActive && myNonBondedOn )
    {
             cout<<__FILE__<<":"<<__LINE__<<" This is fine!"<<endl;
    } else if (myPhysicsWhereYouWantItActive && (!( myNonBondedOn)) ) 
    {
        ErrorManager::instance<<__FILE__<<":"<<__LINE__<<" You turned on Physics Where You Want It but don't have VdW or electrostatic interactions turned on!"<<endl;ErrorManager::instance.treatError();
    } else if ((!(myPhysicsWhereYouWantItActive)) && myNonBondedOn) 
    {
        cout<<__FILE__<<":"<<__LINE__<<" Warning!  you are using non-bonded force field terms, but have not specified a physics zone.  Therefore all atoms in your _system are included in the physics zone.  This could be expensive!"<<endl;
    } else if ((!(myPhysicsWhereYouWantItActive)) && (!( myNonBondedOn)) ) 
    {
             cout<<__FILE__<<":"<<__LINE__<<" This is fine!"<<endl;
    } else 
    {
        ErrorManager::instance<<__FILE__<<":"<<__LINE__<<" Unexplained error!"<<endl;ErrorManager::instance.treatError(); // Just being paranoid..
    }

    //_parameterReader->myBiopolymerClassContainer.printAllIncludedResidues ( _parameterReader->physicsContainer.residueStretchVector ); //includeAllNonBondAtomsInResidueVector);

    ////////////// Reading included water droplet atoms //////////////
    cout<<__FILE__<<":"<<__LINE__<<" About to read water droplet atoms, physicsWhereYouWantIt up to now has included : "<< _dumm.getNumIncludedAtoms () <<" atoms. "<<endl;
    _parameterReader->waterDropletContainer.includeAllAtoms(_dumm);
    // this may be redundant .. should just do it once below:
    cout<<__FILE__<<":"<<__LINE__<<" system.realizeTopology() "<<endl;
    _state = _system.realizeTopology();

    ////////////// Reading included monoAtoms (e.g. ions) //////////////
    cout<<__FILE__<<":"<<__LINE__<<" About to read monoAtoms, physicsWhereYouWantIt up to now has included : "<< _dumm.getNumIncludedAtoms () <<" atoms. "<<endl;
    _parameterReader->myMonoAtomsContainer.includeAllAtoms(_dumm);
    cout<<__FILE__<<":"<<__LINE__<<" system.realizeTopology() "<<endl;
    _state = _system.realizeTopology();
    ////////////// Reading moleculeClassContainer //////////////
    cout<<__FILE__<<":"<<__LINE__<<" About to read moleculeClassContainer, physicsWhereYouWantIt up to now has included : "<< _dumm.getNumIncludedAtoms () <<" atoms. "<<endl;
    _parameterReader->moleculeClassContainer.includeAllAtoms(_dumm);
    cout<<__FILE__<<":"<<__LINE__<<" system.realizeTopology() "<<endl;
    _state = _system.realizeTopology();

    cout<<__FILE__<<":"<<__LINE__<<" read directly specified residues + water droplet atoms, physicsWhereYouWantIt now has included : "<< _dumm.getNumIncludedAtoms () <<" atoms. "<<endl;
    //////////////////////////////////////////////////////////////////
    //DuMM::IncludedAtom 
    if (_parameterReader->verbose){
    cout<<__FILE__<<":"<<__LINE__<<" Dumping "<<_dumm.getNumIncludedAtoms()<<" charged atom types: "<<endl;
    //for (DuMM::NonbondAtomIndex  nbax (0); nbax<_dumm.getNumNonbondAtoms (); nbax++)
    for (DuMM::IncludedAtomIndex iax1 = DuMM::IncludedAtomIndex(0); iax1<_dumm.getNumIncludedAtoms (); iax1++)
    {
                   
        //DuMM::IncludedAtomIndex iax1 = _dumm.getIncludedAtomIndexOfNonbondAtom(nbax);
        cout<<__FILE__<<":"<<__LINE__<<" partial charge:  "<<_dumm.getPartialCharge(iax1)<<" ";//<<endl;
        cout<<__FILE__<<":"<<__LINE__<<" atom type name:  "<<_dumm.getChargedAtomName(iax1) <<endl;
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
            cout<<__FILE__<<":"<<__LINE__<<" : Created Nose-Hoover thermostat."<<endl;
        }
        else if ((_parameterReader->thermostatType).compare("VelocityRescaling") == 0 ) { 
            VelocityRescalingThermostat * myVelocityRescalingThermostat = new VelocityRescalingThermostat(_system,_parameterReader->temperature, _parameterReader->velocityRescalingInterval);
            _system.updDefaultSubsystem().addEventHandler(myVelocityRescalingThermostat);
            ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" : Created velocity rescaling thermostat."<<endl;
        } else {
          ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" : Unrecognized thermostatType = "<<_parameterReader->thermostatType<<endl; 
          ErrorManager::instance.treatError(); 
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
        cout<<__FILE__<<" : "<<__LINE__<<" _parameterReader->setForceAndStericScrubber = "<<_parameterReader->setForceAndStericScrubber<<", now activating scrubber." <<endl;
       
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
        else
            {cout<<"[Repel.h] You must specify an integratorType of one of the supported types!!!"<<endl;assert(0);}

    if (_parameterReader->useFixedStepSize)
        _study->setFixedStepSize(_parameterReader->integratorStepSize);

    _study->setAccuracy(_parameterReader->integratorAccuracy);
    _study->setConstraintTolerance(_parameterReader->constraintTolerance);//_parameterReader->integratorAccuracy);

    ////////////////////////////////////////////

    // have to realize topology since we've added periodic event handlers and reporters and other system objects.
    cout<<__FILE__<<":"<<__LINE__<<" system.realizeTopology() "<<endl;
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
    cout<<"[Repel.h:ConstrainedDynamics] Starting dynamics now."<<endl;
    //cout<<"[Repel.h:ConstrainedDynamics] _study.getPredictedNextStepSize() " << (*_study).getPredictedNextStepSize()<<endl;
    cout<<__FILE__<<":"<<__LINE__<<" size of _matter subsystem Q vector or number of degrees of freedom: "<<(_matter.getQ(_state)).size()<<endl;
    if (_parameterReader->minimize) {
        if (_parameterReader->basePairContainer.numBasePairs() >0) {
            ErrorManager::instance <<__FILE__<<":"<<__LINE__<<"Error! If you want to minimize, you can't have any baseInteraction\'s!"<<endl;
            ErrorManager::instance.treatError();
        } 
        cout<<__FILE__<<":"<<__LINE__<<"Starting minimizer ..."<<endl;
        LocalEnergyMinimizer::minimizeEnergy(_system, _state, .001);
    } 
    _ts->initialize(_state);
    cout<<__FILE__<<":"<<__LINE__<<" Force subsystem topology has been realized: "<<_forces.subsystemTopologyHasBeenRealized()    <<endl;

    cout<<__FILE__<<":"<<__LINE__<<" Matter subsystem mass = "<<_matter.calcSystemMass(_state)<<std::endl;
    
    _nextFrame = 1;
    _previousTimeStep = 0.0;

}

void ConstrainedDynamics::postDynamics(){
    ////////////////////////////////////////////
    /// Post-dynamics tasks                  ///
    ////////////////////////////////////////////

    _output.close();

    cout<<"[Repel.h:ConstrainedDynamics] _study.getPreviousStepSizeTaken() " << (*_study).getPreviousStepSizeTaken()<<endl;
    cout<<"[Repel.h:ConstrainedDynamics] (*_study).getPredictedNextStepSize() " <<(*_study).getPredictedNextStepSize()<<endl;
    cout<<"[Repel.h:ConstrainedDynamics] (*_study).getActualInitialStepSizeTaken() " <<(*_study).getActualInitialStepSizeTaken()<<endl;
    cout<<"[Repel.h:ConstrainedDynamics] (*_study).getNumStepsTaken() " << (*_study).getNumStepsTaken()<<endl;
    cout<<"[Repel.h:ConstrainedDynamics] (*_study).getNumStepsAttempted() " <<(*_study).getNumStepsAttempted()<<endl;
    _state = _ts->getState();
    if (_parameterReader->writeLastFrameFile)
    {
        if ( _parameterReader->useCIFFileFormat )
        {
#ifdef CPP4_MAPS_USAGE
            //======================================== Initialise variables
            mmdb::Manager *mmdb2Manager               = new mmdb::Manager ( );
            mmdb::Model *mmdb2Model                   = new mmdb::Model ( mmdb2Manager, 1 );
            mmdb::io::File cifOutputFile;
            mmdb::mmcif::Data cifOutData;
            int compoundNumber                        = 1;
            int requiredPrecision                     = 19;
            
            //======================================== Create MMDB2 Biomolecule
            mmdb2Manager->MakeBiomolecule             ( 1, 1 );
            
            //======================================== Open file for writing
            cifOutputFile.assign                      ( _parameterReader->lastFrameFileName.c_str ( ) );
            if ( !cifOutputFile.rewrite ( ) )
            {
                ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" Failed to open the file "<< _parameterReader->lastFrameFileName.c_str ( ) << " for writing." << std::endl;
                ErrorManager::instance.treatError     ( );
            }
            
            //======================================== Set mmCIF file head
            std::string strName                       = _parameterReader->lastFrameFileName;
            istringstream iss                         ( strName );
            std::vector<std::string> tokens; std::string token;
            while ( std::getline ( iss, token, '.') ) { if ( !token.empty() ) { tokens.push_back ( token ); } }
            if ( tokens.size() > 2 ) { strName = std::string ( tokens.at(tokens.size()-3) + tokens.at(tokens.size()-2) ); }
            else { strName = "LASTX"; }
            strName.erase                             ( std::remove ( strName.begin(), strName.end(), '/'), strName.end() );
            std::transform                            ( strName.begin(), strName.end(), strName.begin(), ::toupper);
            cifOutData.PutDataName                    ( strName.c_str() );
            
            //======================================== For each compound
            for (SimTK::CompoundSystem::CompoundIndex c(0); c < _system.getNumCompounds(); ++c)
            {
		std::cout <<__FILE__<<":"<<__LINE__<<" c = "<<c<< " compoundNumber = "<<compoundNumber <<std::endl;     
                //==================================== Build the MMDB2 structure
                (_system.getCompound(c)).buildCif     ( _state, mmdb2Model, Transform( Vec3 ( 0 ) ) );
                
                //==================================== Add the model to the data
                 for ( int noModel = 0; noModel < mmdb2Model->GetNumberOfChains(); noModel++ )
                 {
                     if ( mmdb2Model->GetChain(noModel) )
                     {
                         for ( int noRes = 0; noRes < mmdb2Model->GetChain(noModel)->GetNumberOfResidues(); noRes++ )
                         {
                             if ( mmdb2Model->GetChain(noModel)->GetResidue(noRes) )
                             {
                                 for ( int noAt = 0; noAt < mmdb2Model->GetChain(noModel)->GetResidue(noRes)->GetNumberOfAtoms(); noAt++ )
                                 {
                                     if ( mmdb2Model->GetChain(noModel)->GetResidue(noRes)->GetAtom(noAt) )
                                     {
                                         //==================== Initialise variables
                                         mmdb::mmcif::PLoop Loop;
                                         mmdb::AtomName AtName;
                                         mmdb::Element el;
                                         char N[10];
                                         int i,j,RC;
                                         mmdb::PChain chain  = mmdb2Model->GetChain(noModel)->GetResidue(noRes)->chain;
                                         mmdb::PModel model  = mmdb::PModel ( mmdb2Model->GetChain(noModel)->GetModel() );
                                         
                                         //================ Initialise the Loop object
                                         RC           = cifOutData.AddLoop ( mmdb::CIFCAT_ATOM_SITE, Loop );
                                         if ( RC != mmdb::mmcif::CIFRC_Ok )
                                         {
                                             //============ The category was (re)created, provide tags
                                             Loop->AddLoopTag ( mmdb::CIFTAG_GROUP_PDB          ); // ATOM, TER etc.
                                             Loop->AddLoopTag ( mmdb::CIFTAG_ID                 ); // serial number

                                             Loop->AddLoopTag ( mmdb::CIFTAG_TYPE_SYMBOL        ); // element symbol
                                             Loop->AddLoopTag ( mmdb::CIFTAG_LABEL_ATOM_ID      ); // atom name
                                             Loop->AddLoopTag ( mmdb::CIFTAG_LABEL_ALT_ID       ); // alt location
                                             Loop->AddLoopTag ( mmdb::CIFTAG_LABEL_COMP_ID      ); // residue name
                                             Loop->AddLoopTag ( mmdb::CIFTAG_LABEL_ASYM_ID      ); // chain ID
                                             Loop->AddLoopTag ( mmdb::CIFTAG_LABEL_ENTITY_ID    ); // entity ID
                                             Loop->AddLoopTag ( mmdb::CIFTAG_LABEL_SEQ_ID       ); // res seq number
                                             Loop->AddLoopTag ( mmdb::CIFTAG_PDBX_PDB_INS_CODE  ); // insertion code
                                             Loop->AddLoopTag ( mmdb::CIFTAG_SEGMENT_ID         ); // segment ID

                                             Loop->AddLoopTag ( mmdb::CIFTAG_CARTN_X            ); // x-coordinate
                                             Loop->AddLoopTag ( mmdb::CIFTAG_CARTN_Y            ); // y-coordinate
                                             Loop->AddLoopTag ( mmdb::CIFTAG_CARTN_Z            ); // z-coordinate
                                             Loop->AddLoopTag ( mmdb::CIFTAG_OCCUPANCY          ); // occupancy
                                             Loop->AddLoopTag ( mmdb::CIFTAG_B_ISO_OR_EQUIV     ); // temp factor

                                             Loop->AddLoopTag ( mmdb::CIFTAG_CARTN_X_ESD        ); // x-sigma
                                             Loop->AddLoopTag ( mmdb::CIFTAG_CARTN_Y_ESD        ); // y-sigma
                                             Loop->AddLoopTag ( mmdb::CIFTAG_CARTN_Z_ESD        ); // z-sigma
                                             Loop->AddLoopTag ( mmdb::CIFTAG_OCCUPANCY_ESD      ); // occupancy-sigma
                                             Loop->AddLoopTag ( mmdb::CIFTAG_B_ISO_OR_EQUIV_ESD ); // t-factor-sigma

                                             Loop->AddLoopTag ( mmdb::CIFTAG_PDBX_FORMAL_CHARGE ); // charge on atom

                                             Loop->AddLoopTag ( mmdb::CIFTAG_AUTH_SEQ_ID        ); // res seq number
                                             Loop->AddLoopTag ( mmdb::CIFTAG_AUTH_COMP_ID       ); // residue name
                                             Loop->AddLoopTag ( mmdb::CIFTAG_AUTH_ASYM_ID       ); // chain id
                                             Loop->AddLoopTag ( mmdb::CIFTAG_AUTH_ATOM_ID       ); // atom name

                                             Loop->AddLoopTag ( mmdb::CIFTAG_PDBX_PDB_MODEL_NUM ); // model number
                                         }
                                         
                                         //================ Is this a normal atom record?
                                         if ( mmdb2Model->GetChain(noModel)->GetResidue(noRes)->GetAtom(noAt)->WhatIsSet & ( mmdb::ASET_Coordinates | mmdb::ASET_CoordSigma))
                                         {
                                             //============ Yes!
                                             
                                             // group_PDB field
                                             if ( mmdb2Model->GetChain(noModel)->GetResidue(noRes)->GetAtom(noAt)->Het ) { Loop->AddString ( pstr ( "HETATM" ) ); }
                                             else                                                                        { Loop->AddString ( pstr( "ATOM" ) ); }
                                         
                                             // id field
                                             if ( mmdb2Model->GetChain(noModel)->GetResidue(noRes)->GetAtom(noAt)->serNum > 0 ) { Loop->AddInteger ( mmdb2Model->GetChain(noModel)->GetResidue(noRes)->GetAtom(noAt)->serNum ); }
                                             else                                                                               { Loop->AddInteger ( mmdb2Model->GetChain(noModel)->GetResidue(noRes)->GetAtom(noAt)->GetIndex() ); }
                                             
                                             if ( mmdb2Model->GetChain(noModel)->GetResidue(noRes)->GetAtom(noAt)->WhatIsSet & mmdb::ASET_Coordinates)
                                             {
                                                 // type_symbol field
                                                 mmdb::strcpy_css ( el, mmdb2Model->GetChain(noModel)->GetResidue(noRes)->GetAtom(noAt)->element );
                                                 Loop->AddString ( el, true );
                                                 
                                                 // label_atom_id field
                                                 mmdb::strcpy_css ( AtName, mmdb2Model->GetChain(noModel)->GetResidue(noRes)->GetAtom(noAt)->label_atom_id );
                                                 Loop->AddString ( AtName );
                                                 
                                                 // label_alt_id field
                                                 Loop->AddString  ( mmdb2Model->GetChain(noModel)->GetResidue(noRes)->GetAtom(noAt)->altLoc, true );

                                                 // label_comp_id field
                                                 Loop->AddString ( mmdb2Model->GetChain(noModel)->GetResidue(noRes)->label_comp_id );
                                                 
                                                 // label_asym_id field
                                                 Loop->AddString ( mmdb2Model->GetChain(noModel)->GetResidue(noRes)->label_asym_id );
                                                 
                                                 // label_entity_id field
                                                 if ( mmdb2Model->GetChain(noModel)->GetResidue(noRes)->label_entity_id > 0 ) { Loop->AddInteger ( mmdb2Model->GetChain(noModel)->GetResidue(noRes)->label_entity_id ); }
                                                 else                                                                         { Loop->AddNoData  ( mmdb::mmcif::CIF_NODATA_DOT ); }
                                                 
                                                 // label_seq_id field
                                                 if ( mmdb2Model->GetChain(noModel)->GetResidue(noRes)->label_seq_id > mmdb::MinInt4 ) { Loop->AddInteger ( mmdb2Model->GetChain(noModel)->GetResidue(noRes)->label_seq_id ); }
                                                 else                                                                                  { Loop->AddNoData ( mmdb::mmcif::CIF_NODATA_DOT ); }
                                                 
                                                 // pdbx_PDB_ins_code field
                                                 Loop->AddString ( mmdb2Model->GetChain(noModel)->GetResidue(noRes)->insCode, true );

                                                 // segment_id field
                                                 Loop->AddString ( mmdb2Model->GetChain(noModel)->GetResidue(noRes)->GetAtom(noAt)->segID, true );
                                                 
                                                 // Cartn_x, Cartn_y, Cartn_z fields
                                                 Loop->AddReal ( mmdb2Model->GetChain(noModel)->GetResidue(noRes)->GetAtom(noAt)->x, requiredPrecision );
                                                 Loop->AddReal ( mmdb2Model->GetChain(noModel)->GetResidue(noRes)->GetAtom(noAt)->y, requiredPrecision );
                                                 Loop->AddReal ( mmdb2Model->GetChain(noModel)->GetResidue(noRes)->GetAtom(noAt)->z, requiredPrecision );
                                                 
                                                 // occupancy field
                                                 if ( mmdb2Model->GetChain(noModel)->GetResidue(noRes)->GetAtom(noAt)->WhatIsSet & mmdb::ASET_Occupancy ) { Loop->AddReal ( mmdb2Model->GetChain(noModel)->GetResidue(noRes)->GetAtom(noAt)->occupancy, requiredPrecision ); }
                                                 else { Loop->AddNoData ( mmdb::mmcif::CIF_NODATA_QUESTION ); }
                                                 
                                                 // B_iso_or_equiv field
                                                 if ( mmdb2Model->GetChain(noModel)->GetResidue(noRes)->GetAtom(noAt)->WhatIsSet & mmdb::ASET_tempFactor ) { Loop->AddReal ( mmdb2Model->GetChain(noModel)->GetResidue(noRes)->GetAtom(noAt)->tempFactor, requiredPrecision ); }
                                                 else { Loop->AddNoData ( mmdb::mmcif::CIF_NODATA_QUESTION ); }

                                                 // cartn_x_esd, cartn_y_esd, cartn_z_esd fields
                                                 if ( mmdb2Model->GetChain(noModel)->GetResidue(noRes)->GetAtom(noAt)->WhatIsSet & mmdb::ASET_CoordSigma )
                                                 {
                                                   Loop->AddReal ( mmdb2Model->GetChain(noModel)->GetResidue(noRes)->GetAtom(noAt)->sigX, requiredPrecision );
                                                   Loop->AddReal ( mmdb2Model->GetChain(noModel)->GetResidue(noRes)->GetAtom(noAt)->sigY, requiredPrecision );
                                                   Loop->AddReal ( mmdb2Model->GetChain(noModel)->GetResidue(noRes)->GetAtom(noAt)->sigZ, requiredPrecision );
                                                 } else
                                                 {
                                                   Loop->AddNoData ( mmdb::mmcif::CIF_NODATA_QUESTION );
                                                   Loop->AddNoData ( mmdb::mmcif::CIF_NODATA_QUESTION );
                                                   Loop->AddNoData ( mmdb::mmcif::CIF_NODATA_QUESTION );
                                                 }
                                                 
                                                 // occupancy_esd field
                                                 if ( ( mmdb2Model->GetChain(noModel)->GetResidue(noRes)->GetAtom(noAt)->WhatIsSet & mmdb::ASET_OccSigma) &&
                                                      ( mmdb2Model->GetChain(noModel)->GetResidue(noRes)->GetAtom(noAt)->WhatIsSet & mmdb::ASET_Occupancy) )
                                                 {
                                                       Loop->AddReal ( mmdb2Model->GetChain(noModel)->GetResidue(noRes)->GetAtom(noAt)->sigOcc, requiredPrecision );
                                                 }
                                                 else { Loop->AddNoData ( mmdb::mmcif::CIF_NODATA_QUESTION ); }
                                                 
                                               // B_iso_or_equiv_esd field
                                               if ( ( mmdb2Model->GetChain(noModel)->GetResidue(noRes)->GetAtom(noAt)->WhatIsSet & mmdb::ASET_tFacSigma) &&
                                                    ( mmdb2Model->GetChain(noModel)->GetResidue(noRes)->GetAtom(noAt)->WhatIsSet & mmdb::ASET_tempFactor) )
                                               {
                                                     Loop->AddReal ( mmdb2Model->GetChain(noModel)->GetResidue(noRes)->GetAtom(noAt)->sigTemp, requiredPrecision );
                                               }
                                               else { Loop->AddNoData ( mmdb::mmcif::CIF_NODATA_QUESTION ); }
                                             }
                                             else
                                             {
                                               for ( int iter = 0; iter < 18; iter++ )
                                               {
                                                 Loop->AddNoData ( mmdb::mmcif::CIF_NODATA_QUESTION );
                                               }
                                             }
                                             
                                             // pdbx_formal_charge field
                                             if ( mmdb2Model->GetChain(noModel)->GetResidue(noRes)->GetAtom(noAt)->WhatIsSet & mmdb::ASET_Charge)
                                             {
                                               sprintf ( N, "%+2i", mmdb::mround ( mmdb2Model->GetChain(noModel)->GetResidue(noRes)->GetAtom(noAt)->charge ) );
                                               Loop->AddString ( N, true );
                                             }
                                             else { Loop->AddNoData ( mmdb::mmcif::CIF_NODATA_QUESTION ); }

                                             // auth_seq_id field
                                             if ( mmdb2Model->GetChain(noModel)->GetResidue(noRes)->seqNum > mmdb::MinInt4 ) { Loop->AddInteger ( mmdb2Model->GetChain(noModel)->GetResidue(noRes)->seqNum ); }
                                             else                                                                            { Loop->AddNoData ( mmdb::mmcif::CIF_NODATA_DOT ); }
                                             
                                             // auth_comp_id field
                                             Loop->AddString ( mmdb2Model->GetChain(noModel)->GetResidue(noRes)->name );

                                             // auth_asym_id field
                                             if ( mmdb2Model->GetChain(noModel) ) { Loop->AddString ( mmdb2Model->GetChain(noModel)->GetChainID(), true ); }
                                             else { Loop->AddNoData ( mmdb::mmcif::CIF_NODATA_DOT ); }
                                             
                                             // auth_atom_id field
                                             mmdb::strcpy_css ( AtName, mmdb2Model->GetChain(noModel)->GetResidue(noRes)->GetAtom(noAt)->name );
                                             Loop->AddString  ( AtName );
                                             
                                             // pdbx_PDB_model_num field
                                             if ( mmdb2Model->GetSerNum() > 0) { Loop->AddInteger ( mmdb2Model->GetSerNum() ); }
                                             else                              { Loop->AddNoData ( mmdb::mmcif::CIF_NODATA_QUESTION ); }
                                         }
                                     }
                                 }
                             }
                         }
                     }
                 }
                
                //==================================== Write out the CIF file. This is basically a copy of the WriteMMCIF function, but since this function is hardcoded to print the data_ line as a second line and this is not compatible with the current mmCIF formatting, I had to copy and make the single change...
                cifOutputFile.Write                   ( pstr("data_") );
                cifOutputFile.WriteLine               ( cifOutData.GetDataName() );
                
                for ( int i = 0; i < cifOutData.GetNumberOfCategories(); i++ )
                {
                    mmdb::mmcif::Category* cat        = cifOutData.GetCategory(i);
                    if ( cat )
                    {
                        cat->WriteMMCIF               ( cifOutputFile );
                    }
                }
                
                //==================================== Write the entity poly seq loop
                (_system.getCompound(c)).writeEntityPolySeqLoop ( _state, &cifOutputFile, compoundNumber );
		// SCF decided this was being prematurely closed, moved down a few lines, out of the loop
                // cifOutputFile.shut                    ( );
                
                //==================================== Prepare for next compound
                ++compoundNumber;
            }
            cifOutputFile.shut                    ( );
#else
            ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" Error! Requested mmCIF file output, but did not compile with the MMDB2 library. Cannot proceed, if you want to use mmCIF files, please re-compile with the MMDB2 library option allowed." <<endl;
            ErrorManager::instance.treatError         ( );
#endif
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

    cout<<__FILE__<<":"<<__LINE__<<" Stage "<<_parameterReader->currentStage<<" completed successfully! "<<endl;
}

void ConstrainedDynamics::runAllSteps() {
    if(getRemainingFramesNum() <= 0){
        ErrorManager::instance << __FILE__ << " " << __FUNCTION__ << ": Cannot run more steps" << endl;
        ErrorManager::instance.treatError();
    }

    // Normal way to compute dynamics
    //_parameterReader->myBiopolymerClassContainer.printBiopolymerInfo();
    // _ts->stepTo(_parameterReader->numReportingIntervals*_parameterReader->reportingInterval);
    cout<<__FILE__<<":"<<__LINE__<<" About to run time integrator for "<<_parameterReader->numReportingIntervals<< " reporting intervals, each of duration "<<_parameterReader->reportingInterval<<" ps "<<endl;

    // New way, to allow a step by step control later
        
    for(_nextFrame; _nextFrame <= _parameterReader->numReportingIntervals; _nextFrame++){
        //cout<<__FILE__<<":"<<__LINE__<<std::endl; 
        _ts->stepTo(_nextFrame*_parameterReader->reportingInterval);
        _state = _ts->getState();
        //cout<<__FILE__<<":"<<__LINE__<<std::endl; 
        if(_parameterReader->detectConvergence && _parameterReader->converged){
            //cout<<__FILE__<<":"<<__LINE__<<std::endl; 
            _nextFrame = _parameterReader->numReportingIntervals;
        }
    }

}

unsigned int ConstrainedDynamics::runOneStep(){
    if(getRemainingFramesNum() <= 0){
        ErrorManager::instance << __FILE__ << " " << __FUNCTION__ << ": Cannot run more steps" << endl;
        ErrorManager::instance.treatError();
    }

    double nextStep = _previousTimeStep + _parameterReader->reportingInterval;
    _ts->stepTo(nextStep);
    _state = _ts->getState();

    _previousTimeStep = nextStep;
    _nextFrame ++;

    return getRemainingFramesNum();
}

void ConstrainedDynamics::initializeBodies(){
    //cout<<__FILE__<<":"<<__LINE__<<endl;
    initializeMoleculesAndBonds();
    cout<<__FILE__<<":"<<__LINE__<<endl;
    #ifdef USE_OPENMM
    setInterfaceMobilizers();
    #endif
    setMobilizers();
    createMultibodyTree();
}

void ConstrainedDynamics::initializeDynamics(){
    // This is where we call mobilizerContainer.addMobilizerStretchesToVector
    // This is where initializeAtomInfoVectors is called (at least when density map fitting is active): 
    // but note that initializeMoleculesAndBonds() precedes createState()
    cout<<__FILE__<<":"<<__LINE__<<endl;
    //_parameterReader->myBiopolymerClassContainer.printAtomInfoVector(); //   at this point ..
    initializeCustomForcesConstraints();
    cout<<__FILE__<<":"<<__LINE__<<endl;
    createEventHandlers();    
    cout<<__FILE__<<":"<<__LINE__<<endl;
    initializeIntegrator();    
    cout<<__FILE__<<":"<<__LINE__<<" Net charge over "<<_dumm.getNumIncludedAtoms ()<<" atoms included in physics zone = "<< _dumm.getTotalIncludedCharge()<<endl;
}

void ConstrainedDynamics::writeMMBPDB(std::ofstream & filestream){
    CompoundSystem system;
    State state;
    SimbodyMatterSubsystem matter(system);
    DuMMForceFieldSubsystem dumm(system);
    _parameterReader->configureDumm(dumm);
    cout<<__FILE__<<":"<<__LINE__<<endl;
    initializeBiopolymersAndCustomMolecules(system);
    initializeMoleculesAndBonds(system, dumm, matter);
    //cout<<__FILE__<<":"<<__LINE__<<endl;
    //cout<<__FILE__<<":"<<__LINE__<<endl;
    // setInterfaceMobilizers(system, matter, state);
    // setMobilizers();
    createMultibodyTree(system, state);
    _parameterReader->myBiopolymerClassContainer.writePdb(state, system, filestream);
}

void ConstrainedDynamics::runDynamics() {
    cout<<__FILE__<<":"<<__LINE__<<endl;
    if (initializeBiopolymersAndCustomMolecules()){
        cout<<__FILE__<<":"<<__LINE__<<" Warning: Returned an error from initializeBiopolymersAndCustomMolecules"<<std::endl;
        //exit(1) ;
    }
    cout<<__FILE__<<":"<<__LINE__<<endl;
    initializeBodies();
    //cout<<__FILE__<<":"<<__LINE__<<endl;
     
    //cout<<__FILE__<<":"<<__LINE__<<endl;
    //_parameterReader->myBiopolymerClassContainer.printAtomInfoVector(); //  Looks fine at this point ..

    initializeDynamics();

    runAllSteps();

    postDynamics();
};

