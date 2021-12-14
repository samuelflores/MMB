#ifndef Repel_H_
#define Repel_H_

/* -------------------------------------------------------------------------- *
 *                           MMB (MacroMoleculeBuilder)                       *
 * -------------------------------------------------------------------------- *
 *                                                                            *
 * Copyright (c) 2011-12 by the Author.                                       *
 * Author: Samuel Flores                                                      *
 *                                                                            *
 * See RNABuilder.cpp for the copyright and usage agreement.                  *
 * -------------------------------------------------------------------------- */

//const int maxChiWelds = 1000;
//const int maxResidues = 3000;    // was maxChars

#include <iostream>
#include <vector>
#include <stdlib.h>// for MAX_PATH
#include "SimTKmolmodel.h"
#include "PeriodicScrubber.h"  
#include "BiopolymerClass.h"
#include "ParameterReader.h"
#include "PeriodicPdbAndEnergyWriter.h"
#include "BiopolymerClassTwoTransformForces.h"
#ifdef NTC_ENABLED
#include "NtCForces.h"
#endif // NTC_ENABLED
#include "TetherForce.h"
#include "Sterics.h" 
#include "AddNASTForces.h"
#include "WadleyKeatingDuartePyleTorsionForce.h"
#include "SetSingleBondMobility.h" 
#include "DensityForce.h"
#include "ConstraintContainer.h"

//using namespace SimTK;
//using namespace std;


/**
 * /brief This method sets the BondMobility for all bonds within each residue in a certain stretch of residues of any biopolymer   chain
 *
 * added by scf
 */
MMB_EXPORT void setBiopolymerResidueBondMobility (Biopolymer & myChain , BondMobility::Mobility  mobility, ResidueInfo::Index startResidue, ResidueInfo::Index endResidue);


/**
 * /brief This method sets the BondMobility for all bonds in a certain stretch of residues of a Protein        chain
 *
 * added by scf
 */
MMB_EXPORT void setProteinBondMobility (Biopolymer & myProteinChain , BondMobility::Mobility  mobility, ResidueInfo::Index startResidue, ResidueInfo::Index endResidue);



MMB_EXPORT double PointToPlaneDistance (Vec3 Point1, Vec3 Normal1, Vec3 Point2) ;


/**
 * /brief This utility is used to get the name of an atom involved in hydrogen bonding. It was written to help with the problem of pulling together two halves of an A-form double helix. 
 *
 * By convention the first atom of the first residue is a bond donor "H", the second atom of the first residue is an accepter "O".
 * the opposite convention holds for the second residue.
 *
 * /param PdbResidueName is the type of base (AUCG).  Based on this parameter we will determine the name of the appropriate hydrogen-bonding atom.
 * 
 */

// this polymorphism soon to be obsolete

// MMB_EXPORT String ReturnHBondingAtomName(String myPdbResidueName, int firstOrSecondResidue, int firstOrSecondBond );

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
/*
class  SetPolynucleotideChiBondMobility  { public:

    SetPolynucleotideChiBondMobility (
        Biopolymer & myChain,
	int startResidueNumber,
	int endResidueNumber,
	LeontisWesthofClass  & myLeontisWesthofClass ,    
        SimbodyMatterSubsystem & matter,
        State& state,
        Constraint myWeld1[maxChiWelds],
        Constraint myWeld2[maxChiWelds],
        Constraint myWeld3[maxChiWelds]

	) ;  
    SetPolynucleotideChiBondMobility (
        Biopolymer & myChain,
	int startResidueNumber,
	int endResidueNumber,
	LeontisWesthofClass  & myLeontisWesthofClass ,    
	BondMobility::Mobility mobility = BondMobility::Free   
	)   ;
};       
*/
    //static bool compareUpper( const String& param, const char* symbol );


// MMB_EXPORT void setLeontisWesthofBondRowIndex (ParameterReader & myParameterReader,LeontisWesthofClass myLeontisWesthofClass, vector <Biopolymer>  & myMolecule);


class MMB_EXPORT ConstrainedDynamics : public Compound {
public:
    // int setDefaults();  
    void constraintsAndRestraints  (ParameterReader & ,BiopolymerClassContainer & ,GeneralForceSubsystem &  ,SimbodyMatterSubsystem & , State & , CompoundSystem &  );
    // void setBondMobilities(ParameterReader & myParameterReader,vector<Biopolymer> & myMolecule );
    // void addContacts(ParameterReader & myParameterReader,vector<Biopolymer> & myMolecule, GeneralContactSubsystem & contacts, GeneralForceSubsystem & forces, SimbodyMatterSubsystem & matter,CompoundSystem & system );
     
    void runDynamics();

    /**
    * Constructor
    *
    * \param a parameter reader
    */
    ConstrainedDynamics(ParameterReader * myParameterReader);

    /**
    * Destructor
    */
    ~ConstrainedDynamics();

    // Getters
    /**
    * \return the compound system
    */
    CompoundSystem & getCompoundSystem() { return _system; }

    /**
    * \return next frame number
    */
    unsigned int getNextFrameNum() { return _nextFrame; }

    /**
    * \return remaining frames number
    */
    unsigned int getRemainingFramesNum() { return _parameterReader->numReportingIntervals - (_nextFrame - 1); }

    /**
    * \return current state
    */
    State & getCurrentState() { return _state; }



    // Initializers

    void initializeDumm();

    /**
    * Initialize Biopolymers starting position
    */
    int  initializeBiopolymersAndCustomMolecules();
    int  initializeBiopolymersAndCustomMolecules(CompoundSystem & system); 

    /**
    * Initialize one Biopolymer identified by its chain
    * If no input file is provided, initialize with default positions
    */
    void initializeBiopolymer(String chainID, String inputPDBFile="");

    /**
    * Initialize other molecules and Bonds mobilizers
    * Must be called after initializeBiopolymers()
    */
    void initializeMoleculesAndBonds();
    void initializeMoleculesAndBonds(CompoundSystem & system, DuMMForceFieldSubsystem & dumm, SimbodyMatterSubsystem & matter);
    void setInterfaceMobilizers();
    void setInterfaceMobilizers(CompoundSystem & system, SimbodyMatterSubsystem & matter, State & state);
    void setMobilizers();

    /**
    * Create Multibody Tree.
    * Must be called after initializedMoleculesAndBonds()
    */
    void createMultibodyTree();
    void createMultibodyTree(CompoundSystem & system, State & state);


    /**
    * Instantiate Custom and other forces and constraints.
    */
    void initializeCustomForcesConstraints();

    /** 
    * Create periodic event handlers and reporters.
    */
    void createEventHandlers();

    /**
    * Initialize integrator and time stepper.
    */
    void initializeIntegrator();


    /**
    * Initialize other molecules, bonds, 
    * the Mutltibody tree and create the state
    */
    void initializeBodies ();
    void forceAdjustmentsWithFinalMobilizers(); // Some adjustments are made here to forces depending on mobilizers,  the mobilizers in turn are set in initializeBodies().

    /**
    * Initialize everything
    */
    void initializeDynamics();


    // Dynamics operation
    /**
    * Post-dynamics operations
    * Close the trajectory file and delete pointers
    * If required, write the last frame in an independant file
    */
    void postDynamics();


    /**
    * Run dynamics from next frame to numReportingIntervals
    */
    void runAllSteps();

    /**
    * Execute one step of dynamics
    *
    * \return number of remaining frames
    */
    unsigned int runOneStep();

    /**
    * Write original coordinates to the specified file.
    * It uses a dummy CompoundSystem.
    * \param filestream Ouptput file
    */
    void writeMMBPDB(std::ofstream & filestream);

private:
    CompoundSystem          _system;
    SimbodyMatterSubsystem  _matter;
    GeneralForceSubsystem   _forces;
    GeneralContactSubsystem _contacts;
    DuMMForceFieldSubsystem _dumm;
    State                   _state;
    Integrator *            _study;
    TimeStepper*            _ts;

    ofstream                _output;

    ParameterReader *       _parameterReader;

    unsigned int            _nextFrame;  
    double                  _previousTimeStep;  
};   
#endif
