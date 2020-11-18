/* -------------------------------------------------------------------------- *
 *                           MMB (MacroMoleculeBuilder)                       *
 * -------------------------------------------------------------------------- *
 *                                                                            *
 * Copyright (c) 2011-12 by the Author.                                       *
 * Author: Alex Tek                                                           *
 *                                                                            *
 * See RNABuilder.cpp for the copyright and usage agreement.                  *
 * -------------------------------------------------------------------------- */
#include <fstream>
#include <ios>
#include <iostream>
#include <vector>
#include <cstring>
#include <sstream>
#include <string>
#include <utility>

// #include "SimTKmolmodel.h"
 //#include "SimTKsimbody_aux.h"

#include "ParameterReader.h"
#include "MMBLogger.h"
#include "Repel.h"
#include "ParameterReader_wrapper.h"

using namespace std;

extern "C" {

    /***********************************************************************************/
    /**
    * Description of a Biopolymer sequence.
    */
    typedef struct BioPolymerSeq{
        const char * chainID;           ///< chain ID
        const char * sequence;  ///< sequence string
        int firstResNum;        ///< First residue number
        int polyType;           ///< Type of polymer (0:RNA, 1:Protein, 2:DNA, 3:Unassigned)
        const char * pdbFileName; ///< Pdb file containing this sequence
        bool loadFromPdb;
        bool activePhysics;
    } BioPolymerSeq;


    /***********************************************************************************/
    /**
    * Description of a BaseInteraction
    */
    typedef struct BaseInteraction_wrapper{
        int mmbID;
        const char * FirstBPEdge;
        const char * SecondBPEdge;
        const char * FirstBPChain;
        const char * SecondBPChain;
        int FirstBPResidue;
        int SecondBPResidue;
        const char * OrientationBP;
        ////String BasePairIsTwoTransformForce;
        ////int BasePairPriority;
        ////int BasePairTemporary;
        //Rotation rotationCorrection1;
        //Rotation rotationCorrection2;
        //Vec3 translationCorrection1;
        //Vec3 translationCorrection2;
        bool basePairSatisfied;
        int leontisWesthofBondRowIndex;
    }BaseInteraction_wrapper;

    /**
    * Create a BaseInteraction_wrapper from a MMB's BaseInteraction
    */
    BaseInteraction_wrapper createBaseInteraction_wrapper(const BaseInteraction & _struct_){
        BaseInteraction_wrapper _wrap_;
        _wrap_.FirstBPEdge      = strdup(_struct_.FirstBPEdge.c_str());
        _wrap_.SecondBPEdge     = strdup(_struct_.SecondBPEdge.c_str());
        _wrap_.FirstBPChain     = strdup(_struct_.FirstBPChain.c_str());
        _wrap_.SecondBPChain    = strdup(_struct_.SecondBPChain.c_str());
        _wrap_.FirstBPResidue   = _struct_.FirstBPResidue.getResidueNumber();
        _wrap_.SecondBPResidue  = _struct_.SecondBPResidue.getResidueNumber();
        _wrap_.OrientationBP    = strdup(_struct_.OrientationBP.c_str());
    //  _wrap_.BasePairIsTwoTransformForce = _struct_.BasePairIsTwoTransformForce;
    //  _wrap_.BasePairPriority = _struct_.BasePairPriority;
    //  _wrap_.BasePairTemporary = _struct_.BasePairTemporary;
    //  _wrap_.rotationCorrection1 = _struct_.rotationCorrection1;
    //  _wrap_.rotationCorrection2 = _struct_.rotationCorrection2;
    //  _wrap_.translationCorrection1 = _struct_.translationCorrection1;
    //  _wrap_.translationCorrection2 = _struct_.translationCorrection2;
        _wrap_.basePairSatisfied = _struct_.basePairSatisfied;
        _wrap_.leontisWesthofBondRowIndex = _struct_.leontisWesthofBondRowIndex;
        return _wrap_; 
    }
    /***********************************************************************************/

    /***********************************************************************************/   
    /**
    * Description of a Mobilizer.
    */
    typedef struct MobilizerStretch_wrapper{
        int mmbID;              ///< Id in MMB's MobilizerContainer
        const char * chainID;           ///< chain ID
        const char * mobility;  ///< mobility type string
        int resStart;        ///< First residue number
        int resEnd;         ///< Last residue number
    } MobilizerStretch_wrapper;

    /***********************************************************************************/
    /**
    * Description of a MobilizerWithin.
    */
    typedef struct MobilizerWithin_wrapper{
        int mmbID;              ///< Id in MMB's MobilizerContainer
        const char * chainID;           ///< chain ID
        const char * mobility;  ///< mobility type string
        int resID;        ///< residue number
        double radius;         ///< Radius
    } MobilizerWithin_wrapper;

    /***********************************************************************************/   
    /**
    * Description of a Contact.
    */
    typedef struct ContactStretch_wrapper{
        int mmbID;              ///< Id in MMB's ContactContainer
        const char * chainID;           ///< chain ID
        const char * ContactScheme;  ///< mobility type string
        int resStart;        ///< First residue number
        int resEnd;         ///< Last residue number
    } ContactStretch_wrapper;

    /**
    * Create a ContactStretch_wrapper from a MMB's ContactStretch
    */
    ContactStretch_wrapper createContactStretch_wrapper(ContactStretch & _class_){
        ContactStretch_wrapper _wrap_;
        _wrap_.ContactScheme = strdup(_class_.ContactScheme.c_str());
        _wrap_.chainID = strdup(_class_.getChain().c_str());
        _wrap_.resStart = _class_.getStartResidue().getResidueNumber();
        _wrap_.resEnd = _class_.getEndResidue().getResidueNumber();
        return _wrap_; 
    }

    /***********************************************************************************/
    /**
    * Description of a ContactWithin.
    */
    typedef struct ContactWithin_wrapper{
        int mmbID;              ///< Id in MMB's ContactContainer
        const char * chainID;           ///< chain ID
        const char * ContactScheme;  ///< mobility type string
        int resID;        ///< residue number
        double radius;         ///< Radius
    } ContactWithin_wrapper;

    /**
    * Create a ContactWithin_wrapper from a MMB's ContactWithin
    */
    ContactWithin_wrapper createContactWithin_wrapper(ContactWithin & _struct_){
        ContactWithin_wrapper _wrap_;
        _wrap_.ContactScheme = strdup(_struct_.ContactScheme.c_str());
        _wrap_.chainID = strdup(_struct_.Chain.c_str());
        _wrap_.resID = _struct_.Residue.getResidueNumber();
        _wrap_.radius = _struct_.Radius;
        return _wrap_; 
    }


    /***********************************************************************************/
    /**
    * Description of a Constraint.
    */
    typedef struct Constraint_wrapper{
        const char * chain1;           ///< chain ID 1
        const char * chain2;           ///< chain ID 2
        char * atom1;         ///< AtomName 1
        char * atom2;         ///< AtomName 2
        int mmbID;              ///< Id in MMB's ConstraintClassVector
        int res1;        ///< residue number 1
        int res2;        ///< residue number 2
    } Constraint_wrapper;

    /***********************************************************************************/ 

    typedef struct AllResiduesWithin_wrapper{
        int mmbID;
        const char * chain;
        int residue;
        double radius;
    }AllResiduesWithin_wrapper;

    /**
    * Create a AllResiduesWithin_wrapper from a MMB's AllResiduesWithin
    */
    AllResiduesWithin_wrapper createAllResiduesWithin_wrapper(AllResiduesWithin & _struct_){
        AllResiduesWithin_wrapper _wrap_;
        _wrap_.chain    = strdup(_struct_.getChain().c_str());
        _wrap_.residue  = _struct_.getResidue().getResidueNumber();
        _wrap_.radius   = _struct_.getRadius();
        return _wrap_; 
    }

    typedef struct IncludeAllNonBondAtomsInResidue_wrapper{
        int mmbID;
        const char * chain;
        int residue;
    }IncludeAllNonBondAtomsInResidue_wrapper;


    /**
    * Create a IncludeAllNonBondAtomsInResidue_wrapper from a MMB's IncludeAllNonBondAtomsInResidue
    */
    IncludeAllNonBondAtomsInResidue_wrapper createIncludeAllNonBondAtomsInResidue_wrapper(IncludeAllNonBondAtomsInResidue & _struct_){
        IncludeAllNonBondAtomsInResidue_wrapper _wrap_;
        _wrap_.chain    = strdup(_struct_.getChain().c_str());
        _wrap_.residue  = _struct_.getResidue().getResidueNumber();
        return _wrap_; 
    }


    /***********************************************************************************/
    typedef struct ThreadingStruct_wrapper{
        int mmbID;
        const char * chainID1;
        int residueStart1;
        int residueEnd1;
        const char * chainID2;
        int residueStart2;
        int residueEnd2;
        double forceConstant;
        bool backboneOnly;
    }ThreadingStruct_wrapper;

    ThreadingStruct_wrapper createThreadingStruct_wrapper(ThreadingStruct & _struct_){
        ThreadingStruct_wrapper _wrap_;
        _wrap_.chainID1      = strdup(_struct_.chainID1.c_str());
        _wrap_.residueStart1 = _struct_.residueStart1.getResidueNumber();
        _wrap_.residueEnd1   = _struct_.residueEnd1.getResidueNumber();
        _wrap_.chainID2      = strdup(_struct_.chainID2.c_str());
        _wrap_.residueStart2 = _struct_.residueStart2.getResidueNumber();
        _wrap_.residueEnd2   = _struct_.residueEnd2.getResidueNumber();
        _wrap_.forceConstant = _struct_.forceConstant;
        _wrap_.backboneOnly  = _struct_.backboneOnly;
        return _wrap_; 
    }
    ThreadingStruct createThreadingStruct(ThreadingStruct_wrapper * _wrap_){
        ThreadingStruct _struct_;
        _struct_.chainID1       = String(_wrap_->chainID1);
        _struct_.residueStart1  = ResidueID(_wrap_->residueStart1,' ');
        _struct_.residueEnd1    = ResidueID(_wrap_->residueEnd1,' ');
        _struct_.chainID2       = String(_wrap_->chainID2);
        _struct_.residueStart2  = ResidueID(_wrap_->residueStart2,' ');
        _struct_.residueEnd2    = ResidueID(_wrap_->residueEnd2,' ');
        _struct_.forceConstant  = _wrap_->forceConstant;
        _struct_.backboneOnly   = _wrap_->backboneOnly;
        return _struct_;     
    }
    /***********************************************************************************/

    /***********************************************************************************/
    typedef struct AtomSpring_wrapper{
        int mmbID;
        const char * atom1;        
        const char * atom2;        
        int          res1;         
        int          res2;         
        const char * chain1;       
        const char * chain2;       
        bool         toGround;     
        bool         tether;       
        double       forceConstant;
        double       deadLength;   
    }AtomSpring_wrapper;

    AtomSpring_wrapper createAtomSpring_wrapper(AtomSpring & _struct_){
        AtomSpring_wrapper _wrap_;
        _wrap_.atom1           = strdup(_struct_.atom1Name.c_str());
        _wrap_.atom2           = strdup(_struct_.atom2Name.c_str());
        _wrap_.res1            = _struct_.atom1Residue.getResidueNumber();
        _wrap_.res2            = _struct_.atom2Residue.getResidueNumber();
        _wrap_.chain1          = strdup(_struct_.atom1Chain.c_str());
        _wrap_.chain2          = strdup(_struct_.atom2Chain.c_str());
        _wrap_.toGround        = _struct_.toGround;
        _wrap_.tether          = _struct_.tether;
        _wrap_.forceConstant   = _struct_.forceConstant;
        _wrap_.deadLength      = _struct_.deadLength;
        return _wrap_; 
    }
    AtomSpring createAtomSpring(AtomSpring_wrapper * _wrap_){
        AtomSpring _struct_;
        _struct_.atom1Name        = String(_wrap_->atom1 );
        _struct_.atom2Name        = String(_wrap_->atom2 );
        _struct_.atom1Residue     = ResidueID(_wrap_->res1,' ');   
        _struct_.atom2Residue     = ResidueID(_wrap_->res2,' ');   
        _struct_.atom1Chain       = String(_wrap_->chain1);
        _struct_.atom2Chain       = String(_wrap_->chain2);
        _struct_.toGround         = _wrap_->toGround;      
        _struct_.tether           = _wrap_->tether;         
        _struct_.forceConstant    = _wrap_->forceConstant;  
        _struct_.deadLength       = _wrap_->deadLength;     
        return _struct_;     
    }
    /***********************************************************************************/

    /***********************************************************************************/   
    /**
    * Description of a DensityStretch.
    */
    typedef struct DensityStretch_wrapper{
        int mmbID;              ///< Id in MMB's ContactContainer
        const char * chainID;           ///< chain ID
        int resStart;        ///< First residue number
        int resEnd;         ///< Last residue number
    } DensityStretch_wrapper;

    /**
    * Create a DensityStretch_wrapper from a MMB's DensityStretch
    */
    DensityStretch_wrapper createDensityStretch_wrapper(DensityStretch & _class_){
        DensityStretch_wrapper _wrap_;
        _wrap_.chainID  = strdup(_class_.getChain().c_str());
        _wrap_.resStart = _class_.getStartResidue().getResidueNumber();
        _wrap_.resEnd   = _class_.getEndResidue().getResidueNumber();
        return _wrap_; 
    }

    /**
    * Create a MMB's DensityStretch object from a DensityStretch_wrapper struct
    */
    DensityStretch createDensityStretch(DensityStretch_wrapper * _wrap_){
        DensityStretch _class_;
        _class_.setChain(_wrap_->chainID);
        _class_.setStartResidue(ResidueID(_wrap_->resStart,' '));
        _class_.setEndResidue(ResidueID(_wrap_->resEnd,' '));
        return _class_; 
    }

    /***********************************************************************************/
    ParameterReader myParameterReader;  ///< ParameterReader object used to set MMB's parameter
    ConstrainedDynamics * myDynamics = NULL; ///< ConstrainedDynamics object used to initialize and run dynamics

    /***********************************************************************************/
    BioPolymerSeq * sequences = NULL;   ///< Array of current sequences
    MobilizerStretch_wrapper * mobilizers = NULL; ///< Buffer for Mobilizers stretches
    MobilizerWithin_wrapper * mobilizersWithin = NULL; ///< Buffer for MobilizersWithin
    Constraint_wrapper * constraints = NULL; ///< Buffer for Constraints
    ContactStretch_wrapper * contacts = NULL; ///< Buffer for Contacts stretches
    ContactWithin_wrapper * contactsWithin = NULL; ///< Buffer for ContactsWithin
    AllResiduesWithin_wrapper * allResiduesWithins = NULL; ///< Buffer for AllResiduesWithins
    IncludeAllNonBondAtomsInResidue_wrapper * includeAllNonBondAtomsInResidues = NULL; ///< Buffer ncludeAllNonBondAtomsInResidue
    BaseInteraction_wrapper * baseInteractions = NULL; ///< Buffer for BaseInteraction_wrapper
    ThreadingStruct_wrapper * threadings = NULL; ///< Buffer for Threadings
    ThreadingStruct_wrapper * gappedThreadings = NULL; ///< Buffer for Threadings
    AtomSpring_wrapper * atomSprings = NULL; ///< Buffer for AtomSprings
    DensityStretch_wrapper * densityStretches = NULL; ///< Buffor for DensityStretches

    string baseInteractionsStrings; ///< String buffer to temporary store baseInteractions commands
    string residuesWithinString; ///< String buffer to store temporary residues within a cutoff
    

    /***********************************************************************************/
    /**
    * Clear the sequences array
    */
    MMB_EXPORT void clearSequences(char * errorString){
        try{
            if(!sequences)
                return;
            for(int i=0; i < sizeof(sequences)/sizeof(sequences[0]); ++i){
                delete sequences[i].chainID;
                delete sequences[i].sequence;
                delete sequences[i].pdbFileName;
            }
            delete sequences;

            sequences = NULL;

        }
        catch(const MMBException& e){
            strcpy(errorString, e.what());
        }
    }

    /**
    * Clear the baseInteractions array
    */
    MMB_EXPORT void clearBaseInteractions(){
        if(!baseInteractions)
            return;
        for(int i=0; i < sizeof(baseInteractions)/sizeof(baseInteractions[0]); ++i)
        {
            delete baseInteractions[i].FirstBPEdge;
            delete baseInteractions[i].SecondBPEdge;
            delete baseInteractions[i].FirstBPChain;
            delete baseInteractions[i].SecondBPChain;
            delete baseInteractions[i].OrientationBP;
        }
        delete baseInteractions;
        baseInteractions = NULL;
    }

    /**
    * Clear the mobilizers array
    */
    MMB_EXPORT void clearMobilizerStretches(){
        if(!mobilizers)
            return;
        for(int i=0; i < sizeof(mobilizers)/sizeof(mobilizers[0]); ++i)
        {
            delete mobilizers[i].chainID;
            delete mobilizers[i].mobility;
        }
        delete mobilizers;
        mobilizers = NULL;
    }

    /**
    * Clear the mobilizersWithin array
    */
    MMB_EXPORT void clearMobilizersWithin(){
        if(!mobilizersWithin)
            return;
        for(int i=0; i < sizeof(mobilizersWithin)/sizeof(mobilizersWithin[0]); ++i)
        {
            delete mobilizersWithin[i].chainID;
            delete mobilizersWithin[i].mobility;
        }
        delete mobilizersWithin;
        mobilizersWithin = NULL;
    }

    /**
    * Clear the contacts array
    */
    MMB_EXPORT void clearContactStretches(){
        if(!contacts)
            return;
        for(int i=0; i < sizeof(contacts)/sizeof(contacts[0]); ++i)
        {
            delete contacts[i].chainID;
            delete contacts[i].ContactScheme;
        }
        delete contacts;
        contacts = NULL;
    }

    /**
    * Clear the contactsWithin array
    */
    MMB_EXPORT void clearContactsWithin(){
        if(!contactsWithin)
            return;
        for(int i=0; i < sizeof(contactsWithin)/sizeof(contactsWithin[0]); ++i)
        {
            delete contactsWithin[i].chainID;
            delete contactsWithin[i].ContactScheme;
        }
        delete contactsWithin;
        contactsWithin = NULL;
    }

    /**
    * Clear the allResiduesWithins array
    */
    MMB_EXPORT void clearAllResiduesWithins(){
        if(!allResiduesWithins)
            return;
        for(int i=0; i < sizeof(allResiduesWithins)/sizeof(allResiduesWithins[0]); ++i)
            delete allResiduesWithins[i].chain;
        delete allResiduesWithins;
        allResiduesWithins = NULL;
    }

    /**
    * Clear the includeAllNonBondAtomsInResidues array
    */
    MMB_EXPORT void clearIncludeAllNonBondAtomsInResidues(){
        if(!includeAllNonBondAtomsInResidues)
            return;
        for(int i=0; i < sizeof(includeAllNonBondAtomsInResidues)/sizeof(includeAllNonBondAtomsInResidues[0]); ++i)
            delete includeAllNonBondAtomsInResidues[i].chain;
        delete includeAllNonBondAtomsInResidues;
        includeAllNonBondAtomsInResidues = NULL;
    }

    /**
    * Clear the mobilizers array
    */
    MMB_EXPORT void clearConstraints(){
        if(!constraints)
            return;
        for(int i=0; i < sizeof(constraints)/sizeof(constraints[0]); ++i)
        {
            delete constraints[i].chain1;
            delete constraints[i].chain2;
            delete constraints[i].atom1;
            delete constraints[i].atom2;
        }
        delete constraints;
        constraints = NULL;
    }

    /**
    * Remove all forces and constraints from MMB.
    */
    MMB_EXPORT void clearForcesAndConstraints(char * errorString){
        try{
            myParameterReader.clearConstraints();
            myParameterReader.clearForces();
        } 
        catch(const MMBException& e){
            strcpy(errorString, e.what());
        }
    }

    /**
    * Initialize MMB with the defaults parameters
    */
    MMB_EXPORT void init(const char * leontisWesthofFileName, char * errorString){
        try{
            if(myDynamics)
                delete myDynamics;
            myDynamics = new ConstrainedDynamics(&myParameterReader);
            myParameterReader.initializeDefaults(leontisWesthofFileName);
            myParameterReader.readPreviousFrameFile = 0;
        } 
        catch(const MMBException& e){
            strcpy(errorString, e.what());
        }
    }

    /**
    * Postinitialization, to be call after all settings, polymers, forces and constraints are set.
    */
    MMB_EXPORT void postInitialize(char * errorString){
        try{
            myParameterReader.postInitialize();
            myDynamics->initializeDumm();
        } 
        catch(const MMBException& e){
            strcpy(errorString, e.what());
        }
    }

    /**
    * Send a command to MMB. 
    * See MMB's documentation for the available commands.
    */
    MMB_EXPORT void command(const char * comstr, char * errorString){
        try{
            myParameterReader.parameterStringInterpreter(comstr);
        } 
        catch(const MMBException& e){
            strcpy(errorString, e.what());
        }
    }

    /**
    * Write a pdb file of the current state of the polymers.
    * @param outputName the output pdb file.
    */
    MMB_EXPORT void writeDefaultPdb(const char * outputName, char * errorString){
        try{
            // BiopolymerClassContainer & biopolys = myParameterReader.myBiopolymerClassContainer;
            ofstream outfile(outputName);
            // biopolys.writeDefaultPdb(outfile);
            myDynamics->writeMMBPDB(outfile);
            outfile.close();
        } 
        catch(const MMBException& e){
            strcpy(errorString, e.what());
        }
    }

    /**
    * Write a pdb file of the current state of the polymers.
    * @param outputName the output pdb file.
    */
    MMB_EXPORT void writePdb(const char * outputName, char * errorString){
        try{
            BiopolymerClassContainer & biopolys = myParameterReader.myBiopolymerClassContainer;
            ofstream outfile(outputName);
            biopolys.writePdb(myDynamics->getCurrentState(), myDynamics->getCompoundSystem(), outfile, false, 0, 0);
            outfile.close();
        } 
        catch(const MMBException& e){
            strcpy(errorString, e.what());
        }
    }

   /**
    * Initialize other molecules and Bonds mobilizers
    * Must be called after initializeBiopolymers()
    */
    MMB_EXPORT void initializeMoleculesAndBonds(char *errorString){
        try{
            myDynamics->initializeMoleculesAndBonds();
        }
        catch(const MMBException& e){
            strcpy(errorString, e.what());
        }
    }

    /**
    * Initialize bonds, the multibody tree and state.
    * Necessary to get a correct structure.
    */
    MMB_EXPORT void initializeBodies(char *errorString){
        try{
            myDynamics->initializeBodies();
        }
        catch(const MMBException& e){
            strcpy(errorString, e.what());
        }
    }

    /**
    * Initialize the particles coordinates and bonds.
    */
    MMB_EXPORT void initializeMolecules(char * errorString){
        try{
            //myParameterReader.readPreviousFrameFile = 0;
            myDynamics->initializeBiopolymersAndCustomMolecules();
            myDynamics->initializeBodies();
        }
        catch(const MMBException& e){
            strcpy(errorString, e.what());
        }
    }

    /**
    * Return the total number of atoms of MMB's current polymers.
    *@return total number of atoms.
    */
    MMB_EXPORT int getSystemNumAtoms(char *errorString){
        try{
            return myParameterReader.myBiopolymerClassContainer.getTotalNumAtoms();
        }
        catch(const MMBException& e){
            strcpy(errorString, e.what());
        }
    }

    /**
    * Fill a linear float array with the coordinates of all polymers.
    * [x0,y0,z0, x1,y1,z1, ..., xn,yn,zn]
    * @param coordArrayv float array to fill. Must be initialized with at least the number of atoms * 3.
    * @see getSystemNumAtoms
    */
    MMB_EXPORT void getSystemCoordinates(MMB_EXPORT void * coordArrayv, char * errorString){
        try{
            float * coordArray = (float *)coordArrayv;

            BiopolymerClassContainer & bpc = myParameterReader.myBiopolymerClassContainer;
            vector<MMBAtomInfo> concatenatedAtomInfoVector = bpc.getConcatenatedAtomInfoVector(myDynamics->getCurrentState());
            for (int i = 0; i < concatenatedAtomInfoVector.size() ; i++) {
                float * outPos = &(coordArray[i * 3]);
                outPos[0] = concatenatedAtomInfoVector[i].position[0] * 10.0;
                outPos[1] = concatenatedAtomInfoVector[i].position[1] * 10.0;
                outPos[2] = concatenatedAtomInfoVector[i].position[2] * 10.0;
            }
        }
        catch(const MMBException& e){
            strcpy(errorString, e.what());
        }
    }

    /**
    * Update the coordinates target based on a PDB structure.
    *
    */
    MMB_EXPORT void matchCoordinatesFromContent(const char * chainID, const char * pdbFileContent, char * errorString)
    {
        try{
            istringstream fileContent(pdbFileContent);
            BiopolymerClassContainer & bpc = myParameterReader.myBiopolymerClassContainer;
            BiopolymerClass & bp = bpc.updBiopolymerClass(chainID);
            bp.matchCoordinates(fileContent,  
                                myParameterReader.matchExact,  
                                myParameterReader.matchIdealized, 
                                myParameterReader.matchOptimize , 
                                myParameterReader.matchHydrogenAtomLocations, 
                                myParameterReader.matchPurineN1AtomLocations, 
                                myParameterReader.guessCoordinates, 
                                myParameterReader.matchingMinimizerTolerance, 
                                myParameterReader.planarityThreshold );
        }
        catch(const MMBException& e){
            strcpy(errorString, e.what());
        }
    }

    /**
    * Initialize the dynamics.
    * The last frame will be write in a file called last.1.pdb, the trajectory in trajectory.1.pdb and current frame in frame.pdb.
    */
    MMB_EXPORT void initDynamics(char * errorString){
        try{
            // cout << "Fic que je veux: " << myParameterReader.lastFrameFileName << endl;
            myParameterReader.lastFrameFileName = "last.1.pdb";
            myParameterReader.currentStage = 1;
            myParameterReader.removeNonPriorityBasePairs(myParameterReader.currentStage);
            
            myParameterReader.outTrajectoryFileName = "trajectory.1.pdb";

            // clear the trajectory file
            ofstream trajFile(myParameterReader.outTrajectoryFileName);
            
            // myParameterReader.readPreviousFrameFile = 0;
            myParameterReader.readInQVector         = 0;
            if (myParameterReader.setRemoveBasePairsInRigidStretch) 
                myParameterReader.removeBasePairsInRigidStretch();

            myDynamics->initializeDynamics();
        } 
        catch(const MMBException& e){
            strcpy(errorString, e.what());
        }
    }

    /**
    * Make all necessary initialization and run all steps of the dynamics as defined by current parameters.
    * Do not call another initialization if you use runDynamics.
    */
    MMB_EXPORT void runDynamics(char * errorString){
        try{
            myDynamics->runDynamics();
        }
        catch(const MMBException& e){
            strcpy(errorString, e.what());
        }
    }

    /**
    * Run all steps of the dynamics as defined by current parameters.
    */
    MMB_EXPORT void runAllSteps(char * errorString){
        try{
            myDynamics->runAllSteps();
        }
        catch(const MMBException& e){
            strcpy(errorString, e.what());
        }
    }

    /**
    * Run one step of the dynamics as defined by current parameters.
    * @return the remaining number of steps.
    */
    MMB_EXPORT int runOneStep(char * errorString){
        try{
            return myDynamics->runOneStep();
        }
        catch(const MMBException& e){
            strcpy(errorString, e.what());
        }
    }

    /**
    * PostDynamics actions. Write the last frame and delete the integrator.
    */
    MMB_EXPORT void postDynamics(char * errorString){
        try{
            myDynamics->postDynamics();
        } 
        catch(const MMBException& e){
            strcpy(errorString, e.what());
        }
    }

    /**
    * Fill the sequences array with the ones in MMB (in BiopolymerClassContainer)
    * @param seqs a pointer to an array of BioPolymerSeq (for external use)
    */
    MMB_EXPORT int getSequences(BioPolymerSeq ** seqs, char * errorString){
        try{
            clearSequences(errorString);

            BiopolymerClassContainer & biopolys = myParameterReader.myBiopolymerClassContainer;
            int nbPoly = biopolys.getNumBiopolymers();
            sequences = new BioPolymerSeq[nbPoly];

            for(int i=0; i < nbPoly; ++i){
                BiopolymerClass & bp = biopolys.updBiopolymerClass(i);
                sequences[i].chainID = strdup(bp.getChainID().c_str());
                sequences[i].sequence = strdup(bp.getSequence().c_str());
                sequences[i].polyType = bp.getBiopolymerType();
                sequences[i].firstResNum = bp.getFirstResidueID().getResidueNumber();
                sequences[i].pdbFileName = strdup(bp.getPdbFileName().c_str());
                sequences[i].loadFromPdb = bp.getLoadFromPdb();
                sequences[i].activePhysics = bp.getActivePhysics();
            }
            *seqs = sequences;
            return nbPoly;
        }
        catch(const MMBException& e){
            strcpy(errorString, e.what());
        }

    }

    /**
    * MMB is set up with the sequences, creating the corresponding biopolymers.
    * Use initializePolymer to initialize the positions from a pdb file. 
    * @param pdbFileName the name of the pdb file to extract the sequences from
    * @return The number of sequences extracted
    */
    MMB_EXPORT void loadSequencesFromPdb(const char * pdbFileName, char * errorString ){
        try{
            // myParameterReader.previousFrameFileName = pdbFileName;
            // myParameterReader.readPreviousFrameFile = 1;
            myParameterReader.loadSequencesFromPdb(pdbFileName);

            // return getSequences(seqs, errorString);
        }
        catch(const MMBException& e){
            strcpy(errorString, e.what());
        }
    }

    /**
    * MMB is set up with the sequences, creating the corresponding biopolymers.
    * Use initializePolymer to initialize the positions from a pdb file. 
    * @param pdbFileName the name of the pdb file to extract the sequences from
    * @return The number of sequences extracted
    */
    MMB_EXPORT void loadSequencesFromPdbContent(const char * pdbFileContent, char * errorString ){
        try{
            istringstream fileContent(pdbFileContent);
            myParameterReader.loadSequencesFromPdb(pdbFileContent);
        }
        catch(const MMBException& e){
            strcpy(errorString, e.what());
        }
    }

    /**
    * Initialize the coordinates of a polymer identified by its chain ID.
    * Use coordinates frome pdbFileName or generate default coordinates if pdbFileName is empty.
    * The chainID and number of atoms of the polymer must match one of the chain in the pdb file.
    * @param chainID the chain id of the polymer to initialize.
    * @param pdbFileName the pdb file to extract the coordinates from
    */
    MMB_EXPORT void initializePolymer(const char * chainID, const char * pdbFileName, char * errorString){
        try{
            // myParameterReader.previousFrameFileName = pdbFileName;
            // myParameterReader.readPreviousFrameFile = 1;
            cout << "initTopPDB " << pdbFileName << endl;
            myDynamics->initializeBiopolymer(chainID, pdbFileName);
        }
        catch(const MMBException& e){
            strcpy(errorString, e.what());
        }
    }

    /**
    * Delete the polymer identified by its chainID
    */
    MMB_EXPORT void deletePolymer(const char * chainID, char * errorString){
        try{
            BiopolymerClassContainer & biopolys = myParameterReader.myBiopolymerClassContainer;
            biopolys.validateChainID(chainID);
            biopolys.deleteBiopolymerClass(chainID);
        }
        catch(const MMBException& e){
            strcpy(errorString, e.what());
        }
    }

    /**
    * Update a polymer's sequence and/or topology type
    * 
    */
    MMB_EXPORT void updatePolymer(const BioPolymerSeq * wrap, char * errorString){
        try {
            BiopolymerClassContainer & biopolys = myParameterReader.myBiopolymerClassContainer;
            BiopolymerClass & bp = biopolys.updBiopolymerClass(wrap->chainID);
            // bp.changeSequence(wrap->sequence);
            bp.setPdbFileName(wrap->pdbFileName);
            bp.setLoadFromPdb(wrap->loadFromPdb);
            bp.setActivePhysics(wrap->activePhysics);
        }
        catch(const MMBException& e){
            strcpy(errorString, e.what());
        }
    }

    /**
    * Find the residues within the radius (nm) of resID
    * @return a string composed of "chainID resID\n" lines
    */
    MMB_EXPORT const char * getResiduesWithin(const char* chainID, int resID, double radius, char * errorString){
        try{
            BiopolymerClassContainer & bcc = myParameterReader.myBiopolymerClassContainer;
            vector< pair<const BiopolymerClass*, const ResidueID*> > resWithin = bcc.getResiduesWithin(chainID, 
                                                                                         ResidueID(resID,' '), 
                                                                                         radius);//,
                                                                                         // myDynamics->getCurrentState());

            stringstream ss;
            vector< pair<const BiopolymerClass*, const ResidueID*> >::iterator it;
            for(it = resWithin.begin(); it != resWithin.end(); it++){
                // cout << __FUNCTION__ << it->first.getChainID() << " " << it->second.getResidueNumber() << endl;
                ss << it->first->getChainID() << " " << it->second->getResidueNumber() << ",";            
            }
            residuesWithinString = ss.str();
            return residuesWithinString.c_str();
        }
        catch(const MMBException& e){
            strcpy(errorString, e.what());
        }
    }

    /**
    * Fill an array of BaseInteraction_wrapper structs with all BaseInteraction in MMB
    * @param baseInteractionsArray the array to fill for external use
    * @return the number of Base Interactions 
    */
    MMB_EXPORT int getBaseInteractions(BaseInteraction_wrapper ** baseInteractionsArray, char * errorString){
        try{
            BasePairContainer & bc = myParameterReader.basePairContainer;

            clearBaseInteractions();
            baseInteractions = new BaseInteraction_wrapper[bc.numBasePairs()];
            for(size_t i=0; i < bc.numBasePairs(); ++i){
                baseInteractions[i] = createBaseInteraction_wrapper(bc.getBasePair(i));
                baseInteractions[i].mmbID = i;
            }

            *baseInteractionsArray = baseInteractions;
            return bc.numBasePairs();
        }
        catch(const MMBException& e){
            strcpy(errorString, e.what());
        }
    }

    /**
    * Return all the base interactions currently in MMB as an array of strings.
    * MMBid ch1 resN1 edge1 ch2 resN2 edge2 orientation
    * @return array of MMB input strings
    */
    MMB_EXPORT const char * getBaseInteractionsStrings(char * errorString){
        try{
            baseInteractionsStrings = myParameterReader.basePairContainer.getBasePairsStrings();
            return baseInteractionsStrings.c_str();
        }
        catch(const MMBException& e){
            strcpy(errorString, e.what());
        }
    }

    /**
    * Update the base pair 'id'.
    * @param id MMB id of the interaction to update
    * @param ch1 chainID of the first residue
    * @param res1 resID of the first residue
    * @param edge1 LeontisWesthof denomination of the first residue interaction type
    * @param ch2 chainID of the second residue
    * @param res2 resID of the second residue
    * @param edge2 LeontisWesthof denomination of the second residue interaction type
    * @param orient bond orientation (Cis or Trans)
    */
    MMB_EXPORT void updateBasePair(int id, 
                        char* ch1, int res1, char* edge1,
                        char* ch2, int res2, char* edge2,
                        char* orient, char* errorString){
        try{
            myParameterReader.updateBasePair(id, ch1, res1, edge1,
                                             ch2, res2, edge2, orient);
        }
        catch(const MMBException& e){
            strcpy(errorString, e.what());
        }
    }

    /**
    * Remove the base interaction 'id' from MMB.
    * Warning: shift the MMB id of the subsequent interactions from -1.
    * @param id MMB id of the interaction to delete
    */
    MMB_EXPORT void deleteBasePair(int id, char * errorString){
        try{
            myParameterReader.basePairContainer.deleteBasePair(id);
        }
        catch(const MMBException& e){
            strcpy(errorString, e.what());
        }
    }

    MMB_EXPORT int getNumSatisfiedBasePairs(char * errorString){
        try{
            return myParameterReader.satisfiedBasePairs;
        }
        catch(const MMBException& e){
            strcpy(errorString, e.what());
        }
    }

    MMB_EXPORT int getNumUnSatisfiedBasePairs(char * errorString){
        try{
            return myParameterReader.unSatisfiedBasePairs;
        }
        catch(const MMBException& e){
            strcpy(errorString, e.what());
        }
    }

    /**
    * Fill an array of MobilizerStretch_wrapper structs with all MobilizerStretches in MMB
    * @param mobilizers the array to fill for external use
    * @return the number of mobilizers 
    */
    MMB_EXPORT int getMobilizerStretches(MobilizerStretch_wrapper ** mobilArray, char * errorString){
        try{
            MobilizerContainer mc = myParameterReader.mobilizerContainer;
            vector<MobilizerStretch> mmbMobilizers = mc.getResidueStretchVector();

            clearMobilizerStretches();
            mobilizers = new MobilizerStretch_wrapper[mmbMobilizers.size()];
            for(size_t i=0; i < mmbMobilizers.size(); ++i){
                mobilizers[i].mmbID     = i;
                mobilizers[i].chainID   = strdup(mmbMobilizers[i].getChain().c_str());
                mobilizers[i].mobility  = strdup(mmbMobilizers[i].getBondMobilityString().c_str());
                mobilizers[i].resStart  = mmbMobilizers[i].getStartResidue().getResidueNumber();
                mobilizers[i].resEnd    = mmbMobilizers[i].getEndResidue().getResidueNumber();
            }

            *mobilArray = mobilizers;
            return mmbMobilizers.size();
        }
        catch(const MMBException& e){
            strcpy(errorString, e.what());
        }
    }

    /**
    * Update the MobilizerStretch 'id'.
    * @param id index in MMB's MobilizerContainer MobilizerStretch vector
    * @param chainID new chain ID
    * @param startRes new start residue id
    * @param endRes new end residue id
    * @param mobility new bondmobility
    */
    MMB_EXPORT void updateMobilizerStretch(int id, char* chainID,
                               int startRes, int endRes,
                               char* mobility,
                               char* errorString){
        try{
            myParameterReader.updateMobilizerStretch(id, chainID,
                                                    startRes, endRes,
                                                    mobility);
        }
        catch(const MMBException& e){
            strcpy(errorString, e.what());
        }
    }

    /**
    * Remove the MobilizerStretch 'id' from MMB.
    * Warning: shift the MMB id of the subsequent mobilizer stretches from -1.
    * @param id MMB id of the interaction to delete
    */
    MMB_EXPORT void deleteMobilizerStretch(int id, char * errorString){
        try{
            myParameterReader.mobilizerContainer.deleteMobilizerStretch(id);
        }
        catch(const MMBException& e){
            strcpy(errorString, e.what());
        }
    }


    /**
    * Fill an array of MobilizerWithin_wrapper structs with all MobilizersWithin in MMB
    * @param mobilizers the array to fill for external use
    * @return the number of mobilizerWithin 
    */
    MMB_EXPORT int getMobilizersWithin(MobilizerWithin_wrapper ** mobilArray, char * errorString){
        try{
            MobilizerContainer mc = myParameterReader.mobilizerContainer;
            vector<MobilizerWithin> mmbMobilizers = mc.getMobilizerWithinVector();

            clearMobilizersWithin();
            mobilizersWithin = new MobilizerWithin_wrapper[mmbMobilizers.size()];
            for(size_t i=0; i < mmbMobilizers.size(); ++i){
                mobilizersWithin[i].mmbID       = i;
                mobilizersWithin[i].chainID     = strdup(mmbMobilizers[i].getChain().c_str());
                mobilizersWithin[i].mobility    = strdup(mmbMobilizers[i].getBondMobilityString().c_str());
                mobilizersWithin[i].resID       = mmbMobilizers[i].getResidue().getResidueNumber();
                mobilizersWithin[i].radius      = mmbMobilizers[i].getRadius();
            }

            *mobilArray = mobilizersWithin;
            return mmbMobilizers.size();
        }
        catch(const MMBException& e){
            strcpy(errorString, e.what());
        }
    }

    /**
    * Update the MobilizerStretch 'id'.
    * @param id index in MMB's MobilizerContainer MobilizerStretch vector
    * @param chainID new chain ID
    * @param res new residue id
    * @param radius new end residue id
    * @param mobility new bondmobility
    */
    MMB_EXPORT void updateMobilizerWithin(int id, char* chainID,
                               int res,
                               double radius,
                               char* mobility,
                               char* errorString){
        try{
            MobilizerContainer mc = myParameterReader.mobilizerContainer;
            mc.updateMobilizerWithin(id, chainID, ResidueID(res,' '), radius, mobility,
                                     myParameterReader.myBiopolymerClassContainer);
        }
        catch(const MMBException& e){
            strcpy(errorString, e.what());
        }
    }

    /**
    * Remove the MobilizerStretch 'id' from MMB.
    * Warning: shift the MMB id of the subsequent mobilizer stretches from -1.
    * @param id MMB id of the interaction to delete
    */
    MMB_EXPORT void deleteMobilizerWithin(int id, char * errorString){
        try{
            myParameterReader.mobilizerContainer.deleteMobilizerWithin(id);
        }
        catch(const MMBException& e){
            strcpy(errorString, e.what());
        }
    }

    /**
    * Get the rootMobilizer for the chainID
    * @param chainID
    * @param rootMobilizer (output)
    */
    MMB_EXPORT void getRootMobilizer(const char * chainID, char * rootMobilizer, char * errorString){
        try{
            BiopolymerClassContainer & bpc = myParameterReader.myBiopolymerClassContainer;
            strcpy(rootMobilizer, bpc.updBiopolymerClass(chainID).getFirstResidueMobilizerType().c_str());
        }
        catch(const MMBException& e){
            strcpy(errorString, e.what());
        }
    }

    /**
    * Set the rootMobilizer for the chainID
    * @param chainID
    * @param rootMobilizer
    */
    MMB_EXPORT void setRootMobilizer(const char * chainID, const char * rootMobilizer, char * errorString){
        try{
            BiopolymerClassContainer & bpc = myParameterReader.myBiopolymerClassContainer;
            bpc.updBiopolymerClass(chainID).setFirstResidueMobilizerType(rootMobilizer);
        }
        catch(const MMBException& e){
            strcpy(errorString, e.what());
        }
    }

    /**
    * Fill an array of Constraint_wrapper structs with all ConstraintClass in MMB
    * @param constraints the array to fill for external use
    * @return the number of constraints 
    */
    MMB_EXPORT int getConstraints(Constraint_wrapper ** constArray, char * errorString){
        try{
            ConstraintToGroundContainer & cc = myParameterReader.constraintToGroundContainer;
            vector<ConstraintClass> mmbConstraints = cc.getConstraintClassVector();

            clearConstraints();
            constraints = new Constraint_wrapper[mmbConstraints.size()];
            for(size_t i=0; i < mmbConstraints.size(); ++i){
                constraints[i].mmbID    = i;
                constraints[i].chain1   = strdup(mmbConstraints[i].getChain1().c_str());
                constraints[i].res1     = mmbConstraints[i].getResidueID1().getResidueNumber();
                constraints[i].atom1    = strdup(mmbConstraints[i].getAtomName1().c_str());
                constraints[i].chain2   = strdup(mmbConstraints[i].getChain2().c_str());
                constraints[i].res2     = mmbConstraints[i].getResidueID2().getResidueNumber();
                constraints[i].atom2    = strdup(mmbConstraints[i].getAtomName2().c_str());
            }

            *constArray = constraints;
            return mmbConstraints.size();
        }
        catch(const MMBException& e){
            strcpy(errorString, e.what());
        }
    }

    /**
    * Remove the Constraint 'id' from MMB.
    * Warning: shift the MMB id of the subsequent constraints from -1.
    * @param id MMB id of the constraint to delete
    */
    MMB_EXPORT void deleteConstraint(int id, char * errorString){
        try{
            ConstraintToGroundContainer & cc = myParameterReader.constraintToGroundContainer;
            cc.deleteConstraintClass(id);
        }
        catch(const MMBException& e){
            strcpy(errorString, e.what());
        }
    }

    /**
    * Update the COnstraint 'id'.
    * @param id index in MMB's ConstraintToGroundContainer's vector
    * @param chainID1 string new chain ID 1
    * @param res1 int new residue id 1
    * @param atomName1 string new atom name 1
    * @param chainID2 string new chain ID 2
    * @param res2 int new residue id 2
    * @param atomName2 string new atom name 2
    */
    MMB_EXPORT void updateConstraint(int id, 
                          char* chainID1, int res1, char* atomName1,
                          char* chainID2, int res2, char* atomName2,
                          char* errorString){
        try{
            ConstraintToGroundContainer & cc = myParameterReader.constraintToGroundContainer;
            if(strcmp(chainID2,"Ground") == 0)
                cc.updateConstraintToVector(id,
                                            chainID1, ResidueID(res1,' '), atomName1,
                                            myParameterReader.myBiopolymerClassContainer);
            else
                cc.updateConstraintToVector(id,
                                            chainID1, ResidueID(res1,' '), atomName1,
                                            chainID2, ResidueID(res2,' '), atomName2,
                                            myParameterReader.myBiopolymerClassContainer);
        }
        catch(const MMBException& e){
            strcpy(errorString, e.what());
        }
    }

    /**
    * Fill an array of ContactStretch_wrapper structs with all ContactStretches in MMB
    * @param contacts the array to fill for external use
    * @return the number of contacts 
    */
    MMB_EXPORT int getContactStretches(ContactStretch_wrapper ** contactArray, char * errorString){
        try{
            ContactContainer mc = myParameterReader.contactContainer;
            vector<ContactStretch> mmbContacts = mc.getResidueStretchVector();

            clearContactStretches();
            contacts = new ContactStretch_wrapper[mmbContacts.size()];
            for(size_t i=0; i < mmbContacts.size(); ++i){
                contacts[i] = createContactStretch_wrapper(mmbContacts[i]);
                contacts[i].mmbID = i;
            }

            *contactArray = contacts;
            return mmbContacts.size();
        }
        catch(const MMBException& e){
            strcpy(errorString, e.what());
        }
    }

    /**
    * Update the ContactStretch 'id'.
    * @param id index in MMB's ContactContainer ContactStretch vector
    * @param chainID new chain ID
    * @param startRes new start residue id
    * @param endRes new end residue id
    * @param mobility new bondmobility
    */
    MMB_EXPORT void updateContactStretch(int id, char* chainID,
                               int startRes, int endRes,
                               char* contactScheme,
                               char* errorString){
        try{
            ContactContainer cc = myParameterReader.contactContainer;
            cc.updateContact(id, chainID,
                             startRes, endRes,
                             contactScheme, myParameterReader.myBiopolymerClassContainer);
        }
        catch(const MMBException& e){
            strcpy(errorString, e.what());
        }
    }

    /**
    * Remove the ContactStretch 'id' from MMB.
    * Warning: shift the MMB id of the subsequent contact stretches from -1.
    * @param id MMB id of the interaction to delete
    */
    MMB_EXPORT void deleteContactStretch(int id, char * errorString){
        try{
            ContactContainer & cc = myParameterReader.contactContainer;
            cc.deleteContact(id);
        }
        catch(const MMBException& e){
            strcpy(errorString, e.what());
        }
    }


    /**
    * Fill an array of ContactWithin_wrapper structs with all ContactsWithin in MMB
    * @param contacts the array to fill for external use
    * @return the number of contactWithin 
    */
    MMB_EXPORT int getContactsWithin(ContactWithin_wrapper ** mobilArray, char * errorString){
        try{
            ContactContainer & mc = myParameterReader.contactContainer;
            vector<ContactWithin> & mmbContacts = mc.getContactWithinVector();

            clearContactsWithin();
            contactsWithin = new ContactWithin_wrapper[mmbContacts.size()];
            for(size_t i=0; i < mmbContacts.size(); ++i){
                contactsWithin[i] = createContactWithin_wrapper(mmbContacts[i]);
                contactsWithin[i].mmbID = i;
            }

            *mobilArray = contactsWithin;
            return mmbContacts.size();
        }
        catch(const MMBException& e){
            strcpy(errorString, e.what());
        }
    }

    /**
    * Update the ContactStretch 'id'.
    * @param id index in MMB's ContactContainer ContactStretch vector
    * @param chainID new chain ID
    * @param res new residue id
    * @param radius new end residue id
    * @param mobility new bondmobility
    */
    MMB_EXPORT void updateContactWithin(int id, char* chainID,
                               int res,
                               double radius,
                               char* contactScheme,
                               char* errorString){
        try{
            ContactContainer & mc = myParameterReader.contactContainer;
            mc.updateContactWithin(id, chainID, res, radius, contactScheme,
                                   myParameterReader.myBiopolymerClassContainer);
        }
        catch(const MMBException& e){
            strcpy(errorString, e.what());
        }
    }

    /**
    * Remove the ContactStretch 'id' from MMB.
    * Warning: shift the MMB id of the subsequent contact stretches from -1.
    * @param id MMB id of the interaction to delete
    */
    MMB_EXPORT void deleteContactWithin(int id, char * errorString){
        try{
            myParameterReader.contactContainer.deleteContactWithin(id);
        }
        catch(const MMBException& e){
            strcpy(errorString, e.what());
        }
    }


    /**
    * Fill an array of AllResiduesWithin_wrapper structs with all AllResiduesWithins in MMB
    * @param arrayOut the array to fill for external use
    * @return the number of AllResiduesWithins 
    */
    MMB_EXPORT int getAllResiduesWithins(AllResiduesWithin_wrapper ** arrayOut, char * errorString){
        try{
            vector<AllResiduesWithin> mmbVec = myParameterReader.includeAllResiduesWithinVector;

            clearAllResiduesWithins();
            allResiduesWithins = new AllResiduesWithin_wrapper[mmbVec.size()];
            for(size_t i=0; i < mmbVec.size(); ++i){
                allResiduesWithins[i] = createAllResiduesWithin_wrapper(mmbVec[i]);
                allResiduesWithins[i].mmbID = i;
            }

            *arrayOut = allResiduesWithins;
            return mmbVec.size();
        }
        catch(const MMBException& e){
            strcpy(errorString, e.what());
        }
    }

    /**
    * Add an includeAllResiduesWihtin command to MMB
    */
    MMB_EXPORT void addAllResiduesWithin(char * chainID, int resID, double radius, char * errorString){
        try{
            myParameterReader.addAllResiduesWithin(chainID, resID, radius);
        }
        catch(const MMBException& e){
            strcpy(errorString, e.what());
        }
    }

    /**
    * Update an includeAllResiduesWihtin command in MMB
    */
    MMB_EXPORT void updateAllResiduesWithin(int index, char * chainID, int resID, double radius, char * errorString){
        try{
            myParameterReader.updateAllResiduesWithin(index, chainID, resID, radius);
        }
        catch(const MMBException& e){
            strcpy(errorString, e.what());
        }
    }

    /**
    * Remove an includeAllResiduesWihtin command from MMB
    */
    MMB_EXPORT void deleteAllResiduesWithin(int index, char * errorString){
        try{
            myParameterReader.deleteAllResiduesWithin(index);
        }
        catch(const MMBException& e){
            strcpy(errorString, e.what());
        }
    }

    /**
    * Update an includeAllNonBondAtomsInResidue command in MMB
    */
    MMB_EXPORT void updateIncludeAllNonBondAtomsInResidue(int index, char * chainID, int resID, char * errorString){
        try{
            myParameterReader.updateIncludeAllNonBondAtomsInResidue(index, chainID, resID);
        }
        catch(const MMBException& e){
            strcpy(errorString, e.what());
        }
    }

    /**
    * Remove an includeAllNonBonAtomsInResidue comman from MMB
    */
    MMB_EXPORT void deleteIncludeAllNonBondAtomsInResidue(int index, char * errorString){
        try{
            myParameterReader.deleteIncludeAllNonBondAtomsInResidue(index);
        }
        catch(const MMBException& e){
            strcpy(errorString, e.what());
        }
    }

    /**
    * Fill an array of IncludeAllNonBondAtomsInResidue_wrapper structs with all IncludeAllNonBondAtomsInResidues in MMB
    * @param arrayOut the array to fill for external use
    * @return the number of IncludeAllNonBondAtomsInResidue 
    */
    MMB_EXPORT int getIncludeAllNonBondAtomsInResidues(IncludeAllNonBondAtomsInResidue_wrapper ** arrayOut, char * errorString){
        try{
            // TEMP commenting while figuring out how to use myParameterReader.physicsContainer
            vector< IncludeAllNonBondAtomsInResidue > & mmbVec = myParameterReader.physicsContainer.residueStretchVector;

            clearIncludeAllNonBondAtomsInResidues();
            includeAllNonBondAtomsInResidues = new IncludeAllNonBondAtomsInResidue_wrapper[mmbVec.size()];
            for(size_t i=0; i < mmbVec.size(); ++i){
                includeAllNonBondAtomsInResidues[i] = createIncludeAllNonBondAtomsInResidue_wrapper( mmbVec[i] );
                includeAllNonBondAtomsInResidues[i].mmbID = i;
            }

            *arrayOut = includeAllNonBondAtomsInResidues;
            return mmbVec.size();
        }
        catch(const MMBException& e){
            strcpy(errorString, e.what());
        }
    }

    /**
    * Print all setting on the standard output.
    */
    MMB_EXPORT void printAllSettings(char * errorString){
        try{
            myParameterReader.printAllSettings();
        }
        catch(const MMBException& e){
            strcpy(errorString, e.what());
        }
    }

    /**
    * Print coordinates of the polymer identified by chainId.
    * @param chainId the polymer id to print
    */
    MMB_EXPORT void printBiopolymerCoordinates(const char * chainId, char * errorString){
        try{
            BiopolymerClass & bpoly = myParameterReader.myBiopolymerClassContainer.updBiopolymerClass(chainId);
            cout << "Bpoly id: " << bpoly.myBiopolymer.getPdbChainId() << endl;
            PdbChain chain(bpoly.myBiopolymer);
            ofstream pdb("test.pdb");
            int n = 1;
            chain.write(pdb, n);
            pdb.close();
            for(size_t i=0; i<chain.getNumResidues(); ++i){
                PdbResidue res = chain.getResidue(Pdb::ResidueIndex(i));
                cerr << res.getPdbResidueNumber() << " " << res.getName() << endl;
                for(size_t j=0; j<res.getNumAtoms(); ++j){
                    PdbAtom atom = res.getAtom(Pdb::AtomIndex(j));
                    cerr << "\t" << atom.getName() << " " << atom.getLocation() << endl; 
                }
            }
        }
        catch(const MMBException& e){
            strcpy(errorString, e.what());
        }
    }


    /**
    * Clear the threadings array
    */
    MMB_EXPORT void clearThreadings(){
        if(!threadings)
            return;
        for(int i=0; i < sizeof(threadings)/sizeof(threadings[0]); ++i)
        {
            delete threadings[i].chainID1;
            delete threadings[i].chainID2;
        }
        delete threadings;
        threadings = NULL;
    }

    /**
    * Fill an array of ThreadingStruct_wrapper structs with all AllResiduesWithins in MMB
    * @param arrayOut the array to fill for external use
    * @return the number of AllResiduesWithins 
    */
    MMB_EXPORT int getThreadingStructs(ThreadingStruct_wrapper ** arrayOut, char * errorString){
        try{
            vector<ThreadingStruct> mmbVec = myParameterReader.atomSpringContainer.getThreadingVector();

            clearThreadings();
            threadings = new ThreadingStruct_wrapper[mmbVec.size()];
            for(size_t i=0; i < mmbVec.size(); ++i){
                threadings[i] = createThreadingStruct_wrapper(mmbVec[i]);
                threadings[i].mmbID = i;
            }

            *arrayOut = threadings;
            return mmbVec.size();
        }
        catch(const MMBException& e){
            strcpy(errorString, e.what());
        }
    }

    /**
    * Add a Threading command to MMB
    */
    MMB_EXPORT void addThreading(ThreadingStruct_wrapper * wrap, char * errorString){
        try{
            ThreadingStruct thread = createThreadingStruct(wrap);
            myParameterReader.atomSpringContainer.addThreading(thread, 
                                                               myParameterReader.myBiopolymerClassContainer);
        }
        catch(const MMBException& e){
            strcpy(errorString, e.what());
        }
    }

    /**
    * Update a Threading command in MMB
    */
    MMB_EXPORT void updateThreading(int index, ThreadingStruct_wrapper * wrap, char * errorString){
        try{
            ThreadingStruct thread = createThreadingStruct(wrap);
            myParameterReader.atomSpringContainer.updateThreading(index, thread, 
                                                                  myParameterReader.myBiopolymerClassContainer);
        }
        catch(const MMBException& e){
            strcpy(errorString, e.what());
        }
    }

    /**
    * Remove a Threading command from MMB
    */
    MMB_EXPORT void deleteThreading(int index, char * errorString){
        try{
            myParameterReader.atomSpringContainer.deleteThreading(index);
        }
        catch(const MMBException& e){
            strcpy(errorString, e.what());
        }
    }

    /**
    * Clear the gappedThreadings array
    */
    MMB_EXPORT void clearGappedThreadings(){
        if(!gappedThreadings)
            return;
        for(int i=0; i < sizeof(gappedThreadings)/sizeof(gappedThreadings[0]); ++i)
        {
            delete gappedThreadings[i].chainID1;
            delete gappedThreadings[i].chainID2;
        }
        delete gappedThreadings;
        gappedThreadings = NULL;
    }

    /**
    * Fill an array of GappedThreadingStruct_wrapper structs with all AllResiduesWithins in MMB
    * @param arrayOut the array to fill for external use
    * @return the number of AllResiduesWithins 
    */
    MMB_EXPORT int getGappedThreadingStructs(ThreadingStruct_wrapper ** arrayOut, char * errorString){
        try{
            vector<ThreadingStruct> mmbVec = myParameterReader.atomSpringContainer.getGappedThreadingVector();

            clearGappedThreadings();
            gappedThreadings = new ThreadingStruct_wrapper[mmbVec.size()];
            for(size_t i=0; i < mmbVec.size(); ++i){
                gappedThreadings[i] = createThreadingStruct_wrapper(mmbVec[i]);
                gappedThreadings[i].mmbID = i;
            }

            *arrayOut = gappedThreadings;
            return mmbVec.size();
        }
        catch(const MMBException& e){
            strcpy(errorString, e.what());
        }
    }

    /**
    * Add a GappedThreading command to MMB
    */
    MMB_EXPORT void addGappedThreading(ThreadingStruct_wrapper * wrap, char * errorString){
        try{
            ThreadingStruct gappedThread = createThreadingStruct(wrap);
            myParameterReader.atomSpringContainer.addGappedThreading(gappedThread, 
                                                               myParameterReader.myBiopolymerClassContainer);
        }
        catch(const MMBException& e){
            strcpy(errorString, e.what());
        }
    }

    /**
    * Update a GappedThreading command in MMB
    */
    MMB_EXPORT void updateGappedThreading(int index, ThreadingStruct_wrapper * wrap, char * errorString){
        try{
            ThreadingStruct gappedThread = createThreadingStruct(wrap);
            myParameterReader.atomSpringContainer.updateGappedThreading(index, gappedThread, 
                                                                  myParameterReader.myBiopolymerClassContainer);
        }
        catch(const MMBException& e){
            strcpy(errorString, e.what());
        }
    }

    /**
    * Remove a GappedThreading command from MMB
    */
    MMB_EXPORT void deleteGappedThreading(int index, char * errorString){
        try{
            myParameterReader.atomSpringContainer.deleteGappedThreading(index);
        }
        catch(const MMBException& e){
            strcpy(errorString, e.what());
        }
    }


    /**
    * Clear the atomsprings array
    */
    MMB_EXPORT void clearAtomSprings(){
        if(!atomSprings)
            return;
        for(int i=0; i < sizeof(atomSprings)/sizeof(atomSprings[0]); ++i)
        {
            delete atomSprings[i].chain1;
            delete atomSprings[i].chain2;
            delete atomSprings[i].atom1;
            delete atomSprings[i].atom2;
        }
        delete atomSprings;
        atomSprings = NULL;
    }

    /**
    * Fill an array of AtomSpring_wrapper structs with all AtomSprings in MMB
    * @param arrayOut the array to fill for external use
    * @return the number of AtomSprings 
    */
    MMB_EXPORT int getAtomSprings(AtomSpring_wrapper ** arrayOut, char * errorString){
        try{
            vector<AtomSpring> mmbVec = myParameterReader.atomSpringContainer.getAtomSpringVector();

            clearAtomSprings();
            atomSprings = new AtomSpring_wrapper[mmbVec.size()];
            for(size_t i=0; i < mmbVec.size(); ++i){
                atomSprings[i] = createAtomSpring_wrapper(mmbVec[i]);
                atomSprings[i].mmbID = i;
            }

            *arrayOut = atomSprings;
            return mmbVec.size();
        }
        catch(const MMBException& e){
            strcpy(errorString, e.what());
        }
    }

    /**
    * Add a AtomSpring command to MMB
    */
    MMB_EXPORT void addAtomSpring(AtomSpring_wrapper * wrap, char * errorString){
        try{
            AtomSpring atomspring = createAtomSpring(wrap);
            myParameterReader.atomSpringContainer.addAtomSpring(atomspring, 
                                                               myParameterReader.myBiopolymerClassContainer);
        }
        catch(const MMBException& e){
            strcpy(errorString, e.what());
        }
    }

    /**
    * Update a AtomSpring command in MMB
    */
    MMB_EXPORT void updateAtomSpring(int index, AtomSpring_wrapper * wrap, char * errorString){
        try{
            AtomSpring spring = createAtomSpring(wrap);
            myParameterReader.atomSpringContainer.updateAtomSpring(index, spring, 
                                                                  myParameterReader.myBiopolymerClassContainer);
        }
        catch(const MMBException& e){
            strcpy(errorString, e.what());
        }
    }

    /**
    * Remove a AtomSpring command from MMB
    */
    MMB_EXPORT void deleteAtomSpring(int index, char * errorString){
        try{
            myParameterReader.atomSpringContainer.deleteAtomSpring(index);
        }
        catch(const MMBException& e){
            strcpy(errorString, e.what());
        }
    }

    /**
    * Clear the DensityStretches array
    */
    MMB_EXPORT void clearDensityStretches(){
        if(!densityStretches)
            return;
        for(int i=0; i < sizeof(densityStretches)/sizeof(densityStretches[0]); ++i)
        {
            delete densityStretches[i].chainID;
        }
        delete densityStretches;
        densityStretches = NULL;
    }

    /**
    * Fill an array of DensityStretch_wrapper structs with all AllResiduesWithins in MMB
    * @param arrayOut the array to fill for external use
    * @return the number of AllResiduesWithins 
    */
    MMB_EXPORT int getDensityStretches(DensityStretch_wrapper ** arrayOut, char * errorString){
        try{
            vector<DensityStretch> mmbVec = myParameterReader.densityContainer.getDensityStretchVector();

            clearDensityStretches();
            densityStretches = new DensityStretch_wrapper[mmbVec.size()];
            for(size_t i=0; i < mmbVec.size(); ++i){
                densityStretches[i] = createDensityStretch_wrapper(mmbVec[i]);
                densityStretches[i].mmbID = i;
            }

            *arrayOut = densityStretches;
            return mmbVec.size();
        }
        catch(const MMBException& e){
            strcpy(errorString, e.what());
        }
    }

    /**
    * Add a Densitystretch command to MMB
    */
    MMB_EXPORT void addDensityStretch(DensityStretch_wrapper * wrap, char * errorString){
        try{
            DensityStretch stretch = createDensityStretch(wrap);
            myParameterReader.densityContainer.add(stretch, myParameterReader.myBiopolymerClassContainer);
        }
        catch(const MMBException& e){
            strcpy(errorString, e.what());
        }
    }

    /**
    * Update a Densitystretch command in MMB
    */
    MMB_EXPORT void updateDensityStretch(int index, DensityStretch_wrapper * wrap, char * errorString){
        try{
            DensityStretch stretch = createDensityStretch(wrap);
            myParameterReader.densityContainer.updateDensityStretch(index, stretch, 
                                                                  myParameterReader.myBiopolymerClassContainer);
        }
        catch(const MMBException& e){
            strcpy(errorString, e.what());
        }
    }

    /**
    * Remove a Densitystretch command from MMB
    */
    MMB_EXPORT void deleteDensityStretch(int index, char * errorString){
        try{
            myParameterReader.densityContainer.deleteDensityStretch(index);
        }
        catch(const MMBException& e){
            strcpy(errorString, e.what());
        }
    }

    MMB_EXPORT void updateParameterReader_wrapper(ParameterReader_wrapper * _wrap_, char * errorString){
        try{
            updateParameterReader_wrapper(myParameterReader, _wrap_);    
        }catch(const MMBException& e){
            strcpy(errorString, e.what());
        }
    }

    MMB_EXPORT void updateParameterReader(ParameterReader_wrapper * _wrap_, char * errorString){
        try{
            updateParameterReader(_wrap_, myParameterReader);    
        }catch(const MMBException& e){
            strcpy(errorString, e.what());
        }
    }

 }
