/* -------------------------------------------------------------------------- *
 *                           MMB (MacroMoleculeBuilder)                       *
 * -------------------------------------------------------------------------- *
 *                                                                            *
 * Copyright (c) 2011-12 by the Author.                                       *
 * Author: Samuel Flores                                                      *
 *         Alex Tek                                                           *
 *                                                                            *
 * See RNABuilder.cpp for the copyright and usage agreement.                  *
 * -------------------------------------------------------------------------- */

#ifndef ParameterReader_H_
#define ParameterReader_H_

#include <cstdio>
#include <string>
#include <iostream>
#include <vector>

#include "AtomSpringContainer.h"
#include "ContactContainer.h"
#include "ConstraintContainer.h"
#include "DensityMap.h"
#include "WaterDroplet.h"
#include "MonoAtoms.h"
#include "SimTKmolmodel.h"
#include "BaseInteractionParameterReader.h"
#include "BiopolymerClass.h"
#include "ResidueStretchContainer.h"
#include "BasePairContainer.h"
#include "MobilizerContainer.h"
#include "DisplacementContainer.h"
#include "Utils.h"
#include "DensityContainer.h"
#include "MoleculeContainer.h"
#include "Spiral.h"               
//#ifdef BuildNtC
#include "NtC_Class_Container.h"
#include "NTC_FORCE_CLASS.h"
#include "NTC_PARAMETER_READER.h"
//#endif

// Here is where we define GetCurrentDir , an OS dependent function to retrieve the current working directory.
#ifdef _WINDOWS
    #include <direct.h>
    #define GetCurrentDir _getcwd
#else
    #include <unistd.h>
    #define GetCurrentDir getcwd
#endif

//using std::cout;
//using std::endl;

//using namespace SimTK;
//using namespace std  ;
//static double myAtoF(map<const String,double> myUserVariables,const char* value ); 
class MMB_EXPORT ParameterReader { 
private:

    ParameterReader(const ParameterReader &);
    ParameterReader & operator = (const ParameterReader &);

public:

    ParameterReader();
    vector<CovalentBondClass> additionalCovalentBondVector;
    vector<IncludeIntraChainInterface> includeIntraChainInterfaceVector;
    BasePairContainer basePairContainer;    
    map<const ChainResidueIndex, BasePairPartner,twoIndexCmp> basePairPartners;  

       
    // variables previously declared and initialized in Repel.h:
    bool   addAllAtomSterics;
    bool   addAllHeavyAtomSterics;
    bool   addBackboneOxygenForces;
    bool   addProteinBackboneSterics;
    bool   addRNABackboneSterics;
    bool   addSelectedAtoms;          
    bool   addTestSpring;
    bool   useCIFFileFormat;
    std::vector< std::pair < std::string, std::string > > trajectoryFileRemarks;
    std::vector< std::pair < std::string, std::string > > lastFileRemarks;
    bool gemmi_isFirstInStage;
    bool   alignmentForcesIsGapped;
    double alignmentForcesGapPenalty;
    double alignmentForcesDeadLengthFraction;
    double alignmentForcesDeadLength;
    double alignmentForcesDeadLengthIsFractionOfInitialLength;
    double alignmentForcesForceConstant;
    bool   applyC1pSprings;
    vector <BiopolymerModification> biopolymerModificationVector;
    int    calcBaseBodyFramesAtEveryTimeStep;
    bool   calcEnergy      ;
    double totalEnergy;
    double potentialEnergy;
    double kineticEnergy;
    bool   checkSatisfied  ;
    //bool constrainRigidSegments;
    double constraintTolerance;
    InterfaceContainer contactInterfaceContainer;
    bool   guessCoordinates;
    double cutoffRadius    ;
    double cutoffAngle     ;
    double densityAtomFraction;
    double densityNoiseTemperature;
    double densityNoiseScale  ;
    String densityFileName;
    String electroDensityFileName;
    double densityForceConstant;
    bool   densityFitPhosphates; 
    bool   densityNoiseComputeAutocorrelation; 
    bool   densityReportAtEachAtomPosition;
    double electroDensityForceConstant;
    //bool densityMapActivate;
    double excludedVolumeStiffness;
    //String firstResidueMobilizerType;
    int firstStage;
    double fitDefaultTolerance;
    double globalAmberImproperTorsionScaleFactor;
    double globalBondBendScaleFactor ;
    double globalBondStretchScaleFactor ;
    double globalBondTorsionScaleFactor ;
    double globalCoulombScaleFactor;
    double globalGbsaScaleFactor     ;
    double globalVdwScaleFactor;
    double hardSphereStiffnessMultiplier ;
    String  inQVectorFileName;
    double initialSeparation;
    double integratorAccuracy;      
    double integratorStepSize;      
    String integratorType;
    double kbBackboneTorsionGlobalScaleFactor;
    int lastStage;
    String leontisWesthofInFileName;
    bool loadTinkerParameterFile;
    String outQVectorFileName;
    String magnesiumIonChainId;
    double  magnesiumIonRadius ;
    int matchDefaultSkipTopLevelTransform;
    bool matchHydrogenAtomLocations;
    bool matchPurineN1AtomLocations;
    bool matchProteinCarboxylOxygenLocations;
    bool matchExact;    
    bool matchIdealized;
    bool matchOptimize;     
    double matchingMinimizerTolerance;
    int numReportingIntervals ;
    int minimize;
    int monteCarloRun;
    double monteCarloTemperature;
    double monteCarloTemperatureIncrement;
    double nastGlobalBondTorsionScaleFactor;
    double noseHooverTime;
    int numMagnesiumIons;
    String outMonteCarloFileName ;
    String outTrajectoryFileName ;
    //bool physicsWhereYouWantIt;
    float physicsRadius;
    bool  piecewiseRigidify;
    double planarityThreshold; // threshold for considering a bond center to be planar, in rads.
    String potentialType;
    bool prioritize ;
    bool proteinCapping;
    bool useNACappingHydroxyls;
    double excludedVolumeRadius;
    int readInQVector ;
    bool readPreviousFrameFile ;
    int readMagnesiumPositionsFromFile;
    bool removeDensityForcesFromRigidStretches;
    bool removeRigidBodyMomentum;
    double removeMomentumPeriod;
    double reportingInterval  ;
    double restrainingForceConstant;
    double restrainingTorqueConstant;
    bool rigidifyFormedHelices;
    int rigidifyTermini;
    int satisfiedBasePairs ;
    int unSatisfiedBasePairs ;
    double scrubberPeriod  ;
    bool safeParameters;
    int setChiBondMobility;
    int setDefaultStructurePredictionParameters;
    int setDefaultThreadingParameters;

    bool setForceAndStericScrubber    ;
    bool setForceScrubber    ;
    bool setHelicalStacking;
    bool setInitialVelocities;
    bool setLoopBondMobility;//minimize;
    bool setOverallBondMobility;
    bool setRemoveBasePairsInRigidStretch;
    bool setRemoveBasePairsAcrossRigidStretches;
    bool setRepulsiveForce ;
    bool setTemperature;
    double smallGroupInertiaMultiplier;
    Spiral spiral; // This is the helical spiral object, used to fit DNA to icosahedral virus capsids. Defined in Spiral.h
    Vec3   sphericalHelixCenter    ; 
    double sphericalHelixRadius    ; 
    double sphericalHelixStartTheta; 
    double sphericalHelixPhiOffset ;
    double sphericalHelixInterStrandDistance;

    bool stackAllHelicalResidues ;
    String thermostatType;
    String tinkerParameterFileName;
    double twoTransformForceMultiplier;
    bool useFixedStepSize;
    bool useMultithreadedComputation;
    bool useOpenMMAcceleration;
    double vanderWallSphereRadius;
    double velocityRescalingInterval;
    bool verbose;
    int vmdOutput;
    bool waterDropletMake;
    double waterDropletRadius;//Angstroms
    double waterDropletX; //Angstroms
    double waterDropletY; //Angstroms
    double waterDropletZ; //Angstroms

    double waterInertiaMultiplier;
    bool weldToGround;
    double wkdpGlobalBondTorsionScaleFactor;
    bool writeCoordinates;
    bool writeDoublePrecisionTrajectories;
    bool writeFrameFile  ;
    bool writeLastFrameFile  ;
    String workingDirectory    ;

    bool detectConvergence;
    bool converged;
    int convergenceTimeout;
    double convergenceEpsilon;

    BondMobility::Mobility helixBondMobility;
    BondMobility::Mobility loopBondMobility;
    BondMobility::Mobility overallBondMobility;
    BondMobility::Mobility chiBondMobility;

    Vector qVector;
    String lastFrameFileName ;
    String previousFrameFileName ;
    LeontisWesthofClass myLeontisWesthofClass;
    int enforceParallelness  ;
    // end of variables improted from Repel.h

    String sequence;
    String proteinSequence;
    String coarseNucleicAcidSequence;
    int numFirstResidues ;
    int numResetBases    ;
    int numProteinFirstResidues ;
    int numProteinChains;
    int numTemperatures;      
    int numGlobalCoulombScaleFactors;
    int numGlobalVdwScaleFactors;
    double temperature ; 
    double dutyCycle   ; 
    int periodicallyUpdateParameters;
    int currentStage;
    int priority;
    vector<int> residueNumber;
    map<const ChainResidueIndex, int,twoIndexCmp> residueNumberTwo;  

    LeontisWesthofClass     _leontisWesthofClass;  
    //#ifdef BuildNtC 
    NTC_Classes             ntc_classes;
    NTC_PAR_Class           ntc_par_class;
    NTC_FORCE_Class         ntc_force_class;
    NTC_Class_Container ntc_class_container;
    //#endif
    BiopolymerClass         mybiopolymerclass;
    
    mutable map<const String,double> userVariables;
    DensityMap myDensityMap;
    DensityMap myElectroDensityMap;
    MobilizerContainer mobilizerContainer;
    PhysicsContainer physicsContainer;
    ConstraintToGroundContainer constraintToGroundContainer;
    DisplacementContainer displacementContainer;
    AtomSpringContainer atomSpringContainer;
    AtomSpring          dummyAtomSpring; // This is a temporary atomSpring which keeps getting its properties modified by the user prior to being added to AtomSpringContainer. The atomSpringContainer version is not a dummy, it is a real  AtomSpring which will be applied. AtomSpring has a default constructor, so it should not be necessary to initialize.
    BiopolymerClassContainer myBiopolymerClassContainer;
    MoleculeClassContainer   moleculeClassContainer;
    WaterDropletContainer waterDropletContainer;


    map<const String,String> proteinSequences;
    map<const String,String> coarseNucleicAcidSequences;
    map<const String, int> numRigidSegments   ; // scf remove, phased out
    map<const String,int>::iterator firstResidueNumbersIterator;
    //void addRingClosingBond(const String chainID, ResidueID residueID1, String atomName1,String bondCenterName1,  ResidueID residueID2, String atomName2,String bondCenterName2);
    void addC1pSprings (LeontisWesthofClass myLeontisWesthofClass);
    void applyAtomSprings (SimbodyMatterSubsystem & matter, GeneralForceSubsystem & forces, State & state);
    void configureDumm( DuMMForceFieldSubsystem & dumm);
    /*
    static bool aToBool( const char* value ); 
    static bool aToBool( const String& name, const char* value ); 
    */
    //static bool compareUpper( const String& param, const char* symbol );

    vector<BasePair> baseOperationVector;
    ContactContainer contactContainer;
    DensityContainer densityContainer;
    DensityContainer electroDensityContainer;

    vector<SingleBondMobility> singleBondMobilityVector;
    vector<BasePairPartner> basePairPartnerVector;
    //vector<IncludeAllNonBondAtomsInResidue> includeAllNonBondAtomsInResidueVector;
    vector<AllResiduesWithin> includeAllResiduesWithinVector;
    vector<IncludeNonBondAtomInBiopolymerStruct> includeNonBondAtomInBiopolymerVector;
    vector <WaterDropletAboutResidueStruct> waterDropletAboutResidueVector;
    vector<MobilizerDomainsInterface> mobilizerDomainsInterfaceVector;




    void removeBasePairsInRigidStretch ();
    void removeBasePairsAcrossRigidStretches ();

    void printAllSettings (   ostream  & myOstream = std::cout, String remarkString = "") ;
    void printAllSettingsToMMCIF ( std::vector< std::pair < std::string, std::string > > &remarksVec ) ;

    void removeNonPriorityBasePairs (int priorityLevel);
    //int getFirstResidueNumbers(const String myChainId) const ; 

    // int getProteinFirstResidueNumbers(const String myProteinChainId) const ;

    //int getBasePriority(int baseResidueNumber,String baseChain, String basePairingEdge) const ; 
    // int getNumBasePairs() const; 
    void updateBasePair(int index, 
                        String ch1, int res1, String edge1, 
                        String ch2, int res2, String edge2, 
                        String orient);

    void updateMobilizerStretch(int index,
                                String chainId,
                                int startRes,
                                int endRes,
                                String bondMobility);

    void addAllResiduesWithin(String chainID, int resID, double radius);
    void updateAllResiduesWithin(int index, String chainID, int resID, double radius);
    void deleteAllResiduesWithin(int index);

    void updateIncludeAllNonBondAtomsInResidue(int index, const String & chainID, int resID);
    void deleteIncludeAllNonBondAtomsInResidue(int index);

    //int calcHighestPriority();
    //int calcLowestBondingResidue(const String myChainId) ;

    //int calcHighestBondingResidue(const String myChainId);

    void setLeontisWesthofBondRowIndex();

    void parameterStringInterpreter(const String & paramstr);
    void parameterStringInterpreter(const ParameterStringClass & parameterStringClass,
                                   const int readStage = 0,
                                   const bool readAtOneStageOnly = false, 
                                   const bool readOnlyUntilStage = false,
                                   const bool readExcept = false);

    void initializeFromFileOnly(const char * parameterFileName = "./commands.dat" ) ;
    void setFirstAndLastStageAndUseCifFiles(const char * parameterFileName = "./commands.dat" ) ;
    //void setFirstAndLastStage(const char * parameterFileName = "./commands.dat" ) ;

    void loadSequencesFromPdb(const char * pdbFileName, const string & chainsPrefix = "", const bool tempRenumberPdbResidues = 0 );

    //void printRigidSegments();
    // void printBasePairs();
    // void printBaseAssignments(); 

    void postInitialize();

    void clearContainers();
    void clearBiopolymers();
    void clearForces();
    void clearConstraints();
    void initializeDefaults(const char * leontisWesthofInFileName = "./parameters.csv" );

    void initialize(const char * parameterFileName = "./commands.dat" ); 
    MonoAtomsContainer myMonoAtomsContainer;

private: 
    //variables for internal use only:
    int r;  
    int ti;
    int gcsfi;
    int gvsfi;
    int d;
    char * s;
    //int numChains ;
    //int numProteinChains ;
    //int prioritize ;
    //temperature = 300;
    //outQVectorFileName;
    //firstStage = 1;
    //lastStage = 0;// calcHighestPriority();
    //dutyCycle = 1;
    //priority = 0;  //  this will be set in removeNonPriorityBasePairs
};

  
  
// String getResidueName(String chainId,int myResidueNumber,ParameterReader & tempParameterReader, vector<Biopolymer> & tempChain);


#endif                                    
