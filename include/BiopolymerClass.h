/* -------------------------------------------------------------------------- *
 *                           MMB (MacroMoleculeBuilder)                       *
 * -------------------------------------------------------------------------- *
 *                                                                            *
 * Copyright (c) 2011-12 by the Author.                                       *
 * Author: Samuel Flores                                                      *
 *                                                                            *
 * See RNABuilder.cpp for the copyright and usage agreement.                  *
 * -------------------------------------------------------------------------- */

#ifndef BiopolymerClass_H_
#define BiopolymerClass_H_
#define MUTATIONMINORSEPARATOR "-"
#define MUTATIONMAJORSEPARATOR "."

#include "Mutation.h"
#include "SimTKmolmodel.h"
#include <seqan/align.h>
#include "BaseInteractionParameterReader.h"
#include "ConstraintContainer.h"
#include "ReferenceNeighborList.h"

#include <string>
#include <array>
#include <memory>
#include <utility>
#include <vector>

class MobilizerContainer;
class WaterDropletContainer; // forward declaration
class WaterDroplet; // forward declaration
typedef char TChar;                             // character type
typedef seqan::String<TChar> TSequence;                // sequence type
typedef seqan::Align<TSequence, seqan::ArrayGaps> TAlign;
typedef pair <const string, PdbStructure> pdbStructurePairType;
//typedef map pdbStructurePairType PdbStructureMapType ;
typedef map <const string, std::shared_ptr<PdbStructure>> PdbStructureMapType ;


void         printBiopolymerSequenceInfo(const Biopolymer & myBiopolymer) ;
template <class ResidueStretchType>
class ResidueStretchContainer; // Don't you just love forward declarations?


// This data type is for use in atomicPropertyOverrideVector
struct AtomicPropertyOverrideStruct {
    String atomName;
    String property;
    double value;
};

class MMB_EXPORT BiopolymerClass {
private:

    ResidueID   firstResidueID    ;
    String      sequence;
    String      originalSequence;
    String      chainID;
    String      chainPrefix; // This is only used with loadSequencesFromPdb. We have to remember this value for later, when we match coordinates and may need to reopen the input structure file. chainPrefix has to be initialized to "" (we do this in BiopolymerClass::clear()), because loadSequencesFromPdb may never be called.
    String      firstResidueMobilizerType;
    bool        myRenumberPdbResidues;
    bool        proteinCapping;
    vector<MMBAtomInfo> atomInfoVector;
    vector<MMBAtomInfo> ignoreAtomPositionVector;
    vector<ResidueID> residueIDVector; // the element index should match the residue index for fast retrieval
    String      pdbFileName;
    std::shared_ptr<PdbStructure> pdbStructure;
    bool        loadFromPdb;
    void        clear(); // sets all methods to empty values
    //void        validateChainID();
    void        validateChainID() const;
    void        validateResidueNumbersAndInsertionCodes(); // Calls checkResidueNumbersAndInsertionCodes, dies if it finds a problem.
    bool        residueIsPurine (int residueIndex, const String & mySequence);
    bool        residueIsPurine (int residueIndex);
    int         validateSequence() ;  
    int         validateBiopolymerType() const;  
    int         validateProteinCapping();

    // Does not work because of incomplete template class decleration (forward declaration)
    // ResidueStretchContainer<ResidueStretch> inactiveResidueStretches;

    //Temporary
    bool    activePhysics;
    // There is a reason this method is private.  It does not handle the FirstResidue and LastResidue keywords.  Use the BiopolymerClassContainer method, which requires a chain ID.
    ResidueID   residueID(const map<const String,double> &myUserVariables, const char* value) const;
                                       
public:

    Biopolymer  myBiopolymer;

    BiopolymerType::BiopolymerTypeEnum biopolymerType;
    bool    isDNA    () ;
    bool    isRNA    () ;
    String  printOriginalAndRenumberedResidueIDs(const String myPdbId = "XXXX" );
    vector<ResidueID> filterSingleResidueVector (const vector<SingleResidue>  singleResidueVector); // Takes a  vector<SingleResidue> and returns only those ResidueID's belonging to the current BiopolymerClass
    String  getResidueSingleLetterCode(const ResidueID &myResidueID) const {stringstream ss; ss << (myBiopolymer.getResidue(getResidueIndex(myResidueID)).getOneLetterCode ()); return ss.str();}
    void    printTopLevelTransform(){cout<<__FILE__<<":"<<__LINE__<<" Top level transform for chain "<< getChainID() <<" = "<<myBiopolymer.getTopLevelTransform()<<std::endl; }
    int         checkResidueNumbersAndInsertionCodes(); // Same as validateResidueNumbersAndInsertionCodes but if a problem is found it returns 1 rather than dying.
    void    addAtomPositionToIgnore(MMBAtomInfo myAtom){ignoreAtomPositionVector.push_back(myAtom);}
    const ResidueID&   getResidueID    (  int residueIndex) const;
    vector<ResidueID>   getResidueIdVector    (  )  {return residueIDVector;};
    const vector<ResidueID>&   getResidueIdVector (  ) const {return residueIDVector;};
    void    modifyResidue(const BiopolymerModification myBiopolymerModification,Compound  compoundToAdd,  DuMMForceFieldSubsystem & dumm); 
    const String &  getSequence() const {return sequence;}; // gets the sequence //SCF made const 24 Feb 2021
    String  getSequence(const vector <ResidueID> & residueIDVector); // gets the sequence
    String  getSubSequence(const ResidueID &startResidue, const ResidueID &endResidue) const;
    const String&  getChainID() const {return chainID;}
    void    setChainID(String myChainID) {chainID = std::move(myChainID);}
    const String&  getChainPrefix() const {return chainPrefix;}
    void    setChainPrefix(String myChainPrefix ) {chainPrefix = myChainPrefix;} 

    //void    setMutationVector(vector<Mutation> myMutationVector ) { mutationVector = myMutationVector; }
    void    validateMutation( Mutation myMutation);
    void    setOriginalSequence(String); // sets the original sequence
    String  getOriginalSequence() const {cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" originalSequence = "<<originalSequence<<endl; return originalSequence; }; // gets the original (pre-mutation) sequence
    void    setProteinCapping(bool);
    void    setFirstResidueMobilizerType(String myFirstResidueMobilizerType);
    String  getFirstResidueMobilizerType();
    void    setPdbResidueNumbersFromResidueIDVector() ;
    void    renumberPdbResidues(ResidueID firstResidueID);
    void    setBiopolymerType(String); // sets and validates the biopolymerType
    void    setBiopolymerType(BiopolymerType::BiopolymerTypeEnum); // sets and validates the biopolymerType

    void    setPdbFileName(String pdbFileName);
    const String & getPdbFileName() const;

    void    setPdbStructure(const std::shared_ptr<PdbStructure> &structure);
    const std::shared_ptr<PdbStructure> getPdbStructure() const;
    
    void    setLoadFromPdb(bool yesno);
    bool    getLoadFromPdb() const;

    void    setSequence(const String &); // sets and validates the sequence
    void    changeSequence(const String & myNewSequence);
    /*void    renameChain(String newChainID) {
        setChainID(newChainID); myBiopolymer.setPdbChainId(newChainID); 
               
    };*/

    BiopolymerClass(); // sets all properties to empty values.
    BiopolymerClass(String mySequence, String myChainID, ResidueID myFirstResidueNumber, BiopolymerType::BiopolymerTypeEnum myBiopolymerType, bool proteinCapping, bool useNACappingHydroxyls) noexcept;
    BiopolymerClass(const BiopolymerClass &other);
    BiopolymerClass(BiopolymerClass &&other) noexcept;

    BiopolymerClass & operator=(const BiopolymerClass &other);
    BiopolymerClass & operator=(BiopolymerClass &&other) noexcept;
    
    int  initializeBiopolymer(CompoundSystem & system, 
                              bool myProteinCapping, 
                              bool matchExact, bool matchIdealized , 
                              const bool matchOptimize,
                              bool matchHydrogenAtomLocations, 
                              bool matchPurineN1AtomLocations, 
                              bool guessCoordinates,
                              int biopolymerClassIndex, double initialSeparation, 
                              const vector<Displacement> & displacementVector,
                              double matchingMinimizerTolerance, 
                              double myPlanarityThreshold,
                              const vector<SecondaryStructureStretch> &secondaryStructureStretchVector,
                              PdbStructureMapType & pdbStructureMap
                             ) ; //    Should  everything currently done by ConstrainedDynamics::initializeMolecules.  the latter should stop treating biopolymers altogether.  it really should stop treating MonoAtoms as well.

    int     matchCoordinates(const String & inputFileName, 
                             bool matchExact, bool matchIdealized,
                             const bool matchOptimize ,  
                             bool matchHydrogenAtomLocations, 
                             bool matchPurineN1AtomLocations,
                             bool guessCoordinates ,  
                             double matchingMinimizerTolerance, 
                             double myPlanarityThreshold,
			     PdbStructureMapType & pdbStructureMap
			     );   // this parameter sets the out-of-planarity tolerance for identifying planar bonds.  Units: radians.
    int     matchCoordinates(istream & inputFile,
		             PdbStructure::InputType iType,
                             bool matchExact, bool matchIdealized,
                             const bool matchOptimize ,  
                             bool matchHydrogenAtomLocations, 
                             bool matchPurineN1AtomLocations,
                             bool guessCoordinates ,  
                             double matchingMinimizerTolerance, 
                             double myPlanarityThreshold);   // this parameter sets the out-of-planarity tolerance for identifying planar bonds.  Units: radians.
    int     matchCoordinates(const PdbStructure & pdbStructure,  // Return 0 for no issues. Returns 1 for potential problems with maxObservedSinePlaneDeviation
                             bool matchExact, bool matchIdealized,
                             const bool matchOptimize ,  
                             bool matchHydrogenAtomLocations, 
                             bool matchPurineN1AtomLocations,
                             bool guessCoordinates ,  
                             double matchingMinimizerTolerance, 
                             double myPlanarityThreshold);   // this parameter sets the out-of-planarity tolerance for identifying planar bonds.  Units: radians.

    int         getChainLength() const;
    size_t      getNumAtoms();
    ResidueID getFirstResidueID   () const;
    ResidueID getLastResidueID   () const;
    BiopolymerType::BiopolymerTypeEnum  getBiopolymerType() const;
    String  getBiopolymerTypeAsString();
    bool getRenumberPdbResidues (){return myRenumberPdbResidues;}
    void setRenumberPdbResidues (bool tempRenumberPdbResidues);/*{
        myRenumberPdbResidues = tempRenumberPdbResidues;
        std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Just set myRenumberPdbResidues to "<<getRenumberPdbResidues()<<" for chain "<<getChainID()<<std::endl;
        }*/
    bool  getProteinCapping() const {return proteinCapping;}
    //ResidueID   residueID(map<const String,double> myUserVariables, const   char* value);
    ResidueID   residueID(const String &inputString) const;  // this method of converting string to ResidueID has the advantage that it validates against the corresponding biopolymer
    void        validateResidueID(const ResidueID & myResidueID) const;
    void        validateResidueIndex(int myResidueIndex) const;
    void        validateAtomInfoVector();
    bool        hasAtom(  ResidueID myResidueID,   String myAtomName);
    int         validateAtomPathName(  Compound::AtomPathName);
    Compound::AtomPathName atomPathString(const ResidueID &residueID, const String &atomName) const;
    Compound::AtomIndex    atomIndex(const ResidueID &residueID, const String &atomName) const;
    DuMM::AtomIndex    getDuMMAtomIndex(  ResidueID ,   String );
    Vec3        getAtomLocationInMobilizedBodyFrame(  ResidueID myResidueID,   String myAtomName); 
    MobilizedBody & updAtomMobilizedBody(SimbodyMatterSubsystem & matter,   ResidueID myResidueID,   String myAtomName);
    MobilizedBodyIndex getAtomMobilizedBodyIndex( SimbodyMatterSubsystem & matter,   ResidueID myResidueID    ,   String myAtomName);
    Vec3        calcAtomLocationInGroundFrame(const  State & state, const ResidueID &residueID, const String &atomName);
    Vec3        calcDefaultAtomLocationInGroundFrame(const ResidueID &residueID, const String &atomName) const;
    void        loadResidueIDVector();
    void        loadResidueIDVectorAscending(ResidueID firstResidueID);
    ResidueInfo::Index   getResidueIndex(const ResidueID& residueID) const; // called by getPdbResidueName 
    const String&  getPdbResidueName(const ResidueID& residueID) const;
    String      getRepresentativeAtomName() const; // returns the name of an atom which is typically used to represent the entire residue, e.g. CA or C4*.
    double      getRepresentativeAtomMassThreshold() const; // gets a number slightly larger than the maximum expected mass of the representative atom's mobilized body.
    void        setContactParameters(GeneralContactSubsystem & contacts, HuntCrossleyForce & hc, double excludedVolumeStiffness, bool active );
    void        addGeneralSterics (GeneralContactSubsystem & ,ContactSetIndex contactSet, HuntCrossleyForce & hc,SimbodyMatterSubsystem & matter,  double excludedVolumeRadius,double excludedVolumeStiffness    ,    ResidueID startResidue  ,   ResidueID endResidue,   bool endCapsOn ,   bool addHydrogens);
    void        addCustomSterics (GeneralContactSubsystem & contacts, ContactSetIndex contactSet, HuntCrossleyForce & hc,SimbodyMatterSubsystem & matter,LeontisWesthofClass myLeontisWesthofClass,String leontisWesthofInteractionType,   ResidueID startResidue,   ResidueID endResidue,   bool endCapsOn);
    void        setProteinBondMobility (   BondMobility::Mobility  mobility,   ResidueID startResidue,   ResidueID endResidue);
    void        rigidifyTargetedBonds(Compound::AtomTargetLocations  & biopolymerAtomTargets);
    void        setSingleBondMobility(  ResidueID residueID1,  String atomName1,  ResidueID residueID2,   String atomName2,  String mobilityString ); // sets BondMobility for a single bond in the chain.
    Biopolymer & updBiopolymer();
    void        includeNonBondAtom(  ResidueID residueID,   String atomName, State & state, DuMMForceFieldSubsystem & dumm) ;
    ResidueInfo updResidueInfo (  ResidueID residueID) ;
    void        includeAllNonBondAtomsInResidue(  ResidueID residueID, State & state, DuMMForceFieldSubsystem & dumm) ;
    void        includeAllResiduesWithin (  vector<AllResiduesWithin> & includeAllResiduesWithinVector,    vector<IncludeAllNonBondAtomsInResidue> & includeAllNonBondAtomsInResidueVector, const State state);
    void        constrainRigidSegmentsToGround(CompoundSystem & system,  SimbodyMatterSubsystem & matter,State & state, ConstraintToGroundContainer & myConstraintToGroundContainer, bool toGround, ResidueID rootResidue   );
    void        physicsZone(vector<AllResiduesWithin> & myIncludeAllResiduesWithinVector , double radius, SimbodyMatterSubsystem & matter,State & state);
    void        multiplySmallGroupInertia(  ResidueID residueID,   String atomName,   double multiplier, CompoundSystem & system,  SimbodyMatterSubsystem & matter,State & state);
    void        multiplySmallGroupInertia(   double multiplier, CompoundSystem & system, SimbodyMatterSubsystem & matter,State & state) ;
    MMBAtomInfo mmbAtomInfo(const ResidueID &myResidueID, const ResidueInfo::AtomIndex &myResidueInfoAtomIndex, SimbodyMatterSubsystem& matter);
    MMBAtomInfo mmbAtomInfo(const ResidueID &myResidueID, const ResidueInfo::AtomIndex &myResidueInfoAtomIndex, SimbodyMatterSubsystem& matter, DuMMForceFieldSubsystem & dumm);
    //MMBAtomInfo mmbAtomInfo(  ResidueID myResidueID,   ResidueInfo::AtomIndex myResidueInfoAtomIndex,  SimbodyMatterSubsystem& matter, DuMMForceFieldSubsystem & dumm , State & state);
    #ifdef USE_OPENMM
    void        initializeAtomInfoVector(SimbodyMatterSubsystem & matter,  const vector<AtomicPropertyOverrideStruct>  & myAtomicPropertyOverrideVector);
    void        initializeAtomInfoVector(SimbodyMatterSubsystem & matter,DuMMForceFieldSubsystem & dumm, const vector<AtomicPropertyOverrideStruct> & atomicPropertyOverrideVector);
    #endif

    vector<MMBAtomInfo> getAtomInfoVector();
    void printAtomInfoVector(){for (int i = 0 ; i < atomInfoVector.size(); i++) atomInfoVector[i].print(); };
    vector<MMBAtomInfo>  calcAtomInfoVector(ResidueStretch myResidueStretch, SimbodyMatterSubsystem& matter, DuMMForceFieldSubsystem & dumm, const bool includePhosphates = 1);
    void        addRingClosingBond(CovalentBondClass myCovalentBondClass); 
    void        addRingClosingBond( ResidueID residueID1, String atomName1, String bondCenterName1,  ResidueID residueID2, String atomName2,String bondCenterName2, SimTK::BondMobility::Mobility bondMobility) ; 
 
    void        setResidueIDsAndInsertionCodesFromBiopolymer(const Biopolymer & inputBiopolymer, bool endCaps);
    void        setResidueIDsAndInsertionCodesFromBiopolymer(const Biopolymer & inputBiopolymer, Mutation myInsertion, bool endCaps);
    void        setResidueIDsAndInsertionCodesFromBiopolymerWithDeletion(const Biopolymer & oldBiopolymer, ResidueInfo::Index  myDeletionIndex, bool endCaps  );
    void        printBiopolymerInfo() ;
    void        loadMutationVectorFromSequence(); 
    //void        writeMutationFlexibilizers(std::ofstream & output, const int offset,const double radius);
    //void        writeWaterDroplets(std::ofstream & output, const double springConstant, const double radius );
    //void        writeMobilizerWithinMutation(std::ofstream & output,  const double radius  );
    //void        writeMutationBackboneRigidifier(std::ofstream & output, const int offset);
    //void        writePhysicsZones(std::ofstream & output, const int offset);
    //void        writeSubstituteResidueCommands(std::ofstream & output);
    //String      getFormattedMutationsString( String minorSeparator );
    void        setCurrentSequencesFromOriginalSequences();
    //bool    allMutationsDifferFromWildType();
    void        updateMutationResidueTypesFromCurrentSequence() ;
    bool        residueIDLessThanOrEqualTo(ResidueID  residueA, ResidueID  residueB);
    //int         getNumMutationVectorElements() {return mutationVector.size();}; 
    bool        residueIDGreaterThanOrEqualTo(ResidueID  residueA, ResidueID  residueB);
    ResidueID   incrementResidueID(ResidueID & residueID) const;
    ResidueID   decrementResidueID(ResidueID & residueID) const;
    void        setDefaultPeptideDihedralAngle (ResidueID residueID1, ResidueID residueID2, Angle dihedral);
    void    setDefaultPhiAngle (ResidueID residueID, Angle phi);
    void    setDefaultPsiAngle (ResidueID residueID, Angle psi);
    void    setAlphaHelicalDefaultBackboneAngles(ResidueID startResidue, ResidueID endResidue); 
    void        setParallelBetaSheetDefaultBackboneAngles(ResidueID startResidue, ResidueID endResidue); 
    void        setAntiParallelBetaSheetDefaultBackboneAngles(ResidueID startResidue, ResidueID endResidue); 
    
    const int difference(ResidueID  residueA, ResidueID  residueB ) const;
    //ResidueID testSum(ResidueID  oldResidueID, int  increment );
    bool safeSum(ResidueID  inputResidueID, int  increment, ResidueID outputResidueID);
    ResidueID safeSum(ResidueID  inputResidueID, int  increment );
    ResidueID sum(ResidueID  oldResidueID, int  increment ) const;
     
    #ifndef USE_OPENMM
    /** 
    * Get all residues within "distance" (nm) of location
    * @return a vector of ResidueID
    */ 
    std::vector<ResidueID> getResiduesWithin(Vec3 location, double distance);
    #endif
    //void writeDisulphideBridges(std::ofstream & output, SimbodyMatterSubsystem& matter );


    //void lockBiopolymerMobilizedBodies ();

    void        constrainRigidSegmentsToGroundForAllChains(CompoundSystem & system,  SimbodyMatterSubsystem & matter,State & state, ConstraintToGroundContainer & myConstraintToGroundContainer   ) ;
    //void        mutateCurrentSequenceFromMutationVector();

    bool hasResidueStretch(ResidueStretch & residues);

    // Unused for now
    void AddInactiveResidues(ResidueStretch & residues);
    void RemoveInactiveResidues(ResidueStretch & residues);
    vector<ResidueID> getInactiveResiduesVector();

    void setActivePhysics(bool yesno);
    bool getActivePhysics() const;

    template<class ResidueStretchType>
    void selectivelyRemoveResidueStretchFromContainer(const ResidueStretch & residueStretch, ResidueStretchContainer <ResidueStretchType> & residueStretchContainer);

    TAlign    createGappedAlignment(BiopolymerClass otherBiopolymerClass,  double alignmentForcesGapPenalty = -1 );
    int getCorrespondingMutationInCurrentBiopolymer(const BiopolymerClass &otherBiopolymerClass, TAlign align,Mutation mutationInOtherBiopolymer, Mutation & mutationInCurrentBiopolymer);
    int getCorrespondingResidueInCurrentBiopolymer(const BiopolymerClass &otherBiopolymerClass, TAlign align, ResidueID residueIdInOtherBiopolymerClass, ResidueID & correspondingResidueInOtherBiopolymerClass); 
    void sort( vector <ResidueID> & residueIDVector);

};

class MMB_EXPORT BiopolymerClassContainer {
private :
    map <const String, BiopolymerClass> biopolymerClassMap   ;
    vector<Mutation> mutationVector; // this manages the vector of substitution mutations.  Not used for MMB, but breeder uses it extensively.
    PdbStructureMapType  pdbStructureMap;
    //map <const String, PdbStructure> pdbStructureMap;

public:
    #ifdef BuildNtC
    std::vector<std::array<double, 361>> hist;
    std::vector<std::array<double, 361>> prob;
    std::vector<double> counter;
    std::vector<std::array<double, 31>> hist_d;
    std::vector<std::array<double, 31>> prob_d;
    std::vector<double> counter_d;
    #endif
    int     count = 0;

    //BiopolymerClassContainer(){};
    //vector<AtomicPropertyOverrideStruct> atomicPropertyOverrideVector;
    //void        clear(); //: deletes all BiopolymerClass's in biopolymerClassMap, as well as other linked lists
    BiopolymerClassContainer(); //{clear();}; // best to define in .cpp
    const map <const String, BiopolymerClass>& getBiopolymerClassMap () const {return biopolymerClassMap;};
    vector<AtomicPropertyOverrideStruct> atomicPropertyOverrideVector;
    void        clear(); //: deletes all BiopolymerClass's in biopolymerClassMap, as well as other linked lists
    void        writePdbStructureMapDiagnostics(){
        MMBLOG_FILE_FUNC_LINE(INFO,"pdbStructureMap.size() = "<<pdbStructureMap.size()<<std::endl);
        MMBLOG_FILE_FUNC_LINE(INFO,"std::distance(pdbStructureMap.begin(),pdbStructureMap.end()) = "<<std::distance(pdbStructureMap.begin(),pdbStructureMap.end())<<std::endl);
        MMBLOG_FILE_FUNC_LINE(INFO,"pdbStructureMap.empty() = "<<pdbStructureMap.empty()<<std::endl);
    }
    //int         getNumBiopolymers(){return biopolymerClassMap.size();}
    // from HEAD (believe Michal Maly):
    int         getNumBiopolymers() const {return biopolymerClassMap.size();}
    size_t      getTotalNumAtoms();
    vector<SecondaryStructureStretch> secondaryStructureStretchVector;
    void        addBiopolymerClass(String mySequence, 
                                   String myChainID, ResidueID myFirstResidueID, 
                                   BiopolymerType::BiopolymerTypeEnum myBiopolymerType, bool proteinCapping,
                                   String pdbFileName, bool loadFromPdb, 
                                   //String pdbFileName="", bool loadFromPdb=false, 
				   bool useNACappingHydroxyls
                                  ); // validates and sets the listed properties.        
    void        addBiopolymerClass(const String &newChainID, BiopolymerClass newBiopolymerClass) {
        if (hasChainID(newChainID)) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "There is already a biopolymer with chain >"<< newChainID  <<"< !"<<endl);
        }
        biopolymerClassMap.emplace(newChainID, std::move(newBiopolymerClass));
    }
    void        deleteBiopolymerClass(String myChainID);
    void        deleteAllBiopolymerClasses(){biopolymerClassMap.clear();};
    void        deleteAllNonMutatedBiopolymerClasses(); // Deletes all biopolymer classes except those featured in mutationVector
    // Should  everything currently done by ConstrainedDynamics::initializeMolecules.  
    // Should probably split into or call multiple methods for coordinate matching, etc.   
    // The latter should stop treating biopolymers altogether.  it really should stop treating MonoAtoms as well.:
    int         initializeBiopolymers(CompoundSystem & system,
                                      bool myProteinCapping, 
                                      bool matchExact, 
                                      bool matchIdealized ,
                                      const bool matchOptimize ,
                                      bool matchHydrogenAtomLocations, 
                                      bool matchPurineN1AtomLocations,
                                      bool guessCoordinates, 
                                      double initialSeparation,
                                      const vector<Displacement> &displacementVector,
                                      double matchingMinimizerTolerance,
                                      double myPlanarityThreshold
                                     );

    /** 
        Initialize one biopolymer identified by chainID
    */
    int  initializeBiopolymer(const String &chainID, CompoundSystem & system,
                              bool myProteinCapping, bool matchExact, 
                              bool matchIdealized, const bool matchOptimize, 
                              bool matchHydrogenAtomLocations, 
                              bool matchPurineN1AtomLocations,
                              bool guessCoordinates,
                              double initialSeparation, 
                              const vector<Displacement> &displacementVector,
                              double matchingMinimizerTolerance,
                              double myPlanarityThreshold,
                              const vector<SecondaryStructureStretch> &secondaryStructureStretchVector
                             );
    //template<class ResidueStretchTypeA> 
    // This loops through all the MobilzerStretch's in mobilizerContainer, and if these are rigid,  selectively removes them from residueStretchContainer. For example,  residueStretchContainer could be a DensityContainer. If we don't want to apply density forces to Rigid segments, we call this function.
    template<class ResidueStretchType> 
    void selectivelyRemoveRigidMobilizerStretchesFromResidueStretchContainer(MobilizerContainer  & mobilizerContainer, ResidueStretchContainer <ResidueStretchType> & residueStretchContainer);

    String      printOriginalAndRenumberedResidueIDs(const String myPdbId = "XXXX" );
    void        renumberPdbResidues(ResidueID firstResidueID = ResidueID(std::to_string(1))) ; // This renumbers ALL BiopolymerClass's, to start with firstResidueID and increase by consecutive integers from there.
    void        setRenumberPdbResidues (bool myRenumberPdbResidues);

    void        validateAtomInfoVectors();
    const BiopolymerClass &   getBiopolymerClass(const String& myChainID) const;
    BiopolymerClass &   updBiopolymerClass(const String& myChainID);
    const BiopolymerClass &   getBiopolymerClass(int biopolymerClassIndex) const;
    BiopolymerClass &   updBiopolymerClass(int biopolymerClassIndex);
    int                 getBiopolymerClassIndex(String myChainID);

    void        addMutationToVector(Mutation myMutation);
    //void      setBondMobility  (const MobilizerContainer myMobilizerContainer ); 
    void        setBondMobility  (vector<BasePair> & ); 
    void        rigidifyAllChains();
    Vec3        getAtomLocationInMobilizedBodyFrame(String myChainID, ResidueID myResidueID, String myAtomName);
    MobilizedBody &     updAtomMobilizedBody(SimbodyMatterSubsystem &matter, const String &myChainID, const ResidueID &myResidueID, const String &myAtomName);
    //void      addContacts(ContactContainer myContactContainer , GeneralContactSubsystem &,GeneralForceSubsystem &, SimbodyMatterSubsystem &,CompoundSystem & , LeontisWesthofClass & myLeontisWesthofClass ,double excludedVolumeRadius,double excludedVolumeStiffness );
    void        writeDefaultPdb(std::ostream& outputStream);
    void        writePdb(State & state, CompoundSystem & system, std::ostream& outputStream, int modelNumber=1, bool calcEnergy=false, int satisfiedBasePairs=0, int unSatisfiedBasePairs=0);
    bool        hasChainID(const String& chainID) const;
    int         validateChainID(const String& chainID) const;
    Vec3        calcAtomLocationInGroundFrame(const State &state, const String &chainID, const ResidueID &residueId, const String &atomName);
    void        newCalcAxes(const State& state,  LeontisWesthofBondRow myLeontisWesthofBondRow,ResidueID residueID1,ResidueID residueID2,String chain1 , String chain2,Vec3 & xAxisVector1,Vec3 & yAxisVector1, Vec3 & zAxisVector1,Vec3 & xAxisVector2,Vec3 & yAxisVector2 , Vec3 & zAxisVector2,Vec3 & glycosidicNitrogenAtom1LocationInGround,Vec3 & glycosidicNitrogenAtom2LocationInGround, Vec3 & ring1CenterLocationInGround, Vec3 & ring2CenterLocationInGround) ; 
    void        computeCorrection(LeontisWesthofClass &,vector<BaseInteraction> &,State &,SimbodyMatterSubsystem &);
    const String&      getPdbResidueName( const String &chainID, const ResidueID& resID) const;
    void        setSingleBondMobility(const String chainID,const  ResidueID residueID1,const String atomName1,const ResidueID residueID2,const  String atomName2,const String mobilityString ); // sets BondMobility for a single bond in the chain.
    void        setSingleBondMobility(vector<SingleBondMobility>);  
    void        printAllIncludedResidues (const vector<IncludeAllNonBondAtomsInResidue> & includeAllNonBondAtomsInResidueVector);

    #ifdef USE_OPENMM
    std::vector< std::pair<const BiopolymerClass, const ResidueID> > getResiduesWithin(const String & chainID, const ResidueID & resID, double radius, const State & state, OpenMM::NeighborList &  neighborList);
    std::vector< std::pair<const BiopolymerClass, const ResidueID> > getResiduesWithin(const String & chainID, const ResidueID & resID, double radius, OpenMM::NeighborList & neighborList);
    std::vector< std::pair<const BiopolymerClass, const ResidueID> > getResiduesWithin(vector<MMBAtomInfo>& concatenatedAtomInfoVector, const String & chainID, const ResidueID & resID, double radius, OpenMM::NeighborList & neighborList);
    OpenMM::NeighborList getNeighborList(const vector<MMBAtomInfo>& concatenatedAtomInfoVector, double radius);
    void        setNeighborsFromList(vector<MMBAtomInfo>& concatenatedAtomInfoVector, OpenMM::NeighborList& neighborList, double radius);
    template <class type > void findBiopolymerResiduesWithinRadius (const type           & allResiduesWithin      , const State state, vector<SingleResidue> & neighboringResidueVector);
    template <class type2> vector<SingleResidue> findBiopolymerResiduesWithinRadius (const  vector<type2> & allResiduesWithinVector, const State state);
    // Template functions have to be defined in the same compilation unit in which it is used, under many circumstances. Hence the following pass-through or wrapper function:
    vector<SingleResidue>  findBiopolymerResiduesWithinRadius  (const  vector<MobilizerWithin> & allResiduesWithinVector,  const State state);
    void        includeAllResiduesWithin (const vector<AllResiduesWithin> & includeAllResiduesWithinVector,    vector<IncludeAllNonBondAtomsInResidue> & includeAllNonBondAtomsInResidueVector, const State state);
    #endif

    void        includeAllNonBondAtomsInResidues(vector<IncludeAllNonBondAtomsInResidue>  myIncludeAllNonBondAtomsInResidueVector, State & state, DuMMForceFieldSubsystem & dumm) ;
    void        includeNonBondAtoms(  vector<IncludeNonBondAtomInBiopolymerStruct> includeNonBondAtomInBiopolymerVector,  State & state, DuMMForceFieldSubsystem & dumm) ;
    void        includeNonBondAtom(const String chain ,const  ResidueID residue,const  String atomName ,  State & state, DuMMForceFieldSubsystem & dumm) ;
    void        waterDropletAboutResidues (const vector <WaterDropletAboutResidueStruct> waterDropletAboutResidueVector,    WaterDropletContainer & waterDropletContainer );

    void        physicsZone(vector<AllResiduesWithin> & myIncludeAllResiduesWithinVector , double radius, SimbodyMatterSubsystem & matter,State & state);
    void        multiplySmallGroupInertia(const double multiplier, CompoundSystem & system,SimbodyMatterSubsystem & matter,State & state);
    void        initializeAtomInfoVectors(SimbodyMatterSubsystem & matter);
    void        initializeAtomInfoVectors(SimbodyMatterSubsystem & matter,DuMMForceFieldSubsystem & dumm);
    String      extractSequenceFromBiopolymer(const Biopolymer & myBiopolymer, bool endCaps);
    bool        isRNA    (const Biopolymer & inputBiopolymer) ;
    bool        isDNA    (const Biopolymer & inputBiopolymer) ;
    bool        isProtein(const Biopolymer & inputBiopolymer, bool endCaps) ;
    void        loadSequencesFromPdb(const String &inPDBFilename, bool proteinCapping, const String &chainsPrefix, const bool tempRenumberPdbResidues, const bool useNACappingHydroxyls); 
    void        printBiopolymerInfo() ;
    void setResidueIDsAndInsertionCodesFromBiopolymer(const String & chain, const Biopolymer & inputBiopolymer, const bool endCaps);
    ResidueID   residueID(map<const String,double> myUserVariables, const char* value , const String chain); // just like that below, except it can handle user-defined integer variables
    //ResidueID residueID(String inputResidueID, String chain); // this method of converting String to ResidueID has the advantage that it validates against the corresponding biopolymer.  Deprecated!  To be supplanted by the above.
    void constrainAllChainsToEachOther(ConstraintToGroundContainer & constraintToGroundContainer);
    void addConstraintToGround(map<const String,double> myUserVariables,
                               const String inputResidueString,
                               const String chain,
                               const String atomName,
                               ConstraintToGroundContainer & constraintToGroundContainer);
    void addConstraintToGround(map<const String,double> myUserVariables,
                               const String inputResidueString,
                               const String chain,
                               ConstraintToGroundContainer & constraintToGroundContainer);

    void addConstraintToGroundRange(map<const String,double> myUserVariables, const String inputResidueString1, const String inputResidueString2, const String chain, ConstraintToGroundContainer & constraintToGroundContainer);
    void addConstraint(map<const String,double> myUserVariables,
                       const String inputResidueString, const String chain, 
                       const String inputResidueString2, const String chain2, 
                       ConstraintToGroundContainer & constraintToGroundContainer);
    void addConstraint(map<const String,double> myUserVariables,
                       const String atomName1, const String inputResidueString1,const  String chain1, 
                       const String atomName2, const String inputResidueString2,const  String chain2, 
                       ConstraintToGroundContainer & constraintToGroundContainer);
    void addConstraint(map<const String,double> myUserVariables,
               const String atomName1, const String inputResidueString1,const  String chain1, 
               const String atomName2, const String inputResidueString2,const  String chain2, 
               ConstraintType myConstraintType,
               ConstraintToGroundContainer & constraintToGroundContainer);
    void addConstraint(
                   const String chain1,
                   const String chain2,
                   ConstraintToGroundContainer & constraintToGroundContainer);
    void addConstraint(
                   const ResidueID residue1, const String chain1,
                   const ResidueID residue2, const String chain2,
                   ConstraintToGroundContainer & constraintToGroundContainer);

    void constrainRigidSegmentsToGroundForAllChains(CompoundSystem & system,  SimbodyMatterSubsystem & matter,State & state, ConstraintToGroundContainer & myConstraintToGroundContainer  );
    
    void        loadResidueIDVector();
    void        setFirstResidueMobilizerType(const String myFirstResidueMobilizerType);
    void        setContactParameters ( GeneralContactSubsystem & contacts,  HuntCrossleyForce & hc, double excludedVolumeStiffness, bool active );

    void        setOriginalSequence(String myChainID, String myOriginalSequence) {updBiopolymerClass(myChainID).setOriginalSequence(myOriginalSequence);}
    void        setOriginalSequencesFromCurrentSequences()   ; 
    void        replaceBiopolymerWithMutatedBiopolymerClass(const BiopolymerClass & myOldBiopolymerClass,
                                             String & myNewSequence, bool useNACappingHydroxyls = true);
    void        substituteResidue(String myChain , ResidueID myResidue, String mySubstitution, bool proteinCapping);
    void        insertResidue(Mutation myInsertion,  bool proteinCapping) ;
    vector <Mutation> &     updMutationVector(){return mutationVector;}
    vector <Mutation>       getMutationVector(){return mutationVector;}
    void        loadMutationVectorsFromSequence(); 
    void        loadMutationVectorsFromString(String myMutationString); 
    void        writeMutationFlexibilizers(std::ofstream & output, const int offset, const double radius);
    void        writeWaterDroplets(std::ofstream & output, const double springConstant, const double radius );
    void        writeMobilizerWithinMutation(std::ofstream & output,  const double radius  );
    vector<MMBAtomInfo> getConcatenatedAtomInfoVector(bool activeChainsOnly=false);
    vector<MMBAtomInfo> getConcatenatedAtomInfoVector(const State & state,bool activeChainsOnly=false);
    void        printAtomInfoVector();
    void        writeSubstituteResidueCommands(std::ofstream & output);
    int         getNumMutationVectorElements(); 
    String      getFormattedSequencesString();
    String      getFormattedMutationsString( String minorSeparator  );
    String      getFoldxFormattedMutations() ; // This returns a concatenated string of Mutations in the FoldX format
    Mutation    setMutationWildTypeResidueType(Mutation & myMutation);
    Mutation    setMutationWildTypeResidueTypeFromOriginalSequence(Mutation & myMutation);
    void        setCurrentSequencesFromOriginalSequences();
    bool        allMutationsDifferFromWildType();
    void        updateMutationResidueTypesFromCurrentSequence();
    void        substituteResidue(Mutation myMutation, bool safeParameters, bool matchPurineN1AtomLocations, bool proteinCapping );
    void        deleteResidue(Mutation myDeletion,   bool proteinCapping);
    void        setMutationVectorFromString (const std::string mutationString);
    void        addIntraChainInterfaceResidues(String chain, vector<IncludeAllNonBondAtomsInResidue> & myIncludeAllNonBondAtomsInResidueVector , double radius, SimbodyMatterSubsystem & matter,State & state);
    #ifdef USE_OPENMM
    void        createDisulphideBridges(std::ofstream & output);
    void        createDisulphideBridges();
    #endif
    void        loadCysteineAtomInfoVector(vector <MMBAtomInfo> & cysteineAtomInfoVector ) ;
    vector<Mutation> getCompositeMutationVector() {return          mutationVector;}

    void        renameChain(const String &oldChainID, const String &newChainID) {
        BiopolymerClass myBiopolymerClass = getBiopolymerClass(oldChainID); 
        myBiopolymerClass.setChainID(newChainID);
        myBiopolymerClass.myBiopolymer.setPdbChainId(newChainID); 
        // also have to delete and re-add to biopolymerClassMap with the new key
        deleteBiopolymerClass(oldChainID);
        addBiopolymerClass(newChainID, std::move(myBiopolymerClass));
        MMBLOG_FILE_FUNC_LINE(INFO, endl);
    };
    int checkAllResidueNumbersAndInsertionCodes(){for (int i = 0 ; i < getNumBiopolymers(); i++) {if (updBiopolymerClass(i).checkResidueNumbersAndInsertionCodes()) {return 1;}}; return 0; }
    void printTopLevelTransforms(){for (int i = 0 ; i < getNumBiopolymers(); i++) { updBiopolymerClass(i).printTopLevelTransform(); }};  
    /*void setMobilizerTypeForAllChains(const String myMobilizerString, MobilizerContainer & mobilizerContainer){
        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" myMobilizerString = >"<<myMobilizerString<<"< "<<endl;
        for (int i = 0 ; i < myBiopolymerClassContainer.getNumBiopolymers(); i++) {
            String myChainID = myBiopolymerClassContainer.updBiopolymerClass(i).getChainID();
            //String myMobilizerString = parameterStringClass.getString(1);
            cout<<__FILE__<<":"<<__LINE__<<" Adding mobilizer stretch to biopolymer index "<<i<<" , chain "<< myChainID<<endl;
            mobilizerContainer.addMobilizerStretchToVector(
                myChainID,
                myMobilizerString,  
                *this                       
                );
        } // of for
    } */   
};
    
#endif

