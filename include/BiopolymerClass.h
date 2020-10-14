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
#include <string>

#include <array>
#include <utility>
#include <vector>
#include "Mutation.h"
#include "SimTKmolmodel.h"
#include <seqan/align.h>
#include "BaseInteractionParameterReader.h"
#include "ConstraintContainer.h"
#include "ReferenceNeighborList.h"
class WaterDropletContainer; // forward declaration
class WaterDroplet; // forward declaration
typedef char TChar;                             // character type
typedef seqan::String<TChar> TSequence;                // sequence type
typedef seqan::Align<TSequence, seqan::ArrayGaps> TAlign;


void         printBiopolymerSequenceInfo(const Biopolymer & myBiopolymer) ;
bool         letterIsRNA    (String);
bool         letterIsDNA    (String);
bool         letterIsProtein(String);
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
    String      firstResidueMobilizerType;
    bool        myRenumberPdbResidues;
    bool        proteinCapping;
    vector<MMBAtomInfo> atomInfoVector;
    vector<MMBAtomInfo> ignoreAtomPositionVector;
    vector<ResidueID> residueIDVector; // the element index should match the residue index for fast retrieval
    String  pdbFileName;
    PdbStructure pdbStructure;
    bool    loadFromPdb;
    void    clear(); // sets all methods to empty values
    void    validateChainID();
    void    validateResidueNumbersAndInsertionCodes(); // Calls checkResidueNumbersAndInsertionCodes, dies if it finds a problem.
    bool    residueIsPurine (int residueIndex, String mySequence);
    bool    residueIsPurine (int residueIndex);
    int     validateSequence() ;  
    int     validateBiopolymerType() const;  
    int     validateProteinCapping();

    // Does not work because of incomplete template class decleration (forward declaration)
    // ResidueStretchContainer<ResidueStretch> inactiveResidueStretches;

    //Temporary
    bool activePhysics;
    // There is a reason this method is private.  It does not handle the FirstResidue and LastResidue keywords.  Use the BiopolymerClassContainer method, which requires a chain ID.
    ResidueID   residueID(map<const String,double> myUserVariables, const   char* value);
                                       
public:

    Biopolymer  myBiopolymer;

    BiopolymerType::BiopolymerTypeEnum biopolymerType;
    const bool  isDNA    () ;
    const bool  isRNA    () ;
    String  printOriginalAndRenumberedResidueIDs(const String myPdbId = "XXXX" );
    vector<ResidueID> filterSingleResidueVector (const vector<SingleResidue>  singleResidueVector); // Takes a  vector<SingleResidue> and returns only those ResidueID's belonging to the current BiopolymerClass
    String  getResidueSingleLetterCode( ResidueID myResidueID){stringstream ss; ss << (myBiopolymer.getResidue(getResidueIndex(myResidueID)).getOneLetterCode ()); return ss.str();} ;
    void    printTopLevelTransform(){cout<<__FILE__<<":"<<__LINE__<<" Top level transform for chain "<< getChainID() <<" = "<<myBiopolymer.getTopLevelTransform()<<std::endl; }
    int         checkResidueNumbersAndInsertionCodes(); // Same as validateResidueNumbersAndInsertionCodes but if a problem is found it returns 1 rather than dying.
    void    addAtomPositionToIgnore(MMBAtomInfo myAtom){ignoreAtomPositionVector.push_back(myAtom);}
    ResidueID   getResidueID    (  int residueIndex)  ;
    vector<ResidueID>   getResidueIdVector    (  )  {return residueIDVector;};
    void    modifyResidue(const BiopolymerModification myBiopolymerModification,Compound  compoundToAdd,  DuMMForceFieldSubsystem & dumm); 
    String  getSequence(){return sequence;}; // gets the sequence
    String  getSequence(vector <ResidueID> & residueIDVector); // gets the sequence
    String  getSubSequence(const ResidueID startResidue, const ResidueID endResidue); 
    String  getChainID() {return chainID;} 
    void    setChainID(String myChainID ) {chainID = myChainID;cout<<__FILE__<<":"<<__LINE__<<" set chainID to : >"<<getChainID()<<"<"<<endl; } 
    //void    setMutationVector(vector<Mutation> myMutationVector ) { mutationVector = myMutationVector; }
    void    validateMutation( Mutation myMutation);
    void    setOriginalSequence(String); // sets the original sequence
    String  getOriginalSequence(){cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" originalSequence = "<<originalSequence<<endl; return originalSequence; }; // gets the original (pre-mutation) sequence
    void    setProteinCapping(bool);
    void    setFirstResidueMobilizerType(String myFirstResidueMobilizerType);
    String  getFirstResidueMobilizerType();
    void    setPdbResidueNumbersFromResidueIDVector() ;
    void    renumberPdbResidues(ResidueID firstResidueID);
    void    setBiopolymerType(String); // sets and validates the biopolymerType
    void    setBiopolymerType(BiopolymerType::BiopolymerTypeEnum); // sets and validates the biopolymerType

    void    setPdbFileName(String pdbFileName);
    String  getPdbFileName();

    void    setPdbStructure( PdbStructure );
    const   PdbStructure getPdbStructure();
    
    void    setLoadFromPdb(bool yesno);
    bool    getLoadFromPdb();

    void    setSequence(String); // sets and validates the sequence
    void    changeSequence(String myNewSequence);
    /*void    renameChain(String newChainID) {
        setChainID(newChainID); myBiopolymer.setPdbChainId(newChainID); 
               
    };*/

    BiopolymerClass(); // sets all properties to empty values.
    BiopolymerClass(String mySequence, String myChainID, ResidueID myFirstResidueNumber, String myBiopolymerType, bool proteinCapping, bool useNACappingHydroxyls   ); // validates and sets the listed properties.          
    
    int  initializeBiopolymer(CompoundSystem & system, 
                              bool myProteinCapping, 
                              bool matchExact, bool matchIdealized , 
                              const bool matchOptimize,
                              bool matchHydrogenAtomLocations, 
                              bool matchPurineN1AtomLocations, 
                              bool guessCoordinates,
                              int biopolymerClassIndex, double initialSeparation, 
                              const vector<Displacement> displacementVector,
                              double matchingMinimizerTolerance, 
                              double myPlanarityThreshold,
                              vector<SecondaryStructureStretch> secondaryStructureStretchVector 
                             ) ; //    Should  everything currently done by ConstrainedDynamics::initializeMolecules.  the latter should stop treating biopolymers altogether.  it really should stop treating MonoAtoms as well.

    int     matchCoordinates(String inputFileName, 
                             bool matchExact, bool matchIdealized,
                             const bool matchOptimize ,  
                             bool matchHydrogenAtomLocations, 
                             bool matchPurineN1AtomLocations,
                             bool guessCoordinates ,  
                             double matchingMinimizerTolerance, 
                             double myPlanarityThreshold);   // this parameter sets the out-of-planarity tolerance for identifying planar bonds.  Units: radians.
    int     matchCoordinates(istream & inputFile, 
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

    int         getChainLength();
    size_t      getNumAtoms();
    String      getChainID() const{ return chainID;};
    ResidueID   getFirstResidueID   ();
    ResidueID   getLastResidueID   ();
    BiopolymerType::BiopolymerTypeEnum  getBiopolymerType() const;
    String  getBiopolymerTypeAsString();
    bool getRenumberPdbResidues (){return myRenumberPdbResidues;}
    void setRenumberPdbResidues (bool tempRenumberPdbResidues);/*{
        myRenumberPdbResidues = tempRenumberPdbResidues;
        std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Just set myRenumberPdbResidues to "<<getRenumberPdbResidues()<<" for chain "<<getChainID()<<std::endl;
        }*/
    bool  getProteinCapping(){return proteinCapping;}
    //ResidueID   residueID(map<const String,double> myUserVariables, const   char* value);
    ResidueID   residueID(String inputString);  // this method of converting string to ResidueID has the advantage that it validates against the corresponding biopolymer
    void        validateResidueID(ResidueID myResidueID    );
    void        validateResidueIndex(int myResidueIndex);
    void        validateAtomInfoVector();
    bool        hasAtom(  ResidueID myResidueID,   String myAtomName);
    int         validateAtomPathName(  Compound::AtomPathName);
    Compound::AtomPathName atomPathString(  ResidueID myResidueID,   String myAtomName);
    Compound::AtomIndex    atomIndex(  ResidueID ,   String );
    DuMM::AtomIndex    getDuMMAtomIndex(  ResidueID ,   String );
    Vec3        getAtomLocationInMobilizedBodyFrame(  ResidueID myResidueID,   String myAtomName); 
    MobilizedBody & updAtomMobilizedBody(SimbodyMatterSubsystem & matter,   ResidueID myResidueID,   String myAtomName);
    MobilizedBodyIndex getAtomMobilizedBodyIndex( SimbodyMatterSubsystem & matter,   ResidueID myResidueID    ,   String myAtomName);
    Vec3        calcAtomLocationInGroundFrame(const  State & ,    ResidueID residueID,   String atomName);   
    Vec3        calcDefaultAtomLocationInGroundFrame(   ResidueID residueID,   String atomName);
    void        loadResidueIDVector();
    void        loadResidueIDVectorAscending(ResidueID firstResidueID);
    const   ResidueInfo::Index   getResidueIndex(   ResidueID residueID  ); // called by getPdbResidueName 
    String  getPdbResidueName(   ResidueID residueID);
    String      getRepresentativeAtomName(); // returns the name of an atom which is typically used to represent the entire residue, e.g. CA or C4*.
    double      getRepresentativeAtomMassThreshold(); // gets a number slightly larger than the maximum expected mass of the representative atom's mobilized body.
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
    MMBAtomInfo mmbAtomInfo(  ResidueID myResidueID,   ResidueInfo::AtomIndex myResidueInfoAtomIndex,  SimbodyMatterSubsystem& matter  );
    MMBAtomInfo mmbAtomInfo(  ResidueID myResidueID,   ResidueInfo::AtomIndex myResidueInfoAtomIndex,  SimbodyMatterSubsystem& matter, DuMMForceFieldSubsystem & dumm );
    //MMBAtomInfo mmbAtomInfo(  ResidueID myResidueID,   ResidueInfo::AtomIndex myResidueInfoAtomIndex,  SimbodyMatterSubsystem& matter, DuMMForceFieldSubsystem & dumm , State & state);
    #ifdef USE_OPENMM
    void        initializeAtomInfoVector(SimbodyMatterSubsystem & matter,  const vector<AtomicPropertyOverrideStruct>  & myAtomicPropertyOverrideVector);
    void        initializeAtomInfoVector(SimbodyMatterSubsystem & matter,DuMMForceFieldSubsystem & dumm, const vector<AtomicPropertyOverrideStruct> & atomicPropertyOverrideVector);
    #endif

    const       vector<MMBAtomInfo> getAtomInfoVector();
    const       void printAtomInfoVector(){for (int i = 0 ; i < atomInfoVector.size(); i++) atomInfoVector[i].print(); };
    const       vector<MMBAtomInfo>  calcAtomInfoVector(ResidueStretch myResidueStretch, SimbodyMatterSubsystem& matter, DuMMForceFieldSubsystem & dumm, const bool includePhosphates = 1);
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
    ResidueID   incrementResidueID(ResidueID & residueID);
    ResidueID   decrementResidueID(ResidueID & residueID);
    void        setDefaultPeptideDihedralAngle (ResidueID residueID1, ResidueID residueID2, Angle dihedral);
    void    setDefaultPhiAngle (ResidueID residueID, Angle phi);
    void    setDefaultPsiAngle (ResidueID residueID, Angle psi);
    void    setAlphaHelicalDefaultBackboneAngles(ResidueID startResidue, ResidueID endResidue); 
    void        setParallelBetaSheetDefaultBackboneAngles(ResidueID startResidue, ResidueID endResidue); 
    void        setAntiParallelBetaSheetDefaultBackboneAngles(ResidueID startResidue, ResidueID endResidue); 
    
    int difference(ResidueID  residueA, ResidueID  residueB );
    //ResidueID testSum(ResidueID  oldResidueID, int  increment );
    bool safeSum(ResidueID  inputResidueID, int  increment, ResidueID outputResidueID);
    ResidueID safeSum(ResidueID  inputResidueID, int  increment );
    ResidueID sum(ResidueID  oldResidueID, int  increment );
     
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
    void selectivelyRemoveResidueStretchFromContainer(ResidueStretch & residueStretch, ResidueStretchContainer <ResidueStretchType> & residueStretchContainer)
    {    
        // This command crops or deletes residue stretches in the range "residueStretch" from residueStretchVector.  This was intended to cancel any modifications to certain resiude stretches.    
        // We treat three cases:
        // 1. residueStretchVector[i] is a subset of (or is identical to) residueStretch
        //        -> erase residueStretchVector[i]
        // 2. residueStretch is a subset of residueStretchVector[i], with neither endpoint in common, splitting residueStretchVector[i] in two
        //        -> split residueStretchVector[i] into two disjoint residue stretches
        // 3. residueStretch is a subset of residueStretchVector[i], but the start point of residueStretch coincides with that of residueStretchVector[i] .
        //        -> trim  residueStretchVector[i] on left 
        // 4. residueStretch is a subset of residueStretchVector[i], but the end point of residueStretch coincides with that of residueStretchVector[i] .
        //        -> trim  residueStretchVector[i] on right
        // 5. residueStretch and residueStretchVector[i] overlap, with residueStretch starting before residueStretchVector[i].
        //        -> trim  residueStretchVector[i] on left 
        // 6. residueStretch and residueStretchVector[i] overlap, with residueStretchVector[i] starting before residueStretch.
        //        -> trim  residueStretchVector[i] on right
        //const int ResidueStretchContainer::getNumResidueStretches();
        MMBLOG_FILE_FUNC_LINE(INFO, "the Default stretch is :"<<endl);
        residueStretch.printStretch();
        MMBLOG_FILE_FUNC_LINE(INFO, "Now checking "<<residueStretchContainer.getNumResidueStretches()<<" stretches: "<<endl);
        for (int i = 0; i < residueStretchContainer.getNumResidueStretches(); i++) 
        {
            residueStretchContainer.residueStretchVector[i].printStretch();   

            if (residueStretchContainer.residueStretchVector[i].getChain().compare((residueStretch.getChain() )) != 0) {continue;} // in other words, only make modificatiosn to residueStretchContainer if chain ID's match.
            else if ((residueStretch.getStartResidue() <= residueStretchContainer.residueStretchVector[i].getStartResidue()) &&
                (residueStretch.getEndResidue() >= residueStretchContainer.residueStretchVector[i].getEndResidue()))
               {   //case = 1
                   residueStretchContainer.residueStretchVector.erase(residueStretchContainer.residueStretchVector.begin() + i);
                   i--; // vector has been shortened, so make sure we don't skip the next residueStretchContainer.residueStretchVector[i].
                   if (i < -1) {MMBLOG_FILE_FUNC_LINE(CRITICAL, "Unexplained error!"<<endl);}
               }
            else if ((residueStretch.getStartResidue() >  residueStretchContainer.residueStretchVector[i].getStartResidue()) &&
                (residueStretch.getEndResidue() <  residueStretchContainer.residueStretchVector[i].getEndResidue()))
               {   // case = 2 ;
                   MMBLOG_FILE_FUNC_LINE(INFO, endl);
                   MobilizerStretch secondResidueStretch = residueStretchContainer.residueStretchVector[i];
                   ResidueID tempStartResidueID = (residueStretch).getStartResidue(); // getStartResidue() returns a temporary, whereas decrementResidueID expects a reference. can't convert a temporary to a reference.  This is because decrementResidueID might (and will!) try to modify ResidueID (as the name of the function suggests!).
                   //residueStretchContainer.residueStretchVector[i].setEndResidue(decrementResidueID((residueStretch).getStartResidue() ));
                   residueStretchContainer.residueStretchVector[i].setEndResidue(decrementResidueID(tempStartResidueID));//((residueStretch).getStartResidue() )));
                   MMBLOG_FILE_FUNC_LINE(INFO, "Just decreased endpoint of stretch "<<i<<".  New stretch is:"<<endl);
                   residueStretchContainer.residueStretchVector[i].printStretch();
                   ResidueID tempEndResidueID = (residueStretch).getEndResidue();
                   secondResidueStretch.setStartResidue(incrementResidueID(tempEndResidueID));//  residueStretch.getEndResidue()));
                   residueStretchContainer.addResidueStretchToVector(secondResidueStretch);
                   MMBLOG_FILE_FUNC_LINE(INFO, "Just added new  stretch :"<<endl);
                   residueStretchContainer.residueStretchVector[residueStretchContainer.getNumResidueStretches()-1].printStretch();
                   MMBLOG_FILE_FUNC_LINE(INFO, "Moving on to check next stretch. "<<endl);


               }
            else if ((residueStretch.getStartResidue() == residueStretchContainer.residueStretchVector[i].getStartResidue()) &&
                (residueStretch.getEndResidue() <  residueStretchContainer.residueStretchVector[i].getEndResidue()))
               {   // case = 3;
                   MMBLOG_FILE_FUNC_LINE(INFO, "Case 3"<<endl);
                   ResidueID tempEndResidueID = (residueStretch).getEndResidue();
                   residueStretchContainer.residueStretchVector[i].setStartResidue(incrementResidueID(tempEndResidueID));//residueStretch.getEndResidue() ))  ;
               }
            else if ((residueStretch.getEndResidue() == residueStretchContainer.residueStretchVector[i].getEndResidue()) &&
                (residueStretch.getStartResidue() >  residueStretchContainer.residueStretchVector[i].getStartResidue()))
               {   // case = 4;
                   MMBLOG_FILE_FUNC_LINE(INFO, "Case 4"<<endl);
                   
                   ResidueID tempStartResidueID = (residueStretch).getStartResidue();
                   residueStretchContainer.residueStretchVector[i].setEndResidue(decrementResidueID(tempStartResidueID));//residueStretch.getStartResidue()));
               }
            else if ((residueStretch.getStartResidue() <   residueStretchContainer.residueStretchVector[i].getStartResidue()) &&
                (residueStretch.getEndResidue()        >=  residueStretchContainer.residueStretchVector[i].getStartResidue()) &&
                     (residueStretch.getEndResidue()        <   residueStretchContainer.residueStretchVector[i].getEndResidue()))
            {   // case = 5;
                MMBLOG_FILE_FUNC_LINE(INFO, "Case 5"<<endl);
                
                ResidueID tempEndResidueID = (residueStretch).getEndResidue();
                residueStretchContainer.residueStretchVector[i].setStartResidue(incrementResidueID(tempEndResidueID));//residueStretch.getEndResidue()))  ;
            }
            else if ((residueStretch.getEndResidue() >  residueStretchContainer.residueStretchVector[i].getEndResidue()) &&
                     (residueStretch.getStartResidue() >  residueStretchContainer.residueStretchVector[i].getStartResidue())     &&
                     (residueStretch.getStartResidue() <=  residueStretchContainer.residueStretchVector[i].getEndResidue()))
            {    // case = 6;
                MMBLOG_FILE_FUNC_LINE(INFO, "Case 6"<<endl);
                
                ResidueID tempStartResidueID = (residueStretch).getStartResidue();
                residueStretchContainer.residueStretchVector[i].setEndResidue(decrementResidueID(tempStartResidueID));//  residueStretch.getStartResidue()));
            }
            else if (residueStretch.getEndResidue() < residueStretchContainer.residueStretchVector[i].getStartResidue()) {} // do nothing, stretches are disjoint
            else if (residueStretch.getStartResidue() > residueStretchContainer.residueStretchVector[i].getEndResidue()) {} // do nothing, stretches are disjoint
            else {
                // this should never happen
                MMBLOG_FILE_FUNC_LINE(CRITICAL, "Unexplained error!"<<endl);
                }
        }
    }
    TAlign    createGappedAlignment(BiopolymerClass otherBiopolymerClass,  double alignmentForcesGapPenalty = -1 );
    int getCorrespondingMutationInCurrentBiopolymer(BiopolymerClass otherBiopolymerClass, TAlign align,Mutation mutationInOtherBiopolymer, Mutation & mutationInCurrentBiopolymer);
    int getCorrespondingResidueInCurrentBiopolymer(BiopolymerClass otherBiopolymerClass, TAlign align, ResidueID residueIdInOtherBiopolymerClass, ResidueID & correspondingResidueInOtherBiopolymerClass); 
    void sort( vector <ResidueID> & residueIDVector);

};

class MMB_EXPORT BiopolymerClassContainer {
private :
    map <const String, BiopolymerClass> biopolymerClassMap   ;
    vector<Mutation> mutationVector; // this manages the vector of substitution mutations.  Not used for MMB, but breeder uses it extensively.
    map <const String, PdbStructure> pdbStructureMap;

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

    BiopolymerClassContainer(){};
    map <const String, BiopolymerClass> getBiopolymerClassMap () const {return biopolymerClassMap;};
    vector<AtomicPropertyOverrideStruct> atomicPropertyOverrideVector;
    void        clear(); //: deletes all BiopolymerClass's in biopolymerClassMap, as well as other linked lists
    int         getNumBiopolymers(){return biopolymerClassMap.size();}
    size_t      getTotalNumAtoms();
    vector<SecondaryStructureStretch> secondaryStructureStretchVector;
    void        addBiopolymerClass(String mySequence, 
                                   String myChainID, ResidueID myFirstResidueID, 
                                   String myBiopolymerType, bool proteinCapping,
                                   String pdbFileName, bool loadFromPdb, 
                                   //String pdbFileName="", bool loadFromPdb=false, 
				   bool useNACappingHydroxyls
                                  ); // validates and sets the listed properties.        
    void        addBiopolymerClass(String newChainID, BiopolymerClass newBiopolymerClass) {
        if (hasChainID(newChainID)) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "There is already a biopolymer with chain >"<< newChainID  <<"< !"<<endl);
        }
        biopolymerClassMap[newChainID] = newBiopolymerClass;}
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
                                      const vector<Displacement> displacementVector,
                                      double matchingMinimizerTolerance,
                                      double myPlanarityThreshold
                                     );

    /** 
        Initialize one biopolymer identified by chainID
    */
    int  initializeBiopolymer(String chainID, CompoundSystem & system,
                              bool myProteinCapping, bool matchExact, 
                              bool matchIdealized, const bool matchOptimize, 
                              bool matchHydrogenAtomLocations, 
                              bool matchPurineN1AtomLocations,
                              bool guessCoordinates,
                              double initialSeparation, 
                              const vector<Displacement> displacementVector,
                              double matchingMinimizerTolerance,
                              double myPlanarityThreshold,
                              vector<SecondaryStructureStretch> secondaryStructureStretchVector
                             );
    String      printOriginalAndRenumberedResidueIDs(const String myPdbId = "XXXX" );
    void        renumberPdbResidues(ResidueID firstResidueID = ResidueID(std::to_string(1))) ; // This renumbers ALL BiopolymerClass's, to start with firstResidueID and increase by consecutive integers from there.
    void        setRenumberPdbResidues (bool myRenumberPdbResidues);

    void        validateAtomInfoVectors();
    BiopolymerClass &   updBiopolymerClass(String myChainID);
    int                 getBiopolymerClassIndex(String myChainID);
    BiopolymerClass &   updBiopolymerClass(int biopolymerClassIndex);

    void        addMutationToVector(Mutation myMutation);
    //void      setBondMobility  (const MobilizerContainer myMobilizerContainer ); 
    void        setBondMobility  (vector<BasePair> & ); 
    void        rigidifyAllChains();
    Vec3        getAtomLocationInMobilizedBodyFrame(String myChainID, ResidueID myResidueID, String myAtomName);
    MobilizedBody &     updAtomMobilizedBody(SimbodyMatterSubsystem & matter, String myChainID, ResidueID myResidueID, String myAtomName);
    //void      addContacts(ContactContainer myContactContainer , GeneralContactSubsystem &,GeneralForceSubsystem &, SimbodyMatterSubsystem &,CompoundSystem & , LeontisWesthofClass & myLeontisWesthofClass ,double excludedVolumeRadius,double excludedVolumeStiffness );
    void        writeDefaultPdb(std::ostream& outputStream);
    void        writePdb(State & state, CompoundSystem & system, std::ostream& outputStream, int modelNumber=1, bool calcEnergy=false, int satisfiedBasePairs=0, int unSatisfiedBasePairs=0);
    bool        hasChainID(String);
    int         validateChainID(String);
    Vec3        calcAtomLocationInGroundFrame(const State & , String chainID, ResidueID , String );   
    void        newCalcAxes(const State& state,  LeontisWesthofBondRow myLeontisWesthofBondRow,ResidueID residueID1,ResidueID residueID2,String chain1 , String chain2,Vec3 & xAxisVector1,Vec3 & yAxisVector1, Vec3 & zAxisVector1,Vec3 & xAxisVector2,Vec3 & yAxisVector2 , Vec3 & zAxisVector2,Vec3 & glycosidicNitrogenAtom1LocationInGround,Vec3 & glycosidicNitrogenAtom2LocationInGround, Vec3 & ring1CenterLocationInGround, Vec3 & ring2CenterLocationInGround) ; 
    void        computeCorrection(LeontisWesthofClass &,vector<BaseInteraction> &,State &,SimbodyMatterSubsystem &);
    String      getPdbResidueName( const String chainID, ResidueID resID);
    void        setSingleBondMobility(const String chainID,const  ResidueID residueID1,const String atomName1,const ResidueID residueID2,const  String atomName2,const String mobilityString ); // sets BondMobility for a single bond in the chain.
    void        setSingleBondMobility(vector<SingleBondMobility>);  
    void        printAllIncludedResidues (vector<IncludeAllNonBondAtomsInResidue> & includeAllNonBondAtomsInResidueVector );

    #ifdef USE_OPENMM
    std::vector< std::pair<const BiopolymerClass, const ResidueID> > getResiduesWithin(const String & chainID, const ResidueID & resID, double radius, const State & state, OpenMM::NeighborList &  neighborList);
    std::vector< std::pair<const BiopolymerClass, const ResidueID> > getResiduesWithin(const String & chainID, const ResidueID & resID, double radius, OpenMM::NeighborList & neighborList);
    std::vector< std::pair<const BiopolymerClass, const ResidueID> > getResiduesWithin(vector<MMBAtomInfo>& concatenatedAtomInfoVector, const String & chainID, const ResidueID & resID, double radius, OpenMM::NeighborList & neighborList);
    OpenMM::NeighborList getNeighborList(const vector<MMBAtomInfo>& concatenatedAtomInfoVector, double radius);
    void        setNeighborsFromList(vector<MMBAtomInfo>& concatenatedAtomInfoVector, OpenMM::NeighborList& neighborList, double radius);
    template <class type > void findBiopolymerResiduesWithinRadius (const type           & allResiduesWithin      , const State state, vector<SingleResidue> & neighboringResidueVector);
    template <class type2> vector<SingleResidue> findBiopolymerResiduesWithinRadius (const  vector<type2> & allResiduesWithinVector, const State state);
    vector<SingleResidue> findBiopolymerResiduesWithinRadius (const  vector<Mutation> & allResiduesWithinVector, const State state){return findBiopolymerResiduesWithinRadius( allResiduesWithinVector, state); };
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
    const bool  isRNA    (const Biopolymer & inputBiopolymer) ;
    const bool  isDNA    (const Biopolymer & inputBiopolymer) ;
    const bool  isProtein(const Biopolymer & inputBiopolymer, bool endCaps) ;
    void        loadSequencesFromPdb(String inPDBFilename,bool proteinCapping, const String & chainsPrefix , const bool tempRenumberPdbResidues  , const bool useNACappingHydroxyls); 
    const PdbStructure & getPdbStructure(String fileName);
    void        printBiopolymerInfo() ;
    void setResidueIDsAndInsertionCodesFromBiopolymer(const String chain, const Biopolymer & inputBiopolymer,const  bool endCaps);
    ResidueID   residueID(map<const String,double> myUserVariables, const char* value , const String chain); // just like that below, except it can handle user-defined integer variables
    //ResidueID residueID(String inputResidueID, String chain); // this method of converting String to ResidueID has the advantage that it validates against the corresponding biopolymer.  Deprecated!  To be supplanted by the above.
    void constrainAllChainsToEachOther(ConstraintToGroundContainer & constraintToGroundContainer);
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
    void        replaceBiopolymerWithMutatedBiopolymerClass(BiopolymerClass & myOldBiopolymerClass, 
                                             String & myNewSequence, bool useNACappingHydroxyls = true);
    void        substituteResidue(String myChain , ResidueID myResidue, String mySubstitution, bool proteinCapping) ;
    void        insertResidue(Mutation myInsertion,  bool proteinCapping) ;
    vector <Mutation> &     updMutationVector(){return mutationVector;}
    const vector <Mutation> getMutationVector(){return mutationVector;}
    void        loadMutationVectorsFromSequence(); 
    void        loadMutationVectorsFromString(String myMutationString); 
    void        writeMutationFlexibilizers(std::ofstream & output, const int offset, const double radius);
    void        writeWaterDroplets(std::ofstream & output, const double springConstant, const double radius );
    void        writeMobilizerWithinMutation(std::ofstream & output,  const double radius  );
    const       vector<MMBAtomInfo> getConcatenatedAtomInfoVector(bool activeChainsOnly=false);
    const       vector<MMBAtomInfo> getConcatenatedAtomInfoVector(const State & state,bool activeChainsOnly=false);
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
    void        updateMutationResidueTypesFromCurrentSequence() ;
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

    void        renameChain(String oldChainID, String newChainID) {
        BiopolymerClass myBiopolymerClass = updBiopolymerClass(oldChainID); 
        myBiopolymerClass.setChainID(newChainID); myBiopolymerClass.myBiopolymer.setPdbChainId(newChainID); 
        // also have to delete and re-add to biopolymerClassMap with the new key
        deleteBiopolymerClass(oldChainID);
        addBiopolymerClass(newChainID,myBiopolymerClass);
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

