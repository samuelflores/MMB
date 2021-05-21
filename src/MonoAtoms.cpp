/* -------------------------------------------------------------------------- *
 *                           MMB (MacroMoleculeBuilder)                       *
 * -------------------------------------------------------------------------- *
 *                                                                            *
 * Copyright (c) 2011-12 by the Author.                                       *
 * Author: Samuel Flores                                                      *
 *                                                                            *
 * See RNABuilder.cpp for the copyright and usage agreement.                  *
 * -------------------------------------------------------------------------- */

#include "Utils.h"
#include "MonoAtoms.h"
#include "molmodel/internal/Ions.h"
#include <cstdlib>
#include <fstream>

int MonoAtoms::validate() {
	if (chainID.length() >1) {
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "Please do not use chain ID’s longer than one character."<<endl);
    }
	else if (chainID.length() < 1) {
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "chain ID is less than one character long."<<endl);
    }

    //if (numAtoms <1) {
    //    MMBLOG_FILE_FUNC_LINE(CRITICAL, "You must specify more than 0 atoms."<<endl);
    //}
	if (atomName.size() >4) {
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "Atom name must be <= 4 characters long."<<endl);
    }
	else if (atomName.size() < 1) {
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "Atom name must be at least 1 character long."<<endl);
    }
    return 0;
}

MonoAtoms::MonoAtoms () {}

static
Molecule returnIonOfCorrectType(const String &atomName) {
    if (atomName.compare("Mg+2") == 0) {
        return MagnesiumIon{};
    } else if (atomName.compare("Zn+2") == 0 || atomName.compare("ZN+2") == 0) {
        return ZincIon{};
    } else if (atomName.compare("Cl-") == 0) {
        return ChlorideIon{};
    } else if (atomName.compare("Na+") == 0) {
        return SodiumIon{};
    } else if (atomName.compare("K+") == 0) {
        return PotassiumIon{};
    } else if (atomName.compare("Li+") == 0) {
        return LithiumIon{};
    } else if (atomName.compare("Ca+2") == 0) {
        return CalciumIon{};
    } else if (atomName.compare("Cs+" ) == 0) {
        return CesiumIon{};
    } else if (atomName.compare("Rb+") == 0) {
        return RubidiumIon{};
    }

    MMBLOG_FILE_FUNC_LINE(CRITICAL, "You have requested a monoAtoms of an unsupported type: "<<  atomName<<". Currently only the following are supported:  Mg+2, Zn+2, Cl-, Na+, K+, Li+, Ca+2, Cs+, Rb+."<<endl << "Corresponding residue types should be MG, ZN, CL, NA, K, LI, CA, CS, RB."<<endl);
    std::abort();
}



// This method figures out the atom name, chainID, and residue number of the  latest/ highest ion, increments the residue number, adds one ion with that name and the incremented res number to the compound vector, and increments numAtoms.
// Don't forget to renumber residue numbers later!
void MonoAtoms::addMonoAtom ( ) {
    validate();
    Molecule myIon = returnIonOfCorrectType(atomName);
    myIon.setPdbChainId(chainID);
    compoundVector.push_back(myIon);
}
void MonoAtoms::addMonoAtom (Vec3 positionVec3 ) {
    addMonoAtom();
    compoundVector.back().setTopLevelTransform(Transform(positionVec3))    ;
}

MonoAtoms::MonoAtoms (String myChainID,ResidueID myFirstResidueNumber, int   numAtoms, String myAtomName) {
	compoundVector.clear();
	chainID = myChainID;
	firstResidueID     = myFirstResidueNumber;
	//numAtoms = myNumAtoms;	
	atomName = myAtomName;
	validate();
        //Molecule myMolecule;
	for (int i = 0; i < getNumAtoms(); i++) {
            addMonoAtom();		
            /*myMolecule = returnIonOfCorrectType(atomName);
            myMolecule.setPdbChainId(chainID);
            myMolecule.setPdbResidueNumber(i+getFirstResidueID().getResidueNumber());
            compoundVector.push_back(myMolecule);*/
        }
        renumberPdbResidues();
}

String MonoAtoms::getChainID() {
    return chainID;
} 

ResidueID MonoAtoms::getFirstResidueID() {
    return firstResidueID;
    //return ResidueID(compoundVector[0].getPdbResidueNumber(),' ');// compoundVector[0].getPdbInsertionCode());
} 

ResidueID MonoAtoms::getResidueID(const int myResidueIndex) {
    return ResidueID(compoundVector[myResidueIndex         ].getPdbResidueNumber(),' ');
} 

ResidueID MonoAtoms::getLastResidueID() {
	return ResidueID(compoundVector[compoundVector.size()-1].getPdbResidueNumber(),' ');
} 

int MonoAtoms::getResidueIndex(ResidueID myResidueID) {
    for (int i = 0; i < getNumAtoms(); i++) {  
        if (myResidueID == ResidueID(compoundVector[i].getPdbResidueNumber(),' '))//compoundVector[i].getPdbInsertionCode()))
            return i;
    }
    MMBLOG_FILE_FUNC_LINE(CRITICAL, "Requested ResidueID : " <<myResidueID.outString()<<" is invalid."<<endl);
}

int MonoAtoms::getNumAtoms() {
	return compoundVector.size();
	//return numAtoms;
} 

String MonoAtoms::getAtomName() {
	return atomName;
} 

void MonoAtoms::validateResidue(  ResidueID residueNumber) {
	if (residueNumber < getFirstResidueID()) 
		{MMBLOG_FILE_FUNC_LINE(CRITICAL, "Requested a residue lower numberd than the first residue!"<<endl);}
	else if (residueNumber > getLastResidueID()) 
		{MMBLOG_FILE_FUNC_LINE(CRITICAL, "Requested a residue higher numberd than the last residue!"<<endl);}
	else { // do nothing
        };
}

Compound MonoAtoms::getSingleCompound(ResidueID residueNumber){
	validateResidue (residueNumber ); 
	return compoundVector[residueNumber.getResidueNumber() - getFirstResidueID    ().getResidueNumber()];
} 


bool MonoAtoms::hasAtom(ResidueID residueNumber)
{
	validateResidue (residueNumber ); 
        bool myHasAtom = compoundVector[residueNumber.getResidueNumber() - getFirstResidueID()  .getResidueNumber()].hasAtom(Compound::AtomPathName(atomName));
	return myHasAtom;
}

Compound::AtomIndex MonoAtoms::getAtomIndex(ResidueID residueNumber)
{
	validateResidue (residueNumber ); 
	if (hasAtom(residueNumber) == false) 
	    {MMBLOG_FILE_FUNC_LINE(CRITICAL, "Attempted to get a   MonoAtom member that doesn't exist!"<<endl); }
        Compound::AtomIndex myAtomIndex = compoundVector[residueNumber.getResidueNumber() - getFirstResidueID()  .getResidueNumber()].getAtomIndex(Compound::AtomPathName(atomName));
	return myAtomIndex;
}
Vec3 MonoAtoms::getAtomLocationInGroundFrame(const int       residueIndex, const State & state ){
	ResidueID residueNumber = getResidueID(residueIndex);
	validateResidue (residueNumber ); 
	//if (getAtomIndex(residueNumber))
	//    {MMBLOG_FILE_FUNC_LINE(CRITICAL, " Expected an atom index of zero, got:"<<getAtomIndex(residueNumber)<<endl); }
	//atom index should always be zero, but just to be sure I'm getting it explicitly.  
        //return compoundVector[residueIndex   ].calcDefaultAtomLocationInGroundFrame(getAtomIndex(residueNumber));
	MMBLOG_FILE_FUNC_LINE(DEBUG, " About to issue  compoundVector["<<residueIndex<<"].calcAtomLocationInGroundFrame(state,0)"<<endl<<" Got : "<<compoundVector[residueIndex   ].calcAtomLocationInGroundFrame(state,SimTK::Compound::AtomIndex(0))<<endl); 
        return compoundVector[residueIndex   ].calcAtomLocationInGroundFrame(state,SimTK::Compound::AtomIndex(0));
}

Vec3 MonoAtoms::getAtomLocationInMobilizedBodyFrame(ResidueID residueNumber){
	validateResidue (residueNumber ); 
	// it should always be (0,0,0) .. unless someday we figure out how to rigidify the MonoAtoms.
        return compoundVector[residueNumber.getResidueNumber() - getFirstResidueID().getResidueNumber()].getAtomLocationInMobilizedBodyFrame(getAtomIndex(residueNumber));
	//atom index should always be zero, but just to be sure I'm getting it explicitly.  
}
void    MonoAtoms::adoptCompounds(SimTK::CompoundSystem & mySystem, bool readPreviousFrameFile){
	for (int i =  0 ; i<getNumAtoms(); i++) {
                //MMBLOG_FILE_FUNC_LINE( " X,Y,Z position of monoAtom: "<<0.1*sin((double)i)<<" , "<<0.1*cos((double)i)<<", "<<1.0*i/1<<endl;
                /*if (readPreviousFrameFile == 0) {
		    mySystem.adoptCompound(compoundVector[i],Vec3(  0*sin((double)i) ,0*cos((double)i),1.0*i/1));
		}else {*/
		    // 21.05.04 SCF We previously called adoptCompound differently for readPreviousFrameFile == 0. However this was interfering with our ability to generate spherical helices of ions. So now monoatoms are always adopted the same way at this stage. We should create additional commands to create monoAtoms patterns, if future users want to adopt ions in other geometries.
                    mySystem.adoptCompound(compoundVector[i]);
		/*}*/
		//mySystem.adoptCompound(compoundVector[i],Vec3(  0.1*sin((double)i) ,0.1*cos((double)i),1.0*i/1));
	}
}



void MonoAtoms::matchDefaultConfiguration(PdbStructure pdbStructure,bool matchExact, bool matchIdealized)
{
	for (int i = 0; i < (int)compoundVector.size(); i++) {
		Compound::AtomTargetLocations atomTargets = compoundVector[i].createAtomTargets(pdbStructure);
                MMBLOG_FILE_FUNC_LINE(INFO, "contents of atomTargets.size() : "<< atomTargets.size()<<endl);
                //MMBLOG_FILE_FUNC_LINE(" contents of atomTargets : "<<std::endl;
		if (matchExact)
			{
			compoundVector[i].matchDefaultConfiguration(atomTargets,   Compound::Match_Exact );
			}
		if (matchIdealized) 
			{compoundVector[i].matchDefaultConfiguration(atomTargets,   Compound::Match_Idealized );} //planarity tolerance is in Radians, if Sherm's email is to be believed
                MMBLOG_FILE_FUNC_LINE(INFO, "position "<<compoundVector[i].calcDefaultAtomLocationInGroundFrame(compoundVector[i].getAtomName(Compound::AtomIndex(0)))<<endl);
	}
        
}



MobilizedBodyIndex MonoAtoms::getMobilizedBodyIndex(ResidueID residueNumber) {
	validateResidue (residueNumber ); 
	return compoundVector[residueNumber.getResidueNumber()-getFirstResidueID().getResidueNumber()].getAtomMobilizedBodyIndex(getAtomIndex(residueNumber));
 			 //atom index should always be zero, but just to be sure I'm getting it explicitly.  
}

MobilizedBody & MonoAtoms::updMobilizedBody(ResidueID residueNumber, SimbodyMatterSubsystem & myMatter) {
	validateResidue (residueNumber ); 
	MobilizedBodyIndex myMobilizedBodyIndex = getMobilizedBodyIndex(residueNumber);
	return myMatter.updMobilizedBody(myMobilizedBodyIndex);
}

String MonoAtoms::getAtomPathName(ResidueID myResidue) {
	validateResidue(myResidue);
	return String (
		intToString(myResidue.getResidueNumber()) +
                "/" +
		atomName);
}

// here we add all monoAtoms to the physics zone

void MonoAtoms::includeAllAtoms(DuMMForceFieldSubsystem & dumm){
    for (int i = 0; i < (int)compoundVector.size(); i++) {
        MMBLOG_FILE_FUNC_LINE(INFO, "About to get atom named : "<<getAtomName() <<endl);
        Compound::AtomIndex myAtomIndex = compoundVector[i].getAtomIndex(getAtomName() );
        DuMM::AtomIndex myDuMMAtomIndex = compoundVector[i].getDuMMAtomIndex(myAtomIndex);
        dumm.includeNonbondAtom(myDuMMAtomIndex);
    }
}

double dotProduct (const Vec3 vec1 , const Vec3 vec2){
    return (vec1[0]*vec2[0] +vec1[1]*vec2[1] + vec1[2]*vec2[2]);
}
double magnitude  (const Vec3 vec){
    return sqrt(vec[0]*vec[0] +vec[1]*vec[1] + vec[2]*vec[2]);
}

double MonoAtoms::computeCurvatureSquared(const int index, const State & state){
    if (index < 1){
        MMBLOG_FILE_FUNC_LINE(CRITICAL, " Cannot compute curvature for the first or last atom.  "<<endl);
    }
    else if (index > (getNumAtoms()-2)){
        MMBLOG_FILE_FUNC_LINE(CRITICAL, " Cannot compute curvature for the first or last atom.  "<<endl);
    }
    Vec3 priorVec = ( getAtomLocationInGroundFrame(index  ,state) -  getAtomLocationInGroundFrame(index-1,state) );
    Vec3 nextVec = ( getAtomLocationInGroundFrame(index+1,state) -  getAtomLocationInGroundFrame(index  ,state) );
    double angleBetweenVectors = acos(dotProduct(priorVec,nextVec)/magnitude(priorVec)/ magnitude(nextVec));
    MMBLOG_FILE_FUNC_LINE(DEBUG, " For index "<<index<<" computed angle = "<< angleBetweenVectors <<" and squared angle of "<<angleBetweenVectors*angleBetweenVectors <<endl);
    return angleBetweenVectors*angleBetweenVectors;
}

double MonoAtoms::computeTotalCurvatureSquared( const State & state){
    double totalCurvature = 0.;
    for (int i = 1; i < (getNumAtoms()-1); i++){
       totalCurvature += computeCurvatureSquared(i,state);
    }
    if (getNumAtoms() <3){
        MMBLOG_FILE_FUNC_LINE(CRITICAL, " Cannot compute curvature for fewer than 3 atoms.  "<<endl);
    }
    return totalCurvature;
}


MonoAtomsContainer::MonoAtomsContainer() {}

bool MonoAtomsContainer::hasChainID(String myChainID ) {
	if (monoAtomsMap.find(myChainID) == monoAtomsMap.end())
	 	{return false;}
	else
		{return true;}
}

MonoAtoms MonoAtomsContainer::getMonoAtoms(String myChainID) {
	if (hasChainID(myChainID))
		{MonoAtoms myMonoAtoms ( monoAtomsMap[myChainID]);
		return myMonoAtoms;	
                }
	else {MMBLOG_FILE_FUNC_LINE(CRITICAL, "Attempted to get a SingleAtom that doesn't exist!"<<endl);}
}

void MonoAtomsContainer::addMonoAtoms(MonoAtoms myMonoAtoms) {
		if (hasChainID(myMonoAtoms.getChainID()))
 			{MMBLOG_FILE_FUNC_LINE(CRITICAL, "A SingleAtom with that chainID already exists!"<<endl);}
		else {
			monoAtomsMap[myMonoAtoms.getChainID()] = myMonoAtoms;
		}
}	

void    MonoAtomsContainer::adoptCompounds(CompoundSystem & mySystem, bool readPreviousFrameFile){
        map <String,MonoAtoms> :: iterator monoAtomsMapIterator;
	if (monoAtomsMap.size() > 0)
		for (monoAtomsMapIterator = monoAtomsMap.begin();
       			monoAtomsMapIterator != monoAtomsMap.end();
			monoAtomsMapIterator++) 
		{
			(*monoAtomsMapIterator).second.adoptCompounds(mySystem,readPreviousFrameFile);
		}
}

void MonoAtomsContainer::matchDefaultConfiguration(SimTK::PdbStructure pdbStructure,bool matchExact, bool matchIdealized)
{
        map <String,MonoAtoms> :: iterator monoAtomsMapIterator;
	if (monoAtomsMap.size() > 0)
		for (monoAtomsMapIterator = monoAtomsMap.begin();
       			monoAtomsMapIterator != monoAtomsMap.end();
			monoAtomsMapIterator++) 
		{
			(*monoAtomsMapIterator).second.matchDefaultConfiguration(pdbStructure,matchExact, matchIdealized);
		}
}

void MonoAtoms::renumberPdbResidues ( ) {
    for (int i = 0; i < getNumAtoms(); i++) {
        compoundVector[i].setPdbResidueNumber(firstResidueID.getResidueNumber()+i);
    }
}


void MonoAtoms::setPdbChainId (String chainID ) {
    for (int i = 0; i < getNumAtoms(); i++) {
        compoundVector[i].setPdbChainId(chainID);
    }
}

void MonoAtoms::initialize (CompoundSystem & system,  bool readPreviousFrameFile,  String previousFrameFileName, bool matchExact, bool matchIdealized) {
    setPdbChainId(chainID);    
    if (readPreviousFrameFile) { 
        PdbStructure pdbStructure{previousFrameFileName};
        matchDefaultConfiguration(  pdbStructure, matchExact, matchIdealized);
    }
    adoptCompounds(system,readPreviousFrameFile );
}


void MonoAtomsContainer::initialize (CompoundSystem & system , bool readPreviousFrameFile, String previousFrameFileName, bool matchExact, bool matchIdealized ) {
    
        map <String,MonoAtoms> :: iterator monoAtomsMapIterator;
	if (monoAtomsMap.size() > 0)
		for (monoAtomsMapIterator = monoAtomsMap.begin();
       			monoAtomsMapIterator != monoAtomsMap.end();
			monoAtomsMapIterator++) 
		{
			(*monoAtomsMapIterator).second.initialize(system, readPreviousFrameFile,  previousFrameFileName,  matchExact, matchIdealized);
		}
}
void MonoAtomsContainer::remove (String myChainID){
    monoAtomsMap.erase(myChainID); // delete monoAtom element with key = myChainID	
}
String MonoAtomsContainer::getAtomPathName(String myChain,ResidueID myResidue) {
    return monoAtomsMap[myChain].getAtomPathName(myResidue);
}

void MonoAtomsContainer::clear() {
	monoAtomsMap.clear();
}

void MonoAtomsContainer::includeAllAtoms(DuMMForceFieldSubsystem & dumm){


        map <String,MonoAtoms> :: iterator monoAtomsMapIterator;
	if (monoAtomsMap.size() > 0)
		for (monoAtomsMapIterator = monoAtomsMap.begin();
       			monoAtomsMapIterator != monoAtomsMap.end();
			monoAtomsMapIterator++) 
		{
			(*monoAtomsMapIterator).second.includeAllAtoms(dumm);
		}

}

double MonoAtomsContainer::computeTotalCurvatureSquared(const State & state){
        map <String,MonoAtoms> :: iterator monoAtomsMapIterator;
	double totalCurvature = 0.;
	if (monoAtomsMap.size() > 0)
		for (monoAtomsMapIterator = monoAtomsMap.begin();
       			monoAtomsMapIterator != monoAtomsMap.end();
			monoAtomsMapIterator++) 
		{
			// Don't compute if there are fewer than 3 atoms. 
			if ((*monoAtomsMapIterator).second.getNumAtoms()>2) totalCurvature += (*monoAtomsMapIterator).second.computeTotalCurvatureSquared(state);
		}
       MMBLOG_FILE_FUNC_LINE(INFO, " Computed total curvature squared = "<< totalCurvature  <<endl);
       return totalCurvature;

          
}

