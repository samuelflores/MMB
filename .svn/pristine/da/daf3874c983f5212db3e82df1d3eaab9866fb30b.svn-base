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
#include <fstream>

int MonoAtoms::validate() {
	if (chainID.length() >1) {
        ErrorManager::instance <<__FILE__<<":"<<__LINE__<<": Please do not use chain ID’s longer than one character."<<endl; 
        ErrorManager::instance.treatError();
    }
	else if (chainID.length() < 1) {
        ErrorManager::instance <<__FILE__<<":"<<__LINE__<<": chain ID is less than one character long."<<endl; 
        ErrorManager::instance.treatError();
    }

    if (numAtoms <1) {
        ErrorManager::instance <<__FILE__<<":"<<__LINE__<<": You must specify more than 0 atoms."<<endl; 
        ErrorManager::instance.treatError();
    }
	if (atomName.size() >4) {
        ErrorManager::instance <<__FILE__<<":"<<__LINE__<<": Atom name must be <= 4 characters long."<<endl; 
        ErrorManager::instance.treatError();
    }
	else if (atomName.size() < 1) {
        ErrorManager::instance <<__FILE__<<":"<<__LINE__<<": Atom name must be at least 1 character long."<<endl; 
        ErrorManager::instance.treatError();
    }
    return 0;
}

MonoAtoms::MonoAtoms () {}

MonoAtoms::MonoAtoms (String myChainID,ResidueID myFirstResidueNumber, int myNumAtoms, String myAtomName) {
	chainID = myChainID;
	firstResidueID     = myFirstResidueNumber;
	numAtoms = myNumAtoms;	
	atomName = myAtomName;
	validate();

	for (int i = 0; i < myNumAtoms; i++) {
  	    if (atomName.compare("Mg+2") == 0) {
	        MagnesiumIon myMagnesiumIon;
                myMagnesiumIon.setPdbChainId(chainID);
		myMagnesiumIon.setPdbResidueNumber(i+getFirstResidueID().getResidueNumber());
		compoundVector.push_back(myMagnesiumIon);

  	    } else if (atomName.compare("Zn+2") == 0) {
	        ZincIon myIon;
                myIon.setPdbChainId(chainID);
		myIon.setPdbResidueNumber(i+getFirstResidueID().getResidueNumber());
		compoundVector.push_back(myIon);
  	    } else if (atomName.compare("Cl-") == 0) {
	        ChlorideIon myIon;
                myIon.setPdbChainId(chainID);
		myIon.setPdbResidueNumber(i+getFirstResidueID().getResidueNumber());
		compoundVector.push_back(myIon);

  	    } else if (atomName.compare("Na+") == 0) {
	        SodiumIon myIon;
                myIon.setPdbChainId(chainID);
		myIon.setPdbResidueNumber(i+getFirstResidueID().getResidueNumber());
		compoundVector.push_back(myIon);

  	    } else if (atomName.compare("K+") == 0) {
	        PotassiumIon myIon;
                myIon.setPdbChainId(chainID);
		myIon.setPdbResidueNumber(i+getFirstResidueID().getResidueNumber());
		compoundVector.push_back(myIon);

  	    } else if (atomName.compare("ZN+2") == 0) {
	        ZincIon myIon;
                myIon.setPdbChainId(chainID);
		myIon.setPdbResidueNumber(i+getFirstResidueID().getResidueNumber());
		compoundVector.push_back(myIon);

  	    } else if (atomName.compare("Li+") == 0) {
	        LithiumIon myIon;
                myIon.setPdbChainId(chainID);
		myIon.setPdbResidueNumber(i+getFirstResidueID().getResidueNumber());
		compoundVector.push_back(myIon);

  	    } else if (atomName.compare("Ca+2") == 0) {
	        CalciumIon myIon;
                myIon.setPdbChainId(chainID);
		myIon.setPdbResidueNumber(i+getFirstResidueID().getResidueNumber());
		compoundVector.push_back(myIon);

  	    } else if (atomName.compare("Cs+" ) == 0) {
	        CesiumIon myIon;
                myIon.setPdbChainId(chainID);
		myIon.setPdbResidueNumber(i+getFirstResidueID().getResidueNumber());
		compoundVector.push_back(myIon);

  	    } else if (atomName.compare("Rb+") == 0) {
	        RubidiumIon myIon;
                myIon.setPdbChainId(chainID);
		myIon.setPdbResidueNumber(i+getFirstResidueID().getResidueNumber());
		compoundVector.push_back(myIon);

	    } else {
 	        ErrorManager::instance <<__FILE__<<":"<<__LINE__<<": You have requested a monoAtoms of an unsupported type: "<<myAtomName<<". Currently only the following are supported:  Mg+2, Zn+2, Cl-, Na+, K+, Li+, Ca+2, Cs+, Rb+."<<endl << "Corresponding residue types should be MG, ZN, CL, NA, K, LI, CA, CS, RB."<<endl;;
		ErrorManager::instance.treatError();
	    }
        }
        renumberPdbResidues(myFirstResidueNumber);
}

const String MonoAtoms::getChainID() {
	return chainID;
} 

const ResidueID MonoAtoms::getFirstResidueID() {
        return firstResidueID;
	//return ResidueID(compoundVector[0].getPdbResidueNumber(),' ');// compoundVector[0].getPdbInsertionCode());
} 

const ResidueID MonoAtoms::getLastResidueID() {
	return ResidueID(compoundVector[compoundVector.size()-1].getPdbResidueNumber(),' ');
} 

const int MonoAtoms::getResidueIndex(ResidueID myResidueID) {
    for (int i = 0; i < getNumAtoms(); i++) {  
        if (myResidueID == ResidueID(compoundVector[i].getPdbResidueNumber(),' '))//compoundVector[i].getPdbInsertionCode()))
            return i;
    }
    ErrorManager::instance <<__FILE__<<":"<<__LINE__<<"Error!  Requested ResidueID : " <<myResidueID.outString()<<" is invalid."<<endl;
    ErrorManager::instance.treatError();
}

const int MonoAtoms::getNumAtoms() {
	return numAtoms;
} 

const String MonoAtoms::getAtomName() {
	return atomName;
} 

void MonoAtoms::validateResidue(  ResidueID residueNumber) {
	if (residueNumber < getFirstResidueID()) 
		{ErrorManager::instance <<__FILE__<<":"<<__LINE__<<"Error!  Requested a residue lower numberd than the first residue!"<<endl; ErrorManager::instance.treatError();}
	else if (residueNumber > getLastResidueID()) 
		{ErrorManager::instance <<__FILE__<<":"<<__LINE__<<"Error!  Requested a residue higher numberd than the last residue!"<<endl; ErrorManager::instance.treatError();}
	else { // do nothing
        };
}

 Compound MonoAtoms::getSingleCompound(ResidueID residueNumber){
	validateResidue (residueNumber ); 
	return compoundVector[residueNumber.getResidueNumber() - getFirstResidueID    ().getResidueNumber()];
} 


const bool MonoAtoms::hasAtom(ResidueID residueNumber)
{
	validateResidue (residueNumber ); 
        bool myHasAtom = compoundVector[residueNumber.getResidueNumber() - getFirstResidueID()  .getResidueNumber()].hasAtom(Compound::AtomPathName(atomName));
	return myHasAtom;
}

const Compound::AtomIndex MonoAtoms::getAtomIndex(ResidueID residueNumber)
{
	validateResidue (residueNumber ); 
	if (hasAtom(residueNumber) == false) 
	    {ErrorManager::instance <<__FILE__<<":"<<__LINE__<<"Error!  Attempted to get a   MonoAtom member that doesn't exist!"<<endl; ErrorManager::instance.treatError();}
        Compound::AtomIndex myAtomIndex = compoundVector[residueNumber.getResidueNumber() - getFirstResidueID()  .getResidueNumber()].getAtomIndex(Compound::AtomPathName(atomName));
	return myAtomIndex;
}

const Vec3 MonoAtoms::getAtomLocationInMobilizedBodyFrame(ResidueID residueNumber){
	validateResidue (residueNumber ); 
	// it should always be (0,0,0) .. unless someday we figure out how to rigidify the MonoAtoms.
        return compoundVector[residueNumber.getResidueNumber() - getFirstResidueID().getResidueNumber()].getAtomLocationInMobilizedBodyFrame(getAtomIndex(residueNumber));
	//atom index should always be zero, but just to be sure I'm getting it explicitly.  
}
void    MonoAtoms::adoptCompounds(SimTK::CompoundSystem & mySystem, bool readPreviousFrameFile){
	for (int i =  0 ; i<getNumAtoms(); i++) {
                //cout<<__FILE__<<":"<<__LINE__<< " X,Y,Z position of monoAtom: "<<0.1*sin((double)i)<<" , "<<0.1*cos((double)i)<<", "<<1.0*i/1<<endl;
                if (readPreviousFrameFile == 0)
		    mySystem.adoptCompound(compoundVector[i],Vec3(  0*sin((double)i) ,0*cos((double)i),1.0*i/1));
                else mySystem.adoptCompound(compoundVector[i]);
		//mySystem.adoptCompound(compoundVector[i],Vec3(  0.1*sin((double)i) ,0.1*cos((double)i),1.0*i/1));
	}
}



void MonoAtoms::matchDefaultConfiguration(PdbStructure pdbStructure,bool matchExact, bool matchIdealized)
{
	for (int i = 0; i < (int)compoundVector.size(); i++) {
		Compound::AtomTargetLocations atomTargets = compoundVector[i].createAtomTargets(pdbStructure);
                cout<<__FILE__<<":"<<__LINE__<<" contents of atomTargets.size() : "<< atomTargets.size()<<std::endl;
                //cout<<__FILE__<<":"<<__LINE__<<" contents of atomTargets : "<<std::endl;
		if (matchExact)
			{
			compoundVector[i].matchDefaultConfiguration(atomTargets,   Compound::Match_Exact );
			}
		if (matchIdealized) 
			{compoundVector[i].matchDefaultConfiguration(atomTargets,   Compound::Match_Idealized );} //planarity tolerance is in Radians, if Sherm's email is to be believed
                cout<<__FILE__<<":"<<__LINE__<<" position "<<compoundVector[i].calcDefaultAtomLocationInGroundFrame(compoundVector[i].getAtomName(Compound::AtomIndex(0)))<<endl;
	}
        
}



MobilizedBodyIndex MonoAtoms::getMobilizedBodyIndex(ResidueID residueNumber) {
	validateResidue (residueNumber ); 
	return compoundVector[residueNumber.getResidueNumber()-getFirstResidueID().getResidueNumber()].getAtomMobilizedBodyIndex(getAtomIndex(residueNumber));
 			 //atom index should always be zero, but just to be sure I'm getting it explicitly.  
}

MobilizedBody MonoAtoms::updMobilizedBody(ResidueID residueNumber, SimbodyMatterSubsystem & myMatter) {
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
        cout<<__FILE__<<":"<<__LINE__<<" About to get atom named : "<<getAtomName() <<endl;
        Compound::AtomIndex myAtomIndex = compoundVector[i].getAtomIndex(getAtomName() );
        DuMM::AtomIndex myDuMMAtomIndex = compoundVector[i].getDuMMAtomIndex(myAtomIndex);
        dumm.includeNonbondAtom(myDuMMAtomIndex);
    }
}

MonoAtomsContainer::MonoAtomsContainer() {}

const bool MonoAtomsContainer::hasChainID(String myChainID ) {
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
	else {ErrorManager::instance <<__FILE__<<":"<<__LINE__<<"Error!  Attempted to get a SingleAtom that doesn't exist!"<<endl; ErrorManager::instance.treatError();}
}

void MonoAtomsContainer::addMonoAtoms(MonoAtoms myMonoAtoms) {
		if (hasChainID(myMonoAtoms.getChainID()))
 			{ErrorManager::instance <<__FILE__<<":"<<__LINE__<<"Error!  A SingleAtom with that chainID already exists!"<<endl; ErrorManager::instance.treatError();}
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

void MonoAtoms::renumberPdbResidues (ResidueID firstResidueID ) {
    for (int i = 0; i < numAtoms; i++) {
        compoundVector[i].setPdbResidueNumber(firstResidueID.getResidueNumber()+i);
    }
}


void MonoAtoms::setPdbChainId (String chainID ) {
    for (int i = 0; i < numAtoms; i++) {
        compoundVector[i].setPdbChainId(chainID);
    }
}

void MonoAtoms::initialize (CompoundSystem & system,  bool readPreviousFrameFile,  String previousFrameFileName, bool matchExact, bool matchIdealized) {
    //renumberPdbResidues();
    setPdbChainId(chainID);    
    if (readPreviousFrameFile) { 
        std::ifstream inputFile(previousFrameFileName.c_str(), ifstream::in);
        PdbStructure pdbStructure(inputFile);
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

