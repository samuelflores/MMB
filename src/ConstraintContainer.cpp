/*
 *  ConstraintContainer.cpp
 *  MMB
 *
 *  Created by Samuel Flores on 12/13/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#include "ConstraintContainer.h"
#include "BiopolymerClass.h"


void ConstraintToGroundContainer::validateConstraintClass(const ConstraintClass & myConstraintClass, BiopolymerClassContainer & myBiopolymerClassContainer) {
	//cout << __FILE__<<":"<<__LINE__<<" Validating ConstraintClass : "<<endl;
	//myConstraintClass.print();
	if (! myBiopolymerClassContainer.hasChainID(myConstraintClass.getChain1())) 
    { 
        ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" Unable to find chain 1 : >"<<myConstraintClass.getChain1()<<"<"<<endl; 
        ErrorManager::instance.treatError();
    }
	if (myConstraintClass.getConstraintType() != WeldToGround) 
    {	
        // AT: I think this is an erroneous validation
		// if (myConstraintClass.getChain1().compare(myConstraintClass.getChain2()) != 0) {
		// 	ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" Chain 1: "<<myConstraintClass.getChain1()<<" is different from chain 2: "<<myConstraintClass.getChain2()<<endl; 
		// 	ErrorManager::instance.treatError();
		// }
		if (! myBiopolymerClassContainer.hasChainID(myConstraintClass.getChain2())) { 
			ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" Unable to find chain 2 : "<<myConstraintClass.getChain1()<<endl; 
			ErrorManager::instance.treatError();
		}
		if (!(myBiopolymerClassContainer.updBiopolymerClass(myConstraintClass.getChain2()).hasAtom(myConstraintClass.getResidueID2(), myConstraintClass.getAtomName2()))) {
			ErrorManager::instance << __FILE__<<":"<<__LINE__<<" Could not find chain " <<myConstraintClass.getChain2()<<" residue "<<myConstraintClass.getResidueID2().outString()<<", or maybe it has no atom named "<<myConstraintClass.getAtomName2()<<endl; 
			ErrorManager::instance.treatError();
		}
	}
	if (!(myBiopolymerClassContainer.updBiopolymerClass(myConstraintClass.getChain1()).hasAtom(myConstraintClass.getResidueID1(), myConstraintClass.getAtomName1()))) {
		ErrorManager::instance << __FILE__<<":"<<__LINE__<<" Could not find chain " <<myConstraintClass.getChain1()<<" residue "<<myConstraintClass.getResidueID1().outString()<<", or maybe it has no atom named "<<myConstraintClass.getAtomName1()<<endl; 
		ErrorManager::instance.treatError();
	}
	
	 
};

void ConstraintToGroundContainer::validateConstraintClassVector(BiopolymerClassContainer & myBiopolymerClassContainer){
    //printConstraintClasses();
	cout << __FILE__<<":"<<__LINE__<<" About to validate constraintClassVector"<<endl;
	for (int i = 0; i < constraintClassVector.size() ; i++) {
		validateConstraintClass( constraintClassVector[i], myBiopolymerClassContainer);
	}
};

void ConstraintToGroundContainer::pruneCoordinateCouplers(BiopolymerClassContainer & myBiopolymerClassContainer, DuMMForceFieldSubsystem & _dumm){
	/*
	for (int i = 0; i < numConstraintClasses(); i++) {
       		ConstraintClass myConstraintClass = getConstraintClass(i);          
		if ((!(_dumm.updRep().isBondAtom(getDuMMAtomIndex(myConstraintClass.getResidueID1(),myConstraintClass.getAtomName1())))) ||
		    (!(_dumm.updRep().isBondAtom(getDuMMAtomIndex(myConstraintClass.getResidueID2(),myConstraintClass.getAtomName2())))) )
                {
			
			i--;	
		} 	
	} */
}

void ConstraintToGroundContainer::printConstraintClasses() {
	for (int i = 0; i < numConstraintClasses(); i++) printConstraintClass(i);
};

void ConstraintToGroundContainer::applyConstrainChainRigidSegments (BiopolymerClassContainer & biopolymerClassContainer, CompoundSystem & system,  SimbodyMatterSubsystem & matter,State & state){
	for (int i = 0; i < constrainChainRigidSegmentsVector.size(); i++) {
                cout << __FILE__<<":"<<__LINE__<< " constraining rigid segments for chain : "<<constrainChainRigidSegmentsVector[i].chainID<<endl;
		biopolymerClassContainer.updBiopolymerClass(constrainChainRigidSegmentsVector[i].chainID).constrainRigidSegmentsToGround( system, matter, state,  *this , constrainChainRigidSegmentsVector[i].toGround, constrainChainRigidSegmentsVector[i].residueID);
	}
};

void ConstraintToGroundContainer::printConstraintClass(int constraintToGroundIndex) const {
    std::cout<<__FILE__<<":"<<__LINE__<<" Printing constraint with index = "<<constraintToGroundIndex<<std::endl;
    getConstraintClass(constraintToGroundIndex).print();
};

void ConstraintToGroundContainer::addConstraintClassToVector(ConstraintClass myConstraintClass){
    constraintClassVector.push_back (myConstraintClass); 
}

void ConstraintToGroundContainer::addConstraintClassToVector(String myChain, ResidueID myResidueID, String atomName) {
        std::cout<<__FILE__<<":"<<__LINE__<<" About to add constraintToGround for chain ID, ResidueID, and atomName: "<<myChain<<", "<<myResidueID.outString()<<", "<<atomName<<std::endl;
        ConstraintClass myConstraintClass(myChain, myResidueID,atomName); 
        addConstraintClassToVector(myConstraintClass); 
}


void ConstraintToGroundContainer::addConstraintToVector(String myChain, ResidueID myResidueID, String atomName,
                                                        String myChain2, ResidueID myResidueID2, String atomName2
                                                        ){
    ConstraintClass myConstraintClass(myChain, myResidueID,atomName, myChain2, myResidueID2, atomName2, WeldToAtom); 
    addConstraintClassToVector(myConstraintClass); 
}

bool ConstraintToGroundContainer::hasConstraintClass(String myChainID, ResidueID myResidueID) {
    for (int i = 0; i < numConstraintClasses(); i++) {
        if ((getConstraintClass(i).getChain1().compare(myChainID) == 0) &&
            (getConstraintClass(i).getResidueID1() == myResidueID)) {
            return bool(true);
        }
    }
    return bool(false);
}


void ConstraintToGroundContainer::deleteConstraintClass(int index){
    constraintClassVector.erase(constraintClassVector.begin()+index);
}

void ConstraintToGroundContainer::updateConstraintToVector(int index,
                                                           String myChain, ResidueID myResidueID, 
                                                           String atomName,
                                                           BiopolymerClassContainer& biopolymerClassContainer){
    if(index < 0 || index >= constraintClassVector.size())
    {
        ErrorManager::instance << __FILE__<<":"<<__LINE__<<" Could not find constraint with id " << index <<endl; 
        ErrorManager::instance.treatError();
    }
    if(atomName == "") 
        atomName = biopolymerClassContainer.updBiopolymerClass(myChain).getRepresentativeAtomName();
    ConstraintClass newConstraint = ConstraintClass(myChain, myResidueID, atomName);
    validateConstraintClass(newConstraint, biopolymerClassContainer);
    constraintClassVector[index] = newConstraint;
}

void ConstraintToGroundContainer::updateConstraintToVector(int index,
                                                           String myChain, ResidueID myResidueID, 
                                                           String atomName,
                                                           String myChain2, ResidueID myResidueID2, 
                                                           String atomName2,
                                                           BiopolymerClassContainer& biopolymerClassContainer){
    if(index < 0 || index >= constraintClassVector.size())
    {
        ErrorManager::instance << __FILE__<<":"<<__LINE__<<" Could not find constraint with id " << index <<endl; 
        ErrorManager::instance.treatError();
    }
    if(atomName == "") atomName = biopolymerClassContainer.updBiopolymerClass(myChain).getRepresentativeAtomName();
    if(atomName2 == "") atomName2 = biopolymerClassContainer.updBiopolymerClass(myChain2).getRepresentativeAtomName();
    ConstraintClass newConstraint = ConstraintClass(myChain, myResidueID, atomName,
                                 myChain2, myResidueID2, atomName2, WeldToAtom);
    validateConstraintClass(newConstraint, biopolymerClassContainer);
    constraintClassVector[index] = newConstraint;
}






