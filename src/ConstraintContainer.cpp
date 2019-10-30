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
    cout<<__FILE__<<":"<<__LINE__<<endl;
    ConstraintClass myConstraintClass(myChain, myResidueID,atomName, myChain2, myResidueID2, atomName2, WeldToAtom); 
    cout<<__FILE__<<":"<<__LINE__<<endl;
    addConstraintClassToVector(myConstraintClass); 
    cout<<__FILE__<<":"<<__LINE__<<endl;
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


// Determine whether there is any constraint at all defined between two given chain ID's
bool ConstraintToGroundContainer::hasConstraintClass(String myChainID1, String myChainID2) {
    //std::cout<<__FILE__<<":"<<__LINE__<<" Checking whether any constraints at all exist between chains "<<myChainID1<<" and "<<myChainID2<<std::endl;
    for (int i = 0; i < numConstraintClasses(); i++) {
        if ((getConstraintClass(i).getChain1().compare(myChainID1) == 0) &&
            (getConstraintClass(i).getChain2().compare(myChainID2) == 0)) {
            //std::cout<<__FILE__<<":"<<__LINE__<<" TRUE. At least one such constraint exists. It is:"<<std::endl;
            //getConstraintClass(i).print();
            return bool(true);
        }
    }
    std::cout<<__FILE__<<":"<<__LINE__<<" FALSE. No constraint found between chains "<<myChainID1<<" and "<<myChainID2<<std::endl;
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

/*
void ConstraintToGroundContainer::addSingleWeldConstraintPerInterfaceChainPair(   BiopolymerClassContainer & myBiopolymerClassContainer) { // This polymorphism requires that the user specify two sets of chains.  Only residues at the interface between the two sets will be included.  This lets the user leave out other chains (e.g. threading templates) which are in the system but which shouldn't be flexibilized.
        OpenMM::NeighborList neighborList;
        OpenMM::Vec3 boxSize = OpenMM::Vec3(10000,10000,10000);
        vector<MMBAtomInfo> concatenatedAtomInfoVector = myBiopolymerClassContainer.getConcatenatedAtomInfoVector();
        vector<OpenMM::Vec3> particleList(concatenatedAtomInfoVector.size());
        vector<set<int> > exclusions( particleList.size() );
        for (int i = 0; i < concatenatedAtomInfoVector.size() ; i++) {
            particleList[i] = concatenatedAtomInfoVector[i].position;
        }

        cout<<__FILE__<<":"<<__LINE__<<" neighborList size is : "<<neighborList.size()<<endl;
        for (int h = 0 ; h < interfaceContainer.numInterfaces(); h++ ){ // loop through interfaceContainer interfaces ..
            vector<String> referenceChains = interfaceContainer.getInterface(h).getChains();  
            vector<String> partnerChains = interfaceContainer.getInterface(h).getPartnerChains();  
            double         radius        = interfaceContainer.getInterface(h).getDepth();  
            cout<<__FILE__<<":"<<__LINE__<<"Now turning interface "<< h << " to individual constraints between pairs of atoms."<<endl;
            interfaceContainer.getInterface(h).print(); 
            cout<<__FILE__<<":"<<__LINE__<<endl;
            computeNeighborListVoxelHash(neighborList, particleList.size() , particleList, exclusions, &boxSize, false, radius  , 0.0);
            for ( int j = 0 ; j < neighborList.size(); j++) {
                
                if (((( vectorCompare(concatenatedAtomInfoVector[neighborList[j].first].chain , (referenceChains))) == 1) &&
                     (( vectorCompare(concatenatedAtomInfoVector[neighborList[j].second].chain ,(  partnerChains))) == 1))  != //Use an XOR here. This means if the 'partnerChains' evaluation is later set to return 1 when partnerChains is empty, this will still work.
                    ((( vectorCompare(concatenatedAtomInfoVector[neighborList[j].second].chain ,(referenceChains))) == 1) &&
                     (( vectorCompare(concatenatedAtomInfoVector[neighborList[j].first].chain  ,(  partnerChains))) == 1))     //Make sure that exactly one residue is in the 'referenceChains', and the other residue is in the 'partnerChains' .. thus only the desired interface is included
                                                                                                                          )

                {
                    ResidueID residueID1(concatenatedAtomInfoVector[neighborList[j].first].residueID);
                    String chain1(concatenatedAtomInfoVector[neighborList[j].first].chain);
                    String atom1(concatenatedAtomInfoVector[neighborList[j].first].atomName);
                    ResidueID residueID2(concatenatedAtomInfoVector[neighborList[j].second].residueID);
                    String chain2(concatenatedAtomInfoVector[neighborList[j].second].chain);
                    String atom2(concatenatedAtomInfoVector[neighborList[j].second].atomName);
                    if (!(hasConstraintClass(chain1,chain2))) {
                    
                        ConstraintClass myConstraintClass(chain1 ,residueID1,atom1,chain2, residueID2,atom2, WeldToAtom);
                        addConstraintClassToVector(myConstraintClass);
                        cout<<__FILE__<<":"<<__LINE__<<" CONFIRMED that there is not constraint between these two chains. Added ConstraintClass :"<<endl;
                        myConstraintClass.print();
                        cout<<__FILE__<<":"<<__LINE__<<endl;
                        
                        // Right here, should consider deleting particleList  element. But first, are we sure this is kosher for users who go from large Depth to small Depth..? .. in any case, we don't have the index with which to delete from particleList!
                    }
                }
                else {
                    // Right here, should probably delete particleList element.
                }
            }
        }    
    };
*/


#ifdef USE_OPENMM
void ConstraintToGroundContainer::addSingleWeldConstraintPerInterfaceChainPair(   BiopolymerClassContainer & myBiopolymerClassContainer) { // This polymorphism requires that the user specify two sets of chains.  Only residues at the interface between the two sets will be included.  This lets the user leave out other chains (e.g. threading templates) which are in the system but which shouldn't be flexibilized.
        vector<MMBAtomInfo> concatenatedAtomInfoVector = myBiopolymerClassContainer.getConcatenatedAtomInfoVector();
        vector<TwoAtomClass> myTwoAtomClassVector = interfaceContainer.retrieveCloseContactPairs(concatenatedAtomInfoVector);
        for (int i = 0; i < myTwoAtomClassVector.size(); i++){
                    if (!(hasConstraintClass(myTwoAtomClassVector[i]. getChain1(),myTwoAtomClassVector[i].getChain2()))) {
                    
                        ConstraintClass myConstraintClass(myTwoAtomClassVector[i]. getChain1() ,myTwoAtomClassVector[i]. getResidueID1(),myTwoAtomClassVector[i]. getAtomName1(),
                                                          myTwoAtomClassVector[i]. getChain2() ,myTwoAtomClassVector[i]. getResidueID2(),myTwoAtomClassVector[i]. getAtomName2(), WeldToAtom);    

                        addConstraintClassToVector(myConstraintClass);
                        cout<<__FILE__<<":"<<__LINE__<<" CONFIRMED that there is no constraint between these two chains. Added ConstraintClass :"<<endl;
                        myConstraintClass.print();
                        cout<<__FILE__<<":"<<__LINE__<<endl;
                        
                    }
        }
};
#endif



