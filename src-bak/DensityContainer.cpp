/* -------------------------------------------------------------------------- *
 *                           MMB (MacroMoleculeBuilder)                       *
 * -------------------------------------------------------------------------- *
 *                                                                            *
 * Copyright (c) 2011-12 by the Author.                                       *
 * Author: Samuel Flores                                                      *
 *                                                                            *
 * See RNABuilder.cpp for the copyright and usage agreement.                  *
 * -------------------------------------------------------------------------- */

#include "DensityContainer.h"

void DensityContainer::clear(){
    densityStretchVector.clear();
}
void DensityContainer::validate(const DensityStretch & myDensityStretch,BiopolymerClassContainer & myBiopolymerClassContainer){

    /*if (!(myDensityStretch.getStartResidue() == myBiopolymerClassContainer.updBiopolymerClass(myDensityStretch.getChain()).getFirstResidueID ()) ) {
        ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" The start residue (currently "<<myDensityStretch.getStartResidue().outString()<<") is not equal to the first residue number of the chain ("<<myBiopolymerClassContainer.updBiopolymerClass(myDensityStretch.getChain()).getFirstResidueID().outString()<< "). "<<endl; ErrorManager::instance.treatError();    }
    if (!(myDensityStretch.getEndResidue() == myBiopolymerClassContainer.updBiopolymerClass(myDensityStretch.getChain()).getLastResidueID    ()) ) {
        ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" The start residue (currently "<<myDensityStretch.getEndResidue().outString()<<") is not equal to the first residue number of the chain ("<<myBiopolymerClassContainer.updBiopolymerClass(myDensityStretch.getChain()).getLastResidueID    ()    .outString()     << "). "<<endl; ErrorManager::instance.treatError();    } */

    if (myBiopolymerClassContainer.updBiopolymerClass(myDensityStretch.getChain()).difference(myDensityStretch.getEndResidue() , myDensityStretch.getStartResidue()) < 0) {
        ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" The end residue (currently "<<myDensityStretch.getEndResidue().outString()<<") must be greater than or equal to the start residue (currently "<<myDensityStretch.getStartResidue().outString()<<". "<<endl; ErrorManager::instance.treatError();    }
    if ((myDensityStretch.getEndResidue() > myBiopolymerClassContainer.updBiopolymerClass(myDensityStretch.getChain()).getLastResidueID    ()) ) {
        ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" The end residue (currently "<<myDensityStretch.getEndResidue().outString()<<") is greater than the last residue number of the chain."<<endl; ErrorManager::instance.treatError();    }
    if ((myDensityStretch.getStartResidue() < myBiopolymerClassContainer.updBiopolymerClass(myDensityStretch.getChain()).getFirstResidueID    ()) ) {
        ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" The start residue (currently "<<myDensityStretch.getStartResidue().outString()<<") is lesser than the first residue number of the chain."<<endl; ErrorManager::instance.treatError();    }
    if (!(myBiopolymerClassContainer.hasChainID(myDensityStretch.getChain()))){
        ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" Couldn't find chain "<<myDensityStretch.getChain()<<endl; ErrorManager::instance.treatError();
    }   
}

void	DensityContainer::add(const DensityStretch & myDensityStretch,  BiopolymerClassContainer & myBiopolymerClassContainer){
	validate(myDensityStretch,  myBiopolymerClassContainer);
	densityStretchVector.push_back(myDensityStretch);
}

void DensityContainer::updateDensityStretch(int id, const DensityStretch & stretch, BiopolymerClassContainer & myBiopolymerClassContainer){
    if(id < 0 || id >= densityStretchVector.size()){
        ErrorManager::instance <<__FILE__<<":"<<__LINE__<<": you tried to delete a non existing Contact." << endl;
        ErrorManager::instance.treatError();
    }
    validate(stretch, myBiopolymerClassContainer);
    densityStretchVector[id] = stretch;
}

void DensityContainer::deleteDensityStretch(int id){
    if(id < 0 || id >= densityStretchVector.size()){
        ErrorManager::instance <<__FILE__<<":"<<__LINE__<<": you tried to delete a non existing Contact." << endl;
        ErrorManager::instance.treatError();
    }
    densityStretchVector.erase(densityStretchVector.begin()+id);    
}

DensityStretch DensityContainer::getDensityStretch(int densityStretchIndex){
	return densityStretchVector[densityStretchIndex];
}

int DensityContainer::numDensityStretches(){
	return densityStretchVector.size();
}

void DensityContainer::stuffDensityStretchVector( BiopolymerClassContainer & myBiopolymerClassContainer){
	if (numDensityStretches() == 0) {
        	cout<<__FILE__<<":"<<__LINE__<<" All residues in all available biopolymer chains be fitted to the map. "<<endl; 
		for (int i = 0 ; i < myBiopolymerClassContainer.getNumBiopolymers(); i++){
			DensityStretch myDensityStretch;
			myDensityStretch.setChain( myBiopolymerClassContainer.updBiopolymerClass(i).getChainID());
			myDensityStretch.setStartResidue( myBiopolymerClassContainer.updBiopolymerClass(i).getFirstResidueID    ());
			myDensityStretch.setEndResidue( myBiopolymerClassContainer.updBiopolymerClass(i).getLastResidueID    ());
			add(myDensityStretch,myBiopolymerClassContainer);
		}
	} else {
        	ErrorManager::instance <<__FILE__<<":"<<__LINE__<<" The stuffDensityStretchVector can only be used when no chains have been chosen explicitly for fitting to density map. "<<endl; 
                ErrorManager::instance.treatError();
        }
}

void DensityContainer::printDensityStretches(){
    cout<<__FILE__<<":"<<__LINE__<<" About to print all "<< numDensityStretches()<<" density stretches "<<endl;
    for (int i = 0 ; i < numDensityStretches(); i++) {
        DensityStretch myDensityStretch = getDensityStretch(i);
        cout<<__FILE__<<":"<<__LINE__<<" Density Stretch "<<i<<" : chain "<< myDensityStretch.getChain()<<" from residue "<<myDensityStretch.getStartResidue().outString()<<" to "<< myDensityStretch.getEndResidue().outString()<<endl;
    }
}

