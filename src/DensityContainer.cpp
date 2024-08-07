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
    residueStretchVector.clear();
}

void DensityContainer::updateDensityStretch(int id, const DensityStretch & stretch, BiopolymerClassContainer & myBiopolymerClassContainer){
    if(id < 0 || id >= residueStretchVector.size()){
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "you tried to delete a non existing Contact." << endl);
    }
    validateResidueStretch(stretch, myBiopolymerClassContainer);
    residueStretchVector[id] = stretch;
}
/*
void DensityContainer::deleteDensityStretch(int id){
    if(id < 0 || id >= densityStretchVector.size()){
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "you tried to delete a non existing Contact." << endl);
    }
    densityStretchVector.erase(densityStretchVector.begin()+id);    
}*/
/*
DensityStretch DensityContainer::getDensityStretch(int densityStretchIndex){
	return densityStretchVector[densityStretchIndex];
}

int DensityContainer::numDensityStretches(){
	return densityStretchVector.size();
}*/

void DensityContainer::stuffDensityStretchVector( BiopolymerClassContainer & myBiopolymerClassContainer){
	if (numDensityStretches() == 0) {
        MMBLOG_FILE_FUNC_LINE(INFO, "All residues in all available biopolymer chains be fitted to the map. "<<endl);
		for (auto i = 0 ; i < myBiopolymerClassContainer.getNumBiopolymers(); i++){
			DensityStretch myDensityStretch;
			myDensityStretch.setChain( myBiopolymerClassContainer.updBiopolymerClass(i).getChainID());
			myDensityStretch.setStartResidue( myBiopolymerClassContainer.updBiopolymerClass(i).getFirstResidueID    ());
			myDensityStretch.setEndResidue( myBiopolymerClassContainer.updBiopolymerClass(i).getLastResidueID    ());
			myDensityStretch.printStretch();
                        //MMBLOG_FILE_FUNC_LINE(DEBUG, "Density Stretch  chain "<< myDensityStretch.getChain()<<" from residue "<<myDensityStretch.getStartResidue().outString()<<" to "<< myDensityStretch.getEndResidue().outString()<<endl);
			addStretch(myDensityStretch);//,myBiopolymerClassContainer);
		}
	} else {
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "The stuffDensityStretchVector can only be used when no chains have been chosen explicitly for fitting to density map. "<<endl);
        }
}
/*
void DensityContainer::printDensityStretches(){
    MMBLOG_FILE_FUNC_LINE(INFO, "About to print all "<< numDensityStretches()<<" density stretches "<<endl);
    for (auto i = 0 ; i < numDensityStretches(); i++) {
        DensityStretch myDensityStretch = getDensityStretch(i);
        MMBLOG_FILE_FUNC_LINE(INFO, "Density Stretch "<<i<<" : chain "<< myDensityStretch.getChain()<<" from residue "<<myDensityStretch.getStartResidue().outString()<<" to "<< myDensityStretch.getEndResidue().outString()<<endl);
    }
}*/

