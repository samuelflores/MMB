/* -------------------------------------------------------------------------- *
 *                           MMB (MacroMoleculeBuilder)                       *
 * -------------------------------------------------------------------------- *
 *                                                                            *
 * Copyright (c) 2011-12 by the Author.                                       *
 * Author: Samuel Flores                                                      *
 *                                                                            *
 * See RNABuilder.cpp for the copyright and usage agreement.                  *
 * -------------------------------------------------------------------------- */


#include "DisplacementContainer.h"


void DisplacementContainer::validateDisplacement(const Displacement displacement,  BiopolymerClassContainer & myBiopolymerContainer ){
 
    myBiopolymerContainer.validateChainID(displacement.chain); 
    ValidateVec3(displacement.displacement);
    if (hasChain(displacement.chain)){
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "You cannot create more than one displacement per chain"<<endl);
    }
}

bool DisplacementContainer::hasChain(String chain){
    for (int i = 0; i < numDisplacements(); i++){
	if (chain.compare(updDisplacement(i).chain) == 0)  
	    return true;
    }   
    return false;
}  

Displacement & DisplacementContainer::updDisplacement(String chain) {
    for (int i = 0; i < numDisplacements(); i++){
	if (chain.compare(updDisplacement(i).chain) == 0)  
	    return updDisplacement(i);
    }
    // if requested chain was not found, return an error
    MMBLOG_FILE_FUNC_LINE(CRITICAL, "You have requested a chain for which no Displacement is available : "<<chain<<". You will first need to call the initialDisplacement command, even if it is only to issue a displacement of 0,0,0.  "<<endl);
}


Vec3 DisplacementContainer::getInitialDisplacementVec3(String chain){
   if (! hasChain(chain)) {
       return Vec3(0);
   } else {
       return updDisplacement( chain).displacement;
   }
}
