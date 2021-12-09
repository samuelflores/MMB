/* -------------------------------------------------------------------------- *
 *                           MMB (MacroMoleculeBuilder)                       *
 * -------------------------------------------------------------------------- *
 *                                                                            *
 * Copyright (c) 2011-12 by the Author.                                       *
 * Author: Samuel Flores                                                      *
 *                                                                            *
 * See RNABuilder.cpp for the copyright and usage agreement.                  *
 * -------------------------------------------------------------------------- */

#include "BasePairContainer.h"
#include "BiopolymerClass.h"
#include "Utils.h"          
#include "NtC_Class_Container.h"
#include <sstream>

void BasePairContainer::clear(){
    myBasePairVector.clear();
};

void BasePairContainer::addBasePair(BiopolymerClassContainer & myBiopolymerClassContainer, 
                                    const LeontisWesthofClass & lhClass, 
                                    BaseInteraction myBasePair, bool helicalStacking){
    validateBasePair(myBiopolymerClassContainer,  lhClass, myBasePair, helicalStacking);
    myBasePairVector.push_back(myBasePair);
}

void BasePairContainer::deleteBasePair(int basePairIndex) {
    MMBLOG_FILE_FUNC_LINE(INFO, "about to delete BASE PAIR " << basePairIndex <<" of "<<numBasePairs()<<endl);
    (myBasePairVector).erase((myBasePairVector).begin()+basePairIndex);          
}

void BasePairContainer::updateBasePair(int index, String ch1, int res1, String edge1, 
                                                  String ch2, int res2, String edge2, 
                                       String orient,
                                       BiopolymerClassContainer& myBiopolymerClassContainer,
                                       const LeontisWesthofClass& lhClass,
                                       bool helicalStacking){
    BaseInteraction original;
    vector<BaseInteraction>::iterator it;

    try{
        // validate index
        myBasePairVector.at(index);

        // copy of the base pair
        original = BaseInteraction(myBasePairVector[index]);
        // remove the current base pair to avoid validation errors
        it = myBasePairVector.erase(myBasePairVector.begin()+index);

        // new BaseInteraction with the new values
        BaseInteraction bi;
        bi.FirstBPChain = ch1;
        bi.SecondBPChain = ch2;
        bi.FirstBPResidue = ResidueID(res1, ' ');
        bi.SecondBPResidue = ResidueID(res2, ' ');
        bi.FirstBPEdge = edge1;
        bi.SecondBPEdge = edge2;
        bi.OrientationBP = orient;

        validateBasePair(myBiopolymerClassContainer, lhClass, bi, helicalStacking);
        
        myBasePairVector.insert(it, bi);
    }
    catch(const out_of_range& oor){
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "BaseInteraction index " << index << " invalid for update." << endl);
    }
    catch(const MMBException& mmbError){
        myBasePairVector.insert(it, original);
        throw mmbError;
    }
}

String BasePairContainer::getBasePairsStrings(){
    stringstream returnStr;
    for(size_t i = 0; i < myBasePairVector.size(); i++){
        returnStr << i << " " 
                     << myBasePairVector[i].FirstBPChain << " "
                     << myBasePairVector[i].FirstBPResidue.getResidueNumber() << " "
                     << myBasePairVector[i].FirstBPEdge << " "
                     << myBasePairVector[i].SecondBPChain << " "
                     << myBasePairVector[i].SecondBPResidue.getResidueNumber() << " "
                     << myBasePairVector[i].SecondBPEdge << " "
                     << myBasePairVector[i].OrientationBP << endl;
    }

    return returnStr.str();
}

void BasePairContainer::validateBasePair(BiopolymerClassContainer & myBiopolymerClassContainer, 
                                         const LeontisWesthofClass & lhClass,
                                         BaseInteraction & myBasePair, bool helicalStacking){
    myBiopolymerClassContainer.updBiopolymerClass(myBasePair.FirstBPChain).validateResidueID(myBasePair.FirstBPResidue);
    myBiopolymerClassContainer.updBiopolymerClass(myBasePair.SecondBPChain).validateResidueID(myBasePair.SecondBPResidue);
    if (!((myBiopolymerClassContainer.updBiopolymerClass(myBasePair.FirstBPChain).biopolymerType ==  BiopolymerType::RNA) ||
          (myBiopolymerClassContainer.updBiopolymerClass(myBasePair.FirstBPChain).biopolymerType == BiopolymerType::DNA))) {
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "You have attempted to create a baseInteraction for a chain which is of type : "<<myBiopolymerClassContainer.updBiopolymerClass(myBasePair.FirstBPChain).biopolymerType<<" rather than of DNA or RNA."<<endl);
    }

    if (!((myBiopolymerClassContainer.updBiopolymerClass(myBasePair.SecondBPChain).biopolymerType == BiopolymerType::RNA) ||
          (myBiopolymerClassContainer.updBiopolymerClass(myBasePair.SecondBPChain).biopolymerType == BiopolymerType::DNA))) {
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "You have attempted to create a baseInteraction for a chain which is of type : "<<myBiopolymerClassContainer.updBiopolymerClass(myBasePair.SecondBPChain).biopolymerType<<" rather than of DNA or RNA."<<endl);
    }
    if ((myBasePair.FirstBPResidue == myBasePair.SecondBPResidue) &&
        (myBasePair.FirstBPChain == myBasePair.SecondBPChain)) {
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "You have attempted to create a baseInteraction between chain "<< myBasePair.FirstBPChain<< ", residue "<<myBasePair.FirstBPResidue.outString()<<" and itself.  This is not allowed!"<<endl);
    }

    if(helicalStacking && 
       myBasePair.FirstBPEdge=="WatsonCrick" && 
       myBasePair.SecondBPEdge=="WatsonCrick" &&
       myBasePair.OrientationBP=="Cis")
    {
        if(hasWatsonCrickCisPair(myBasePair.FirstBPChain, myBasePair.FirstBPResidue)){
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "There is already one WatsonCrick/WatsonCrick/Cis base pair assigned to chain "<<myBasePair.FirstBPChain<<", residue "<<myBasePair.FirstBPResidue.outString()<<".  If you insist on this, please turn off setHelicalStacking . "<<endl);
        }
        if(hasWatsonCrickCisPair(myBasePair.SecondBPChain, myBasePair.SecondBPResidue)){
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "There is already one WatsonCrick/WatsonCrick/Cis base pair assigned to chain "<<myBasePair.SecondBPChain<<", residue "<<myBasePair.SecondBPResidue.outString()<<".  If you insist on this, please turn off setHelicalStacking . "<<endl);
        }
    }

    // Check Leontis and Westhof coherence
    String resName1 = myBiopolymerClassContainer.getPdbResidueName(myBasePair.FirstBPChain, myBasePair.FirstBPResidue);
    String resName2 = myBiopolymerClassContainer.getPdbResidueName(myBasePair.SecondBPChain, myBasePair.SecondBPResidue);
    // if getLeontisWesthofBondRow doesn't throw an exception, it's good
    LeontisWesthofBondRow br = lhClass.getLeontisWesthofBondRow(
                                                           myBasePair.FirstBPResidue, myBasePair.SecondBPResidue,
                                                           resName1, myBasePair.FirstBPEdge,
                                                           resName2, myBasePair.SecondBPEdge,
                                                           myBasePair.OrientationBP,
                                                           "baseInteraction"
                                                          );
    myBasePair.leontisWesthofBondRowIndex = lhClass.getLeontisWesthofBondRowIndex(
                                                           resName1, resName2,
                                                           myBasePair.FirstBPEdge,
                                                           myBasePair.SecondBPEdge,
                                                           myBasePair.OrientationBP,
                                                           "baseInteraction"
                                                          ) ;

}


const BaseInteraction & BasePairContainer::getBasePair(int basePairIndex) {
    return myBasePairVector[basePairIndex];
}
/*BasePair & BasePairContainer::updBasePair(int basePairIndex)
    return myBasePairVector[basePairIndex] ;
}*/

int BasePairContainer::numBasePairs() {
    return myBasePairVector.size();
}

bool BasePairContainer::hasWatsonCrickCisPair(String chainID ,ResidueID residueNumber) { //,String secondChain,int secondResidueNumber,String orientation){
    // should also ensure that redundant WatsonCrick's don't exist.
    int numWatsonCrickPairs = 0;
    for (int i = 0; i < numBasePairs(); i++) {
        if (  (getBasePair(i).FirstBPEdge.compare("WatsonCrick") == 0)  &&
              (getBasePair(i).SecondBPEdge.compare("WatsonCrick") == 0) &&
              (getBasePair(i).OrientationBP.compare("Cis") == 0)          &&
            (((getBasePair(i).FirstBPResidue == residueNumber) && (getBasePair(i).FirstBPChain.compare(chainID) == 0  )) || 
             ((getBasePair(i).SecondBPResidue == residueNumber) && (getBasePair(i).SecondBPChain.compare(chainID) == 0)))   
           ) {
            numWatsonCrickPairs++;
            //return True;
        }

    }
    if (numWatsonCrickPairs == 0) {
        return false ; // no WC pairs found for this chainID and residueNumber.
    } else if (numWatsonCrickPairs == 1) {
        return true; // one WC pair found for this chainID and residueNumber.
    } else {
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "There is more than one WatsonCrick/WatsonCrick/Cis base pair assigned to chain "<<chainID<<", residue "<<residueNumber.outString()<<".  If you insist on this, please turn off setHelicalStacking . "<<endl);
    }
}

const BaseInteraction BasePairContainer::getWatsonCrickCisPair(String chainID ,ResidueID residueNumber) { 
    // this also ensures that redundant WatsonCrick's don't exist.
    if (! hasWatsonCrickCisPair(chainID, residueNumber)) {
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "You have attempted to retrieve a baseInteraction for chain "<<chainID<<" residue "<<residueNumber.outString()<<" which doesn't exist. "<<endl);
    } else {
        for (int i = 0; i < numBasePairs(); i++) {
            if (  (getBasePair(i).FirstBPEdge.compare("WatsonCrick") == 0)  &&
                  (getBasePair(i).SecondBPEdge.compare("WatsonCrick") == 0) &&
                  (getBasePair(i).OrientationBP.compare("Cis") == 0)          &&
                (((getBasePair(i).FirstBPResidue == residueNumber) && (getBasePair(i).FirstBPChain.compare(chainID) == 0  )) || 
                 ((getBasePair(i).SecondBPResidue == residueNumber) && (getBasePair(i).SecondBPChain.compare(chainID) == 0)))   
               ) 
            {
                return getBasePair(i); 
            }
        } // of for i
    } // of else
    MMBLOG_FILE_FUNC_LINE(CRITICAL, "UNKNOWN ERROR -- shouldn't get here" << endl);
} // of method




const ResidueID BasePairContainer::getWatsonCrickCisPairingResidue(String chainID ,ResidueID residueNumber) { 
    BaseInteraction myBasePair = getWatsonCrickCisPair( chainID , residueNumber);
    if ((myBasePair.FirstBPResidue == residueNumber) &&
        (myBasePair.FirstBPChain == chainID))  {
        return myBasePair.SecondBPResidue;
    } else if  ((myBasePair.SecondBPResidue == residueNumber) &&
              (myBasePair.SecondBPChain == chainID))  {
        return myBasePair.FirstBPResidue;
    } else {
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "Unexplained error! "<<endl);
    }
}

const String BasePairContainer::getWatsonCrickCisPairingChain(String chainID ,ResidueID residueNumber) {
    BaseInteraction myBasePair = getWatsonCrickCisPair( chainID , residueNumber);
    if ((myBasePair.FirstBPResidue == residueNumber) &&
        (myBasePair.FirstBPChain == chainID))  {
        return myBasePair.SecondBPChain;
    } else if  ((myBasePair.SecondBPResidue == residueNumber) &&
              (myBasePair.SecondBPChain == chainID))  {
        return myBasePair.FirstBPChain;
    } else {
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "Unexplained error! "<<endl);
    }
}

// Given a residue and chain, this gets the last residue of the same chain, that is part of the same helix.
// This means that it must interact with consecutive residues of the same or different chain.
const ResidueID  BasePairContainer::getLastWatsonCrickCisPairingResidueOfRun(String chainID ,ResidueID firstResidueNumberInStack, BiopolymerClassContainer & myBiopolymerClassContainer) { 
    if (! hasWatsonCrickCisPair(  chainID,firstResidueNumberInStack) ){ 
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "The residue number you provided, "<< firstResidueNumberInStack.outString()<<" is not engaged in any WatsonCrick Cis interactions!  "<<endl);
        // return firstResidueNumberInStack;
    }
            

    ResidueID    firstInteractingResidue =  getWatsonCrickCisPairingResidue(  chainID,firstResidueNumberInStack);
    String firstInteractingChain   =  getWatsonCrickCisPairingChain  (  chainID, firstResidueNumberInStack);
    BiopolymerClass  myBiopolymerClass = myBiopolymerClassContainer.updBiopolymerClass(chainID); 
    BiopolymerClass  myInteractingBiopolymerClass = myBiopolymerClassContainer.updBiopolymerClass(firstInteractingChain); 
    for (ResidueID lastResidueNumberInStack = firstResidueNumberInStack ; lastResidueNumberInStack <=  myBiopolymerClass.getLastResidueID(); myBiopolymerClass.incrementResidueID (lastResidueNumberInStack) ) {
                if (hasWatsonCrickCisPair(chainID, lastResidueNumberInStack)) {
                    if (lastResidueNumberInStack == myBiopolymerClass.getLastResidueID()) {
                        // the last residue number in the stack cannot be any higher than the last residue number.
                        return lastResidueNumberInStack;
                    }
                    ResidueID    lastInteractingResidue =  getWatsonCrickCisPairingResidue(chainID, lastResidueNumberInStack);
                    String lastInteractingChain   =  getWatsonCrickCisPairingChain  (chainID, lastResidueNumberInStack);

                    MMBLOG_FILE_FUNC_LINE(INFO, "firstInteractingChain = >"<<firstInteractingChain<<"< lastInteractingChain = >"<<lastInteractingChain<<"< "<<endl);
                    if (firstInteractingChain.compare(lastInteractingChain) != 0) {
                        MMBLOG_FILE_FUNC_LINE(INFO, "Residue "<<lastResidueNumberInStack.outString()<<" of chain "<<chainID<<" is interacting with chain "<<lastInteractingChain<<" Residue "<<lastInteractingResidue.outString()<<endl);
                        MMBLOG_FILE_FUNC_LINE(INFO, "Whereas  "<<firstResidueNumberInStack.outString()<<" of chain "<<chainID<<" is interacting with chain "<<firstInteractingChain<<" Residue "<<firstInteractingResidue.outString()<<endl);
                        MMBLOG_FILE_FUNC_LINE(INFO, "Whereas the first residue of the stack, "<<firstResidueNumberInStack.outString()<<" of chain "<<chainID<<" was interacting with residue "<<firstInteractingResidue.outString()<<" of chain "<<firstInteractingChain<<endl);
                        MMBLOG_FILE_FUNC_LINE(INFO, "Since the interacting chains are different, residue "<<lastResidueNumberInStack.outString() <<" of chain "<<chainID<<" is not part of the stack"<<endl);
                        return (myBiopolymerClass.decrementResidueID(lastResidueNumberInStack) );
                    }

                    // lastInteractingResidue is progressively decrementing as we increment lastResidueNumberInStack, if it is still part of the growing stack.
                    // therefore both sides of the comparison below should be >0 if lastInteractingResidue is still in the stack.
                    if ((myInteractingBiopolymerClass.difference(firstInteractingResidue , lastInteractingResidue) == myBiopolymerClass.difference (lastResidueNumberInStack , firstResidueNumberInStack)) &&
                        (lastInteractingChain.compare(firstInteractingChain)==0)){
                            // lastResidueNumberInStack is indeed still part of the stack.  Do nothing.
                            MMBLOG_FILE_FUNC_LINE(INFO, "firstInteractingChain = "<<firstInteractingChain<<" lastInteractingChain = "<<lastInteractingChain<<endl);
                        } else {
                            // lastResidueNumberInStack is not part of the stack.  
                            return (myBiopolymerClass.decrementResidueID(lastResidueNumberInStack) );
                        }
                }
                else {
                    // lastResidueNumberInStack is not part of the stack.  
                    return myBiopolymerClass.decrementResidueID(lastResidueNumberInStack );
                }
    } // of for lastResidueNumberInStack
    MMBLOG_FILE_FUNC_LINE(CRITICAL, "Unexplained error! "<<endl);
} // of getLastWatsonCrickCisPairingResidueOfRun

//void BasePairContainer::     addHelicalStackingInteractions(String chainID, ResidueID firstResidue, ResidueID lastResidue,BiopolymerClassContainer & myBiopolymerClassContainer, const LeontisWesthofClass & lhClass){
/*void BasePairContainer::generateAorBFormNtCs(String chainID, ResidueID firstResidue, ResidueID lastResidue,BiopolymerClassContainer & myBiopolymerClassContainer, const LeontisWesthofClass & lhClass){
    if ( myBiopolymerClassContainer.updBiopolymerClass(chainID).difference(lastResidue , firstResidue) <1) {
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "It's not possible to apply helical stacking interactions to a run of fewer than 2 residues! "<<endl);
    }
}*/

void BasePairContainer::generateHelicalStackingInteractions(String chainID, ResidueID firstResidue, ResidueID lastResidue,BiopolymerClassContainer & myBiopolymerClassContainer, const LeontisWesthofClass & lhClass){

    if ( myBiopolymerClassContainer.updBiopolymerClass(chainID).difference(lastResidue , firstResidue) <1) {
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "It's not possible to apply helical stacking interactions to a run of fewer than 2 residues! "<<endl);
    }
    for (ResidueID i = firstResidue; i < lastResidue;myBiopolymerClassContainer.updBiopolymerClass(chainID).incrementResidueID( i) ) {
        BaseInteraction myBasePair;
        myBasePair.rotationCorrection1 = Rotation(0.0,UnitVec3(0,0,1));
        myBasePair.rotationCorrection2 = Rotation(0.0,UnitVec3(0,0,1));
        myBasePair.FirstBPChain    = chainID;
        myBasePair.SecondBPChain   = chainID;
        myBasePair.FirstBPResidue  =  i;
        myBasePair.SecondBPResidue =  myBiopolymerClassContainer.updBiopolymerClass(chainID).incrementResidueID( i);
        myBiopolymerClassContainer.updBiopolymerClass(chainID).decrementResidueID( i);
        if      (myBiopolymerClassContainer.updBiopolymerClass(chainID).getBiopolymerType() == BiopolymerType::RNA) 
        {
            myBasePair.FirstBPEdge     = "HelicalStackingA3";
            myBasePair.SecondBPEdge    = "HelicalStackingA5";
        } else if (myBiopolymerClassContainer.updBiopolymerClass(chainID).getBiopolymerType() == BiopolymerType::DNA )
        {
            myBasePair.FirstBPEdge     = "HelicalStackingB3";
            myBasePair.SecondBPEdge    = "HelicalStackingB5";
        } else {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "Invalid biopolymerType: "<< myBiopolymerClassContainer.updBiopolymerClass(chainID).getBiopolymerType() <<endl);
        }
        myBasePair.OrientationBP = "Cis";
        addBasePair( myBiopolymerClassContainer,lhClass,myBasePair, true); 
    }
    
}

// Add helical stacking interactions for all runs of Watson-Crick base pairs, beyond some minimum lenght

void BasePairContainer::addHelicalStacking(BiopolymerClassContainer & myBiopolymerClassContainer, const LeontisWesthofClass & lhClass, const  NTC_PAR_Class & ntc_par_class , NTC_Class_Container  & my_ntc_class_container){
    MMBLOG_FILE_FUNC_LINE(INFO, "Running addHelicalStacking  "<<endl);
    for (auto i = 0; i < myBiopolymerClassContainer.getNumBiopolymers(); i++) {
        //MMBLOG_FILE_FUNC_LINE(" So far we have "<<numBasePairs() <<" baseInteraction's.  They are: "<<endl;
        BiopolymerClass & myBiopolymerClass = myBiopolymerClassContainer.updBiopolymerClass(i);
        String myChainID = myBiopolymerClass.getChainID();
        MMBLOG_FILE_FUNC_LINE(INFO, "Running addHelicalStacking for chain "<<myChainID<<endl);
        //myResidueNumber counts down chain i from the first to last residue number, looking for stretches of chain i to which it can apply helical stacking runs.  at the end of each stacking run, it increments to one residue after the last residue of the run.
        ResidueID myResidue = myBiopolymerClass.getFirstResidueID(); 
        while  ( myResidue !=  myBiopolymerClass.getLastResidueID()){
            MMBLOG_FILE_FUNC_LINE(INFO, "Checking myResidue "<<myResidue.outString()<<" vs. myBiopolymerClass.getLastResidueID() = "<<myBiopolymerClass.getLastResidueID().outString()<<endl);
            if ( hasWatsonCrickCisPair(myChainID,myResidue) ){ 
                ResidueID lastPairingResidue = getLastWatsonCrickCisPairingResidueOfRun(myChainID,myResidue,myBiopolymerClassContainer);
                MMBLOG_FILE_FUNC_LINE(INFO, "Generating stacking interactions from "<<myResidue.outString()<<" to "<<lastPairingResidue.outString()<<endl);
                if ( myBiopolymerClass.difference (lastPairingResidue , myResidue) > 1){ // if we have at least 3BP in the helix.  increase this if you wish, but minimum is zero.
                    generateHelicalStackingInteractions(myChainID,myResidue,lastPairingResidue,myBiopolymerClassContainer, lhClass) ;
                    #ifdef Build_NtC
                    my_ntc_class_container.generateAorBFormNtCs( myBiopolymerClassContainer, myChainID, myResidue, lastPairingResidue, 0.5 , ntc_par_class);
                    #endif
                }
                if (lastPairingResidue <  myBiopolymerClass.getLastResidueID()) {
                    myResidue = myBiopolymerClass.incrementResidueID(lastPairingResidue) ; // continue the outer loop after the just-discovered stacking run.
                    if (myResidue ==  myBiopolymerClass.getLastResidueID()) break ; // If we've just incremented to the end of the chain, we're done here.
                }
                else if (lastPairingResidue ==  myBiopolymerClass.getLastResidueID()) break ; // we've already set all helical stackings for this chain.
                else { // lastPairingResidue > myBiopolymerClass.getLastResidueID() .. can't happen!
                    MMBLOG_FILE_FUNC_LINE(CRITICAL, "Unexplained error!"<<endl);
                }
            } else if (myResidue ==  myBiopolymerClass.getLastResidueID()) {
                break; // we're at the last residue of the chain
            // Catchall :
            } else {//if (myResidue < myBiopolymerClass.getLastResidueID() ){
                myResidue = myBiopolymerClass.incrementResidueID(myResidue); // increment myResidue
                if (myResidue ==  myBiopolymerClass.getLastResidueID()) break ; // If we've just incremented to the end of the chain, we're done here.
            } 

            if ( myResidue ==  myBiopolymerClass.getLastResidueID()) { // this should have been caught above.
                MMBLOG_FILE_FUNC_LINE(CRITICAL, "Unexplained error!"<<endl);
            }
            
        } // of for myResidue
    } // of for i
    MMBLOG_FILE_FUNC_LINE(INFO, "Done adding helical stackings. So far we have "<<numBasePairs() <<" baseInteraction's.  They are: "<<endl);
    printBasePairs();
} // of addHelicalStacking


void BasePairContainer::printBasePairs() {
    MMBLOG_FILE_FUNC_LINE(INFO, "Printing all base pairs:"<<endl);
    for (auto i = 0; i < numBasePairs(); i++) {
        MMBLOG_FILE_FUNC_LINE(INFO, getBasePair(i).FirstBPChain<<" "<<getBasePair(i).FirstBPResidue.outString()<<" "<<   getBasePair(i).FirstBPEdge<<" "<<getBasePair(i).SecondBPChain<<" "<<getBasePair(i).SecondBPResidue.outString()<<" "<<   getBasePair(i).SecondBPEdge<<" "<<  getBasePair(i).OrientationBP<<endl);
    }
}

void BasePairContainer::setBasePairSatisfied(int basePairIndex ,bool isSatisfied) {
    myBasePairVector[basePairIndex].basePairSatisfied = isSatisfied;    
}

vector<int> BasePairContainer::getSatisfiedBasePairs(){
    vector<int> bps;
    for (auto i = 0; i < numBasePairs(); i++) {
        if(getBasePair(i).basePairSatisfied){
            bps.push_back(i);
        }
    }

    return bps;
}

