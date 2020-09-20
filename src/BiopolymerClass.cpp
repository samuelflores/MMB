/* -------------------------------------------------------------------------- *
 *                           MMB (MacroMoleculeBuilder)                       *
 * -------------------------------------------------------------------------- *
 *                                                                            *
 * Copyright (c) 2011-12 by the Author.                                       *
 * Author: Samuel Flores                                                      *
 *                                                                            *
 * See RNABuilder.cpp for the copyright and usage agreement.                  *
 * -------------------------------------------------------------------------- */

//#include <boost/algorithm/string.hpp>
//#include <algorithm>
#include "BiopolymerClass.h"
#include "Utils.h"
#include "ContactContainer.h"
#include "ResidueStretchContainer.h"
//#include <string>
#include <map>
#include <set>
#include <cstdlib>
//#include <stdlib.h>
#include "ReferenceNeighborList.h"
#include  <utility>

// #define  _DEBUG_FLAGS_ON_

using namespace std;
using namespace SimTK;

    //template <class ResidueStretchType>
    //class ResidueStretchContainer<ResidueStretchType>;
    //template <class ResidueStretchType>
    //const int ResidueStretchContainer<ResidueStretchType>::getNumResidueStretches();

void printBiopolymerSequenceInfo(const Biopolymer & myBiopolymer) {
    for (int i = 0; i < myBiopolymer.getNumResidues(); i++) {
        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Residue type, number, and insertion code: "<<myBiopolymer.getResidue(ResidueInfo::Index(i)).getOneLetterCode() <<", "<<myBiopolymer.getResidue(ResidueInfo::Index(i)).getPdbResidueNumber()<<", "<<myBiopolymer.getResidue(ResidueInfo::Index(i)).getPdbInsertionCode()<<endl;
    } 
};

bool letterIsPurine(String myLetter) {
    if ((myLetter.compare("A") == 0 ) ||  
        (myLetter.compare("G") == 0 ))    
            {return true; }
    else return false;
}

bool letterIsRNA(String myLetter) {
    if (myLetter.compare("A") == 0)  {return true; }
    else if (myLetter.compare("C") == 0) {return true;}
    else if (myLetter.compare("G") == 0) {return true;}
    else if (myLetter.compare("U") == 0) {return true;}
    else {
        cout <<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": You have specified a non-RNA residue, single letter code = "<<myLetter<<endl;
        //ErrorManager::instance <<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": You have specified a non-RNA residue, single letter code = "<<myLetter<<endl;
        //ErrorManager::instance.treatError(); 
        return false;}
}

bool letterIsDNA(String myLetter) {
    if (myLetter.compare("A") == 0)  {return true; }
    else if (myLetter.compare("C") == 0) {return true;}
    else if (myLetter.compare("G") == 0) {return true;}
    else if (myLetter.compare("T") == 0) {return true;}
    else {
        cout <<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": You have specified a non-DNA residue, single letter code = "<<myLetter<<endl;
        //ErrorManager::instance <<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": You have specified a non-DNA residue, single letter code = "<<myLetter<<endl;
        //ErrorManager::instance.treatError(); 
        return false;}
}


bool letterIsProtein(String   myLetter) {
    if      (myLetter.compare("C") == 0) {return true;}
    else if (myLetter.compare("X") == 0) {return true;}
    else if (myLetter.compare("H") == 0) {return true;}
    else if (myLetter.compare("I") == 0) {return true;}
    else if (myLetter.compare("M") == 0) {return true;}
    else if (myLetter.compare("S") == 0) {return true;}
    else if (myLetter.compare("V") == 0) {return true;}
    else if (myLetter.compare("A") == 0) {return true;}
    else if (myLetter.compare("G") == 0) {return true;}
    else if (myLetter.compare("L") == 0) {return true;}
    else if (myLetter.compare("P") == 0) {return true;}
    else if (myLetter.compare("T") == 0) {return true;}
    else if (myLetter.compare("R") == 0) {return true;}
    else if (myLetter.compare("F") == 0) {return true;}
    else if (myLetter.compare("Y") == 0) {return true;}
    else if (myLetter.compare("W") == 0) {return true;}
    else if (myLetter.compare("D") == 0) {return true;}
    else if (myLetter.compare("N") == 0) {return true;}
    else if (myLetter.compare("E") == 0) {return true;}
    else if (myLetter.compare("Q") == 0) {return true;}
    else if (myLetter.compare("K") == 0) {return true;}
    else {
        //ErrorManager::instance <<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<< ": The symbol " << myLetter << " is not in the protein alphabet\n";
        cout <<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<< ": The symbol " << myLetter << " is not in the protein alphabet"<<endl;;
                //ErrorManager::instance.treatError();
                return false;
    } 

}


void   BiopolymerClass:: modifyResidue( const BiopolymerModification myBiopolymerModification,Compound  compoundToAdd,  DuMMForceFieldSubsystem & dumm){//DuMMForceFieldSubsystem & dumm) { 

        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<endl;
        Compound addedAtom = UnivalentAtom("HG", Element::Hydrogen());
        DuMM::ChargedAtomTypeIndex      myChargedAtomTypeIndex = dumm.getBiotypeChargedAtomType( compoundToAdd.getAtomBiotypeIndex(Compound::AtomIndex( 0)));
            //dumm.getNextUnusedChargedAtomTypeIndex (); 
        String residueName ("Cysteine (-SH)");
        String specificAtomName = myBiopolymerModification.getAtomOnAddedCompound();
        stringstream ss; ss<< getResidueIndex( myBiopolymerModification.getResidueToModify ());
        String residueIndexString = ss.str();
        String longAtomName = residueIndexString +String("/") + specificAtomName ;
        String chargedAtomTypeName("blah");
        Ordinality::Residue myOrdinality(SimTK::Ordinality::Any);
        /*Biotype::defineBiotype (
                Element::getBySymbol("H"),
                1,                           // valence
                residueName,
                specificAtomName             //should this be specificAtomName or genericAtomName?
            );  */
        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" compoundToAdd.getNumAtoms() "<< compoundToAdd.getNumAtoms()<<endl;
        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" compoundToAdd.getAtomBiotypeIndex(0) "<< compoundToAdd.getAtomBiotypeIndex(Compound::AtomIndex( 0)) <<endl;
        String bondName =  "0/SG/bond2";
        double bondLength = .14;
        Angle myDihedral = 180*Deg2Rad;
        BondMobility::Mobility myBondMobility = stringToBondMobility("Rigid");
        Compound myCompound(UnivalentAtom(specificAtomName, Element::Hydrogen()));
        myCompound.setPdbResidueNumber(1);
        myCompound.setPdbChainId("A");
        myCompound.setPdbResidueName("CYX");
        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<endl;
        myBiopolymer.bondCompound(
            residueIndexString , // This should be the residue index.  Will be prepended with / to atom name
            myCompound ,
            bondName,     
            bondLength,
            myDihedral);
        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<endl;
        Compound::AtomIndex myAtomIndex = myBiopolymer.getAtomIndex( Compound::AtomPathName (longAtomName) );
        Compound::AtomPathName tempAtomName("HG");
        ResidueInfo myResidueInfo = myBiopolymer.updResidue(ResidueInfo::Index(0));
        myBiopolymer.updResidue(ResidueInfo::Index(0)) .addAtom(
		 myAtomIndex,
                 specificAtomName);            
        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" myBiopolymer.getAtomBiotypeIndex(Compound::AtomIndex( myBiopolymer.getAtomIndex(longAtomName))) : "<<  myBiopolymer.getAtomBiotypeIndex(Compound::AtomIndex( myBiopolymer.getAtomIndex(longAtomName))) <<endl;
        myBiopolymer.setBiotypeIndex(longAtomName, compoundToAdd.getAtomBiotypeIndex(Compound::AtomIndex( 0)));
        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" myBiopolymer.getAtomBiotypeIndex(Compound::AtomIndex( myBiopolymer.getAtomIndex(longAtomName))) : "<<  myBiopolymer.getAtomBiotypeIndex(Compound::AtomIndex( myBiopolymer.getAtomIndex(longAtomName))) <<endl;
        /*myBiopolymer.nameAtom( Compound::AtomName ("0/HG"),  Compound::AtomName ("HG"));
        PdbChain myPdbChain = PdbChain(myBiopolymer,Transform());
        PdbResidue  myPdbResidue = myPdbChain.updResidue(Pdb::ResidueIndex(0));
        myPdbResidue = PdbResidue(myBiopolymer, int (0), Transform()); 
        PdbAtom myPdbAtom(myCompound,String(specificAtomName),Transform());
        myPdbResidue.addAtom(myPdbAtom);
        PdbAtom myPdbAtom2 = myPdbResidue.getAtom(String("HG"));*/
}


String  BiopolymerClass::getSubSequence(const ResidueID startResidue, const ResidueID endResidue)
{
    validateResidueID(startResidue);
    validateResidueID(  endResidue);
    String subSequence = sequence.substr(getResidueIndex(startResidue ), (difference (endResidue, startResidue) +1)) ; // gets the subsequence bounded by startResidue, endResidue inclusive. If we omit the +1, then we would not get endResidue.	
    //std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" You have requested the subsequence of chain "<<getChainID()<<", which is : "<<sequence<<" bounded by "<<startResidue.outString()     <<", "<<endResidue.outString()<<". Result is: "<<subSequence<<std::endl;
    return subSequence;
}

// prints out cleaned residue numbers (staring at 1) and correspondence to 

String BiopolymerClass::printOriginalAndRenumberedResidueIDs(const String myPdbId) {
    //stringstream returnStringStream; 
    String myQuery = "" ; 
    for (int i = 0 ; i < myBiopolymer.getNumResidues(); i++) {
        myQuery += " insert into cleanedNonCleanedResidueId  (cleanedResidueNumber,originalResidueNumber,originalInsertionCode,chainId,pdbId) VALUES ( "+ std::to_string(i + 1 - proteinCapping) + " , "  + std::to_string(getResidueID(i).getResidueNumber())+" , '" +  getResidueID(i).getInsertionCode() + "' , '"+ getChainID()  + "' , '"    + myPdbId + "') ; \n" ; 
        //String myQuery = " update cleanedNonCleanedResidueId set cleanedResidueNumber = "+ std::to_string(i + 1 - proteinCapping)+" where originalResidueNumber = "+std::to_string(getResidueID(i).getResidueNumber())+" and originalInsertionCode = '"+getResidueID(i).getInsertionCode()+"' and pdbId = '" + myPdbId + "' ; \n" ; 
        //std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" update cleanedNonCleanedResidueId set cleanedResidueNumber = "<< i + 1 - proteinCapping<<" where originalResidueNumber = "<<getResidueID(i).getResidueNumber()<<" and originalInsertionCode = '"<<getResidueID(i).getInsertionCode()<<"' and pdbId = 'XXXX';" <<endl; 
        //returnStringStream<<myQuery;  
    }
    cout<<myQuery;
    return myQuery;
    //return returnStringStream;
}


void BiopolymerClass::clear() {
    myBiopolymer =  Biopolymer();
    std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" sizeof(PdbStructure) = "<< sizeof(PdbStructure) <<std::endl;
    std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" sizeof(myBiopolymer) = "<< sizeof(myBiopolymer) <<std::endl;
    
    //firstResidueNumber = 0;
    biopolymerType =    BiopolymerType::Unassigned;
    setProteinCapping ( false);
    sequence =    "";
    originalSequence = "";
    chainID  =    "";
    atomInfoVector.clear();
    ignoreAtomPositionVector.clear();
    residueIDVector.clear();
    setFirstResidueMobilizerType(String("Free")); // set the default for this variable.  This means the root atom is connected to ground by a Free mobilizer conferring 6 DOFs.  The alternative is a Weld mobilizer, conferring 0 DOFs.
    setActivePhysics(true);
    //pdbStructure = NULL;
}

void  BiopolymerClass::validateChainID(){
    //if (chainID.length() == 1) {
    //    return 1;
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" validating chain ID >"<<chainID<<"<"<<endl;
    if (chainID.length() == 0) {
        ErrorManager::instance <<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" : Chain ID must be at least one character long.  Yours has length "<<chainID.length()<<endl; 
        ErrorManager::instance.treatError();
    } else if (chainID.compare(" ") == 0) {
        ErrorManager::instance <<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" : You are not allowed to set the chain ID to \' \' in MMB, even though that's probably kosher by the PDB."<<endl; 
        ErrorManager::instance.treatError();
    }else if (chainID.length() >= 2) 
    {
        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": In the PDB format chain ID's must be a exactly one character long. Yours has length "<< chainID.length() <<endl;
        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": MMB actually permits this!  But make sure you know how the REMARK-SimTK-long-chainId tag works. "<<endl; //And make sure you don't turn on loadSequencesFromPdb .. that won't work with long chain IDs."<<endl;
    }
}


//  Returns 0 if all is OK, otherwise returns 1
int  BiopolymerClass::checkResidueNumbersAndInsertionCodes(){
    
    int myResidueIndex = getResidueIndex(getFirstResidueID());
    ResidueID myResidueID = getFirstResidueID();
    cout <<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" : "<<getChainID() <<endl;
    while (myResidueIndex < getResidueIndex(getLastResidueID())){
        cout <<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" : "<<getResidueID(myResidueIndex).outString()<<"."<<flush <<endl;
        // Check that at least the integer part of the ResidueID is nondecreasing
        if ( getResidueID(myResidueIndex).getResidueNumber() <= getResidueID(myResidueIndex + 1).getResidueNumber()) {
            // all is well. Even if insertion codes are wonky, we can deal with that.
        } else {
           cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<"  The residue ID's are problematic! Specifically, " <<getResidueID(myResidueIndex).outString() <<" !<= "<<getResidueID(myResidueIndex + 1).outString()<<" . We can tolerate wonky insertion code ordering, but the integer part of the residueID cannot decrease, otherwise our structure matching algorithm pukes. Kindly follow non-bizarre numbering conventions."<<endl;
           return 1; // Non zero return value indicates error
           //ErrorManager::instance <<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": The residue ID's are problematic! Specifically, " <<getResidueID(myResidueIndex).outString() <<" !<= "<<getResidueID(myResidueIndex + 1).outString()<<" . We can tolerate wonky insertion code ordering, but the integer part of the residueID cannot decrease, otherwise our structure matching algorithm pukes. Kindly follow non-bizarre numbering conventions."<<endl;
           //ErrorManager::instance.treatError();
        }          

        myResidueIndex ++;
    }
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Done checking all residue numbers for chain "<<getChainID()<<" . The simple check passed, returning 0."<<endl;
    return 0;
}
void BiopolymerClass::validateResidueNumbersAndInsertionCodes(){
    /*cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" This validation step is currently not active."<<endl;
    
    int myResidueIndex = getResidueIndex(getFirstResidueID());
    ResidueID myResidueID = getFirstResidueID();
    while (myResidueIndex < getResidueIndex(getLastResidueID())){
        // cout<<getResidueID(myResidueIndex).outString()<<"."<<flush;
        // Check that at least the integer part of the ResidueID is nondecreasing
        if ( getResidueID(myResidueIndex).getResidueNumber() <= getResidueID(myResidueIndex + 1).getResidueNumber()) {
            // all is well. Even if insertion codes are wonky, we can deal with that.
        } else {
           ErrorManager::instance <<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": The residue ID's are problematic! Specifically, " <<getResidueID(myResidueIndex).outString() <<" !<= "<<getResidueID(myResidueIndex + 1).outString()<<" . We can tolerate wonky insertion code ordering, but the integer part of the residueID cannot decrease, otherwise our structure matching algorithm pukes. Kindly follow non-bizarre numbering conventions."<<endl;
           ErrorManager::instance.treatError();
        }          

        myResidueIndex ++;
    }
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<endl;*/
    if (checkResidueNumbersAndInsertionCodes()) { // returns 1 in case of problems
           ErrorManager::instance <<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": The residue ID's are problematic for chain "<<getChainID()<<" .. see message above."<<std::endl;
           ErrorManager::instance.treatError();
    }
    
}
bool BiopolymerClass::residueIsPurine (int residueIndex, String mySequence) {
    //cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<endl;
    if ((biopolymerType == BiopolymerType::RNA)  || (biopolymerType == BiopolymerType::DNA))
        {
        //cout<<residueIndex<<":"<<mySequence.substr(residueIndex,1)<<":"<<letterIsPurine(mySequence.substr(residueIndex,1))<<"."<<flush;     
        return letterIsPurine(mySequence.substr(residueIndex,1)) ;       
        } 
    else {cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": This function is intended only for nucleic acids!"<<endl; exit(1); return false;}
}

bool BiopolymerClass::residueIsPurine (int residueIndex) {
    return residueIsPurine(residueIndex,sequence);
}

int BiopolymerClass::validateSequence() {
   if (biopolymerType == BiopolymerType::Unassigned) { 
       //cout<<Unassigned<<endl;
       if (sequence.length() >0){
           ErrorManager::instance <<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": Your biopolymerType is "<<biopolymerType<<".  sequence must be zero length."          <<endl;
           ErrorManager::instance.treatError();
       } else return 0;
   } else if (biopolymerType == BiopolymerType::RNA) {
       if (sequence.length() <1){
           ErrorManager::instance <<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": Your biopolymerType is "<<biopolymerType<<".  sequence must be of length > 0." <<endl;
           ErrorManager::instance.treatError();
       } 
       for (int i = 0; i < (int)sequence.length(); i++) {
           if (!
               letterIsRNA((sequence.substr(i,1)))
              ) {
                   ErrorManager::instance <<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": The provided sequence contains a residue : "<<sequence.substr(i,1)<< " which is not a canonical RNA residue type." <<endl;
                   ErrorManager::instance.treatError();
               }
       }

   } else if (biopolymerType == BiopolymerType::DNA) {
       if (sequence.length() <1){
           ErrorManager::instance <<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": Your biopolymerType is "<<biopolymerType<<".  sequence must be of length > 0." <<endl;
           ErrorManager::instance.treatError();
       } 
       for (int i = 0; i < (int)sequence.length(); i++) {
           if (!
               letterIsDNA((sequence.substr(i,1)))
              ) {
                   ErrorManager::instance <<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": The provided sequence contains a residue : "<<sequence.substr(i,1)<< " which is not a canonical DNA residue type." <<endl;
                   ErrorManager::instance.treatError();
               }

       }
   } else if (biopolymerType == BiopolymerType::Protein ) {

       if (sequence.length() <1){
           ErrorManager::instance <<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": Your biopolymerType is "<<biopolymerType<<".  sequence must be of length > 0." <<endl;
           ErrorManager::instance.treatError();
       } 
       for (int i = 0; i < (int)sequence.length(); i++) {
           if (!
               letterIsProtein((sequence.substr(i,1)))
              ) {
                   ErrorManager::instance <<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": The provided sequence contains a residue : "<<sequence.substr(i,1)<< " which is not a canonical Protein residue type." <<endl;
                   ErrorManager::instance.treatError();
               }
       }
       if ((proteinCapping) && ((sequence.substr((strlen(sequence.c_str())-1),1)).compare("P") == 0 )) {
           ErrorManager::instance <<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" The last residue on a specified protein chain is a proline.  RNABuilder cannot add an end cap to a proline.  Please set \"proteinCapping False\" somewhere prior to specifying this chain.  If you can't do this then you might want to delete or change the residue type for this residue. "<<endl;
           ErrorManager::instance.treatError();
       }

   } else {
       ErrorManager::instance <<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": Your have requested an unsupported biopolymerType : " <<biopolymerType<<".  " <<endl;
       ErrorManager::instance.treatError();
   }
   return 0;
};
 
int BiopolymerClass::validateBiopolymerType () const {
        if (biopolymerType == BiopolymerType::RNA) {
           return 0;
        } else if (biopolymerType == BiopolymerType::DNA ) {
           return 0;
        } else if (biopolymerType == BiopolymerType::Protein ) {
           return 0;

        } else {
             ErrorManager::instance <<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": Your have requested an unsupported biopolymerType : " <<biopolymerType<<".  " <<endl;
             ErrorManager::instance.treatError();
        }

}

void BiopolymerClass::validateAtomInfoVector(){
    if (atomInfoVector.size() == 0) {
             ErrorManager::instance <<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": Your atomInfoVector has no elements! "  <<endl;
             ErrorManager::instance.treatError();

    } else {
        // Everything OK, do nothing
        //cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" For chain "<<getChainID()<<" Your atomInfoVector has "<<atomInfoVector.size()<<" elements"<<endl;     
    }
}


//Return the sequence corresponding to a list of ResidueID's from the current BiopolymerClass
String BiopolymerClass::getSequence(vector <ResidueID> & residueIDVector){
    String mySequence = "";
    for (int i = 0; i < residueIDVector.size(); i++){
        if ((i>0) && (getResidueIndex(residueIDVector[i-1]) > getResidueIndex(residueIDVector[i]))){std::cout <<__FILE__<<":"<<__LINE__<< " Error! Found two consecutive residues : "<< residueIDVector[i-1].outString() <<" and "<< residueIDVector[i].outString()  <<" which are not in increasing order."<<std::endl; exit(1); }
        mySequence +=  getResidueSingleLetterCode(residueIDVector[i]);
    }
    return mySequence;
}

void BiopolymerClass::validateMutation( Mutation myMutation) {
    myMutation.validate(); 
    if ( myMutation.getChain().compare(chainID) != 0) {
        ErrorManager::instance<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Bad mutant"<<endl; 
        ErrorManager::instance.treatError(); 
    }; 
    validateResidueID(myMutation.getResidue());
}




void BiopolymerClass::setSequence(String mySequence) {
    //cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<endl;
    this->sequence = mySequence; 
    //cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<endl;
    if (  biopolymerType == BiopolymerType::RNA) {
        //cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<endl;
        this->myBiopolymer = SimTK::RNA(mySequence,1);
        //cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<endl;
    } else if (  biopolymerType == BiopolymerType::DNA) {
        this->myBiopolymer = SimTK::DNA(mySequence,1);
    } else if (  biopolymerType == BiopolymerType::Protein) {
        this->myBiopolymer = SimTK::Protein(mySequence,BondMobility::Rigid,proteinCapping);
    }
    validateSequence();
    validateBiopolymerType();
    // renumberPdbResidues( this->firstResidueID ); //SCF
    //setPdbResidueNumbersFromResidueIDVector();
    // Note that this will not have the right PDB residue numbering. 
}


void BiopolymerClass::changeSequence(String myNewSequence) {
    String myOldSequence = getSequence();
    //String myNewSequence = myOldSequence;
    if (myNewSequence.length() != myOldSequence.length()) {
        ErrorManager::instance <<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" The new sequence : "<<myNewSequence<<" is of a different length than the old :"<<myOldSequence<<std::endl; 
        ErrorManager::instance.treatError();
    }
    // ResidueID myFirstResidueNumber = myOldBiopolymerClass.getFirstResidueID();
    //myNewSequence[myOldBiopolymerClass.getResidueIndex( myResidue) ] = *(mySubstitution.c_str()); // careful! getResidueIndex would potentially be wrong .. here we want the first letter of the sequence to correspond to position zero, with no regard to proteinCapping.
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": old sequence = "<<myOldSequence<<endl;
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": new sequence = "<<myNewSequence<<endl;
    setSequence(myNewSequence);  // Note that this will not have the right PDB residue numbering. Hence the next line:
    setPdbResidueNumbersFromResidueIDVector();
}

/*void BiopolymerClass::renameChain(String newChainID) {
    setChainID(newChainID);
    myBiopolymer.setPdbChainID(newChainID);
}*/

void BiopolymerClass::setOriginalSequence(String mySequence) {
    originalSequence = mySequence; 
    //validateSequence();
}

void BiopolymerClass::setBiopolymerType(String myBiopolymerType) {
    if ( myBiopolymerType.compare ("RNA") == 0) {
        biopolymerType = BiopolymerType::RNA ;}
    else if (myBiopolymerType.compare("DNA") == 0) {
        biopolymerType = BiopolymerType::DNA;}
    else if (myBiopolymerType.compare("Protein") == 0) {
        biopolymerType = BiopolymerType::Protein;}
    else {
             ErrorManager::instance <<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": Your have requested an unsupported biopolymerType : " << myBiopolymerType<<".  " <<endl;
             ErrorManager::instance.treatError();
    }
    validateBiopolymerType();
}

void BiopolymerClass::setBiopolymerType(BiopolymerType::BiopolymerTypeEnum myBiopolymerType) {
    biopolymerType = myBiopolymerType;
    validateBiopolymerType();
}
int  BiopolymerClass::validateProteinCapping () {
    if (proteinCapping && (biopolymerType == BiopolymerType::RNA)) {
        ErrorManager::instance <<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<"Error! You have set proteinCapping TRUE for an RNA.  This is not allowed. "<<endl;
        ErrorManager::instance.treatError();
    }
    else if ((!proteinCapping) && (biopolymerType == BiopolymerType::RNA)) {} // this is OK
    else if (proteinCapping && (biopolymerType == BiopolymerType::DNA)) {
        ErrorManager::instance <<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<"Error! You have set proteinCapping TRUE for an DNA.  This is not allowed. "<<endl;
        ErrorManager::instance.treatError();
    } 
    else if ((!proteinCapping) && (biopolymerType == BiopolymerType::DNA)) {} // this is OK
    else if (biopolymerType == BiopolymerType::Protein) {} // this is always OK.
    else if (proteinCapping && (biopolymerType == BiopolymerType::Unassigned)) {
        ErrorManager::instance <<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<"Error! You have set proteinCapping TRUE for a bipolymer of type Unassigned.  This is not allowed. "<<endl;
        ErrorManager::instance.treatError();
    }
    else if ((!proteinCapping) && (biopolymerType == BiopolymerType::Unassigned)) {} // this is OK.
    else {
        ErrorManager::instance <<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<"Error! You have tried to set proteinCapping for biopolymerType "<<biopolymerType<<", which is not of a type which can have proteinCapping set. "<<endl;
        ErrorManager::instance.treatError();

    }
    return 0;
}

void BiopolymerClass::setProteinCapping(bool myProteinCapping) {
    proteinCapping = myProteinCapping;
    validateProteinCapping();
}

void BiopolymerClass::setFirstResidueMobilizerType(String myFirstResidueMobilizerType){
    if ((myFirstResidueMobilizerType.compare("Weld") == 0 ) || (myFirstResidueMobilizerType.compare("Free") == 0)) 
        firstResidueMobilizerType = myFirstResidueMobilizerType;
    else {
        ErrorManager::instance <<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<"You have tried to specify an invalid firstResidueMobilizerType = "<<myFirstResidueMobilizerType<<endl; ErrorManager::instance.treatError();
    }
}

String BiopolymerClass::getFirstResidueMobilizerType(){
    return firstResidueMobilizerType;
}

BiopolymerClass::BiopolymerClass() {
    clear();
    std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" sizeof(PdbStructure) = "<< sizeof(PdbStructure) <<std::endl;
    std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" sizeof(myBiopolymer) = "<< sizeof(myBiopolymer) <<std::endl;
}

BiopolymerClass::BiopolymerClass(String mySequence, String myChainID, ResidueID myFirstResidueNumber, String myBiopolymerType, bool proteinCapping , bool useNACappingHydroxyls ){
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<endl;
    clear();
    setBiopolymerType(     myBiopolymerType);
    //cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<endl;
    setProteinCapping (proteinCapping); 
    //cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<endl;
    setChainID(myChainID);
    setSequence(mySequence);  // Note that this will not have the right PDB residue numbering. Hence the next line:
    renumberPdbResidues( myFirstResidueNumber );
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<endl;
    //setChainID(myChainID);
    pdbFileName = "";
    //cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<endl;
    loadFromPdb = false;
    
    validateChainID();
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<endl;
    if (  biopolymerType == BiopolymerType::RNA) {
	// useNACappingHydroxyls , when true (default) replaces the 5' phosphorus with an H5T.      
        myBiopolymer = SimTK::RNA(mySequence,useNACappingHydroxyls);
    } else if (  biopolymerType == BiopolymerType::DNA) {
        myBiopolymer = SimTK::DNA(mySequence,useNACappingHydroxyls);
    } else if (  biopolymerType == BiopolymerType::Protein) {
        //cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" "<<mySequence<<", "<<proteinCapping<<endl;
        myBiopolymer = SimTK::Protein   (mySequence,BondMobility::Rigid,proteinCapping);
        //myBiopolymer = SimTK::Protein   (mySequence,BondMobility::Rigid,proteinCapping);
        //cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<endl;
    } else { 
        ErrorManager::instance <<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<"Error! You have tried to make a Biopolymer of an unsupported type: "<< myBiopolymerType<<endl;
        ErrorManager::instance.treatError();
    }
    validateSequence();
    validateBiopolymerType();
    validateProteinCapping();
    renumberPdbResidues((myFirstResidueNumber ) );
    //cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<endl;
}


void BiopolymerClass::setPdbResidueNumbersFromResidueIDVector() {
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Inside setPdbResidueNumbersFromResidueIDVector() for chain "<< getChainID()<< endl;
    ResidueID myFirstResidueID = getFirstResidueID() ;
    validateProteinCapping();
    int numResidues = -1;
    int capping = -1111;
    if (proteinCapping && (biopolymerType == BiopolymerType::Protein)) {
        capping = 1;
        ResidueID myResidueID = getFirstResidueID();
        myBiopolymer.updResidue(ResidueInfo::Index(0)).setPdbResidueNumber((decrementResidueID(myResidueID)).getResidueNumber());
        myBiopolymer.updResidue(ResidueInfo::Index(0)).setPdbInsertionCode((decrementResidueID(myResidueID)).getInsertionCode());
    } else if ((biopolymerType == BiopolymerType::RNA) || (biopolymerType == BiopolymerType::Protein) || (biopolymerType == BiopolymerType::DNA) ) {
        capping = 0;
    } else { 
        ErrorManager::instance <<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<"Error! You have tried to renumber a Biopolymer of an unsupported type: "<<   biopolymerType<<endl;
        ErrorManager::instance.treatError();
    }
    for (int i = capping ; i < (residueIDVector.size() + capping) ; i++) {
        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" ResidueInfo::Index(i) = "<<ResidueInfo::Index(i)<<" and residueIDVector[i-capping] = "<<residueIDVector[i-capping].outString()<<endl;
	myBiopolymer.updResidue(ResidueInfo::Index(i)).setPdbResidueNumber(residueIDVector[i-capping].getResidueNumber()); // Remember, residue with index 0 will be the ACE N-terminal cap, if this is a protein and proteinCapping == true
	myBiopolymer.updResidue(ResidueInfo::Index(i)).setPdbInsertionCode(residueIDVector[i-capping].getInsertionCode()); 
    } 
    if ( getFirstResidueID() != residueIDVector[0]) {
	ErrorManager::instance <<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" getFirstResidueID() returned "<< getFirstResidueID().outString() <<" while residueIDVector[0] returned "<<residueIDVector[0].outString()<<" . This makes no sense! "<<endl;
	ErrorManager::instance.treatError();
    }
    if ( getLastResidueID() != residueIDVector[residueIDVector.size()-1]) {
	ErrorManager::instance <<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" getLastResidueID() returned "<< getLastResidueID().outString() <<" while residueIDVector[residueIDVector.size()-1] returned "<<residueIDVector[residueIDVector.size()-1].outString()<<" . This makes no sense! "<<endl;
	ErrorManager::instance.treatError();
    }
    
}
 
// This renumbers myBiopolymer . You still need to setResidueIDsAndInsertionCodesFromBiopolymer
void BiopolymerClass::renumberPdbResidues(ResidueID firstResidueID ) {
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" firstResidueID = >"<<firstResidueID.outString()<<"< "<<endl; 
    this->firstResidueID = firstResidueID;
    validateProteinCapping();
    int myNewResidueNumberWithoutInsertionCode = firstResidueID.getResidueNumber() ;
    if (proteinCapping && (biopolymerType == BiopolymerType::Protein)) {
        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<endl;
        myBiopolymer.renumberPdbResidues(decrementResidueID(firstResidueID).getResidueNumber() ); 
        //myNewResidueNumberWithoutInsertionCode -= 1; // proteins have an ACE residue added to the N-terminus.  Therefore that residue must be given the first residue number minus one!
    } else if ((biopolymerType == BiopolymerType::RNA) || (biopolymerType == BiopolymerType::Protein) || (biopolymerType == BiopolymerType::DNA) ) {
        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<endl;
        // myNewResidueNumberWithoutInsertionCode = (firstResidueID.getResidueNumber() ); // proteins have an ACE residue added to the N-terminus.  Therefore that residue must be given the first residue number minus one!
        myBiopolymer.renumberPdbResidues(firstResidueID.getResidueNumber()); // this is now done below.
    } else { 
        ErrorManager::instance <<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<"Error! You have tried to renumber a Biopolymer of an unsupported type: "<<   biopolymerType<<endl;
        ErrorManager::instance.treatError();
    }
    //myBiopolymer.renumberPdbResidues((myNewResidueNumberWithoutInsertionCode));
    //residueIDVector.clear();
    //setResidueIDsAndInsertionCodesFromBiopolymer(myBiopolymer,proteinCapping);//loadResidueIDVector();
}

ResidueID BiopolymerClass::getResidueID(const int residueIndex)  {
    if (residueIDVector.size() > 0) {
        return residueIDVector[residueIndex];
    }
    //else
    ErrorManager::instance<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Why are you even contemplating this ridiculous means of getting a residueID?  Just use the residueIDVector!"<<endl;
    ErrorManager::instance.treatError();

    ResidueID myResidueID;
    myResidueID.setResidueNumber(myBiopolymer.getResidue(ResidueInfo::Index(residueIndex)).getPdbResidueNumber());
    myResidueID.setInsertionCode(myBiopolymer.getResidue(ResidueInfo::Index(residueIndex)).getPdbInsertionCode());
    return myResidueID;
}

int  BiopolymerClass::matchCoordinates(String inputFileName, 
                                       bool matchExact, bool matchIdealized,
                                       const bool matchOptimize ,  
                                       bool matchHydrogenAtomLocations, 
                                       bool matchPurineN1AtomLocations,
                                       bool guessCoordinates ,  
                                       double matchingMinimizerTolerance, 
                                       double myPlanarityThreshold   // this parameter sets the out-of-planarity tolerance for identifying planar bonds.  Units: radians.

    ) {
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" about to match chain \""<< getChainID()<<"\" to file name : "<<inputFileName<<endl;
    //if(pdbStructure == NULL)
    {
        PdbStructure myPdbStructure;
        
        //============================================ Read in PDB or CIF
        if ( inputFileName.length() > 4 )
        {
            if ( inputFileName.substr ( inputFileName.length() - 4, inputFileName.length() - 1) == ".pdb" )
            {
                std::ifstream inputFile               ( inputFileName.c_str(), ifstream::in );
                
                if ( !inputFile.good() )
                {
                    ErrorManager::instance << "!!! Error !!! The file " << inputFileName << " could not be opened. If this is not the file you wanted to open, please supply the requested file name after the loadSequencesFromPdb command. Note that the supported file extensions currently are \".pdb\", \".cif\" and \".cif.gz\"." << std::endl;
                    ErrorManager::instance.treatError ( );
                }
                else
                {
                    myPdbStructure                    = PdbStructure (inputFile);
                }
            }
            else if ( inputFileName.substr ( inputFileName.length() - 4, inputFileName.length() - 1) == ".cif" )
            {
#ifdef GEMMI_USAGE
                std::ifstream testOpen                ( inputFileName.c_str() );
                if ( testOpen.good() )
                {
                    myPdbStructure                    = PdbStructure (inputFileName);
                }
                else
                {
                    std::string pdbFileHlp            = inputFileName;
                    pdbFileHlp.append                 ( ".gz" );
                    
                    std::ifstream testOpen2           ( pdbFileHlp.c_str() );
                    if ( testOpen2.good() )
                    {
                        myPdbStructure                = PdbStructure ( pdbFileHlp );
                    }
                    else
                    {
                        ErrorManager::instance << "!!! Error !!! The file " << inputFileName << " could not be opened. If this is not the file you wanted to open, please supply the requested file name after the loadSequencesFromPdb command. Note that the supported file extensions currently are \".pdb\", \".cif\" and \".cif.gz\"." << std::endl;
                        ErrorManager::instance.treatError ( );
                    }
                    testOpen2.close                   ( );
                }
                testOpen.close                        ( );
#else
                ErrorManager::instance << "!!! Error !!! MMB was not compiled with the Gemmi library required for mmCIF support. Cannot proceed, if you want to use mmCIF files, please re-compile with the Gemmi library option allowed." << std::endl;
                ErrorManager::instance.treatError     ( );
#endif
            }
            else if ( inputFileName.length() > 7 )
            {
                if ( inputFileName.substr ( inputFileName.length() - 7, inputFileName.length() - 1) == ".cif.gz" )
                {
#ifdef GEMMI_USAGE
                myPdbStructure                        = PdbStructure (inputFileName);
#else
                ErrorManager::instance << "!!! Error !!! MMB was not compiled with the Gemmi library required for mmCIF support. Cannot proceed, if you want to use mmCIF files, please re-compile with the Gemmi library option allowed." << std::endl;
                ErrorManager::instance.treatError     ( );
#endif
                }
                else
                {
                    ErrorManager::instance << "!!! Error !!! The file " << inputFileName << " could not be opened. If this is not the file you wanted to open, please supply the requested file name after the loadSequencesFromPdb command. Note that the supported file extensions currently are \".pdb\", \".cif\" and \".cif.gz\"." << std::endl;
                    ErrorManager::instance.treatError ( );
                }
            }
            else
            {
                ErrorManager::instance << "!!! Error !!! The file " << inputFileName << " could not be opened. If this is not the file you wanted to open, please supply the requested file name after the loadSequencesFromPdb command. Note that the supported file extensions currently are \".pdb\", \".cif\" and \".cif.gz\"." << std::endl;
                ErrorManager::instance.treatError     ( );
            }
        }
        else
        {
            ErrorManager::instance << "!!! Error !!! The file " << inputFileName << " could not be opened. If this is not the file you wanted to open, please supply the requested file name after the loadSequencesFromPdb command. Note that the supported file extensions currently are \".pdb\", \".cif\" and \".cif.gz\"." << std::endl;
            ErrorManager::instance.treatError         ( );
        }
        
        return matchCoordinates( myPdbStructure, matchExact, matchIdealized, matchOptimize,
                         matchHydrogenAtomLocations, matchPurineN1AtomLocations,
                         guessCoordinates, matchingMinimizerTolerance, myPlanarityThreshold);
    }
    /*
    else
    {
        cout << "using cached PdbStructure" << endl;
        return matchCoordinates(*pdbStructure, matchExact, matchIdealized, matchOptimize,
                         matchHydrogenAtomLocations, matchPurineN1AtomLocations,
                         guessCoordinates, matchingMinimizerTolerance, myPlanarityThreshold);
    }*/
}

int  BiopolymerClass::matchCoordinates(istream & inputFile,
                                       bool matchExact, bool matchIdealized,
                                       const bool matchOptimize ,  
                                       bool matchHydrogenAtomLocations, 
                                       bool matchPurineN1AtomLocations,
                                       bool guessCoordinates ,  
                                       double matchingMinimizerTolerance, 
                                       double myPlanarityThreshold   // this parameter sets the out-of-planarity tolerance for identifying planar bonds.  Units: radians.

    ) {
    PdbStructure myPdbStructure(inputFile);
    // inputFile.close();
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" PdbStructure done for chain " << getChainID() << endl;
    return matchCoordinates(myPdbStructure, matchExact, matchIdealized, matchOptimize,
                     matchHydrogenAtomLocations, matchPurineN1AtomLocations,
                     guessCoordinates, matchingMinimizerTolerance, myPlanarityThreshold);
    //cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Coordinates matched for chain " << getChainID() << endl;
}


// Return zero if successful, 1 for error.
int  BiopolymerClass::matchCoordinates(const PdbStructure & myPdbStructure, 
                                       bool matchExact, bool matchIdealized,
                                       const bool matchOptimize ,  
                                       bool matchHydrogenAtomLocations, 
                                       bool matchPurineN1AtomLocations,
                                       bool guessCoordinates ,  
                                       double matchingMinimizerTolerance, 
                                       double myPlanarityThreshold   // this parameter sets the out-of-planarity tolerance for identifying planar bonds.  Units: radians.

    ) {
    double maxObservedSinePlaneDeviation = 0;
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<endl;
    Compound::AtomTargetLocations biopolymerAtomTargets = myBiopolymer.createAtomTargets(myPdbStructure,guessCoordinates); 

    bool matchProteinCarboxylOxygenLocations = false;
    bool matchNucleotideSideGroups = false;

    Compound::AtomTargetLocations::iterator biopolymerAtomTargetsIterator = biopolymerAtomTargets.begin();
    map<Compound::AtomIndex, Vec3>::iterator it;
    map<Compound::AtomIndex, Vec3>::iterator next;
    next = biopolymerAtomTargets.begin();
    while (next != biopolymerAtomTargets.end())
    {
        it = next;
        Compound::AtomIndex m = (*it).first;
        Element myAtomElement = myBiopolymer.getAtomElement(m);
        Compound::AtomName myAtomName = myBiopolymer.getAtomName(m);
        size_t pos = myAtomName.find("/");
        String myAtomNameSubstr = myAtomName.substr(pos);
        String myAtomNameOnly = myAtomName.substr(pos+1);
        ResidueInfo::Index  myResidueIndex (    atoi(myAtomName.substr(0,pos).c_str()) );
        next++;

        if( !matchHydrogenAtomLocations && ((myAtomElement.getName()).compare("hydrogen") == 0) )
        {
            biopolymerAtomTargets.erase(it);
            continue;
        } 
	#ifdef _DEBUG_FLAGS_ON_
	    //cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<endl;
	    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" biopolymerType = >"<<biopolymerType<<"< compare to protein: "<<BiopolymerType::Protein<<endl;
	    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" matchProteinCarboxylOxygenLocations = >"<<matchProteinCarboxylOxygenLocations<<"< "<<endl;
	    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" myAtomName = >"<<myAtomName<<"< "<<endl;
	    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" myAtomNameSubstr = >"<<myAtomNameSubstr<<"< "<<endl;
	    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" myResidueIndex = >"<<myResidueIndex<<"< "<<endl;
	    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" getResidue = >"<<getResidueID(myResidueIndex).outString()<<"< "<<endl;
	#endif
        if( (biopolymerType==BiopolymerType::Protein) && (!matchProteinCarboxylOxygenLocations) && ((myAtomNameSubstr).compare("/O")==0) )
        {
	    #ifdef _DEBUG_FLAGS_ON_
		cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" erasing atom from biopolymerAtomTargets, with myAtomNameSubstr = >"<<myAtomNameSubstr<<"< "<<endl;
	    #endif
            biopolymerAtomTargets.erase(it);
            continue;
        } 
        for (int i = 0; i < ignoreAtomPositionVector.size(); i++) {
            //MMBAtomInfo myMMBAtomInfo = ignoreAtomPositionVector(i);
            //if (myMMBAtomInfo.getResidueIndex() == myResidueIndex) (
            if (ignoreAtomPositionVector[i].getResidueIndex() == myResidueIndex) {
                if (ignoreAtomPositionVector[i].getAtomName().compare(myAtomNameOnly) == 0) {
		    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" erasing atom from biopolymerAtomTargets, with myAtomNameSubstr = >"<<myAtomNameSubstr<<"< "<<endl;
		    biopolymerAtomTargets.erase(it);
		    continue;
                }
            }     
        }
        if( (((biopolymerType == BiopolymerType::RNA) || (biopolymerType == BiopolymerType::DNA)) && (!(matchPurineN1AtomLocations ))) && residueIsPurine(myResidueIndex,getOriginalSequence()) )
        {
            // have to check residue type of ORIGINAL sequence, otherwise the new residue might be a pyrimidine and we would not delete the N1, defeating the whole point.
            if(myAtomNameSubstr.compare("/N1") == 0 ) 
            {
                // cout<<"-DEL!-"<<flush;
                biopolymerAtomTargets.erase(it);
                continue;
            }
        }
        if( (((biopolymerType == BiopolymerType::RNA)|| (biopolymerType == BiopolymerType::DNA )) && (! matchNucleotideSideGroups )) && ((myAtomNameSubstr.compare("/N2") == 0) ||
                   (myAtomNameSubstr.compare("/N6") == 0)) )
        {
            biopolymerAtomTargets.erase(it);
            continue;
        }
    }
    if (matchExact) {
            // low tolerance breaks planarity just about everywhere
            //std::ofstream tempStream(String("match.382.pdb").c_str(),ios_base::out);
            //myBiopolymer.writeDefaultPdb(tempStream,Transform(Vec3(0)));
            cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<"about to myBiopolymer.matchDefaultAtomChirality for chain  "<<getChainID()<<endl;
            //int matchDefaultAtomChiralityReturnValue = 0;

            myBiopolymer.matchDefaultAtomChirality(biopolymerAtomTargets, maxObservedSinePlaneDeviation,  myPlanarityThreshold, false);

            cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" maxObservedSinePlaneDeviation = "<<maxObservedSinePlaneDeviation<<std::endl;
            if (maxObservedSinePlaneDeviation > .2 ) {
                cout<<__FILE__<<":"<<__LINE__<<" maxObservedSinePlaneDeviation = "<<maxObservedSinePlaneDeviation<<" . Warning! This might be too high!"<<std::endl;  
                //return 1;
            } // We will use this as our error threshold
            
            cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<"done with myBiopolymer.matchDefaultAtomChirality"<<endl;
            //std::ofstream tempStream2(String("match.385.pdb").c_str(),ios_base::out);
            //myBiopolymer.writeDefaultPdb(tempStream2,Transform(Vec3(0)));
            myBiopolymer.matchDefaultBondLengths(biopolymerAtomTargets);
            cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<"done with myBiopolymer.matchDefaultBondLengths"<<endl;
            //std::ofstream tempStream3(String("match.388.pdb").c_str(),ios_base::out);
            //myBiopolymer.writeDefaultPdb(tempStream3,Transform(Vec3(0)));
            myBiopolymer.matchDefaultBondAngles(biopolymerAtomTargets);
            cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<"done with myBiopolymer.matchDefaultBondAngles   "<<endl;
            //std::ofstream tempStream6(String("match.393.pdb").c_str(),ios_base::out);
            //myBiopolymer.writeDefaultPdb(tempStream6,Transform(Vec3(0)));
            // Set dihedral angles even when bonded atoms are planar
            myBiopolymer.matchDefaultDihedralAngles(biopolymerAtomTargets, Compound::DistortPlanarBonds);//KeepPlanarBonds); //was DistortPlanarBonds

            cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<"done with myBiopolymer.matchDefaultDihedralAngles   "<<endl;
            #ifdef _DEBUG_FLAGS_ON_
                cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<endl;
		//std::ofstream tempStream4(String("match.397.pdb").c_str(),ios_base::out);
                cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<endl;
		//#myBiopolymer.writeDefaultPdb(tempStream4,Transform(Vec3(0)));
		cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" just wrote match.397.pdb"<<endl;
            #endif
            cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" starting matchDefaultTopLevelTransform:"<<endl;
            myBiopolymer.matchDefaultTopLevelTransform(biopolymerAtomTargets);
            cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" confirming transform: "<<endl;
            cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" myBiopolymer.getTopLevelTransform() = "<<  myBiopolymer.getTopLevelTransform()<<endl;
            cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" about to writeDefaultPdb:"<<endl;
            std::ofstream tempStream5(String("match.400.pdb").c_str(),ios_base::out);
            myBiopolymer.writeDefaultPdb(tempStream5,Transform(Vec3(0)));
            cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Now writing with top level transform provided:"<<endl;
            std::ofstream tempStream6(String("match.401.pdb").c_str(),ios_base::out);
            myBiopolymer.writeDefaultPdb(tempStream6,myBiopolymer.getTopLevelTransform()  );
            cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<endl;
    }

    if (matchIdealized) {
        rigidifyTargetedBonds(biopolymerAtomTargets);
        PdbAtom::setWriteFullPrecisionLocation(true); // PDB stucts are used to set the default coordinates in the final steps of this method.  Let's use higher precision.
        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" About to optimize chain "<<getChainID()<<" using ObservedPointFitter, with matchingMinimizerTolerance = "<<matchingMinimizerTolerance<<".  If this fails to converge, try increasing this parameter. If it converges but is not sufficiently close to your input structure file, try decreasing the parameter."<< endl;
        // third parameter defaults to useObservedPointFitter = true. This means that LocalEnergyMinimizer will NOT be called.       
        bool myUseObservedPointFitter = true ;
        myBiopolymer.matchDefaultConfiguration(biopolymerAtomTargets,Compound::Match_Idealized, myUseObservedPointFitter, matchingMinimizerTolerance ); 
        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<endl;
    } 
    // std::ofstream tempStream6(String("match.416.pdb").c_str(),ios_base::out);
    // myBiopolymer.writeDefaultPdb(tempStream6,Transform(Vec3(0)));
    if (matchOptimize){
 // found that the following is redundant:
        bool myUseObservedPointFitter = false;
        myUseObservedPointFitter = true ;
        
        myBiopolymer.matchDefaultConfiguration(biopolymerAtomTargets,Compound::Match_Idealized, myUseObservedPointFitter, matchingMinimizerTolerance ); 

        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" A small tolerance means more accuracy but takes longer .. the default is typically accurate enough, but you might experiment with a larger one."<<endl;
    } 
    // std::ofstream tempStream7(String("match.428.pdb").c_str(),ios_base::out);
    // myBiopolymer.writeDefaultPdb(tempStream7,Transform(Vec3(0)));
    if (getRenumberPdbResidues()){ // Once matching is done, it should be safe to renumber
        //SCF
        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<endl;
        renumberPdbResidues(ResidueID("1"));
    }
    return 0; // If we got this far, all is well.
}

void BiopolymerClass::rigidifyTargetedBonds(Compound::AtomTargetLocations & biopolymerAtomTargets) {

    map<Compound::AtomIndex, Vec3>::iterator it;
    map<Compound::AtomIndex, Vec3>::iterator next;
    next = biopolymerAtomTargets.begin();
    while (next != biopolymerAtomTargets.end())
    {
       it = next;
       Compound::AtomIndex m = (*it).first;
       Element myAtomElement = myBiopolymer.getAtomElement(m);
       next++;
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<(myAtomElement.getName())<<endl;
       cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" "<<(myAtomElement.getName())<<", "<<  biopolymerAtomTargets[m]  <<endl;
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" "<<m<<","<<(myBiopolymer.getAtomName(m))<<endl;
    }
/*
    ResidueInfo& residue = updResidue(r);

    // Create a list of all atoms in the residue
    std::set<Compound::AtomIndex> residueAtoms;
    for (ResidueInfo::AtomIndex a(0); a < residue.getNumAtoms(); ++a)
        residueAtoms.insert(residue.getAtomIndex(a));

    // Set mobility on bonds that are within residue
    // atoms
    std::set<Compound::AtomIndex>::const_iterator aI;
    for (aI = residueAtoms.begin(); aI !=residueAtoms.end(); ++aI) {
        const CompoundAtom& atom = getImpl().getAtom(*aI);
        // atoms->bondCenters
        for (CompoundAtom::BondCenterIndex bc(0); bc < atom.getNumBondCenters(); ++bc) {
            BondCenterInfo::AtomKey key(*aI, bc);
            const BondCenterInfo& bondCenter = getImpl().getBondCenterInfo(key);
            if (! bondCenter.isBonded()) continue; // skip unbonded centers
            const BondCenterInfo& otherBondCenter =
                getImpl().getBondCenterInfo(bondCenter.getBondPartnerBondCenterIndex());
            Compound::AtomIndex otherAtomIndex = otherBondCenter.getAtomIndex();
            // skip bonds that go outside of this residue
            if (residueAtoms.find(otherAtomIndex) == residueAtoms.end()) continue;
            BondInfo& bondInfo = updImpl().updBondInfo(bondCenter.getBondIndex());
            updImpl().updBond(bondInfo).setMobility(mobility);
        }
    }

    //return *this;


        for (int q=0;q<(int)baseOperationVector.size();q++)  

           if (((baseOperationVector[q]).BasePairIsTwoTransformForce).compare("mobilizer") == 0){
              cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Setting mobilizer of type "<<(baseOperationVector[q]).FirstBPEdge<<" for chain "<<(baseOperationVector[q]).FirstBPChain<<" from residue "<<(baseOperationVector[q]).FirstBPResidue.outString()<<" to residue "<<(baseOperationVector[q]).SecondBPResidue.outString()<<endl;
               MobilizerStretch dummyMobilizerStretch;
               BondMobility::Mobility myBondMobility = dummyMobilizerStretch.setBondMobility(baseOperationVector[q].FirstBPEdge ) ;
               BiopolymerClass & myBiopolymerClass ( updBiopolymerClass((baseOperationVector[q]).FirstBPChain));
                   if (biopolymerClassMap[baseOperationVector[q].FirstBPChain].biopolymerType == BiopolymerType::RNA){
                        (static_cast<RNA&>( myBiopolymerClass.myBiopolymer)).setRNABondMobility(myBondMobility,
                            SimTK::ResidueInfo::Index (myBiopolymerClass.getResidueIndex((baseOperationVector[q]).FirstBPResidue)), 
                            SimTK::ResidueInfo::Index (myBiopolymerClass.getResidueIndex((baseOperationVector[q]).SecondBPResidue))); 
                   } else if (biopolymerClassMap[baseOperationVector[q].FirstBPChain].biopolymerType == BiopolymerType::DNA){
                        //ErrorManager::instance <<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" At this time DNA is not a supported biopolymer type. "<< baseOperationVector[q].FirstBPEdge <<". "<<endl;
                        //ErrorManager::instance.treatError();
                        (static_cast<DNA&>( myBiopolymerClass.myBiopolymer)).setDNABondMobility(myBondMobility,
                            SimTK::ResidueInfo::Index (myBiopolymerClass.getResidueIndex ((baseOperationVector[q]).FirstBPResidue)), 
                            SimTK::ResidueInfo::Index (myBiopolymerClass.getResidueIndex ((baseOperationVector[q]).SecondBPResidue))); 
                   } else if (biopolymerClassMap[(baseOperationVector[q]).FirstBPChain].biopolymerType == BiopolymerType::Protein) {
                       myBiopolymerClass.setProteinBondMobility(
                           myBondMobility,
                           (baseOperationVector[q]).FirstBPResidue,
                           (baseOperationVector[q]).SecondBPResidue
                           );
                   } else exit(1);
}
*/
}

void BiopolymerClass::setSingleBondMobility(ResidueID residueID1,  String atomName1,ResidueID residueID2, String atomName2, String mobilityString ) {
    SimTK::BondMobility::Mobility myBondMobility ;  
    if (mobilityString.compare("Rigid") ==0) { myBondMobility = BondMobility::Rigid;}
    else if (mobilityString.compare("Torsion") ==0) { myBondMobility = BondMobility::Torsion;}
    else if (mobilityString.compare("Free") ==0) { myBondMobility = BondMobility::Free;}
    else {ErrorManager::instance <<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Unrecognied bond mobility: "<<mobilityString<<endl; ErrorManager::instance.treatError();}
    myBiopolymer.setBondMobility(myBondMobility,atomPathString(residueID1,atomName1),atomPathString(residueID2,atomName2 ));         
}
    /**
     * \brief Set chain ID, renumber residues, match coordinates, and adopt the compound.
     *
     */
int  BiopolymerClass::initializeBiopolymer(CompoundSystem & system, 
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
                                          ) 
{
    //int returnValue = 0;
    if (biopolymerType == BiopolymerType::Protein) {
        setProteinCapping (myProteinCapping);
    }
    else{
        setProteinCapping (false);                           
    }
    //setRenumberPdbResidues(0); // Assume we will not be renumbering. This may change later depending on user input.
    myBiopolymer.setPdbChainId(((chainID.c_str())));
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Just issued myBiopolymer.setPdbChainId("<<chainID.c_str()<<") based on chainID = >"<<chainID<<"<"<<endl;
    Vec3 initialDisplacementVec3 = Vec3(0); // essential to initialize, since displacementVector may be empty
    Rotation myRotation;
    myRotation.setRotationToIdentityMatrix ();

    for (int i = 0; i < (int)displacementVector.size(); i++){
        if (displacementVector[i].chain.compare(getChainID()) == 0) {
            cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Displacement vector index "<<i<<" chain "<<displacementVector[i].chain<<" matches current chain "<<getChainID()<<".  Applying displacement from input structure file of : "<< displacementVector[i].displacement <<" Å "<<endl;
            myRotation              = displacementVector[i].rotation;
            initialDisplacementVec3 = displacementVector[i].displacement;
            break;
        }
        initialDisplacementVec3 = Vec3(0); // this is redundant 
    }
    if (this->loadFromPdb) {
        //returnValue = 
        if (matchCoordinates(this->pdbFileName, matchExact, matchIdealized,matchOptimize ,matchHydrogenAtomLocations,matchPurineN1AtomLocations, guessCoordinates, matchingMinimizerTolerance,myPlanarityThreshold )) {
            cout<<__FILE__<<":"<<__LINE__<<" Warning: Returned an error from matchCoordinates"<<std::endl;
            //return 1;
        }

        ////////////////////
        SecondaryStructureStretch mySecondaryStructureStretch  ;
        for (int i = 0; i <   secondaryStructureStretchVector.size(); i++) 
        {
            mySecondaryStructureStretch =   secondaryStructureStretchVector[i]; 
            if (mySecondaryStructureStretch.getChain().compare( getChainID()) == 0) 
            {     
                cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Applying secondary structure default phi, psi angles for stretch : "<<i<<endl;
                mySecondaryStructureStretch.printStretch();
                if (mySecondaryStructureStretch.getSecondaryStructureType() == Alpha) 
                {
                    setAlphaHelicalDefaultBackboneAngles(mySecondaryStructureStretch.getStartResidue(), mySecondaryStructureStretch.getEndResidue());
                } else if (mySecondaryStructureStretch.getSecondaryStructureType() == ParallelBeta) 
                {
                    setParallelBetaSheetDefaultBackboneAngles(mySecondaryStructureStretch.getStartResidue(), mySecondaryStructureStretch.getEndResidue());
                } else if (mySecondaryStructureStretch.getSecondaryStructureType() == AntiParallelBeta) 
                {
                    setAntiParallelBetaSheetDefaultBackboneAngles(mySecondaryStructureStretch.getStartResidue(), mySecondaryStructureStretch.getEndResidue());
                } else 
                {
                    ErrorManager::instance <<"["<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<"] : At the moment we can only enforce secondary structure types Alpha, ParallelBeta, and AntiParallelBeta."<<endl;
                    ErrorManager::instance.treatError();
                }
            } // of if chain
        }
        ////////////////////
    }
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<endl;   
    if (getRenumberPdbResidues()){
        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<endl;   
        renumberPdbResidues(ResidueID("1")); 
        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<endl;   
        // second argument evaluates to True if proteinCapping is requested, and this is a protein 
        setResidueIDsAndInsertionCodesFromBiopolymer(myBiopolymer, (getProteinCapping() && (biopolymerType == BiopolymerType::Protein) ));
        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<endl;   
    }
 //895     if (biopolymerType == BiopolymerType::Protein) {
 //896         setProteinCapping (myProteinCapping);


    if (this->loadFromPdb) {
        cout << __FILE__ << " " << __FUNCTION__ << " :" << firstResidueID.getResidueNumber() << endl; 
        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Adopting chain "<<getChainID()<<" with displacement from input structure file of : "<<initialDisplacementVec3<<" Å "<<getSequence()<<endl;
        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Current rotation : "<<myRotation<<endl;  
        system.adoptCompound(myBiopolymer ,Transform(myRotation, (initialDisplacementVec3/1)) );} // used to convert to nm, now using nm directly. 
    else {
        system.adoptCompound(myBiopolymer ,Vec3(biopolymerClassIndex,biopolymerClassIndex,biopolymerClassIndex  )*initialSeparation/1);  // used to convert to nm, now using nm directly

    }



    validateResidueNumbersAndInsertionCodes();
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<std::endl;
    return 0;//returnValue;
}
int BiopolymerClass::getChainLength() {
    return sequence.length();
};

size_t BiopolymerClass::getNumAtoms() {
    return myBiopolymer.getNumAtoms();
}

ResidueID BiopolymerClass::getFirstResidueID(){
    int myResidueNumber = myBiopolymer.getResidue(ResidueInfo::Index(int(proteinCapping))).getPdbResidueNumber();    
    int myInsertionCode = myBiopolymer.getResidue(ResidueInfo::Index(int(proteinCapping))).getPdbInsertionCode();    
    return ResidueID(myResidueNumber, myInsertionCode);
};
ResidueID BiopolymerClass::getLastResidueID(){
    int myResidueNumber = myBiopolymer.getResidue(ResidueInfo::Index(myBiopolymer.getNumResidues() - int(proteinCapping) -1)).getPdbResidueNumber();    
    int myInsertionCode = myBiopolymer.getResidue  ( ResidueInfo::Index( myBiopolymer.getNumResidues() - int(proteinCapping) -1)).getPdbInsertionCode();    
    return ResidueID(myResidueNumber, myInsertionCode);
};
BiopolymerType::BiopolymerTypeEnum BiopolymerClass::getBiopolymerType() const {
    validateBiopolymerType();
    return biopolymerType;
}

String BiopolymerClass::getBiopolymerTypeAsString() {
    if (   biopolymerType == BiopolymerType::RNA) return "RNA";
    else if (  biopolymerType == BiopolymerType::DNA) return "DNA";
    else if (  biopolymerType == BiopolymerType::Protein) return "Protein";
    else {
        ErrorManager::instance<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": Your have an unsupported biopolymerType : " << biopolymerType<<".  " <<endl;
        ErrorManager::instance.treatError();
    }
    validateBiopolymerType();
}

bool BiopolymerClass::setRenumberPdbResidues (bool tempRenumberPdbResidues){
    myRenumberPdbResidues = tempRenumberPdbResidues;
    std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Just set myRenumberPdbResidues to "<<getRenumberPdbResidues()<<" for chain "<<getChainID()<<std::endl;
    }


/**
 * \brief This polymorphism supports e.g. FirstResidue, LastResidue.
 *
 */

ResidueID BiopolymerClass::residueID(map<const String,double> myUserVariables,  const char* value){
        String tempString(value);
        if ((tempString.substr(0,1)).compare("@") ==0) { // if the String starts with '@' , this is a user-defined integer variable.  Note that insertion codes cannot be specified with this method.
            ResidueID myResidueID( myUserVariables, value);
            validateResidueID(myResidueID );
            return myResidueID;
        } else { // if the residue ID is supplied as a String literal, just validate and return it.  This String can contain insertion codes.
            return residueID(String(value));
        } 
};

    /**
     * \brief Make sure residue number is in proper range.
     * This polymorphism does NOT support e.g. FirstResidue, LastResidue.
     */
ResidueID BiopolymerClass::residueID(String inputString){
        
        ResidueID myResidueID(inputString/*, false*/); // set validate=false, because BiopolymerClass has its own validation.
        validateResidueID(myResidueID );
        return myResidueID;
};

void BiopolymerClass::validateResidueID(ResidueID myResidueID){
        //cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Validating requested residue ID "<<myResidueID.outString()<<endl;
        int myResidueIndex = getResidueIndex(myResidueID);
        validateResidueIndex(myResidueIndex);
}


void BiopolymerClass::validateResidueIndex(int myResidueIndex){
    if (myResidueIndex < 0) {
        ErrorManager::instance <<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<"Error! Residue index  "<<myResidueIndex<<" for chain "<<getChainID()<<" is less than zero"<<endl;
        ErrorManager::instance.treatError();
    }
    if (proteinCapping && (biopolymerType == BiopolymerType::Protein)) {
        if (myResidueIndex >  getChainLength()){
            ErrorManager::instance <<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<"Error! Residue index  "<<myResidueIndex<< " for chain "<<getChainID()<<  " is too big."<<endl;
            ErrorManager::instance.treatError();
        }
    } else if ((biopolymerType == BiopolymerType::RNA) || (biopolymerType == BiopolymerType::Protein) || (biopolymerType == BiopolymerType::DNA) ) {
        if (myResidueIndex >  (getChainLength() - 1)){
            ErrorManager::instance <<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<"Error! Residue index  "<<myResidueIndex<<" for chain "<<getChainID()<<" is too big."<<endl;
            ErrorManager::instance.treatError();
        }
    } else { 
        ErrorManager::instance <<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<"Unexplained Error!  "<< endl;
        ErrorManager::instance.treatError();
    }
}

////////////////////////////////////////////////////////////////////////////////////////
///calls myBiopolymer's hasAtom method to make sure the atom exists, dies otherwise. //
////////////////////////////////////////////////////////////////////////////////////////
int BiopolymerClass::validateAtomPathName(Compound::AtomPathName myAtomPathName){
    if (! myBiopolymer.hasAtom(myAtomPathName)){
        ErrorManager::instance <<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<"Error! Atom Path  "<<myAtomPathName<<" specifies an atom which does not exist"<<endl;
        ErrorManager::instance.treatError();
    } else return 0;
}

bool BiopolymerClass::hasAtom(ResidueID myResidueID, String myAtomName) {
         int myResidueIndex = getResidueIndex(myResidueID);
     Compound::AtomPathName myAtomPathName = atomPathString(myResidueID   , myAtomName);
     if ( myBiopolymer.hasAtom(myAtomPathName)){
         return bool(1);
     } else return bool(0);
}

////////////////////////////////////////////////////////////////////////////////////////
///generates an atom path e.g. 3/C1* .                                               //
////////////////////////////////////////////////////////////////////////////////////////

Compound::AtomPathName BiopolymerClass::atomPathString(ResidueID myResidueID, String myAtomName) {
         int myResidueIndex = getResidueIndex(myResidueID);
     Compound::AtomPathName myAtomPathName =  Compound::AtomPathName
         (String(
            intToString(myResidueIndex) +  // does this properly correct for proteinCapping? confirm empirically later.
            String("/") +
            myAtomName
            )
         );
   // validateAtomPathName(myAtomPathName);
   return myAtomPathName;
}

Compound::AtomIndex    BiopolymerClass::atomIndex(ResidueID myResidueID, String myAtomName) {
   Compound::AtomPathName myAtomPathString = atomPathString(myResidueID,  myAtomName);
   return myBiopolymer.getAtomIndex(myAtomPathString);
}


DuMM::AtomIndex    BiopolymerClass::getDuMMAtomIndex(ResidueID myResidueID, String myAtomName) {
   Compound::AtomPathName myAtomPathString = atomPathString(myResidueID,  myAtomName);
   Compound::AtomIndex myAtomIndex = myBiopolymer.getAtomIndex(myAtomPathString);
   return myBiopolymer.getDuMMAtomIndex(myAtomIndex); 
}


//////////////////////////////////////////////////////////////////////
/// Retrieves the location of the atom in its mobilized body frame ///
//////////////////////////////////////////////////////////////////////

Vec3 BiopolymerClass::getAtomLocationInMobilizedBodyFrame(ResidueID myResidueID, String myAtomName){
   // this will call atomPathString, which will validate the residue number and name.
   Compound::AtomIndex myAtomIndex = atomIndex(myResidueID,  myAtomName );
   return myBiopolymer.getAtomLocationInMobilizedBodyFrame(myAtomIndex);   
}

// mmbAtomInfo WITHOUT dumm, doesn't set mass, atomicNumber, mobilizedBody, or mobilizedBodyIndex.

MMBAtomInfo BiopolymerClass::mmbAtomInfo(ResidueID myResidueID, ResidueInfo::AtomIndex myResidueInfoAtomIndex,  SimbodyMatterSubsystem& matter ) {
        ResidueInfo myResidueInfo = myBiopolymer.updResidue(getResidueIndex(myResidueID ));    
        Compound::AtomIndex myAtomIndex = myResidueInfo.getAtomIndex(myResidueInfoAtomIndex );
        Compound::AtomName  myAtomName  = myResidueInfo.getAtomName(myResidueInfoAtomIndex );
        MMBAtomInfo myMMBAtomInfo;
        myMMBAtomInfo.compoundAtomIndex = myAtomIndex;
            myMMBAtomInfo.atomName = myAtomName; 
            myMMBAtomInfo.residueID = myResidueID;
            myMMBAtomInfo.chain = getChainID();
        return myMMBAtomInfo;
}
// mmbAtomInfo WITH dumm, adds mass, atomicNumber, mobilizedBody, and mobilizedBodyIndex.
MMBAtomInfo BiopolymerClass::mmbAtomInfo(ResidueID myResidueID, ResidueInfo::AtomIndex myResidueInfoAtomIndex,  SimbodyMatterSubsystem& matter, DuMMForceFieldSubsystem & dumm) {
        ResidueInfo myResidueInfo = myBiopolymer.updResidue(getResidueIndex(myResidueID ));    
        Compound::AtomIndex myAtomIndex = myResidueInfo.getAtomIndex(myResidueInfoAtomIndex );
        Compound::AtomName  myAtomName  = myResidueInfo.getAtomName(myResidueInfoAtomIndex );
        MMBAtomInfo myMMBAtomInfo = mmbAtomInfo(myResidueID, myResidueInfoAtomIndex, matter);
        DuMM::AtomIndex myDuMMAtomIndex = myBiopolymer.getDuMMAtomIndex(myAtomIndex);
        myMMBAtomInfo.mobilizedBody = updAtomMobilizedBody(matter, myResidueID ,myAtomName);
        myMMBAtomInfo.mobilizedBodyIndex =  myMMBAtomInfo.mobilizedBody.getMobilizedBodyIndex();
        //cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" initializing  myMMBAtomInfo.mobilizedBodyIndex = "<<  myMMBAtomInfo.mobilizedBodyIndex  <<endl;
        myMMBAtomInfo.atomName = myAtomName;
        myMMBAtomInfo.mass = dumm.getAtomMass(myDuMMAtomIndex); 
        myMMBAtomInfo.atomicNumber = dumm.getAtomElement(myDuMMAtomIndex); 
        myMMBAtomInfo.mobilizedBody = updAtomMobilizedBody(matter, myResidueID ,myAtomName);
        myMMBAtomInfo.mobilizedBodyIndex =  myMMBAtomInfo.mobilizedBody.getMobilizedBodyIndex();
        myMMBAtomInfo.partialCharge = dumm.getPartialCharge(myDuMMAtomIndex);
        return myMMBAtomInfo;
}
// Overrides default atom position, uses that from State instead. Not sure we'll ever need this. 
/*MMBAtomInfo BiopolymerClass::mmbAtomInfo(ResidueID myResidueID, ResidueInfo::AtomIndex myResidueInfoAtomIndex,  SimbodyMatterSubsystem& matter, DuMMForceFieldSubsystem & dumm, State & state) {
        ResidueInfo myResidueInfo = myBiopolymer.updResidue(getResidueIndex(myResidueID ));    
        Compound::AtomName  myAtomName  = myResidueInfo.getAtomName(myResidueInfoAtomIndex );
        MMBAtomInfo myMMBAtomInfo = mmbAtomInfo(myResidueID, myResidueInfoAtomIndex, matter, dumm);
            Vec3 myPositionVec3 = calcAtomLocationInGroundFrame(state, myResidueID, myAtomName);
            myMMBAtomInfo.position =openmmVecType(myPositionVec3[0],myPositionVec3[1], myPositionVec3[2]);
        return myMMBAtomInfo;
}*/

// This function loops through a provided vector<MMBAtomInfo> and for all phosphate atoms, sets atomicNumber to zero. This is intended to mask phosphates for density map fitting.
void overrideAtomInfoVectorProperties(BiopolymerClass & myBiopolymerClass, vector<MMBAtomInfo> & subjectAtomInfoVector, const vector<AtomicPropertyOverrideStruct> & myAtomicPropertyOverrideVector){
    //if (myBiopolymerClass.isRNA() || myBiopolymerClass.isDNA() ) {
        // Actually we will change atomicNumber to zero so it is inactive in density map fitting.
        for (int i = 0; i < subjectAtomInfoVector.size(); i++){
            //cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": About to check "<<subjectAtomInfoVector[i].atomName<<" in subjectAtomInfoVector["<<i<<"] "<<std::endl;
            for (int overrideVectorIndex = 0; overrideVectorIndex < myAtomicPropertyOverrideVector.size() ; overrideVectorIndex++){
                if (subjectAtomInfoVector[i].atomName == myAtomicPropertyOverrideVector[overrideVectorIndex].atomName){
                    if (myAtomicPropertyOverrideVector[overrideVectorIndex].property == "atomicNumber") {
                        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": For atom # "<<i<<", with name "<<subjectAtomInfoVector[i].atomName <<",  property atomicNumber is currently set to "<< subjectAtomInfoVector[i].atomicNumber <<std::endl; 
                        // NOte we are casting double as int:
                        subjectAtomInfoVector[i].atomicNumber = myAtomicPropertyOverrideVector[overrideVectorIndex].value;
                        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": For atom # "<<i<<", with name "<<subjectAtomInfoVector[i].atomName <<", just overrode property atomicNumber to "<< subjectAtomInfoVector[i].atomicNumber <<std::endl; 
                    } else {
	                ErrorManager::instance<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": You have tried to change the property : >" <<myAtomicPropertyOverrideVector[overrideVectorIndex].property << "< for atoms of name : >"<<  myAtomicPropertyOverrideVector[overrideVectorIndex].atomName <<"< .  This is not supported!"<<endl; 
                        ErrorManager::instance.treatError();
                    } // of if atomicNumber
                } // of if atomName match
            } // of for overrideVectorIndex
            /*
            if ((subjectAtomInfoVector[i].atomName == "P") ||
               (subjectAtomInfoVector[i].atomName == "OP1") ||
               (subjectAtomInfoVector[i].atomName == "OP2") ||
               (subjectAtomInfoVector[i].atomName == "O5*") ||
               (subjectAtomInfoVector[i].atomName == "O5'") ||
               (subjectAtomInfoVector[i].atomName == "O3*") ||
               (subjectAtomInfoVector[i].atomName == "O3'") ) {
                cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": About to change atomicNumber  from "<<subjectAtomInfoVector[i].atomicNumber <<" in subjectAtomInfoVector["<<i<<"] "<<std::endl;
                subjectAtomInfoVector[i].atomicNumber = 0; // Artificially set the atomic number to zero so it becomes inactive in density map fitting.
                cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<":  atomicNumber  is now  "<<subjectAtomInfoVector[i].atomicNumber <<" in subjectAtomInfoVector["<<i<<"] "<<std::endl;
                //cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": About to delete "<<subjectAtomInfoVector[i].atomName<<" in subjectAtomInfoVector["<<i<<"] "<<std::endl;
                //subjectAtomInfoVector.erase(subjectAtomInfoVector.begin() + i);
                //i--; // Now we will need to revisit the current i, since the vector has been shortened at this position.
            } // of if atomName
            */
        } // of for i
    //} // of if RNA/DNA
    //else { 
        // We are a protein, so includePhosphates should not affect us. Do nothing.
    //}

} // of overrideAtomInfoVectorProperties

#ifdef USE_OPENMM
// Without dumm, doesn't load certain properties..
void BiopolymerClass::initializeAtomInfoVector(SimbodyMatterSubsystem& matter,  const vector<AtomicPropertyOverrideStruct>  & myAtomicPropertyOverrideVector ) {
    if (atomInfoVector.size() > 0 ) {
	  ErrorManager::instance<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": initializeAtomInfoVector has already been called!"<<endl; 
	  ErrorManager::instance.treatError();
    }
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Pre-loading atomInfoVector for chain >"<<getChainID()<<"<"<<endl;
    PdbChain myPdbChain = PdbChain(myBiopolymer,myBiopolymer.getTopLevelTransform());

    for (ResidueID j = getFirstResidueID(); j <= getLastResidueID() ; incrementResidueID(j)) { 
        ResidueInfo myResidueInfo = myBiopolymer.updResidue(getResidueIndex(j));
        for (ResidueInfo::AtomIndex k (0) ;k < myResidueInfo.getNumAtoms() ; k++) {
            //cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" j,k "<<j.outString()<<", "<<k<<endl;
            MMBAtomInfo myAtomInfo = mmbAtomInfo(j,k,matter);
            myAtomInfo.setResidueIndex(getResidueIndex(j));
            myAtomInfo.setChain(getChainID());
            myAtomInfo.setResidueID(j);	
            //cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<endl;
           const PdbAtom& myPdbAtom = myPdbChain.getAtom(myAtomInfo.atomName, PdbResidueId(j.getResidueNumber(), j.getInsertionCode()));
            //cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<endl;
            Vec3 myPositionVec3 = myPdbAtom.getCoordinates();
            //cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<endl;
            myAtomInfo.position =openmmVecType (myPositionVec3[0],myPositionVec3[1], myPositionVec3[2]);
            //myAtomInfo.position =OpenMM::Vec3(myPositionVec3[0],myPositionVec3[1], myPositionVec3[2]);
            //cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<endl;
            //myAtomInfo.print();
            //cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<endl;
            atomInfoVector.push_back(myAtomInfo);   
            //cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<endl;
        } // of for k
        if (j == getLastResidueID() ) break;
    } // of for j
    // now if   maskPhosphates is true,  we set the corresponding atomic numbers to zero.
    //if ((maskPhosphates)){
    overrideAtomInfoVectorProperties(*this, atomInfoVector,myAtomicPropertyOverrideVector);
        //otherwise, do nothing. Phosphates on nucleic acids will get treated just like all other atoms for density map fitting purposes.
    //}



} // of initializeAtomInfoVector

void BiopolymerClass::initializeAtomInfoVector(SimbodyMatterSubsystem& matter, DuMMForceFieldSubsystem & dumm, const vector<AtomicPropertyOverrideStruct>  & myAtomicPropertyOverrideVector) {
    // have to allow rerun actually, since the dumm version can't be called much earlier.
    /*if (atomInfoVector.size() > 0 ) {
      ErrorManager::instance<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": initializeAtomInfoVector has already been called!"<<endl; 
      ErrorManager::instance.treatError();
    }*/
          atomInfoVector.clear();
          ignoreAtomPositionVector.clear();
          PdbChain myPdbChain = PdbChain(myBiopolymer,myBiopolymer.getTopLevelTransform());
          //cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<endl;
          for (ResidueID j = getFirstResidueID(); j <= getLastResidueID() ; incrementResidueID(j)) { 
            ResidueInfo myResidueInfo = myBiopolymer.updResidue(getResidueIndex(j));
                   for (ResidueInfo::AtomIndex k (0) ;k < myResidueInfo.getNumAtoms() ; k++) {
                        MMBAtomInfo myAtomInfo = mmbAtomInfo(j,k,matter,dumm);
                        myAtomInfo.setResidueIndex(getResidueIndex(j));
		        myAtomInfo.setChain(getChainID());
		        myAtomInfo.setResidueID(j);	
                        // Added this as a fix:
			const PdbAtom& myPdbAtom = myPdbChain.getAtom(myAtomInfo.atomName, PdbResidueId(j.getResidueNumber(), j.getInsertionCode()));
			Vec3 myPositionVec3 = myPdbAtom.getCoordinates();
			myAtomInfo.position =openmmVecType(myPositionVec3[0],myPositionVec3[1], myPositionVec3[2]);
                        ////

                        atomInfoVector.push_back(myAtomInfo);   
                   } // of for k
            if (j == getLastResidueID() ) break;
    } // of for j
    // now if   maskPhosphates is true,  we set the corresponding atomic numbers to zero.
    //if ((maskPhosphates)){
    overrideAtomInfoVectorProperties(*this,atomInfoVector, myAtomicPropertyOverrideVector);
        //otherwise, do nothing. Phosphates on nucleic acids will get treated just like all other atoms for density map fitting purposes.
    //}
} // of initializeAtomInfoVector
#endif

const vector<MMBAtomInfo>  BiopolymerClass::getAtomInfoVector(){
    validateAtomInfoVector();
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Inside getAtomInfoVector(). Chain "<<getChainID()<<" has an atomInfoVector of length : "<<atomInfoVector.size()<<endl;
    return atomInfoVector;
}


const vector<MMBAtomInfo>  BiopolymerClass::calcAtomInfoVector(ResidueStretch myResidueStretch, SimbodyMatterSubsystem& matter, DuMMForceFieldSubsystem & dumm, const bool includePhosphates ) {


    vector<MMBAtomInfo> returnAtomInfoVector;
    if ((myResidueStretch.getStartResidue() == getFirstResidueID()) && 
        (myResidueStretch.getEndResidue()   == getLastResidueID()  )) {
        validateAtomInfoVector(); //return atomInfoVector;
        returnAtomInfoVector = atomInfoVector;
    } // just return the precomputed atomInfoVector
    else {
          vector<MMBAtomInfo>::iterator startAtomInfoIterator;
          vector<MMBAtomInfo>::iterator endAtomInfoIterator;
          ResidueInfo myEndResidueInfo = myBiopolymer.updResidue(getResidueIndex(  myResidueStretch.getEndResidue() ));
          MMBAtomInfo   myStartAtomInfo =  mmbAtomInfo(myResidueStretch.getStartResidue(), ResidueInfo::AtomIndex(0), matter,dumm ) ;
          MMBAtomInfo   myEndAtomInfo   =  mmbAtomInfo(myResidueStretch.getEndResidue(), ResidueInfo::AtomIndex(myEndResidueInfo.getNumAtoms()-1), matter,dumm ) ;

          startAtomInfoIterator =   atomInfoVector.begin();
          ResidueID indexResidueID = getFirstResidueID();
          while ( indexResidueID < myResidueStretch.getStartResidue()) { 
              startAtomInfoIterator += myBiopolymer.updResidue(getResidueIndex(indexResidueID)).getNumAtoms();
              if (indexResidueID <  getLastResidueID() ) incrementResidueID(indexResidueID); else break; // make sure we don't increment past the last residue
       
          }
          endAtomInfoIterator =   startAtomInfoIterator ;
          ResidueID indexResidueID2 = myResidueStretch.getStartResidue();      
          while ( indexResidueID2 <= myResidueStretch.getEndResidue()) { 
              endAtomInfoIterator += myBiopolymer.updResidue( getResidueIndex(indexResidueID2) ).getNumAtoms();
              if (indexResidueID2 <  getLastResidueID() ) incrementResidueID(indexResidueID2); else break; // make sure we don't increment past the last residue
          }
          endAtomInfoIterator -= 1;
          returnAtomInfoVector = vector<MMBAtomInfo>  (startAtomInfoIterator, endAtomInfoIterator+1);
          //return vector<MMBAtomInfo>  (startAtomInfoIterator, endAtomInfoIterator+1);
    }
 
        /* 
    if (includePhosphates){
        // do nothing. Phosphates on nucleic acids will get treated just like all other atoms for density map fitting purposes.
    } else if (isRNA() || isDNA() ) {
        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": This section is obsolete!"<<std::endl;
        // Now we need to delete the phosphates from returnAtomInfoVector.
        for (int i = 0; i < returnAtomInfoVector.size(); i++){
            //cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": About to check "<<returnAtomInfoVector[i].atomName<<" in returnAtomInfoVector["<<i<<"] "<<std::endl;
            if ((returnAtomInfoVector[i].atomName == "P") ||
               (returnAtomInfoVector[i].atomName == "OP1") ||
               (returnAtomInfoVector[i].atomName == "OP2") ||
               (returnAtomInfoVector[i].atomName == "O5*") ||
               (returnAtomInfoVector[i].atomName == "O5'") ||
               (returnAtomInfoVector[i].atomName == "O3*") ||
               (returnAtomInfoVector[i].atomName == "O3'") ) {
                //cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": About to delete "<<returnAtomInfoVector[i].atomName<<" in returnAtomInfoVector["<<i<<"] "<<std::endl;
                returnAtomInfoVector.erase(returnAtomInfoVector.begin() + i);
                i--; // Now we will need to revisit the current i, since the vector has been shortened at this position.
            } // of if atomName
        } // of for i
    } // of if RNA/DNA
    else { 
        // We are a protein, so includePhosphates should not affect us. Do nothing.
    }
        */
    return returnAtomInfoVector;
    //cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": Unexplained error! "<<endl; exit (0);
}

void BiopolymerClass::addRingClosingBond( CovalentBondClass myCovalentBondClass) {
        addRingClosingBond( myCovalentBondClass.getResidueID1(), myCovalentBondClass.getAtomName1(), myCovalentBondClass.getBondCenterName1(), 
            myCovalentBondClass.getResidueID2(), myCovalentBondClass.getAtomName2(), myCovalentBondClass.getBondCenterName2(), myCovalentBondClass.getBondMobility());
}

void BiopolymerClass::addRingClosingBond( ResidueID residueID1, String atomName1, String bondCenterName1,  ResidueID residueID2, String atomName2,String bondCenterName2, SimTK::BondMobility::Mobility bondMobility ){
        const Compound::BondCenterPathName & centerName1 =  String(getResidueIndex(residueID1))+String('/')+String(atomName1)+String('/')+String(bondCenterName1);
        const Compound::BondCenterPathName & centerName2 = String(getResidueIndex(residueID2))+String('/')+String(atomName2)+String('/')+String(bondCenterName2);
        if (!(myBiopolymer.hasBondCenter(centerName1))){
            ErrorManager::instance<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Unable to find bond center "<<centerName1<<std::endl;        
            ErrorManager::instance.treatError();
    }
        if (!(myBiopolymer.hasBondCenter(centerName2))){
            ErrorManager::instance<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Unable to find bond center "<<centerName2<<std::endl;        
            ErrorManager::instance.treatError();
    }
    double bondLength = 111.111 ; // This doesn't matter, so I set to an absurd value
    double dihedralAngle = 0.0; // This doesn't matter either.
    myBiopolymer.addRingClosingBond( centerName1,    centerName2 , bondLength, dihedralAngle, bondMobility); 
}

MobilizedBody & BiopolymerClass::updAtomMobilizedBody(SimbodyMatterSubsystem & matter, ResidueID myResidueID    , String myAtomName){ 
    //Compound::AtomIndex myAtomIndex = atomIndex (myResidueID,myAtomName ); 
    MobilizedBodyIndex myAtomMobilizedBodyIndex = getAtomMobilizedBodyIndex(matter,myResidueID,myAtomName ); 
    return matter.updMobilizedBody(myAtomMobilizedBodyIndex);
}

MobilizedBodyIndex BiopolymerClass::getAtomMobilizedBodyIndex(SimbodyMatterSubsystem & matter, ResidueID myResidueID    , String myAtomName){ 
    Compound::AtomIndex myAtomIndex = atomIndex (myResidueID,myAtomName ); 
    MobilizedBodyIndex myAtomMobilizedBodyIndex = myBiopolymer.getAtomMobilizedBodyIndex(myAtomIndex); 
    return myAtomMobilizedBodyIndex;
    //return matter.updMobilizedBody(myAtomMobilizedBodyIndex);
}

Vec3 BiopolymerClass::calcDefaultAtomLocationInGroundFrame(ResidueID myResidueID, String atomName ){
    Compound::AtomIndex myAtomIndex = atomIndex(myResidueID,atomName);
    return myBiopolymer.calcDefaultAtomLocationInGroundFrame (atomPathString(myResidueID,atomName));    
}

Vec3 BiopolymerClass::calcAtomLocationInGroundFrame(const State & state,  ResidueID myResidueID, String atomName ){
    Compound::AtomIndex myAtomIndex = atomIndex(myResidueID,atomName);
    //cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": "<< myBiopolymer.calcAtomLocationInGroundFrame(state,myAtomIndex)<<endl;    
    return myBiopolymer.calcAtomLocationInGroundFrame(state,myAtomIndex);    
}


void BiopolymerClass::loadResidueIDVector() {
    ErrorManager::instance<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": You should not be doing this at this stage!  This is being done in BiopolymerClass::setResidueIDsAndInsertionCodesFromBiopolymer."<<endl; 
    ErrorManager::instance.treatError();

    for (int i = 0; i < getChainLength(); i++) {
        ResidueID myResidueID;    
        ResidueInfo myResidueInfo ( myBiopolymer.getResidue(ResidueInfo::Index( i)));
        myResidueID.setResidueNumber  (myResidueInfo.getPdbResidueNumber());
        myResidueID.setInsertionCode  (myResidueInfo.getPdbInsertionCode());
        residueIDVector.push_back(myResidueID);
    }
}

void BiopolymerClass::loadResidueIDVectorAscending(ResidueID firstResidueID ){
    if (residueIDVector.size() > 0) {
        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" "<<std::endl;
        for (int i = 0; i < residueIDVector.size() ; i++){
            cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": "<<residueIDVector[i].outString()<<std::endl;
        }
        ErrorManager::instance<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Why does residueIDVector have something in it already?"<<endl; 
        ErrorManager::instance.treatError();
    }
    ResidueID myResidueID = firstResidueID;
    for (int i = 0; i < getChainLength(); i++) {
        residueIDVector.push_back(myResidueID);
        myResidueID.setResidueNumber(myResidueID.getResidueNumber()+1);
        myResidueID.setInsertionCode(' ');
    }
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Just loaded residueIDVector with defaults starting with "<<firstResidueID.outString()<<endl;
    // for (int i = 0; i < getChainLength(); i++) {
    //     cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<residueIDVector[i].outString()<<endl;
    // }
}

const ResidueInfo::Index BiopolymerClass::getResidueIndex(ResidueID residueID){
    int residueIndex;
    if (residueIDVector.size() >0){
    vector<ResidueID>::iterator residueIDVectorIterator ;
    //ResidueID tempRes =(* (residueIDVector.begin()) );          
    residueIDVectorIterator = find(residueIDVector.begin(), residueIDVector.end(), residueID);
    int residueIDVectorPosition = residueIDVectorIterator-residueIDVector.begin();
        
        residueIndex = ResidueInfo::Index(residueIDVectorPosition);
        //std::cout <<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Your residue ID: "<<residueID.outString() << " has a corresponding residue index : "<<residueIndex<<std::endl;  
        if ((residueIndex < 0 ) || (residueIndex >= getChainLength())) {
            ErrorManager::instance<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Encountered a problem with residue ID "<<residueID.outString()<<" of chain "<<getChainID() <<" . This returned and index of "<<residueIndex<<". The residue ID should lie in the closed interval "<<getFirstResidueID().outString()<<" , "<<getLastResidueID().outString()<<". If you are performing an arithmetic (+/-) operation on a residue number, the leftmost term correspond to an existing residue number, while the rest of the terms are increments in sequence to be added or subtracted from that residue number. Or, you maybe you issued loadSequencesFromPdb and there is no residue numbered "<<residueID.outString() <<endl;
            ErrorManager::instance <<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" The computed index : "<<residueIndex<< " is unreasonable and would be expected to be in the range : 0 to "<<(getChainLength()-1)<<endl;    
            //ErrorManager::instance <<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Tried to extract residue index corresponding to residue ID "<<residueID.outString()<<" . The residue indices should lie in the closed interval "<<getFirstResidueID().outString()<<" , "<<getLastResidueID().outString()<<" . The resulting index : "<<residueIndex<< " is out of the corresponding range : 0 to "<<(getChainLength()-1)<<endl;
            ErrorManager::instance.treatError();
        }
        return ResidueInfo::Index(residueIndex);
    } else { // this "else" block is really inefficient .. for when residueIDVector is empty.  Try to avoid going into this!
        // ResidueInfo myResidueInfo ( myBiopolymer.getResidue(ResidueInfo::Index( 0)));
        for (int i = 0; i < getChainLength(); i++) {
            ResidueInfo myResidueInfo = ResidueInfo ( myBiopolymer.getResidue(ResidueInfo::Index( i)));
            if (
                (myResidueInfo.getPdbResidueNumber() ==   residueID.getResidueNumber()) &&
                (myResidueInfo.getPdbInsertionCode() ==   residueID.getInsertionCode()) 
                )
                {
                    residueIndex = i;
                    validateResidueIndex(residueIndex);
                    return ResidueInfo::Index(residueIndex);
                    //continue;
                }
        }     
        ErrorManager::instance<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Encountered a problem with residue ID "<<residueID.outString()<<" . The residue ID should lie in the closed interval "<<getFirstResidueID().outString()<<" , "<<getLastResidueID().outString()<<endl;
        ErrorManager::instance <<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" The computed index : "<<residueIndex<< " is unreasonable and would be expected to be in the range : 0 to "<<(getChainLength()-1)<<endl;
        ErrorManager::instance.treatError();
    }
        
    ErrorManager::instance <<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<"Unexplained error"<<endl;  
    ErrorManager::instance.treatError();
}

/////////////////////////////////////////////////////////////////////////////
/// calls getResidueIndex, which validates residue number. returns name.  ///
/////////////////////////////////////////////////////////////////////////////
String  BiopolymerClass::getPdbResidueName( ResidueID residueID){
    return myBiopolymer.getResidue(getResidueIndex(residueID)).getPdbResidueName(); 
}

String BiopolymerClass::getRepresentativeAtomName(){
    if (biopolymerType     ==  BiopolymerType::RNA)     {
        return  "C4*";
    }
    else if (biopolymerType     ==  BiopolymerType::DNA)     {
        return  "C4*";
    }
    else if (biopolymerType == BiopolymerType::Protein) {
        return  "CA";
    }
    else {
        ErrorManager::instance << __FUNCTION__ << ": biopolymerType " << biopolymerType << " unknown" << endl;
        ErrorManager::instance.treatError();
    }
}


double BiopolymerClass::getRepresentativeAtomMassThreshold(){
    if (biopolymerType     ==  BiopolymerType::RNA)     {
        return  30.   ;
    }
    else if (biopolymerType     ==  BiopolymerType::DNA)     {
        return  30.   ;
    }
    else if (biopolymerType == BiopolymerType::Protein) {
        return  18.  ;
    }
    else {
        ErrorManager::instance << __FUNCTION__ << ": biopolymerType " << biopolymerType << " unknown" << endl;
        ErrorManager::instance.treatError();
    }
}

// Set active = false to effectively turn off sterics (e.g. when using the "scrubber").  Note that the compute cost of evaluating the forces may not change. 
void BiopolymerClass::setContactParameters(GeneralContactSubsystem & contacts,  HuntCrossleyForce & hc, double excludedVolumeStiffness, bool active ){
    double huntCrossleyDissipation  =.0;
    if (!(active)) excludedVolumeStiffness = .0; // set force constants to zero
    hc.setBodyParameters(SimTK::ContactSurfaceIndex(contacts.getNumBodies(hc.getContactSetIndex()  )-1),excludedVolumeStiffness ,huntCrossleyDissipation, 0., 0., 0.);
};


void BiopolymerClass::addGeneralSterics(GeneralContactSubsystem & contacts, ContactSetIndex contactSet, HuntCrossleyForce & hc,SimbodyMatterSubsystem & matter,double excludedVolumeRadius,double excludedVolumeStiffness,  ResidueID startResidue, ResidueID endResidue, bool endCapsOn, bool addHydrogens) {

// start adding spheres
    double huntCrossleyDissipation  =.0;
    
    for (int q=(getResidueIndex(startResidue)); q<=(getResidueIndex(endResidue));q++) 
    {
        ResidueInfo myResidueInfo = myBiopolymer.updResidue((ResidueInfo::Index(q)));
        String myPdbResidueName1 =  getPdbResidueName(getResidueID(q));
        for (ResidueInfo::AtomIndex r(0) ; r<ResidueInfo::AtomIndex(  myResidueInfo.getNumAtoms()); r++)
        { //loop over four possible interacting pairs of atoms.
            { //if the atom name field is not blank, do this
                stringstream ss3;
                ss3<<q<<"/"<<myResidueInfo.getAtomName(r);
                // excluded volume radius WAS 1.25.  but now using much smaller radii.                        // following:
                //Whitford PC, Noel JK, Gosavi S, Schug A, Sanbonmatsu KY & Onuchic JN, "An All-atom Structure-Based Potential for Proteins: Bridging Minimal Models with All-atom Empirical Forcefields" PROTEINS (2008) DOI: 10.1002/prot.22253.


                SimTK_ERRCHK_ALWAYS(
                    (myBiopolymer.hasAtom(ss3.str()) ),
                    __FILE__,
                    "Failed to attach steric sphere.  Could not find specified atom"
                );//: %s .",String(ss3.str()));

                if (myBiopolymer.hasAtom(ss3.str()))
                    if (addHydrogens || (((myBiopolymer.getAtomElement(myResidueInfo.getAtomIndex( r  ))).getSymbol()).compare("H") != 0))
                    {                           
                        contacts.addBody(contactSet,
                                         (matter.updMobilizedBody(myBiopolymer.getAtomMobilizedBodyIndex(Compound::AtomIndex(myBiopolymer.getAtomIndex(ss3.str()))))),
                                          ContactGeometry::Sphere(excludedVolumeRadius/1), //used to convert from angstrom to nanometers.
                                          (myBiopolymer.getAtomLocationInMobilizedBodyFrame(myBiopolymer.getAtomIndex(ss3.str())))                                        
                                        );
                        setContactParameters(contacts, hc, excludedVolumeStiffness, true);
                        // hc.setBodyParameters(SimTK::ContactSurfaceIndex(contacts.getNumBodies(contactSet)-1),excludedVolumeStiffness ,huntCrossleyDissipation, 0., 0., 0.);
                    } //if
            }//
        }// for
    }// for   
}

void BiopolymerClass::addCustomSterics(GeneralContactSubsystem & contacts, ContactSetIndex contactSet, HuntCrossleyForce & hc,SimbodyMatterSubsystem & matter,LeontisWesthofClass myLeontisWesthofClass,String leontisWesthofInteractionType, ResidueID startResidue,ResidueID endResidue, bool endCapsOn) {
// start adding spheres
     double huntCrossleyDissipation  =.0;
     //for (int q=(getResidueIndex(startResidue)); q<=(getResidueIndex(endResidue));q++) {
     for (ResidueID q=(startResidue); q<=(endResidue);incrementResidueID(q)) {
            ResidueInfo myResidueInfo = myBiopolymer.updResidue(ResidueInfo::Index(getResidueIndex(q)));
            String myPdbResidueName1 = (myResidueInfo).getPdbResidueName();
            LeontisWesthofBondRow myLeontisWesthofBondRow = myLeontisWesthofClass.getLeontisWesthofBondRow(
                ResidueID((myResidueInfo).getPdbResidueNumber(), (myResidueInfo).getPdbInsertionCode() ),
                ResidueID((myResidueInfo).getPdbResidueNumber(), (myResidueInfo).getPdbInsertionCode()  ) ,
                myPdbResidueName1,
                leontisWesthofInteractionType,"",leontisWesthofInteractionType,"X","contact");
            for (int r =0; r<4; r++) { //loop over four possible interacting pairs of atoms.
                if ((myLeontisWesthofBondRow.residue1Atom[r]).compare("") != 0) { //if the atom name field is not blank, do this

                        stringstream ss4;
                        ss4<<q.outString()<<"/"<<myLeontisWesthofBondRow.residue1Atom[r];                  
                        if ((q == getFirstResidueID()  ) && 
                            (
                            ((myLeontisWesthofBondRow.residue1Atom[r]).compare("P") == 0) ||
                            ((myLeontisWesthofBondRow.residue1Atom[r]).compare("OP1") == 0) ||
                            ((myLeontisWesthofBondRow.residue1Atom[r]).compare("OP2") == 0) 
                            )
                           ) {
                            cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Warning:  You are attempting to place a contact sphere on an atom that doesn't exist, in this case an omitted 5' Phosphate atom "<<myLeontisWesthofBondRow.residue1Atom[r]<<" at residue ID "<<q.outString()<<"; this is pretty harmless.                 "<<endl;
                        }     
                        else {
                           String ss3 = atomPathString(q,myLeontisWesthofBondRow.residue1Atom[r]); // can't have this statement earlier because of the P, OP1, OP2 problem.
                           contacts.addBody(contactSet,
                                        (matter.updMobilizedBody(myBiopolymer.getAtomMobilizedBodyIndex(Compound::AtomIndex(myBiopolymer.getAtomIndex(ss3))))),
                                        ContactGeometry::Sphere(myLeontisWesthofBondRow.bondLength[r]/10),                                        (myBiopolymer.getAtomLocationInMobilizedBodyFrame(myBiopolymer.getAtomIndex(ss3)))
                                        );
                           hc.setBodyParameters(SimTK::ContactSurfaceIndex(contacts.getNumBodies(contactSet)-1),myLeontisWesthofBondRow.springConstant
[r] ,huntCrossleyDissipation, 0., 0., 0.);
                       } 

                }
            }



         if (q==(endResidue)) break;
     }   
}

void BiopolymerClass::setProteinBondMobility ( BondMobility::Mobility  mobility, ResidueID startResidue, ResidueID endResidue){
    if (!(biopolymerType == BiopolymerType::Protein)) {
        ErrorManager::instance <<"setProteinBondMobility can only be used on Protein, type "<< biopolymerType<<" setProteinBondMobility.  You have tried to apply it to "<<biopolymerType<<endl;
        ErrorManager::instance.treatError();
    }

    validateResidueID(startResidue);
    validateResidueID(  endResidue);
        ResidueID i = startResidue;
        while ( i <= endResidue){
            myBiopolymer.setResidueBondMobility(getResidueIndex(i), mobility);
            if (i>(startResidue)) {
                ResidueID tempResidueID = i;
                myBiopolymer.setBondMobility(mobility ,atomPathString(decrementResidueID(tempResidueID),String("C")) ,atomPathString(i ,String("N"))  );
            }
            if (i == endResidue) {
                break;
            }
            incrementResidueID(i);
        }
        }

Biopolymer & BiopolymerClass::updBiopolymer() {
    return  myBiopolymer;
}


void BiopolymerClass::includeNonBondAtom(ResidueID residueID, String atomName, State & state, DuMMForceFieldSubsystem & dumm) {
        // cout<< __FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" About to add residue : "<<residueID.outString()<<" atom : "<<atomName<<endl;
    DuMM::AtomIndex myDuMMAtomIndex = getDuMMAtomIndex(residueID, atomName);
    dumm.includeNonbondAtom(myDuMMAtomIndex);
}

ResidueInfo BiopolymerClass::updResidueInfo (ResidueID residueID) {
        validateResidueID(residueID);
    return myBiopolymer.updResidue((getResidueIndex(residueID)));
}

void BiopolymerClass::includeAllNonBondAtomsInResidue(ResidueID residueID, State & state, DuMMForceFieldSubsystem & dumm) {
    ResidueInfo myResidueInfo = updResidueInfo( residueID )  ;
    for (SimTK::ResidueInfo::AtomIndex i  = SimTK::ResidueInfo::AtomIndex(0) ; i <myResidueInfo.getNumAtoms(); i++) {
        includeNonBondAtom ( residueID, myResidueInfo.getAtomName(i) , state,dumm)  ;
    }

}


void BiopolymerClass::constrainRigidSegmentsToGround(CompoundSystem & system,  SimbodyMatterSubsystem & matter,State & state, ConstraintToGroundContainer & myConstraintToGroundContainer , bool toGround = true, ResidueID baseResidue = ResidueID() ) {
    MobilizedBody myOldBody = updAtomMobilizedBody(matter,getFirstResidueID(),getRepresentativeAtomName());
    MobilizedBody myBody = myOldBody;
    ResidueID myResidueID = getFirstResidueID() ;
    ResidueID myOldResidueID = myResidueID;
    String baseChain = getChainID();
    
    
    while (myResidueID <= getLastResidueID()) {
            myBody = updAtomMobilizedBody(matter,myResidueID,getRepresentativeAtomName());
            ConstraintClass  myConstraintClass; 
            if (toGround) {
                // do nothing WRT setting chain2, residueID2, atomName2
            } else if (! (toGround)){
                // then we should specify what chain, residue, and atom name to constrain to:
                if (!(hasAtom(baseResidue, getRepresentativeAtomName()))) {
                    ErrorManager::instance << __FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Could not find residue "<<baseResidue.outString()<<", or maybe it has no atom named "<<baseResidue.outString()<<endl; ErrorManager::instance.treatError();
                }
                myConstraintClass.setChain2(baseChain);
                myConstraintClass.setResidueID2(baseResidue);
                myConstraintClass.setAtomName2(getRepresentativeAtomName());
                myConstraintClass.setConstraintType(WeldToAtom);//setToGround(false);
                
                // if the rigid segment to be constrained contains chain2, residueID2, atomName2, then we should increment residue numbers until that's over:               
                MobilizedBody body2 = updAtomMobilizedBody(matter, baseResidue, getRepresentativeAtomName());
                while (myBody.getMobilizedBodyIndex() == body2.getMobilizedBodyIndex() &&
                       (myResidueID < getLastResidueID()) )  {
                    incrementResidueID(myResidueID);
                    myBody = updAtomMobilizedBody(matter,myResidueID,getRepresentativeAtomName());
                    if (myResidueID == getLastResidueID()) {goto stop;}
                }
            }
            
            myOldBody      = myBody;
            myOldResidueID = myResidueID;
            MassProperties myBodyMassProperties = myBody.getBodyMassProperties(state);
                    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<"Checking out chain, residueID, atomName: "<<getChainID()<<", "<< myResidueID.outString() <<" , "<<getRepresentativeAtomName()<<endl;
                    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<"got mass: "<< myBody.getBodyMassProperties(state).getMass()<<endl;
                    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<"got MobilizedBodyIndex: "<< myBody.getMobilizedBodyIndex()               <<endl;
                    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<"got getFirstResidueMobilizerType() : "<< getFirstResidueMobilizerType()               <<endl;
                    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" ((myResidueID == getFirstResidueID()) evaluates to : "<< (myResidueID == getFirstResidueID())<<endl;
                    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" (getFirstResidueMobilizerType().compare(Weld)==0  ) evaluates to : "<<  (getFirstResidueMobilizerType().compare("Weld")==0  )<<endl;
            if (myBodyMassProperties.getMass() > getRepresentativeAtomMassThreshold()){ 
      
                while ((myBody.getMobilizedBodyIndex() == myOldBody.getMobilizedBodyIndex()) &&
                       (myResidueID < getLastResidueID()))
                {// figure out where rigid segment ends
                    incrementResidueID(myResidueID);
                    myBody = updAtomMobilizedBody(matter,myResidueID,getRepresentativeAtomName());
                }  
                myBody = updAtomMobilizedBody(matter,myResidueID,getRepresentativeAtomName());
                if ((myResidueID == getLastResidueID()) &&
                    (myOldResidueID ==getLastResidueID() ) &&
                    ( myBody.getBodyMassProperties(state).getMass() >getRepresentativeAtomMassThreshold()  )
                    ) 
                { 
                    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<"About to add constraint to chain, residueID, atomName: "<<getChainID()<<", "<< getLastResidueID().outString() <<" , "<<getRepresentativeAtomName()<<endl;
                    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<"got mass: "<< myBody.getBodyMassProperties(state).getMass()<<endl;
                    if (!( (getFirstResidueMobilizerType().compare("Weld")==0 ) && (getLastResidueID() == getFirstResidueID()))) 
                    {
                        myConstraintClass.setChain1(getChainID());
                        myConstraintClass.setResidueID1(getLastResidueID());
                        myConstraintClass.setAtomName1(getRepresentativeAtomName());
                        myConstraintToGroundContainer.addConstraintClassToVector(myConstraintClass); // getChainID(),getLastResidueID() ,getRepresentativeAtomName() );
                    }
                    break; // Our work is done.  Break out of while loop; if we don't we'll be in this loop infinitely.

                } else if ((myResidueID == getLastResidueID()) &&
                    (myOldResidueID < getLastResidueID() ) &&
                    (!( myBody.getBodyMassProperties(state).getMass() >getRepresentativeAtomMassThreshold()  ))
                    ) 
                { 
                    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<"About to add constraint to chain, residueID, atomName: "<<getChainID()<<", "<< getLastResidueID().outString() <<" , "<<getRepresentativeAtomName()<<endl;
                    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<"got mass: "<< myBody.getBodyMassProperties(state).getMass()<<endl;
                    if (!( (getFirstResidueMobilizerType().compare("Weld")==0 ) && (myOldResidueID == getFirstResidueID()))) {

                        myConstraintClass.setChain1(getChainID());
                        myConstraintClass.setResidueID1(myOldResidueID);
                        myConstraintClass.setAtomName1(getRepresentativeAtomName());
                        myConstraintToGroundContainer.addConstraintClassToVector(myConstraintClass );//getChainID(),myOldResidueID,getRepresentativeAtomName() );
                    }
                    break; // Our work is done.  Break out of while loop; if we don't we'll be in this loop infinitely.

                } else if ((myResidueID == getLastResidueID()) &&
                    (myOldResidueID == getLastResidueID() ) &&
                    (!( myBody.getBodyMassProperties(state).getMass() >getRepresentativeAtomMassThreshold()  ))
                    ) 
                { 
                    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<"About to add constraint to chain, residueID, atomName: "<<getChainID()<<", "<< getLastResidueID().outString() <<" , "<<getRepresentativeAtomName()<<endl;
                    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<"got mass: "<< myBody.getBodyMassProperties(state).getMass()<<endl;
                    if (!( (getFirstResidueMobilizerType().compare("Weld")==0 ) && (myOldResidueID == getFirstResidueID()))){
                        myConstraintClass.setChain1(getChainID());
                        myConstraintClass.setResidueID1(myOldResidueID);
                        myConstraintClass.setAtomName1(getRepresentativeAtomName());
                        myConstraintToGroundContainer.addConstraintClassToVector(myConstraintClass );//getChainID(),myOldResidueID,getRepresentativeAtomName() );
                    }
                    break; // Our work is done.  Break out of while loop; if we don't we'll be in this loop infinitely.

                } else if ((myResidueID == getLastResidueID()) &&
                    (myOldResidueID <  getLastResidueID() ) &&
                    (( myBody.getBodyMassProperties(state).getMass() > getRepresentativeAtomMassThreshold() ))
                    ) 
                { 
                    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<"About to add constraint to chain, residueID, atomName: "<<getChainID()<<", "<< getLastResidueID().outString() <<" , "<<getRepresentativeAtomName()<<endl;
                    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<"got mass: "<< myBody.getBodyMassProperties(state).getMass()<<endl;
                    if (!( (getFirstResidueMobilizerType().compare("Weld")==0 ) && (myOldResidueID == getFirstResidueID()))){
                        myConstraintClass.setChain1(getChainID());
                        myConstraintClass.setResidueID1(myOldResidueID);
                        myConstraintClass.setAtomName1(getRepresentativeAtomName());
                        myConstraintClass.setConstraintType(WeldToGround);
                        myConstraintToGroundContainer.addConstraintClassToVector(myConstraintClass);///getChainID(),myOldResidueID,getRepresentativeAtomName() );
                    }
                    break; // Our work is done.  Break out of while loop; if we don't we'll be in this loop infinitely.
                /*} else if ((myResidueID == getFirstResidueID()) && (getFirstResidueMobilizerType().compare("Weld")==0  )) {
                    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<endl;
                    break; // No point in adding a constraint to the first residue, if the root mobilizer is Weld!  That means it is already welded to ground!
                */
                } else if (myResidueID <  getLastResidueID()) 
                {
                    if ((myOldResidueID == getFirstResidueID()) && (getFirstResidueMobilizerType().compare("Weld")==0  )) {
                        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<endl;
                        continue; // No point in adding a constraint to the first residue, if the root mobilizer is Weld!  That means it is already welded to ground!
                    } 
                    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<"myResidueID ="<<myResidueID.outString()<<endl; 
                    cout<<" myBody.getBodyMassProperties(state).getMass() > getRepresentativeAtomMassThreshold() "<<endl;
                    cout<< myBody.getBodyMassProperties(state).getMass() <<", "<<getRepresentativeAtomMassThreshold() <<endl;
                    if ( myOldBody.getBodyMassProperties(state).getMass() >getRepresentativeAtomMassThreshold()  ) 
                    {
                        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<"About to add constraint to chain, residueID, atomName: "<<getChainID()<<", "<< myOldResidueID.outString() <<" , "<<getRepresentativeAtomName()<<endl;
                        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<"got mass: "<< myOldBody.getBodyMassProperties(state).getMass()<<endl;
                        if (!( (getFirstResidueMobilizerType().compare("Weld")==0 ) && (myOldResidueID == getFirstResidueID())))
                        myConstraintClass.setChain1(getChainID());
                        myConstraintClass.setResidueID1(myOldResidueID);
                        myConstraintClass.setAtomName1(getRepresentativeAtomName());
                        myConstraintToGroundContainer.addConstraintClassToVector(myConstraintClass);
                    }else {
                        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Skipped this constraint -- it's already immobile since getFirstResidueMobilizerType() = "<<getFirstResidueMobilizerType()<<endl;
                    }
                } else {
                    ErrorManager::instance <<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Unexplained error! "<<endl;
                    ErrorManager::instance.treatError();
                }
            } else if (myResidueID <  getLastResidueID()) 
                incrementResidueID(myResidueID); // check next residue's mass
            else if (myResidueID == getLastResidueID()) break;
    } // of while
    stop: cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Ending constrainRigidSegmentsToGround."<<endl;
    //myConstraintClass.validate
}


// This takes the vector<AllResiduesWithin> , and adds residues to it, if those residues are flexible/mobile as deterimined by checking their atom-associated masses.

void BiopolymerClass::physicsZone(vector<AllResiduesWithin> & myIncludeAllResiduesWithinVector , double radius, SimbodyMatterSubsystem & matter,State & state) {
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Creating physics zone "<<radius<<" in size about flexible atoms in chain "<<getChainID()<<endl;
    AllResiduesWithin myAllResiduesWithin;
    myAllResiduesWithin.setResidue ( ResidueID(-11111,' '));
    for (int i = 0; i < atomInfoVector.size() ; i++) {
            if (atomInfoVector[i].residueID  == myAllResiduesWithin.getResidue()) {continue;} // Each residue needs be added only once, no matter how many flexible atoms it has.
            MobilizedBody myBody = updAtomMobilizedBody(matter,atomInfoVector[i].residueID,atomInfoVector[i].atomName);
            MassProperties myBodyMassProperties = myBody.getBodyMassProperties(state);
            if (myBodyMassProperties.getMass() < 40){
                myAllResiduesWithin.setRadius   ( radius);
                myAllResiduesWithin.setChain   ( getChainID());
                myAllResiduesWithin.setResidue   ( atomInfoVector[i].residueID);
                myIncludeAllResiduesWithinVector.push_back(myAllResiduesWithin); 
               
                cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Implementng physicsRadius. Include around :"<< endl;
                myAllResiduesWithin.print();
            }
        
    } 
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" All AllResiduesWithin records have been printed for chain "<<getChainID()<<" . If none were printed, then MMB thinks this chain has no flexible residues." <<endl;
}

void BiopolymerClass::multiplySmallGroupInertia(ResidueID residueID, String atomName, double multiplier, CompoundSystem & system,  SimbodyMatterSubsystem & matter,State & state) {
            MobilizedBody myBody = updAtomMobilizedBody(matter,residueID,atomName);
            MassProperties myBodyMassProperties = myBody.getBodyMassProperties(state);
            if (myBodyMassProperties.getMass() < 40){
                myBody.setDefaultMassProperties (MassProperties(myBodyMassProperties.getMass(), myBodyMassProperties.getMassCenter(),  myBodyMassProperties.getInertia ()* (double)multiplier));
        state = system.realizeTopology();
                system.realize(state,Stage::Position);
            }
}

void BiopolymerClass::multiplySmallGroupInertia(double multiplier, CompoundSystem & system, SimbodyMatterSubsystem & matter,State & state) {

    if (multiplier != 1.0) // if unity, no need to waste time and risk rounding error.
        for (ResidueID i = getFirstResidueID()    ; i < getLastResidueID(); incrementResidueID(i)) {
            //cout<<i<<","<< getPdbResidueName(i)<<endl;
            if ( getPdbResidueName(i).compare("THR") == 0) {
                multiplySmallGroupInertia(i,"CG2",multiplier,system,matter,state); 
                multiplySmallGroupInertia(i,"OG1",multiplier,system,matter,state); 
            }

            else if ( getPdbResidueName(i).compare("ILE") == 0) {
                multiplySmallGroupInertia(i,"CG2",multiplier,system,matter,state); 
                multiplySmallGroupInertia(i,"CD" ,multiplier,system,matter,state); 
            }

            else if ( getPdbResidueName(i).compare("LEU") == 0) {
                multiplySmallGroupInertia(i,"CD1",multiplier,system,matter,state); 
                multiplySmallGroupInertia(i,"CD2" ,multiplier,system,matter,state); 
            }

            else if ( getPdbResidueName(i).compare("SER") == 0) {
                multiplySmallGroupInertia(i,"OG" ,multiplier,system,matter,state); 
            } 


            else if ( getPdbResidueName(i).compare("ASN") == 0) {
                multiplySmallGroupInertia(i,"ND2" ,multiplier,system,matter,state); 
            }

            else if ( getPdbResidueName(i).compare("GLN") == 0) {
                multiplySmallGroupInertia(i,"NE2" ,multiplier,system,matter,state); 
            }

            else if ( getPdbResidueName(i).compare("LYS") == 0) {
                multiplySmallGroupInertia(i,"NZ" ,multiplier,system,matter,state); 
            }

            else if ( getPdbResidueName(i).compare("GLU") == 0) {
                multiplySmallGroupInertia(i,"CD" ,multiplier,system,matter,state); 
            }

            else if ( getPdbResidueName(i).compare("MET") == 0) {
                multiplySmallGroupInertia(i,"CE" ,multiplier,system,matter,state); 
            }
            else if ( getPdbResidueName(i).compare("ALA") == 0) {
                multiplySmallGroupInertia(i,"CB" ,multiplier,system,matter,state); 
            }
            else if ( getPdbResidueName(i).compare("ARG") == 0) {
                multiplySmallGroupInertia(i,"NH1" ,multiplier,system,matter,state); 
                multiplySmallGroupInertia(i,"NH2" ,multiplier,system,matter,state); 
            }

            else if ( getPdbResidueName(i).compare("VAL") == 0) {
                multiplySmallGroupInertia(i,"CG1" ,multiplier,system,matter,state); 
                multiplySmallGroupInertia(i,"CG2" ,multiplier,system,matter,state); 
            }
            if (i ==getLastResidueID()) break;
    }
}


void BiopolymerClass::setResidueIDsAndInsertionCodesFromBiopolymer(const Biopolymer & inputBiopolymer, bool endCaps = 0      ) {
    // if (residueIDVector.size() > 0) 
    // {
    //     ErrorManager::instance<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Why does residueIDVector have something in it already?"<<endl; 
    //     ErrorManager::instance.treatError();
    // }
    residueIDVector.clear();
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Setting residue numbers and insertion codes from biopolymer in input structure file, for chain "<<getChainID()<<endl;
    for (int inputResidueIndex = ( 0 + endCaps ); inputResidueIndex < (inputBiopolymer.getNumResidues() - endCaps) ; inputResidueIndex ++) {
        int myResidueIndex = inputResidueIndex - endCaps + proteinCapping;
        const ResidueInfo inputResidueInfo = inputBiopolymer.getResidue(ResidueInfo::Index( inputResidueIndex ));
        const int  inputResidueNumber  = (inputResidueInfo).getPdbResidueNumber();
        const char inputInsertionCode  = (inputResidueInfo).getPdbInsertionCode();
        //const char inputOneLetterCode  = (inputResidueInfo).getOneLetterCode();
        //ResidueInfo myResidueInfo = myBiopolymer.updResidue(ResidueInfo::Index( myResidueIndex ));
        myBiopolymer.updResidue(ResidueInfo::Index( myResidueIndex )).setPdbResidueNumber(inputResidueNumber);
        myBiopolymer.updResidue(ResidueInfo::Index( myResidueIndex )).setPdbInsertionCode(inputInsertionCode);
        ResidueInfo myResidueInfo = myBiopolymer.updResidue(ResidueInfo::Index( myResidueIndex ));
        // cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Residue index, number, insertion code, and residue type: "<<myResidueIndex<<"," <<myResidueInfo.getPdbResidueNumber()   <<","<<myResidueInfo.getPdbInsertionCode()   <<", "<<myResidueInfo.getOneLetterCode()   <<endl;
        
        ////// Tbis was previously done in loadResidueIDVector().  However most of the time is spent in the myBiopolymer.updResidue step.  So for efficiency I am now doing it here: //// 
        ResidueID myResidueID;    
        myResidueID.setResidueNumber  (myResidueInfo.getPdbResidueNumber());
        myResidueID.setInsertionCode  (myResidueInfo.getPdbInsertionCode());
        residueIDVector.push_back(myResidueID);
        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" myResidueIndex "<<myResidueIndex<<" inputResidueNumber "<<inputResidueNumber<<" inputInsertionCode "<< inputInsertionCode << " myResidueID " <<myResidueID.outString()<<endl;
        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" myBiopolymer.updResidue(ResidueInfo::Index( myResidueIndex )) : "<< myBiopolymer.updResidue(ResidueInfo::Index( myResidueIndex )).getPdbResidueNumber()<< myBiopolymer.updResidue(ResidueInfo::Index( myResidueIndex )).getPdbInsertionCode()<<endl;   
        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" getResidueID(myResidueIndex) = >"<<getResidueID(myResidueIndex).outString()<<"< "<<endl; 
    }

    //cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" About to validate residue numbers and insertion codes  "<<endl;
    //validateResidueNumbersAndInsertionCodes();
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Just finished setting residueID's and insertion codes for chain "<<getChainID()<<" from biopolymer in input structure file"<<endl;
    //printBiopolymerInfo();
    
};


void BiopolymerClass::setResidueIDsAndInsertionCodesFromBiopolymer(const Biopolymer & inputBiopolymer, Mutation myInsertion, bool endCaps = 0 ){
    bool residueInserted = false;
    //cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Seeking to renumber, after inserting a new residue numbered: >"<<myInsertion.getResidueID().outString()<<"< "<<endl;
    if (residueIDVector.size() > 0) 
    {
        ErrorManager::instance<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Why does residueIDVector have something in it already?"<<endl; 
        ErrorManager::instance.treatError();
    } 
    // Should make sure this works for insertions at the beginning of the chain, esp. for the case of endCaps = true.
    //cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Setting residue numbers and insertion codes from biopolymer in input structure file, for chain "<<getChainID()<<endl;
    //int myResidueIndex = -1111;
    int firstInputResidueIndex = 0 + endCaps;
    int myFirstResidueIndex = firstInputResidueIndex - endCaps + proteinCapping - 1 ; // This is the index that counts over the NEW biopolymer, WITH the insertion. The -1 is just to let the ++ operation below get it to its starting value.
    int myResidueIndex = myFirstResidueIndex; 

    //ResidueInfo myResidueInfo = myBiopolymer.updResidue(ResidueInfo::Index( myResidueIndex + 1));
    for (int inputResidueIndex = ( firstInputResidueIndex ); inputResidueIndex < (inputBiopolymer.getNumResidues() - endCaps) ; inputResidueIndex ++) {
        myResidueIndex++;
        const ResidueInfo inputResidueInfo = inputBiopolymer.getResidue(ResidueInfo::Index( inputResidueIndex ));
        const int  inputResidueNumber    = (inputResidueInfo).getPdbResidueNumber();
        const char inputInsertionCode    = (inputResidueInfo).getPdbInsertionCode();
        ResidueID inputResidueID((inputResidueInfo).getPdbResidueNumber(), (inputResidueInfo).getPdbInsertionCode());
        //const char inputOneLetterCode  = (inputResidueInfo).getOneLetterCode();
        ResidueInfo myResidueInfo = myBiopolymer.updResidue(ResidueInfo::Index( myResidueIndex ));
        
        ////// Tbis was previously done in loadResidueIDVector().  However most of the time is spent in the myBiopolymer.updResidue step.  So for efficiency I am now doing it here: //// 
        ResidueID myResidueID;    
        myResidueID.setResidueNumber  (myResidueInfo.getPdbResidueNumber());
        myResidueID.setInsertionCode  (myResidueInfo.getPdbInsertionCode());
        //cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<"  (myInsertion.getResidueID() <  inputResidueID), for "<<myInsertion.getResidueID().outString()<<" and "<<inputResidueID.outString()<<"  = "<< (myInsertion.getResidueID() <  inputResidueID) <<endl;
        if ((myInsertion.getResidue() <  inputResidueID ) && (! (residueInserted)))
        {
            // Note that an insertion will shift residue indices in myBiopolymer(with insertion) compared to inputBiopolymer (insertionless).  
            // In this case, the insertion will affect indices in the NEXT inputResidueIndex, not this one.
            cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Inserting residue number: >"<<myInsertion.getResidue().getResidueNumber() <<"< "<<endl;
            cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Inserting insertion code: >"<<myInsertion.getResidue().getInsertionCode() <<"< "<<endl;
        myBiopolymer.updResidue(ResidueInfo::Index( myResidueIndex )).setPdbResidueNumber(myInsertion.getResidue().getResidueNumber() );
        myBiopolymer.updResidue(ResidueInfo::Index( myResidueIndex )).setPdbInsertionCode(myInsertion.getResidue().getInsertionCode() );
            ResidueInfo myResidueInfo2 = myBiopolymer.updResidue(ResidueInfo::Index( myResidueIndex ));
        //cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" NEW chain Residue index , number, insertion code, and residue type: "<<myResidueIndex<<", "<<myResidueInfo2.getPdbResidueNumber()   <<","<<myResidueInfo2.getPdbInsertionCode()   <<", "<<myResidueInfo2.getOneLetterCode()   <<endl;
            //cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" About to push back >"<<myInsertion.getResidueID().outString()<<"< "<< endl;
            residueIDVector.push_back(myInsertion.getResidue()); 
            residueInserted = true; myResidueIndex ++; // Note we increment here AFTER adding the INSERTED residue.
        }
        // set residue ID and insertion code in myBiopolymer, then add to residueIDVector.
        myBiopolymer.updResidue(ResidueInfo::Index( myResidueIndex  )).setPdbResidueNumber(inputResidueNumber);
        myBiopolymer.updResidue(ResidueInfo::Index( myResidueIndex  )).setPdbInsertionCode(inputInsertionCode);
        //ResidueInfo myResidueInfo3 = myBiopolymer.updResidue(ResidueInfo::Index( myResidueIndex ));
        //cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" NEW chain Residue index , number, insertion code, and residue type: "<<myResidueIndex<<"," <<myResidueInfo.getPdbResidueNumber()   <<","<<myResidueInfo.getPdbInsertionCode()   <<", "<<myResidueInfo.getOneLetterCode()   <<endl;
        
        //cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" About to push back >"<<inputResidueID.outString()<<"< "<<endl;
        residueIDVector.push_back(inputResidueID); 

    //for (int i = 0; i < residueIDVector.size(); i++) {
    //    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" residue ID : "<<residueIDVector[i].outString()<<endl;
    //}

    }
    if (myResidueIndex == myFirstResidueIndex) { // This would indicate myResidueIndex never changed!
        ErrorManager::instance<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Unexplained error!"<<endl; 
        ErrorManager::instance.treatError();
    }
    if (! (residueInserted)) { // It would appear that the insertion is to be added to the end of the chain:
            residueInserted = true; myResidueIndex++;
        myBiopolymer.updResidue(ResidueInfo::Index( myResidueIndex  )).setPdbResidueNumber(myInsertion.getResidue().getResidueNumber() );
        myBiopolymer.updResidue(ResidueInfo::Index( myResidueIndex  )).setPdbInsertionCode(myInsertion.getResidue().getInsertionCode() );
            //myResidueInfo = myBiopolymer.updResidue(ResidueInfo::Index( myResidueIndex ));
            //cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" NEW chain Residue index , number, insertion code, and residue type: "<<myResidueIndex<<"," <<myResidueInfo.getPdbResidueNumber()   <<","<<myResidueInfo.getPdbInsertionCode()   <<", "<<myResidueInfo.getOneLetterCode()   <<endl;
            cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" About to push back >"<<myInsertion.getResidue().outString()<<"< "<<endl;
            residueIDVector.push_back(myInsertion.getResidue()); 
    }
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" About to list the residue IDs for this chain, having  completed the insertion operation : "<<endl;
    for (int i = 0; i < residueIDVector.size(); i++) {
        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" residue ID : "<<residueIDVector[i].outString()<<endl;
    }
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" About to validate residue numbers and insertion codes  "<<endl;
    validateResidueNumbersAndInsertionCodes();
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Just finished setting residueID's and insertion codes for chain "<<getChainID()<<" from biopolymer in input structure "<<endl;
    
};

void BiopolymerClass::setResidueIDsAndInsertionCodesFromBiopolymerWithDeletion(const Biopolymer & oldBiopolymer, ResidueInfo::Index  myDeletionIndex, bool endCaps = 0 ){
    if (residueIDVector.size() > 0) 
    {
        ErrorManager::instance<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Why does residueIDVector have something in it already?"<<endl; 
        ErrorManager::instance.treatError();
    } 
    // Should make sure this works for insertions at the beginning of the chain, esp. for the case of endCaps = true.
    int firstInputResidueIndex = 0 + endCaps;
    int myFirstResidueIndex = firstInputResidueIndex - endCaps + proteinCapping - 1 ; // This is the index that counts over the NEW biopolymer, WITH the insertion. The -1 is just to let the ++ operation below get it to its starting value.
    int myResidueIndex = myFirstResidueIndex;  // this counts over residues in the NEW chain

    //ResidueInfo myResidueInfo = myBiopolymer.updResidue(ResidueInfo::Index( myResidueIndex + 1));
    for (int oldResidueIndex = ( firstInputResidueIndex ); oldResidueIndex < (oldBiopolymer.getNumResidues() - endCaps) ; oldResidueIndex ++) {
        if (oldResidueIndex == myDeletionIndex) { 
            // Do nothing.  This is the residue to be deleted, and therefore ignored here.
        } else {
            myResidueIndex++; 
        const ResidueInfo oldResidueInfo = oldBiopolymer.getResidue(ResidueInfo::Index( oldResidueIndex )); // This is from the OLD chain
        const int  oldResidueNumber    = (oldResidueInfo).getPdbResidueNumber();
        const char oldInsertionCode    = (oldResidueInfo).getPdbInsertionCode();
        ResidueID oldResidueID((oldResidueInfo).getPdbResidueNumber(), (oldResidueInfo).getPdbInsertionCode());
        ResidueInfo myResidueInfo = myBiopolymer.updResidue(ResidueInfo::Index( myResidueIndex ));
        ////// Tbis was previously done in loadResidueIDVector().  However most of the time is spent in the myBiopolymer.updResidue step.  So for efficiency I am now doing it here: //// 
        ResidueID myResidueID;    
        myResidueID.setResidueNumber  (myResidueInfo.getPdbResidueNumber());
        myResidueID.setInsertionCode  (myResidueInfo.getPdbInsertionCode());
        // set residue ID and insertion code in myBiopolymer, then add to residueIDVector.
        myBiopolymer.updResidue(ResidueInfo::Index( myResidueIndex  )).setPdbResidueNumber(oldResidueNumber);
        myBiopolymer.updResidue(ResidueInfo::Index( myResidueIndex  )).setPdbInsertionCode(oldInsertionCode);
        //myResidueInfo = myBiopolymer.updResidue(ResidueInfo::Index( myResidueIndex ));
        residueIDVector.push_back(oldResidueID); 
        }

    }
    if (myResidueIndex == myFirstResidueIndex) { // This would indicate myResidueIndex never changed!
        ErrorManager::instance<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Unexplained error! Did you delete the last residue in the chain?"<<endl; 
        ErrorManager::instance.treatError();
    }
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" About to list the residue IDs for this chain, having  completed the insertion operation : "<<endl;
    for (int i = 0; i < residueIDVector.size(); i++) {
        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" residue ID : "<<residueIDVector[i].outString()<<endl;
    }
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" About to validate residue numbers and insertion codes  "<<endl;
    validateResidueNumbersAndInsertionCodes();
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Just finished setting residueID's and insertion codes for chain "<<getChainID()<<" from biopolymer in old structure "<<endl;
};

void BiopolymerClass::printBiopolymerInfo() {
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" printing biopolymer information for chain "<< getChainID()<<" first residue : "<<getFirstResidueID().outString()<< " and last residue : " << getLastResidueID().outString()<<endl;
    for (int i = 0; i < getChainLength(); i++) {
        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Chain ID, Residue type, number, and insertion code: "<<getChainID()<<", "  <<myBiopolymer.getResidue(ResidueInfo::Index(i)).getOneLetterCode() <<", "<<myBiopolymer.getResidue(ResidueInfo::Index(i)).getPdbResidueNumber()<<", "<<myBiopolymer.getResidue(ResidueInfo::Index(i)).getPdbInsertionCode()<<endl;
    } 
}





/*void BiopolymerClass::writeMutationBackboneRigidifier(std::ofstream & output, const int offset) {
                int leftFlexibleOffset = offset;
                int rightFlexibleOffset = offset;
                for (int i = 0 ; i < getNumMutationVectorElements(); i++) {
                        output <<"mobilizer Default "<<mutationVector[i].getChainID()<<" " ;
                        output<<safeSum(mutationVector[i].getResidueID(),(- leftFlexibleOffset)).outString()<<" ";
            output <<safeSum(mutationVector[i].getResidueID(),rightFlexibleOffset).outString()<<std::endl;
                        //for (ResidueID myResidueID = (safeSum(mutationVector[i].getResidueID(),(- leftFlexibleOffset))) ; 
              //myResidueID <= safeSum(mutationVector[i].getResidueID(),rightFlexibleOffset);
              //incrementResidueID(myResidueID))   
                        ResidueID myResidueID = (safeSum(mutationVector[i].getResidueID(),(- leftFlexibleOffset))) ; 

            while(  myResidueID <= safeSum(mutationVector[i].getResidueID(),rightFlexibleOffset))
            {
                output<<"singleBondMobility "<<getChainID()<<" "<<myResidueID.outString()<<" N Rigid "<<getChainID()<<" "<<myResidueID.outString()<<" CA "<<endl;
                output<<"singleBondMobility "<<getChainID()<<" "<<myResidueID.outString()<<" CA Rigid "<<getChainID()<<" "<<myResidueID.outString()<<" C "<<endl;
                if (myResidueID < safeSum(mutationVector[i].getResidueID(),rightFlexibleOffset)){ // We have to be careful not to increment myResidueID past the end of the chain
                    incrementResidueID(myResidueID) ;  }
                else if (myResidueID == safeSum(mutationVector[i].getResidueID(),rightFlexibleOffset)) {
                    break;
                }
            }
                }
}*/


/*void BiopolymerClass::writePhysicsZones(std::ofstream & output, const int offset) {
                int leftFlexibleOffset = offset;
                int rightFlexibleOffset = offset;
                for (int i = 0 ; i < getNumMutationVectorElements(); i++) {
                        output <<"includeAllResiduesWithin 1.2 "<<mutationVector[i].getChainID()<<" " ;
                        output <<safeSum(mutationVector[i].getResidueID(),(- leftFlexibleOffset)).outString() <<endl;
                        output <<"includeAllResiduesWithin 1.2 "<<mutationVector[i].getChainID()<<" " ;
            output <<safeSum(mutationVector[i].getResidueID(),rightFlexibleOffset).outString()<<std::endl;
                        output <<"includeAllResiduesWithin 1.2 "<<mutationVector[i].getChainID()<<" "<<mutationVector[i].getResidueID().outString() <<std::endl;
                }
}*/



bool BiopolymerClass::residueIDLessThanOrEqualTo(ResidueID  residueA, ResidueID  residueB){
    return (getResidueIndex(residueA) <= getResidueIndex(residueB));
}

bool BiopolymerClass::residueIDGreaterThanOrEqualTo(ResidueID  residueA, ResidueID  residueB){
    return (getResidueIndex(residueA) >= getResidueIndex(residueB));
};

ResidueID BiopolymerClass::incrementResidueID(ResidueID  & residueID){
    if (residueID == getLastResidueID()) {
        ErrorManager::instance <<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" You cannot increment the last residue ID!"<<endl;
        ErrorManager::instance.treatError();
    }
    residueID = getResidueID( getResidueIndex(residueID) + 1);   
    validateResidueID(residueID); // getResidueIndex (residueID) should get the index directly from myBiopolymer, but this is one higher .. needs validation.
    return residueID;
};

ResidueID BiopolymerClass::decrementResidueID(ResidueID & residueID){

    if (residueID == getFirstResidueID()) {
        ErrorManager::instance <<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" You cannot decrement the first residue ID!"<<endl;
        ErrorManager::instance.treatError();
    }
    residueID = getResidueID( getResidueIndex(residueID) - 1);   
    validateResidueID(residueID);
    return residueID;
};

void        BiopolymerClass::setDefaultPhiAngle (ResidueID residueID, Angle phi) {
    //myBiopolymer.updResidue(getResidueIndex(residueID)).setDefaultPhiAngle(phi);
        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Setting angle connecting C of residue "<<((sum(residueID , -1)).outString())<<", and N,CA,C of residue "<<(residueID.outString())<<endl;
        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Residue ID : "<<residueID.outString()<<" atom N has path name:  >"<< (atomPathString((residueID), String("N")))<<"< "<<endl;
        myBiopolymer.setDefaultDihedralAngle(phi, 
        atomPathString(sum(residueID , -1), String("C")),
        atomPathString((residueID), String("N")),
        atomPathString((residueID), String("CA")),
        atomPathString((residueID), String("C"))
                );
};

void        BiopolymerClass::setDefaultPsiAngle (ResidueID residueID, Angle psi){
    //myBiopolymer.updResidue(getResidueIndex(residueID)).setDefaultPsiAngle(psi);
        myBiopolymer.setDefaultDihedralAngle(psi, 
        atomPathString((residueID), String("N")),
        atomPathString((residueID), String("CA")),
        atomPathString((residueID), String("C")),
        atomPathString(sum(residueID , 1), String("N"))
        );
};
void        BiopolymerClass::setDefaultPeptideDihedralAngle (ResidueID residueID1, ResidueID residueID2, Angle dihedral){
    //myBiopolymer.updResidue(getResidueIndex(residueID)).setDefaultPsiAngle(psi);
    if (!( sum(residueID1,1) == residueID2)){
            ErrorManager::instance <<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" You must use consecutive increaseing residue ID's here"<<endl;
            ErrorManager::instance.treatError();
        }
        // Prolines don't have HN.  So we have to measure dihedral with respect to CA, which is 180 degrees from HN.
        Angle myDihedral = dihedral - (180*Deg2Rad);
        myBiopolymer.setDefaultDihedralAngle(myDihedral,
        atomPathString((residueID1), String("O"  )),
        atomPathString((residueID1), String("C")),
        atomPathString((residueID2), String("N")),
        atomPathString((residueID2), String("CA"))
        );
};

void        BiopolymerClass::setAlphaHelicalDefaultBackboneAngles(ResidueID startResidue, ResidueID endResidue){
        validateResidueID(sum(startResidue, -1));
        validateResidueID(sum(endResidue , 1));
        // these defaults are from (http://www.chembio.uoguelph.ca/educmat/phy456/456lec01.htm)
    Angle alphaHelicalPhi = -60*Deg2Rad; //Angle is in radians
    Angle alphaHelicalPsi = -60*Deg2Rad;
    Angle transAngle = 180*Deg2Rad ;
    for (ResidueID myResidueID = startResidue; myResidueID <= endResidue; incrementResidueID(myResidueID)){
                setDefaultPeptideDihedralAngle (sum(myResidueID,-1),myResidueID, transAngle );
        setDefaultPhiAngle (myResidueID, alphaHelicalPhi);
        setDefaultPsiAngle (myResidueID, alphaHelicalPsi);
                setDefaultPeptideDihedralAngle (myResidueID,sum(myResidueID,1),transAngle );
    }
};
void        BiopolymerClass::setParallelBetaSheetDefaultBackboneAngles(ResidueID startResidue, ResidueID endResidue){
        validateResidueID(sum(startResidue, -1));
        validateResidueID(sum(endResidue , 1));
        // these defaults come from JE Wampler, 1996 (http://www.bmb.uga.edu/wampler/tutorial/prot2.html)
    Angle myPhi = -119*Deg2Rad; //Angle is in radians
    Angle myPsi = 113*Deg2Rad;
    Angle transAngle = 180*Deg2Rad ;
    for (ResidueID myResidueID = startResidue; myResidueID <= endResidue; incrementResidueID(myResidueID)){
                setDefaultPeptideDihedralAngle (sum(myResidueID,-1),myResidueID, transAngle );
        setDefaultPhiAngle (myResidueID, myPhi);
        setDefaultPsiAngle (myResidueID, myPsi);
                setDefaultPeptideDihedralAngle (myResidueID,sum(myResidueID,1),transAngle );
    }
};
void        BiopolymerClass::setAntiParallelBetaSheetDefaultBackboneAngles(ResidueID startResidue, ResidueID endResidue){
        validateResidueID(sum(startResidue, -1));
        validateResidueID(sum(endResidue , 1));
        // these defaults come from JE Wampler, 1996 (http://www.bmb.uga.edu/wampler/tutorial/prot2.html)
    Angle myPhi = -139*Deg2Rad; //Angle is in radians
    Angle myPsi = 136*Deg2Rad;
    Angle transAngle = 180*Deg2Rad ;
    for (ResidueID myResidueID = startResidue; myResidueID <= endResidue; incrementResidueID(myResidueID)){
                setDefaultPeptideDihedralAngle (sum(myResidueID,-1),myResidueID, transAngle );
        setDefaultPhiAngle (myResidueID, myPhi);
        setDefaultPsiAngle (myResidueID, myPsi);
                setDefaultPeptideDihedralAngle (myResidueID,sum(myResidueID,1),transAngle );
    }
};

int BiopolymerClass::difference(ResidueID  residueA, ResidueID  residueB ){
    return (getResidueIndex(residueA) - getResidueIndex(residueB));
};

// This test to make sure that it is possible to sum increment to inputResidueID and not go out of range.
bool BiopolymerClass::safeSum(ResidueID  inputResidueID, int  increment, ResidueID outputResidueID){
    if ((getResidueIndex(inputResidueID) + increment) >  getResidueIndex(getLastResidueID())) {
    outputResidueID = getLastResidueID();
    return false;
    } else if ((getResidueIndex(inputResidueID) + increment) < getResidueIndex(getFirstResidueID())) {
    outputResidueID = getFirstResidueID();
    return false;
    } else {
        outputResidueID = residueIDVector[(getResidueIndex(inputResidueID) + increment)];
    return true;
    }
}

// This version of safeSum returns a ResidueID which is bounded, i.e. at least FirstResidue and at most LastResidue.
ResidueID BiopolymerClass::safeSum(ResidueID  inputResidueID, int  increment){
    if ((getResidueIndex(inputResidueID) + increment) >  getResidueIndex(getLastResidueID())) {
    return getLastResidueID();
    } else if ((getResidueIndex(inputResidueID) + increment) < getResidueIndex(getFirstResidueID())) {
    return getFirstResidueID();
    } else {
        return residueIDVector[(getResidueIndex(inputResidueID) + increment)];
    }
}

void BiopolymerClass::setCurrentSequencesFromOriginalSequences() {
    setSequence(getOriginalSequence()); // Note that this will not have the right PDB residue numbering. Hence the next line:
    setPdbResidueNumbersFromResidueIDVector();
}




ResidueID BiopolymerClass::sum(ResidueID  oldResidueID, int  increment ){
    ResidueID newResidueID = oldResidueID;
    if (residueIDVector.size() > 0) {
        int oldResidueIndex = getResidueIndex(oldResidueID);
        int newResidueIndex = oldResidueIndex + increment;
        if ((newResidueIndex < 0) || (newResidueIndex >= residueIDVector.size())) {
            ErrorManager::instance <<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" You have tried to add "<<increment<<" to residue ID "<<oldResidueID.outString()<<".  The result is out of range."<<endl; ErrorManager::instance.treatError();
        } else {
            newResidueID = getResidueID(newResidueIndex);
            return(newResidueID);
        }
            
    } else {
        ErrorManager::instance <<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" This really should be obsolete now.  If you get this message, something has gone wrong."<<endl; ErrorManager::instance.treatError();
        if (increment >=0)
            for (int i = 0 ; i < increment; i++)
                newResidueID = incrementResidueID(newResidueID);
        else
            for (int i = 0 ; i > increment; i--)
                newResidueID = decrementResidueID(newResidueID);
    }   
    return newResidueID;
};

void BiopolymerClass::setPdbFileName(String pdbFileName){
    this->pdbFileName = pdbFileName;
}
String BiopolymerClass::getPdbFileName(){
    return pdbFileName;
}

void  BiopolymerClass::setPdbStructure(const PdbStructure myPdbStructure)
{
    this->pdbStructure = myPdbStructure;
}
const PdbStructure BiopolymerClass::getPdbStructure()
{
    return this->pdbStructure;
}

void BiopolymerClass::setLoadFromPdb(bool yesno){
    this->loadFromPdb = yesno;
}
bool BiopolymerClass::getLoadFromPdb(){
    return loadFromPdb;
}

bool BiopolymerClass::hasResidueStretch(ResidueStretch & residues)
{
    if(residues.getChain() != getChainID())
    {
        residues.printStretch();
        cout <<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" : The ResidueStretch above does not belong to Biopolymer "<< getChainID()<<endl; 
        return false;
    }

    if(residues.getStartResidue() < getFirstResidueID() || residues.getStartResidue() > getLastResidueID() ||
       residues.getEndResidue() < getFirstResidueID() || residues.getEndResidue() > getLastResidueID() )
    {
        residues.printStretch();
        cout <<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" : The ResidueID range of the ResidueStretch above should be included in ["<< getFirstResidueID().outString()<< ":"<< getLastResidueID().outString()<<"]."<<endl; 
        return false;
    }

    return true;
}


// void BiopolymerClass::AddInactiveResidues(ResidueStretch & residues)
// {

//     if(!hasResidueStretch(residues))
//     {
//         ErrorManager::instance << __FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" : AddInactiveResidues : You tried to add an invalid residue stretch." << endl;
//         ErrorManager::instance.treatError();
//     }
//     inactiveResidueStretches.addResidueStretchToVector(residues);
// }


// void BiopolymerClass::RemoveInactiveResidues(ResidueStretch & residues)
// {
//     if(!hasResidueStretch(residues))
//     {
//         ErrorManager::instance << __FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" : RemoveInactiveResidues : You tried to remove an invalid residue stretch." << endl;
//         ErrorManager::instance.treatError();
//     }
//     inactiveResidueStretches.removeResidueStretchFromVector(residues);
// }

// vector<ResidueID> BiopolymerClass::getInactiveResiduesVector()
// {
//     vector<ResidueID> vec;
//     vector<ResidueStretch> residueStretchVector = inactiveResidueStretches.getResidueStretchVector();
//     vector<ResidueStretch>::iterator it;
//     for(it = residueStretchVector.begin(); it != residueStretchVector.end(); it++)
//     {
//         ResidueID indexResidueID = getFirstResidueID();
//         while(indexResidueID < it->getStartResidue()) { 
//             vec.push_back(updResidue(getResidueIndex(indexResidueID)));
//             if(indexResidueID <  getLastResidueID()) 
//                 incrementResidueID(indexResidueID); 
//             else 
//                 break; // make sure we don't increment past the last residue
//         }
//     }

//     return vec;
// }

void BiopolymerClass::setActivePhysics(bool yesno){
    activePhysics = yesno;
}
bool BiopolymerClass::getActivePhysics() const{
    return activePhysics;
}

#ifndef USE_OPENMM
// Outdated, it shuld'nt be used
vector<ResidueID> BiopolymerClass::getResiduesWithin(Vec3 location, double distance){
    vector<ResidueID> residuesWithin;
    for (ResidueID j = getFirstResidueID(); j <= getLastResidueID(); incrementResidueID(j)) {
        double myDistance = (double)(calcDefaultAtomLocationInGroundFrame(j, getRepresentativeAtomName()) - location).norm();
        if(myDistance <= distance){
            residuesWithin.push_back(j);
        }

        if(j == getLastResidueID()) break;
    }
    return residuesWithin;
}
#endif

/**
 * /brief This method locks all MobilizedBody's in a biopolymer. It is equivalent to using BondMobility::Rigid, but with constraints rather than mobilizers. It is intended to be used for adaptive dynamics, because we will be able to monitor the reaction forces required to maintain the constraints. These forces can be the criterion for "melting" a DOF.
 *
 */
/*
    void BiopolymerClass::lockBiopolymerMobilizedBodies (){
        for (int i = 0; i < atomInfoVector.size(); i++) {
            atomInfoVector[i].mobilizedBody.lockMobilizer();  
        }
    }
 */




/*
template<class ResidueStretchType>
void BiopolymerClass::selectivelyRemoveResidueStretchFromContainer(ResidueStretch & residueStretch, ResidueStretchContainer <ResidueStretchType> & residueStretchContainer){
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
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" the Default stretch is :"<<endl;
    residueStretch.printStretch();
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Now checking "<<residueStretchContainer.getNumResidueStretches()<<" stretches: "<<endl;
    for (int i = 0; i < residueStretchContainer.getNumResidueStretches(); i++)
    {
        residueStretchContainer.residueStretchVector[i].printStretch();

        if (residueStretchContainer.residueStretchVector[i].getChain().compare((residueStretch.getChain() )) != 0) {continue;} // in other words, only make modificatiosn to residueStretchContainer if chain ID's match.
        else if ((residueStretch.getStartResidue() <= residueStretchContainer.residueStretchVector[i].getStartResidue()) &&
            (residueStretch.getEndResidue() >= residueStretchContainer.residueStretchVector[i].getEndResidue()))
           {   //case = 1
               residueStretchContainer.residueStretchVector.erase(residueStretchContainer.residueStretchVector.begin() + i);
               i--; // vector has been shortened, so make sure we don't skip the next residueStretchContainer.residueStretchVector[i].
               if (i < -1) {ErrorManager::instance <<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<"Unexplained error!"<<endl; ErrorManager::instance.treatError();}
           }
        else if ((residueStretch.getStartResidue() >  residueStretchContainer.residueStretchVector[i].getStartResidue()) &&
            (residueStretch.getEndResidue() <  residueStretchContainer.residueStretchVector[i].getEndResidue()))
           {   // case = 2 ;
               cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<"  "<<endl;
               MobilizerStretch secondResidueStretch = residueStretchContainer.residueStretchVector[i];
               ResidueID tempStartResidueID = (residueStretch).getStartResidue(); // getStartResidue() returns a temporary, whereas decrementResidueID expects a reference. can't convert a temporary to a reference.  This is because decrementResidueID might (and will!) try to modify ResidueID (as the name of the function suggests!).
               //residueStretchContainer.residueStretchVector[i].setEndResidue(decrementResidueID((residueStretch).getStartResidue() ));
               residueStretchContainer.residueStretchVector[i].setEndResidue(decrementResidueID(tempStartResidueID));//((residueStretch).getStartResidue() )));
               cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Just decreased endpoint of stretch "<<i<<".  New stretch is:"<<endl;
               residueStretchContainer.residueStretchVector[i].printStretch();
               ResidueID tempEndResidueID = (residueStretch).getEndResidue();
               secondResidueStretch.setStartResidue(incrementResidueID(tempEndResidueID));//  residueStretch.getEndResidue()));
               residueStretchContainer.addResidueStretchToVector(secondResidueStretch);
               cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Just added new  stretch :"<<endl;
               residueStretchContainer.residueStretchVector[residueStretchContainer.getNumResidueStretches()-1].printStretch();
               cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Moving on to check next stretch. "<<endl;


           }
        else if ((residueStretch.getStartResidue() == residueStretchContainer.residueStretchVector[i].getStartResidue()) &&
            (residueStretch.getEndResidue() <  residueStretchContainer.residueStretchVector[i].getEndResidue()))
           {   // case = 3;

               cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<"  Case 3"<<endl;
               ResidueID tempEndResidueID = (residueStretch).getEndResidue();
               residueStretchContainer.residueStretchVector[i].setStartResidue(incrementResidueID(tempEndResidueID));//residueStretch.getEndResidue() ))  ;
           }
        else if ((residueStretch.getEndResidue() == residueStretchContainer.residueStretchVector[i].getEndResidue()) &&
            (residueStretch.getStartResidue() >  residueStretchContainer.residueStretchVector[i].getStartResidue()))
           {   // case = 4;
               cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<"  Case 4"<<endl;

               ResidueID tempStartResidueID = (residueStretch).getStartResidue();
               residueStretchContainer.residueStretchVector[i].setEndResidue(decrementResidueID(tempStartResidueID));//residueStretch.getStartResidue()));
           }
        else if ((residueStretch.getStartResidue() <   residueStretchContainer.residueStretchVector[i].getStartResidue()) &&
            (residueStretch.getEndResidue()        >=  residueStretchContainer.residueStretchVector[i].getStartResidue()) &&
                 (residueStretch.getEndResidue()        <   residueStretchContainer.residueStretchVector[i].getEndResidue()))
        {   // case = 5;
            cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<"  Case 5"<<endl;

            ResidueID tempEndResidueID = (residueStretch).getEndResidue();
            residueStretchContainer.residueStretchVector[i].setStartResidue(incrementResidueID(tempEndResidueID));//residueStretch.getEndResidue()))  ;
        }
        else if ((residueStretch.getEndResidue() >  residueStretchContainer.residueStretchVector[i].getEndResidue()) &&
                 (residueStretch.getStartResidue() >  residueStretchContainer.residueStretchVector[i].getStartResidue())     &&
                 (residueStretch.getStartResidue() <=  residueStretchContainer.residueStretchVector[i].getEndResidue()))
        {    // case = 6;
            cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<"  Case 6"<<endl;

            ResidueID tempStartResidueID = (residueStretch).getStartResidue();
            residueStretchContainer.residueStretchVector[i].setEndResidue(decrementResidueID(tempStartResidueID));//  residueStretch.getStartResidue()));
        }
        else if (residueStretch.getEndResidue() < residueStretchContainer.residueStretchVector[i].getStartResidue()) {} // do nothing, stretches are disjoint
        else if (residueStretch.getStartResidue() > residueStretchContainer.residueStretchVector[i].getEndResidue()) {} // do nothing, stretches are disjoint
        else {
            // this should never happen
                    ErrorManager::instance <<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<"Unexplained error!"<<endl; ErrorManager::instance.treatError();
            }
    }
}*/


TAlign BiopolymerClass::createGappedAlignment(BiopolymerClass otherBiopolymerClass, double alignmentForcesGapPenalty ){ // Set a default value of -1 for the gap penalty to allow gaps. For ungapped, do a big value e.g. -1000
    //String chainA = thread.chainID1;
    //String chainB = thread.chainID2;
    //BiopolymerClass & bpA = myBiopolymerClassContainer.updBiopolymerClass(chainA);
    //BiopolymerClass & bpB = myBiopolymerClassContainer.updBiopolymerClass(chainB);
    //typedef seqan::String<char> TSequence;                 // sequence type
    //typedef seqan::Align<TSequence,seqan::ArrayGaps> TAlign;      // align type
    TSequence seqA = getSubSequence(getFirstResidueID(),getLastResidueID()  ).c_str();  // Need a new BiopolymerClass method which retrieves subsequences.!
    TSequence seqB = otherBiopolymerClass.getSubSequence(otherBiopolymerClass.getFirstResidueID(), otherBiopolymerClass.getLastResidueID() ).c_str();
    TAlign align;
    seqan::resize(rows(align), 2);
    assignSource(row(align,0),seqA);
    assignSource(row(align,1),seqB);
    // simple alignment:
    int score = globalAlignment(align, seqan::Score<int,seqan::Simple>(0,-1, alignmentForcesGapPenalty )); // ..signature:Score<TValue, Simple>(match, mismatch, gap [, gap_open])
    return align;
}


int BiopolymerClass::getCorrespondingMutationInCurrentBiopolymer(BiopolymerClass otherBiopolymerClass, TAlign align,Mutation mutationInOtherBiopolymer, Mutation & mutationInCurrentBiopolymer){
    //Mutation mutationInCurrentBiopolymer;
    std::string chainInCurrentBiopolymer = getChainID();
    ResidueID residueIdInCurrentBiopolymer ;
    if (getCorrespondingResidueInCurrentBiopolymer(otherBiopolymerClass, align, mutationInOtherBiopolymer.getResidue(), residueIdInCurrentBiopolymer)) { // Return value of 0 indicates success.
       cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Failed to  translate the residue  "<<mutationInOtherBiopolymer.getResidue().outString() <<" from other biopolymer, with chain = "<< otherBiopolymerClass.getChainID()  <<" to the current biopolymer, with chain ="<< chainInCurrentBiopolymer <<  endl;
       return 1; // return value of 1    indicates failure, .
    } else {
       cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Successfully translated the other biopolymer, with chain = "<< otherBiopolymerClass.getChainID()  <<" residue  "<<mutationInOtherBiopolymer.getResidue().outString() <<" to residue ID: " << residueIdInCurrentBiopolymer.outString()<<" in the current biopolymer, with chain ="<<chainInCurrentBiopolymer  <<  endl;
    };
    std::string substitutedResidueTypeInCurrentBiopolymer = mutationInOtherBiopolymer. getSubstitutedResidueType();
    mutationInCurrentBiopolymer.setChain    (chainInCurrentBiopolymer)    ;
    mutationInCurrentBiopolymer.setResidue  (residueIdInCurrentBiopolymer);
    mutationInCurrentBiopolymer.setSubstitutedResidueType(substitutedResidueTypeInCurrentBiopolymer)    ;
    //return mutationInCurrentBiopolymer    ;
    return 0 ; // Return value of 0 to indicate normal (successful) function. Calling program should monitor this, a return value of 1 or nonzero indicates failure.

}

// The return value of this method is 0 for failure, 1 for success.
int  BiopolymerClass::getCorrespondingResidueInCurrentBiopolymer(BiopolymerClass otherBiopolymerClass, TAlign align, ResidueID residueIdInOtherBiopolymerClass , ResidueID & correspondingResidueIdInCurrentBiopolymerClass  ){


        int status = 1    ; // set to 1 to indicate failure
        int aIndex = 0; int bIndex = 0; // Indices which count over residues in chains A and B.

        int i  = 0; // counts over columns in alignment
        correspondingResidueIdInCurrentBiopolymerClass = ResidueID(-11111,' '); // Set to a ridiculous value to make error easier to detect
        //ResidueID correspondingResidueIdInCurrentBiopolymerClass("-11111"," ");
        ResidueID tempResidueIdInOtherBiopolymerClass;
        //while ((aIndex < getSubSequence(thread.residueStart1,thread.residueEnd1 ).length()) &&
        //       (bIndex < otherBiopolymerClass.getSubSequence(thread.residueStart2, thread.residueEnd2 ).length()   ))
        while ((aIndex < getSubSequence(getFirstResidueID(),getLastResidueID() ).length()) &&
               (bIndex < otherBiopolymerClass.getSubSequence(otherBiopolymerClass.getFirstResidueID() , otherBiopolymerClass.getLastResidueID () ).length()   ))
        {
            tempResidueIdInOtherBiopolymerClass = otherBiopolymerClass.sum(otherBiopolymerClass.getFirstResidueID(), bIndex);
            if (tempResidueIdInOtherBiopolymerClass == residueIdInOtherBiopolymerClass) // Remember, the inequalities <,> are probably not reliable now that we allow non sequential residueID's
            {
                 if (String(seqan::row (align,0)[i]).compare("-")  == 0 ) {
                     cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Error! The residue ID in the other biopolymer: "<<residueIdInOtherBiopolymerClass.outString()<<" appears to be inserted with respect to the current biopolymer. Unable to provide a corresponding residue ID"<<std::endl; //exit(1);
                     status = 1;  // The calling procedure should be monitoring this.
                     return status; // Failed! success  would be 0
                 } else if (String(seqan::row (align,1)[i]).compare("-")  == 0  ) {
                      cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Error! The residue ID in the other biopolymer: "<<residueIdInOtherBiopolymerClass.outString()<<" returned a gap in that other biopolymer.  This is a very odd error!  "<<std::endl; //exit(1);
                      status = 1;  // The calling procedure should be monitoring this.
                      return status; // Failed! success  would be 0
                      //return correspondingResidueIdInCurrentBiopolymerClass;
                 } else { // There is no error. We can now return the residue ID in the current biopolymer
                     correspondingResidueIdInCurrentBiopolymerClass = sum(getFirstResidueID()   , aIndex);
                     cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Successfully found a counterpart to the other biopolymer's residue "<<residueIdInOtherBiopolymerClass.outString()<<" . The counterpart in the current biopolymer is : "<<correspondingResidueIdInCurrentBiopolymerClass.outString()<<std::endl;
                     status = 0; // The calling procedure should be monitoring this. We have succeeded, so return 0.
                     return status;
                     //return correspondingResidueIdInCurrentBiopolymerClass;
                 }
            }


            if (String(seqan::row (align,0)[i]).compare("-")  != 0  ) {
                aIndex ++;
            }
            if (String(seqan::row (align,1)[i]).compare("-")  != 0  ) {
                bIndex ++;
            }
            i++;
        } // End While
        // If we got this far, we failed to find a match.
        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" : Failed to find a  residue in the current biopolymer, corresponding to the other biopolymer's residue :"<< residueIdInOtherBiopolymerClass.outString()<< endl;
        return 1;
         //ErrorManager::instance <<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": Failed to find a  residue in the current biopolymer, corresponding to the other biopolymer's residue :"<< residueIdInOtherBiopolymerClass.outString()<< endl;
         //ErrorManager::instance.treatError();
}

//template <class ResidueStretchType>

//template<class ResidueStretchType> 
//void BiopolymerClass::selectivelyRemoveResidueStretchFromContainer(ResidueStretch & residueStretch, ResidueStretchContainer <ResidueStretch> & residueStretchContainer){}







/////////////////////////////////////////////////////////////////////////////
/// gets rid of all the BiopolymerClass objects in its biopolymerClassMap.///
/////////////////////////////////////////////////////////////////////////////

void BiopolymerClassContainer::clear(){
    mutationVector.clear();
    pdbStructureMap.clear();
    biopolymerClassMap.clear();
    atomicPropertyOverrideVector.clear();
}




/////////////////////////////////////////////////////////////////////////////
/// initialize all biopolymers in container.                              ///
/////////////////////////////////////////////////////////////////////////////

int  BiopolymerClassContainer::initializeBiopolymers(CompoundSystem & system,
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
                                                     double myPlanarityThreshold)
{
    map<const String, BiopolymerClass>::iterator biopolymerClassMapIterator = biopolymerClassMap.begin();
    int n = 0;
    int returnValue = 0;
    for(biopolymerClassMapIterator = biopolymerClassMap.begin(); biopolymerClassMapIterator != biopolymerClassMap.end(); biopolymerClassMapIterator++) {
        returnValue = (biopolymerClassMapIterator->second).initializeBiopolymer(system, 
                                                                  myProteinCapping, matchExact, 
                                                                  matchIdealized ,matchOptimize, 
                                                                  matchHydrogenAtomLocations, 
                                                                  matchPurineN1AtomLocations, 
                                                                  guessCoordinates,
                                                                  n,
                                                                  initialSeparation, 
                                                                  displacementVector,
                                                                  matchingMinimizerTolerance,
                                                                  myPlanarityThreshold,
                                                                  secondaryStructureStretchVector
                                                                  );  
        if (returnValue){
            std::cout<<__FILE__<<":"<<__LINE__<<" Warning: Returned an error from initializeBiopolymer"<<std::endl;
            //returnValue = 1;
        }; 
        n++;
    }
    return returnValue;
}

int  BiopolymerClassContainer::initializeBiopolymer(String chainID, CompoundSystem & system,
                                                    bool myProteinCapping, bool matchExact, 
                                                    bool matchIdealized, const bool matchOptimize, 
                                                    bool matchHydrogenAtomLocations, 
                                                    bool matchPurineN1AtomLocations,
                                                    bool guessCoordinates,
                                                    double initialSeparation, 
                                                    const vector<Displacement> displacementVector,
                                                    double matchingMinimizerTolerance,
                                                    double myPlanarityThreshold,
                                                    vector<SecondaryStructureStretch> secondaryStructureStretchVector)
{
    BiopolymerClass & bpc = updBiopolymerClass(chainID);
    if (bpc.initializeBiopolymer(system, myProteinCapping, matchExact, 
                             matchIdealized, matchOptimize,
                             matchHydrogenAtomLocations, matchPurineN1AtomLocations,
                             guessCoordinates, 
                             getBiopolymerClassIndex(chainID),
                             initialSeparation, displacementVector,
                             matchingMinimizerTolerance,
                             myPlanarityThreshold,
                             secondaryStructureStretchVector)) {
            std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Returned an error from initializeBiopolymer"<<std::endl;
            return 1;
    };

}

String BiopolymerClassContainer::printOriginalAndRenumberedResidueIDs(const String myPdbId) {
    String myQuery = "";
    for (auto  biopolymerClassMapIterator = biopolymerClassMap.begin() ; biopolymerClassMapIterator != biopolymerClassMap.end(); biopolymerClassMapIterator++) {
        // argument returns the chain ID
        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Writing correspondence between original and renumbered residue IDs for chain >"<<biopolymerClassMapIterator->first<<"< ."<<endl;
        myQuery += updBiopolymerClass(biopolymerClassMapIterator->first).printOriginalAndRenumberedResidueIDs(myPdbId);
    }
    return myQuery;
}

void BiopolymerClassContainer::renumberPdbResidues(ResidueID firstResidueID) {
    for (auto  biopolymerClassMapIterator = biopolymerClassMap.begin() ; biopolymerClassMapIterator != biopolymerClassMap.end(); biopolymerClassMapIterator++) {
        // This returns the chain ID
        updBiopolymerClass(biopolymerClassMap.begin()->first).renumberPdbResidues(firstResidueID);
    }
}

void BiopolymerClassContainer::validateAtomInfoVectors(){
    map<const String, BiopolymerClass>::iterator biopolymerClassMapIterator = biopolymerClassMap.begin();
    for(biopolymerClassMapIterator = biopolymerClassMap.begin(); biopolymerClassMapIterator != biopolymerClassMap.end(); biopolymerClassMapIterator++) {
        (biopolymerClassMapIterator->second).validateAtomInfoVector(); 
    }
}



////////////////////////////////////////////////////
/// Accessor method to add a new BiopolymerClass ///
////////////////////////////////////////////////////


void BiopolymerClassContainer::addBiopolymerClass(String mySequence, String myChainID, ResidueID myFirstResidueNumber, 
                                                  String myBiopolymerType, bool proteinCapping, String pdbFileName, bool loadFromPdb, bool useNACappingHydroxyls)
{
    BiopolymerClass bp(mySequence, myChainID, myFirstResidueNumber, myBiopolymerType, proteinCapping,useNACappingHydroxyls);
    bp.setRenumberPdbResidues(0); // Default value
    biopolymerClassMap[myChainID] = bp;
    biopolymerClassMap.at(myChainID).setPdbFileName(pdbFileName);
    biopolymerClassMap.at(myChainID).setLoadFromPdb(loadFromPdb);
}
/*void BiopolymerClassContainer::addBiopolymerClass(String newChainID, BiopolymerClass newBiopolymerClass) 
{
    biopolymerClassMap[newChainID] = newBiopolymerClass;
}*/

void BiopolymerClassContainer::deleteBiopolymerClass(String myChainID ) {
    //int myIndex = getBiopolymerClassIndex(String myChainID);
    validateChainID(myChainID);
    biopolymerClassMap.erase (myChainID);
}

void BiopolymerClassContainer::deleteAllNonMutatedBiopolymerClasses(){
    map<const String, BiopolymerClass>::iterator it;
    map<const String, BiopolymerClass>::iterator next;
    next = biopolymerClassMap.begin();
    
    while (next != biopolymerClassMap.end())
    {
        bool match = false;   
	it = next;
	for (int i = 0; i < mutationVector.size(); i++) {
		if ((it->second).getChainID().compare( mutationVector[i].getChain()) == 0){
                    match = true;
                    //biopolymerClassMap.erase ((it->second).getChainID());
                } 
	}
        if (! match) {
                cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Chain "<< (it->second).getChainID() << " is never mutated. Deleting this chain."<<endl;
                biopolymerClassMap.erase ((it->second).getChainID());
        }
        next++;
    }
}


////////////////////////////////////////////
/// Fetches a non-const BiopolymerClass  ///
////////////////////////////////////////////
BiopolymerClass &   BiopolymerClassContainer::updBiopolymerClass(String myChainID) {
    validateChainID(myChainID);
    return biopolymerClassMap[myChainID];
}

int   BiopolymerClassContainer::getBiopolymerClassIndex(String myChainID){
    validateChainID(myChainID);
    map<const String, BiopolymerClass>::iterator it;
    map<const String, BiopolymerClass>::iterator next;
    int i = 0;
    next = biopolymerClassMap.begin();
    while (next != biopolymerClassMap.end())
    {
       it = next;
       if ((it->second).getChainID().compare(myChainID) == 0) return i ;
       next++;
       i++;
    }
    ErrorManager::instance <<"["<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<"] : unable to retrieve index for requested chain ID: "<<  myChainID   <<endl;
    ErrorManager::instance.treatError();
                
}

BiopolymerClass &   BiopolymerClassContainer::updBiopolymerClass(int biopolymerClassIndex){
    if ((biopolymerClassIndex <0) || (biopolymerClassIndex >= getNumBiopolymers())) {
        ErrorManager::instance <<"["<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<"] : biopolymerClassIndex out of range:"<<biopolymerClassIndex<<endl;
        ErrorManager::instance.treatError();
    }
    map<const String,BiopolymerClass>::iterator biopolymerClassMapIterator = biopolymerClassMap.begin();
    int i = 0;
    for (i = 0; i<= biopolymerClassIndex; i++) {
        if (i == biopolymerClassIndex)
            return biopolymerClassMapIterator->second;
        else biopolymerClassMapIterator++;
    }
    ErrorManager::instance <<"["<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<"] : unable to retrieve BiopolymerClass for requested index : "<<  biopolymerClassIndex   <<endl;
    ErrorManager::instance.treatError();

}

size_t BiopolymerClassContainer::getTotalNumAtoms(){
    size_t numAtoms = 0;
    map<const String, BiopolymerClass>::iterator it;
    for(it = biopolymerClassMap.begin(); it != biopolymerClassMap.end(); it++){
        numAtoms += it->second.getNumAtoms();
    }

    return numAtoms;
}

////////////////////////////////////////////
/// applies all the mobilizer commands.  ///
////////////////////////////////////////////


void BiopolymerClassContainer::setBondMobility ( vector<BasePair> & baseOperationVector) {
    for (int q=0;q<(int)baseOperationVector.size();q++){

       if (((baseOperationVector[q]).BasePairIsTwoTransformForce).compare("mobilizer") == 0){
            cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Setting mobilizer of type "<<(baseOperationVector[q]).FirstBPEdge<<" for chain "<<(baseOperationVector[q]).FirstBPChain<<" from residue "<<(baseOperationVector[q]).FirstBPResidue.outString()<<" to residue "<<(baseOperationVector[q]).SecondBPResidue.outString()<<endl;
            MobilizerStretch dummyMobilizerStretch;
            BondMobility::Mobility myBondMobility = dummyMobilizerStretch.setBondMobility(baseOperationVector[q].FirstBPEdge ) ;
            BiopolymerClass & myBiopolymerClass ( updBiopolymerClass((baseOperationVector[q]).FirstBPChain));

            BiopolymerType::BiopolymerTypeEnum btype = biopolymerClassMap[baseOperationVector[q].FirstBPChain].biopolymerType;
            
            if (btype == BiopolymerType::RNA){
                (static_cast<RNA&>( myBiopolymerClass.myBiopolymer)).setRNABondMobility(myBondMobility,
                    SimTK::ResidueInfo::Index (myBiopolymerClass.getResidueIndex((baseOperationVector[q]).FirstBPResidue)), 
                    SimTK::ResidueInfo::Index (myBiopolymerClass.getResidueIndex((baseOperationVector[q]).SecondBPResidue))); 
            } 
            else if (btype == BiopolymerType::DNA){
                //ErrorManager::instance <<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" At this time DNA is not a supported biopolymer type. "<< baseOperationVector[q].FirstBPEdge <<". "<<endl;
                //ErrorManager::instance.treatError();
                (static_cast<DNA&>( myBiopolymerClass.myBiopolymer)).setDNABondMobility(myBondMobility,
                    SimTK::ResidueInfo::Index (myBiopolymerClass.getResidueIndex ((baseOperationVector[q]).FirstBPResidue)), 
                    SimTK::ResidueInfo::Index (myBiopolymerClass.getResidueIndex ((baseOperationVector[q]).SecondBPResidue))); 
            } 
            else if (btype == BiopolymerType::Protein) {
               myBiopolymerClass.setProteinBondMobility(
                   myBondMobility,
                   (baseOperationVector[q]).FirstBPResidue,
                   (baseOperationVector[q]).SecondBPResidue
                   );
            } 
            else {
                ErrorManager::instance << __FUNCTION__ << ": biopolymerType " << btype << " unknown" << endl;
                ErrorManager::instance.treatError();
            }
        }
    }
}

void BiopolymerClassContainer::rigidifyAllChains() {

    map<const String,BiopolymerClass>::iterator it;
    map<const String,BiopolymerClass>::iterator next;
    next =   biopolymerClassMap.begin();
    while (next != biopolymerClassMap.end())
    {
        it = next;
        //(it->second).myBiopolymer.writePdb(state,outputStream,Transform(Vec3(0)));
        BiopolymerClass & myBiopolymerClass ( (it->second));//.myBiopolymer);
        if (myBiopolymerClass.biopolymerType == BiopolymerType::RNA){
            (static_cast<RNA&>( myBiopolymerClass.myBiopolymer)).setRNABondMobility(
                    BondMobility::Rigid,
                    SimTK::ResidueInfo::Index ( myBiopolymerClass.getResidueIndex(myBiopolymerClass.getFirstResidueID())),
                    SimTK::ResidueInfo::Index ( myBiopolymerClass.getResidueIndex(myBiopolymerClass.getLastResidueID())  )
            );
        } else if (myBiopolymerClass.biopolymerType == BiopolymerType::Protein) {
            myBiopolymerClass.setProteinBondMobility(
                   BondMobility::Rigid,
                   myBiopolymerClass.getFirstResidueID(),  
                   myBiopolymerClass.getLastResidueID()     
            );
        } 
        else {
            ErrorManager::instance << __FUNCTION__ << ": biopolymerType " << myBiopolymerClass.biopolymerType << " unknown" << endl;
            ErrorManager::instance.treatError();
        }
        next++;
    }
}

Vec3 BiopolymerClassContainer::getAtomLocationInMobilizedBodyFrame(String myChainID, ResidueID myResidueID, String myAtomName){
    return updBiopolymerClass(myChainID).getAtomLocationInMobilizedBodyFrame(myResidueID,  myAtomName);
}

MobilizedBody & BiopolymerClassContainer::updAtomMobilizedBody(SimbodyMatterSubsystem & matter,String chainID, ResidueID myResidueID, String myAtomName){ 
    return updBiopolymerClass(chainID).updAtomMobilizedBody(matter,myResidueID,myAtomName);
}


void BiopolymerClassContainer::writeDefaultPdb(std::ostream& outputStream)
{
    map<const String,BiopolymerClass>::iterator it;
    map<const String,BiopolymerClass>::iterator next;
    next =   biopolymerClassMap.begin();
    while (next != biopolymerClassMap.end())
    {
        it = next;
        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<"Calling myBiopolymer.writeDefaultPdb for chain "<<it->first<<endl;
        (it->second).myBiopolymer.writeDefaultPdb(outputStream, (it->second).myBiopolymer.getTopLevelTransform()  );
        //(it->second).myBiopolymer.writeDefaultPdb(outputStream,Transform(Vec3(0)));
        next++;
    }
}


void BiopolymerClassContainer::writePdb(State & state, CompoundSystem & system, std::ostream& outputStream, int modelNumber, bool calcEnergy, int satisfiedBasePairs, int unSatisfiedBasePairs) // get the latter two from ParameterReader
{
    // static int modelNumber = 1; // increments by one at each reporting step

    // to calculate potential energy we need to be at least a Dynamics stage
    system.realize(state, Stage::Dynamics);
    outputStream << "MODEL     " << std::setw(4) << modelNumber << std::endl;
     
    map<const String,BiopolymerClass>::iterator it;
    map<const String,BiopolymerClass>::iterator next;
    next =   biopolymerClassMap.begin();
    while (next != biopolymerClassMap.end())
    {
        it = next;
        (it->second).myBiopolymer.writePdb(state,outputStream,Transform(Vec3(0)));
        //Element myAtomElement = myBiopolymer.getAtomElement(m);
        next++;
    }

    //scf added time reporting
    time_t rawtime;
    struct tm * timeinfo;
    time ( &rawtime );
    timeinfo = localtime ( &rawtime );

    outputStream << "ENDMDL" << std::endl;

    outputStream <<"REMARK seconds since January 1st 1970: "<<time ( NULL     )<<std::endl; //<<"REMARK elapsed time: "<<(clock()/CLOCKS_PER_SEC)<<std::endl;
    outputStream <<"REMARK Current time is: "<<asctime (timeinfo) <<"REMARK elapsed time: "<<(clock()/CLOCKS_PER_SEC)<<std::endl;
    outputStream.setf(ios::fixed, ios::floatfield); // set output to fixed rather than scientific format
    if (calcEnergy) 
        outputStream <<"REMARK Energy = "<<system.calcEnergy(state) <<std::endl;
    outputStream <<"REMARK Angular, Linear Momentum = "<<system.calcSystemRigidBodyMomentum(state)<<endl;
    //system.removeSystemRigidBodyMomentum(state,false);
    //outputStream<<"REMARK zeroing out rigid body momentum "<<endl;
    //outputStream <<"REMARK Angular, Linear Momentum = "<<system.calcSystemRigidBodyMomentum(state)<<endl;

    outputStream <<"REMARK ["<< __FILE__<<"] state.getNU()    = "<<state.getNU()            <<std::endl;
    outputStream <<"REMARK ["<< __FILE__<<"]Satisfied contacts : "<<satisfiedBasePairs<<endl;
    outputStream <<"REMARK ["<< __FILE__<<"]Unsatisfied contacts : "<<unSatisfiedBasePairs<<endl;

    cout<<"Just wrote structure for reporting interval # "<<modelNumber<<std::endl;
    cout <<"Satisfied contacts : "<<satisfiedBasePairs<<endl;
    // ++modelNumber;
}

bool BiopolymerClassContainer::hasChainID(String chainID){
        if (biopolymerClassMap.find(chainID) == biopolymerClassMap.end())
                {return false;}
        else
                {return true;}
}


int  BiopolymerClassContainer::validateChainID(String chainID){
    if   (! hasChainID(chainID))
    {
        ErrorManager::instance <<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": You have requested a chain ID which does not correspond to any instantiated Biopolymer: "<<chainID<<endl;
        ErrorManager::instance.treatError();
    } else return 0   ;
}

Vec3  BiopolymerClassContainer::calcAtomLocationInGroundFrame(const State & state, String chainID, ResidueID residueID, String atomName){
    validateChainID(chainID);
    return updBiopolymerClass(chainID).calcAtomLocationInGroundFrame(state,residueID,atomName);
};


//void BiopolymerClassContainer::newCalcAxes(State & state, String chain1, ResidueID residueID1, String chain2, ResidueID residueID2, LeontisWesthofBondRow myLeontisWesthofBondRow ) {
void BiopolymerClassContainer::newCalcAxes(const State& state,  LeontisWesthofBondRow myLeontisWesthofBondRow,ResidueID residueID1,ResidueID residueID2,String chain1 , String chain2,Vec3 & xAxisVector1,Vec3 & yAxisVector1, Vec3 & zAxisVector1,Vec3 & xAxisVector2,Vec3 & yAxisVector2 , Vec3 & zAxisVector2,Vec3 & glycosidicNitrogenAtom1LocationInGround,Vec3 & glycosidicNitrogenAtom2LocationInGround, Vec3 & ring1CenterLocationInGround, Vec3 & ring2CenterLocationInGround)  {

            glycosidicNitrogenAtom1LocationInGround = calcAtomLocationInGroundFrame(state, chain1, residueID1, myLeontisWesthofBondRow.residue1Atom[0]);
            glycosidicNitrogenAtom2LocationInGround = calcAtomLocationInGroundFrame(state, chain2, residueID2, myLeontisWesthofBondRow.residue2Atom[0]);

            Vec3 firstRingAtomvector1 = calcAtomLocationInGroundFrame(state,chain1,residueID1,myLeontisWesthofBondRow.residue1Atom[1])  - glycosidicNitrogenAtom1LocationInGround;
            Vec3 secondRingAtomvector1 = calcAtomLocationInGroundFrame(state,chain1,residueID1,myLeontisWesthofBondRow.residue1Atom[2])  - glycosidicNitrogenAtom1LocationInGround;
            Vec3 firstRingAtomvector2 = calcAtomLocationInGroundFrame(state,chain2,residueID2,myLeontisWesthofBondRow.residue2Atom[1])  - glycosidicNitrogenAtom2LocationInGround;
            Vec3 secondRingAtomvector2 = calcAtomLocationInGroundFrame(state,chain2,residueID2,myLeontisWesthofBondRow.residue2Atom[2])  - glycosidicNitrogenAtom2LocationInGround;
            String resName1 = myLeontisWesthofBondRow.pdbResidueName1;
            
            if ((resName1.compare("A  ") == 0) || (resName1.compare("G  ") == 0) ||
                (resName1.compare("DA ") == 0) || (resName1.compare("DG ") == 0)) 
            { //if purine

                xAxisVector1 =  -5.88327 * firstRingAtomvector1 - 6.13617 * secondRingAtomvector1;
                ring1CenterLocationInGround = (calcAtomLocationInGroundFrame(state,chain1,residueID1,String("N3")) 
                                              +calcAtomLocationInGroundFrame(state,chain1,residueID1, String("C6")))/2; 
            }
            else if ((resName1.compare("C  ") == 0) || (resName1.compare("DC ") == 0)) 
            {
                xAxisVector1 = -7.83435 * firstRingAtomvector1 -6.99265          *secondRingAtomvector1;
                ring1CenterLocationInGround = (calcAtomLocationInGroundFrame(state,chain1,residueID1,String("N1")) 
                                              +calcAtomLocationInGroundFrame(state,chain1,residueID1, String("C4")))/2; 
            }
            else if ((resName1.compare("U  ") == 0) || (resName1.compare("DT ") == 0)) 
            {
                xAxisVector1 = -7.3491 * firstRingAtomvector1 -6.47606 *secondRingAtomvector1;    
                ring1CenterLocationInGround = (calcAtomLocationInGroundFrame(state,chain1,residueID1,String("N1")) 
                                              +calcAtomLocationInGroundFrame(state,chain1,residueID1, String("C4")))/2; 

            } else 
            {
                ErrorManager::instance <<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<"Unsupported residue type : >"<<resName1<<"<"<<  endl; 
                ErrorManager::instance.treatError();
            }

            String resName2 = myLeontisWesthofBondRow.pdbResidueName2;
            if ((resName2.compare("A  ") == 0) || (resName2.compare("G  ") == 0) ||
                (resName2.compare("DA ") == 0) || (resName2.compare("DG ") == 0))
            { //if purine
                xAxisVector2 = -5.88327 * firstRingAtomvector2 -6.13617 *secondRingAtomvector2;
                ring2CenterLocationInGround = (calcAtomLocationInGroundFrame(state,chain2,residueID2,String("N3")) 
                                              +calcAtomLocationInGroundFrame(state,chain2,residueID2, String("C6")))/2; 
            }  
            else if ((resName2.compare("C  ") == 0) || (resName2.compare("DC ") == 0))
            {
                ring2CenterLocationInGround = (calcAtomLocationInGroundFrame(state,chain2,residueID2,String("N1")) 
                                              +calcAtomLocationInGroundFrame(state,chain2,residueID2, String("C4")))/2; 
                xAxisVector2 = -7.83435 * firstRingAtomvector2 -6.99265 *secondRingAtomvector2;
            }
            else if ((resName2.compare("U  ") == 0) || (resName2.compare("DT ") == 0)) 
            {
                xAxisVector2 = -7.3491  * firstRingAtomvector2 -6.47606 *secondRingAtomvector2;
                ring2CenterLocationInGround = (calcAtomLocationInGroundFrame(state,chain2,residueID2,String("N1")) 
                                              +calcAtomLocationInGroundFrame(state,chain2,residueID2, String("C4")))/2; 
            }
            else { 
                ErrorManager::instance <<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Unrecognized residue type : >"<< resName2<<"<"<< endl; 
                ErrorManager::instance.treatError();
            } // trap errors

            zAxisVector1 = (firstRingAtomvector1%secondRingAtomvector1);
            zAxisVector1 = zAxisVector1/zAxisVector1.norm();
            zAxisVector2 = (firstRingAtomvector2%secondRingAtomvector2);
            zAxisVector2 = zAxisVector2/zAxisVector2.norm();
            yAxisVector1 = zAxisVector1%xAxisVector1;
            yAxisVector1= yAxisVector1/yAxisVector1.norm();
            yAxisVector2 = zAxisVector2%xAxisVector2;
            yAxisVector2= yAxisVector2/yAxisVector2.norm();
}

void BiopolymerClassContainer::computeCorrection(LeontisWesthofClass & myLeontisWesthofClass, vector<BaseInteraction>& baseInteractionVector   ,State & state,SimbodyMatterSubsystem &  matter) {
        int i;
        String chain1, chain2;
        ResidueID residue1, residue2;
        for ( i = 0; i <(int)baseInteractionVector.size(); i ++) {
            chain1   = baseInteractionVector[i].FirstBPChain;
            chain2   = baseInteractionVector[i].SecondBPChain;
            residue1 = baseInteractionVector[i].FirstBPResidue;
            residue2 = baseInteractionVector[i].SecondBPResidue;
            if (
                (( 

                ((updBiopolymerClass((baseInteractionVector[i].FirstBPChain)).getBiopolymerType()  == BiopolymerType::RNA)) ||
                ((updBiopolymerClass((baseInteractionVector[i].FirstBPChain)).getBiopolymerType()  == BiopolymerType::DNA)) 
                ) && (
                ((updBiopolymerClass((baseInteractionVector[i].SecondBPChain)).getBiopolymerType()  == BiopolymerType::RNA)) ||
                ((updBiopolymerClass((baseInteractionVector[i].SecondBPChain)).getBiopolymerType()  == BiopolymerType::DNA)) 
                ))  
               )   
            {   
                String myResidueName1 = getPdbResidueName(chain1,residue1);
                String myResidueName2 =    
                    getPdbResidueName(chain2,residue2);
                LeontisWesthofBondRow myLeontisWesthofBondRow = myLeontisWesthofClass.getLeontisWesthofBondRow(
                    residue1, 
                    residue2, 
                    myResidueName1, 
                    baseInteractionVector[i].FirstBPEdge,   
                    myResidueName2, 
                    baseInteractionVector[i].SecondBPEdge ,
                    baseInteractionVector[i].OrientationBP ,   
                    "baseInteraction"//Vector[i].BasePairIsTwoTransformForce
                );  
                Vec3 xAxisVector1 ;
                Vec3 yAxisVector1;
                Vec3 zAxisVector1;
                Vec3 xAxisVector2;
                Vec3 yAxisVector2;
                Vec3 zAxisVector2;
                Vec3 glycosidicNitrogenAtom1LocationInGround;
                Vec3 glycosidicNitrogenAtom2LocationInGround;
                Vec3 ring1CenterLocationInGround;
                Vec3 ring2CenterLocationInGround;
                MobilizedBody body1 = updAtomMobilizedBody(matter,chain1,residue1, myLeontisWesthofBondRow.residue1Atom[0])  ;

                MobilizedBody body2 = updAtomMobilizedBody(matter,chain2,residue2, myLeontisWesthofBondRow.residue2Atom[0])  ;

              
                newCalcAxes(state,
                   myLeontisWesthofBondRow, 
                   residue1,
                   residue2,
                   chain1,
                   chain2 ,
                   xAxisVector1,yAxisVector1,zAxisVector1,xAxisVector2,yAxisVector2,zAxisVector2,
                   glycosidicNitrogenAtom1LocationInGround,
                   glycosidicNitrogenAtom2LocationInGround,
                   ring1CenterLocationInGround,ring2CenterLocationInGround);

                Rotation rotation1FromRingAtoms(Mat33(xAxisVector1,yAxisVector1,zAxisVector1));
                Rotation rotation2FromRingAtoms(Mat33(xAxisVector2,yAxisVector2,zAxisVector2));
                Rotation    myRotationCorrection1 = ~rotation1FromRingAtoms * ( matter.getMobilizedBody(body1).getBodyTransform(state)).R();
                Rotation myRotationCorrection2 = (~rotation2FromRingAtoms * ( matter.getMobilizedBody(body2).getBodyTransform(state)).R());
                Vec3 myTranslationCorrection1 = (~( matter.getMobilizedBody(body1).getBodyTransform(state)).R()*(glycosidicNitrogenAtom1LocationInGround - ( matter.getMobilizedBody(body1).getBodyTransform(state)).T()  ));
                Vec3 myTranslationCorrection2 = (~( matter.getMobilizedBody(body2).getBodyTransform(state)).R()*(glycosidicNitrogenAtom2LocationInGround - ( matter.getMobilizedBody(body2).getBodyTransform(state)).T()  ));
                (baseInteractionVector[i]).rotationCorrection1 =myRotationCorrection1;
                (baseInteractionVector[i]).rotationCorrection2 = myRotationCorrection2;
                (baseInteractionVector[i]).translationCorrection1 = myTranslationCorrection1;
                (baseInteractionVector[i]).translationCorrection2 = myTranslationCorrection2;
}
} // of for i
}

String BiopolymerClassContainer::getPdbResidueName( String chainID, ResidueID residueNumber){
    validateChainID(chainID);
    return updBiopolymerClass(chainID).getPdbResidueName(residueNumber);
}


void        BiopolymerClassContainer::setSingleBondMobility(String chainID, ResidueID residueID1,  String atomName1,ResidueID residueID2, String atomName2, String mobilityString ) {
    updBiopolymerClass(chainID).setSingleBondMobility(residueID1,  atomName1, residueID2, atomName2, mobilityString );
}


void BiopolymerClassContainer::setSingleBondMobility( vector<SingleBondMobility> mySingleBondMobilityVector) {
    

    for (int q=0;q<(int)mySingleBondMobilityVector.size();q++) {

        if ((mySingleBondMobilityVector[q].chain1).compare(mySingleBondMobilityVector[q].chain2) != 0) {ErrorManager::instance <<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<"for singleBondMobility, both atoms must be on same chain."<<endl; ErrorManager::instance.treatError();}

        setSingleBondMobility(mySingleBondMobilityVector[q].chain1, mySingleBondMobilityVector[q].residue1 , mySingleBondMobilityVector[q].atom1,  mySingleBondMobilityVector[q].residue2 , mySingleBondMobilityVector[q].atom2, mySingleBondMobilityVector[q].mobility );

    }
}

void BiopolymerClassContainer::printAllIncludedResidues (vector<IncludeAllNonBondAtomsInResidue> & includeAllNonBondAtomsInResidueVector ) {
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": Listing all residues to be included in physics zone:"<<endl; 
    for (int i = 0 ; i < (int)includeAllNonBondAtomsInResidueVector.size(); i++) {
        IncludeAllNonBondAtomsInResidue myIncludeAllNonBondAtomsInResidue = includeAllNonBondAtomsInResidueVector[i];
        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": chain = "<<myIncludeAllNonBondAtomsInResidue.getChain()<<", residue = "<<myIncludeAllNonBondAtomsInResidue.getResidue().outString()<<endl;
    }
}

#ifdef USE_OPENMM
vector< pair<const BiopolymerClass*, const ResidueID*> > BiopolymerClassContainer::getResiduesWithin(const String & chainID, const ResidueID & resID, double radius, const State & state, OpenMM::NeighborList * neighborList){
    vector<MMBAtomInfo> concatenatedAtomInfoVector = getConcatenatedAtomInfoVector(state);
    return getResiduesWithin(concatenatedAtomInfoVector, chainID, resID, radius, neighborList); // calls two below.
}

vector< pair<const BiopolymerClass*, const ResidueID*> > BiopolymerClassContainer::getResiduesWithin(const String & chainID, const ResidueID & resID, double radius, OpenMM::NeighborList * neighborList){
    vector<MMBAtomInfo> concatenatedAtomInfoVector = getConcatenatedAtomInfoVector();
    return getResiduesWithin(concatenatedAtomInfoVector, chainID, resID, radius, neighborList); // calls one below
}

vector< pair<const BiopolymerClass*, const ResidueID*> > BiopolymerClassContainer::getResiduesWithin(vector<MMBAtomInfo>& concatenatedAtomInfoVector, const String & chainID, const ResidueID & resID, double radius, OpenMM::NeighborList * neighborList){

    vector< pair<const BiopolymerClass*, const ResidueID*> > residuesWithin;
    BiopolymerClass & primaryBiopolymerClass = updBiopolymerClass(chainID);

    // We add the given residue first
    residuesWithin.push_back(make_pair(&primaryBiopolymerClass,&resID));

    // Get the neighborlist
    if(neighborList == NULL)
    {
        OpenMM::NeighborList nl = getNeighborList(concatenatedAtomInfoVector, radius);
        neighborList = &nl;
    }


    cout << "Going through the neighbors" << endl;
    // Go through the list
    for ( int j = 0 ; j < neighborList->size(); j++) 
    {
        if(j % 1000000 == 0)
            cout << "NeighborList; read " << j << " neighbors" << endl;
        unsigned int id1 = (*neighborList)[j].first;
        unsigned int id2 = (*neighborList)[j].second;

        MMBAtomInfo & atom1 = concatenatedAtomInfoVector[id1];
        MMBAtomInfo & atom2 = concatenatedAtomInfoVector[id2];

        BiopolymerClass & bpc1 = updBiopolymerClass(atom1.chain);
        BiopolymerClass & bpc2 = updBiopolymerClass(atom2.chain);

        double dist = atom1.distance(atom2);

        // // we care only about representative atoms
        // if(atom1.atomName != bpc1.getRepresentativeAtomName() || atom2.atomName != bpc2.getRepresentativeAtomName())
        //     continue;

        // if residue 1 is the given residue we add residue 2
        if(atom1.chain == chainID && atom1.residueID == resID && dist <= radius)
        {
            residuesWithin.push_back(make_pair(&(bpc2),&(atom2.residueID)));
        }
        // if residue 2 is the given residue we add residue 1
        else if(atom2.chain == chainID && atom2.residueID == resID && dist <= radius)
        {
            residuesWithin.push_back(make_pair(&(bpc1),&(atom1.residueID)));
        }
    }
    return residuesWithin;
}

OpenMM::NeighborList BiopolymerClassContainer::getNeighborList(const vector<MMBAtomInfo>& concatenatedAtomInfoVector, double radius)
{
    // Generate particle list for OpenMM
    vector<openmmVecType> particleList(concatenatedAtomInfoVector.size());
    for (int i = 0; i < concatenatedAtomInfoVector.size() ; i++) 
    {
        particleList[i] = concatenatedAtomInfoVector[i].position;
    }
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<endl;

    // Now the neighbors list
    vector<set<int> > exclusions( particleList.size() );
    OpenMM::NeighborList neighborList;
    //OpenMM::Vec3 * boxSize;
    openmmVecType boxSize = openmmVecType (10000,10000,10000);
    computeNeighborListVoxelHash(neighborList, particleList.size() , particleList, exclusions, &boxSize, false, radius, 0.0);
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<< ": NeighborList computed. Size: " << neighborList.size() << " - Radius: " << radius <<endl;
    return neighborList;
}

void BiopolymerClassContainer::setNeighborsFromList(vector<MMBAtomInfo>& concatenatedAtomInfoVector, OpenMM::NeighborList& neighborList, double radius)
{
    for (int i = 0; i < concatenatedAtomInfoVector.size() ; i++) 
    {
        concatenatedAtomInfoVector[i].clearNeighbors();
    }
    cout << "Going through the neighbors" << endl;
    // Go through the list
    for ( int j = 0 ; j < neighborList.size(); j++) 
    {
        if(j % 1000000 == 0)
            cout << "NeighborList; read " << j << " neighbors" << endl;
        unsigned int id1 = neighborList[j].first;
        unsigned int id2 = neighborList[j].second;

        MMBAtomInfo & atom1 = concatenatedAtomInfoVector[id1];
        MMBAtomInfo & atom2 = concatenatedAtomInfoVector[id2];

        double dist = atom1.distance(atom2);
        if(dist <= radius)
        {
            atom1.addNeighbor(&atom2);
            atom2.addNeighbor(&atom1); // For applyContactsWithin, ContactContainer.cpp : ContactContainer::createContactsWithin loops over a list of atoms and looks at its neighbor. since the neighborlist is on order of larger on left, smaller on right, then we have to do both forward and reverse linking. This of course doubles the number of neighbors. But I think it's more efficient than trying to do double looping.  Hopefully memory won't be an issue. 

            //cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Added neighbor atom 1 = "<<atom1.getChain()<<atom1.getResidueIndex()<<atom1.getAtomName() << " , atom 2 = "<<atom2.getChain()<<atom2.getResidueIndex()<<atom2.getAtomName()<<std::endl;
            //cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Added neighbor atom 2 = "<<atom2.getChain()<<atom2.getResidueIndex()<<atom2.getAtomName() << " , atom 1 = "<<atom1.getChain()<<atom1.getResidueIndex()<<atom1.getAtomName()<<std::endl;
	     
        }
    }    
}



//vector<ResidueID> BiopolymerClassContainer::findBiopolymerResiduesWithinRadius (const AllResiduesWithin & allResiduesWithin, const String targetChainID, const State state) {
    //vector <SingleResidue> neighboringSingleResidueVector = findBiopolymerResiduesWithinRadius(allResiduesWithin, state);
    //return neighboringSingleResidueVector;
//}

// Returns only those elements of vector<SingleResidue>  singleResidueVector that belong to the current chain.
vector<ResidueID> BiopolymerClass::filterSingleResidueVector          (const vector<SingleResidue>  singleResidueVector) {
    vector<ResidueID> myResidueIDVector; myResidueIDVector.clear();
    for (int i = 0 ; i < singleResidueVector.size(); i++){
        if (singleResidueVector[i].getChain() == getChainID()){
           myResidueIDVector.push_back(singleResidueVector[i].getResidue()); 
        }
    }
    sort (myResidueIDVector); // No reason not to return this nicely sorted. We will depend on that in ThreadingStruct.
    return myResidueIDVector;
}

// SCF : should create a generic method which takes a vector of residues, and produces another vector of residues within a distance of the first.  Then we can use that method for the applyMobilizersWithin command.

template <class type>
// here again, target type is AllResiduesWithin, but can also be e.g. MobilizerWithin
void BiopolymerClassContainer::findBiopolymerResiduesWithinRadius (const type & allResiduesWithin, const State state, vector<SingleResidue> & neighboringResidueVector) {
    // Find the max radius and add the requested residues to the physics vector
    //double maxRadius = 0.0;
    // In order to avoid duplicates, we are taking neighboringResidueVector as an argument now, and appending to it.
    //vector<SingleResidue> neighboringResidueVector; 
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" allResiduesWithin = "<< endl; allResiduesWithin.print();  //OK to here
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" to now, return vector has  " <<neighboringResidueVector.size() <<" residues"<<endl; 
    
    // Add the residue to  vector
    SingleResidue neighboringResidue;
    // Wait, is AllResiduesWithin a descendant of SingleResidue?
    neighboringResidue.setChain( allResiduesWithin.getChain());
    neighboringResidue.setResidue( allResiduesWithin.getResidue());
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<"  "<<endl;
    auto it0 = std::find(std::begin(neighboringResidueVector), std::end(neighboringResidueVector), neighboringResidue);            
    // If neighboringResidue is not already in neighborVector, add it:
    if(it0 == neighboringResidueVector.end() && updBiopolymerClass(neighboringResidue.getChain()).getActivePhysics()) {
        neighboringResidueVector.push_back(neighboringResidue);
    }
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" to now, return vector has  " <<neighboringResidueVector.size() <<" residues"<<endl; 
    
    cout << __FILE__ <<":"<<__LINE__<<": allResiduesWithin.getRadius() " <<allResiduesWithin.getRadius() << endl; // OK to here
    // allResiduesWithinVector was empty
    // No residues to add so no need to compute the neighbor list
    if(allResiduesWithin.getRadius() <= 1E-14)
    {
        if (neighboringResidueVector.size() > 1) {
            //ErrorManager::instance 
            cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" WARNING: neighboringResidueVector should have size <= 1, instead size = "<<neighboringResidueVector.size() <<" If this corresponds to the accumulated number of residues given in allResiduesWithin over several calls to this function, these are self-matches, this is probably fine. "<<endl;
            //ErrorManager::instance.treatError();
        }
        //return neighboringResidueVector;
    }

    // Generate particle list for OpenMM
    cout << __FILE__ <<":"<<__LINE__<<endl;
    vector<MMBAtomInfo> concatenatedAtomInfoVector = getConcatenatedAtomInfoVector(); // Why was this not necessary before?
    vector<openmmVecType> particleList(concatenatedAtomInfoVector.size());
    for (int i = 0; i < concatenatedAtomInfoVector.size() ; i++) 
    {
        #ifdef NEIGHBORLISTCAONLY
        // This alternate compilation means only CA atoms will be taken into account when figuring out the physics zone from a given radius. 
        if (concatenatedAtomInfoVector[i].getAtomName( ) == updBiopolymerClass(concatenatedAtomInfoVector[i].getChain()).getRepresentativeAtomName()) {
        #endif
            
            particleList[i] = concatenatedAtomInfoVector[i].position;
        #ifdef NEIGHBORLISTCAONLY
        } else {
            particleList.erase(particleList.begin() + i); 
            concatenatedAtomInfoVector.erase(concatenatedAtomInfoVector.begin() + i); 
            i--;  
            cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" i,  concatenatedAtomInfoVector.size() = "<< i <<", " << concatenatedAtomInfoVector.size() <<endl;
        }
        #endif
    }
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<endl;

    // Now the neighbors list
    vector<set<int> > exclusions( particleList.size() );
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<endl;
    OpenMM::NeighborList neighborList;
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<endl;
    openmmVecType boxSize (10000,10000,10000);
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<endl;
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" neighborList size is : "<<neighborList.size()<<endl;
    computeNeighborListVoxelHash(neighborList, particleList.size() , particleList, exclusions, &boxSize, false, allResiduesWithin.getRadius(), 0.0);
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" neighborList size is : "<<neighborList.size()<<endl;

    // Go through the list
    for ( int j = 0 ; j < neighborList.size(); j++) 
    {
        unsigned int id1 = neighborList[j].first;
        unsigned int id2 = neighborList[j].second;

        String name1 = concatenatedAtomInfoVector[id1].atomName;
        String name2 = concatenatedAtomInfoVector[id2].atomName;

        SingleResidue incl1; incl1.setChain(concatenatedAtomInfoVector[id1].chain); incl1.setResidue(concatenatedAtomInfoVector[id1].residueID);
        SingleResidue incl2; incl2.setChain(concatenatedAtomInfoVector[id2].chain); incl2.setResidue(concatenatedAtomInfoVector[id2].residueID); 

        double dist = concatenatedAtomInfoVector[id1].distance(concatenatedAtomInfoVector[id2]);

        if ((incl1 == allResiduesWithin) && updBiopolymerClass(incl1.getChain()).getActivePhysics())
        {
            if(dist <= allResiduesWithin.getRadius() &&
                //continue;
               (incl2 != SingleResidue(allResiduesWithin)) && 
               updBiopolymerClass(incl2.getChain()).getActivePhysics())
            {

                auto it = std::find(std::begin(neighboringResidueVector), std::end(neighboringResidueVector), incl2);            
		//vector<type>::iterator it = find (neighboringResidueVector.begin(), neighboringResidueVector.end(), incl2); 
		if(it == neighboringResidueVector.end() && updBiopolymerClass(incl2.getChain()).getActivePhysics()) {
                    #ifdef NEIGHBORLISTCAONLY
                    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Confirm that these are representative atom names : "<<endl;// name1 << " , "<<name2<<endl;
                    concatenatedAtomInfoVector[id1].print(); 
                    concatenatedAtomInfoVector[id2].print(); 
                    #endif
		    neighboringResidueVector.push_back(incl2);
                    incl2.printStretch();
		}
                //neighboringResidueVector.push_back(incl2);
            }
        }
        else if ((incl2 == allResiduesWithin) && updBiopolymerClass(incl2.getChain()).getActivePhysics())
        {
            if ((dist <= allResiduesWithin.getRadius()) &&
            //    continue;
               (incl1 != SingleResidue(allResiduesWithin)) && 
               updBiopolymerClass(incl1.getChain()).getActivePhysics())
            {
		vector<SingleResidue>::iterator it = find (neighboringResidueVector.begin(), neighboringResidueVector.end(), incl1); 
		if(it == neighboringResidueVector.end() && updBiopolymerClass(incl1.getChain()).getActivePhysics()) {
                    #ifdef NEIGHBORLISTCAONLY
                    //cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Confirm that these are representative atom names : "<< name1 << " , "<<name2<<endl;
                    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Confirm that these are representative atom names : "<<endl;// name1 << " , "<<name2<<endl;
                    concatenatedAtomInfoVector[id1].print(); 
                    concatenatedAtomInfoVector[id2].print(); 
                    #endif
		    neighboringResidueVector.push_back(incl1);
                    incl1.printStretch();
		}
            }
        }
    }
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" returning vector with " <<neighboringResidueVector.size() <<" residues"<<endl; 
    //return neighboringResidueVector;
}

//template <type2>
// target type is AllResiduesWithin. But can also be e.g. MobilizerWithin
template <class type2> vector<SingleResidue> BiopolymerClassContainer::findBiopolymerResiduesWithinRadius (const vector<type2> & allResiduesWithinVector, const State state) {
    //vector<SingleResidue> neighborVector; neighborVector.clear();
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" allResiduesWithinVector.size() = "<<allResiduesWithinVector.size() <<std::endl;
    vector<SingleResidue> neighborVector ; neighborVector.clear();
    for (int i = 0 ; i < allResiduesWithinVector.size() ; i++) {
        // This will progressively extend neighborVector, avoiding duplicates:  
        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" i = "<<i<<" "<<endl; // OK to now
        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" allResiduesWithinVector["<<i<<"] = "<<endl; allResiduesWithinVector[i].print();   // This is returning correct radius 
        findBiopolymerResiduesWithinRadius<type2>(allResiduesWithinVector[i],state, neighborVector );
        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<"  "<<endl;
        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" neighborVector now has "<<neighborVector.size()<<" elements."<<endl;
    }
   cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" neighborVector now has "<<neighborVector.size()<<" elements."<<endl;
   return neighborVector;
}
vector<SingleResidue> BiopolymerClassContainer::findBiopolymerResiduesWithinRadius (const vector<MobilizerWithin> & allResiduesWithinVector, const State state) {
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" allResiduesWithinVector.size() = "<<allResiduesWithinVector.size()  <<endl;
    return findBiopolymerResiduesWithinRadius<MobilizerWithin>(allResiduesWithinVector,state);
}

// Keep the code below, it is more efficient.  Doesn't rebuild the hash each time, and also doesn't make redundant additions to the return vector.
/*
vector<SingleResidue> BiopolymerClassContainer::findBiopolymerResiduesWithinRadius (const vector<AllResiduesWithin> & allResiduesWithinVector, const State state) {
    // Find the max radius and add the requested residues to the physics vector
    double maxRadius = 0.0;
    vector<AllResiduesWithin>::const_iterator itARW;
    vector<SingleResidue>   neigboringResidueVector ;
    neigboringResidueVector.clear();
    for(itARW = allResiduesWithinVector.begin(); itARW != allResiduesWithinVector.end(); itARW++)
    {
        if(itARW->getRadius() > maxRadius)
            maxRadius = itARW->getRadius();
        if(updBiopolymerClass(itARW->getChain()).getActivePhysics() == false)
            continue;

        // Add the query residue to neighboring residue vector
        SingleResidue neighboringResidue;
        neighboringResidue.setChain( itARW->getChain());
        neighboringResidue.setResidue( itARW->getResidue());
        vector<SingleResidue>::iterator it = find (neigboringResidueVector.begin(), neigboringResidueVector.end(), neighboringResidue); 
        if(it == neigboringResidueVector.end()) {
            neigboringResidueVector.push_back(neighboringResidue);
            cout << __FILE__ <<":"<<__LINE__<<endl;
            neighboringResidue.printStretch();
        }
    }

    cout << __FILE__ <<":"<<__LINE__<<": maxRadius " << maxRadius << endl;
    // allResiduesWithinVector was empty
    // No residues to add so no need to compute the neighbor list
    if(maxRadius <= 1E-14)
    {
        if (neigboringResidueVector.size() > 0) {
            ErrorManager::instance <<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" : Something odd happened.  neigboringResidueVector should have size 0, instead size = "<<neigboringResidueVector.size() <<endl;
            ErrorManager::instance.treatError();
        }
        return neigboringResidueVector;
    }

    // Generate particle list for OpenMM
    cout << __FILE__ <<":"<<__LINE__<<endl;
    vector<MMBAtomInfo> concatenatedAtomInfoVector = getConcatenatedAtomInfoVector(); // Why was this not necessary before?
    vector<openmmVecType> particleList(concatenatedAtomInfoVector.size());
    for (int i = 0; i < concatenatedAtomInfoVector.size() ; i++) 
    {
        particleList[i] = concatenatedAtomInfoVector[i].position;
    }
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<endl;

    // Now the neighbors list
    vector<set<int> > exclusions( particleList.size() );
    OpenMM::NeighborList neighborList;
    openmmVecType boxSize (10000,10000,10000);
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" neighborList size is : "<<neighborList.size()<<endl;
    computeNeighborListVoxelHash(neighborList, particleList.size() , particleList, exclusions, &boxSize, false, maxRadius, 0.0);
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" neighborList size is : "<<neighborList.size()<<endl;

    // Go through the list
    for ( int j = 0 ; j < neighborList.size(); j++) 
    {
        unsigned int id1 = neighborList[j].first;
        unsigned int id2 = neighborList[j].second;

        String name1 = concatenatedAtomInfoVector[id1].atomName;
        String name2 = concatenatedAtomInfoVector[id2].atomName;

        SingleResidue incl1; incl1.setChain(concatenatedAtomInfoVector[id1].chain); incl1.setResidue(concatenatedAtomInfoVector[id1].residueID);
        SingleResidue incl2; incl2.setChain(concatenatedAtomInfoVector[id2].chain); incl2.setResidue(concatenatedAtomInfoVector[id2].residueID); 

        double dist = concatenatedAtomInfoVector[id1].distance(concatenatedAtomInfoVector[id2]);


        // if neighbor 1 is requested, add neighbor 2 if active and not already added.
        itARW = find(allResiduesWithinVector.begin(), allResiduesWithinVector.end(), incl1);
        if(itARW!= allResiduesWithinVector.end() && updBiopolymerClass(incl1.getChain()).getActivePhysics())
        {
            if(dist > itARW->getRadius())
                continue;
            vector<SingleResidue>::iterator it = find (neigboringResidueVector.begin(), neigboringResidueVector.end(), incl2); 
            if(it == neigboringResidueVector.end() && updBiopolymerClass(incl2.getChain()).getActivePhysics()) {
                neigboringResidueVector.push_back(incl2);
            }
        }

        // if neighbor 2 is requested, add neighbor 1 if active and not already added.
        itARW = find(allResiduesWithinVector.begin(), allResiduesWithinVector.end(), incl2);
        if(itARW != allResiduesWithinVector.end() && updBiopolymerClass(incl1.getChain()).getActivePhysics())
        {
            if(dist > itARW->getRadius())
                continue;
            vector<SingleResidue>::iterator it = find (neigboringResidueVector.begin(), neigboringResidueVector.end(), incl1); 
            if(it == neigboringResidueVector.end() && updBiopolymerClass(incl1.getChain()).getActivePhysics()) {
                neigboringResidueVector.push_back(incl1);


            }
        }
    }
    return neigboringResidueVector;
} */

void BiopolymerClassContainer::includeAllResiduesWithin (const vector<AllResiduesWithin> & includeAllResiduesWithinVector, 
                            vector<IncludeAllNonBondAtomsInResidue> & includeAllNonBondAtomsInResidueVector, 
                            const State state) {
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" includeAllNonBondAtomsInResidueVector  size is now : "<<includeAllNonBondAtomsInResidueVector.size()<<endl;
    vector<IncludeAllNonBondAtomsInResidue> tempIncludeAllNonBondAtomsInResidueVector = findBiopolymerResiduesWithinRadius(includeAllResiduesWithinVector, state);
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" tempIncludeAllNonBondAtomsInResidueVector  size is now : "<<tempIncludeAllNonBondAtomsInResidueVector.size()<<endl; // This is returning 0 .. check/
    includeAllNonBondAtomsInResidueVector.reserve(includeAllNonBondAtomsInResidueVector.size() + tempIncludeAllNonBondAtomsInResidueVector.size());
    includeAllNonBondAtomsInResidueVector.insert(includeAllNonBondAtomsInResidueVector.end(), tempIncludeAllNonBondAtomsInResidueVector.begin(), tempIncludeAllNonBondAtomsInResidueVector.end());
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" includeAllNonBondAtomsInResidueVector  size is now : "<<includeAllNonBondAtomsInResidueVector.size()<<endl;
    
}
#endif

void BiopolymerClassContainer::includeAllNonBondAtomsInResidues(vector<IncludeAllNonBondAtomsInResidue>  myIncludeAllNonBondAtomsInResidueVector, State & state, DuMMForceFieldSubsystem & dumm) {
    for (int i = 0; i < (int)myIncludeAllNonBondAtomsInResidueVector.size(); i++){
        // Skip residues in non active chains
        if(updBiopolymerClass(myIncludeAllNonBondAtomsInResidueVector[i].getChain()).getActivePhysics()==false)
            continue;
        updBiopolymerClass(myIncludeAllNonBondAtomsInResidueVector[i].getChain()).includeAllNonBondAtomsInResidue(myIncludeAllNonBondAtomsInResidueVector[i].getResidue(),state,dumm);
    }
}

/*void BiopolymerClassContainer::includeAllResiduesWithin (const vector<AllResiduesWithin> & includeAllResiduesWithinVector, 
                            vector<IncludeAllNonBondAtomsInResidue> & includeAllNonBondAtomsInResidueVector, 
                            const State state) {
    // Find the max radius and add the requested residues to the physics vector
    double maxRadius = 0.0;
    vector<AllResiduesWithin>::const_iterator itARW;
    for(itARW = includeAllResiduesWithinVector.begin(); itARW != includeAllResiduesWithinVector.end(); itARW++)
    {
        if(itARW->getRadius() > maxRadius)
            maxRadius = itARW->getRadius();

        if(updBiopolymerClass(itARW->getChain()).getActivePhysics() == false)
            continue;
        // Add the residue to physics vector
        IncludeAllNonBondAtomsInResidue myIncludeAllNonBondAtomsInResidue;
        myIncludeAllNonBondAtomsInResidue.setChain( itARW->getChain());
        myIncludeAllNonBondAtomsInResidue.setResidue( itARW->getResidue());
        vector<IncludeAllNonBondAtomsInResidue>::iterator it = find (includeAllNonBondAtomsInResidueVector.begin(), includeAllNonBondAtomsInResidueVector.end(), myIncludeAllNonBondAtomsInResidue); 
        if(it == includeAllNonBondAtomsInResidueVector.end()) {
            includeAllNonBondAtomsInResidueVector.push_back(myIncludeAllNonBondAtomsInResidue);
            cout << __FILE__ <<":"<<__LINE__<<endl;
            myIncludeAllNonBondAtomsInResidue.printStretch();
        }
    }

    cout << __FILE__ <<":"<<__LINE__<<": maxRadius " << maxRadius << endl;
    // includeAllResiduesWithinVector was empty
    // No residues to add so no need to compute the neighbor list
    if(maxRadius <= 0)
    {
        return;
    }

    // Generate particle list for OpenMM
    cout << __FILE__ <<":"<<__LINE__<<endl;
    //printAtomInfoVector(); // NOT fine at this point
    vector<MMBAtomInfo> concatenatedAtomInfoVector = getConcatenatedAtomInfoVector(); // Why was this not necessary before?
    vector<openmmVecType> particleList(concatenatedAtomInfoVector.size());
    for (int i = 0; i < concatenatedAtomInfoVector.size() ; i++) 
    {
        particleList[i] = concatenatedAtomInfoVector[i].position;
    }
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<endl;

    // Now the neighbors list
    vector<set<int> > exclusions( particleList.size() );
    OpenMM::NeighborList neighborList;
    openmmVecType boxSize (10000,10000,10000);
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" neighborList size is : "<<neighborList.size()<<endl;
    computeNeighborListVoxelHash(neighborList, particleList.size() , particleList, exclusions, &boxSize, false, maxRadius, 0.0);
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" neighborList size is : "<<neighborList.size()<<endl;

    // Go through the list
    for ( int j = 0 ; j < neighborList.size(); j++) 
    {
        unsigned int id1 = neighborList[j].first;
        unsigned int id2 = neighborList[j].second;

        String name1 = concatenatedAtomInfoVector[id1].atomName;
        String name2 = concatenatedAtomInfoVector[id2].atomName;

        IncludeAllNonBondAtomsInResidue incl1; incl1.setChain(concatenatedAtomInfoVector[id1].chain); incl1.setResidue(concatenatedAtomInfoVector[id1].residueID);
        IncludeAllNonBondAtomsInResidue incl2; incl2.setChain(concatenatedAtomInfoVector[id2].chain); incl2.setResidue(concatenatedAtomInfoVector[id2].residueID); 

        double dist = concatenatedAtomInfoVector[id1].distance(concatenatedAtomInfoVector[id2]);


        // if neighbor 1 is requested, add neighbor 2 if active and not already added.
        itARW = find(includeAllResiduesWithinVector.begin(), includeAllResiduesWithinVector.end(), incl1);
        if(itARW!= includeAllResiduesWithinVector.end() && updBiopolymerClass(incl1.getChain()).getActivePhysics())
        {
            if(dist > itARW->getRadius())
                continue;
            vector<IncludeAllNonBondAtomsInResidue>::iterator it = find (includeAllNonBondAtomsInResidueVector.begin(), includeAllNonBondAtomsInResidueVector.end(), incl2); 
            if(it == includeAllNonBondAtomsInResidueVector.end() && updBiopolymerClass(incl2.getChain()).getActivePhysics()) {
                includeAllNonBondAtomsInResidueVector.push_back(incl2);
            }
        }

        // if neighbor 2 is requested, add neighbor 1 if active and not already added.
        itARW = find(includeAllResiduesWithinVector.begin(), includeAllResiduesWithinVector.end(), incl2);
        if(itARW != includeAllResiduesWithinVector.end() && updBiopolymerClass(incl1.getChain()).getActivePhysics())
        {
            if(dist > itARW->getRadius())
                continue;
            vector<IncludeAllNonBondAtomsInResidue>::iterator it = find (includeAllNonBondAtomsInResidueVector.begin(), includeAllNonBondAtomsInResidueVector.end(), incl1); 
            if(it == includeAllNonBondAtomsInResidueVector.end() && updBiopolymerClass(incl1.getChain()).getActivePhysics()) {
                includeAllNonBondAtomsInResidueVector.push_back(incl1);


            }
        }
    }
}*/

void BiopolymerClassContainer::includeNonBondAtoms(  vector<IncludeNonBondAtomInBiopolymerStruct> includeNonBondAtomInBiopolymerVector,  State & state, DuMMForceFieldSubsystem & dumm) {
    for (int i = 0 ; i < (int)includeNonBondAtomInBiopolymerVector.size(); i++) {
        includeNonBondAtom(includeNonBondAtomInBiopolymerVector[i].chain,  includeNonBondAtomInBiopolymerVector[i].residue, includeNonBondAtomInBiopolymerVector[i].atomName, state,dumm);
    }
}

void BiopolymerClassContainer::includeNonBondAtom(String chain , ResidueID residue, String atomName ,  State & state, DuMMForceFieldSubsystem & dumm) {
    // Skip atoms in non active chains
    if(updBiopolymerClass(chain).getActivePhysics()==false)
        return;
    updBiopolymerClass(chain).includeNonBondAtom(residue, atomName, state, dumm);
}

void BiopolymerClassContainer::waterDropletAboutResidues (const vector <WaterDropletAboutResidueStruct> waterDropletAboutResidueVector,    WaterDropletContainer & waterDropletContainer  )     {
        for (int i = 0; i < (int)waterDropletAboutResidueVector.size(); i++) {
                 BiopolymerClass  primaryBiopolymerClass = updBiopolymerClass(waterDropletAboutResidueVector[i]. biopolymerChainID );
                 cout<<__FILE__<<" : "<<__LINE__<<" : "<<primaryBiopolymerClass.getRepresentativeAtomName()<<endl;
                 Vec3 myLocation = (primaryBiopolymerClass.calcDefaultAtomLocationInGroundFrame(waterDropletAboutResidueVector[i].residue, primaryBiopolymerClass.getRepresentativeAtomName()))*(1.0); // used to convert to Å, now using nm

                 /*
                 WaterDroplet myWaterDroplet;
                 myWaterDroplet.chainID   =waterDropletAboutResidueVector[i].waterDropletChainID ;

                 myWaterDroplet.center = myLocation;
                 myWaterDroplet.setRadius ( waterDropletAboutResidueVector[i].radius);
                 myWaterDroplet.tetherStrength =  waterDropletAboutResidueVector[i].tetherStrength;
                 waterDropletContainer.add(myWaterDroplet); 
                 */
                 // WaterDropletContainer is now an incomplete type. So this no longer works. Leaving broken for now. Should fix before using. Will exit. Do not try to #include WaterDroplet .. this won't work ever since we are #include Threading.h  
                 cout<<__FILE__<<" : "<<__LINE__<<" : "<<endl;
                 exit(1);
   } // of for i
}

void BiopolymerClassContainer::physicsZone(vector<AllResiduesWithin> & myIncludeAllResiduesWithinVector , double radius, SimbodyMatterSubsystem & matter,State & state) {
    map<const String,BiopolymerClass>::iterator biopolymerClassMapIterator = biopolymerClassMap.begin();
    for(biopolymerClassMapIterator = biopolymerClassMap.begin(); biopolymerClassMapIterator != biopolymerClassMap.end(); biopolymerClassMapIterator++) {
        if(biopolymerClassMapIterator->second.getActivePhysics())
            (biopolymerClassMapIterator->second).physicsZone(myIncludeAllResiduesWithinVector, radius, matter, state);
    }  
}

void BiopolymerClassContainer::multiplySmallGroupInertia( double multiplier, CompoundSystem & system,SimbodyMatterSubsystem & matter,State & state) {
    map<const String,BiopolymerClass>::iterator biopolymerClassMapIterator = biopolymerClassMap.begin();
    for(biopolymerClassMapIterator = biopolymerClassMap.begin(); biopolymerClassMapIterator != biopolymerClassMap.end(); biopolymerClassMapIterator++) {
        (biopolymerClassMapIterator->second). multiplySmallGroupInertia(multiplier,system,matter,state);
    }  

}

// parameter endCaps, when set to True, tells us that myBiopolymer has end caps which should be ignored when extracting the sequence.

String BiopolymerClassContainer::extractSequenceFromBiopolymer(const Biopolymer & myBiopolymer, bool endCaps = 0      ){
    stringstream mySequence;
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<endl;
    for (int i = (0 + endCaps) ; i < (myBiopolymer.getNumResidues() - endCaps); i++) {
        // cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" "<<i<<" "<<myBiopolymer.getResidue      (ResidueInfo::Index(i)).getOneLetterCode()<<endl;
        mySequence<<myBiopolymer.getResidue      (ResidueInfo::Index(i)).getOneLetterCode(); 
    }
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Extracted sequence: "<< mySequence.str() << endl;
    return mySequence.str();
};

#ifdef USE_OPENMM
void BiopolymerClassContainer::initializeAtomInfoVectors(SimbodyMatterSubsystem& matter ) {
    map<const String,BiopolymerClass>::iterator biopolymerClassMapIterator = biopolymerClassMap.begin();
    for(biopolymerClassMapIterator = biopolymerClassMap.begin(); biopolymerClassMapIterator != biopolymerClassMap.end(); biopolymerClassMapIterator++) {
        (biopolymerClassMapIterator->second).initializeAtomInfoVector(matter,  atomicPropertyOverrideVector);
    }  
};


void BiopolymerClassContainer::initializeAtomInfoVectors(SimbodyMatterSubsystem& matter, DuMMForceFieldSubsystem & dumm) {
    map<const String,BiopolymerClass>::iterator biopolymerClassMapIterator = biopolymerClassMap.begin();
    for(biopolymerClassMapIterator = biopolymerClassMap.begin(); biopolymerClassMapIterator != biopolymerClassMap.end(); biopolymerClassMapIterator++) {
        (biopolymerClassMapIterator->second).initializeAtomInfoVector(matter, dumm,   atomicPropertyOverrideVector);
    }  
};
#endif

const bool isRNAtest(const Biopolymer & inputBiopolymer){
    for (int i = 0; i < inputBiopolymer.getNumResidues(); i++) {
        const ResidueInfo myResidueInfo = inputBiopolymer.getResidue(ResidueInfo::Index(i));
        const char myOneLetterCode = myResidueInfo.getOneLetterCode();
        if (! letterIsRNA(String(myOneLetterCode))) {
            return false;    
        }
        if (! inputBiopolymer.hasAtom("0/O2*")) {
            //cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" No O2* atom found on first residue.  This is not RNA! "<<endl;
            return false;
        }
    }
    return true;

}


const bool BiopolymerClassContainer::isRNA(const Biopolymer & inputBiopolymer)  {
   /*
    for (int i = 0; i < inputBiopolymer.getNumResidues(); i++) {
        const ResidueInfo myResidueInfo = inputBiopolymer.getResidue(ResidueInfo::Index(i));
        const char myOneLetterCode = myResidueInfo.getOneLetterCode();
        if (! letterIsRNA(String(myOneLetterCode))) {
            return false;    
        }
        if (! inputBiopolymer.hasAtom("0/O2*")) {
            cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" No O2* atom found on first residue.  This is not RNA! "<<endl;
            return false;
        }
    }
    return true;
*/
    return isRNAtest(inputBiopolymer);
};

const bool BiopolymerClass::isRNA()  {
    /*
    for (int i = 0; i < this->updBiopolymer().getNumResidues(); i++) {
        const ResidueInfo myResidueInfo = this->updBiopolymer().getResidue(ResidueInfo::Index(i));
        const char myOneLetterCode = myResidueInfo.getOneLetterCode();
        if (! letterIsRNA(String(myOneLetterCode))) {
            return false;    
        }
        if (! this->updBiopolymer().hasAtom("0/O2*")) {
            cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" No O2* atom found on first residue.  This is not RNA! "<<endl;
            return false;
        }
    }
    return true;
 */
    return isRNAtest(this->updBiopolymer());
};

const bool isDNAtest(const Biopolymer & inputBiopolymer)  {
    for (int i = 0; i < inputBiopolymer.getNumResidues(); i++) {
        const ResidueInfo myResidueInfo = inputBiopolymer.getResidue(ResidueInfo::Index(i));
        const char myOneLetterCode = myResidueInfo.getOneLetterCode();
        if (! letterIsDNA(String(myOneLetterCode))) {
            return false;    
        }
        if ( inputBiopolymer.hasAtom("0/O2*")) {
            cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" O2* atom found on first residue.  This is not DNA! "<<endl;
            return false;
        }
    }
    return true;
};

const bool BiopolymerClass::isDNA()  {
    return isDNAtest(this->updBiopolymer());
    /*
    for (int i = 0; i < this->updBiopolymer().getNumResidues(); i++) {
        const ResidueInfo myResidueInfo = this->updBiopolymer().getResidue(ResidueInfo::Index(i));
        const char myOneLetterCode = myResidueInfo.getOneLetterCode();
        if (! letterIsDNA(String(myOneLetterCode))) {
            return false;    
        }
        if ( this->updBiopolymer().hasAtom("0/O2*")) {
            cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" O2* atom found on first residue.  This is not DNA! "<<endl;
            return false;
        }
    }
    return true;*/
};
const bool BiopolymerClassContainer::isDNA(const Biopolymer & inputBiopolymer)  {
    return isDNAtest(inputBiopolymer);
    /*
    for (int i = 0; i < inputBiopolymer.getNumResidues(); i++) {
        const ResidueInfo myResidueInfo = inputBiopolymer.getResidue(ResidueInfo::Index(i));
        const char myOneLetterCode = myResidueInfo.getOneLetterCode();
        if (! letterIsDNA(String(myOneLetterCode))) {
            return false;    
        }
        if ( inputBiopolymer.hasAtom("0/O2*")) {
            cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" O2* atom found on first residue.  This is not DNA! "<<endl;
            return false;
        }
    }
    return true;
    */
};

const bool BiopolymerClassContainer::isProtein(const Biopolymer & inputBiopolymer, bool endCaps = true)  {
    
    for (int i = (0+ endCaps) ; i < (inputBiopolymer.getNumResidues() - endCaps ); i++) {
        const ResidueInfo myResidueInfo = inputBiopolymer.getResidue(ResidueInfo::Index(i));
        const char myOneLetterCode = myResidueInfo.getOneLetterCode();
        if (! letterIsProtein(String(myOneLetterCode))) {
            return false;    
        }

        if (! inputBiopolymer.hasAtom("0/CA")) {
            return false;
        }

    }
    return true;
};

void BiopolymerClassContainer::loadSequencesFromPdb(const String inPDBFileName,const bool proteinCapping, const String & chainsPrefix, const bool tempRenumberPdbResidues, bool useNACappingHydroxyls ){
    //std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" >"<< deletedResidueVector.size() <<"<"<<std::endl;
    std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<std::endl;
    std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" About to load sequences from file : "<<inPDBFileName<<std::endl;
    //CheckFile myCheckFile(inPDBFileName);
    //myCheckFile.validateExists();
    //myCheckFile.validateNonZeroSize(); 
    struct stat st;
    // Just querying the members of st is not a good idea. First, check to make sure stat succeeded at all: 
    if ((stat(inPDBFileName.c_str(), &st)) == -1){
        ErrorManager::instance <<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__ <<" ERROR: Tried to determine status of file >"<< inPDBFileName << "< but stat returned a failure code. Perhaps the file does not exist, or path permissions are not correct."<<std::endl;
        ErrorManager::instance.treatError();
    }

    stat(inPDBFileName.c_str(), &st);

    std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__ <<" About to check that "<<inPDBFileName<<" has nonzero size.."<<std::endl;
    if ( st.st_size == 0){
        ErrorManager::instance <<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__ <<" ERROR: Apparently "<<inPDBFileName<<" has size "<<st.st_size <<" . Dying now."<<std::endl;
        ErrorManager::instance.treatError();
    } else {
        std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Apparently "<<inPDBFileName<<" has size "<<st.st_size <<" . This seems OK."<<std::endl;
    }
    PDBReader myPDBReader ( inPDBFileName ); //, deletedResidueVector);
    std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<std::endl;
    CompoundSystem system;
    std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<std::endl;
    SimbodyMatterSubsystem  matter(system);
    std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<std::endl;
    GeneralForceSubsystem forces(system);
    DuMMForceFieldSubsystem dumm(system);
    dumm.loadAmber99Parameters();
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" About to issue myPDBReader.createCompounds( system)"<<endl;
    myPDBReader.createCompounds( system, chainsPrefix );
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Done with myPDBReader.createCompounds( system)"<<endl;
    cout<<__FILE__<<":"<<__LINE__<<" system.getNumCompounds() = "<<system.getNumCompounds() <<endl;
    
    //================================================ Use PDB reader or CIF reader depending on the extension.
    PdbStructure myPdbStructure;
    if ( inPDBFileName.substr ( inPDBFileName.length() - 4, inPDBFileName.length() - 1) == ".pdb" )
    {
        //============================================ No problem, continue as usual
        std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Filename " << inPDBFileName << " suggests PDB file. Using the PDB file reader ... " << std::endl;
        ifstream pdbfile                              ( inPDBFileName.c_str()) ;
        myPdbStructure                                = PdbStructure ( pdbfile, chainsPrefix );
        pdbfile.close                                 ( );
    }
    else
    {
        //============================================ This should be a CIF file, read it using MMDB
        std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Caching the PdbStructure from CIF file " << inPDBFileName << std::endl;
        myPdbStructure                                = PdbStructure ( inPDBFileName );
    }
    
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<endl;
    pdbStructureMap.insert(pair<String, PdbStructure>(inPDBFileName, myPdbStructure) );
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" myPdbStructure.getNumModels() "<<myPdbStructure.getNumModels()<<endl;
    int myNumChains =  myPdbStructure.getModel(Pdb::ModelIndex(0)).getNumChains();
    // PdbStructure can sometimes come up with a higher chain count. Maybe it puts in some HETATOM's or HOH's as extra chains. So we will use this one which we get more from CompoundSystem:
    int myNumChainsFromSystem =  system.getNumCompounds() /  myPdbStructure.getNumModels() ;
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" myNumChains "<<myNumChains<<endl;
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<endl;

    cout<<__FILE__<<":"<<__LINE__<<" system.getNumCompounds() = "<<system.getNumCompounds() <<endl;
    // let's check which chains we have already:    
    cout<<__FILE__<<":"<<__LINE__<<" This BiopolymerClassContainer already has getNumBiopolymers() = "<< getNumBiopolymers() <<endl;

    // system.getNumCompounds() returns the number of chains * number of models!  This is too many chains. We only want the chain in model 0 by arbitrary convention.
    //for (SimTK::CompoundSystem::CompoundIndex c(0); c < system.getNumCompounds(); ++c) 
    // We will instead use myNumChains (myPdbStructure.getModel(Pdb::ModelIndex(0)).getNumChains()) which is only the number of chains in model 0
        
    // for (SimTK::CompoundSystem::CompoundIndex c(0); c < myNumChains; ++c) // This way got us some extra chains for some reason.
    for (SimTK::CompoundSystem::CompoundIndex c(0); c < myNumChainsFromSystem; ++c) 
    {
        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Processing chain >"<<system.getCompound(c).getPdbChainId()<<"< ."<<endl;
        if (Molecule::isInstanceOf(system.getCompound(c) )) 
        {
            const Molecule   & myMolecule   = Molecule::downcast(system.getCompound(c));                
            if (Biopolymer::isInstanceOf(myMolecule ))
            {
                const Biopolymer & myBiopolymer = Biopolymer::downcast(myMolecule);

                stringstream myChainId;
                myChainId << myBiopolymer.getPdbChainId();
                String myChainIdString = myChainId.str();
                cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Chain ID :"<<myChainId.str()<<endl;
                //printBiopolymerSequenceInfo(myBiopolymer);

                if(hasChainID(myChainIdString))
                {
                    ErrorManager::instance <<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" The chain Id "<< myChainIdString << " found in PDB file " << inPDBFileName << " is already assigned to a chain in MMB. We suggest you to rename or remove the chain in the pdb file or you can use the 'deleteChain " << myChainIdString << "' command before loading the sequences in this file." << endl;
                    ErrorManager::instance.treatError();
                }

                if (isRNA(myBiopolymer) ) // if RNA 
                {
                    const RNA& myRNA=static_cast<const RNA&>(myBiopolymer); 
                    String myBiopolymerType("RNA"); 
                    String mySequence = extractSequenceFromBiopolymer(myRNA);
                    ResidueID myFirstResidueNumber ( 
                        myRNA.getResidue( ResidueInfo::Index(0)).getPdbResidueNumber(),
                        myRNA.getResidue( ResidueInfo::Index(0)).getPdbInsertionCode()
                        ); // retrieve first residue in chain. don't forget, we haven't actually called getPdbInsertionCode  () !  also, i think this doesn't take into account proteinCapping, if that's an issue. 
                    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" "<<mySequence<<endl;
                    addBiopolymerClass(
                     mySequence, 
                     myChainIdString,
                     myFirstResidueNumber, 
                     myBiopolymerType, 
                     proteinCapping,
                     inPDBFileName,
                     true, useNACappingHydroxyls);
                    //updBiopolymerClass(myChainIdString).setResidueNumbersAndInsertionCodesFromBiopolymer(myRNA, 'FALSE'     ); // provided biopolymer has end caps because it's a protein.
                    setOriginalSequence(myChainIdString,mySequence); 
                    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" mySequence "<<mySequence<<endl;
                    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" myChainIdString "<<myChainIdString<<endl;
                    //if (tempRenumberPdbResidues){
                    //   cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" "<<endl;
                    //   renumberPdbResidues(ResidueID("1"));
                    //}
                    setResidueIDsAndInsertionCodesFromBiopolymer(myChainIdString,myRNA,0); // provided biopolymer has no end caps because it's an RNA.
                    updBiopolymerClass(myChainIdString).setRenumberPdbResidues(tempRenumberPdbResidues);
                    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<endl;
                }
                else if (isDNA(myBiopolymer) ) // if DNA 
                {
                    const DNA& myDNA=static_cast<const DNA&>(myBiopolymer); 
                    String myBiopolymerType("DNA"); 
                    String mySequence = extractSequenceFromBiopolymer(myDNA);
                    ResidueID myFirstResidueNumber ( 
                        myDNA.getResidue( ResidueInfo::Index(0)).getPdbResidueNumber(),
                        myDNA.getResidue( ResidueInfo::Index(0)).getPdbInsertionCode()
                        ); // retrieve first residue in chain. don't forget, we haven't actually called getPdbInsertionCode  () !  also, i think this doesn't take into account proteinCapping, if that's an issue. 
                    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" "<<mySequence<<endl;
                    addBiopolymerClass(
                     mySequence, 
                     myChainIdString,
                     myFirstResidueNumber, 
                     myBiopolymerType, 
                     proteinCapping,
                     inPDBFileName,
                     true,useNACappingHydroxyls);
                    //if (tempRenumberPdbResidues){
                    //   renumberPdbResidues(ResidueID("1"));
                    //}
                    setResidueIDsAndInsertionCodesFromBiopolymer(myChainIdString,myDNA,0); // provided biopolymer has no end caps because it's a DNA.
                    updBiopolymerClass(myChainIdString).setRenumberPdbResidues(tempRenumberPdbResidues);
                }
                else if (isProtein(myBiopolymer)) {

                    const Protein& myProtein=static_cast<const Protein&>(myBiopolymer); 
                    String myBiopolymerType("Protein"); 
                    String mySequence = extractSequenceFromBiopolymer(myProtein, true);
                    ResidueID myFirstResidueNumber ( 
                        myProtein.getResidue( ResidueInfo::Index(1)).getPdbResidueNumber(),
                        myProtein.getResidue( ResidueInfo::Index(1)).getPdbInsertionCode()
                        ); // retrieve first residue in chain. don't forget, we haven't actually called getPdbInsertionCode  () !  also, i think this doesn't take into account proteinCapping, if that's an issue. 
                    //printBiopolymerSequenceInfo(myProtein);
                    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" "<<mySequence<<endl;
                    addBiopolymerClass(
                        mySequence, 
                        myChainIdString,
                        myFirstResidueNumber,
                        myBiopolymerType, 
                        proteinCapping,
                        inPDBFileName,
                        true, // this will setRenumberPdbResidues(-1)
			useNACappingHydroxyls); 
                    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<endl;
                    setResidueIDsAndInsertionCodesFromBiopolymer(myChainIdString,myProtein,1); // provided biopolymer has end caps because it's a protein.
                    updBiopolymerClass(myChainIdString).setRenumberPdbResidues(tempRenumberPdbResidues);
                    /*if (tempRenumberPdbResidues){
                       cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<endl;
                       updBiopolymerClass(myChainIdString).renumberPdbResidues(ResidueID("1"));
                    }*/
                    //cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<endl;
                }
                else {
                   ErrorManager::instance <<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" loadSequencesFromPdb can only be used with Protein, RNA and DNA.  You have one or more atoms in the file "<< inPDBFileName<< " which belong to none of these. Please get rid of these atoms and try again."<<endl;
                   ErrorManager::instance.treatError();
                }
                BiopolymerClass & myBiopolymerClass = updBiopolymerClass(myChainIdString);
                myBiopolymerClass.setPdbStructure((pdbStructureMap.at(inPDBFileName)));
            } // of if Biopolymer
        } // of if Molecule
        cout<<__FILE__<<":"<<__LINE__<<" This BiopolymerClassContainer now     has getNumBiopolymers() = "<< getNumBiopolymers() <<endl;
    } // for
    cout<<__FILE__<<":"<<__LINE__<<" Done adding compounds for now. This BiopolymerClassContainer now     has getNumBiopolymers() = "<< getNumBiopolymers() <<endl;

    //printBiopolymerSequenceInfo(updBiopolymerClass("g").myBiopolymer);
};

void BiopolymerClassContainer::resetAllPdbFileNames ( String newPdbFileNames )
{
    //================================================ Iterate through all biopolymers
    map<const String, BiopolymerClass>::iterator it;
    map<const String, BiopolymerClass>::iterator next;

    next                                              = this->biopolymerClassMap.begin();
    while (next != this->biopolymerClassMap.end())
    {
        it = next;
        (it->second).setPdbFileName                   ( newPdbFileNames );
        next++;
    }
}

void BiopolymerClassContainer::printBiopolymerInfo() {
    map<const String,BiopolymerClass>::iterator biopolymerClassMapIterator = biopolymerClassMap.begin();
    for(biopolymerClassMapIterator = biopolymerClassMap.begin(); biopolymerClassMapIterator != biopolymerClassMap.end(); biopolymerClassMapIterator++) {
        (biopolymerClassMapIterator->second).printBiopolymerInfo();
    }  
};

void BiopolymerClassContainer::setResidueIDsAndInsertionCodesFromBiopolymer(String chain, const Biopolymer & inputBiopolymer, bool endCaps){
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<endl;
    updBiopolymerClass(chain).setResidueIDsAndInsertionCodesFromBiopolymer(   inputBiopolymer,  endCaps);
}
ResidueID BiopolymerClassContainer::residueID(map<const String,double> myUserVariables,  const char* value , String chain) {
    String inputResidueID(value);
    size_t plusPosition = inputResidueID.find('+',1); 
    size_t minusPosition = inputResidueID.find('-',1); 
    size_t ampersandPosition = inputResidueID.find('@',0 ); // Look for "@" anywhere in "value" 
    //cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<endl;
    if (inputResidueID.compare("FirstResidue") == 0) {
        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" You have requested the first residue of chain "<<chain<<", which is : "<<updBiopolymerClass(chain).getFirstResidueID().outString()<<endl;
        return updBiopolymerClass(chain).getFirstResidueID();
    } 
    else if (inputResidueID.compare("LastResidue") == 0) {
        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" You have requested the last residue of chain "<<chain<<", which is : "<<updBiopolymerClass(chain).getLastResidueID().outString()<<endl;
        return updBiopolymerClass(chain).getLastResidueID();
    }
    //else if ((inputResidueID.substr(0,1)).compare("@") ==0) { // if the String starts with '@' , this is a user-defined integer variable.  Note that insertion codes cannot be specified with this method.
    //    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<endl;
    //    return ResidueID(myUserVariables, value  ); // This method knows how to do arithmetic involving literal integers and user variables, in any order. Or if it is just a user variable with no arithmetic, it can handle that also.
    //} 
    else if ((plusPosition != String::npos) || (minusPosition != String::npos) ){ // This is the case that there is a +/- operation to do
        size_t leftMostPlusMinus = min(plusPosition,minusPosition) ;
        String myResidueIDString = inputResidueID.substr(0, (leftMostPlusMinus + 0) ); // The part before the first +/- is assumed to be the residue ID.
        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" You have specified an arithmetic operation '+/-' be performed on a residue ID: "<<inputResidueID<<endl;
        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" The starting residue ID is taken to be: "<<myResidueIDString <<endl;
        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<endl;
        // I am NOT adding +1 to leftMostPlusMinus .. this means that the +/1 sign goes with the myResidueIncrementString.
        String myResidueIncrementString = (inputResidueID.substr(leftMostPlusMinus ,1000)); // the second parameter is ridiculously large, but will be truncated at the end of the input String.
        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<endl;
        stringstream myResidueIncrementStringStream(myResidueIncrementString);
        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<endl;
        int myResidueIncrement = -1111;
        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" myResidueIncrementString = "<<myResidueIncrementString<<endl; 
        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Note that the above should NOT include any insertion codes."<<endl; 
        myResidueIncrement = myAtoI(myUserVariables, myResidueIncrementString ); // This can handle additional +/- as well as user variables in any order or position
        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" myResidueIncrement = >"<<myResidueIncrement<<"<"<<endl;
	cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" chain: "<<chain<<endl;
        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" myResidueIDString >"<<myResidueIDString<<"<"<<endl;
        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" myResidueIDString.c_str() >"<<myResidueIDString.c_str()<<"<"<<endl;
        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" You wish to add the following increment : "<<myResidueIncrement<<" to the following residue ID: ";
        //cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<endl;
        ResidueID myResidueID = residueID(myUserVariables, myResidueIDString.c_str(), chain); // this is recursive, calls self. Should be able to handle LastResidue, FirstResidue, @ variables, and literal strings.
        cout<<myResidueID.outString()<<endl;
        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" The result is: "<<updBiopolymerClass(chain).sum(myResidueID,myResidueIncrement).outString()<<endl;
        return updBiopolymerClass(chain).sum(myResidueID,myResidueIncrement);
        
    }
    else if (ampersandPosition == 0           ){
        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<endl;
        return ResidueID(myUserVariables, value  ); // This method knows how to do arithmetic involving literal integers and user variables, in any order.  Or if it is just a user variable with no arithmetic, it can handle that also.  However if we are calling this, "value" does not have any +/- in this leaf. Any +/- is being done one layer up in the recursion. 
    }
    else { // must be a plain integer, or integer followed by insertion code .. 
        return updBiopolymerClass(chain).residueID(inputResidueID); }
    ErrorManager::instance <<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Unexplained error! Could not parse residue ID: >"<<inputResidueID<<"< of chain >"<<chain<<"<" <<endl;  //. Is it possible you tried to use a literal residue ID and insertion code AFTER (to the right of an arithmetic operator?" <<endl; // if we get to this point, something is wrong.
    // if we get to this point, something is wrong.
    // Is it possible you tried to use a literal residue ID and insertion code AFTER (to the right of an arithmetic operator?" <<endl; // if we get to this point, something is wrong.
    ErrorManager::instance.treatError();
}


void BiopolymerClassContainer::constrainAllChainsToEachOther(ConstraintToGroundContainer & constraintToGroundContainer){
    for (auto  biopolymerClassMapIterator = biopolymerClassMap.begin() ; biopolymerClassMapIterator != biopolymerClassMap.end(); biopolymerClassMapIterator++) {
        if (biopolymerClassMapIterator != biopolymerClassMap.begin()) { // For all but the first biopolymer
            addConstraint( ((biopolymerClassMap.begin())->first) , // chain ID of first BiopolymerClass
                           (biopolymerClassMapIterator->first)   , // chain ID of other BiopolymerClass
                           constraintToGroundContainer);           // Parameter reader has a ConstraintToGroundContainer member.
        }
    }
}

void BiopolymerClassContainer::addConstraintToGround(map<const String,double> myUserVariables, const String inputResidueString, const String chain, ConstraintToGroundContainer & constraintToGroundContainer){
    constraintToGroundContainer.addConstraintClassToVector( 
        chain,
        residueID(myUserVariables,
                  inputResidueString,
                  chain),
        updBiopolymerClass(chain).getRepresentativeAtomName() 
        ); 
}

void BiopolymerClassContainer::addConstraint(map<const String,double> myUserVariables,
                   const String inputResidueString, const String chain1, 
                   const String inputResidueString2, const String chain2, 
                   ConstraintToGroundContainer & constraintToGroundContainer)
{
    ResidueID residue1 = residueID(myUserVariables,inputResidueString,chain1);
    ResidueID residue2 = residueID(myUserVariables,inputResidueString2,chain2);
    addConstraint(residue1,chain1, residue2,chain2,constraintToGroundContainer);        
    /* constraintToGroundContainer.addConstraintToVector(
        chain,
        residueID(myUserVariables,inputResidueString,chain),
        updBiopolymerClass(chain).getRepresentativeAtomName(),
        chain2,
        residueID(myUserVariables,inputResidueString2,chain2),
        updBiopolymerClass(chain2).getRepresentativeAtomName() 
    ); */
}
void BiopolymerClassContainer::addConstraint(
                   const String chain1,
                   const String chain2,
                   ConstraintToGroundContainer & constraintToGroundContainer)
{
    constraintToGroundContainer.addConstraintToVector(
        chain1,
        updBiopolymerClass(chain1). getFirstResidueID(),
        updBiopolymerClass(chain1).getRepresentativeAtomName(),
        chain2,
        updBiopolymerClass(chain2). getFirstResidueID(),
        updBiopolymerClass(chain2).getRepresentativeAtomName()
    );
}

void BiopolymerClassContainer::addConstraint(
                   const ResidueID residue1, const String chain1,
                   const ResidueID residue2, const String chain2,
                   ConstraintToGroundContainer & constraintToGroundContainer)
{
    constraintToGroundContainer.addConstraintToVector(
        chain1,
        residue1,
        updBiopolymerClass(chain1).getRepresentativeAtomName(),
        chain2,
        residue2,
        updBiopolymerClass(chain2).getRepresentativeAtomName()
    );
}

// To add a constraint which specifies chain, residue, and atom name:
void BiopolymerClassContainer::addConstraint(map<const String,double> myUserVariables,
                   const String atomName, const String inputResidueString,const  String chain, 
                   const String atomName2, const String inputResidueString2,const  String chain2, 
                   ConstraintToGroundContainer & constraintToGroundContainer)
{
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<endl;
    ResidueID residueID1;
    if (hasChainID(chain)){ residueID1 = residueID(myUserVariables,inputResidueString,chain);}
    else {
        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Warning! chain >"<<chain<<"< is not a biopolymer! Constraints to non-biopolymers are a new feature still in beta!"<<endl;
        residueID1 = ResidueID(inputResidueString);}

    ResidueID residueID2;
    if (hasChainID(chain2)){ residueID2 = residueID(myUserVariables,inputResidueString2,chain2);}
    else { 
        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Warning! chain >"<<chain2<<"< is not a biopolymer! Constraints to non-biopolymers are a new feature still in beta!"<<endl;
        residueID2 = ResidueID(inputResidueString2);}

    constraintToGroundContainer.addConstraintToVector(
        chain, 
        residueID1,  
        atomName,
        chain2,
        residueID2,  
        atomName2 
    ); 
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<endl;
} 

void BiopolymerClassContainer::addConstraint(map<const String,double> myUserVariables,
                   const String atomName1, const String inputResidueString1,const  String chain1, 
                   const String atomName2, const String inputResidueString2,const  String chain2, 
                   ConstraintType myConstraintType,
                   ConstraintToGroundContainer & constraintToGroundContainer)
{
    ConstraintClass myConstraintClass;
    myConstraintClass.setChain1(chain1 );
    myConstraintClass.setChain2(chain2);
    myConstraintClass.setResidueID1( residueID(myUserVariables,inputResidueString1, myConstraintClass.getChain1() ));
    myConstraintClass.setResidueID2( residueID(myUserVariables,inputResidueString2, myConstraintClass.getChain2() ));
    myConstraintClass.setAtomName1(atomName1);
    myConstraintClass.setAtomName2(atomName2);

    // These two lines are to validate the atom name:
    updBiopolymerClass(myConstraintClass.getChain1()).atomPathString(myConstraintClass.getResidueID1(),myConstraintClass.getAtomName1());
    updBiopolymerClass(myConstraintClass.getChain2()).atomPathString(myConstraintClass.getResidueID2(),myConstraintClass.getAtomName2());
    //   
    myConstraintClass.setConstraintType(CoupledCoordinate);
        constraintToGroundContainer.addConstraintClassToVector(myConstraintClass);


    /*constraintToGroundContainer.addConstraintToVector(
        chain, 
        residueID(myUserVariables,inputResidueString,chain),  
        atomName,
        chain2,
        residueID(myUserVariables,inputResidueString2,chain2),  
        atomName2 
    );*/ 
} 

void BiopolymerClassContainer::constrainRigidSegmentsToGroundForAllChains(CompoundSystem & system,  SimbodyMatterSubsystem & matter,State & state, ConstraintToGroundContainer & myConstraintToGroundContainer  ) {

    map<const String, BiopolymerClass>::iterator biopolymerClassMapIterator = biopolymerClassMap.begin();
    for(biopolymerClassMapIterator = biopolymerClassMap.begin(); biopolymerClassMapIterator != biopolymerClassMap.end(); biopolymerClassMapIterator++) 
    {
        BiopolymerClass       myBiopolymerClass = (biopolymerClassMapIterator->second);
        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Constraining chain "<<myBiopolymerClass.getChainID()<<" to ground."<<endl;
        myBiopolymerClass.constrainRigidSegmentsToGround(system,  matter,state,  myConstraintToGroundContainer   );
    }
    //cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<"At the end of constrainRigidSegmentsToGroundForAllChains, running validateConstraintClassVector:"<<endl;
    myConstraintToGroundContainer.validateConstraintClassVector(*this); 
};



void BiopolymerClassContainer::loadResidueIDVector(){
    map<const String, BiopolymerClass>::iterator biopolymerClassMapIterator = biopolymerClassMap.begin();
    for(biopolymerClassMapIterator = biopolymerClassMap.begin(); biopolymerClassMapIterator != biopolymerClassMap.end(); biopolymerClassMapIterator++) {
        (biopolymerClassMapIterator->second).loadResidueIDVector();    
    }
}


void BiopolymerClassContainer::setFirstResidueMobilizerType(String myFirstResidueMobilizerType) {
    map<const String, BiopolymerClass>::iterator biopolymerClassMapIterator = biopolymerClassMap.begin();
    for(biopolymerClassMapIterator = biopolymerClassMap.begin(); biopolymerClassMapIterator != biopolymerClassMap.end(); biopolymerClassMapIterator++) {
        (biopolymerClassMapIterator->second).setFirstResidueMobilizerType(myFirstResidueMobilizerType);
    }
}

void BiopolymerClassContainer::setContactParameters ( GeneralContactSubsystem & contacts, HuntCrossleyForce & hc, double excludedVolumeStiffness, bool active ) {
    map<const String, BiopolymerClass>::iterator biopolymerClassMapIterator = biopolymerClassMap.begin();
    for(biopolymerClassMapIterator = biopolymerClassMap.begin(); biopolymerClassMapIterator != biopolymerClassMap.end(); biopolymerClassMapIterator++) {
        (biopolymerClassMapIterator->second).setContactParameters( contacts,    hc,  excludedVolumeStiffness, active );
    }
}

void BiopolymerClassContainer::setOriginalSequencesFromCurrentSequences() {
    map<const String, BiopolymerClass>::iterator biopolymerClassMapIterator = biopolymerClassMap.begin();
    for(biopolymerClassMapIterator = biopolymerClassMap.begin(); biopolymerClassMapIterator != biopolymerClassMap.end(); biopolymerClassMapIterator++) {
        (biopolymerClassMapIterator->second).setOriginalSequence((biopolymerClassMapIterator->second).getSequence() );
    }
}



void BiopolymerClassContainer::substituteResidue(String myChain , ResidueID myResidue, String mySubstitution, bool proteinCapping) {
    BiopolymerClass myOldBiopolymerClass = updBiopolymerClass(myChain);

    String myOldSequence = myOldBiopolymerClass.getSequence();
    String myOriginalSequence = myOldBiopolymerClass.getOriginalSequence();
    String myNewSequence = myOldSequence;
    ResidueID myFirstResidueNumber = myOldBiopolymerClass.getFirstResidueID();
    myNewSequence[myOldBiopolymerClass.getResidueIndex( myResidue) ] = *(mySubstitution.c_str()); // careful! getResidueIndex would potentially be wrong .. here we want the first letter of the sequence to correspond to position zero, with no regard to proteinCapping.
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": old sequence = "<<myOldSequence<<endl;
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": new sequence = "<<myNewSequence<<endl;

    Biopolymer tempBiopolymer = myOldBiopolymerClass.myBiopolymer;
    replaceBiopolymerWithMutatedBiopolymerClass(myOldBiopolymerClass, myNewSequence);
    updBiopolymerClass(myChain).setResidueIDsAndInsertionCodesFromBiopolymer(tempBiopolymer, proteinCapping);
}

void BiopolymerClassContainer::replaceBiopolymerWithMutatedBiopolymerClass(BiopolymerClass & myOldBiopolymerClass, 
                                                            String & myNewSequence, bool useNACappingHydroxyls )
{
    String myChain = myOldBiopolymerClass.getChainID();
    ResidueID myFirstResidueNumber = myOldBiopolymerClass.getFirstResidueID();
    bool proteinCapping = myOldBiopolymerClass.getProteinCapping();
    String myOriginalSequence = myOldBiopolymerClass.getOriginalSequence();
    String oldBiopolymerClassBiopolymerType = myOldBiopolymerClass.getBiopolymerTypeAsString();
    String oldBiopolymerClassPdbFileName = myOldBiopolymerClass.getPdbFileName();
    bool oldBiopolymerClassLoadFromPdb = myOldBiopolymerClass.getLoadFromPdb();
    bool oldActivePhysics = myOldBiopolymerClass.getActivePhysics();
    deleteBiopolymerClass(myChain);
    if (hasChainID(myChain)){
        ErrorManager::instance <<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": Unexplained error!"<<endl;
        ErrorManager::instance.treatError();
    }
    addBiopolymerClass(myNewSequence,myChain, myFirstResidueNumber ,
                       oldBiopolymerClassBiopolymerType  ,proteinCapping,
                       oldBiopolymerClassPdbFileName, oldBiopolymerClassLoadFromPdb, useNACappingHydroxyls);
    updBiopolymerClass(myChain).setActivePhysics(oldActivePhysics);
    setOriginalSequence(myChain,myOriginalSequence);
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": Restoring residue numbers and insertion codes after mutating.. "<<endl;
}


void BiopolymerClassContainer::loadMutationVectorsFromSequence() {
    if (mutationVector.size() > 0) {
        //cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" There are already some mutations in mutationVector !  "<<endl; exit(1);
        ErrorManager::instance<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" There are already some mutations in mutationVector ! "<<endl; 
        ErrorManager::instance.treatError();
    }
    map<const String, BiopolymerClass>::iterator biopolymerClassMapIterator = biopolymerClassMap.begin();
    for(biopolymerClassMapIterator = biopolymerClassMap.begin(); biopolymerClassMapIterator != biopolymerClassMap.end(); biopolymerClassMapIterator++) 
    {
	    for (unsigned int i = 0 ; i < (biopolymerClassMapIterator->second).getSequence().size() ; i ++) {
		if ((biopolymerClassMapIterator->second).getSequence().substr(i,1).compare((biopolymerClassMapIterator->second).getOriginalSequence().substr(i,1)) != 0 ) { 
		    std::cout<<"Found a mutation at residue index "<<i<<std::endl; 
		    Mutation myMutation;
		    myMutation.setChain   ((biopolymerClassMapIterator->second).getChainID());
		    myMutation.setResidue ((biopolymerClassMapIterator->second).getResidueID(i)  );
		    myMutation.setSubstitutedResidueType((biopolymerClassMapIterator->second).getSequence().substr(i,1));
		    addMutationToVector(myMutation);
		    //mutationVector.push_back(myMutation); 
		    //myNumMutations++;
		} else { // do nothing
		}   
	    }
    //cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Done with loadMutationVectorFromSequence .. mutations look like: "<< getFormattedMutationsString(MUTATIONMINORSEPARATOR )<<endl;
    }
}

/*void BiopolymerClassContainer::writeMutationFlexibilizers(std::ofstream & output, const int offset , const double radius = 0.0) {
    map<const String, BiopolymerClass>::iterator biopolymerClassMapIterator = biopolymerClassMap.begin();
    for(biopolymerClassMapIterator = biopolymerClassMap.begin(); biopolymerClassMapIterator != biopolymerClassMap.end(); biopolymerClassMapIterator++) {
        (biopolymerClassMapIterator->second).writeMutationFlexibilizers( output,  offset, radius);
    }
}*/
void BiopolymerClassContainer::writeMutationFlexibilizers(std::ofstream & output, const int offset, const double radius = 0.0 ) {
                int leftFlexibleOffset = offset;
                int rightFlexibleOffset = offset;
                for (int i = 0 ; i <       mutationVector.size(); i++) {
                        output <<"mobilizer Default "<<mutationVector[i].getChain()<<" " ;
                        output<<updBiopolymerClass(mutationVector[i].getChain()).safeSum(mutationVector[i].getResidue(),(- leftFlexibleOffset)).outString()<<" ";
                        output <<updBiopolymerClass(mutationVector[i].getChain()).safeSum(mutationVector[i].getResidue(),rightFlexibleOffset).outString()<<std::endl;
                        output <<"applyMobilizersWithin Default "<<radius<<" "<<mutationVector[i].getChain()<<" "<<mutationVector[i].getResidue().outString()<<std::endl;
                }
}

/*void BiopolymerClassContainer::writeWaterDroplets(std::ofstream & output, const double springConstant = 300, const double radius = 0.0 ) {
    map<const String, BiopolymerClass>::iterator biopolymerClassMapIterator = biopolymerClassMap.begin();
    for(biopolymerClassMapIterator = biopolymerClassMap.begin(); biopolymerClassMapIterator != biopolymerClassMap.end(); biopolymerClassMapIterator++) {
        (biopolymerClassMapIterator->second).writeWaterDroplets( output,  springConstant, radius);
    }
}*/
void BiopolymerClassContainer::writeWaterDroplets(std::ofstream & output, const double springConstant = 300, const double radius = 0.0 ) {
                for (int i = 0 ; i <       mutationVector.size   (); i++) {
                        output <<"waterDropletAboutResidue "<<mutationVector[i].getChain()<<" " ;
                        output<<mutationVector[i].getResidue().outString()<<" ";
                        output <<radius<<" ";
                        output <<springConstant<<" "<< mutationVector[i].getMutationAsString()  <<std::endl; // getMutationAsString() will return a multi-character string.  MMB knows how to handle this. However in the trajectory file, the chain ID will be " ", with remarks telling MMB how to parse.
                }
}

/*void BiopolymerClassContainer::writeMobilizerWithinMutation(std::ofstream & output,  const double radius = 0.0 ) {
    map<const String, BiopolymerClass>::iterator biopolymerClassMapIterator = biopolymerClassMap.begin();
    for(biopolymerClassMapIterator = biopolymerClassMap.begin(); biopolymerClassMapIterator != biopolymerClassMap.end(); biopolymerClassMapIterator++) {
        (biopolymerClassMapIterator->second).writeMobilizerWithinMutation( output,  radius);
    }
}*/

void BiopolymerClassContainer::writeMobilizerWithinMutation(std::ofstream & output,  const double radius = 0.0 ) {
                for (int i = 0 ; i < mutationVector.size(); i++) {
                        output <<"applyMobilizersWithin Default  "<<radius<<" "<<mutationVector[i].getChain()<<" " ;
                        output <<mutationVector[i].getResidue().outString()<<std::endl;
                }
}



/*void BiopolymerClassContainer::writeMutationBackboneRigidifier (std::ofstream & output, const int offset) {
    map<const String, BiopolymerClass>::iterator biopolymerClassMapIterator = biopolymerClassMap.begin();
    for(biopolymerClassMapIterator = biopolymerClassMap.begin(); biopolymerClassMapIterator != biopolymerClassMap.end(); biopolymerClassMapIterator++) {
        (biopolymerClassMapIterator->second).writeMutationBackboneRigidifier ( output,  offset);
    }
}*/

/*void BiopolymerClassContainer::writePhysicsZones(std::ofstream & output, const int offset) {
    map<const String, BiopolymerClass>::iterator biopolymerClassMapIterator = biopolymerClassMap.begin();
    for(biopolymerClassMapIterator = biopolymerClassMap.begin(); biopolymerClassMapIterator != biopolymerClassMap.end(); biopolymerClassMapIterator++) {
        (biopolymerClassMapIterator->second).writePhysicsZones( output, offset);
    }
}*/

const       vector<MMBAtomInfo> BiopolymerClassContainer::getConcatenatedAtomInfoVector(bool activeChainsOnly) {
    vector<MMBAtomInfo> myAtomInfoVector;
    vector<MMBAtomInfo> tempAtomInfoVector;
    myAtomInfoVector.clear();
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<endl;
    map<const String, BiopolymerClass>::iterator biopolymerClassMapIterator = biopolymerClassMap.begin();
    for (biopolymerClassMapIterator = biopolymerClassMap.begin(); biopolymerClassMapIterator != biopolymerClassMap.end(); biopolymerClassMapIterator++) {
        if(activeChainsOnly && !biopolymerClassMapIterator->second.getActivePhysics())
            continue;
        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Inside getConcatenatedAtomInfoVector(). Doing chain : "<<(biopolymerClassMapIterator->second).getChainID()<<endl;
        tempAtomInfoVector = (biopolymerClassMapIterator->second).getAtomInfoVector(); 
    myAtomInfoVector.insert(myAtomInfoVector.end(),tempAtomInfoVector.begin(), tempAtomInfoVector.end());
    }
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" At the end of getConcatenatedAtomInfoVector(). The returned vector has length : "<<myAtomInfoVector.size()<<endl;
    return myAtomInfoVector;
}

const       vector<MMBAtomInfo> BiopolymerClassContainer::getConcatenatedAtomInfoVector(const State & state,bool activeChainsOnly) {
    vector<MMBAtomInfo> myAtomInfoVector;
    vector<MMBAtomInfo> tempAtomInfoVector;
    myAtomInfoVector.clear();
    map<const String, BiopolymerClass>::iterator biopolymerClassMapIterator = biopolymerClassMap.begin();
    for (biopolymerClassMapIterator = biopolymerClassMap.begin(); biopolymerClassMapIterator != biopolymerClassMap.end(); biopolymerClassMapIterator++) {
        if(activeChainsOnly && !biopolymerClassMapIterator->second.getActivePhysics())
            continue;
        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Inside getConcatenatedAtomInfoVector(). Doing chain : "<<(biopolymerClassMapIterator->second).getChainID()<<endl;
        tempAtomInfoVector = (biopolymerClassMapIterator->second).getAtomInfoVector();
        for (int m = 0; m < (int)tempAtomInfoVector.size(); m++) 
        {
            MMBAtomInfo & tempAtomInfo = tempAtomInfoVector[m];
            Vec3 v = (biopolymerClassMapIterator->second).myBiopolymer.calcAtomLocationInGroundFrame(state, tempAtomInfo.compoundAtomIndex);
            tempAtomInfo.position[0] = v[0];
            tempAtomInfo.position[1] = v[1];
            tempAtomInfo.position[2] = v[2];
        } 
    myAtomInfoVector.insert(myAtomInfoVector.end(),tempAtomInfoVector.begin(), tempAtomInfoVector.end());
    }
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" At the end of getConcatenatedAtomInfoVector(). The returned vector has length : "<<myAtomInfoVector.size()<<endl;
    return myAtomInfoVector;
}
 

void BiopolymerClassContainer::printAtomInfoVector() {
    map<const String, BiopolymerClass>::iterator biopolymerClassMapIterator = biopolymerClassMap.begin();
    for(biopolymerClassMapIterator = biopolymerClassMap.begin(); biopolymerClassMapIterator != biopolymerClassMap.end(); biopolymerClassMapIterator++) {
        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Printing atomInfoVector for chain "<<(biopolymerClassMapIterator->second).getChainID()<<endl;
        (biopolymerClassMapIterator->second).printAtomInfoVector();
    }
}
/*
void BiopolymerClassContainer::writeSubstituteResidueCommands(std::ofstream & output) {
    map<const String, BiopolymerClass>::iterator biopolymerClassMapIterator = biopolymerClassMap.begin();
    for(biopolymerClassMapIterator = biopolymerClassMap.begin(); biopolymerClassMapIterator != biopolymerClassMap.end(); biopolymerClassMapIterator++) {
        (biopolymerClassMapIterator->second).writeSubstituteResidueCommands( output);
    }
}*/
void BiopolymerClassContainer::writeSubstituteResidueCommands(std::ofstream & output) {
    for (int i = 0 ; i <       mutationVector.size   (); i++) {
    output <<"substituteResidue "<<mutationVector[i].getChain()<<" "<<mutationVector[i].getResidue().outString()<<" "<<mutationVector[i].getSubstitutedResidueType()<<std::endl;
    }
}

int BiopolymerClassContainer::getNumMutationVectorElements() {
    /*int numMutations  = 0;
    map<const String, BiopolymerClass>::iterator biopolymerClassMapIterator = biopolymerClassMap.begin();
    for(biopolymerClassMapIterator = biopolymerClassMap.begin(); biopolymerClassMapIterator != biopolymerClassMap.end(); biopolymerClassMapIterator++) {
        numMutations += (biopolymerClassMapIterator->second).getNumMutationVectorElements();
    }*/
    return mutationVector.size();//numMutations;
}
String BiopolymerClassContainer::getFormattedSequencesString() {
    String sequencesString = "";    
    map<const String, BiopolymerClass>::iterator biopolymerClassMapIterator = biopolymerClassMap.begin();
    for (biopolymerClassMapIterator = biopolymerClassMap.begin(); biopolymerClassMapIterator != biopolymerClassMap.end(); biopolymerClassMapIterator++) {
         sequencesString += (biopolymerClassMapIterator->second).getChainID();
         sequencesString += ":"; // separates chain ID from sequence
         sequencesString += (biopolymerClassMapIterator->second).getSequence();
         biopolymerClassMapIterator ++ ;
         if ((biopolymerClassMapIterator) != biopolymerClassMap.end())  // if this is not the last chain
             sequencesString += "."; // a "." connects chains
         biopolymerClassMapIterator -- ; // decrement biopolymerClassMapIterator again
    }
    if (sequencesString.length() == 0) {
        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Unexpectedly, the getFormattedSequencesString() string is of length zero."<<endl; exit(1);
    } 
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" sequencesString = "<<sequencesString<<endl;
    return sequencesString;
}

String BiopolymerClassContainer::getFormattedMutationsString(String minorSeparator =  MUTATIONMINORSEPARATOR) {
    String sequencesString = "";             
    for (int i = 0 ; i < getNumMutationVectorElements() ; i++) {
        sequencesString += mutationVector[i].getChain();
        sequencesString += minorSeparator; // separates chain ID from residueID
        sequencesString +=  mutationVector[i].getResidue().outString(); 
        sequencesString += minorSeparator; // separates residueID from mutation
    sequencesString += mutationVector[i].getSubstitutedResidueType();
        if (i < (getNumMutationVectorElements() -1)) {
        sequencesString += MUTATIONMAJORSEPARATOR ; // connect mutants // was "."
        }
    }    
    return sequencesString;
}

String BiopolymerClassContainer::getFoldxFormattedMutations() {

    String combinedMutationString = "";

    for (int i = 0 ; i < getNumMutationVectorElements() ; i++) {
        Mutation myMutation;
        myMutation = mutationVector[i];
        //setMutationWildTypeResidueType(myMutation);
        setMutationWildTypeResidueTypeFromOriginalSequence(myMutation);
        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" myMutation.getWildTypeResidueType : >"<<myMutation.getWildTypeResidueType()<<"< "<<endl;  
        combinedMutationString += myMutation.getMutationAsFoldxString();
        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<"  myMutation.getMutationAsFoldxString()  : >"<< myMutation.getMutationAsFoldxString() <<"< "<<endl;  
        if (i < (getNumMutationVectorElements() -1)) {
        combinedMutationString += FOLDXSEPARATOR ; // connect single mutants with ","
        }
    }
    return combinedMutationString;
}

Mutation  BiopolymerClassContainer::setMutationWildTypeResidueType(Mutation & myMutation){
    BiopolymerClass myBiopolymerClass = updBiopolymerClass(myMutation. getChain());
    ResidueID myResidue = myBiopolymerClass.residueID(myMutation.getResidue().outString()); // This BiopolymerClass method has a validation step. Requires a String.
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" myBiopolymerClass.residueID(myMutation.getResidue().outString()) returns >"<<myBiopolymerClass.residueID(myMutation.getResidue().outString()).outString()<<"< "<<endl; 
    String myWildTypeResidueType = myBiopolymerClass.getResidueSingleLetterCode(myResidue);
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" myBiopolymerClass.getResidueSingleLetterCode(myResidue) = >"<<myBiopolymerClass.getResidueSingleLetterCode(myResidue)<<"< "<<endl; 
    myMutation. setWildTypeResidueType(myWildTypeResidueType);
    return myMutation;
} ;

// This is a variation of setMutationWildTypeResidueType. Sometimes the current sequence is mutated, so the original residue type is lost. This is a way of recovering it. 
Mutation  BiopolymerClassContainer::setMutationWildTypeResidueTypeFromOriginalSequence(Mutation & myMutation){
    BiopolymerClass myBiopolymerClass = updBiopolymerClass(myMutation. getChain());
    ResidueInfo::Index myResidueIndex = myBiopolymerClass.getResidueIndex(myMutation.getResidue());
    String myOriginalWildTypeResidueType = myBiopolymerClass.getOriginalSequence().substr(myResidueIndex,1);
    myMutation. setWildTypeResidueType(myOriginalWildTypeResidueType);
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<endl;
    myMutation.print();
    /*
    ResidueID myResidue = myBiopolymerClass.residueID(myMutation.getResidue().outString()); // This BiopolymerClass method has a validation step. Requires a String.
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" myBiopolymerClass.residueID(myMutation.getResidue().outString()) returns >"<<myBiopolymerClass.residueID(myMutation.getResidue().outString()).outString()<<"< "<<endl; 
    String myWildTypeResidueType = myBiopolymerClass.getResidueSingleLetterCode(myResidue);
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" myBiopolymerClass.getResidueSingleLetterCode(myResidue) = >"<<myBiopolymerClass.getResidueSingleLetterCode(myResidue)<<"< "<<endl; 
    myMutation. setWildTypeResidueType(myWildTypeResidueType);
    */ 
    return myMutation;
} ;

/*String BiopolymerClassContainer::getFormattedMutationsString(String minorSeparator = MUTATIONMINORSEPARATOR) {
    String mutationString = ""; 
    map<const String, BiopolymerClass>::iterator biopolymerClassMapIterator = biopolymerClassMap.begin();
    for (biopolymerClassMapIterator = biopolymerClassMap.begin(); biopolymerClassMapIterator != biopolymerClassMap.end(); biopolymerClassMapIterator++) {
     if (mutationString.substr(0,1).compare(".") == 0)
	 mutationString = mutationString.substr(1,(mutationString.length()-1)); // get rid of any leading periods left in the last cycle of this loop.
         cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" chain ID = >"<<(biopolymerClassMapIterator->second).getChainID()<<"<"<<endl;
         cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" current mutationString = >"<< mutationString <<"<"<<endl;
         if ((biopolymerClassMapIterator->second).getFormattedMutationsString(minorSeparator ).size() > 0) 
         {
             cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" adding formatted mutation string  = >"<<(biopolymerClassMapIterator->second).getFormattedMutationsString(minorSeparator )<<  "<"<<endl;;
             if (mutationString.length() > 0 ) {
                 cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" plus string = >"<<mutationString<<"<"<<endl;
                 mutationString =  (biopolymerClassMapIterator->second).getFormattedMutationsString(minorSeparator )+String(".")+mutationString  ; // += seems to do current + new , rather than vice versa.
             } else {
                 cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" .. to no other string = >"<<mutationString<<"<"<< endl;
                 mutationString =  (biopolymerClassMapIterator->second).getFormattedMutationsString(minorSeparator );
             }
         } else {
             cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" mutation string to be added had no length.. doing nothing."<<endl;
         }
         //mutationString += (biopolymerClassMapIterator->second).getFormattedMutationsString(minorSeparator ); // This was reversing the order of the mutations in the string
         cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" mutationString is now >"<<mutationString<<"<"<<endl;
         //biopolymerClassMapIterator ++ ;
         //if ((biopolymerClassMapIterator) != biopolymerClassMap.end())  // if this is not the last chain
         //mutationString += "."; // a "." connects mutations
         //biopolymerClassMapIterator -- ; // decrement biopolymerClassMapIterator again
    }
    for (int i = 0; i < getNumBiopolymers() ; i++) // there could be one "." added for each chain
            // This may have become redundant due to more careful placement of major separators.  But safer to leave in..
            if (mutationString.length() > 0)
            if (mutationString.substr((mutationString.length()-1),1).compare(".") == 0) {
                cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Dropping one trailing \'.\' "<<endl;
                mutationString = mutationString.substr(0,(mutationString.length()-1)); // get rid of trailing periods.
            }
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" mutationString is now >"<<mutationString<<"<"<<endl;
    return mutationString;
}*/

void BiopolymerClassContainer::setCurrentSequencesFromOriginalSequences() {
    map<const String, BiopolymerClass>::iterator biopolymerClassMapIterator = biopolymerClassMap.begin();
    for(biopolymerClassMapIterator = biopolymerClassMap.begin(); biopolymerClassMapIterator != biopolymerClassMap.end(); biopolymerClassMapIterator++) {
        (biopolymerClassMapIterator->second).setCurrentSequencesFromOriginalSequences(); //setSequence((biopolymerClassMapIterator->second).getOriginalSequence() );
    }
}

/*
bool BiopolymerClassContainer::allMutationsDifferFromWildType(){
    map<const String, BiopolymerClass>::iterator biopolymerClassMapIterator = biopolymerClassMap.begin();
    for(biopolymerClassMapIterator = biopolymerClassMap.begin(); biopolymerClassMapIterator != biopolymerClassMap.end(); biopolymerClassMapIterator++) {
        if (!(biopolymerClassMapIterator->second).allMutationsDifferFromWildType()) return false;
    }
    return true;
}*/

bool BiopolymerClassContainer::allMutationsDifferFromWildType() { // This tells us whether any of the proposed mutants actually do not  change the residue type at the specified position.
    for (unsigned int i = 0 ; i <       mutationVector.size   () ; i ++) {
    String updatedResidueType = updBiopolymerClass(mutationVector[i].getChain()).getOriginalSequence().substr(updBiopolymerClass(mutationVector[i].getChain()). getResidueIndex(mutationVector[i].getResidue()) ,1);
    if (updatedResidueType.compare(mutationVector[i].getSubstitutedResidueType()) == 0) {
                cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" The substituted residue type : >"<<mutationVector[i].getSubstitutedResidueType()<<"< is the same as the existing residue type : >"<<updatedResidueType<<endl;
                return false;
    }
    }
    return true; 
};
/*
void BiopolymerClassContainer::updateMutationResidueTypesFromCurrentSequence() {
    map<const String, BiopolymerClass>::iterator biopolymerClassMapIterator = biopolymerClassMap.begin();
    for(biopolymerClassMapIterator = biopolymerClassMap.begin(); biopolymerClassMapIterator != biopolymerClassMap.end(); biopolymerClassMapIterator++) {
        (biopolymerClassMapIterator->second).updateMutationResidueTypesFromCurrentSequence();
    }
}    */
void BiopolymerClassContainer::updateMutationResidueTypesFromCurrentSequence() {
    for (unsigned int i = 0 ; i < getNumMutationVectorElements() ; i ++) {
        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<endl;
        mutationVector[i].print(); 
        String updatedResidueType = updBiopolymerClass(mutationVector[i].getChain()).getSequence().substr(updBiopolymerClass(mutationVector[i].getChain()).getResidueIndex(mutationVector[i].getResidue()) ,1);
        if (updatedResidueType.compare(mutationVector[i].getSubstitutedResidueType()) == 0) {
            ErrorManager::instance<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" The substitution has not changed!  This is suspicious!"<<std::endl; ErrorManager::instance.treatError();
        }
        mutationVector[i].setSubstitutedResidueType( updatedResidueType );
        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<endl;
        mutationVector[i].print(); 
        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Done with updateMutationResidueTypesFromCurrentSequence .. mutations look like: "<< getFormattedMutationsString()<<endl;
    }
};


bool BiopolymerClassContainer::setRenumberPdbResidues (bool myRenumberPdbResidues){ 
    for (auto  biopolymerClassMapIterator = biopolymerClassMap.begin() ; biopolymerClassMapIterator != biopolymerClassMap.end(); biopolymerClassMapIterator++) { 
        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" About to setRenumberPdbResidues("<<myRenumberPdbResidues<<") for chain "<<biopolymerClassMapIterator->first      <<endl;
        updBiopolymerClass(biopolymerClassMapIterator->first).setRenumberPdbResidues(myRenumberPdbResidues);  
    }
}


void BiopolymerClassContainer::addMutationToVector(Mutation myMutation) { 
    updBiopolymerClass(myMutation.getChain()).validateMutation(myMutation);
    mutationVector.push_back(myMutation);
}

void BiopolymerClassContainer::substituteResidue(Mutation myMutation, 
                                                 bool safeParameters, 
                                                 bool matchPurineN1AtomLocations, 
                                                 bool proteinCapping) 
{
    String myChain = myMutation.getChain();
    ResidueID myResidue = myMutation.getResidue();
    //cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" >"<<myResidue.getInsertionCode()<<endl;
    String mySubstitution = myMutation.getSubstitutedResidueType();
    BiopolymerClass myOldBiopolymerClass = updBiopolymerClass(myChain);
    if (safeParameters) if  (myOldBiopolymerClass.getBiopolymerType() != BiopolymerType::Protein ) if (matchPurineN1AtomLocations) {
        ErrorManager::instance<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": In order to substitute a nucleic acid residue, you must first set matchPurineN1AtomLocations FALSE.  Otherwise you might mutate a purine to pyrmidine, and the N1 atom of the watson-crick edge would be taken as the glycosidic nitrogen of the pyrimidine, generating a physically irrational structure in the mutant."<<endl;
        ErrorManager::instance.treatError();
    }
    String myOldSequence = myOldBiopolymerClass.getSequence();
    String myOriginalSequence = myOldBiopolymerClass.getOriginalSequence();
    String myNewSequence = myOldSequence;
    ResidueID myFirstResidueNumber = myOldBiopolymerClass.getFirstResidueID();
    myNewSequence[myOldBiopolymerClass.getResidueIndex( myResidue) ] = *(mySubstitution.c_str()); // careful! getResidueIndex would potentially be wrong .. here we want the first letter of the sequence to correspond to position zero, with no regard to proteinCapping.
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": old sequence = "<<myOldSequence<<endl;
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": new sequence = "<<myNewSequence<<endl;    

    Biopolymer tempBiopolymer = myOldBiopolymerClass.myBiopolymer;
    replaceBiopolymerWithMutatedBiopolymerClass(myOldBiopolymerClass, myNewSequence);
    updBiopolymerClass(myChain).setResidueIDsAndInsertionCodesFromBiopolymer(tempBiopolymer, proteinCapping);
}

// this deletes a given residue from a given chain
void BiopolymerClassContainer::deleteResidue(Mutation myDeletion,   bool proteinCapping) {
    String myChain = myDeletion.getChain();
    ResidueID myDeletedResidueID = myDeletion.getResidue();
    String mySubstitution = myDeletion.getSubstitutedResidueType();
    BiopolymerClass myOldBiopolymerClass = updBiopolymerClass(myChain);
    String myOldSequence = myOldBiopolymerClass.getSequence();
    String myOriginalSequence = myOldBiopolymerClass.getOriginalSequence();
    String myNewSequence = myOldSequence;
    ResidueID myFirstResidueNumber = myOldBiopolymerClass.getFirstResidueID();;
    if (myDeletedResidueID == myOldBiopolymerClass.getFirstResidueID())  { myOldBiopolymerClass.incrementResidueID(myFirstResidueNumber);} // If we're deleting the first residue of the chain, the second residue becomes the first.

    ResidueInfo::Index deletedResidueIndex = myOldBiopolymerClass.getResidueIndex(myDeletedResidueID);  // Residue index of the residue to be deleted.  As numbered in the "old" biopolymer, of course.

    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": old sequence = "<<myOldSequence<<endl;
    myNewSequence.erase(deletedResidueIndex,1);
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": new sequence = "<<myNewSequence<<endl;

    Biopolymer tempBiopolymer = myOldBiopolymerClass.myBiopolymer;
    replaceBiopolymerWithMutatedBiopolymerClass(myOldBiopolymerClass, myNewSequence);
    updBiopolymerClass(myChain).setResidueIDsAndInsertionCodesFromBiopolymerWithDeletion(tempBiopolymer, deletedResidueIndex, proteinCapping);      
}

// this inserts a residue into an existing chain. myInsertion contains the chain ID and residue ID of the residue to insert.  The location of the insertion will be deduced by the existing residue IDs of the chain, under the assumption that alphabetical PDB ordering is to be respected.
 
void BiopolymerClassContainer::insertResidue(Mutation myInsertion,   bool proteinCapping) {
    String myChain = myInsertion.getChain();
    ResidueID myInsertedResidueID = myInsertion.getResidue();
    //cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" >"<<myInsertedResidueID.getInsertionCode()<<endl;
    String mySubstitution = myInsertion.getSubstitutedResidueType();
    BiopolymerClass myOldBiopolymerClass = updBiopolymerClass(myChain);
    String myOldSequence = myOldBiopolymerClass.getSequence();
    String myOriginalSequence = myOldBiopolymerClass.getOriginalSequence();
    String myNewSequence = myOldSequence;
    ResidueID myFirstResidueNumber = myOldBiopolymerClass.getFirstResidueID();
    ResidueInfo::Index insertedResidueIndex ; //= ResidueInfo::Index(-1111);
    // Now we will deduce the residue index for the insertion:
    for (ResidueID tempResidueID = myFirstResidueNumber; tempResidueID <= myOldBiopolymerClass.getLastResidueID(); myOldBiopolymerClass.incrementResidueID(tempResidueID)) {
        if (myInsertedResidueID == tempResidueID) {
            ErrorManager::instance<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": You have tried to insert residue >"<< myInsertedResidueID.outString()<<"< but this residue already exists!"<<endl;
            ErrorManager::instance.treatError();
        }
        if (myInsertedResidueID < tempResidueID) { // compare on the basis of residue number and insertion code (alphabetical)
            insertedResidueIndex  = myOldBiopolymerClass.getResidueIndex(tempResidueID);
            break; // exit for loop
        } else if (tempResidueID == myOldBiopolymerClass.getLastResidueID()) {
            insertedResidueIndex = myOldBiopolymerClass.getResidueIndex(tempResidueID) ; insertedResidueIndex++; // Append the insertion after the end of the chain.
            break; // exit for loop
        } else {} // Do nothing; we are not at the end of the loop so increment tempResidueID and keep looking.
    }

    if (mySubstitution.length() != 1) {
        ErrorManager::instance << __FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Substitution should be in single letter code!  Yours : >"<<mySubstitution<<"< is not exactly one character long." << endl;
        ErrorManager::instance.treatError();
    }
     
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": old sequence = "<<myOldSequence<<endl;
    myNewSequence.insert(insertedResidueIndex,mySubstitution)  ; 
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": new sequence = "<<myNewSequence<<endl;  

    Biopolymer tempBiopolymer = myOldBiopolymerClass.myBiopolymer;
    replaceBiopolymerWithMutatedBiopolymerClass(myOldBiopolymerClass, myNewSequence);
    updBiopolymerClass(myChain).setResidueIDsAndInsertionCodesFromBiopolymer(tempBiopolymer, myInsertion, proteinCapping);  
}

/*void        BiopolymerClassContainer::deleteResidue(String chainID, Residue startResidue, Residue endResidue) {

            String myChain = myMutation.getChainID();

            if (safeParameters) if  (myOldBiopolymerClass.getBiopolymerType() != BiopolymerType::Protein ) if (matchPurineN1AtomLocations) {
                ErrorManager::instance<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": In order to substitute a nucleic acid residue, you must first set matchPurineN1AtomLocations FALSE.  Otherwise you might mutate a purine to pyrmidine, and the N1 atom of the watson-crick edge would be taken as the glycosidic nitrogen of the pyrimidine, generating a physically irrational structure in the mutant."<<endl;
                ErrorManager::instance.treatError();
            }

            for (ResidueID myResidue = startResidue; myResidue <= endResidue;  updBiopolymerClass(chainID).incrementResidue(myResidue)) {
                BiopolymerClass myOldBiopolymerClass = updBiopolymerClass(myChain);
        String myOldSequence = myOldBiopolymerClass.getSequence();
        String myOriginalSequence = myOldBiopolymerClass.getOriginalSequence();
        String myNewSequence = myOldSequence;
                myNewSequence[myOldBiopolymerClass.getResidueIndex( myResidue) ] = *(mySubstitution.c_str()); // careful! getResidueIndex would potentially be wrong .. here we want the first letter of the sequence to correspond to position zero, with no regard to proteinCapping.
        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": old sequence = "<<myOldSequence<<endl;
        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": new sequence = "<<myNewSequence<<endl;
            
        }

            ResidueID myFirstResidueNumber = myOldBiopolymerClass.getFirstResidueID();
            Biopolymer tempBiopolymer = myOldBiopolymerClass. myBiopolymer;
            String oldBiopolymerClassBiopolymerType = myOldBiopolymerClass.getBiopolymerTypeAsString();
            deleteBiopolymerClass(myChain);
            if (hasChainID(myChain)){
                ErrorManager::instance <<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": Unexplained error!"<<endl;
                ErrorManager::instance.treatError();
            }
            addBiopolymerClass(myNewSequence,myChain, myFirstResidueNumber ,oldBiopolymerClassBiopolymerType  ,proteinCapping);
            setOriginalSequence(myChain,myOriginalSequence);
            cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": Restoring residue numbers and insertion codes after mutating.. "<<endl;
            updBiopolymerClass(myChain).setResidueIDsAndInsertionCodesFromBiopolymer(tempBiopolymer, proteinCapping);
            updBiopolymerClass(myChain).setMutationVector(myOldBiopolymerClass.getMutationVector());
            if (updBiopolymerClass(myChain).getNumMutationVectorElements() != myOldBiopolymerClass.getNumMutationVectorElements()){
            cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Unexplained error!"<<std::endl;
        }       
}*/

/*vector<Mutation> BiopolymerClassContainer::getCompositeMutationVector() {
    std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" This is obsolete! Just call getMutationVector(). It will return a const vector <Mutation> ."<<std::endl;
    return          mutationVector;
};*/    

    // have  Chromosome::setMutationVectorFromString call this also.
void BiopolymerClassContainer::setMutationVectorFromString (const std::string mutationString) {
    /* // Old way:
    std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Parsing "<<mutationString<<std::endl;
    size_t minorSeparatorPosition1 = mutationString.find(MUTATIONMINORSEPARATOR); // find the start position of residue ID
    size_t minorSeparatorPosition2 = mutationString.find(MUTATIONMINORSEPARATOR, (minorSeparatorPosition1 + 1)); // find the end position of residue ID
    size_t dotPosition = mutationString.find('.'); // find the end position of the mutant record
    if (minorSeparatorPosition2 > dotPosition) {std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Ill-formatted mutation string! Remember to use the format X"<<MUTATIONMINORSEPARATOR<<"NNI"<<MUTATIONMINORSEPARATOR<<"S .. where X is the chain ID, NNI is the residue ID and insertion code (if any) and S is the substituted residue type. Separate the mutants with " <<"."<< mutationString<<std::endl; exit(1);}
    string myChainID = mutationString.substr(0,minorSeparatorPosition1);
    string residueIDString = mutationString.substr(minorSeparatorPosition1 + 1, ( minorSeparatorPosition2 - minorSeparatorPosition1 - 1)); // get the string fragment between minorSeparatorPosition1 and minorSeparatorPosition2.
    string mutatedResidueTypeString = mutationString.substr(minorSeparatorPosition2 + 1, (dotPosition - minorSeparatorPosition2 - 1)); // get the string fragment between  minorSeparatorPosition2 and dotPosition -- this is the mutant residue type.
    std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Parsing "<<mutationString<<std::endl;
    std::cout<<". first chain = >"<<myChainID<<"<"<<std::endl;
    std::cout<<" with residue ID >"<<residueIDString <<"<"<< std::endl;
    std::cout<< " and mutant residue type = >"<<mutatedResidueTypeString<<"< "<<std::endl;
    std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<std::endl;
    if (!(hasChainID(myChainID))) {
            std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Could not find a chain  "<<myChainID<<std::endl; exit(1);
    }
    std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<std::endl;
    ResidueID mutantResidueID = updBiopolymerClass(myChainID).residueID( residueIDString ) ; // convert residueIDString to a real ResidueID
    Mutation myMutation;
    myMutation.setResidue(mutantResidueID);
    myMutation.setChain(myChainID);
    myMutation.setSubstitutedResidueType(mutatedResidueTypeString);
    std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" adding the following mutation to the mutation vector of chain >"<<myChainID<<"< :"<<std::endl;
    myMutation.print();
    addMutationToVector(myMutation );
    std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" getNumMutationVectorElements() "<< getNumMutationVectorElements() << std::endl;
    std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" getFormattedMutationsString() "<< getFormattedMutationsString(MUTATIONMINORSEPARATOR )<<std::endl;
    std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" getFormattedMutationsString(MUTATIONMINORSEPARATOR) >"<< getFormattedMutationsString(MUTATIONMINORSEPARATOR)<<"< "<<std::endl;
    if (dotPosition == 0) {std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Bad format! "<<std::endl; exit(1); }
    if ((dotPosition < (mutationString.length()-1))  &&
            (dotPosition != string::npos)) // dotPosition == string::npos indicates no dot was found
    { // recursively send rest of string to the same setMutationVectorFromString function to take care of remaining mutations.
            std::string newMutationString = mutationString.substr((dotPosition+1), string::npos);
            setMutationVectorFromString(newMutationString);
    }

    */ 


///////////////////// New way, which detects and tolerates FoldX formatted mutation strings:


                std::cout<<__FILE__<<":"<<__LINE__<<" Parsing "<<mutationString<<std::endl;

                size_t dotPosition = mutationString.find(','); // If we find a comma, then this is a SKEMPI formatted mutation string.
                if (dotPosition == std::string::npos){dotPosition = mutationString.find(MUTATIONMAJORSEPARATOR); } // If there was no comma, then we are looking for a dot '.', because this could be a bree
                //der-formatted mutation string. Or it could be a single-substitution mutant in SKEMPI format.

                string mySingleMutationString = mutationString.substr(0,dotPosition);
                std::cout<<__FILE__<<":"<<__LINE__<<std::endl;
                Mutation myMutation;
                std::cout<<__FILE__<<":"<<__LINE__<<std::endl;
                myMutation.setChainSubstitutionFromSingleMutationString(mySingleMutationString); // This method automatically detects whether we are using breeder or SKEMPI formatted mutation string, and 
                // parses accordingly. 
                if (!(hasChainID(myMutation.getChain()))) {
                        std::cout<<__FILE__<<":"<<__LINE__<<" BiopolymerClass does not have a chain  "<<myMutation.getChain()<<std::endl; exit(1);
                }
                std::cout<<__FILE__<<":"<<__LINE__<<" adding the following mutation to the mutation vector of chain >"<<myMutation.getChain()<<"< :"<<std::endl;
                myMutation.print();
                addMutationToVector(myMutation );
                std::cout<<__FILE__<<":"<<__LINE__<<" biopolymerClassContainer.getNumMutationVectorElements() "<< getNumMutationVectorElements() << std::endl;
                std::cout<<__FILE__<<":"<<__LINE__<<" biopolymerClassContainer.getFormattedMutationsString() "<< getFormattedMutationsString(MUTATIONMINORSEPARATOR )<<std::endl;
                std::cout<<__FILE__<<":"<<__LINE__<<" biopolymerClassContainer.getFormattedMutationsString(MUTATIONMINORSEPARATOR) >"<< getFormattedMutationsString(MUTATIONMINORSEPARATOR)<<"< "<<std::endl; 
                if (dotPosition == 0) {std::cout<<__FILE__<<":"<<__LINE__<<" Bad format! "<<std::endl; exit(1); }
                if ((dotPosition < (mutationString.length()-1))  &&
                        (dotPosition != string::npos)) // dotPosition == string::npos indicates no dot was found
                { // recursively send rest of string to the same setMutationVectorFromString function to take care of remaining mutations.
                        std::string newMutationString = mutationString.substr((dotPosition+1), string::npos);
                        setMutationVectorFromString(newMutationString);
                }
/////////////////////




};

void BiopolymerClassContainer::addIntraChainInterfaceResidues(String chain, vector<IncludeAllNonBondAtomsInResidue> & myIncludeAllNonBondAtomsInResidueVector , double radius, SimbodyMatterSubsystem & matter,State & state) {
    ResidueStretchContainer <SingleResidue> myResidueStretchContainer;
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" myResidueStretchContainer.getNumResidueStretches() = "<<myResidueStretchContainer.getNumResidueStretches()<< endl;
    #ifdef USE_OPENMM
    myResidueStretchContainer.addIntraChainInterfaceResidues( radius, chain, *this );
    #endif
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" myResidueStretchContainer.getNumResidueStretches() = "<<myResidueStretchContainer.getNumResidueStretches()<< endl;
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<endl;
    IncludeAllNonBondAtomsInResidue myIncludeAllNonBondAtomsInResidue;
    for (int i = 0 ; i < myResidueStretchContainer.getNumResidueStretches(); i++) {
        //cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<endl;
        myIncludeAllNonBondAtomsInResidue.setChain(   myResidueStretchContainer.getResidueStretch(i).getChain());
        if (myIncludeAllNonBondAtomsInResidue.getChain().compare(chain) != 0){
            ErrorManager::instance <<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": Error! The chain ID in myResidueStretchContainer["<<i<<"] is not "<<chain<<endl; ErrorManager::instance.treatError();
        }
        myIncludeAllNonBondAtomsInResidue.setResidue (  myResidueStretchContainer.getResidueStretch(i).getStartResidue() );
        myIncludeAllNonBondAtomsInResidueVector.push_back(myIncludeAllNonBondAtomsInResidue);
    }
};

#ifdef USE_OPENMM
void BiopolymerClassContainer::createDisulphideBridges(std::ofstream & output) {
    vector <MMBAtomInfo> cysteineAtomInfoVector; cysteineAtomInfoVector.clear();
    loadCysteineAtomInfoVector(cysteineAtomInfoVector);   

    OpenMM::NeighborList neighborList;
    openmmVecType boxSize = openmmVecType(10000,10000,10000);
    vector<openmmVecType> particleList(cysteineAtomInfoVector.size());
    vector<set<int> > exclusions( particleList.size() );
    for (int i = 0; i < cysteineAtomInfoVector.size() ; i++) {
	particleList[i] = cysteineAtomInfoVector[i].position;
    }
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" neighborList size is : "<<neighborList.size()<<endl;
    double         radius = .27;
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<endl;
    computeNeighborListVoxelHash(neighborList, particleList.size() , particleList, exclusions, &boxSize, false, radius  , 0.0);
    for ( int j = 0 ; j < neighborList.size(); j++) {
	    //cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<endl;
	    ResidueID residueID1(cysteineAtomInfoVector[neighborList[j].first].residueID);
	    String chain1(cysteineAtomInfoVector[neighborList[j].first].chain);
	    String atom1(cysteineAtomInfoVector[neighborList[j].first].atomName);
	    ResidueID residueID2(cysteineAtomInfoVector[neighborList[j].second].residueID);
	    String chain2(cysteineAtomInfoVector[neighborList[j].second].chain);
	    String atom2(cysteineAtomInfoVector[neighborList[j].second].atomName);
            if (chain1.compare(chain2) == 0) {
		if ( atom1.compare("SG"  ) != 0){
		    ErrorManager::instance<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": Unexpectedly trying to form disulphide bridge between non-SG  atoms "<<atom1<<endl; 
		    ErrorManager::instance.treatError();
		} 
		if ( atom2.compare("SG"  ) != 0){
		    ErrorManager::instance<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": Unexpectedly trying to form disulphide bridge between non-SG  atoms "<<atom2<<endl; 
		    ErrorManager::instance.treatError();
		} 
                output<<"substituteResidue "<<chain1<<" "<<residueID1.outString()<<" X"<<endl;
                output<<"substituteResidue "<<chain2<<" "<<residueID2.outString()<<" X"<<endl;
                output<<"addRingClosingBond "<<chain1<<" "<<residueID1.outString()<<" SG bond2 "<<residueID2.outString()<<" SG bond2 "<<endl;
                //substituteResidue(chain1, residueID1, String("X"),updBiopolymerClass(chain1).getProteinCapping() );		
                //substituteResidue(chain2, residueID2, String("X"),updBiopolymerClass(chain2).getProteinCapping() );		
                //updBiopolymerClass(chain1).addRingClosingBond( residueID1,  atom1,  String("bond2"),   residueID2,  atom2,String("bond2") , SimTK::BondMobility::Free                  );
	    }
	}
     
    }

void BiopolymerClassContainer::createDisulphideBridges() {
    vector <MMBAtomInfo> cysteineAtomInfoVector; cysteineAtomInfoVector.clear();
    loadCysteineAtomInfoVector(cysteineAtomInfoVector);   

    OpenMM::NeighborList neighborList;
    openmmVecType boxSize = openmmVecType(10000,10000,10000);
    vector<openmmVecType> particleList(cysteineAtomInfoVector.size());
    vector<set<int> > exclusions( particleList.size() );
    for (int i = 0; i < cysteineAtomInfoVector.size() ; i++) {
	particleList[i] = cysteineAtomInfoVector[i].position;
    }
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" neighborList size is : "<<neighborList.size()<<endl;
    double         radius = .27;
    cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<endl;
    computeNeighborListVoxelHash(neighborList, particleList.size() , particleList, exclusions, &boxSize, false, radius  , 0.0);
    for ( int j = 0 ; j < neighborList.size(); j++) {
	    //cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<endl;
	    ResidueID residueID1(cysteineAtomInfoVector[neighborList[j].first].residueID);
	    String chain1(cysteineAtomInfoVector[neighborList[j].first].chain);
	    String atom1(cysteineAtomInfoVector[neighborList[j].first].atomName);
	    ResidueID residueID2(cysteineAtomInfoVector[neighborList[j].second].residueID);
	    String chain2(cysteineAtomInfoVector[neighborList[j].second].chain);
	    String atom2(cysteineAtomInfoVector[neighborList[j].second].atomName);
            if (chain1.compare(chain2) == 0) {
		if ( atom1.compare("SG"  ) != 0){
		    ErrorManager::instance<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": Unexpectedly trying to form disulphide bridge between non-SG  atoms "<<atom1<<endl; 
		    ErrorManager::instance.treatError();
		} 
		if ( atom2.compare("SG"  ) != 0){
		    ErrorManager::instance<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": Unexpectedly trying to form disulphide bridge between non-SG  atoms "<<atom2<<endl; 
		    ErrorManager::instance.treatError();
		} 
                
                substituteResidue(chain1, residueID1, String("X"),updBiopolymerClass(chain1).getProteinCapping() );		
                substituteResidue(chain2, residueID2, String("X"),updBiopolymerClass(chain2).getProteinCapping() );		
                updBiopolymerClass(chain1).addRingClosingBond( residueID1,  atom1,  String("bond2"),   residueID2,  atom2,String("bond2") , SimTK::BondMobility::Free                  );
	    }
	}
     
    }
#endif

void BiopolymerClassContainer::loadCysteineAtomInfoVector(vector <MMBAtomInfo> & cysteineAtomInfoVector ) {
    vector <MMBAtomInfo> myConcatenatedAtomInfoVector;
    myConcatenatedAtomInfoVector = getConcatenatedAtomInfoVector();
    if (myConcatenatedAtomInfoVector.size() ==0 ) {
      ErrorManager::instance<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": initializeAtomInfoVector has not been called!"<<endl; 
      ErrorManager::instance.treatError();
    }
    if (cysteineAtomInfoVector.size() > 0 ) {
      ErrorManager::instance<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<":cysteineAtomInfoVector is not empty!"<<endl; 
      ErrorManager::instance.treatError();
    }
    //cysteineAtomInfoVector.clear();
    for  (int i = 0 ; i < myConcatenatedAtomInfoVector.size(); i++){
        myConcatenatedAtomInfoVector[i].print();
        if(myConcatenatedAtomInfoVector[i].atomName.compare("SG") ==0) {
            cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" Found an SG.."<<endl; 
            //cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" PDB residue name 0 "<<updBiopolymerClass(myConcatenatedAtomInfoVector[i].chain ).updResidueInfo(myConcatenatedAtomInfoVector[i].residueID).getPdbResidueName ()<<endl;
            //if (updBiopolymerClass(myConcatenatedAtomInfoVector[i].chain ).updResidueInfo(myConcatenatedAtomInfoVector[i].residueID).getPdbResidueName () .compare("CYS")){
                cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" adding a cysteine SG .."<<endl;  
                cysteineAtomInfoVector.push_back(myConcatenatedAtomInfoVector[i]);
            //}
        }
    }
    
}

void BiopolymerClass::sort( vector <ResidueID> & residueIDVector){
    // Use an insertion sort
    for (int i = 1; i < (residueIDVector.size() - 1); i++){
        int j = i;
        while ((j > 0) && (getResidueIndex(residueIDVector[j-1]) > getResidueIndex(residueIDVector[j]))){
            // swap elements j and j-1:
            swap(residueIDVector[j], residueIDVector[j-1]);
            //ResidueID tempj =     residueIDVector[j];
            //residueIDVector[j] = residueIDVector[j-1];
            //residueIDVector[j-1] = tempj;

            j--;
        } // of while
    } // of for
}; // of method 
