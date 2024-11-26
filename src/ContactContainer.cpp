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
#include "BiopolymerClass.h"
#include "ContactContainer.h"
#include "SimTKmolmodel.h"


void ContactContainer::clear(){
    residueStretchVector.clear();
    MMBLOG_FILE_FUNC_LINE(INFO, "Just cleared residueStretchVector .. this now contains "<<numContacts()<<" contacts "<<endl);
    contactWithinVector.clear();

}

void ContactContainer::validateContact(ContactStretch myContactStretch, BiopolymerClassContainer & myBiopolymerClassContainer){
    validateResidueStretch(myContactStretch, myBiopolymerClassContainer);
    bool myEndCapsOn;
    if (
	(
	(!(myBiopolymerClassContainer.updBiopolymerClass((myContactStretch.getChain())).getBiopolymerType()  == BiopolymerType::RNA))
	&&
	(!(myBiopolymerClassContainer.updBiopolymerClass((myContactStretch.getChain())).getBiopolymerType()  == BiopolymerType::DNA))
	) &&
       (((myContactStretch).ContactScheme).compare("SelectedAtoms") == 0))  {
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "The SelectedAtoms sterics can only be applied to RNA or DNA chains.  You are attempting to apply them to a "<< myBiopolymerClassContainer.updBiopolymerClass((myContactStretch.getChain())).getBiopolymerType() <<" chain. "<<endl);
    }


	if ((myBiopolymerClassContainer.updBiopolymerClass((myContactStretch.getChain())).getBiopolymerType()  == BiopolymerType::Protein) && myBiopolymerClassContainer.updBiopolymerClass((myContactStretch.getChain())).getProteinCapping())
	    {    myEndCapsOn = true;}
	else {myEndCapsOn = false;} // importantly, myEndCapsOn is always false for RNA

	if (((myContactStretch).ContactScheme).compare("ExcludedVolume") ==0 ){
	    MMBLOG_FILE_FUNC_LINE(CRITICAL, "ExcludedVolume is no longer a supported contact type."<<endl);
	}




	 else if (((myContactStretch).ContactScheme).compare("RNABackboneSterics") ==0 ){
	     MMBLOG_FILE_FUNC_LINE(CRITICAL, "RNABackboneSterics is no longer supported. "<<endl);

	     if (!(myBiopolymerClassContainer.updBiopolymerClass((myContactStretch.getChain())).getBiopolymerType()  == BiopolymerType::RNA)) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "The RNABackboneSterics can only be applied to RNA chains.  You are attempting to apply them to a non-RNA chain. "<<endl);
		 }
	} else {
		    // everything is OK!
	}
}

void ContactContainer::validateContactWithin(ContactWithin contactWithin ,  BiopolymerClassContainer & myBiopolymerClassContainer){
  	myBiopolymerClassContainer.updBiopolymerClass(contactWithin.Chain).validateResidueID(contactWithin.Residue);
}

void ContactContainer::pushContactWithin ( ContactWithin contactWithin, BiopolymerClassContainer & myBiopolymerClassContainer){
    validateContactWithin(contactWithin,myBiopolymerClassContainer);
    contactWithinVector.push_back(contactWithin);
}

void ContactContainer::deleteContactWithin(int id){
    if(id < 0 || id >= contactWithinVector.size()){
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "you tried to delete a non existing Contact." << endl);
    }
    contactWithinVector.erase(contactWithinVector.begin()+id);
}

void ContactContainer::updateContactWithin(int id, String chainID, int resID, double radius, String contactScheme, BiopolymerClassContainer & myBiopolymerClassContainer){
    if(id < 0 || id >= contactWithinVector.size()){
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "you tried to update a non existing Contact." << endl);
    }
    ContactWithin & newContactWithin = contactWithinVector[id];
    newContactWithin.Chain = chainID;
    newContactWithin.Residue = ResidueID(resID, ' ');
    newContactWithin.Radius =  radius;
    newContactWithin.ContactScheme = contactScheme;

    validateContactWithin(newContactWithin, myBiopolymerClassContainer);
    // contactWithinVector[id] = newContactWithin;
}


void ContactContainer::listDistances ( BiopolymerClassContainer & myBiopolymerClassContainer, State & state ){
    for (int g = 0; g < myBiopolymerClassContainer.getNumBiopolymers(); g++) 
    {
	BiopolymerClass  primaryBiopolymerClass = myBiopolymerClassContainer.updBiopolymerClass(g);
	
	for (ResidueID h = primaryBiopolymerClass.getFirstResidueID(); h < primaryBiopolymerClass.getLastResidueID(); primaryBiopolymerClass.incrementResidueID( h)) {
            
	    for (int i = g+1; i < myBiopolymerClassContainer.getNumBiopolymers(); i++) { // start from h+1 to prevent double counting and calculating distances to self.
		BiopolymerClass  partnerBiopolymerClass = myBiopolymerClassContainer.updBiopolymerClass(i);

		for (ResidueID j = partnerBiopolymerClass.getFirstResidueID(); j <= partnerBiopolymerClass.getLastResidueID();   partnerBiopolymerClass.incrementResidueID( j)) {
		    double myDistance = (double)(
			partnerBiopolymerClass.calcAtomLocationInGroundFrame  (state, j, partnerBiopolymerClass.getRepresentativeAtomName())
			- primaryBiopolymerClass.calcAtomLocationInGroundFrame(state, h, primaryBiopolymerClass.getRepresentativeAtomName())
			).norm();
		    MMBLOG_FILE_FUNC_LINE(INFO, "Distance between chain "<<primaryBiopolymerClass.getChainID()<<" , residue "<<h.outString()<<" , atom "<<primaryBiopolymerClass.getRepresentativeAtomName()<<" and chain "<<partnerBiopolymerClass.getChainID()<<" residue "<<j.outString()<<" atom "<<partnerBiopolymerClass.getRepresentativeAtomName() <<" = "<<myDistance<<endl);
		} // of for j
	    } // of for i
      } // of for h
    } // of for g
} // of method

void ContactContainer::createContactsWithin ( BiopolymerClassContainer & myBiopolymerClassContainer, State & state ){
    double maxRadius = 0.0;
    for (int h = 0; h < (int)contactWithinVector.size(); h++)
    {
        if(contactWithinVector[h].Radius > maxRadius)
            maxRadius = contactWithinVector[h].Radius;
    }
    if(maxRadius <= 0)
    {
        return;
    }
    vector<MMBAtomInfo> concatenatedAtomInfoVector = myBiopolymerClassContainer.getConcatenatedAtomInfoVector(state,true);
    MMBLOG_FILE_FUNC_LINE(INFO, "ContactContainer: caching neighborList"<<endl);
    OpenMM::NeighborList neighborList = myBiopolymerClassContainer.getNeighborList(concatenatedAtomInfoVector, maxRadius); 
    myBiopolymerClassContainer.setNeighborsFromList(concatenatedAtomInfoVector, neighborList, maxRadius);
    for (int h = 0; h < (int)contactWithinVector.size(); h++)
    {
        MMBLOG_FILE_FUNC_LINE(DEBUG, "Starting to apply contacts within "<<contactWithinVector[h].Radius<<" nm of chain "<<contactWithinVector[h].Chain<<" residue "<<contactWithinVector[h].Residue.outString()<<endl);
        MMBLOG_FILE_FUNC_LINE(DEBUG, "Note that in MMB 2.10 and earlier, we took this radius in Å.  We are going back to nm for consistency, with apologies for the confusion."<<endl);
        // vector< pair<const BiopolymerClass*, const ResidueID*> > residuesWithin = myBiopolymerClassContainer.getResiduesWithin(contactWithinVector[h].Chain, contactWithinVector[h].Residue, contactWithinVector[h].Radius, state, &neighborList);

        // vector< pair<const BiopolymerClass*, const ResidueID*> >::iterator it;
        // for(it = residuesWithin.begin(); it != residuesWithin.end(); it++){


        vector<MMBAtomInfo>::iterator itAtom;
        for(itAtom = concatenatedAtomInfoVector.begin(); itAtom != concatenatedAtomInfoVector.end(); itAtom++)
        {
            if(itAtom->getChain() == contactWithinVector[h].Chain && itAtom->getResidueID() == contactWithinVector[h].Residue)
            {
                vector<MMBAtomInfo*>::iterator itAtomPointer;
                for(itAtomPointer=itAtom->neighbors.begin(); itAtomPointer!=itAtom->neighbors.end(); itAtomPointer++)
                {
                    BiopolymerClass & myBiopolymerClass = myBiopolymerClassContainer.updBiopolymerClass((*itAtomPointer)->getChain());
                    // Skip residues in unactive chains
                    if(myBiopolymerClass.getActivePhysics() == false)
                        continue;
                    MMBLOG_FILE_FUNC_LINE(DEBUG, "distance is less than "<<contactWithinVector[h].Radius<<", checking if contact already exists : ");
                    ContactStretch myContact;
                    myContact.setChain (myBiopolymerClass.getChainID());
                    myContact.setStartResidue((*itAtomPointer)->getResidueID());
                    myContact.setEndResidue ((*itAtomPointer)->getResidueID());
                    myContact.ContactScheme = contactWithinVector[h].ContactScheme;
                    MMBLOG_FILE_FUNC_LINE(DEBUG, "result : "<<hasSharedContact(myContact)<<endl);
                    if (! (hasSharedContact(myContact))) {
                        MMBLOG_FILE_FUNC_LINE(DEBUG, "Adding contact : ");
                        printContact(myContact);
                        addContactToVector(myContact, myBiopolymerClassContainer);
                    } else if (hasSharedContact(myContact) ) {
                        MMBLOG_FILE_FUNC_LINE(DEBUG, "Not adding contact : ");
                        printContact(myContact);
                    } else {  
                        MMBLOG_FILE_FUNC_LINE(CRITICAL, "Unexplained error!"<<endl);
                    }
                    MMBLOG_FILE_FUNC_LINE(DEBUG, endl);//" result : "<<hasSharedContact(myContact)<<endl;
                }
            }
        }
    }
    MMBLOG_FILE_FUNC_LINE(INFO, "done with createContactsWithin .. at this stage, the list of contacts is:"<<endl);
    //MMBLOG_FILE_FUNC_LINE(" Preparing to print all "<<numContacts()<<" contacts "<<endl;
    printContacts();
    MMBLOG_FILE_FUNC_LINE(INFO, "done with createContactsWithin .. moving on .."<<endl);
} // of method

void ContactContainer::printContact(ContactStretch contactStretch)  {
        MMBLOG_FILE_FUNC_LINE(INFO, "Contact scheme = "<<contactStretch.ContactScheme<<" chain= "<<contactStretch.getChain()<<" from residue "<<contactStretch.getStartResidue().outString()<<" to "<<contactStretch.getEndResidue().outString()<<endl);
}
void ContactContainer::printContact(int contactIndex) {
	ContactStretch contactStretch = getResidueStretch(contactIndex );
          
        MMBLOG_FILE_FUNC_LINE(INFO, "Contact "<<contactIndex<<" ");//scheme = "<<contactStretch.ContactScheme<<" chain= "<<contactStretch.getChain()<<" from residue "<<contactStretch.getStartResidue().outString()<<" to "<<contactStretch.getEndResidue().outString()<<endl;
        printContact(contactStretch);
}


void ContactContainer::printContacts(){
    MMBLOG_FILE_FUNC_LINE(INFO, "Preparing to print all "<<numContacts()<<" contacts "<<endl);
    for (int i = 0; i < numContacts(); i++) {
        printContact(i);
    }
}

/*void ContactContainer::createHuntCrossleyForce(forces, contacts, contactSetLargeSpheres) {//  HuntCrossleyForce & myHuntCrossleyForce ){
    huntCrossleyForce = HuntCrossleyForce (forces, contacts, contactSetLargeSpheres ) ;
};
HuntCrossleyForce & ContactContainer::getHuntCrossleyForce(){
    return huntCrossleyForce;
};*/

    
/*void ContactContainer::createHuntCrossleyForce(forces, contacts, contactSetLargeSpheres) {//  HuntCrossleyForce & myHuntCrossleyForce ){
    huntCrossleyForce = HuntCrossleyForce (forces, contacts, contactSetLargeSpheres ) ;
};
HuntCrossleyForce & ContactContainer::getHuntCrossleyForce(){
    return huntCrossleyForce;
};*/

void ContactContainer::applyContactsToBiopolymers(BiopolymerClassContainer & myBiopolymerClassContainer,GeneralContactSubsystem &  contacts,GeneralForceSubsystem & forces,SimbodyMatterSubsystem &  matter, LeontisWesthofClass & myLeontisWesthofClass,double excludedVolumeRadius, double excludedVolumeStiffness ){
ContactSetIndex contactSetLargeSpheres = contacts.createContactSet();

HuntCrossleyForce hcLargeSpheres(forces, contacts, contactSetLargeSpheres);
for (int q=0;q<numContacts();q++)  
{    
    ContactStretch myContactStretch = getResidueStretch(q);  
    bool myEndCapsOn;
    if ( 
	(    
	(!(myBiopolymerClassContainer .updBiopolymerClass((myContactStretch.getChain())).getBiopolymerType()  == BiopolymerType::RNA))                &&
	(!(myBiopolymerClassContainer .updBiopolymerClass((myContactStretch.getChain())).getBiopolymerType()  == BiopolymerType::DNA))
	) &&                (((myContactStretch).ContactScheme).compare("SelectedAtoms") == 0))  {
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "The SelectedAtoms sterics can only be applied to RNA or DNA chains. You are attempting to apply them to a "<< myBiopolymerClassContainer .updBiopolymerClass((myContactStretch.getChain())).getBiopolymerType() <<" chain. "<<endl);
    };


    if ((myBiopolymerClassContainer .updBiopolymerClass((myContactStretch.getChain())).getBiopolymerType()  == BiopolymerType::Protein) && myBiopolymerClassContainer .updBiopolymerClass((myContactStretch.getChain())).getProteinCapping())
	{    myEndCapsOn = true;}
    else {myEndCapsOn = false;} // importantly, myEndCapsOn is always false for RNA

    if ((myContactStretch.ContactScheme).compare("ExcludedVolume") ==0 ){
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "ExcludedVolume is no longer a supported contact type."<<endl);
    }

     else if ((myContactStretch.ContactScheme).compare("RNABackboneSterics") ==0 ){
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "RNABackboneSterics is no longer supported. "<<endl);//s.  You are attempting to apply them to a non-RNA chain. "<<endl;
     }
     else if ((myContactStretch.ContactScheme).compare("AllHeavyAtomSterics") ==0 )
	myBiopolymerClassContainer.updBiopolymerClass(myContactStretch.getChain()).addGeneralSterics(contacts,contactSetLargeSpheres ,hcLargeSpheres, matter,excludedVolumeRadius, excludedVolumeStiffness , (myContactStretch).getStartResidue(),(myContactStretch).getEndResidue(), myEndCapsOn, false );
    else if (((myContactStretch).ContactScheme).compare("AllAtomSterics") ==0 ){
	myBiopolymerClassContainer.updBiopolymerClass(myContactStretch.getChain()).addGeneralSterics(contacts,contactSetLargeSpheres ,hcLargeSpheres, matter, excludedVolumeRadius, excludedVolumeStiffness  , (myContactStretch).getStartResidue(),(myContactStretch).getEndResidue(),myEndCapsOn, true );
    }

    else  { // any contact types defined in the parameter file are applied here. 
	    myBiopolymerClassContainer.updBiopolymerClass (myContactStretch.getChain() ).addCustomSterics(contacts,contactSetLargeSpheres, hcLargeSpheres, matter,myLeontisWesthofClass,(myContactStretch).ContactScheme, (myContactStretch).getStartResidue(),(myContactStretch).getEndResidue(),myEndCapsOn);
    }   

} // of for q  
    MMBLOG_FILE_FUNC_LINE(INFO, "Have applied "<<contacts.getNumBodies(contactSetLargeSpheres)<<" contact spheres to atoms."<<endl);
} // of method




void ContactContainer::addContactToVector(ContactStretch myContactStretch, BiopolymerClassContainer & myBiopolymerClassContainer) {
    validateContact(myContactStretch, myBiopolymerClassContainer); 
    ResidueStretchContainer<ContactStretch>::addStretch  (myContactStretch );
    printContact(myContactStretch);
 //    getNumResidueStretches();
 //    MMBLOG_FILE_FUNC_LINE(" Just added a contact to the residue stretch vector.  Now have "<<numContacts()<<" contacts "<<endl;
 //    int myInt = getNumResidueStretches();
 //    myInt --;

	// ContactStretch contactStretch = getResidueStretch(myInt );
 //        MMBLOG_FILE_FUNC_LINE( " from residue "<<contactStretch.getStartResidue().outString()<<" to "<<contactStretch.getEndResidue().outString()<<endl;
 //        MMBLOG_FILE_FUNC_LINE(   " Contact "<<myInt<<" scheme = "<<contactStretch.ContactScheme;//<<" chain= "<<contactStretch.getChain()<<" from residue "<<contactStretch.getStartResidue().outString()<<" to "<<contactStretch.getEndResidue().outString()<<endl;
 //        printContact(contactStretch);

 //    printContact(myInt);
    
}

void ContactContainer::deleteContact(int id){
    if(id < 0 || id >= residueStretchVector.size()){
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "you tried to delete a non existing Contact." << endl);
    }
    residueStretchVector.erase(residueStretchVector.begin()+id);
}
void ContactContainer::updateContact(int id, string myChain, int myStartResidue, 
                                     int myEndResidue, string myContactScheme, 
                                     BiopolymerClassContainer & myBiopolymerClassContainer){
    if(id < 0 || id >= residueStretchVector.size()){
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "you tried to update a non existing Contact." << endl);
    }
    ContactStretch newStretch;
    newStretch.setChain(myChain);
    newStretch.setStartResidue(ResidueID(myStartResidue, ' '));
    newStretch.setEndResidue(ResidueID(myEndResidue, ' '));
    newStretch.ContactScheme = myContactScheme;
    validateContact(newStretch, myBiopolymerClassContainer);
    residueStretchVector[id] = newStretch;
}

ContactStretch ContactContainer::getContact(int myContactIndex ) {
    printContact(ResidueStretchContainer<ContactStretch>::getResidueStretch(myContactIndex)); // this printed out fine!
    return ResidueStretchContainer<ContactStretch>::getResidueStretch(myContactIndex);
}


int ContactContainer::numContacts(){
    //MMBLOG_FILE_FUNC_LINE(" num contacts = "<<getNumResidueStretches()<<endl;
    return getNumResidueStretches();
    //return residueStretchVector.size();           
}

bool  ContactContainer::hasSharedContact(String chainID, ResidueID startResidueID, ResidueID endResidueID, String contactScheme )  {
    for (int i = 0; i < numContacts(); i++){
        ContactStretch tempContactStretch = getResidueStretch(i);
        if (
            (tempContactStretch.ContactScheme.compare(contactScheme ) == 0 ) && 
            (tempContactStretch.getChain().compare(chainID ) == 0 ) && 
            ((tempContactStretch.getStartResidue() <= endResidueID)   && (tempContactStretch.getEndResidue() >= endResidueID) ||   // the supplied end residue is in the range of this ContactStretch (with index $i).
             (tempContactStretch.getStartResidue() <= startResidueID) && (tempContactStretch.getEndResidue() >= startResidueID)) // the supplied start residue is in the range of this ContactStretch (with index $i).
            ) {
            MMBLOG_FILE_FUNC_LINE(DEBUG, "The stretch of residues you proposed (chain "<<chainID<<" residues "<<startResidueID.outString()<<" to "<< endResidueID.outString()<<" overlaps with an existing stretch, from "<<tempContactStretch.getStartResidue().outString()<<" to "<<tempContactStretch.getEndResidue().outString()<<". "<<endl);
            return true  ;
            MMBLOG_FILE_FUNC_LINE(CRITICAL, " Unexplained error!"<<endl);
        }
    }
    MMBLOG_FILE_FUNC_LINE(DEBUG, "The stretch of residues you proposed (chain "<<chainID<<" residues "<<startResidueID.outString()<<" to "<< endResidueID.outString()<<" does not overlap with any existing stretch"<<endl);
    return (false);
    MMBLOG_FILE_FUNC_LINE(CRITICAL, "Unexplained error!"<<endl);
}

bool  ContactContainer::hasSharedContact(ContactStretch contactStretch)  {
    return hasSharedContact(contactStretch.getChain(), contactStretch.getStartResidue(), contactStretch.getEndResidue(), contactStretch.ContactScheme);
}
