/* -------------------------------------------------------------------------- *
 *                           MMB (MacroMoleculeBuilder)                       *
 * -------------------------------------------------------------------------- *
 *                                                                            *
 * Copyright (c) 2011-12 by the Author.                                       *
 * Author: Samuel Flores                                                      *
 *                                                                            *
 * See RNABuilder.cpp for the copyright and usage agreement.                  *
 * -------------------------------------------------------------------------- */

#ifndef ContactContainer_H_
#define ContactContainer_H_

#include "BiopolymerClass.h"
#include "ResidueStretchContainer.h"

class MMB_EXPORT ContactContainer : public ResidueStretchContainer <ContactStretch> {
private:
	//vector <ContactStretch> contactStretchVector;
	vector <ContactWithin> contactWithinVector;
    //HuntCrossleyForce huntCrossleyForce;
public:
  	void clear();
  	void validateContact(ContactStretch, BiopolymerClassContainer & myBiopolymerClassContainer);
	void addContactToVector(ContactStretch myContactStretch) ;
	void addContactToVector(ContactStretch myContactStretch, BiopolymerClassContainer & myBiopolymerClassContainer) ;
	void addContactToVector(string myChain, int myStartResidue, int myEndResidue, string myContactScheme, BiopolymerClassContainer & myBiopolymerClassContainer) ;

    void deleteContact(int id);
    void updateContact(int id, string myChain, int myStartResidue, int myEndResidue, string myContactScheme, BiopolymerClassContainer & myBiopolymerClassContainer);

    vector <ContactWithin> & getContactWithinVector() { return contactWithinVector; }
	void validateContactWithin(ContactWithin contactWithin ,  BiopolymerClassContainer & myBiopolymerClassContainer);
	void pushContactWithin(ContactWithin contactWithin, BiopolymerClassContainer & myBiopolymerClassContainer);
	void createContactsWithin( BiopolymerClassContainer & myBiopolymerClassContainer, State & state );
    void deleteContactWithin(int id);
    void updateContactWithin(int id, String chainID, int resID, double radius, String contactScheme, BiopolymerClassContainer & myBiopolymerClassContainer);

	void listDistances( BiopolymerClassContainer & myBiopolymerClassContainer, State & state );


	void printContact(ContactStretch contactStretch);
    void printContact(int contactIndex) ;
    void printContacts() ;
    //void createHuntCrossleyForce(GeneralForceSubsystem & forces,GeneralContactSubsystem & contacts,ContactSetIndex contactSetLargeSpheres)
    //HuntCrossleyForce & getHuntCrossleyForce();  
	void applyContactsToBiopolymers	(BiopolymerClassContainer & myBiopolymerClassContainer,GeneralContactSubsystem &  contacts,GeneralForceSubsystem & forces,SimbodyMatterSubsystem & matter, LeontisWesthofClass & myLeontisWesthofClass,double excludedVolumeRadius, double excludedVolumeStiffness);
	ContactStretch getContact(int contactIndex);
	int numContacts() ;
    bool hasSharedContact(String chainID, ResidueID startResidueID, ResidueID endResidueID,  String contactScheme) ;
    bool hasSharedContact(ContactStretch contactStretch) ;
};

#endif
