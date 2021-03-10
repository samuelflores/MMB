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
#include "MobilizerContainer.h"
#include "SimTKmolmodel.h"
#include "ReferenceNeighborList.h"


void MobilizerContainer::clear(){
    residueStretchVector.clear();
    mobilizerWithinVector.clear();
    interfaceContainer.clear();
};
bool MobilizerContainer::isEmpty() {
    if ((residueStretchVector.size() >0) ||
        (mobilizerWithinVector.size() >0) ||
        (interfaceContainer.numInterfaces() >0) )
        {return false;}
    else {return true;}
}



void MobilizerContainer::validateMobilizerStretch(MobilizerStretch & myMobilizerStretch, BiopolymerClassContainer & myBiopolymerClassContainer){
    
    if (!(myBiopolymerClassContainer.hasChainID(myMobilizerStretch.getChain()))){
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "Couldn't find chain "<<myMobilizerStretch.getChain()<<endl);
    }
    if ( myBiopolymerClassContainer.updBiopolymerClass(myMobilizerStretch.getChain()).difference(myMobilizerStretch.getEndResidue() , myMobilizerStretch.getStartResidue()) < 0) {
    }
    myBiopolymerClassContainer.updBiopolymerClass(myMobilizerStretch.getChain()).validateResidueID(myMobilizerStretch.getEndResidue());
    if ((myMobilizerStretch.getStartResidue() < myBiopolymerClassContainer.updBiopolymerClass(myMobilizerStretch.getChain()).getFirstResidueID()) ) {
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "The start residue (currently "<<myMobilizerStretch.getStartResidue().outString()<<") is lesser than the first residue number of the chain."<<endl);
    }
    validateResidueStretch(myMobilizerStretch,myBiopolymerClassContainer); // This method from the parent class has basic checks, e.g. making sure EndResidue > StartResidue.
};

void MobilizerContainer::validateMobilizerStretch(int  mobilizerStretchIndex, BiopolymerClassContainer & myBiopolymerClassContainer){
    MobilizerStretch myMobilizerStretch = getResidueStretch(mobilizerStretchIndex);
    validateMobilizerStretch(myMobilizerStretch,myBiopolymerClassContainer);
};

void MobilizerContainer::printMobilizerStretch(int mobilizerStretchIndex){
    MMBLOG_FILE_FUNC_LINE(INFO, "Mobilizer stretch "<<mobilizerStretchIndex<<" BondMobility = "<<getResidueStretch(mobilizerStretchIndex ).getBondMobilityString()<<endl);
    MMBLOG_FILE_FUNC_LINE(INFO, "chain= "<<getResidueStretch(mobilizerStretchIndex ).getChain()  <<" from residue "<<getResidueStretch(mobilizerStretchIndex ).getStartResidue().outString()<<" to "<<getResidueStretch(mobilizerStretchIndex ).getEndResidue().outString()<<endl);
};

void MobilizerContainer::printMobilizerStretches(){
    for (int i = 0; i < getNumResidueStretches(); i++) {
        printMobilizerStretch(i);
    }
};
    




void MobilizerContainer::addMobilizerStretchToVector(MobilizerStretch myMobilizerStretch, BiopolymerClassContainer & myBiopolymerClassContainer) {
    validateMobilizerStretch(myMobilizerStretch, myBiopolymerClassContainer); 
    MMBLOG_FILE_FUNC_LINE(INFO, "Adding mobilizer stretch to vector:"<<endl);
    myMobilizerStretch.print();
    addStretch(myMobilizerStretch );
    //printMobilizerStretch(getNumResidueStretches()-1);
};

void MobilizerContainer::addMobilizerStretchToVector(String myChain, ResidueID myStartResidue, ResidueID myEndResidue, String bondMobilityString, BiopolymerClassContainer & myBiopolymerClassContainer) {
    //cout<<__FILE__<<":"<<__LINE__<<endl;
    MobilizerStretch myMobilizerStretch;
    //cout<<__FILE__<<":"<<__LINE__<<endl;
    myMobilizerStretch.setChain ( myChain);
    myMobilizerStretch.setStartResidue ( myStartResidue);
    myMobilizerStretch.setEndResidue ( myEndResidue);
    //cout<<__FILE__<<":"<<__LINE__<<endl;
    myMobilizerStretch.setBondMobility(bondMobilityString );
    //cout<<__FILE__<<":"<<__LINE__<<" bondMobilityString = >"<<bondMobilityString<<"<"<<endl;
    //cout<<__FILE__<<":"<<__LINE__<<endl;
    myMobilizerStretch.print();
    addMobilizerStretchToVector(myMobilizerStretch,myBiopolymerClassContainer);
};


void MobilizerContainer::addMobilizerStretchToVector(String myChain, String bondMobilityString, BiopolymerClassContainer & myBiopolymerClassContainer) {
    
    addMobilizerStretchToVector(myChain, myBiopolymerClassContainer.updBiopolymerClass(myChain).getFirstResidueID(), myBiopolymerClassContainer.updBiopolymerClass(myChain).getLastResidueID(),bondMobilityString, myBiopolymerClassContainer);
};

void MobilizerContainer::addMobilizerStretchToVector(MobilizerStretch residueStretch, String  bondMobilityString , BiopolymerClassContainer & myBiopolymerClassContainer) {
        addMobilizerStretchToVector(MobilizerStretch(residueStretch, bondMobilityString) ,  myBiopolymerClassContainer);  
};

void MobilizerContainer::updateMobilizerStretch(int id, String myChain, 
                                                ResidueID myStartResidue, 
                                                ResidueID myEndResidue, 
                                                String bondMobilityString, 
                                                BiopolymerClassContainer & myBiopolymerClassContainer){
    MobilizerStretch & ms = residueStretchVector[id];
    ms.setChain(myChain);
    ms.setStartResidue(myStartResidue);
    ms.setEndResidue(myEndResidue);
    ms.setBondMobility(bondMobilityString);
    validateMobilizerStretch(ms, myBiopolymerClassContainer);
}

void MobilizerContainer::deleteMobilizerStretch(int id){
    residueStretchVector.erase(residueStretchVector.begin()+id);
}

/*
void MobilizerContainer::addMobilizerStretchesToVector(vector <MobilizerStretch> residueStretchVector, String  bondMobilityString , BiopolymerClassContainer & myBiopolymerClassContainer) {
    for (int i = 0; i < residueStretchVector.size()  ; i++) {
        cout<<__FILE__<<":"<<__LINE__<<" Inside addMobilizerStretchesToVector, adding stretch "<<i<<endl;
        cout<<__FILE__<<":"<<__LINE__<<" index "<<i<<" chain "<<residueStretchVector[i].getChain()<<", residue = "<<residueStretchVector[i].getStartResidue().outString()<<" to "<<residueStretchVector[i].getEndResidue().outString()<<endl;
        addMobilizerStretchToVector(residueStretchVector[i], bondMobilityString,  myBiopolymerClassContainer);
    }
};

// loop through all specified Interfaces, add interface residues to residueStretchVector

void MobilizerContainer::addMobilizerStretchesToVector(BiopolymerClassContainer  & myBiopolymerClassContainer) {
	for (int i = 0; i < interfaceContainer.numInterfaces(); i++) {
                // create temporary ResidueStretchContainer which will figure out which residues are in the interface.  We then  add these to residueStretchVector.
		ResidueStretchContainer <MobilizerStretch> tempResidueStretchContainer;
		Interface tempMobilizerInterface = interfaceContainer.getInterface(i);
		tempResidueStretchContainer.addAllMutualChainResidues(tempMobilizerInterface.Depth, tempMobilizerInterface.Chains, tempMobilizerInterface.PartnerChains, myBiopolymerClassContainer);
                cout<<__FILE__<<":"<<__LINE__<<endl;
		addMobilizerStretchesToVector(tempResidueStretchContainer.getResidueStretchVector(), tempMobilizerInterface.MobilizerString, myBiopolymerClassContainer);   
                cout<<__FILE__<<":"<<__LINE__<<endl;
	}
};



void MobilizerContainer::addMobilizerStretchesToVector(vector <MobilizerStretch> residueStretchVector,  BiopolymerClassContainer & myBiopolymerClassContainer) {
    for (int i = 0; i < (int)residueStretchVector.size(); i++) {
        addMobilizerStretchToVector(residueStretchVector[i],  myBiopolymerClassContainer);
    }
};*/

/*void MobilizerContainer::addMobilizerStretchesToVector(vector <MobilizerStretch> residueStretchVector,  BiopolymerClassContainer & myBiopolymerClassContainer) {
    for (int i = 0; i < (int)residueStretchVector.size(); i++) {
        addMobilizerStretchToVector(residueStretchVector[i],  state, compoundSystem);
    }
};*/

#ifdef USE_OPENMM
void MobilizerContainer::addMobilizerDomainsInterfacesToVector(const vector<MobilizerDomainsInterface> & mDIVector, BiopolymerClassContainer & myBiopolymerClassContainer)
{
    vector<MMBAtomInfo> atomInfoVector = myBiopolymerClassContainer.getConcatenatedAtomInfoVector();

    // Find the maximum range to compute only one neighbor list
    double maxRange = 0;
    vector<MobilizerDomainsInterface>::const_iterator it;
    for(it=mDIVector.begin(); it!=mDIVector.end(); it++)
    {
        if(it->range > maxRange)
            maxRange = it->range;
    }

    if(maxRange <= 0)
        return;

    vector<MobilizerStretch> addedResidues;
    vector<MobilizerStretch>::iterator msIt;

    OpenMM::NeighborList neighborList = myBiopolymerClassContainer.getNeighborList(atomInfoVector, maxRange);
    // Go through the list
    for ( size_t j = 0 ; j < neighborList.size(); j++)
    {
        unsigned int id1 = neighborList[j].first;
        unsigned int id2 = neighborList[j].second;

        String chain1 = atomInfoVector[id1].chain;
        ResidueID res1 = atomInfoVector[id1].residueID;
        String chain2 = atomInfoVector[id2].chain;
        ResidueID res2 = atomInfoVector[id2].residueID;

        double dist = atomInfoVector[id1].distance(atomInfoVector[id2]);

        // Check for all the MobilizerDomainsInterface
        for(it=mDIVector.begin(); it!=mDIVector.end(); it++)
        {
            // check wether neighbors are across the interface and within the requested range
            if( (dist <= it->range) &&
                (it->domain1.contains(chain1, res1) && it->domain2.contains(chain2, res2)) ||
                (it->domain2.contains(chain1, res1) && it->domain1.contains(chain2, res2)))
            {
                MobilizerStretch rs1 = MobilizerStretch(chain1, res1, res1, it->MobilizerString);
                msIt = find(addedResidues.begin(), addedResidues.end(), rs1);
                if(msIt == addedResidues.end())
                {
                    addMobilizerStretchToVector(rs1, myBiopolymerClassContainer);
                    addedResidues.push_back(residueStretchVector.back());
                    if(it->rigidBackbone)
                        addPhiPsiMobility(chain1, res1, res1, "Rigid", myBiopolymerClassContainer);
                }
                MobilizerStretch rs2 = MobilizerStretch(chain2, res2, res2, it->MobilizerString);
                msIt = find(addedResidues.begin(), addedResidues.end(), rs2);
                if(msIt == addedResidues.end())
                {
                    addMobilizerStretchToVector(rs2, myBiopolymerClassContainer);
                    addedResidues.push_back(residueStretchVector.back());
                    if(it->rigidBackbone)
                        addPhiPsiMobility(chain2, res2, res2, "Rigid", myBiopolymerClassContainer);
                }
            }
        }

    }

    for(msIt=addedResidues.begin(); msIt!=addedResidues.end(); msIt++)
        cerr << msIt->getChain() << " " << msIt->getStartResidue().outString() << endl;
}
#endif

String MobilizerContainer::getChain(int mobilizerStretchIndex){
    return getResidueStretch(mobilizerStretchIndex).getChain();
};
ResidueID MobilizerContainer::getStartResidue(int mobilizerStretchIndex){
    return getResidueStretch(mobilizerStretchIndex).getStartResidue();
};

ResidueID MobilizerContainer::getEndResidue(int mobilizerStretchIndex){
    return getResidueStretch(mobilizerStretchIndex).getEndResidue();
};


void MobilizerContainer::setBiopolymerBondMobility (BiopolymerClassContainer & myBiopolymerClassContainer) {
    for (int q=0;q<getNumResidueStretches();q++)
    {               
        //BondMobility myBondMobility = getBondMobility(q);
        BiopolymerClass & myBiopolymerClass ( myBiopolymerClassContainer.updBiopolymerClass((residueStretchVector[q]).getChain()));
        BiopolymerType::BiopolymerTypeEnum btype = myBiopolymerClassContainer.updBiopolymerClass(getChain(q)).biopolymerType; 

        if (btype == BiopolymerType::RNA){
            (static_cast<RNA&>( myBiopolymerClass.myBiopolymer)).setRNABondMobility(getResidueStretch(q).getBondMobility(),
            SimTK::ResidueInfo::Index (myBiopolymerClass.getResidueIndex(getStartResidue(q) )),
            SimTK::ResidueInfo::Index ( myBiopolymerClass.getResidueIndex( getEndResidue(q))));                   

        } else if (btype == BiopolymerType::DNA){
            (static_cast<DNA&>( myBiopolymerClass.myBiopolymer)).setDNABondMobility(getResidueStretch(q).getBondMobility(),
            SimTK::ResidueInfo::Index (myBiopolymerClass.getResidueIndex  ( getStartResidue(q) )),
            SimTK::ResidueInfo::Index (myBiopolymerClass.getResidueIndex   ( getEndResidue(q) )));

        } else if (btype == BiopolymerType::Protein) {
            myBiopolymerClass.setProteinBondMobility(
                getResidueStretch(q).getBondMobility(),
                getStartResidue(q),
                getEndResidue(q)
            );
        } 
        else {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "biopolymerType " << btype << " unknown" << endl);
        }
    }
};

#ifdef USE_OPENMM
void MobilizerContainer::createMobilizersWithin ( BiopolymerClassContainer & myBiopolymerClassContainer, State & state ){
    MMBLOG_FILE_FUNC_LINE(DEBUG, " calling f for mobilizerWithinVector  of size  "<<mobilizerWithinVector.size()<< endl);
    vector <MobilizerWithin> tempSingleMobilizerWithin; 
    for (size_t h = 0 ; h < mobilizerWithinVector.size() ; h++) { 
         //  it is possible to pass the entire vector and let findBiopolymerResiduesWithinRadius loop over it. however we might then be unable to control the mobilizer types separately
        tempSingleMobilizerWithin.clear(); tempSingleMobilizerWithin.push_back(mobilizerWithinVector[h]);
        MMBLOG_FILE_FUNC_LINE(DEBUG, " calling f for mobilizerWithinVector  of size  "<<mobilizerWithinVector.size()<< endl);
        vector<SingleResidue> myMobilizerResidueVector = myBiopolymerClassContainer.findBiopolymerResiduesWithinRadius(tempSingleMobilizerWithin ,state);
        for (size_t i = 0 ; i < myMobilizerResidueVector.size() ; i ++) 
         {
            MobilizerStretch myMobilizer;// =myMobilizerResidueVector[i];
            myMobilizer.setChain (myMobilizerResidueVector[i].getChain());
            myMobilizer.setStartResidue ( myMobilizerResidueVector[i].getResidue() );
            myMobilizer.setEndResidue (myMobilizerResidueVector[i].getResidue() );
            myMobilizer.setBondMobility(tempSingleMobilizerWithin[0].getBondMobilityString() );
            addMobilizerStretchToVector(myMobilizer, myBiopolymerClassContainer);
            MMBLOG_FILE_FUNC_LINE(DEBUG, endl);
            myMobilizer.printStretch();
         } 
    } // of for h
}; // of method
#endif

void MobilizerContainer::pushMobilizerWithin ( MobilizerWithin mobilizerWithin, BiopolymerClassContainer & myBiopolymerClassContainer){
    validateMobilizerWithin(mobilizerWithin,myBiopolymerClassContainer);
    mobilizerWithinVector.push_back(mobilizerWithin);
};


void MobilizerContainer::validateMobilizerWithin(MobilizerWithin mobilizerWithin ,  BiopolymerClassContainer & myBiopolymerClassContainer){
        myBiopolymerClassContainer.updBiopolymerClass(mobilizerWithin.getChain()).validateResidueID(mobilizerWithin.getResidue()); // not really necessary; this was validated in ParameterReader.cpp

};

int MobilizerContainer::numMobilizerWithin() {
    return mobilizerWithinVector.size() ;
};

void MobilizerContainer::addPhiPsiMobility(String chain, ResidueID startResidue, ResidueID endResidue, String bondMobilityString , BiopolymerClassContainer& myBiopolymerClassContainer ) {
    SingleBondMobility mySingleBondMobility;
    mySingleBondMobility.residue1 = startResidue;	
    mySingleBondMobility.residue2 = startResidue;	
    myBiopolymerClassContainer.validateChainID(chain);
    myBiopolymerClassContainer.updBiopolymerClass(chain).validateResidueID(startResidue);   
    myBiopolymerClassContainer.updBiopolymerClass(chain).validateResidueID(endResidue);   
    mySingleBondMobility.chain1 = chain;
    mySingleBondMobility.chain2 = chain;
    mySingleBondMobility.mobility = bondMobilityString;
    while(  mySingleBondMobility.residue1 <=  endResidue) {
                 // First do the N-CA bond:
                mySingleBondMobility.atom1    = String("N"); 
                mySingleBondMobility.atom2    = String("CA"); 
                // atomPathString used to self-validate.  However it no longer does that, because occasionally we want to create paths that will be validated only later.  In any event, we need to validate explicitly here:
                Compound::AtomPathName myAtomPathName1 = myBiopolymerClassContainer.updBiopolymerClass(mySingleBondMobility.chain1).atomPathString(mySingleBondMobility.residue1, mySingleBondMobility.atom1);
                myBiopolymerClassContainer.updBiopolymerClass(mySingleBondMobility.chain1).validateAtomPathName(myAtomPathName1);
                Compound::AtomPathName myAtomPathName2 = myBiopolymerClassContainer.updBiopolymerClass(mySingleBondMobility.chain2).atomPathString(mySingleBondMobility.residue2, mySingleBondMobility.atom2);
                myBiopolymerClassContainer.updBiopolymerClass(mySingleBondMobility.chain2).validateAtomPathName(myAtomPathName2);
                singleBondMobilityVector.push_back(mySingleBondMobility);

                // Now do the C-CA bond:
                // this might become necessary:
                SingleBondMobility mySingleBondMobility2 = mySingleBondMobility ;
                mySingleBondMobility2.atom1    = String("C"); 
                myAtomPathName1 = myBiopolymerClassContainer.updBiopolymerClass(mySingleBondMobility2.chain1).atomPathString(mySingleBondMobility2.residue1, mySingleBondMobility2.atom1);
                myBiopolymerClassContainer.updBiopolymerClass(mySingleBondMobility2.chain1).validateAtomPathName(myAtomPathName1);
                singleBondMobilityVector.push_back(mySingleBondMobility2);

                if (mySingleBondMobility2.residue1 < endResidue){ // We have to be careful not to increment myResidueID past the end of the chain
                        mySingleBondMobility.residue1 = myBiopolymerClassContainer.updBiopolymerClass(mySingleBondMobility.chain1).incrementResidueID(mySingleBondMobility.residue1);
                        mySingleBondMobility.residue2 = myBiopolymerClassContainer.updBiopolymerClass(mySingleBondMobility.chain2).incrementResidueID(mySingleBondMobility.residue2); 
			if (!(mySingleBondMobility.residue1 == mySingleBondMobility.residue2)) {
			    MMBLOG_FILE_FUNC_LINE(CRITICAL, "residue1 and residue2 are not equal ."<<endl);}
		}
                else if (mySingleBondMobility.residue1  == endResidue) {
                        break;
                }
    }

}

void MobilizerContainer::addPhiPsiMobility(String bondMobilityString , BiopolymerClassContainer& myBiopolymerClassContainer ) {
	
    for (int i = 0; i < myBiopolymerClassContainer.getNumBiopolymers(); i++) {
	BiopolymerClass       myBiopolymerClass = myBiopolymerClassContainer.updBiopolymerClass(i);
        if (myBiopolymerClass.getBiopolymerType() == BiopolymerType::Protein){
	    //String chain = myBiopolymerClass.getChainID();
            addPhiPsiMobility(myBiopolymerClass.getChainID() , myBiopolymerClass.getFirstResidueID(), myBiopolymerClass.getLastResidueID(), bondMobilityString, myBiopolymerClassContainer);
	    
        }
    }
};

void MobilizerContainer::deleteMobilizerWithin(int id){
    mobilizerWithinVector.erase(mobilizerWithinVector.begin()+id);
}

void MobilizerContainer::updateMobilizerWithin(int id, String myChain, ResidueID myRes, double myRadius, String bondMobilityString, BiopolymerClassContainer & myBiopolymerClassContainer){
    MobilizerWithin & mw = mobilizerWithinVector[id];
    mw.setChain ( myChain);
    mw.setResidue  (myRes);
    mw.setRadius  (myRadius);
    mw.setBondMobilityString ( bondMobilityString);
    validateMobilizerWithin(mw, myBiopolymerClassContainer);
}

void MobilizerContainer::setMobilizerTypeForAllChains(const String myMobilizerString, BiopolymerClassContainer & myBiopolymerClassContainer){
    MMBLOG_FILE_FUNC_LINE(INFO, "myMobilizerString = >"<<myMobilizerString<<"< "<<endl);
    for (int i = 0 ; i < myBiopolymerClassContainer.getNumBiopolymers(); i++) {
        String myChainID = myBiopolymerClassContainer.updBiopolymerClass(i).getChainID();
        //String myMobilizerString = parameterStringClass.getString(1);
        MMBLOG_FILE_FUNC_LINE(INFO, "Adding mobilizer stretch to biopolymer index "<<i<<" , chain "<< myChainID<<endl);
        addMobilizerStretchToVector(
            myChainID,
            myMobilizerString,
            myBiopolymerClassContainer
            );
    } // of for
}
