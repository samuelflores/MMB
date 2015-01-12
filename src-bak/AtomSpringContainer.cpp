/* -------------------------------------------------------------------------- *
 *                           MMB (MacroMoleculeBuilder)                       *
 * -------------------------------------------------------------------------- *
 *                                                                            *
 * Copyright (c) 2011-12 by the Author.                                       *
 * Author: Samuel Flores                                                      *
 *                                                                            *
 * See RNABuilder.cpp for the copyright and usage agreement.                  *
 * -------------------------------------------------------------------------- */


#include "AtomSpringContainer.h"
#include "BiopolymerClass.h"
#include <seqan/align.h>

using namespace std;
using namespace SimTK;

AtomSpring & AtomSpringContainer::initializeAtomSpring(AtomSpring & atomSpring) {
   atomSpring.atom1Name = ""           ;   
   atomSpring.atom2Name = ""           ;   
   atomSpring.atom1Residue =  ResidueID(0, ' ')        ;   
   atomSpring.atom2Residue =  ResidueID(0, ' ')        ;   
   atomSpring.atom1Chain = "";
   atomSpring.atom2Chain = "";
   atomSpring.toGround   = false;
   atomSpring.tether     = false;
   atomSpring.groundLocation = Vec3(0);
   atomSpring.forceConstant  = 0.0 ;
   atomSpring.deadLength     = 0.0 ;   
   return atomSpring;
};


void AtomSpringContainer::printAtomSpring(const AtomSpring atomSpring){
    cout<<__FILE__<<":"<<__LINE__
        <<" atom1Chain     = " << atomSpring.atom1Chain    
        <<" atom1Residue   = " << atomSpring.atom1Residue.outString()   
	<<" atom1Name      = " << atomSpring.atom1Name     
        <<" atom2Chain     = " << atomSpring.atom2Chain    
        <<" atom2Residue   = " << atomSpring.atom2Residue.outString()   
        <<" atom2Name      = " << atomSpring.atom2Name      
        <<" toGround       = " << atomSpring.toGround      
        <<" tether         = " << atomSpring.tether        
        <<" groundLocation = " << atomSpring.groundLocation <<" (nm,nm,nm) "
        <<" forceConstant  = " << atomSpring.forceConstant <<" (kJ/mol/nm/nm) "
        <<" deadLength     = " << atomSpring.deadLength <<" (nm) "
        <<endl;     
};

void AtomSpringContainer::printAtomSpring(int atomSpringIndex){
    AtomSpring myAtomSpring = getAtomSpring (atomSpringIndex);
    printAtomSpring(myAtomSpring);
};

void AtomSpringContainer::printAtomSprings(){
    for (int i = 0 ; i < numAtomSprings(); i++) 
        printAtomSpring(i); 
};


void AtomSpringContainer::validateAtomSpring(const AtomSpring & atomSpring){//,  BiopolymerClassContainer & myBiopolymerContainer ){
    // a significant amount of validation is already being done in ParameterReader.cpp, when the atomSpring or related command is read.
    // the following two calls will ensure that the desired atoms exist:
    //myBiopolymerContainer.updBiopolymerClass(atomSpring.atom1Chain).atomPathString(atomSpring.atom1Residue,atomSpring.atom1Name);
    //if (! atomSpring.toGround) // if it's .toGround, then there is no second atom
    //  myBiopolymerContainer.updBiopolymerClass(atomSpring.atom2Chain).atomPathString(atomSpring.atom2Residue,atomSpring.atom2Name);
    // make sure groundLocation is not NaN or Inf.  This applies even if .toGround is false, because it should anyhow be well formed:
    //cout<<__FILE__<<":"<<__LINE__<<" Validating atom spring : "<<endl;
    //printAtomSpring(  atomSpring);
    ValidateVec3(atomSpring.groundLocation);
    ValidateReal(atomSpring.forceConstant);
    ValidateReal(atomSpring.deadLength);
};

void AtomSpringContainer::validateAtomSpring(const AtomSpring & atomSpring,  BiopolymerClassContainer & myBiopolymerContainer ){
    // a significant amount of validation is already being done in ParameterReader.cpp, when the atomSpring or related command is read.
    // the following two calls will ensure that the desired atoms exist:
    myBiopolymerContainer.updBiopolymerClass(atomSpring.atom1Chain).atomPathString(atomSpring.atom1Residue,atomSpring.atom1Name);
    if (! atomSpring.toGround) // if it's .toGround, then there is no second atom
    	myBiopolymerContainer.updBiopolymerClass(atomSpring.atom2Chain).atomPathString(atomSpring.atom2Residue,atomSpring.atom2Name);
    // make sure groundLocation is not NaN or Inf.  This applies even if .toGround is false, because it should anyhow be well formed:
    //cout<<__FILE__<<":"<<__LINE__<<" Validating atom spring : "<<endl;
    //printAtomSpring(  atomSpring);
    ValidateVec3(atomSpring.groundLocation);
    ValidateReal(atomSpring.forceConstant);
    ValidateReal(atomSpring.deadLength);
};

void AtomSpringContainer::addAtomSpring(const AtomSpring & atomSpring, BiopolymerClassContainer & myBiopolymerClassContainer){
    validateAtomSpring(atomSpring, myBiopolymerClassContainer);
    atomSpringVector.push_back(atomSpring);
}

void AtomSpringContainer::deleteAtomSpring(int id){
    if(id < 0 || id >= atomSpringVector.size()){
        ErrorManager::instance <<__FILE__<<":"<<__LINE__<<": you tried to delete a non existing AtomSpring." << endl;
        ErrorManager::instance.treatError();
    }
    atomSpringVector.erase(atomSpringVector.begin()+id);
}

void AtomSpringContainer::updateAtomSpring(const int id, const AtomSpring & newSpring, BiopolymerClassContainer & myBiopolymerClassContainer){
    if(id < 0 || id >= atomSpringVector.size()){
        ErrorManager::instance <<__FILE__<<":"<<__LINE__<<": you tried to update a non existing AtomSpring." << endl;
        ErrorManager::instance.treatError();
    }
    validateAtomSpring(newSpring, myBiopolymerClassContainer);
    atomSpringVector[id] = newSpring;
}

void AtomSpringContainer::clear()
{
    atomSpringVector.clear();
    clearThreading();
    clearGappedThreading();
}

void AtomSpringContainer::clearThreading(){
    threadingStructVector.clear();
}

void AtomSpringContainer::validateThreading(const ThreadingStruct & thread, BiopolymerClassContainer & myBiopolymerClassContainer){
    BiopolymerClass & bpc1 = myBiopolymerClassContainer.updBiopolymerClass(thread.chainID1);
    BiopolymerClass & bpc2 = myBiopolymerClassContainer.updBiopolymerClass(thread.chainID2);

    if(bpc1.getBiopolymerType() != bpc2.getBiopolymerType())
    {
       ErrorManager::instance << __FILE__ << " " << __LINE__ << ": In the threading command, botch chains must be of the same type." << endl;
       ErrorManager::instance.treatError(); 
    }

    if( thread.residueStart1 > thread.residueEnd1)
    {
        ErrorManager::instance << __FILE__ << " " << __LINE__ << ": In the threading command, the end residue must be greater than or equal to the start residue for each chain." << endl;
        ErrorManager::instance.treatError();
    }
    if( thread.residueStart2 > thread.residueEnd2)
    {
        ErrorManager::instance << __FILE__ << " " << __LINE__ << ": In the threading command, the end residue must be greater than or equal to the start residue for each chain." << endl;
        ErrorManager::instance.treatError();
    }
    if( bpc1.difference(thread.residueEnd1,thread.residueStart1) != bpc2.difference(thread.residueEnd2,thread.residueStart2))
    {
        ErrorManager::instance << __FILE__ << " " << __LINE__ << ": In the threading command, the two threaded segments must be of the same length." << endl;
        ErrorManager::instance.treatError();
    }
}

void AtomSpringContainer::addThreading(const ThreadingStruct & threadingStruct, 
                                       BiopolymerClassContainer & myBiopolymerClassContainer){
    validateThreading(threadingStruct, myBiopolymerClassContainer);
    threadingStructVector.push_back(threadingStruct);
}

void AtomSpringContainer::addThreading(String chain1, ResidueID resStart1, ResidueID resEnd1, 
                                       String chain2, ResidueID resStart2, ResidueID resEnd2, 
                                       double forceConstant, bool backboneOnly, 
                                       BiopolymerClassContainer & myBiopolymerClassContainer){
    ThreadingStruct thread(chain1, resStart1, resEnd1, chain2, resStart2, resEnd2, forceConstant, backboneOnly);
    addThreading(thread, myBiopolymerClassContainer);
}

void AtomSpringContainer::deleteThreading(int id){
    if(id < 0 || id >= threadingStructVector.size()){
        ErrorManager::instance <<__FILE__<<":"<<__LINE__<<": you tried to delete a non existing Threading." << endl;
        ErrorManager::instance.treatError();
    }
    threadingStructVector.erase(threadingStructVector.begin()+id);
}

void AtomSpringContainer::updateThreading(int id, const ThreadingStruct & newThread, 
                                          BiopolymerClassContainer & myBiopolymerClassContainer){
    if(id < 0 || id >= threadingStructVector.size()){
        ErrorManager::instance <<__FILE__<<":"<<__LINE__<<": you tried to update a non existing Threading." << endl;
        ErrorManager::instance.treatError();
    }
    validateThreading(newThread, myBiopolymerClassContainer);
    threadingStructVector[id] = newThread;
}

void AtomSpringContainer::updateThreading(int id, String chain1, ResidueID resStart1, ResidueID resEnd1, 
                                          String chain2, ResidueID resStart2, ResidueID resEnd2, 
                                          double forceConstant, bool backboneOnly, 
                                          BiopolymerClassContainer & myBiopolymerClassContainer){
    if(id < 0 || id >= threadingStructVector.size()){
        ErrorManager::instance <<__FILE__<<":"<<__LINE__<<": you tried to update a non existing Threading." << endl;
        ErrorManager::instance.treatError();
    }
    ThreadingStruct newThread(chain1, resStart1, resEnd1, chain2, resStart2, resEnd2, forceConstant, backboneOnly);
    this->updateThreading(id, newThread, myBiopolymerClassContainer);
}

void AtomSpringContainer::createSpringsFromThreading(BiopolymerClassContainer & myBiopolymerClassContainer)
{
    vector<ThreadingStruct>::iterator it;
    for(it = threadingStructVector.begin(); it != threadingStructVector.end(); it++){
        ThreadingStruct & thread = *it;
        BiopolymerClass & bp1 = myBiopolymerClassContainer.updBiopolymerClass(thread.chainID1);
        int threadedLength = bp1.difference (thread.residueEnd1 , thread.residueStart1) + 1;

        if(bp1.getBiopolymerType() == BiopolymerType::Protein && thread.backboneOnly)
        {
            for (int i = 0; i < threadedLength; i++) 
            {
                AtomSpring myAtomSpring1(thread.chainID1, myBiopolymerClassContainer.updBiopolymerClass(thread.chainID1).sum(thread.residueStart1 , i), String("N"),
                                         thread.chainID2, myBiopolymerClassContainer.updBiopolymerClass(thread.chainID2).sum(thread.residueStart2 , i), String("N"),
                                         thread.forceConstant
                                        );

                AtomSpring myAtomSpring2(thread.chainID1, myBiopolymerClassContainer.updBiopolymerClass(thread.chainID1).sum(thread.residueStart1 , i), String("CA"),
                                         thread.chainID2, myBiopolymerClassContainer.updBiopolymerClass(thread.chainID2).sum(thread.residueStart2 , i), String("CA"),
                                         thread.forceConstant
                                        );

                AtomSpring myAtomSpring3(thread.chainID1, myBiopolymerClassContainer.updBiopolymerClass(thread.chainID1).sum(thread.residueStart1 , i), String("C"),
                                         thread.chainID2, myBiopolymerClassContainer.updBiopolymerClass(thread.chainID2).sum(thread.residueStart2 , i), String("C"),
                                         thread.forceConstant
                                        );
                this->add   (myAtomSpring1);
                this->add   (myAtomSpring2);
                this->add   (myAtomSpring3);

            }        
        }
        else
        {
            for (int i = 0; i < threadedLength; i++) 
            {
                ResidueInfo myResidueInfoA = myBiopolymerClassContainer.updBiopolymerClass(thread.chainID1).updResidueInfo(myBiopolymerClassContainer.updBiopolymerClass(thread.chainID1).sum( thread.residueStart1 , i));

                ResidueInfo myResidueInfoB = myBiopolymerClassContainer.updBiopolymerClass(thread.chainID2).updResidueInfo(myBiopolymerClassContainer.updBiopolymerClass(thread.chainID2).sum(thread.residueStart2 , i));
                for (int j = 0; j < (int)myResidueInfoA.getNumAtoms(); j++) {
                    String atomNameA =  myResidueInfoA.getAtomName(ResidueInfo::AtomIndex (j));
                    if (
                        ((
                          (atomNameA.substr(0,1).compare("0") == 0) || 
                          (atomNameA.substr(0,1).compare("1") == 0) || 
                          (atomNameA.substr(0,1).compare("2") == 0) || 
                          (atomNameA.substr(0,1).compare("3") == 0) || 
                          (atomNameA.substr(0,1).compare("4") == 0) || 
                          (atomNameA.substr(0,1).compare("5") == 0) || 
                          (atomNameA.substr(0,1).compare("6") == 0) || 
                          (atomNameA.substr(0,1).compare("7") == 0) || 
                          (atomNameA.substr(0,1).compare("8") == 0) || 
                          (atomNameA.substr(0,1).compare("9") == 0) 
                         )
                         &&
                         (atomNameA.substr(1,1 ).compare("H") == 0)) ||
                        (atomNameA.substr(0,1 ).compare("H") == 0) 
                       ) 
                    { // do nothing; leaving out hydrogens
                    } else {
                        if (myBiopolymerClassContainer.updBiopolymerClass(thread.chainID2).hasAtom( myBiopolymerClassContainer.updBiopolymerClass(thread.chainID2).sum (thread.residueStart2 , i),atomNameA)) {

                            AtomSpring myAtomSpring1(thread.chainID1, 
                                                     myBiopolymerClassContainer.updBiopolymerClass(thread.chainID1).sum(thread.residueStart1 , i), 
                                                     atomNameA,
                                                     thread.chainID2, 
                                                     myBiopolymerClassContainer.updBiopolymerClass(thread.chainID2).sum(thread.residueStart2 , i), 
                                                     atomNameA,
                                                     thread.forceConstant
                                                    );
                            //cout<<__FILE__<<":"<<__LINE__<<" Created atomSpring for proteinThreading: atomNameA, thread.chainID1, residueA, thread.chainID2, residueB : >"<<atomNameA<<"< " <<thread.chainID1<<", "<<thread.residueStart1 + i<<", "<<thread.chainID2<<", "<<thread.residueStart2 + i  <<endl;
                            this->add(myAtomSpring1);

                        }
                    } //of if not H
                } // of for numatoms

            } // of for residues
        }
    }
}

void AtomSpringContainer::clearGappedThreading(){
    gappedThreadingStructVector.clear();
}

ThreadingStruct AtomSpringContainer::createGappedThreading(String chain1, String chain2, double forceConstant, bool backboneOnly, BiopolymerClassContainer & myBiopolymerClassContainer)
{
    ThreadingStruct thread;
    thread.chainID1 = chain1;
    thread.chainID2 = chain2;
    thread.residueStart1 = myBiopolymerClassContainer.updBiopolymerClass(thread.chainID1).getFirstResidueID();
    thread.residueStart2 = myBiopolymerClassContainer.updBiopolymerClass(thread.chainID2).getFirstResidueID();
    thread.residueEnd1 = myBiopolymerClassContainer.updBiopolymerClass(thread.chainID1).getLastResidueID();
    thread.residueEnd2 = myBiopolymerClassContainer.updBiopolymerClass(thread.chainID2).getLastResidueID();
    thread.forceConstant = forceConstant;
    thread.backboneOnly = backboneOnly;

    return thread;
}

void AtomSpringContainer::addGappedThreading(const ThreadingStruct & threadingStruct, 
                                       BiopolymerClassContainer & myBiopolymerClassContainer){
    // validateThreading(threadingStruct, myBiopolymerClassContainer);
    gappedThreadingStructVector.push_back(threadingStruct);
}

void AtomSpringContainer::addGappedThreading(String chain1, String chain2, double forceConstant, bool backboneOnly, 
                                       BiopolymerClassContainer & myBiopolymerClassContainer){
    ThreadingStruct thread = createGappedThreading(chain1, chain2, forceConstant, backboneOnly, myBiopolymerClassContainer);
    addGappedThreading(thread, myBiopolymerClassContainer);
}

void AtomSpringContainer::deleteGappedThreading(int id){
    if(id < 0 || id >= gappedThreadingStructVector.size()){
        ErrorManager::instance <<__FILE__<<":"<<__LINE__<<": you tried to delete a non existing Threading." << endl;
        ErrorManager::instance.treatError();
    }
    gappedThreadingStructVector.erase(gappedThreadingStructVector.begin()+id);
}

void AtomSpringContainer::updateGappedThreading(int id, const ThreadingStruct & newThread, 
                                          BiopolymerClassContainer & myBiopolymerClassContainer){
    if(id < 0 || id >= gappedThreadingStructVector.size()){
        ErrorManager::instance <<__FILE__<<":"<<__LINE__<<": you tried to update a non existing Threading." << endl;
        ErrorManager::instance.treatError();
    }
    // validateThreading(newThread, myBiopolymerClassContainer);
    gappedThreadingStructVector[id] = newThread;
}

void AtomSpringContainer::updateGappedThreading(int id, String chain1, String chain2, double forceConstant, bool backboneOnly, 
                                          BiopolymerClassContainer & myBiopolymerClassContainer){
    if(id < 0 || id >= gappedThreadingStructVector.size()){
        ErrorManager::instance <<__FILE__<<":"<<__LINE__<<": you tried to update a non existing Threading." << endl;
        ErrorManager::instance.treatError();
    }
    ThreadingStruct newThread = createGappedThreading(chain1, chain2, forceConstant, backboneOnly, myBiopolymerClassContainer);
    this->updateGappedThreading(id, newThread, myBiopolymerClassContainer);
}

void AtomSpringContainer::createSpringsFromGappedThreading(BiopolymerClassContainer & myBiopolymerClassContainer)
{
    vector<ThreadingStruct>::iterator it;
    for(it = gappedThreadingStructVector.begin(); it != gappedThreadingStructVector.end(); it++)
    {
        ThreadingStruct & thread = *it;
        cout << "ThreadForceConstant " << thread.forceConstant << endl;
        String chainA = thread.chainID1;
        String chainB = thread.chainID2;
        BiopolymerClass & bpA = myBiopolymerClassContainer.updBiopolymerClass(chainA);
        BiopolymerClass & bpB = myBiopolymerClassContainer.updBiopolymerClass(chainB);

        typedef seqan::String<char> TSequence;                 // sequence type
        typedef seqan::Align<TSequence,seqan::ArrayGaps> TAlign;      // align type           
        TSequence seqA = bpA.getSequence().c_str(); 
        TSequence seqB = bpB.getSequence().c_str();
        TAlign align;

        seqan::resize(rows(align), 2);
        assignSource(row(align,0),seqA);
        assignSource(row(align,1),seqB); 
        int score = globalAlignment(align, seqan::Score<int,seqan::Simple>(0,-1,-1));
        std::cout << "Score: " << score << ::std::endl;
        std::cout << align << ::std::endl;
        cout<<__FILE__<<":"<<__LINE__<<" : "<< seqan::row(align,0)<<endl;
        cout<<__FILE__<<":"<<__LINE__<<" : "<< seqan::row(align,1)<<endl;
        //String row0 = seqan::row(align,0);
        int sequenceALength = bpA.getSequence().length();
        int alignmentLength = sizeof(seqan::row(align,0) ); cout<<alignmentLength<<endl; 
        alignmentLength = sizeof(seqan::row(align,0)[300] ) ; cout<<alignmentLength<<endl; 
        alignmentLength = sizeof(seqan::row(align,0) ) / sizeof(seqan::row(align,0)[0] ) ; cout<<alignmentLength<<endl; 
        int aIndex = 0; int bIndex = 0; // Indices which count over residues in chains A and B.
        ResidueID startResidueA = bpA.getFirstResidueID();
        ResidueID startResidueB = bpB.getFirstResidueID();
        int i  = 0; // counts over columns in alignment
        while ((aIndex < bpA.getSequence().length()) && (bIndex < bpB.getSequence().length()   )) 
        {
            //cout<<__FILE__<<":"<<__LINE__<<" : "<< seqan::row (align,0)[i];
            //cout<<" : "<< seqan::row (align,1)[i]<<endl;
            if ((String(seqan::row (align,0)[i]).compare("-")  != 0  )  &&
                (String(seqan::row (align,1)[i]).compare("-")  != 0  )) 
            { 
                //cout<<__FILE__<<":"<<__LINE__<<" : These residues "<<bpA.sum(startResidueA  , aIndex  ).outString() <<" and "<<bpB.sum(startResidueB, bIndex).outString()<<" are aligned. Will thread them."<<endl;
                //cout<<__FILE__<<":"<<__LINE__<<" : Indices are "<<aIndex<<" and "<<bIndex<<endl;
                ResidueInfo myResidueInfoA = bpA.updResidueInfo(bpA.sum( startResidueA  , aIndex  )) ;
                ResidueInfo myResidueInfoB = bpB.updResidueInfo(bpB.sum(startResidueB, bIndex)) ;
                for (int j = 0; j < (int)myResidueInfoA.getNumAtoms(); j++) 
                {
                    String atomNameA =  myResidueInfoA.getAtomName(ResidueInfo::AtomIndex (j));
                    if (
                            ((
                              (atomNameA.substr(0,1).compare("0") == 0) || 
                              (atomNameA.substr(0,1).compare("1") == 0) || 
                              (atomNameA.substr(0,1).compare("2") == 0) || 
                              (atomNameA.substr(0,1).compare("3") == 0) || 
                              (atomNameA.substr(0,1).compare("4") == 0) || 
                              (atomNameA.substr(0,1).compare("5") == 0) || 
                              (atomNameA.substr(0,1).compare("6") == 0) || 
                              (atomNameA.substr(0,1).compare("7") == 0) || 
                              (atomNameA.substr(0,1).compare("8") == 0) || 
                              (atomNameA.substr(0,1).compare("9") == 0) 
                             )
                             &&
                             (atomNameA.substr(1,1 ).compare("H") == 0)) ||
                            (atomNameA.substr(0,1 ).compare("H") == 0) 
                       ) 
                    { // do nothing; leaving out hydrogens
                    } else 
                    {
                        if (bpB.hasAtom(bpB.sum(startResidueB , bIndex),atomNameA)) 
                        {

                            AtomSpring myAtomSpring1(chainA, bpA.sum(startResidueA , aIndex), atomNameA,
                                                     chainB, bpB.sum(startResidueB , bIndex), atomNameA,
                                                     thread.forceConstant
                                                    );

                            // AtomSpring myAtomSpring1;

                            // myAtomSpring1.atom1Chain   = chainA     ;
                            // myAtomSpring1.atom1Residue = bpA.sum(startResidueA , aIndex); // this is in PDB residue numbering
                            // myAtomSpring1.atom1Name    = atomNameA;  
                            // myAtomSpring1.atom2Chain   = chainB     ;
                            // myAtomSpring1.atom2Residue = bpB.sum(startResidueB , bIndex);
                            // myAtomSpring1.atom2Name    = atomNameA;   
                            // myAtomSpring1.groundLocation[0]    = 0.0;
                            // myAtomSpring1.groundLocation[1]    = 0.0;
                            // myAtomSpring1.groundLocation[2]    = 0.0;
                            // myAtomSpring1.toGround     = false;
                            // myAtomSpring1.tether       = false;
                            // myAtomSpring1.deadLength   = 0.0;  
                            // myAtomSpring1.forceConstant= myForceConstant;
                            //cout<<__FILE__<<":"<<__LINE__<<" Created atomSpring for proteinThreading: atomNameA, chainA, residueA, chainB, residueB : >"<<atomNameA<<"< " <<chainA<<", "<<startResidueA + i<<", "<<chainB<<", "<<startResidueB + i  <<endl;
                            this->add(myAtomSpring1);

                        }
                    } //of if not H
                } // of for numatoms
            } //End if seqan...
            if (String(seqan::row (align,0)[i]).compare("-")  != 0  ) {
                aIndex ++;
            }
            if (String(seqan::row (align,1)[i]).compare("-")  != 0  ) {
                bIndex ++;
            }
            i++;
        } // End While
    } // End for
}

