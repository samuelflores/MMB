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
//#include <seqan/align.h>
#include <cmath> // contains sqrt
#include <ctgmath>
//#include "/home/sam/svn/molmodel/include/molmodel/internal/Superpose.h"
//#include "/home/sam/svn/molmodel/include/molmodel/internal/Superpose.h"


//#include "superimposer.h"

using namespace std;
using namespace SimTK;


// This function returns the Vec3 connecting the endpoints of a single AtomSpring


Vec3 getLeftAtomLocationFromAtomSpring(State & state,  BiopolymerClassContainer & biopolymerClassContainer, AtomSpring & atomSpring){
    MMBLOG_FILE_FUNC_LINE(INFO, "atomSpring.atom1Residue, atomSpring.atom1Name  = "<<atomSpring.atom1Residue.outString() <<", "<<   atomSpring.atom1Name<<endl);
    return biopolymerClassContainer.updBiopolymerClass(atomSpring. atom1Chain).calcAtomLocationInGroundFrame(state,    
        atomSpring.atom1Residue ,   atomSpring.atom1Name);
    
}
Vec3 getRightAtomLocationFromAtomSpring(State & state,  BiopolymerClassContainer & biopolymerClassContainer, AtomSpring & atomSpring){
    MMBLOG_FILE_FUNC_LINE(INFO, "atomSpring.atom2Residue ,   atomSpring.atom2Name  = "<<atomSpring.atom2Residue.outString() <<", "<<   atomSpring.atom2Name<<endl);
    return biopolymerClassContainer.updBiopolymerClass(atomSpring. atom2Chain).calcAtomLocationInGroundFrame(state,    
        atomSpring.atom2Residue ,   atomSpring.atom2Name);
}

Vec3 getDisplacementVec(State & state,  BiopolymerClassContainer & biopolymerClassContainer, AtomSpring atomSpring){
    Vec3 atom1Location = getLeftAtomLocationFromAtomSpring(state , biopolymerClassContainer , atomSpring);
    //biopolymerClassContainer.updBiopolymerClass(atomSpring. atom1Chain).calcAtomLocationInGroundFrame(state,    
        //atomSpring.atom1Residue ,   atomSpring.atom1Name);
    Vec3 atom2Location = getRightAtomLocationFromAtomSpring(state , biopolymerClassContainer , atomSpring);
    //biopolymerClassContainer.updBiopolymerClass(atomSpring. atom2Chain).calcAtomLocationInGroundFrame(state,    
        //atomSpring.atom2Residue ,   atomSpring.atom2Name);
    Vec3 displacementVec = atom2Location - atom1Location;
    return displacementVec;
}

 // This function returns the extension of a single AtomSpring

double getExtension(State & state,  BiopolymerClassContainer & biopolymerClassContainer, AtomSpring atomSpring) {
    Vec3 displacementVec = getDisplacementVec(state, biopolymerClassContainer, atomSpring);
    double myDotProduct = DotProduct(displacementVec,displacementVec); 
    double displacementVecMagnitude = sqrt(myDotProduct); // square root of dot product gives us vector modulus
    return displacementVecMagnitude;
}

// this function computes RMSD where "D" is the extension of each AtomSpring in atomSpringVector.

double AtomSpringContainer::calcRmsd(State & state, BiopolymerClassContainer & biopolymerClassContainer){
    double sumSquareExtension = 0.0;
    double myExtension = 0.0;
    int numSprings = 0;
    for (size_t i = 0; i <  atomSpringVector.size() ; i++) { // Count over all AtomSpring's
        Vec3 myDisplacementVec = (getDisplacementVec(state, biopolymerClassContainer, atomSpringVector[i]));
        sumSquareExtension += DotProduct(myDisplacementVec,myDisplacementVec)   ; // Dot product of myDisplacementVec is the square of the extension
        numSprings++;
    }
    double myRmsd = sqrt(sumSquareExtension/numSprings);
    return myRmsd;
}

// This function calls the Kabsch algorith, as implemented in superimposer.cpp . This will return the minimum (optimum) attainable RMSD under a rigid translation-rotation assumption. It assumes all atoms on the left side of the AtomSpringContainer are to be superimposed rigidly on all atoms on the right side. This will of course not hold if there is any flexibility between the atoms on one or the other side.
float AtomSpringContainer::calcKabschRmsd(State & state, BiopolymerClassContainer & biopolymerClassContainer){
    // Algorithm takes 2D vector of the form [ [x], [y], [z] ]
    std::vector< std::vector<float> > coord0; 
    std::vector< std::vector<float> > coord1;
    //for (int left = 0 ; left < atomSpringVector.size() ; left ++)
    Kabsch78::VectorSet myVectorSet; myVectorSet.clear();
    double totalExtension = 0.;
    for (size_t i = 0; i <  atomSpringVector.size() ; i++) { // Count over all AtomSpring's
        MMBLOG_FILE_FUNC_LINE(INFO, endl);
        // getLeftAtomLocationFromAtomSpring, getRightAtomLocationFromAtomSpring both return Vec3 but we need vector<float>
        Vec3 leftVec3  = getLeftAtomLocationFromAtomSpring (state , biopolymerClassContainer , atomSpringVector[i] )  ;
        Vec3 rightVec3(-9999.9,-9999.9,-9999.9);
        if (atomSpringVector[i].toGround){
            MMBLOG_FILE_FUNC_LINE(INFO, "Detected that you have some sort of spring or tether to ground. Please remove these if you want to call calcKabschRmsd. If you are not intereted in Kabsch you can ignore this message"<<endl);
            rightVec3=atomSpringVector[i].groundLocation;
            // This is not the sort of spring we want to do a Kabsch algorithm calc on .
            return -9999.9;
        } else {
            MMBLOG_FILE_FUNC_LINE(INFO, endl);
            rightVec3 = getRightAtomLocationFromAtomSpring(state , biopolymerClassContainer , atomSpringVector[i] )  ;}
        double myExtension = getExtension(state, biopolymerClassContainer, atomSpringVector[i]);
        totalExtension += (myExtension*myExtension);
        MMBLOG_FILE_FUNC_LINE(INFO, "Comparing left location : "<< leftVec3<< " vs. right location : "<<rightVec3<<" . This gives extension of : "<<myExtension<<endl);
        Vec3Pair myVec3Pair(leftVec3, rightVec3);
        myVectorSet.push_back(myVec3Pair);
         
    }
     
    MMBLOG_FILE_FUNC_LINE(INFO, "atomSpringVector.size() = "<<atomSpringVector.size()<<endl);
    TransformAndResidual myTransformAndResidual = Kabsch78::superpose(myVectorSet);
    if (atomSpringVector.size() > 0) {
        totalExtension = sqrt(totalExtension/ atomSpringVector.size());
        MMBLOG_FILE_FUNC_LINE(INFO, "Detected <atomSpringVector.size() > 0. Kabsch algorithm gives RMSD for atomSpring's of : "<<myTransformAndResidual.residual<<" pre-Kabsch extension was "<<totalExtension<<endl);}
    return myTransformAndResidual.residual;
        
}

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
   atomSpring.deadLengthIsFractionOfInitialLength = 0;
   return atomSpring;
};


void AtomSpringContainer::printAtomSpring(const AtomSpring atomSpring){
    MMBLOG_FILE_FUNC_LINE(INFO,
          "atom1Chain     = " << atomSpring.atom1Chain    
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
        <<endl);
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
    ValidateDouble(atomSpring.forceConstant);
    ValidateDouble(atomSpring.deadLength);
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
    ValidateDouble(atomSpring.forceConstant);
    ValidateDouble(atomSpring.deadLength);
};

void AtomSpringContainer::addAtomSpring(const AtomSpring & atomSpring, BiopolymerClassContainer & myBiopolymerClassContainer){
    validateAtomSpring(atomSpring, myBiopolymerClassContainer);
    atomSpringVector.push_back(atomSpring);
}

void AtomSpringContainer::deleteAtomSpring(int id){
    if(id < 0 || id >= atomSpringVector.size()){
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "you tried to delete a non existing AtomSpring." << endl);
    }
    atomSpringVector.erase(atomSpringVector.begin()+id);
}

void AtomSpringContainer::updateAtomSpring(const int id, const AtomSpring & newSpring, BiopolymerClassContainer & myBiopolymerClassContainer){
    if(id < 0 || id >= atomSpringVector.size()){
        MMBLOG_FILE_FUNC_LINE(INFO, "you tried to update a non existing AtomSpring." << endl);
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
    BiopolymerClass  bpc1 = (thread.getThreadingPartner(0).biopolymerClass);// myBiopolymerClassContainer.updBiopolymerClass(thread.updThreadingPartner(0).biopolymerClass.getChainID());
    BiopolymerClass  bpc2 = (thread.getThreadingPartner(1).biopolymerClass);//myBiopolymerClassContainer.updBiopolymerClass(thread.updThreadingPartner(1).biopolymerClass.getChainID());
    //BiopolymerClass & bpc1 = myBiopolymerClassContainer.updBiopolymerClass(thread.updThreadingPartner(0).biopolymerClass.getChainID());
    //BiopolymerClass & bpc2 = myBiopolymerClassContainer.updBiopolymerClass(thread.updThreadingPartner(1).biopolymerClass.getChainID());
    //BiopolymerClass & bpc2 = myBiopolymerClassContainer.updBiopolymerClass(thread.);

    if(bpc1.getBiopolymerType() != bpc2.getBiopolymerType())
    {
       MMBLOG_FILE_FUNC_LINE(CRITICAL, "In the threading command, both chains must be of the same type." << endl);
    }

    if( thread.getThreadingPartner(0).startResidue > thread.getThreadingPartner(0).endResidue)
    {
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "In the threading command, the end residue must be greater than or equal to the start residue for each chain." << endl);
    }
    if( thread.getThreadingPartner(1).startResidue  > thread.getThreadingPartner(1).endResidue )
    {
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "In the threading command, the end residue must be greater than or equal to the start residue for each chain." << endl);
    }
    if( bpc1.difference(thread.getThreadingPartner(0).endResidue,thread.getThreadingPartner(0).startResidue) != bpc2.difference(thread.getThreadingPartner(1).endResidue,thread.getThreadingPartner(1).startResidue))
    {
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "In the threading command, the two threaded segments must be of the same length." << endl);
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
    ThreadingStruct thread(myBiopolymerClassContainer.updBiopolymerClass(chain1) , resStart1, resEnd1,myBiopolymerClassContainer.updBiopolymerClass( chain2), resStart2, resEnd2, forceConstant, backboneOnly);
    addThreading(thread, myBiopolymerClassContainer);
}

void AtomSpringContainer::deleteThreading(int id){
    if(id < 0 || id >= threadingStructVector.size()){
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "you tried to delete a non existing Threading." << endl);
    }
    threadingStructVector.erase(threadingStructVector.begin()+id);
}

void AtomSpringContainer::updateThreading(int id, const ThreadingStruct & newThread, 
                                          BiopolymerClassContainer & myBiopolymerClassContainer){
    if(id < 0 || id >= threadingStructVector.size()){
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "you tried to update a non existing Threading." << endl);
    }
    validateThreading(newThread, myBiopolymerClassContainer);
    threadingStructVector[id] = std::move(newThread);
}

void AtomSpringContainer::updateThreading(int id, String chain1, ResidueID resStart1, ResidueID resEnd1, 
                                          String chain2, ResidueID resStart2, ResidueID resEnd2, 
                                          double forceConstant, bool backboneOnly, 
                                          BiopolymerClassContainer & myBiopolymerClassContainer){
    if(id < 0 || id >= threadingStructVector.size()){
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "you tried to update a non existing Threading." << endl);
    }
    ThreadingStruct newThread(myBiopolymerClassContainer.updBiopolymerClass(chain1), resStart1, resEnd1,myBiopolymerClassContainer.updBiopolymerClass( chain2), resStart2, resEnd2, forceConstant, backboneOnly);
    this->updateThreading(id, newThread, myBiopolymerClassContainer);
}

void AtomSpringContainer::createSpringsFromThreading(BiopolymerClassContainer & myBiopolymerClassContainer)
{
    vector<ThreadingStruct>::iterator it;
    for(it = threadingStructVector.begin(); it != threadingStructVector.end(); it++){
        ThreadingStruct & thread = *it;
        BiopolymerClass  bp1 = (thread.getThreadingPartner(0).biopolymerClass) ; //myBiopolymerClassContainer.updBiopolymerClass(thread.getThreadingPartner(0).biopolymerClass.getChainID()  );
        int threadedLength = bp1.difference (thread.getThreadingPartner(0).endResidue  ,thread.getThreadingPartner(0).startResidue  ) + 1;

        if(bp1.getBiopolymerType() == BiopolymerType::Protein && thread.backboneOnly)
        {
            for (int i = 0; i < threadedLength; i++) 
            {
                AtomSpring myAtomSpring1(thread.getThreadingPartner(0).biopolymerClass.getChainID(), thread.getThreadingPartner(0).biopolymerClass.sum(thread.getThreadingPartner(0).startResidue , i), String("N"),
                                         thread.getThreadingPartner(1).biopolymerClass.getChainID(), thread.getThreadingPartner(1).biopolymerClass.sum(thread.getThreadingPartner(1).startResidue , i), String("N"),
                                         thread.forceConstant,thread.deadLength,thread.deadLengthFraction
                                        );

                AtomSpring myAtomSpring2(thread.getThreadingPartner(0).biopolymerClass.getChainID(), thread.getThreadingPartner(0).biopolymerClass.sum(thread.getThreadingPartner(0).startResidue , i), String("CA"),
                                         thread.getThreadingPartner(1).biopolymerClass.getChainID(), thread.getThreadingPartner(1).biopolymerClass.sum(thread.getThreadingPartner(1).startResidue , i), String("CA"),
                                         thread.forceConstant,thread.deadLength,thread.deadLengthFraction
                                        );

                AtomSpring myAtomSpring3(thread.getThreadingPartner(0).biopolymerClass.getChainID(), thread.getThreadingPartner(0).biopolymerClass.sum(thread.getThreadingPartner(0).startResidue , i), String("C"),
                                         thread.getThreadingPartner(1).biopolymerClass.getChainID(), thread.getThreadingPartner(1).biopolymerClass.sum(thread.getThreadingPartner(1).startResidue , i), String("C"),
                                         thread.forceConstant,thread.deadLength,thread.deadLengthFraction
                                        );
                //myAtomSpring1.deadLength = thread.deadLength;
                //myAtomSpring2.deadLength = thread.deadLength;
                //myAtomSpring3.deadLength = thread.deadLength;
                this->add   (myAtomSpring1);
                this->add   (myAtomSpring2);
                this->add   (myAtomSpring3);

            }        
        }
        else
        {
            for (int i = 0; i < threadedLength; i++) 
            {
                ResidueInfo myResidueInfoA = thread.getThreadingPartner(0).biopolymerClass.updResidueInfo(thread.getThreadingPartner(0).biopolymerClass.sum(thread.getThreadingPartner(0).startResidue  , i));
                //ResidueInfo myResidueInfoA = myBiopolymerClassContainer.updBiopolymerClass(thread.chainID1).updResidueInfo(myBiopolymerClassContainer.updBiopolymerClass(thread.chainID1).sum( thread.residueStart1 , i));

                ResidueInfo myResidueInfoB = thread.getThreadingPartner(1).biopolymerClass. updResidueInfo(thread.getThreadingPartner(1).biopolymerClass. sum(thread.getThreadingPartner(1).startResidue , i));
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
                        if (thread.getThreadingPartner(1).biopolymerClass.hasAtom( thread.getThreadingPartner(1).biopolymerClass.sum (thread.getThreadingPartner(1).startResidue , i),atomNameA)) {

                            AtomSpring myAtomSpring1(thread.getThreadingPartner(0).biopolymerClass. getChainID(), 
                                                     thread.getThreadingPartner(0).biopolymerClass. sum(thread.getThreadingPartner(0).startResidue , i), 
                                                     atomNameA,
                                                     thread.getThreadingPartner(1).biopolymerClass. getChainID(), 
                                                     thread.getThreadingPartner(1).biopolymerClass. sum(thread.getThreadingPartner(1).startResidue , i), 
                                                     atomNameA,
                                                     thread.forceConstant,
                                                     thread.deadLength,
                                                     thread.deadLengthFraction
                                                    );
                            MMBLOG_FILE_FUNC_LINE(DEBUG," Created atomSpring for proteinThreading: "<<endl);
			    //atomNameA, thread.chainID1, residueA, thread.chainID2, residueB : >"<<atomNameA<<"< " << thread.getThreadingPartner(0).biopolymerClass. getChainID()  <<", "<<  <<", "<<thread.chainID2<<", "<<thread.residueStart2 + i  <<endl);
                            //MMBLOG_FILE_FUNC_LINE(DEBUG," Created atomSpring for proteinThreading: atomNameA, thread.chainID1, residueA, thread.chainID2, residueB : >"<<atomNameA<<"< " << thread.getThreadingPartner(0).biopolymerClass. getChainID()  <<", "<<  <<", "<<thread.chainID2<<", "<<thread.residueStart2 + i  <<endl);
                            MMBLOG_FILE_FUNC_LINE(DEBUG, " Adding AtomSpring from ThreadingStruct: "<<endl);
			    myAtomSpring1.printDebug();

                            //"Created a ThreadingStruct connecting chain "<<chain1<<" residue "<<thread.updThreadingPartner(0).startResidue.outString() <<" to "<<thread.updThreadingPartner(0).endResidue.outString()<<" . "<<endl);
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

// Creates a gapped alignment using only explicitly specified stretches of residues in the two aligned chains.
ThreadingStruct AtomSpringContainer::createGappedThreading(String chain1, ResidueID startResidue1,  ResidueID endResidue1, String chain2,  ResidueID startResidue2,  ResidueID endResidue2,  double forceConstant, bool backboneOnly, BiopolymerClassContainer & myBiopolymerClassContainer)
{
    ThreadingStruct thread;
    thread.updThreadingPartner(0).biopolymerClass = myBiopolymerClassContainer.updBiopolymerClass(chain1);
    thread.updThreadingPartner(1).biopolymerClass = myBiopolymerClassContainer.updBiopolymerClass(chain2);
    //thread.chainID1 = chain1;
    //thread.chainID2 = chain2;
    MMBLOG_FILE_FUNC_LINE(INFO, endl);
    ResidueStretch  myResidueStretch1 = ResidueStretch(chain1, startResidue1, endResidue1);
    ResidueStretch  myResidueStretch2 = ResidueStretch(chain2, startResidue2, endResidue2);
    if (!(myBiopolymerClassContainer.updBiopolymerClass(chain1).hasResidueStretch(myResidueStretch1))) {
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "You have specified an invalid residue stretch: Chain "<<chain1<<" from residue "<<startResidue1.outString()<<" to "<<endResidue1.outString() << endl);
    }
    if (!(myBiopolymerClassContainer.updBiopolymerClass(chain2).hasResidueStretch(myResidueStretch2))) {
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "You have specified an invalid residue stretch: Chain "<<chain2<<" from residue "<<startResidue2.outString()<<" to "<<endResidue2.outString() << endl);
    }
    thread.updThreadingPartner(0).startResidue = startResidue1;//  myBiopolymerClassContainer.updBiopolymerClass(thread.chainID1).getFirstResidueID();
    thread.updThreadingPartner(1).startResidue = startResidue2; //myBiopolymerClassContainer.updBiopolymerClass(thread.chainID2).getFirstResidueID();
    thread.updThreadingPartner(0).endResidue   = endResidue1; // myBiopolymerClassContainer.updBiopolymerClass(thread.chainID1).getLastResidueID();
    thread.updThreadingPartner(1).endResidue   = endResidue2; //myBiopolymerClassContainer.updBiopolymerClass(thread.chainID2).getLastResidueID();
    thread.forceConstant = forceConstant;
    thread.backboneOnly = backboneOnly;
    MMBLOG_FILE_FUNC_LINE(INFO, "Created a ThreadingStruct connecting chain "<<chain1<<" residue "<<thread.updThreadingPartner(0).startResidue.outString() <<" to "<<thread.updThreadingPartner(0).endResidue.outString()<<" . "<<endl);
    MMBLOG_FILE_FUNC_LINE(INFO, ":                                   to chain "<<chain2<<" residue "<<thread.updThreadingPartner(1).startResidue.outString() <<" to "<<thread.updThreadingPartner(1).endResidue.outString()<<" . "<<endl);
    return thread;
}

// Creates a gapped alignment using all residues in the two aligned chains.
ThreadingStruct AtomSpringContainer::createGappedThreading(String chain1, String chain2, double forceConstant, bool backboneOnly, BiopolymerClassContainer & myBiopolymerClassContainer){
    MMBLOG_FILE_FUNC_LINE(INFO, endl);
    ResidueID myStartResidue1 = myBiopolymerClassContainer.updBiopolymerClass(chain1).getFirstResidueID();
    ResidueID myStartResidue2 = myBiopolymerClassContainer.updBiopolymerClass(chain2).getFirstResidueID();
    ResidueID myEndResidue1   = myBiopolymerClassContainer.updBiopolymerClass(chain1).getLastResidueID();
    ResidueID myEndResidue2   = myBiopolymerClassContainer.updBiopolymerClass(chain2).getLastResidueID();
    return createGappedThreading((chain1), myStartResidue1, myEndResidue1 ,(chain2) , myStartResidue2,  myEndResidue2, forceConstant, backboneOnly, myBiopolymerClassContainer);
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
void AtomSpringContainer::addGappedThreading(String chain1, ResidueID startResidue1,  ResidueID endResidue1,   String chain2, ResidueID startResidue2,  ResidueID endResidue2,  double forceConstant, bool backboneOnly, 
                                       BiopolymerClassContainer & myBiopolymerClassContainer){
    ThreadingStruct thread = createGappedThreading(chain1, startResidue1,endResidue1, chain2,startResidue2,endResidue2, forceConstant, backboneOnly, myBiopolymerClassContainer);
    addGappedThreading(thread, myBiopolymerClassContainer);
}

void AtomSpringContainer::deleteGappedThreading(int id){
    if(id < 0 || id >= gappedThreadingStructVector.size()){
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "you tried to delete a non existing Threading." << endl);
    }
    gappedThreadingStructVector.erase(gappedThreadingStructVector.begin()+id);
}

void AtomSpringContainer::updateGappedThreading(int id, const ThreadingStruct & newThread, 
                                          BiopolymerClassContainer & myBiopolymerClassContainer){
    if(id < 0 || id >= gappedThreadingStructVector.size()){
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "you tried to update a non existing Threading." << endl);
    }
    // validateThreading(newThread, myBiopolymerClassContainer);
    gappedThreadingStructVector[id] = newThread;
}

void AtomSpringContainer::updateGappedThreading(int id, String chain1, String chain2, double forceConstant, bool backboneOnly, 
                                          BiopolymerClassContainer & myBiopolymerClassContainer){
    if(id < 0 || id >= gappedThreadingStructVector.size()){
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "you tried to update a non existing Threading." << endl);
    }
    ThreadingStruct newThread = createGappedThreading(chain1, chain2, forceConstant, backboneOnly, myBiopolymerClassContainer);
    this->updateGappedThreading(id, newThread, myBiopolymerClassContainer);
}

//seqan::AlignmentStats 
void AtomSpringContainer::createSpringsFromGappedThreading(BiopolymerClassContainer & myBiopolymerClassContainer)
{
    vector<ThreadingStruct>::iterator it;
    for(it = gappedThreadingStructVector.begin(); it != gappedThreadingStructVector.end(); it++)
    {
        ThreadingStruct & thread = *it;
        cout << "ThreadForceConstant " << thread.forceConstant << endl;
        String chainA = thread.updThreadingPartner(0).biopolymerClass. getChainID();//thread.chainID1;
        String chainB = thread.updThreadingPartner(1).biopolymerClass. getChainID();
        BiopolymerClass & bpA = myBiopolymerClassContainer.updBiopolymerClass(chainA);
        BiopolymerClass & bpB = myBiopolymerClassContainer.updBiopolymerClass(chainB);

        //TSequence seqA = bpA.getSubSequence(thread.updThreadingPartner(0).startResidue,thread.updThreadingPartner(0).endResidue).c_str();  // Need a new BiopolymerClass method which retrieves subsequences.!
        //TSequence seqB = bpB.getSubSequence(thread.updThreadingPartner(1).startResidue, thread.updThreadingPartner(1).endResidue ).c_str();
        // the above is now done by:
        thread. setLongSequences();
        TAlign align = thread.computeAlign(); //setLongSequences also calls computeAlign, but we are being paranoid. Plus this gives us a convenient return value.
        //TAlign align;
        /*
        seqan::resize(rows(align), 2);
        assignSource(row(align,0),seqA);
        assignSource(row(align,1),seqB); 
        // simple alignment: 
        seqan::Blosum62 scoringScheme(-1, -12);
        int score = globalAlignment(align,scoringScheme ); // ..signature:Score<TValue, Simple>(match, mismatch, gap [, gap_open])

        std::cout << "Score: " << score << ::endl;
        std::cout << align << ::endl;

        seqan::AlignmentStats stats;
        computeAlignmentStats(stats, align, scoringScheme);
        thread.setAlignmentStats(stats); 
        thread.printAlignmentStats(); // This is done in computeAlign 
        cout<<__FILE__<<":"<<__LINE__<<" : "<< seqan::row(align,0)<<endl;
        cout<<__FILE__<<":"<<__LINE__<<" : "<< seqan::row(align,1)<<endl;
        */
        int aIndex = 0; int bIndex = 0; // Indices which count over residues in chains A and B.

        int i  = 0; // counts over columns in alignment
        //while ((aIndex < bpA.getSubSequence(thread.residueStart1,thread.residueEnd1 ).length()) && 
        //       (bIndex < bpB.getSubSequence(thread.residueStart2, thread.residueEnd2 ).length()   )) 
        while ((aIndex < length(thread.updThreadingPartner(0).sequence)) && 
               (bIndex < length(thread.updThreadingPartner(1).sequence))) 
        {
            if ((String(seqan::row (align,0)[i]).compare("-")  != 0  )  &&
                (String(seqan::row (align,1)[i]).compare("-")  != 0  )) 
            { 
                MMBLOG_FILE_FUNC_LINE(INFO, "Applying threading forces to "<<thread.updThreadingPartner(0).biopolymerClass. getChainID()<<" : "<<bpA.sum(thread.updThreadingPartner(0).startResidue   , aIndex).outString()<<" to "<<thread.updThreadingPartner(1).biopolymerClass. getChainID()<<" : "<<bpB.sum(thread.updThreadingPartner(1).startResidue ,bIndex).outString()<<endl);
                ResidueInfo myResidueInfoA = bpA.updResidueInfo(bpA.sum(thread.updThreadingPartner(0).startResidue   , aIndex  )) ;
                ResidueInfo myResidueInfoB = bpB.updResidueInfo(bpB.sum(thread.updThreadingPartner(1).startResidue   , bIndex)) ;
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
                        if (bpB.hasAtom(bpB.sum(thread.updThreadingPartner(1).startResidue   , bIndex),atomNameA)) 
                        {

                            AtomSpring myAtomSpring1(chainA, bpA.sum(thread.updThreadingPartner(0).startResidue   , aIndex), atomNameA,
                                                     chainB, bpB.sum(thread.updThreadingPartner(1).startResidue   , bIndex), atomNameA,
                                                     thread.forceConstant,thread.deadLength,thread.deadLengthFraction
                                                    );

                            //myAtomSpring1.deadLengthIsFractionOfInitialLength = thread.deadLengthIsFractionOfInitialLength;
                            //myAtomSpring1.deadLengthFraction = thread.deadLengthFraction;
                            //myAtomSpring1.deadLength = thread.deadLength;

                            MMBLOG_FILE_FUNC_LINE(DEBUG, " Adding AtomSpring from ThreadingStruct: "<<endl);
                            myAtomSpring1.printDebug(); 

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
    //return stats;
}

