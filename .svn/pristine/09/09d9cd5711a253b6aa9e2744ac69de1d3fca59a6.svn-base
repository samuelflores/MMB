/* -------------------------------------------------------------------------- *
 *                           MMB (MacroMoleculeBuilder)                       *
 * -------------------------------------------------------------------------- *
 *                                                                            *
 * Copyright (c) 2011-12 by the Author.                                       *
 * Author: Samuel Flores                                                      *
 *                                                                            *
 * See RNABuilder.cpp for the copyright and usage agreement.                  *
 * -------------------------------------------------------------------------- */

#ifndef MobilizerContainer_H_
#define MobilizerContainer_H_

#include "BiopolymerClass.h"
#include "ResidueStretchContainer.h"
 
class MMB_EXPORT MobilizerContainer : public ResidueStretchContainer <MobilizerStretch> {
private:
    vector <MobilizerWithin> mobilizerWithinVector;

public:
    void clear();
    bool isEmpty();
    MobilizerContainer(){clear();}
    vector<SingleBondMobility> singleBondMobilityVector;	
    //InterfaceContainer mobilizerInterfaceContainer;

    void validateMobilizerStretch(int mobilizerStretchIndex, BiopolymerClassContainer & myBiopolymerClassContainer);
    void validateMobilizerStretch(MobilizerStretch & ms, BiopolymerClassContainer & myBiopolymerClassContainer);
    void addMobilizerStretchToVector(MobilizerStretch myMobilizerStretch, BiopolymerClassContainer & myBiopolymerClassContainer) ;
    void addMobilizerStretchToVector(String myChain, ResidueID myStartResidue, ResidueID myEndResidue, String bondMobilityString, BiopolymerClassContainer & myBiopolymerClassContainer) ;
    void addMobilizerStretchToVector(String myChain, String bondMobilityString, BiopolymerClassContainer & myBiopolymerClassContainer) ;

    void updateMobilizerStretch(int id, String myChain, ResidueID myStartResidue, ResidueID myEndResidue, String bondMobilityString, BiopolymerClassContainer & myBiopolymerClassContainer);

    void deleteMobilizerStretch(int id);

    void printMobilizerStretch(int mobilizerStretchIndex);
    void printMobilizerStretches();
    int numMobilizerStretches();
    String getChain(int mobilizerStretchIndex);
    ResidueID getStartResidue(int mobilizerStretchIndex);
    ResidueID getEndResidue(int mobilizerStretchIndex);
    void setBiopolymerBondMobility (BiopolymerClassContainer & myBiopolymerClassContainer);

    vector<MobilizerWithin> & getMobilizerWithinVector() { return mobilizerWithinVector; }
    void createMobilizersWithin ( BiopolymerClassContainer & myBiopolymerClassContainer);//, State & state );
    void pushMobilizerWithin ( MobilizerWithin mobilizerWithin, BiopolymerClassContainer & myBiopolymerClassContainer);
    void validateMobilizerWithin(MobilizerWithin mobilizerWithin ,  BiopolymerClassContainer & myBiopolymerClassContainer);
    int numMobilizerWithin();

    void deleteMobilizerWithin(int id);
    void updateMobilizerWithin(int id, String myChain, ResidueID myRes, double myRadius, String bondMobilityString, BiopolymerClassContainer & myBiopolymerClassContainer);

    void addPhiPsiMobility(String chain, ResidueID startResidue, ResidueID endResidue, String bondMobilityString , BiopolymerClassContainer& biopolymerClassContainer );
    void addPhiPsiMobility(String bondMobilityString , BiopolymerClassContainer& biopolymerClassContainer );
    void addMobilizerStretchToVector(MobilizerStretch mobilizerStretch, String  bondMobilityString , BiopolymerClassContainer & myBiopolymerClassContainer);
    void addMobilizerStretchesToVector(vector <MobilizerStretch> residueStretchVector, String  bondMobilityString , BiopolymerClassContainer & myBiopolymerClassContainer);
    void addMobilizerStretchesToVector(BiopolymerClassContainer&);
    //void addMobilizerStretchesToVector(BiopolymerClassContainer&, CompoundSystem&, State&);
    void addMobilizerStretchesToVector(vector<MobilizerStretch>, BiopolymerClassContainer&);
    void addSingleBondMobilityToAllChains(String bondMobilityString , BiopolymerClassContainer& biopolymerClassContainer );
    void addMobilizerDomainsInterfacesToVector(const vector<MobilizerDomainsInterface> & mDIVector, BiopolymerClassContainer & myBiopolymerClassContainer);

};

#endif
