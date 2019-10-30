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
    #ifdef USE_OPENMM
    void createMobilizersWithin ( BiopolymerClassContainer & myBiopolymerClassContainer, State & state );
    #endif
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
    #ifdef USE_OPENMM
    void addMobilizerDomainsInterfacesToVector(const vector<MobilizerDomainsInterface> & mDIVector, BiopolymerClassContainer & myBiopolymerClassContainer);
    #endif
    void setMobilizerTypeForAllChains(const String myMobilizerString, BiopolymerClassContainer & myBiopolymerClassContainer);/*{
        cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" myMobilizerString = >"<<myMobilizerString<<"< "<<endl;
        for (int i = 0 ; i < myBiopolymerClassContainer.getNumBiopolymers(); i++) {
            String myChainID = myBiopolymerClassContainer.updBiopolymerClass(i).getChainID();
            //String myMobilizerString = parameterStringClass.getString(1);
                          cout<<__FILE__<<":"<<__LINE__<<" Adding mobilizer stretch to biopolymer index "<<i<<" , chain "<< myChainID<<endl;
            addMobilizerStretchToVector(
                myChainID,
                myMobilizerString,
                myBiopolymerClassContainer
                );
        } // of for
    }*/

};

#endif
