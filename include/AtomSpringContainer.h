/* -------------------------------------------------------------------------- *
 *                           MMB (MacroMoleculeBuilder)                       *
 * -------------------------------------------------------------------------- *
 *                                                                            *
 * Copyright (c) 2011-12 by the Author.                                       *
 * Author: Samuel Flores                                                      *
 *                                                                            *
 * See RNABuilder.cpp for the copyright and usage agreement.                  *
 * -------------------------------------------------------------------------- */

#ifndef AtomSpringContainer_H_
#define AtomSpringContainer_H_

#include "Utils.h"
#include "Threading.h"
#include <seqan/align.h>
//#include <Superpose.h>



class BiopolymerClassContainer;
typedef seqan::String<char> TSequence;                 // sequence type
typedef seqan::Align<TSequence,seqan::ArrayGaps> TAlign;

class MMB_EXPORT AtomSpringContainer {
private:
        std::vector <AtomSpring> atomSpringVector;
        std::vector <ThreadingStruct> threadingStructVector;
        std::vector <ThreadingStruct> gappedThreadingStructVector;
public:
        void clear();
        AtomSpringContainer() {clear(); };
        double calcRmsd      (State & state, BiopolymerClassContainer & biopolymerClassContainer);
        float  calcKabschRmsd(State & state, BiopolymerClassContainer & biopolymerClassContainer);
        //void printAtomSpring  (const AtomSpring   atomSpring);
        void printAtomSpring  (int atomSpringIndex);   
        void printAtomSprings ();   
        void validateAtomSpring(const AtomSpring & atomSpring );// ,   BiopolymerClassContainer & myBiopolymerContainer);
        void validateAtomSpring(const AtomSpring & atomSpring,   BiopolymerClassContainer & myBiopolymerContainer);
        AtomSpring & initializeAtomSpring(AtomSpring & atomSpring);
        void add(AtomSpring & atomSpring ){validateAtomSpring(atomSpring); atomSpringVector.push_back(atomSpring); };
        void addAtomSpring(const AtomSpring & atomSpring, BiopolymerClassContainer & myBiopolymerContainer);
        AtomSpring & updAtomSpring(int atomSpringIndex) {return atomSpringVector[atomSpringIndex]; };
        const AtomSpring getAtomSpring(int atomSpringIndex) {return atomSpringVector[atomSpringIndex]; };
        int numAtomSprings() {return atomSpringVector.size(); };
        const std::vector <AtomSpring> getAtomSpringVector( ){return atomSpringVector; };
        void deleteAtomSpring(const int atomSpringIndex);
        void updateAtomSpring(const int atomSpringIndex, const AtomSpring & newSpring, BiopolymerClassContainer & myBiopolymerClassContainer);

        void clearThreading();
        std::vector<ThreadingStruct> & getThreadingVector() { return threadingStructVector; }
        void validateThreading(const ThreadingStruct & threadingStruct, BiopolymerClassContainer & myBiopolymerClassContainer);
        void addThreading(const ThreadingStruct & threadingStruct, BiopolymerClassContainer & myBiopolymerClassContainer);
        void addThreading(String chain1, ResidueID resStart1, ResidueID resEnd1, String chain2, ResidueID resStart2, ResidueID resEnd2, double forceConstant, bool backboneOnly, BiopolymerClassContainer & myBiopolymerClassContainer);
        void deleteThreading(int id);
        void updateThreading(int id, const ThreadingStruct & threadingStruct, BiopolymerClassContainer & myBiopolymerClassContainer);
        void updateThreading(int id, String chain1, ResidueID resStart1, ResidueID resEnd1, String chain2, ResidueID resStart2, ResidueID resEnd2, double forceConstant, bool backboneOnly, BiopolymerClassContainer & myBiopolymerClassContainer);
        void createSpringsFromThreading(BiopolymerClassContainer & myBiopolymerClassContainer);

        void clearGappedThreading();
        std::vector<ThreadingStruct> & getGappedThreadingVector() { return gappedThreadingStructVector; }
        void createSpringsFromGappedThreading(BiopolymerClassContainer & myBiopolymerClassContainer);
        void addGappedThreading(const ThreadingStruct & threadingStruct, BiopolymerClassContainer & myBiopolymerClassContainer);
        void addGappedThreading(String chain1, String chain2, double forceConstant, bool backboneOnly, BiopolymerClassContainer & myBiopolymerClassContainer);
        void addGappedThreading(String chain1, ResidueID startResidue1,  ResidueID endResidue1,   String chain2, ResidueID startResidue2,  ResidueID endResidue2,  double forceConstant, bool backboneOnly,
                                       BiopolymerClassContainer & myBiopolymerClassContainer);
        void deleteGappedThreading(int id);
        void updateGappedThreading(int id, const ThreadingStruct & threadingStruct, BiopolymerClassContainer & myBiopolymerClassContainer);
        void updateGappedThreading(int id, String chain1, String chain2, double forceConstant, bool backboneOnly, BiopolymerClassContainer & myBiopolymerClassContainer);

        static ThreadingStruct createGappedThreading(String chain1, String chain2, double forceConstant, bool backboneOnly, BiopolymerClassContainer & myBiopolymerClassContainer1,  BiopolymerClassContainer & myBiopolymerClassContainer2);
        static ThreadingStruct createGappedThreading(String chain1, String chain2, double forceConstant, bool backboneOnly, BiopolymerClassContainer & myBiopolymerClassContainer) {return createGappedThreading( chain1,  chain2,  forceConstant, backboneOnly,  myBiopolymerClassContainer,   myBiopolymerClassContainer);}; // For backward compatibility. Just calls its two-BiopolymerClassContainer counterpart with both BiopolymerClassContainer's set to the same object.
        static ThreadingStruct createGappedThreading(String chain1, ResidueID startResidue1,  ResidueID endResidue1, String chain2,  ResidueID startResidue2,  ResidueID endResidue2,  double forceConstant, bool backboneOnly, BiopolymerClassContainer & myBiopolymerClassContainer1,BiopolymerClassContainer & myBiopolymerClassContainer2   ) ;
        static ThreadingStruct createGappedThreading(String chain1, ResidueID startResidue1,  ResidueID endResidue1, String chain2,  ResidueID startResidue2,  ResidueID endResidue2,  double forceConstant, bool backboneOnly, BiopolymerClassContainer & myBiopolymerClassContainer) {return createGappedThreading(chain1, startResidue1,  endResidue1, chain2,  startResidue2,   endResidue2,  forceConstant, backboneOnly,  myBiopolymerClassContainer ,  myBiopolymerClassContainer);}; // For backward compatibility. Just calls its two-BiopolymerClassContainer counterpart with both BiopolymerClassContainer's set to the same object.
        void printAllAlignmentStats(){for (int i = 0; i < threadingStructVector.size(); i++){threadingStructVector[i].printAlignmentStats();}}

};




#endif

