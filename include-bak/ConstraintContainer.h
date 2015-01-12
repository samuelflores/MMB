/* -------------------------------------------------------------------------- *
 *                           MMB (MacroMoleculeBuilder)                       *
 * -------------------------------------------------------------------------- *
 *                                                                            *
 * Copyright (c) 2011-12 by the Author.                                       *
 * Author: Samuel Flores                                                      *
 *                                                                            *
 * See RNABuilder.cpp for the copyright and usage agreement.                  *
 * -------------------------------------------------------------------------- */

#ifndef ConstraintContainer_H_
#define ConstraintContainer_H_
//#include "BiopolymerClass.h"
#include "Utils.h"
#include <vector>
#include <iostream>

class BiopolymerClassContainer; // trying a forward declaration

class MMB_EXPORT ConstraintToGroundContainer {
private:
	std::vector <ConstraintClass>      constraintClassVector;
	std::vector <ChainResidueToGround> constrainChainRigidSegmentsVector;
public:
  	void clear() {constraintClassVector.clear(); constrainChainRigidSegmentsVector .clear();};

	const ConstraintClass getConstraintClass(int constraintToGroundIndex) const  {
            return constraintClassVector[constraintToGroundIndex]; };

    std::vector<ConstraintClass> & getConstraintClassVector() { return constraintClassVector;}
	// Constrain rigid segments in chain to either ground or a given residue on same chain.  First step:  storing in constrainChainRigidSegmentsVector to queue them:
	void queueConstrainChainRigidSegments (ChainResidueToGround chainResidue){constrainChainRigidSegmentsVector.push_back(chainResidue); };
	// This will take the ChainResidueToGround's in constrainChainRigidSegmentsVector and apply them.  
	void applyConstrainChainRigidSegments (BiopolymerClassContainer & biopolymerClassContainer, CompoundSystem & system,  SimbodyMatterSubsystem & matter,State & state);

	void printConstraintClass(int constraintToGroundIndex) const;
	void validateConstraintClass(const ConstraintClass & myConstraintClass, BiopolymerClassContainer & myBiopolymerClassContainer) ;

					
	void validateConstraintClassVector(BiopolymerClassContainer & myBiopolymerClassContainer);
	void pruneCoordinateCouplers(BiopolymerClassContainer & myBiopolymerClassContainer, DuMMForceFieldSubsystem & _dumm);
	void addConstraintClassToVector(ConstraintClass myConstraintClass);
	void addConstraintClassToVector(String myChain, ResidueID myResidueID, String atomName);
    void addConstraintToVector(String myChain, ResidueID myResidueID, String atomName,
	                           String myChain2, ResidueID myResidueID2, String atomName2);

    void deleteConstraintClass(int index);
    void updateConstraintToVector(
                    int index,
                    String myChain, ResidueID myResidueID, String atomName,
                    BiopolymerClassContainer& biopolymerClassContainer);
    void updateConstraintToVector(
                    int index,
                    String myChain, ResidueID myResidueID, String atomName,
                    String myChain2, ResidueID myResidueID2, String atomName2,
                    BiopolymerClassContainer& biopolymerClassContainer);

	const int numConstraintClasses(){return constraintClassVector.size();};
	void printConstraintClasses();
    bool hasConstraintClass(String myChainID, ResidueID myResidueID);
};

#endif
