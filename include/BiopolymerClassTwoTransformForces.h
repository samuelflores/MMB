
/* -------------------------------------------------------------------------- *
 *                           MMB (MacroMoleculeBuilder)                       *
 * -------------------------------------------------------------------------- *
 *                                                                            *
 * Copyright (c) 2011-12 by the Author.                                       *
 * Author: Samuel Flores                                                      *
 *                                                                            *
 * See RNABuilder.cpp for the copyright and usage agreement.                  *
 * -------------------------------------------------------------------------- */

//#include "SimTKsimbody_aux.h"
#include <iostream>
#include <ostream>
#include <stdlib.h>// for MAX_PATH
#include "SimTKmolmodel.h"
#include "ParameterReader.h"
#include "BiopolymerClass.h"
using namespace SimTK;
using namespace std; 


inline std::ostream& operator<<(std::ostream& o, const ParameterReader&) {
    assert(false);
    return o;
};


class AllTwoTransformLinearSprings : public Force::Custom::Implementation { 

private: 
    SimbodyMatterSubsystem& matter;
    ParameterReader& myParameterReader;	 
    LeontisWesthofClass& myLeontisWesthofClass;
    BiopolymerClassContainer & myBiopolymerClassContainer;
    //Biopolymer * myChain;
    mutable int parameterReaderIndex; 
    std::ostream& outputStream; 
public: 

    AllTwoTransformLinearSprings (SimbodyMatterSubsystem& matter,ParameterReader& myParameterReader,  LeontisWesthofClass& myLeontisWesthofClass, BiopolymerClassContainer & myBiopolymerClassContainer, std::ostream& outputStream ) ; // : matter(matter),myParameterReader(myParameterReader), myLeontisWesthofClass (myLeontisWesthofClass), myBiopolymerClassContainer(myBiopolymerClassContainer), outputStream(outputStream);

    void         calcAxes (const State& state ,LeontisWesthofBondRow myLeontisWesthofBondRow,ResidueID residueNumber1,ResidueID residueNumber2,String,String,Vec3 & xAxisVector1,Vec3 & yAxisVector1, Vec3 & zAxisVector1,Vec3 & xAxisVector2,Vec3 & yAxisVector2 , Vec3 & zAxisVector2,Vec3 & glycosidicNitrogenAtom1LocationInGround,Vec3 & glycosidicNitrogenAtom2LocationInGround, Vec3 & ring1CenterLocationInGround, Vec3 & ring2CenterLocationInGround) const ;  
    int isThisATwoTransformForce(String myBPEdge) const;
    void calcForce(const State& state, Vector_<SpatialVec>& bodyForces,  
            Vector_<Vec3>& particleForces, Vector& mobilityForces) const ; 
    Real calcPotentialEnergy(const State& state) const; 
    bool dependsOnlyOnPositions() const;
};
