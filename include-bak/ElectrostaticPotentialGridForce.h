/* -------------------------------------------------------------------------- *
 *                           MMB (MacroMoleculeBuilder)                       *
 * -------------------------------------------------------------------------- *
 *                                                                            *
 * Copyright (c) 2011-12 by the Author.                                       *
 * Author: Samuel Flores                                                      *
 *                                                                            *
 * See RNABuilder.cpp for the copyright and usage agreement.                  *
 * -------------------------------------------------------------------------- */

#ifndef ElectrostaticPotentialGridForce_H_
#define ElectrostaticPotentialGridForce_H_

#include <math.h>
//#include "SimTKsimbody_aux.h"
#include <iostream>
#include <ostream>
#include <stdlib.h>// for MAX_PATH
#include "SimTKmolmodel.h"
#include "ParameterReader.h"
#include "BiopolymerClass.h"
#include "DensityMap.h"
using namespace SimTK;
using namespace std; 


/*inline std::ostream& operator<<(std::ostream& o, const ParameterReader&) {
    assert(false);
    return o;
};*/


class ElectrostaticPotentialGridForce : public Force::Custom::Implementation { 

protected: 
    SimbodyMatterSubsystem& matter;
   ParameterReader& myParameterReader;   
   DensityMap     & myDensityMap     ;
   DuMMForceFieldSubsystem & dumm;   
    //LeontisWesthofClass& myLeontisWesthofClass;
    BiopolymerClassContainer & myBiopolymerClassContainer;
    Biopolymer * myChain;
    mutable int parameterReaderIndex; 
    std::ostream& outputStream; 
public: 

    ElectrostaticPotentialGridForce (SimbodyMatterSubsystem& matter,ParameterReader& myParameterReader, DensityMap & myDensityMap, DuMMForceFieldSubsystem & dumm,  BiopolymerClassContainer & myBiopolymerClassContainer, std::ostream& outputStream ) ; 

    void calcForce(const State& state, Vector_<SpatialVec>& bodyForces,  
            Vector_<Vec3>& particleForces, Vector& mobilityForces) const ; 
    Real calcPotentialEnergy(const State& state) const; 
    bool dependsOnlyOnPositions() const;
};
#endif