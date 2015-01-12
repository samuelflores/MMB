/* -------------------------------------------------------------------------- *
 *                           MMB (MacroMoleculeBuilder)                       *
 * -------------------------------------------------------------------------- *
 *                                                                            *
 * Copyright (c) 2011-12 by the Author.                                       *
 * Author: Samuel Flores                                                      *
 *                                                                            *
 * See RNABuilder.cpp for the copyright and usage agreement.                  *
 * -------------------------------------------------------------------------- */
#ifndef _TetherForce_H_
#define _TetherForce_H_


#include <math.h>
//#include "SimTKsimbody_aux.h"
#include <iostream>
#include <ostream>
#include <stdlib.h>// for MAX_PATH
#include "SimTKmolmodel.h"
#include "ParameterReader.h"
using namespace SimTK;
using namespace std; 

class TetherForce : public Force::Custom::Implementation { 
public: 

    TetherForce (
                SimbodyMatterSubsystem &   	matter, 
                GeneralForceSubsystem  &   	 forces,
		const MobilizedBodyIndex   	body1Index,
		const Vec3  & 	station1,
		const MobilizedBodyIndex   	body2Index,
		const Vec3  & 	station2,
		Real  	k,
		Real  	x0
                );

    void calcForce(const State& state, Vector_<SpatialVec>& bodyForces,  
            Vector_<Vec3>& particleForces, Vector& mobilityForces) const ; 
    Real calcPotentialEnergy(const State& state) const; 
    bool dependsOnlyOnPositions() const;
private:
    
                SimbodyMatterSubsystem &   	matter; 
                GeneralForceSubsystem  & 	forces; 
		const MobilizedBodyIndex     	body1Index;
		const Vec3               	station1;
		const MobilizedBodyIndex     	body2Index;
		const Vec3               	station2;
		Real  	                 	k;
		Real                     	x0;
};

#endif
