
/* -------------------------------------------------------------------------- *
 *                           MMB (MacroMoleculeBuilder)                       *
 * -------------------------------------------------------------------------- *
 *                                                                            *
 * Copyright (c) 2011-12 by the Author.                                       *
 * Author: Samuel Flores                                                      *
 *                                                                            *
 * See RNABuilder.cpp for the copyright and usage agreement.                  *
 * -------------------------------------------------------------------------- */

#include <math.h>
//#include "SimTKsimbody_aux.h"
#include <iostream>
#include <ostream>
#include <stdlib.h>// for MAX_PATH
#include <string.h>
#include "SimTKmolmodel.h"
#include "ParameterReader.h"
#include "BiopolymerClass.h"
#include "Utils.h"
using namespace SimTK;
using namespace std; 


//inline std::ostream& operator<<(std::ostream& o, const ParameterReader&) {
//    assert(false);
//    return o;
//};


class NTC_Torque : public Force::Custom::Implementation { 

private: 
    SimbodyMatterSubsystem& matter;
    ParameterReader& myParameterReader;	 
    NTC_PAR_Class& myNTC_PAR_Class;
    BiopolymerClassContainer & myBiopolymerClassContainer;
    //Biopolymer * myChain;
    mutable int parameterReaderIndex; 
    std::ostream& outputStream; 
public: 

    NTC_Torque (SimbodyMatterSubsystem& matter,ParameterReader& myParameterReader,  NTC_PAR_Class& myNTC_PAR_Class, BiopolymerClassContainer & myBiopolymerClassContainer, std::ostream& outputStream ) ; // : matter(matter),myParameterReader(myParameterReader), myNTC_PAR_Class (myNTC_PAR_Class), myBiopolymerClassContainer(myBiopolymerClassContainer), outputStream(outputStream);

    void  calcAxes (const State& state ,NTC_PAR_BondRow myNTC_PAR_BondRow,ResidueID residueNumber1,ResidueID residueNumber2,String chain1,String chain2,Vec3 & xAxisVector1,Vec3 & yAxisVector1, Vec3 & zAxisVector1,Vec3 & xAxisVector2,Vec3 & yAxisVector2 , Vec3 & zAxisVector2,
        Vec3 & glycosidicNitrogenAtom1LocationInGround,Vec3 & glycosidicNitrogenAtom2LocationInGround, Vec3 & ring1CenterLocationInGround, Vec3 & ring2CenterLocationInGround) const ;  
    int isThisATwoTransformForce(String myBPEdge) const;
    void calcForce(const State& state, Vector_<SpatialVec>& bodyForces,  
            Vector_<Vec3>& particleForces, Vector& mobilityForces) const ; 
    Real calcPotentialEnergy(const State& state) const; 
    bool dependsOnlyOnPositions() const;
    
    Real return_dist_ang(double angle,double rotationAngle) const;
    Real return_angle(Vec3 cross_1,Vec3 cross_2,Vec3 cross_3,Vec3 d_d2) const; 
};
