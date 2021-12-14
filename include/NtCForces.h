
/* -------------------------------------------------------------------------- *
 *                           MMB (MacroMoleculeBuilder)                       *
 * -------------------------------------------------------------------------- *
 *                                                                            *
 * Copyright (c) 2011-12 by the Author.                                       *
 * Author: Samuel Flores                                                      *
 *                                                                            *
 * See RNABuilder.cpp for the copyright and usage agreement.                  *
 * -------------------------------------------------------------------------- */

#include "SimTKmolmodel.h"
#include "ParameterReader.h"
#include "BiopolymerClass.h"
#include "Utils.h"

using namespace SimTK;
using namespace std;

class NTC_Torque : public Force::Custom::Implementation {
private:
    SimbodyMatterSubsystem& matter;
    ParameterReader& myParameterReader;
    NTC_PAR_Class& myNTC_PAR_Class;
    BiopolymerClassContainer & myBiopolymerClassContainer;

public:
    NTC_Torque(SimbodyMatterSubsystem &matter,
               ParameterReader &myParameterReader,
               NTC_PAR_Class &myNTC_PAR_Class,
               BiopolymerClassContainer &myBiopolymerClassContainer);

    void calcForce(const State &state, Vector_<SpatialVec> &bodyForces, Vector_<Vec3> &particleForces, Vector &mobilityForces) const;
    Real calcPotentialEnergy(const State& state) const;
    bool dependsOnlyOnPositions() const;
    Real return_dist_ang(double angle,double rotationAngle) const;
    Real return_angle(const Vec3 &cross_1, const Vec3 &cross_2, const Vec3 &cross_3, const Vec3 &d_d2) const;
};
