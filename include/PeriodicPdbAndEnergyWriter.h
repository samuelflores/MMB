/* -------------------------------------------------------------------------- *
 *                           MMB (MacroMoleculeBuilder)                       *
 * -------------------------------------------------------------------------- *
 *                                                                            *
 * Copyright (c) 2011-12 by the Author.                                       *
 * Author: Samuel Flores                                                      *
 *                                                                            *
 * See RNABuilder.cpp for the copyright and usage agreement.                  *
 * -------------------------------------------------------------------------- */

#ifndef SimTK_MOLMODEL_PERIODICPDBANDENERGYWRITER_H_
#define SimTK_MOLMODEL_PERIODICPDBANDENERGYWRITER_H_
#include "ParameterReader.h"
#include "molmodel/internal/Compound.h"
#include "BiopolymerClass.h"
#include "DensityForce.h"
//#include "DistanceMatrix.h"
#include <vector>

namespace SimTK 
{

/// Writes atomic coordinates in PDB format to a file stream at
/// specified intervals during a simulation.
class PeriodicPdbAndEnergyWriter : public PeriodicEventHandler {
public:
    PeriodicPdbAndEnergyWriter(
        const CompoundSystem& system, 
        const DuMMForceFieldSubsystem& dumm  , 
        const DensityForce           & densityForce,
        std::ostream& outputStream,
        double interval,
        ParameterReader & myParameterReader,
        BiopolymerClassContainer & myBiopolymerClassContainer
        //vector<MagnesiumIon> myMagnesiumIonVec

    );
    
    void handleEvent(State& state, Real accuracy, bool& shouldTerminate) const;

private:
    const CompoundSystem& system;
    const DuMMForceFieldSubsystem& dumm        ; 
    const DensityForce           & densityForce; 
    std::ostream& outputStream;
    ParameterReader & myParameterReader;
    BiopolymerClassContainer & myBiopolymerClassContainer;
    static std::vector<double> myEnergies;

};

} // namespace SimTK

#endif // SimTK_MOLMODEL_PERIODICPDBANDENERGYWRITER_H_
