
#ifndef CifOutput_h
#define CifOutput_h

#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include "molmodel/internal/Compound.h"

#include "PeriodicPdbAndEnergyWriter.h"
#include "SimTKsimbody.h"

#ifdef GEMMI_USAGE
    #include <gemmi/align.hpp>
#endif

namespace SimTK
{
    namespace CIFOut
    {
#ifdef GEMMI_USAGE
        void assignEntities ( gemmi::Structure &outStruc, const map <const String, BiopolymerClass>& biopolymers, const CompoundSystem &system );
        void buildModel     ( const State& state, gemmi::Model& gModel, const map<const String, BiopolymerClass>& biopolymers, const CompoundSystem& system, int precision );
        void writeOutCif    ( const gemmi::Structure& outStruct, const std::string& fileName, const std::vector < std::pair < std::string, std::string > >& remarks );
        void reWriteOutCif  ( const gemmi::Model& gModel, const std::string& modelName, const std::string& fileName, ParameterReader& myParameterReader, const CompoundSystem& system, bool firstInStage );
#endif
    }                                                 // End namespace CIFOut
}                                                     // End namespace SimTK

#endif /* CifOutput_h */
