
#ifndef CifOutput_h
#define CifOutput_h

#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include "molmodel/internal/Compound.h"


namespace SimTK
{
    namespace CIFOut
    {
#ifdef GEMMI_USAGE
        void writeOutCif ( gemmi::Structure outStruct, std::string fileName, std::vector < std::pair < std::string, std::string > > remarks );
#endif
    }                                                 // End namespace CIFOut
}                                                     // End namespace SimTK

#endif /* CifOutput_h */
