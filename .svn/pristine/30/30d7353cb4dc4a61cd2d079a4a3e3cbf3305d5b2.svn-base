/* -------------------------------------------------------------------------- *
 *                           MMB (MacroMoleculeBuilder)                       *
 * -------------------------------------------------------------------------- *
 *                                                                            *
 * Copyright (c) 2011-12 by the Author.                                       *
 * Author: Samuel Flores                                                      *
 *                                                                            *
 * See RNABuilder.cpp for the copyright and usage agreement.                  *
 * -------------------------------------------------------------------------- */

//#ifndef SimTK_MOLMODEL_PERIODICPDBWRITER_H_
//#define SimTK_MOLMODEL_PERIODICPDBWRITER_H_
#include <cstdio>
#include "SimTKsimbody.h"
#include "molmodel/internal/Compound.h"
//#include "DistanceMatrix.h"
#include <iostream>
#include <iomanip>
#include <vector>

namespace SimTK 
{

/// Writes atomic coordinates in PDB format to a file stream at
/// specified intervals during a simulation.
class PeriodicPdbAndEnergySingleFrameWriter : public PeriodicEventReporter {
public:
    PeriodicPdbAndEnergySingleFrameWriter(
        const CompoundSystem& system, 
        char* outputFrameFileName, //std::ostream& outputStream,
        Real interval) 
        : PeriodicEventReporter(interval), 
          system(system), 
        outputFrameFileName(outputFrameFileName) 
    {}

    void handleEvent(const State& state) const {
		static int modelNumber = 1; // increments by one at each reporting step
		int nextAtomSerialNumber = 1; // atom serial number for each compound picks up where previous compound left off

        system.realize(state, Stage::Position);	
	        std::ofstream outputStream;
		outputStream.open(  outputFrameFileName);
	        //std::ostream outputStream(outputFrameFileName);
		outputStream << "MODEL     " << std::setw(4) << modelNumber << std::endl;
		//cout<<"[PeriodicPdbAndEnergySingleFrameWriter.h] check 1"<<endl;
		for (SimTK::CompoundSystem::CompoundIndex c(0); c < system.getNumCompounds(); ++c)
			system.getCompound(c).writePdb(state, outputStream, nextAtomSerialNumber);

		outputStream << "ENDMDL" << std::endl;
	
		//scf added time reporting 
                time_t rawtime;
                struct tm * timeinfo;
                time ( &rawtime );
                timeinfo = localtime ( &rawtime );
                outputStream <<"REMARK Current time is: "<<asctime (timeinfo) <<"REMARK elapsed time: "<<(clock()/CLOCKS_PER_SEC)<<std::endl;
                outputStream.setf(ios::fixed, ios::floatfield); // set output to fixed rather than scientific format
                outputStream <<"REMARK Energy = "<<system.calcEnergy(state) <<std::endl;
	        outputStream.close();
		++modelNumber;
    }

private:
    const CompoundSystem& system;
    //std::ostream& outputStream;
    char* outputFrameFileName;// outputStream(outputStream) 
};

} // namespace SimTK

//#endif // SimTK_MOLMODEL_PERIODICPDBWRITER_H_
