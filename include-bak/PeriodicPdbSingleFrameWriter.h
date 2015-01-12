#ifndef PeriodicPdbSingleFrameWriter_H_
#define PeriodicPdbSingleFrameWriter_H_
/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Molmodel                            *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
 * Authors: Christopher Bruns                                                 *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

//#ifndef SimTK_MOLMODEL_PERIODICPDBWRITER_H_
//#define SimTK_MOLMODEL_PERIODICPDBWRITER_H_
#include <cstdio>
#include "SimTKsimbody.h"
#include "molmodel/internal/Compound.h"
//#include "DistanceMatrix.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
//using namespace std;
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
                //outputStream <<"REMARK Energy = "<<system.calcEnergy(state) <<std::endl;
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
#endif
