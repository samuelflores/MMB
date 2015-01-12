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

#ifndef SimTK_MOLMODEL_PERIODICPDBANDCOMWRITER_H_
#define SimTK_MOLMODEL_PERIODICPDBANDCOMWRITER_H_

#include "SimTKsimbody.h"
#include "molmodel/internal/Compound.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <list>
#include <stdio.h>
#include "RNATools/include/CenterOfMass.h"

namespace SimTK {

/// Writes atomic coordinates in PDB format to a file stream at
/// specified intervals during a simulation.

class PeriodicPdbAndCOMWriter : public PeriodicEventReporter {
public:
    PeriodicPdbAndCOMWriter(
        const CompoundSystem& system, 
        std::ostream& outputStream,
        Real interval,
        SimbodyMatterSubsystem & matter,//(system)
        Biopolymer & myMolecule,
        TinkerDuMMForceFieldSubsystem & dumm,
        list<int> residueList 

	) 
        : PeriodicEventReporter(interval), 
          system(system), 
          outputStream(outputStream)
          ,matter(matter),
          myMolecule(myMolecule), 
          dumm(dumm),
          residueList(residueList)
    {}
    void handleEvent(const State& state) const {
		static int modelNumber = 1; // increments by one at each reporting step
		int nextAtomSerialNumber = 1; // atom serial number for each compound picks up where previous compound left off

                system.realize(state, Stage::Position);

		outputStream << "MODEL     " << std::setw(4) << modelNumber << std::endl;

		for (SimTK::Compound::Index c(0); c < system.getNumCompounds(); ++c)
			system.getCompound(c).writePdb(state, outputStream, nextAtomSerialNumber);

	
		//scf added time reporting 
                time_t rawtime;
                struct tm * timeinfo;
                time ( &rawtime );
                timeinfo = localtime ( &rawtime );
                outputStream <<"REMARK Current time is: "<<asctime (timeinfo) <<"REMARK elapsed time: "<<(clock()/CLOCKS_PER_SEC)<<std::endl;
		//
                //scf added Geometric Center reporting
                Vec3 helix2GC = CenterOfMass ( myMolecule , residueList, matter,  state,dumm);
                cout<<"REMARK Geometric Center ="<<helix2GC<<endl;
                cout<<"REMARK Energy ="<<system.calcEnergy(state)<<endl;
                printf (        "ATOM      1  GC  GC  C   1    %8.3f%8.3f%8.3f",helix2GC[0]*10,helix2GC[1]*10,helix2GC[2]*10); 
                char buffer [70];
                sprintf (buffer,"ATOM      1  GC  GC  C   1    %8.3f%8.3f%8.3f",helix2GC[0]*10,helix2GC[1]*10,helix2GC[2]*10); 
                outputStream<<(string(buffer))<<endl; 
 	        //
		outputStream << "ENDMDL" << std::endl;


		++modelNumber;
    }

private:
    const CompoundSystem& system;
    std::ostream& outputStream;
    SimbodyMatterSubsystem & matter;
    Biopolymer & myMolecule; 
    TinkerDuMMForceFieldSubsystem & dumm;
    list<int> residueList ;
};
/*
*/
} // namespace SimTK

#endif // SimTK_MOLMODEL_PERIODICPDBWRITER_H_
