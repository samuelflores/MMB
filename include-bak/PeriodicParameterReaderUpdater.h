#ifndef PeriodicParameterReaderUpdater_H_
#define PeriodicParameterReaderUpdater_H_
/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Molmodel                            *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
 * Authors: Samuel Flores                                                     *
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
#include <iomanip>
#include <vector>

namespace SimTK 
{

/// Writes atomic coordinates in PDB format to a file stream at
/// specified intervals during a simulation.
class PeriodicParameterReaderUpdater : public PeriodicEventHandler  {
public:
    PeriodicParameterReaderUpdater(
        const CompoundSystem& system, 
        AllTwoTransformLinearSprings & myAllTwoTransformLinearSprings,
        Real interval) 
        : PeriodicEventHandler (interval), 
          system(system), myAllTwoTransformLinearSprings(myAllTwoTransformLinearSprings), interval(interval)
    {}

    void handleEvent( State& state ,Real accuracy, const Vector& yWeights, const Vector& 
        ooConstraintTols, Stage& lowestModified, bool& shouldTerminate )  const {

        cout<<"[PeriodicParameterReaderUpdater.h] updating ParameterReader object"<<endl; 
 //       cout<<"[PeriodicParameterReaderUpdater.h] interval"<<interval<<endl;//updating ParameterReader object"<<endl; 
        ParameterReader myParameterReader = myAllTwoTransformLinearSprings.getParameterReader( state);
        myParameterReader.initializeDefaults();
	myParameterReader.initializeFromFileOnly("commands.dat");        
	myParameterReader.initializeFromFileOnly("output.txt"  );        
        myParameterReader.postInitialize();
	assert(myParameterReader.currentStage>0);
        myParameterReader.removeNonPriorityBasePairs(myParameterReader.currentStage);
        myParameterReader.printBasePairs();
        myAllTwoTransformLinearSprings.setParameterReader( state,myParameterReader);
/*
*/
    }
Real getNextEventTime(const State& state, bool includeCurrentTime) const { 
    Real currentTime = state.getTime();
    long count = (long)std::floor(currentTime/interval);// impl->eventInterval);
    //long count = (long)std::floor(currentTime/ impl->eventInterval);
    volatile Real eventTime = count*interval;//mpl->eventInterval;
    //volatile Real eventTime = count*impl->eventInterval;
    while (eventTime < currentTime || (eventTime == currentTime && !includeCurrentTime)) {
        count++;
        eventTime = count*interval;//mpl->eventInterval;
        //eventTime = count*impl->eventInterval;
    }   
    return eventTime;


} 
private:
    const CompoundSystem& system;
    AllTwoTransformLinearSprings & myAllTwoTransformLinearSprings;
    //std::ostream& outputStream;
    Real interval;
    //ParameterReader myParameterReader;
};

} // namespace SimTK

//#endif // SimTK_MOLMODEL_PERIODICPDBWRITER_H_
#endif
