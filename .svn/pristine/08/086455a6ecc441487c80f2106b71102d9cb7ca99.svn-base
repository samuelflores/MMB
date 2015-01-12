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
//#include "ParameterReader.h"
#include "molmodel/internal/Compound.h"
//#include "DistanceMatrix.h"
#include <iostream>
#include <iomanip>
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
        std::ostream& outputStream,
        Real interval,
        ParameterReader & myParameterReader//,
        //vector<MagnesiumIon> myMagnesiumIonVec

        ) 
        : PeriodicEventHandler(interval), 
          system(system), 
      dumm(dumm), // the second dumm comes from the parameters fed to the constructor.  the first is the private member.
          outputStream(outputStream) ,
          myParameterReader(myParameterReader)//,
          //myMagnesiumIonVec(myMagnesiumIonVec)
    {}
void handleEvent(State& state, Real accuracy, bool& shouldTerminate) const  {

        static int modelNumber = 1; // increments by one at each reporting step

        system.realize(state, Stage::Dynamics);
        outputStream << "MODEL     " << std::setw(4) << modelNumber << std::endl;
        PdbAtom::setWriteFullPrecisionLocation(myParameterReader.writeDoublePrecisionTrajectories);// false by default to save a little disk space here
        for (SimTK::CompoundSystem::CompoundIndex c(0); c < system.getNumCompounds(); ++c)
            (system.getCompound(c)).writePdb(state, outputStream,Transform(Vec3(0)));//, nextAtomSerialNumber);
        filebuf fb;
        fb.open("frame.pdb",ios::out);
        std::ostream  fbstream (&fb);
        PdbAtom::setWriteFullPrecisionLocation(true); // get higher precision from file that might be reused for reading.
        for (SimTK::CompoundSystem::CompoundIndex c(0); c < system.getNumCompounds(); ++c)
            (system.getCompound(c)).writePdb(state, fbstream,Transform(Vec3(0)));

        //scf added time reporting 
        time_t rawtime;
        struct tm * timeinfo;
        time ( &rawtime );
        timeinfo = localtime ( &rawtime );

        outputStream << "ENDMDL" << std::endl;
                     
        outputStream <<"REMARK seconds since January 1st 1970: "<<time ( NULL     )<<std::endl; //<<"REMARK elapsed time: "<<(clock()/CLOCKS_PER_SEC)<<std::endl;
        outputStream <<"REMARK Current time is: "<<asctime (timeinfo) <<"REMARK elapsed time: "<<(clock()/CLOCKS_PER_SEC)<<std::endl;
        outputStream.setf(ios::fixed, ios::floatfield); // set output to fixed rather than scientific format
        if (myParameterReader.calcEnergy) {
            double myPotentialEnergy = system.calcPotentialEnergy(state);
            myParameterReader.potentialEnergy = myPotentialEnergy;
            outputStream <<"REMARK Potential Energy = "<<myPotentialEnergy <<" kJ/mol, "<<myPotentialEnergy/4.184<<" kcal/mol"<<std::endl;
            cout <<"REMARK Potential Energy = "<<myPotentialEnergy <<" kJ/mol, "<<myPotentialEnergy/4.184<<" kcal/mol"<<std::endl;

            double myKineticEnergy = system.calcKineticEnergy(state);
            myParameterReader.kineticEnergy = myKineticEnergy;
            outputStream <<"REMARK Kinetic Energy = "<< myKineticEnergy <<" kJ/mol, "<<myKineticEnergy/4.184<<" kcal/mol"<<std::endl;
            cout <<"REMARK Kinetic Energy = "<< myKineticEnergy <<" kJ/mol, "<<myKineticEnergy/4.184<<" kcal/mol"<<std::endl;

            double myEnergy = system.calcEnergy(state);
            myParameterReader.totalEnergy = myEnergy;
            outputStream <<"REMARK Energy = "<< myEnergy <<" kJ/mol "<<std::endl;
            cout <<"REMARK Energy = "<< myEnergy <<" kJ/mol "<<std::endl;
            myEnergies.push_back(myEnergy);
        }
        outputStream <<"REMARK Angular, Linear Momentum = "<<system.calcSystemRigidBodyMomentum(state)<<endl;

        //cout<<__FILE__<<" : "<<__LINE__<<" "<<dumm.getVdwMixingRuleName (dumm.getVdwMixingRule())<<endl;

        outputStream <<"REMARK ["<< __FILE__<<"] state.getNU()    = "<<state.getNU()            <<std::endl;
        //outputStream <<"REMARK ["<< __FILE__<<"]Satisfied contacts : "<<myParameterReader.satisfiedBasePairs<<endl;
        //outputStream <<"REMARK ["<< __FILE__<<"]Unsatisfied contacts : "<<myParameterReader.unSatisfiedBasePairs<<endl;

        cout<<"Just wrote structure for reporting interval # "<<modelNumber<<std::endl; 
        //cout <<"Satisfied base pairs : "<<myParameterReader.satisfiedBasePairs<<" out of "<<myParameterReader.satisfiedBasePairs+myParameterReader.unSatisfiedBasePairs<<endl;
        //cout <<"Unsatisfied contacts : "<<myParameterReader.unSatisfiedBasePairs<<endl;
        ++modelNumber;

        // Check if converged or not
        if (myParameterReader.detectConvergence) {
            if(myEnergies.size() > myParameterReader.convergenceTimeout) {
                double lastEnergy = myEnergies.back();
                double energyDiffMean = 0.0;
                vector<double>::iterator it;
                for(it = myEnergies.end()-myParameterReader.convergenceTimeout; it != myEnergies.end(); it++)
                {
                    double diffEnergy = fabs(*it - *(it-1));
                    // cout << diffEnergy << " " << endl;
                    energyDiffMean += diffEnergy;
                }
                energyDiffMean /= myParameterReader.convergenceTimeout;
                // cout << "Last " << myParameterReader.convergenceTimeout << " energies diff mean: " << energyDiffMean << endl;
                if(energyDiffMean < myParameterReader.convergenceEpsilon)
                {
                    myParameterReader.converged = true;   
                    cout << "Converged! Energy difference between two reporting intervals has been less than "<< myParameterReader.convergenceEpsilon << " kJ/mol for the last " << myParameterReader.convergenceTimeout << " frames." << endl;
                }
            }
        }
    }

private:
    const CompoundSystem& system;
    const DuMMForceFieldSubsystem& dumm  ; 
    std::ostream& outputStream;
    ParameterReader & myParameterReader;
    static std::vector<double> myEnergies; 

};

std::vector<double> PeriodicPdbAndEnergyWriter::myEnergies;

} // namespace SimTK

//#endif // SimTK_MOLMODEL_PERIODICPDBWRITER_H_
