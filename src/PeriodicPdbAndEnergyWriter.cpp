/* -------------------------------------------------------------------------- *
 *                           MMB (MacroMoleculeBuilder)                       *
 * -------------------------------------------------------------------------- *
 *                                                                            *
 * Copyright (c) 2011-12 by the Author.                                       *
 * Author: Samuel Flores                                                      *
 *         Alex Tek                                                           *  
 *                                                                            *
 * See RNABuilder.cpp for the copyright and usage agreement.                  *
 * -------------------------------------------------------------------------- */
#include <cstdio>
#include <iostream>
#include <iomanip>

#include <cassert>

#include "MMBLogger.h"
#include "PeriodicPdbAndEnergyWriter.h"
#include "SimTKsimbody.h"
#include "CifOutput.h"
#include "ProgressWriter.h"

using namespace std;

std::vector<double> PeriodicPdbAndEnergyWriter::myEnergies;

SimTK::PeriodicPdbAndEnergyWriter::PeriodicPdbAndEnergyWriter( 	const CompoundSystem& system,
    															const DuMMForceFieldSubsystem& dumm  , 
    															std::ostream& outputStream,
    															double interval,
    															ParameterReader & myParameterReader,
    															BiopolymerClassContainer & myBiopolymerClassContainer
    															//vector<MagnesiumIon> myMagnesiumIonVec
    														 ) : PeriodicEventHandler(interval),
																 system(system), 
																 dumm(dumm), // the second dumm comes from the parameters fed to the constructor.  the first is the private member.
																 outputStream(outputStream) ,
																 myParameterReader(myParameterReader),
                                                                                                                                 myBiopolymerClassContainer(myBiopolymerClassContainer)
																 //myMagnesiumIonVec(myMagnesiumIonVec)
{
    myParameterReader.trajectoryFileRemarks.clear     ( );
    myParameterReader.gemmi_isFirstInStage            = true;
}

void SimTK::PeriodicPdbAndEnergyWriter::handleEvent(State& state, Real accuracy, bool& shouldTerminate) const  {

    static int modelNumber = 1; // increments by one at each reporting step
    int compoundNumber                                = 1;
    int requiredPrecision                             = 12;
    
    system.realize(state, Stage::Dynamics);

    //================================================ Create new gemmi model for this position
    gemmi::Model gModel                               ( std::to_string( modelNumber ) );
    
    //================================================ Solve the data name
    std::string strName                               = myParameterReader.outTrajectoryFileName;
    istringstream iss                                 ( strName );
    std::vector<std::string> tokens; std::string token;
    while ( std::getline ( iss, token, '.') )         { if ( !token.empty() ) { tokens.push_back ( token ); } }
    if ( tokens.size() > 2 )                          { strName = std::string ( tokens.at(tokens.size()-3) + tokens.at(tokens.size()-2) ); }
    else                                              { strName = "TRAJECTORYX"; }
    if ( strName[0] == '/' )                          { strName.erase(0, 1); }
    
    if (! myParameterReader.useCIFFileFormat )
    {
        outputStream << "MODEL     " << std::setw(4) << modelNumber << std::endl;
    }
    PdbAtom::setWriteFullPrecisionLocation(myParameterReader.writeDoublePrecisionTrajectories);// false by default to save a little disk space here

    if (myParameterReader.contactInterfaceContainer.numInterfaces() > 0) {
	vector<MMBAtomInfo> concatenatedAtomInfoVector = myBiopolymerClassContainer.getConcatenatedAtomInfoVector(state);
	vector<TwoAtomClass> myTwoAtomClassVector = myParameterReader.contactInterfaceContainer.retrieveCloseContactPairs(concatenatedAtomInfoVector); // Just calling this method is enough to get the close contacts written to stdout.
    }

    if ( myParameterReader.useCIFFileFormat )
    {
        const CIFOut::Data data{myParameterReader.myBiopolymerClassContainer.getBiopolymerClassMap(), myParameterReader.myMonoAtomsContainer};
        CIFOut::buildModel( state, gModel, data, system, 3 );
    }
    else
    {
        for (SimTK::CompoundSystem::CompoundIndex c(0); c < system.getNumCompounds(); ++c)
        {
            (system.getCompound(c)).writePdb(state, outputStream,Transform(Vec3(0)));//, nextAtomSerialNumber);
        }
    }

    if ( myParameterReader.useCIFFileFormat )
    {
        SimTK::CIFOut::reWriteOutCif                  ( gModel, strName, "frame.cif", myParameterReader, system, true );
    }
    else
    {
        filebuf fb;
        fb.open("frame.pdb",ios::out);
        std::ostream  fbstream (&fb);
        PdbAtom::setWriteFullPrecisionLocation(true); // get higher precision from file that might be reused for reading.
        for (SimTK::CompoundSystem::CompoundIndex c(0); c < system.getNumCompounds(); ++c)
            (system.getCompound(c)).writePdb(state, fbstream,Transform(Vec3(0)));
    }

    //scf added time reporting 
    time_t rawtime;
    struct tm * timeinfo;
    time ( &rawtime );
    timeinfo = localtime ( &rawtime );

    if ( myParameterReader.useCIFFileFormat )
    {
        //============================================ Generate trajectory remarks
        myParameterReader.trajectoryFileRemarks.push_back ( std::pair < std::string, std::string > ( "3", "Trajectory " + std::to_string ( modelNumber ) + ": Seconds since January 1st 1970           : " + std::to_string ( time ( NULL ) ) ) );
        std::string curTimeHlp                        ( asctime (timeinfo) );
        curTimeHlp.erase                              ( std::remove ( curTimeHlp.begin(), curTimeHlp.end(), '\n' ), curTimeHlp.end() );
        myParameterReader.trajectoryFileRemarks.push_back ( std::pair < std::string, std::string > ( "3", "Trajectory " + std::to_string ( modelNumber ) + ": Current time is                          : " + curTimeHlp ) );
        myParameterReader.trajectoryFileRemarks.push_back ( std::pair < std::string, std::string > ( "3", "Trajectory " + std::to_string ( modelNumber ) + ": Elapsed time                             : " + std::to_string ( clock()/CLOCKS_PER_SEC ) ) );
        
        if (myParameterReader.calcEnergy)
        {
            double myPotentialEnergy                  = system.calcPotentialEnergy(state);
            Real dummPotentialEnergy                  = dumm.calcPotentialEnergy(state);
            double myKineticEnergy                    = system.calcKineticEnergy(state);
            double myEnergy                           = system.calcEnergy(state);
            myParameterReader.kineticEnergy           = myKineticEnergy;
            myParameterReader.totalEnergy             = myEnergy;
            myEnergies.push_back                      ( myEnergy );
            
            myParameterReader.trajectoryFileRemarks.push_back ( std::pair < std::string, std::string > ( "3", "Trajectory " + std::to_string ( modelNumber ) + ": Total Potential Energy                   = " + std::to_string ( myPotentialEnergy ) + " kJ/mol, " + std::to_string( myPotentialEnergy / 4.184 ) + " kcal/mol" ) );
            myParameterReader.trajectoryFileRemarks.push_back ( std::pair < std::string, std::string > ( "3", "Trajectory " + std::to_string ( modelNumber ) + ": DuMM Molecular Dynamics Potential Energy = " + std::to_string ( dummPotentialEnergy ) + " kJ/mol, " + std::to_string( dummPotentialEnergy / 4.184 ) + " kcal/mol" ) );
            myParameterReader.trajectoryFileRemarks.push_back ( std::pair < std::string, std::string > ( "3", "Trajectory " + std::to_string ( modelNumber ) + ": Kinetic Energy                           = " + std::to_string ( myKineticEnergy ) + " kJ/mol, " + std::to_string( myKineticEnergy / 4.184 ) + " kcal/mol" ) );
            myParameterReader.trajectoryFileRemarks.push_back ( std::pair < std::string, std::string > ( "3", "Trajectory " + std::to_string ( modelNumber ) + ": Energy                                   = " + std::to_string ( myEnergy ) + " kJ/mol" ) );

            MMBLOG_PLAIN_NOSEV(ALWAYS, "Total Potential Energy =                   " << myPotentialEnergy << " kJ/mol, " << myPotentialEnergy/4.184 << " kcal/mol" << std::endl);
            MMBLOG_PLAIN_NOSEV(ALWAYS, "DuMM Molecular Dynamics Potential Energy = " << dummPotentialEnergy << " kJ/mol, " << dummPotentialEnergy/4.184 << " kcal/mol"<<std::endl);
            MMBLOG_PLAIN_NOSEV(ALWAYS, "Kinetic Energy =  " << myKineticEnergy << " kJ/mol, " << myKineticEnergy/4.184 << " kcal/mol" << std::endl);
            MMBLOG_PLAIN_NOSEV(ALWAYS, "REMARK Energy = " << myEnergy << " kJ/mol " << std::endl);
        }
        
        std::stringstream alMomHlp; alMomHlp << system.calcSystemRigidBodyMomentum(state);
        myParameterReader.trajectoryFileRemarks.push_back ( std::pair < std::string, std::string > ( "3", "Trajectory " + std::to_string ( modelNumber ) + ": Angular, Linear Momentum                 = " + alMomHlp.str() ) );
        std::stringstream getNuHlp; getNuHlp << state.getNU();
        myParameterReader.trajectoryFileRemarks.push_back ( std::pair < std::string, std::string > ( "3", "Trajectory " + std::to_string ( modelNumber ) + ": [" + __FILE__ + "] state.getNU()    = " + getNuHlp.str() ) );
        myParameterReader.trajectoryFileRemarks.push_back ( std::pair < std::string, std::string > ( "3", "" ) );
        
        //============================================ Write out CIF
        SimTK::CIFOut::reWriteOutCif                  ( gModel, strName, myParameterReader.outTrajectoryFileName, myParameterReader, system, myParameterReader.gemmi_isFirstInStage );
        if ( myParameterReader.gemmi_isFirstInStage ) { myParameterReader.gemmi_isFirstInStage = false; }
    }
    else
    {
        outputStream << "ENDMDL" << std::endl;
        
        outputStream <<"REMARK seconds since January 1st 1970: "<<time ( NULL     )<<std::endl; //<<"REMARK elapsed time: "<<(clock()/CLOCKS_PER_SEC)<<std::endl;
        outputStream <<"REMARK Current time is: "<<asctime (timeinfo) <<"REMARK elapsed time: "<<(clock()/CLOCKS_PER_SEC)<<std::endl;
        outputStream.setf(ios::fixed, ios::floatfield); // set output to fixed rather than scientific format
        if (myParameterReader.calcEnergy) {
            double myPotentialEnergy = system.calcPotentialEnergy(state);
             //DuMMForceFieldSubsystem& tempDumm = dumm;
            Real dummPotentialEnergy =     dumm.calcPotentialEnergy(state);
            outputStream <<"REMARK Total Potential Energy = "<<myPotentialEnergy <<" kJ/mol, "<<myPotentialEnergy/4.184<<" kcal/mol"<<std::endl;
            outputStream <<"REMARK DuMM Molecular Dynamics Potential Energy = "<<dummPotentialEnergy <<" kJ/mol, "<<dummPotentialEnergy/4.184<<" kcal/mol"<<std::endl;
            MMBLOG_PLAIN_NOSEV(ALWAYS, "Total Potential Energy =                   " << myPotentialEnergy << " kJ/mol, " << myPotentialEnergy/4.184 << " kcal/mol" << std::endl);
            MMBLOG_PLAIN_NOSEV(ALWAYS, "DuMM Molecular Dynamics Potential Energy = " << dummPotentialEnergy << " kJ/mol, " << dummPotentialEnergy/4.184<<" kcal/mol"<<std::endl);

            double myKineticEnergy = system.calcKineticEnergy(state);
            myParameterReader.kineticEnergy = myKineticEnergy;
            outputStream <<"REMARK Kinetic Energy = "<< myKineticEnergy <<" kJ/mol, "<<myKineticEnergy/4.184<<" kcal/mol"<<std::endl;
            MMBLOG_PLAIN_NOSEV(ALWAYS, "Kinetic Energy = " << myKineticEnergy << " kJ/mol, " << myKineticEnergy/4.184 << " kcal/mol" << std::endl);

            double myEnergy = system.calcEnergy(state);
            myParameterReader.totalEnergy = myEnergy;
            outputStream <<"REMARK Energy = "<< myEnergy <<" kJ/mol "<<std::endl;
            MMBLOG_PLAIN_NOSEV(ALWAYS, "REMARK Energy = " << myEnergy << " kJ/mol " << std::endl);
            myEnergies.push_back(myEnergy);
        }
        outputStream <<"REMARK Angular, Linear Momentum = "<<system.calcSystemRigidBodyMomentum(state)<<endl;

        //cout<<__FILE__<<" : "<<__LINE__<<" "<<dumm.getVdwMixingRuleName (dumm.getVdwMixingRule())<<endl;

        outputStream <<"REMARK ["<< __FILE__<<"] state.getNU()    = "<<state.getNU()            <<std::endl;
        //outputStream <<"REMARK ["<< __FILE__<<"]Satisfied contacts : "<<myParameterReader.satisfiedBasePairs<<endl;
        //outputStream <<"REMARK ["<< __FILE__<<"]Unsatisfied contacts : "<<myParameterReader.unSatisfiedBasePairs<<endl;
    }

    //cout<<__FILE__<<" : "<<__LINE__<<" "<<dumm.getVdwMixingRuleName (dumm.getVdwMixingRule())<<endl;

//    outputStream <<"REMARK ["<< __FILE__<<"] state.getNU()    = "<<state.getNU()            <<std::endl;
    //outputStream <<"REMARK ["<< __FILE__<<"]Satisfied contacts : "<<myParameterReader.satisfiedBasePairs<<endl;
    //outputStream <<"REMARK ["<< __FILE__<<"]Unsatisfied contacts : "<<myParameterReader.unSatisfiedBasePairs<<endl;


    GlobalProgressWriter::get().update(ProgressWriter::State::RUNNING, modelNumber);
    MMBLOG_PLAIN_NOSEV(ALWAYS, "Just wrote structure for reporting interval # " << modelNumber << std::endl);
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
                MMBLOG_PLAIN_NOSEV(ALWAYS, "Converged! Energy difference between two reporting intervals has been less than "<< myParameterReader.convergenceEpsilon << " kJ/mol for the last " << myParameterReader.convergenceTimeout << " frames." << endl);
            }
        }
    }

    MMBLogger::instance().flush();
}
