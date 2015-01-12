#ifndef PeriodicScrubber_H_
#define PeriodicScrubber_H_

/* -------------------------------------------------------------------------- *
 *                           MMB (MacroMoleculeBuilder)                       *
 * -------------------------------------------------------------------------- *
 *                                                                            *
 * Copyright (c) 2011-12 by the Author.                                       *
 * Author: Samuel Flores                                                      *
 *                                                                            *
 * See RNABuilder.cpp for the copyright and usage agreement.                  *
 * -------------------------------------------------------------------------- */

#include "SimTKsimbody.h"
#include <stdio.h>
#include "ParameterReader.h"
//#include 
using namespace std;
namespace SimTK {

/// Writes atomic coordinates in PDB format to a file stream at
/// specified intervals during a simulation.

class PeriodicScrubber : public PeriodicEventHandler  {
public:
    PeriodicScrubber(
        const CompoundSystem& system, 
        GeneralForceSubsystem& forces,
        DuMMForceFieldSubsystem& dumm,
        ParameterReader& myParameterReader,
        Real halfPeriod//,
        /*GeneralContactSubsystem & contacts,
        ContactSetIndex contactSet,
        HuntCrossleyForce & hc*/
	) 
        : PeriodicEventHandler (halfPeriod/10), 
          system(system), 
	  forces(forces),
          dumm(dumm),
          myParameterReader(myParameterReader),
	  halfPeriod(halfPeriod)//,
          /*contacts(contacts),
	  contactSet(contactSet),
          hc(hc)*/
    {}
    void handleEvent( State& state, Real accuracy, bool& shouldTerminate) const {
        
        double dutyCycle  = myParameterReader.dutyCycle;
 
        if (((fmod(state.getTime(),  2*halfPeriod)/(2*halfPeriod)) -  (1-dutyCycle))    > -.0000001)  {
            std::cout<<__FILE__<<":"<<__LINE__<<": setting setForceIsDisabled to 0"<<std::endl;
            for (int i = 0; i < forces.getNumForces(); i++) { 
                forces.setForceIsDisabled(state,SimTK::ForceIndex(i),0);
	    }
	    dumm.setCoulombGlobalScaleFactor(myParameterReader.globalCoulombScaleFactor);
	    dumm.setBondTorsionGlobalScaleFactor(myParameterReader.globalBondTorsionScaleFactor);
	    dumm.setGbsaGlobalScaleFactor(myParameterReader.globalGbsaScaleFactor);
	    dumm.setVdwGlobalScaleFactor(myParameterReader.globalVdwScaleFactor);
	    dumm.setBondStretchGlobalScaleFactor(myParameterReader.globalBondStretchScaleFactor);
	    dumm.setBondBendGlobalScaleFactor(myParameterReader.globalBondBendScaleFactor);
	    dumm.setAmberImproperTorsionGlobalScaleFactor(myParameterReader.globalAmberImproperTorsionScaleFactor);
	    dumm.setCustomBondStretchGlobalScaleFactor(0);
	    dumm.setCustomBondBendGlobalScaleFactor(0);
            //myParameterReader.myBiopolymerClassContainer.setContactParameters(contacts,  hc,  myParameterReader.excludedVolumeStiffness, false /* set contact forces to zero */      ); 
        }
	else {
            std::cout<<__FILE__<<":"<<__LINE__<<": setting setForceIsDisabled to 1 for "<< forces.getNumForces()<<" forces."<<  std::endl;
            for (int i = 0; i <   forces.getNumForces(); i++) { 
                forces.setForceIsDisabled(state,SimTK::ForceIndex(i),1);
	    }
            dumm.setAllGlobalScaleFactors(0);
            //myParameterReader.myBiopolymerClassContainer.setContactParameters(contacts,  hc,  myParameterReader.excludedVolumeStiffness, true  /* set contact forces to default strength */      ); 
        }
        // we need to realize to Topology stage, since we just modified the forces:
        // for some reason, this was giving bad results:
        //state = system.realizeTopology();
        // removed, hopefully topology will be realized when needed.

    }

private:
    const CompoundSystem& system;
    GeneralForceSubsystem& forces;
    DuMMForceFieldSubsystem& dumm;
    ParameterReader& myParameterReader;
    Real halfPeriod;
    /*GeneralContactSubsystem & contacts;
    ContactSetIndex contactSet;
    HuntCrossleyForce & hc;*/
};
} // namespace SimTK
#endif
