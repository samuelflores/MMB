#include "SimTKmolmodel.h"
#include "SimTKsimbody_aux.h"

#include <fstream>
#include <ios>
#include <iostream>
#include <vector>
#include <stdlib.h>
//#include "SPRNA.h"
using namespace SimTK;
using namespace std;

int main () {
        CompoundSystem system;
        SimbodyMatterSubsystem  matter(system);
        GeneralForceSubsystem forces(system);

        DuMMForceFieldSubsystem dumm(system);
        filebuf fb;
        fb.open ("/Users/Sam/svn/RNAToolbox/include/resources/tinker_amber99_sam.prm",ios::in);
        istream is(&fb);

        dumm.populateFromTinkerParameterFile(is );
        //dumm.loadAmber99Parameters();
        dumm.setCoulombGlobalScaleFactor(0);//meterReader.globalCoulombScaleFactor);
        dumm.setBondTorsionGlobalScaleFactor(0);//meterReader.globalBondTorsionScaleFactor);
        dumm.setGbsaGlobalScaleFactor(0);//meterReader.globalGbsaScaleFactor);
        dumm.setVdwGlobalScaleFactor(0);//meterReader.globalVdwScaleFactor);
        dumm.setBondStretchGlobalScaleFactor(0);//meterReader.globalBondStretchScaleFactor);
        dumm.setBondBendGlobalScaleFactor(0);//meterReader.globalBondBendScaleFactor);
        dumm.setAmberImproperTorsionGlobalScaleFactor(0);//meterReader.globalAmberImproperTorsionScaleFactor);
        dumm.setAllGlobalScaleFactors(0);        
        Biopolymer myMolecule;
        Biopolymer myMolecule2;
        myMolecule = RNA(    "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"); //100 //ok
        myMolecule2 = RNA(    "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"); //100 //ok
        //myMolecule = RNA(    "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"); //200 //ok
        //myMolecule = RNA(    "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"); //250 //ok
        //myMolecule = RNA(    "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"); //300  //fail   
        //myMolecule2 = RNA(    "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"); //300  //fail   
        //myMolecule = RNA(    "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"); //400 //fail


        system.adoptCompound(myMolecule,Vec3(  0,0,0));
        system.adoptCompound(myMolecule2,Vec3(  0,0,0));
        cout<<"about to modelCompounds"<<endl;
        time_t rawtime;
        struct tm * timeinfo;
        time ( &rawtime );
        timeinfo = localtime ( &rawtime );
        cout<<" time:"<< asctime (timeinfo)   <<endl;
        system.modelCompounds();
        time ( &rawtime );
        timeinfo = localtime ( &rawtime );
        cout<<" time:"<< asctime (timeinfo)   <<endl;
        cout<<"done with modelCompounds"<<endl;
        State &  state = system.updDefaultState();
        myMolecule.setCompoundBondMobility(BondMobility::Rigid);
        myMolecule2.setCompoundBondMobility(BondMobility::Rigid);
        state = system.realizeTopology();
        system.realize(state,Stage::Position);
        VerletIntegrator study(system,.002);
        study.setFixedStepSize(.001);

        string myFileName5 = "check.5.pdb";
        ofstream  myOfstream5(  myFileName5.c_str());


        if (0) Force::TwoPointLinearSpring mySpring(forces,  
            matter.updMobilizedBody(myMolecule.getAtomMobilizedBodyIndex(Compound::AtomIndex(myMolecule.getAtomIndex(string("0/N9") )))),
            Vec3(0,0,0),
            matter.updMobilizedBody(myMolecule2.getAtomMobilizedBodyIndex(Compound::AtomIndex(myMolecule2.getAtomIndex(string("0/N9") )))),
            Vec3(0),
            1,0);
        state = system.realizeTopology();
        system.realize(state,Stage::Position);
             
        //myMolecule.writePdb(state,myOfstream5,Vec3(0));
        
        TimeStepper ts(system,study);
        state.updU() = 1;
        ts.initialize(state);
        //# putting writeDefaultPdb here was fine to 2200 residues
        cout <<"1 about to write default pdb coords to  test.pdb"<<endl;
        //myMolecule.writeDefaultPdb("test.pdb",Vec3(0));
        //myMolecule2.writeDefaultPdb("test.pdb",Vec3(0));
        cout <<"1 just wrote default pdb coords to  test.pdb"<<endl;
        cout<<"[Repel.h:ConstrainedDynamics] Starting dynamics now."<<endl;
        time ( &rawtime );
        timeinfo = localtime ( &rawtime );
        cout<<" time:"<< asctime (timeinfo)   <<endl;

        //myMolecule.writePdb(state,myOfstream5,Vec3(0));
        //myMolecule2.writePdb(state,myOfstream5,Vec3(0));
        //cout <<"2 just wrote default pdb coords to  test.pdb"<<endl;
        ts.stepTo(.10);
        //system.realize(state,Stage::Dynamics);
        //myMolecule.writePdb(state,myOfstream5,Vec3(0));
        //myMolecule2.writePdb(state,myOfstream5,Vec3(0));
        state = ts.getState();
        system.realize(state,Stage::Dynamics);
        cout         <<"REMARK Angular, Linear Momentum = "<<system.calcSystemRigidBodyMomentum(state)<<endl;
        time ( &rawtime );
        timeinfo = localtime ( &rawtime );
        cout<<" time:"<< asctime (timeinfo)   <<endl;


       cout<<"done."<<endl;
};   

