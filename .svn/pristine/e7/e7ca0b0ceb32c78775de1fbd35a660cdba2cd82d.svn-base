#include "SimTKmolmodel.h"
#include "SimTKsimbody_aux.h"
#include "CoarseNucleicAcid.h"
#include <fstream>
#include <ios>
#include <iostream>
#include <vector>
#include <stdlib.h>
using namespace SimTK;
using namespace std;

int main () {
        CompoundSystem system;
        SimbodyMatterSubsystem  matter(system);
        GeneralForceSubsystem forces(system);

        DuMMForceFieldSubsystem dumm(system);
	

        //dumm.loadAmber99Parameters();
        filebuf fb;
        fb.open (("/Users/Sam/svn/RNAToolbox/include/resources/tinker_amber99_sam.prm"),ios::in);
        istream is(&fb);
        //if (myParameterReader.loadTinkerParameterFile) 
        {
            //cout<<"You have specified tinkerParameterFileName = "<<myParameterReader.tinkerParameterFileName<<".  Checking this file.."<<endl;
            SimTK_ERRCHK_ALWAYS(fb.is_open(),"[Repel.h]", "The Tinker parameter file you specified could not be opened.  Please check your tinkerParameterFileName parameter, or set \"loadTinkerParameterFile 0\" to use the hard-coded Tinker parameters instead.");
            dumm.populateFromTinkerParameterFile (is);}

        dumm.setCoulombGlobalScaleFactor(0);//meterReader.globalCoulombScaleFactor);
        dumm.setBondTorsionGlobalScaleFactor(1);//meterReader.globalBondTorsionScaleFactor);
        dumm.setGbsaGlobalScaleFactor(0);//meterReader.globalGbsaScaleFactor);
        dumm.setVdwGlobalScaleFactor(0);//meterReader.globalVdwScaleFactor);
        dumm.setBondStretchGlobalScaleFactor(0);//meterReader.globalBondStretchScaleFactor);
        dumm.setBondBendGlobalScaleFactor(0);//meterReader.globalBondBendScaleFactor);
        dumm.setAmberImproperTorsionGlobalScaleFactor(0);//meterReader.globalAmberImproperTorsionScaleFactor);
         
        Biopolymer myMolecule;
        //myMolecule = Protein(    "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"); //50 //ok
        myMolecule = CoarseNucleicAcid(    "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"); //50 //ok

        myMolecule.assignBiotypes(); 
        system.adoptCompound(myMolecule,Vec3(  0,0,0));
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
        state = system.realizeTopology();
        system.realize(state,Stage::Position);
        VerletIntegrator study(system,.002);
        
        TimeStepper ts(system,study);
        ts.initialize(state);
        //# putting writeDefaultPdb here was fine to 2200 residues
        cout <<"1 about to write default pdb coords to  test.pdb"<<endl;
        myMolecule.writeDefaultPdb("test.pdb",Vec3(0));
        cout <<"1 just wrote default pdb coords to  test.pdb"<<endl;
        cout<<"[Repel.h:ConstrainedDynamics] Starting dynamics now."<<endl;
        time ( &rawtime );
        timeinfo = localtime ( &rawtime );
        cout<<" time:"<< asctime (timeinfo)   <<endl;
        //cout <<"2 about to write default pdb coords to  test.pdb"<<endl;
        //myMolecule.writeDefaultPdb("test.pdb",Vec3(0));
        //cout <<"2 just wrote default pdb coords to  test.pdb"<<endl;
        ts.stepTo(.1);
 
        time ( &rawtime );
        timeinfo = localtime ( &rawtime );
        cout<<" time:"<< asctime (timeinfo)   <<endl;


       cout<<"done."<<endl;
};   

