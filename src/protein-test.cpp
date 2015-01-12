#include "SimTKmolmodel.h"
#include "SimTKsimbody_aux.h"
#include "ParameterReader.h"
//#include "PeriodicPdbAndEnergyWriter.h"
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
	

        dumm.loadAmber99Parameters();
        dumm.setCoulombGlobalScaleFactor(0);//meterReader.globalCoulombScaleFactor);
        dumm.setBondTorsionGlobalScaleFactor(1);//meterReader.globalBondTorsionScaleFactor);
        dumm.setGbsaGlobalScaleFactor(0);//meterReader.globalGbsaScaleFactor);
        dumm.setVdwGlobalScaleFactor(0);//meterReader.globalVdwScaleFactor);
        dumm.setBondStretchGlobalScaleFactor(1);//meterReader.globalBondStretchScaleFactor);
        dumm.setBondBendGlobalScaleFactor(0);//meterReader.globalBondBendScaleFactor);
        dumm.setAmberImproperTorsionGlobalScaleFactor(0);//meterReader.globalAmberImproperTorsionScaleFactor);
         
        vector <Biopolymer> myMolecule;
        //ParameterReader myParameterReader("commands.dat");
       
        //myMolecule = Protein(    "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"); //50 //ok
        Biopolymer tempMolecule = RNA    (    "UACGUAAGGA"); //50 //ok
        //Biopolymer tempMolecule = Protein ("FKVLQEPTCVSDYMSISTCEWKMNGPTNCSTELRLLYQLVFLLSEAHTCIPENNGGAGCVCHLLMDDVVSADNYTLDLWAGQQLLWKGSFKPSEHVKPRAPGNLTVHDTLLLTWSNPYPPDNYLYNHLTYAVNIWSENDPADFRIYNVTYLEPSLRI");
        tempMolecule.setPdbChainId(*("A"));
	ifstream inputFile("/Users/Sam/svn/RNAToolbox/trunk/src/in.pdb");

	PdbStructure pdbStructure(inputFile);
	inputFile.close();
        
        bool guessCoords = true ;
	Compound::AtomTargetLocations biopolymerAtomTargets = tempMolecule.createAtomTargets(pdbStructure,guessCoords);
	bool  matchHydrogenAtomLocations = false;
	if (! matchHydrogenAtomLocations)
	    {
	       map<Compound::AtomIndex, Vec3>::iterator it;
	       map<Compound::AtomIndex, Vec3>::iterator next;
	       next = biopolymerAtomTargets.begin();
	       while (next != biopolymerAtomTargets.end())
	       {
		   it = next;
		   Compound::AtomIndex m = (*it).first;
		   Element myAtomElement = tempMolecule.getAtomElement(m);
		   next++;
		   //cout<<__FILE__<<":"<<__LINE__<<(myAtomElement.getName())<<endl;
		   if  ((myAtomElement.getName()).compare("hydrogen") == 0)
			{
			biopolymerAtomTargets.erase(it);
			}
		   }
	    }


        cout<<__FILE__<<":"<<__LINE__<<" "<<biopolymerAtomTargets.size()<<endl;

        //tempMolecule.matchDefaultConfiguration(biopolymerAtomTargets,Compound::Match_Exact);

        guessCoords = false  ;
	biopolymerAtomTargets = tempMolecule.createAtomTargets(pdbStructure,guessCoords);

        tempMolecule.writeDefaultPdb("test.1.pdb",Vec3(0));
	tempMolecule.matchDefaultAtomChirality(biopolymerAtomTargets, 0.01, false);
        tempMolecule.writeDefaultPdb("test.2.pdb",Vec3(0));

	tempMolecule.matchDefaultBondLengths(biopolymerAtomTargets);
        tempMolecule.writeDefaultPdb("test.3.pdb",Vec3(0));
	tempMolecule.matchDefaultBondAngles(biopolymerAtomTargets);
        tempMolecule.writeDefaultPdb("test.3.pdb",Vec3(0));

	// Set dihedral angles even when bonded atoms are planar
	tempMolecule.matchDefaultDihedralAngles(biopolymerAtomTargets, Compound::DistortPlanarBonds);
        tempMolecule.writeDefaultPdb("test.4.pdb",Vec3(0));

	tempMolecule.matchDefaultTopLevelTransform(biopolymerAtomTargets);
        tempMolecule.writeDefaultPdb("test.5.pdb",Vec3(0));

        tempMolecule.fitDefaultConfiguration(biopolymerAtomTargets, 0.005);
        tempMolecule.writeDefaultPdb("test.6.pdb",Vec3(0));

/*
	tempMolecule.matchDefaultAtomChirality(biopolymerAtomTargets, 0.01, false);
	tempMolecule.matchDefaultBondLengths(biopolymerAtomTargets);
	tempMolecule.matchDefaultBondAngles(biopolymerAtomTargets);
	tempMolecule.matchDefaultDihedralAngles(biopolymerAtomTargets, Compound::DistortPlanarBonds);

	tempMolecule.matchDefaultTopLevelTransform(biopolymerAtomTargets);

        tempMolecule.fitDefaultConfiguration(biopolymerAtomTargets, 0.00001);
        tempMolecule.writeDefaultPdb("test.7.pdb",Vec3(0));
*/

        cout<<__FILE__<<":"<<__LINE__<<" "<<biopolymerAtomTargets.size()<<endl;
        //tempMolecule.matchDefaultConfiguration(biopolymerAtomTargets,Compound::Match_Idealized);
        //tempMolecule.matchDefaultTopLevelTransform(biopolymerAtomTargets);
        //tempMolecule.writeDefaultPdb("test.1.pdb",Vec3(0));

        tempMolecule.assignBiotypes(); 
        //myMolecule[0].assignBiotypes(); 
        
        system.adoptCompound(tempMolecule,Vec3(  0,0,0));
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

	ofstream outFile5("/Users/Sam/svn/RNAToolbox/trunk/src/test.9.pdb");
        tempMolecule.writePdb(state,outFile5 ,Vec3(0));
        //PeriodicPdbAndEnergyWriter * myPeriodicPdbWriter = new PeriodicPdbAndEnergyWriter(system,cout  ,  .1, myParameterReader,myMolecule );
        //system.updDefaultSubsystem().addEventReporter(myPeriodicPdbWriter);

        TimeStepper ts(system,study);
        ts.initialize(state);
        //# putting writeDefaultPdb here was fine to 2200 residues
        cout <<"1 about to write default pdb coords to  test.pdb"<<endl;
        tempMolecule.writeDefaultPdb("test.pdb",Vec3(0));
        cout <<"1 just wrote default pdb coords to  test.pdb"<<endl;
        //cout <<"REMARK Angular, Linear Momentum = "<<system.calcSystemRigidBodyMomentum(state)<<endl;
        cout<<"[Repel.h:ConstrainedDynamics] Starting dynamics now."<<endl;
        time ( &rawtime );
        timeinfo = localtime ( &rawtime );
        cout<<" time:"<< asctime (timeinfo)   <<endl;
        //cout <<"2 about to write default pdb coords to  test.pdb"<<endl;
        //myMolecule[0].writeDefaultPdb("test.pdb",Vec3(0));
        //cout <<"2 just wrote default pdb coords to  test.pdb"<<endl;
        ts.stepTo(.1);
        system.realize(state,Stage::Dynamics);

        cout <<"REMARK Angular, Linear Momentum = "<<system.calcSystemRigidBodyMomentum(state)<<endl;
        time ( &rawtime );
        timeinfo = localtime ( &rawtime );
        cout<<" time:"<< asctime (timeinfo)   <<endl;


       cout<<"done."<<endl;
};   

