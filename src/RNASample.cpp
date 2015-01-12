#include "SimTKmolmodel.h"
#include "SimTKsimbody_aux.h"

//#include "/Users/samuelflores/rna-dynamics/Ligands.h"  
//#include "/Users/samuelflores/rna-dynamics/Water.h"    
//#include "/Users/samuelflores/rna-dynamics/WaterDroplet.h"
#include <list>
#include "Ligands.h"  
#include "RNATools/include/PeriodicPdbAndCOMWriter.h"  
#include "Water.h"    
//#include "molmodel/internal/GeometricCenter.h"
//#include "/Users/samuelflores/svn/molmodel/internal/molmodel/GeometricCenter.h"    
//#include "/home/scflores/rna-scratch/Water.h"    
#include "WaterDroplet.h"
//#include "/home/scflores/rna-scratch/WaterDroplet.h"
//#include "/Users/samuelflores/svn/molmodel/include/molmodel/internal/Ligands.h"
#include <iostream>
#include <vector>

using std::cout;
using std::endl;

using namespace SimTK;
using namespace std;
int main(int argc, char *argv[]) {

    // user configurable defaults
    int firstres=17;
    char* inFileName="1ARJ.simple.pdb";
    char * outTrajectoryName = "out.movie.pdb";
    char* moleculeChainId = "B";
    char*  sequence = "GGCAGAUCUGAGCCUGGGAGCUCUCUGCC";
    char* myBondMobilityType;
    float myCoulombGlobalScaleFactor = 0.0;
    float myGbsaGlobalScaleFactor = 0.0;
    float myVdwGlobalScaleFactor = 1.0;
    float myReportingInterval    =.1;
    float maxTime                =10000;//ps
    float thermostatTimeConstant = 0.0;
    int myUseMultithreadedComputation =0;
    float myTemperature =10000.;
    float waterRadius   = 1.50;
    int numWater         =0;//133 ;
    int numDivalents     = 0  ;
    int makeWaterDroplet=0;
    // end user configurable defaults

 for  (int q =1; q<argc;q+=2)
    {   
        cout<<"argv["<<q<<"]  = "<<argv[q]<<endl;
        string key = argv[q];
        if      (key=="-T") {myTemperature = atof(argv[q+1]);}
        else if (key     == "-I")   {myReportingInterval= atof(argv[q+1]);}
        else if (key     == "-S")   {inFileName = argv[q+1];}
        else if (key     == "-TFN") {outTrajectoryName = argv[q+1];}
        else if (key     == "-V")   {myVdwGlobalScaleFactor        = atof(argv[q+1]);}
        else if (key     == "-G")   {myGbsaGlobalScaleFactor       = atof(argv[q+1]);}
        else if (key     == "-C")   {myCoulombGlobalScaleFactor    = atof(argv[q+1]);}
        else if (key     == "-M")   {maxTime                       = atof(argv[q+1]);}
        else if (key     == "-TC")  {thermostatTimeConstant       = atof(argv[q+1]);}
        else if (key     == "-MT")  {myUseMultithreadedComputation= atoi(argv[q+1]);}
        else if (key     == "-ND")  {numDivalents                 = atoi(argv[q+1]);}
        else if (key     == "-NW")  {numWater                     = atoi(argv[q+1]);}
        else if (key     == "-WD")  {makeWaterDroplet             = atoi(argv[q+1]);}
        else if (key     == "-WR")  {waterRadius                  = atof(argv[q+1]);}
        else if (key     == "-BM")  {myBondMobilityType           = (argv[q+1]);}
        else if (key     == "-SEQ") {sequence           =(argv[q+1]);}
        else if (key     == "-CID") {moleculeChainId           =(argv[q+1]);}
    } 
    

    CompoundSystem system;
    SimbodyMatterSubsystem  matter(system);
    TinkerDuMMForceFieldSubsystem dumm(system);
    GeneralForceSubsystem forces; //this is needed only for the VanderWallSphere in WaterDroplet

    dumm.setCoulombGlobalScaleFactor(myCoulombGlobalScaleFactor);
    dumm.setGbsaGlobalScaleFactor(myGbsaGlobalScaleFactor);
    dumm.setVdwGlobalScaleFactor(myVdwGlobalScaleFactor);
    dumm.setUseMultithreadedComputation(myUseMultithreadedComputation);


    ifstream tinkerStream("./tinker_amber99_sam.prm");
    //ifstream tinkerStream("/Users/samuelflores/svn/molmodel/resources/tinker_amber99_sam.prm");
    dumm.populateFromTinkerParameterFile(tinkerStream);
    tinkerStream.close();
    RNA myMolecule(sequence,0);//
    MagnesiumIon myGCHelix2;

    MagnesiumIon myMagnesiumIonVec[numDivalents]; 
    for (int i = 0; i < numDivalents; i++) {
        (myMagnesiumIonVec[i]).setAmberLikeParameters(dumm);
        myMagnesiumIonVec[i].setPdbResidueNumber(i);
        myMagnesiumIonVec[i].setPdbChainId('C');    
        myMagnesiumIonVec[i].setPdbResidueName("MG2");
    
    } 

    myMolecule.assignBiotypes();
    
        (myGCHelix2).setAmberLikeParameters(dumm);
        myGCHelix2.setPdbResidueNumber(1);
        myGCHelix2.setPdbChainId('C');
        myGCHelix2.setPdbResidueName("MG2");


    myMolecule.setPdbChainId(*moleculeChainId);
    for (int q = 0; q <  (myMolecule.getNResidues()); q++)
         {
          //myMolecule.updResidue(Compound::Index(q)).setPdbChainId('B');
          cout<<"myMolecule.updResidue(Compound::Index(q)).setPdbResidueNumber(firstres+q);"<<std::endl;
          myMolecule.updResidue(Compound::Index(q)).setPdbResidueNumber(firstres+q);
        }

    std::ifstream inFileStream(inFileName,ifstream::in);
    PdbStructure pdbStructure(inFileStream);
    Compound::AtomTargetLocations atomTargets = myMolecule.createAtomTargets(pdbStructure);
    std::cout<<"atomtargest.szie "<<atomTargets.size()<<"versus atoms in myMolecule = "<<myMolecule.getNAtoms()<<std::endl;


    // Four steps to a perfect match
    myMolecule.matchDefaultAtomChirality(atomTargets);
    myMolecule.matchDefaultBondLengths(atomTargets);
    myMolecule.matchDefaultBondAngles(atomTargets);
    myMolecule.matchDefaultDihedralAngles(atomTargets);
    myMolecule.matchDefaultTopLevelTransform(atomTargets);
    Real residual = myMolecule.getTransformAndResidual(atomTargets).residual;
    //system.adoptCompound(myGCHelix2);
    for (int i = 1; i<numDivalents; i++)
    {
        system.adoptCompound(myMagnesiumIonVec[i],Vec3(i,0,0));
    }
    system.adoptCompound(myMolecule);
    if (makeWaterDroplet) WaterDroplet myReturnInt(system,dumm,forces);

    myMolecule.setCompoundBondMobility(BondMobility::Free);
    myMolecule.setRNABondMobility(BondMobility::Rigid,(17-firstres),(22-firstres));
    myMolecule.setRNABondMobility(BondMobility::Rigid,(26-firstres),(39-firstres));
    myMolecule.setRNABondMobility(BondMobility::Rigid,(40-firstres),(45-firstres));
     //39-17=22, 40-17=23
    myMolecule.setBondMobility(BondMobility::Free,"23/P","22/O3*");
    myMolecule.setBondMobility(BondMobility::Free,"22/O3*","22/C3*");
    myMolecule.setBondMobility(BondMobility::Free,"23/P"  ,"23/O5*");
    myMolecule.setBondMobility(BondMobility::Free,"23/O5*","23/C5*");
    

    system.modelCompounds();


    State & state = system.updDefaultState();
    system.realizeTopology();
    system.realize(state,Stage::Position);
    
    const SimTK::Transform& mytransform2 = matter.getMobilizedBody(myMolecule.getAtomMobilizedBodyIndex(Compound::AtomIndex(0))).findBodyTransformInAnotherBody(state,matter.getMobilizedBody(myMolecule.getAtomMobilizedBodyIndex(Compound::AtomIndex(myMolecule.getNAtoms()-1))));
    Constraint::Weld myWeld(matter.updMobilizedBody(
        myMolecule.getAtomMobilizedBodyIndex(Compound::AtomIndex(0))),Transform(Vec3(0)),
        matter.updMobilizedBody( myMolecule.getAtomMobilizedBodyIndex(Compound::AtomIndex(myMolecule.getNAtoms()-1))),mytransform2
        );
    system.realizeTopology();//realizeTopolgy is special.  invoke explicitly.
    system.realize(state,Stage::Position);
    Constraint::Weld myWeld2(matter.updMobilizedBody(
            myMolecule.getAtomMobilizedBodyIndex(Compound::AtomIndex(0))),Transform(Vec3(0)),
            matter.Ground(),
            matter.getMobilizedBody(myMolecule.getAtomMobilizedBodyIndex(Compound::AtomIndex(0))).findBodyTransformInAnotherBody(state,matter.Ground())
            //~Transform((matter.updMobilizedBody(myMolecule.getAtomMobilizedBodyIndex(Compound::AtomIndex(0)))).getBodyTransform(state)*myMolecule.getAtomLocationInMobilizedBodyFrame(Compound::AtomIndex(0)))
            );


    system.realizeTopology();//realizeTopolgy is special.  invoke explicitly.
    system.realize(state,Stage::Position);

   list<int> residueList ;
   residueList.push_back  ( 26  -  firstres  );
   residueList.push_back  ( 27  -  firstres  );
   residueList.push_back  ( 28  -  firstres  );
   residueList.push_back  ( 29  -  firstres  );
   residueList.push_back  ( 36  -  firstres  );
   residueList.push_back  ( 37  -  firstres  );
   residueList.push_back  ( 38  -  firstres  );
   residueList.push_back  ( 39  -  firstres  );

   cout<<"residueList begin, end:" << (residueList.front()) << ","<< (residueList.back())<<endl;
 
    //Vec3 helix2GC = GeometricCenter ( myMolecule , residueList, matter,  state);
    /*
*/
    //(matter.updMobilizedBody(((myGCHelix2.getAtomMobilizedBodyIndex(Compound::AtomIndex(0)))))).setQFromVector(state,Vector((0,0,0,0)));///helix2GC));
    //system.realizeTopology();//realizeTopolgy is special.  invoke explicitly.
    system.realize(state,Stage::Position);
    //const SimTK::Transform& mytransform3 = matter.getMobilizedBody((myMagnesiumIonVec[0]).getAtomMobilizedBodyIndex(Compound::AtomIndex(0))).findBodyTransformInAnotherBody(state,matter.getMobilizedBody((myMolecule.updResidue(Compound::Index(26 -firstres))).getAtomMobilizedBodyIndex(Compound::AtomIndex(0))));
/*
    Constraint::Weld myGCWeld(
        matter.updMobilizedBody( (myMolecule.updResidue(Compound::Index(26 -firstres))).getAtomMobilizedBodyIndex(Compound::AtomIndex(0))),
        ~((matter.updMobilizedBody( (myMolecule.updResidue(Compound::Index(26 -firstres))).getAtomMobilizedBodyIndex(Compound::AtomIndex(0)))).getBodyTransform(state)) *  Transform(helix2GC),
        matter.updMobilizedBody(((myGCHelix2)).getAtomMobilizedBodyIndex(Compound::AtomIndex(0))),
        Transform(Vec3(0))
        );
*/
    system.realizeTopology(); // realizeTopolgy is special.  invoke explicitly.
    system.realize(state,Stage::Position);


    VelocityRescalingThermostat * myVelocityRescalingThermostat = new VelocityRescalingThermostat(system,  myTemperature, myReportingInterval);
    system.updDefaultSubsystem().addEventHandler(myVelocityRescalingThermostat);
    cout<<"check 7.0"<<endl;
    RungeKuttaMersonIntegrator study(system);
    cout<<"check 8.0"<<endl;
    //system.realize(state,Stage::Position);
    study.initialize(state);
    cout<<"check 9.0"<<endl;
/*
    time_t rawtime;
    struct tm * timeinfo;
*/
        stringstream ss1;
        ss1<<"out.movie.pdb";//<<p<<".pdb";
        ofstream outputFrame(outTrajectoryName);
        //ofstream outputFrame(ss1.str().c_str());
              /*  list<int> residueList ;
                int firstres = 17; 
                residueList.push_back ( 26 - firstres );  
                residueList.push_back ( 27 - firstres );  
                residueList.push_back ( 28 - firstres );  
                residueList.push_back ( 29 - firstres );  
                residueList.push_back ( 36 - firstres );  
                residueList.push_back ( 37 - firstres );  
                residueList.push_back ( 38 - firstres );  
                residueList.push_back ( 39 - firstres );  
*/
        PeriodicPdbAndCOMWriter * myPeriodicPdbWriter = new PeriodicPdbAndCOMWriter(system,outputFrame,myReportingInterval,matter,myMolecule,dumm,residueList);
        system.updDefaultSubsystem().addEventReporter(myPeriodicPdbWriter);
        TimeStepper ts(system,study);
        ts.initialize(state);
        ts.stepTo(maxTime);//10000 *myReportingInterval);
}
