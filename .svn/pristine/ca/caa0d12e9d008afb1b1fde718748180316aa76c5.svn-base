#include "SimTKmolmodel.h"
#include "SimTKsimbody_aux.h"

#include <fstream>
#include <ios>
#include <iostream>
#include <vector>

//#include "/Users/samuelflores/svn/RNAToolbox/include/Repel.h"
using namespace SimTK;
using namespace std;

int main(int argc, char *argv[] ) {


 char*  inputFileName = "/Users/samuelflores/svn/tar-dynamics/parameterization/in.pdb";	
 char*  outputFileName = "/Users/samuelflores/svn/tar-dynamics/parameterization/out.pdb";
 char firstChainID = ' ';	
 char secondChainID = ' ';
 char* firstResidueType = "G";
 char* secondResidueType = "G";
 int firstResidueNumber = 0;	
 int secondResidueNumber = 0;	
 BondMobility::Mobility overallBondMobility = BondMobility::Free;   
/* cout<<  "-FRN  firstResidueNumber"<<endl;
 cout<<  "-SRN  secondResidueNumber"<<endl;
 cout<<  "-FRT  firstResidueType"<<endl;
 cout<<  "-SRT  secondResidueType"<<endl;
 cout<<  "-FCID firstChainID"<<endl;
 cout<<  "-SCID secondChainID"<<endl;
 cout<<  "-FIN  inputFileNam"<<endl;
 cout<<  "-FOUT outputFileName"<<endl;
*/

 for  (int q =1; q<argc;q+=2)
    {   
        cout<<"argv["<<q<<"]  ="<<argv[q]<<endl;
        string key = argv[q];
        if      (key == "-FRN")  {firstResidueNumber  = atoi(argv[q+1]);}
        else if (key == "-SRN")  {secondResidueNumber = atoi(argv[q+1]);}
        else if (key == "-FRT")  {firstResidueType    = (argv[q+1]);}
        else if (key == "-SRT")  {secondResidueType   = (argv[q+1]);}
        else if (key == "-FCID") {firstChainID        = *(argv[q+1]);}
        else if (key == "-SCID") {secondChainID       = *(argv[q+1]);}
        else if (key == "-FIN")  {inputFileName       = (argv[q+1]);}
        else if (key == "-FOUT") {outputFileName      = (argv[q+1]);}
        //else if (key == "-BM"  ) {overallBondMobility = String(argv[q+1]);}
    } 
    if   (!(firstResidueNumber>0)) {cout<<"ERROR : invalid firstResidueNumber (-FRN) "<<endl;  }
    assert(firstResidueNumber>0);  
    if   (!(secondResidueNumber>0)) {cout<<"ERROR : invalid secondResidueNumber (-SRN) "<<endl; } 
    assert(secondResidueNumber>0);  
        cout<<endl<<"You have specified or defaulted to the following parameters:"<<endl<<endl;
        cout<<"-FRN  firstResidueNumber  ="<<firstResidueNumber<<endl;
        cout<<"-SRN  secondResidueNumber ="<<secondResidueNumber<<endl; 
        cout<<"-FRT  firstResidueType    ="<<firstResidueType  <<endl;
        cout<<"-SRT  secondResidueType   ="<<secondResidueType  <<endl; 
        cout<<"-FCID firstChainID        =>"<<firstChainID <<"<"<<endl;
        cout<<"-SCID secondChainID       =>"<< secondChainID<<"<"<<endl;
        cout<<"-FIN  inputFileName       ="<< inputFileName<<endl;
        //cout<<"-FOUT outputFileName      ="<< outputFileName<<endl;
        cout<<"-BM   overallBondMobility =" <<  overallBondMobility <<endl<<endl; 


    if   (!(firstResidueNumber>0)) {cout<<"ERROR : invalid firstResidueNumber (-FRN) "<<endl;  }
    assert(firstResidueNumber>0);  
    if   (!(secondResidueNumber>0)) {cout<<"ERROR : invalid secondResidueNumber (-SRN) "<<endl; } 
    assert(secondResidueNumber>0);  

    String atomPathNameA , atomPathNameB;
    //if (((string(firstResidueType)).compare("A") == 0 ) || ((string(firstResidueType)).compare("G") == 0 )) {atomPathNameA = "0/C1*";}
    if (((string(firstResidueType)).compare("A") == 0 ) || ((string(firstResidueType)).compare("G") == 0 )) {atomPathNameA = "0/N9";}
    else if (((string(firstResidueType)).compare("C") == 0 ) || ((string(firstResidueType)).compare("U") == 0 )) {atomPathNameA = "0/N1";}
    if (((string(secondResidueType)).compare("A") == 0 ) || ((string(secondResidueType)).compare("G") == 0 )) {atomPathNameB = "0/N9";}
    else if (((string(secondResidueType)).compare("C") == 0 ) || ((string(secondResidueType)).compare("U") == 0 )) {atomPathNameB = "0/N1";}
    RNA myMoleculeA(firstResidueType);
    //RNA myMoleculeA("G");
    RNA myMoleculeB(secondResidueType);
    myMoleculeA.setCompoundBondMobility(overallBondMobility);	
    myMoleculeB.setCompoundBondMobility(overallBondMobility);	
    //RNA myMoleculeB("G");
    (myMoleculeA.updResidue(ResidueInfo::Index(0))).setPdbResidueNumber(firstResidueNumber);
    //(myMoleculeA.updResidue(ResidueInfo::Index(1))).setPdbResidueNumber(2);
    (myMoleculeB.updResidue(ResidueInfo::Index(0))).setPdbResidueNumber(secondResidueNumber);
    myMoleculeA.setPdbChainId(firstChainID);    
    myMoleculeB.setPdbChainId(secondChainID);    
    //myMoleculeA.setPdbChainId('A');    
    //myMoleculeB.setPdbChainId('B');    
    //myMoleculeB.setPdbChainId('B');    
    //don't forget to remove the terminal H3 from the  input file:
    //string basePairFileName = "/Users/samuelflores/svn/tar-dynamics/parameterization/A-WC-U-WC.pdb";	
   
    // End user configurable parameters


    CompoundSystem system;
    GeneralForceSubsystem forces(system);
    DuMMForceFieldSubsystem dumm(system);
    HuntCrossleyContact myHuntCrossleyContact(system);

    dumm.loadAmber99Parameters();

    SimbodyMatterSubsystem  matter(system);

    //string basePairFileName = "/Users/samuelflores/svn/tar-dynamics/parameterization/A-WC-U-WC.pdb";	
    ifstream basePairFile (inputFileName);//basePairFileName.c_str());	 
    assert(basePairFile.good());
    PdbStructure pdbStructure(basePairFile);
    Compound::AtomTargetLocations atomTargetsA = myMoleculeA.createAtomTargets(pdbStructure); 
    Compound::AtomTargetLocations atomTargetsB = myMoleculeB.createAtomTargets(pdbStructure); 
    cout<<"size of atomTargets ="<<atomTargetsA.size()<<endl;    
    //cout<<"number of atoms in myMolecule ="<<myMoleculeA.getNAtoms()<<endl;

    //myMoleculeA.fitDefaultConfiguration(atomTargetsA,0.01);

    myMoleculeA.matchDefaultAtomChirality(atomTargetsA,.010);
    myMoleculeA.matchDefaultBondLengths(atomTargetsA);
    myMoleculeA.matchDefaultBondAngles(atomTargetsA);
    myMoleculeA.matchDefaultDihedralAngles(atomTargetsA);
    myMoleculeA.matchDefaultTopLevelTransform(atomTargetsA);
    cout<<"check 1"<<endl;
    //myMoleculeB.fitDefaultConfiguration(atomTargetsB,0.01);
    myMoleculeB.matchDefaultAtomChirality(atomTargetsB,.010);
    myMoleculeB.matchDefaultBondLengths(atomTargetsB);
    myMoleculeB.matchDefaultBondAngles(atomTargetsB);
    myMoleculeB.matchDefaultDihedralAngles(atomTargetsB);
    myMoleculeB.matchDefaultTopLevelTransform(atomTargetsB);
    cout<<"check 2"<<endl;
    basePairFile.close(); 
    cout<<"check 3"<<endl;
    system.adoptCompound(myMoleculeA);
    cout<<"check 4"<<endl;
    system.adoptCompound(myMoleculeB);
    system.modelCompounds();
    cout<<"check 5"<<endl;
    State state = system.realizeTopology();
    cout<<"check 6"<<endl;
    system.realize(state,Stage::Position);
    cout<<"check 7"<<endl;
    /*ofstream outputFileStream(outputFileName);
    myMoleculeA.writePdb(state,outputFileName,Transform(Vec3(0)));
    myMoleculeB.writePdb(state,outputFileName,Transform(Vec3(0)));
    outputFileStream.close();*/
    myMoleculeA.writeDefaultPdb(("firstResidueOut.pdb"),Transform(Vec3(0)));
    myMoleculeB.writeDefaultPdb(("secondResidueOut.pdb"),Transform(Vec3(0)));
    cout<<"check 7.5"<<endl;
    //cout<<"[generate-base-pair-transform.cpp]   :"<< myMoleculeA.getAtomIndex("0/N1")<<","<<myMoleculeA.getAtomIndex("0/N9")  <<endl;
    cout<<"check 7.6"<<endl;
    //cout<<"[generate-base-pair-transform.cpp]   :"<<myMoleculeB.getAtomMobilizedBodyIndex( myMoleculeB.getAtomIndex("0/N1"))<<endl;
    cout<<","<<myMoleculeA.getAtomMobilizedBodyIndex(Compound::AtomIndex(0));//eculeA.getAtomIndex("0/N9"))  <<endl;
    //cout<<","<<myMoleculeA.getAtomMobilizedBodyIndex(myMoleculeA.getAtomIndex("0/N9"))  <<endl;
    cout<<"check 7.8"<<endl;
    cout <<"[generate-base-pair-transform.cpp] chain A location of "<<atomPathNameA<<" = "<<myMoleculeA.calcAtomLocationInGroundFrame(state,myMoleculeA.getAtomIndex(atomPathNameA))<<endl;
    cout << "[generate-base-pair-transform.cpp] chain B location of "<<atomPathNameB<<" = "<<myMoleculeB.calcAtomLocationInGroundFrame(state,myMoleculeB.getAtomIndex(atomPathNameB))<<endl;
    Transform bodyTransform2InFrameOf1 =matter.getMobilizedBody(myMoleculeB.getAtomMobilizedBodyIndex(myMoleculeB.getAtomIndex(atomPathNameB))).findBodyTransformInAnotherBody(state,matter.getMobilizedBody(myMoleculeA.getAtomMobilizedBodyIndex(myMoleculeA.getAtomIndex(atomPathNameA))));
    cout<<"check 8"<<endl;
    cout<<"[generate-base-pair-transform.cpp] bodyTransform2InFrameOf1.T() :"<< bodyTransform2InFrameOf1.T()<<endl;
    //cout<<"[generate-base-pair-transform.cpp] bodyTransform2InFrameOf1.T().norm :"<< (bodyTransform2InFrameOf1.T()).norm()<<endl;
    cout<<"[generate-base-pair-transform.cpp] bodyTransform2InFrameOf1.R() :"<< (bodyTransform2InFrameOf1.R()).convertRotationToAngleAxis()<<endl;
}

