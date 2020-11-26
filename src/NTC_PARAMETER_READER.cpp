/* -------------------------------------------------------------------------- *
 *                           MMB (MacroMoleculeBuilder)                       *
 * -------------------------------------------------------------------------- *
 *                                                                            *
 * Copyright (c) 2011-12 by the Author.                                       *
 * Author: Samuel Flores                                                      *
 *                                                                            *
 * See RNABuilder.cpp for the copyright and usage agreement.                  *
 * -------------------------------------------------------------------------- */

//#include <ostream>
#include <iostream>
#include <istream>
#include <fstream>
#include "SimTKsimbody.h"
#include "NtC_Class_Container.h"
#include "NTC_FORCE_CLASS.h"
#include "NTC_PARAMETER_READER.h"

const int numNTC_PAR_BondMatrixRows=(39600);  // yes, global constants are bad.  Couldn't think of an elegant way around this one though.  This number should be exactly equal to the number of rows in the leontisWesthofBondMatrix
//const int maxParallelTorques = 1000; //max number of parallel torques to be applied.  This can be huge, minimal cost for doing that.

using namespace SimTK;
using namespace std;

/**
 * 
 * 
 * /param 
 * myPdbResidueName1,2 must be one of "A","C","G","U".
 * bondingEdge1,2 must be one of "WatsonCrick","Hoogsteen","Sugar","Bifurcated".
 * dihedraltype must be either "Cis" or "Trans".
 *
 */

    NTC_PAR_BondKey::NTC_PAR_BondKey(String myPdbResidueName1, String myPdbResidueName2,String myBondingEdge1, String myBondingEdge2, String dihedraltype, String myIsTwoTransformForce) :
        pdbResidueName1(std::move(myPdbResidueName1)),
        pdbResidueName2(std::move(myPdbResidueName2)),
        bondingEdge1(std::move(myBondingEdge1)),
        bondingEdge2(std::move(myBondingEdge2)),
        dihedraltype(std::move(dihedraltype)),
        isTwoTransformForce(std::move(myIsTwoTransformForce)) {
    }

    NTC_PAR_BondKey::NTC_PAR_BondKey(NTC_PAR_BondRow myNTC_PAR_BondRow) {
        pdbResidueName1 = myNTC_PAR_BondRow.pdbResidueName1;
        pdbResidueName2 = myNTC_PAR_BondRow.pdbResidueName2;
        bondingEdge1 = myNTC_PAR_BondRow.bondingEdge1;
        bondingEdge2 =myNTC_PAR_BondRow.bondingEdge2;
        dihedraltype = myNTC_PAR_BondRow.dihedraltype;
        isTwoTransformForce = myNTC_PAR_BondRow.isTwoTransformForce; 
    }; 	
//};

    int  NTC_PAR_Class::initialize( String inFileName) {
        NTC_PAR_Map.clear();
        myNTC_PAR_BondMatrix.myNTC_PAR_BondRow.clear();
        ifstream inFile(inFileName.c_str(),ifstream::in);
        cout<<"Now checking for existence of "<<inFileName<<endl;

        if (!(inFile.good())) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "Unable to open parameter file "<<inFileName<<endl);
        }
        int q=0;
	//char * s; 
        string s;
        //s = new char[500];
        while (inFile.good()) {
            std::getline(inFile,s,',');
            //inFile.getline(s,500,',');
             

            if ((String(s)).compare("NTCRECORD") == 0)  { //if this is a RECORD entry
                std::getline(inFile,s,',');
            	//inFile.getline( s, 100,',' );
                NTC_PAR_BondRow tempRow;
                myNTC_PAR_BondMatrix.myNTC_PAR_BondRow.push_back(tempRow);
	        (myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).pdbResidueName1 = String(s); // resname A, G, T, C 
                std::getline(inFile,s,',');
                //inFile.getline( s, 100,',' );
	        (myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).pdbResidueName2 = String(s); // resname A, G, T, C column 3
                std::getline(inFile,s,',');
                //inFile.getline( s, 100,',' );
	        (myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).bondingEdge1    = String(s); // type NTC Class C 4 
                std::getline(inFile,s,',');
                //inFile.getline( s, 100,',' );
	        (myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).bondingEdge2    = String(s); // type NTC Class C 5 
                std::getline(inFile,s,',');
                //inFile.getline( s, 100,',' );
	        (myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).dihedraltype = String(s); // type dihedral C 6
;
                for (int r=0; r<4; r++) {
                    std::getline(inFile,s,',');
	            //inFile.getline( s, 100,',' );
	            (myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).residue1Atom[r] = String(s); // res atom 1 C7-C10
                } 
                for (int r=0; r<4; r++) {
                    std::getline(inFile,s,',');
	            //inFile.getline( s, 100,',' );
	            (myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).atom_shift[r] = String(s); // atom shift C11-C14
                } 
                for (int r=0; r<4; r++) {
                    std::getline(inFile,s,',');
	            //inFile.getline( s, 100,',' );
	            (myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).bondLength[r] =   (double)atof(s.c_str()); // bond length C15-C18
                } 
                for (int r=0; r<4; r++) {
                    std::getline(inFile,s,',');
	            //inFile.getline( s, 100,',' );
	            (myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).springConstant[r] = (double)atof(s.c_str()); // 
                } 
                std::getline(inFile,s,',');
	        //inFile.getline(s,100,',');
	        (myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).torqueConstant =  (double)atof(s.c_str()); // torque Constant
		assert ((myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).torqueConstant >= 0);
                for (int r=0; r<3; r++) {
                    std::getline(inFile,s,',');
	            //inFile.getline( s, 100,',' );
	            (myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).attachmentPoint[r] = atof(s.c_str());
                }
                std::getline(inFile,s,',');
                //inFile.getline( s, 100,',' );
                (myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).rotationAngle = (double)atof(s.c_str());
                for (int r=0; r<2; r++) {
                    std::getline(inFile,s,',');
	            //inFile.getline( s, 100,',' );
	            (myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).rotationAxis[r] = atof(s.c_str());
                }
                
                std::getline(inFile,s,',');
                (myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).CONFALVALUE = atof(s.c_str());
                
                std::getline(inFile,s,',');
                //inFile.getline( s, 100,',' );
                (myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).isTwoTransformForce= String(s);

                std::getline(inFile,s,',');
                //inFile.getline( s, 100,',' );
                (myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).distanceC1pC1p = (double)atof(s.c_str());

                NTC_PAR_Map[NTC_PAR_BondKey(myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q])] = myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q];
                
              //  cout << (myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).pdbResidueName1 << (myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).pdbResidueName2 << (myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).dihedraltype << endl;
              //  cout << (myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).atom_shift[0] << (myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).atom_shift[1] << (myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).atom_shift[2] << (myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).atom_shift[3] << endl; 
              //  cout << (myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).isTwoTransformForce << endl;
                
              q++;
            }
            }
        SimTK_ERRCHK_ALWAYS(NTC_PAR_Map.size() ==  myNTC_PAR_BondMatrix.myNTC_PAR_BondRow.size() ,"[BaseInteractionParameterReader.cpp]","Inconsistency in number of NTC parameter rows. This probably means that your parameter file tried to specify parameters for the same interaction twice!"); 
         
      //  cout << myNTC_PAR_BondMatrix.myNTC_PAR_BondRow.size() << "my NTC PAR bond row size " << endl; 
        
        //delete[] s;	
        inFile.close();
        cout<<"done initializing myNTC_PAR_BondMatrix"<<endl;

        return(0);
        };


    void NTC_PAR_Class::printNTC_PAR_BondRows () {    
        for   (int q =0; q< (int)myNTC_PAR_BondMatrix.myNTC_PAR_BondRow.size(); q++) 
            cout<<"[NTCParameterReader.cpp] 269: "<<(myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).pdbResidueName1 
              <<(myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).pdbResidueName2
                <<(myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).bondingEdge1
                <<(myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).bondingEdge2
                <<(myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).dihedraltype << endl;
    //            <<(myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).isTwoTransformForce<<endl;
    };


    int  NTC_PAR_Class::getNTC_PAR_BondRowIndex(String myPdbResidueName1, String myPdbResidueName2,String Classtype,String dihedraltype,String myBasePairIsTwoTransformForce, NTC_Classes NTC) const {
        //if (0) { //!((myBasePairIsTwoTransformForce.compare("aromatic") == 0) || (myBasePairIsTwoTransformForce.compare("baseInteraction") == 0))) {
        //    return -11111;
        //} else  {
      
        cout << (int)myNTC_PAR_BondMatrix.myNTC_PAR_BondRow.size() << " br size " << " " << Classtype << " " << dihedraltype << " " << myPdbResidueName1 << " "<< myPdbResidueName2 << endl;
        
        for   (int q =0; q< (int)myNTC_PAR_BondMatrix.myNTC_PAR_BondRow.size(); q++) {
            
            
            if (
                ((((myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).pdbResidueName1).compare(myPdbResidueName1)) ==0)  &&
                ((((myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).pdbResidueName2).compare(myPdbResidueName2)) ==0)&&
                ((((myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).bondingEdge1).compare(Classtype))==0) &&
                ((((myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).bondingEdge2).compare(Classtype))==0) &&
                ((((myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).dihedraltype).compare(dihedraltype)) == 0) //&&
                //((((myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).isTwoTransformForce).compare(myBasePairIsTwoTransformForce) == (0) )  )
                )
                {
                    
                    (NTC.NtC_atom_type1)=((myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).residue1Atom[0]);
                    (NTC.NtC_atom_type2)=((myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).residue1Atom[1]);
                    (NTC.NtC_atom_type3)=((myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).residue1Atom[2]);
                    (NTC.NtC_atom_type4)=((myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).residue1Atom[3]);
                    (NTC.NtC_dihedraltype)=((dihedraltype));
                    (NTC.Harmonic_pot_constant)=((myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).springConstant[0]);
                    (NTC.Residue_shift_atom1)=((myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).atom_shift[0]);
                    (NTC.Residue_shift_atom2)=((myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).atom_shift[1]);
                    (NTC.Residue_shift_atom3)=((myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).atom_shift[2]);
                    (NTC.Residue_shift_atom4)=((myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).atom_shift[3]);                    
                    (NTC.Rotation_angle)=((myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).rotationAngle);
                     NTC.NTC_PAR_BondRowIndex = q;
                    
        //       cout <<  (myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).pdbResidueName1 <<  (myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).pdbResidueName2 <<  (myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).bondingEdge1 << (myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).bondingEdge2 << (myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).dihedraltype << endl;  
                    
		     if(0)  cout<<"[BaseInteractionParameterReader.cpp] found the right NTC_PAR_BondRow. residue1Atom[0], residue1Atom[1], residue1Atom[2] ,residue1Atom[3] residue2Atom[0], residue2Atom[1], residue2Atom[2] ,residue2Atom[3]  ="<<
                    (myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).residue1Atom[0]<<","<<
                    (myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).residue1Atom[1]<<","<<
                    (myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).residue1Atom[2]<<","<<
                    (myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).residue1Atom[3]<<","<<
                    (myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).atom_shift[0]<<","<<
                    (myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).atom_shift[1]<<","<<
                    (myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).atom_shift[2]<<","<<
                    (myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).atom_shift[3]<<
	            endl;
           //         if (myBasePairIsTwoTransformForce.compare("contact") != 0) SimTK_ERRCHK_ALWAYS( fabs(((myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).rotationAxis).norm()-1.0) <.001 ,"[BaseInteractionParameterReader.cpp]","The desired interaction was found but the norm of its rotationAxis is not unity within tolerance of .001.  The interaction may be blank or incorrect.  "); 
                    
                    return q ; //myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q];
                }	 	
                  
              if(0) cout<<"[BaseInteractionParameterReader.cpp] looking at NTC_PAR_BondRow. residue1Atom[0], residue1Atom[1], residue1Atom[2] ,residue1Atom[3] residue2Atom[0], residue2Atom[1], residue2Atom[2] ,residue2Atom[3], bondingEdge1, bondingEdge2, dihedraltype  ="<<
                    
                    ((((myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).pdbResidueName1)))<<","<<  
                    ((((myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).pdbResidueName2)))<<","<<
                    (myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).residue1Atom[0]<<","<<
                    (myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).residue1Atom[1]<<","<<
                    (myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).residue1Atom[2]<<","<<
                    (myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).residue1Atom[3]<<","<<
                    (myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).atom_shift[0]<<","<<
                    (myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).atom_shift[1]<<","<<
                    (myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).atom_shift[2]<<","<<
                    (myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).atom_shift[3]<<","<<

                    (myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).bondingEdge1   <<","<<
                    (myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).bondingEdge2   <<","<<
                    (myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).dihedraltype<<","<<
	            endl; 
                    if (0) cout<<"[BaseInteractionParameterReader.cpp] pdbResidueName1, pdbResidueName2, ="<<(myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).pdbResidueName1<<","<<(myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).pdbResidueName2<<endl;
                   
        }
	    /*cout<<"[BaseInteractionParameterReader.cpp] failed to match :"<<endl<<
               ","<<(myBasePairIsTwoTransformForce)<<
                ","<<(myPdbResidueName1)<<
                ","<<(myPdbResidueName2)<<
                ","<<(myBondingEdge1)<<
                ","<<(myBondingEdge2)<<
                ","<<(dihedraltype)<<endl;  */
                //","<<(dihedraltype)<<endl;  
            SimTK_ERRCHK_ALWAYS(0,"[BaseInteractionParameterReader.cpp]","Found no match for the above user-specified interaction.  Either add this interaction type to the parameter file, or check your spelling, syntax, or semantics.");

    //}
    };


    
    NTC_PAR_BondRow NTC_PAR_Class::getNTC_PAR_BondRow(ResidueID myResidueNumber1,ResidueID myResidueNumber2, String myPdbResidueName1, String myBondingEdge1, String myPdbResidueName2,String myBondingEdge2, String dihedraltype,String myBasePairIsTwoTransformForce) const {

        static map <const NTC_PAR_BondKey, NTC_PAR_BondRow, NTC_PAR_BondKeyCmp>::iterator iter = NTC_PAR_Map.begin();

        iter = NTC_PAR_Map.find(NTC_PAR_BondKey(myPdbResidueName1, myPdbResidueName2, myBondingEdge1,  myBondingEdge2,  dihedraltype,  myBasePairIsTwoTransformForce));
       
        NTC_PAR_BondRow myReturnNTC_PAR_BondRow;


        if 
            (!((myPdbResidueName1.compare(myReturnNTC_PAR_BondRow.pdbResidueName1) == 0) &&
            ( myPdbResidueName2.compare(myReturnNTC_PAR_BondRow.pdbResidueName2) == 0) &&            (myBondingEdge1.compare(myReturnNTC_PAR_BondRow.bondingEdge1) == 0)        &&
            (myBondingEdge2.compare(myReturnNTC_PAR_BondRow.bondingEdge2) == 0)        &&            (dihedraltype.compare(myReturnNTC_PAR_BondRow.dihedraltype) == 0) &&
            (myBasePairIsTwoTransformForce.compare(myReturnNTC_PAR_BondRow.isTwoTransformForce) == 0)))
                {

                cout<<"[BaseInteractionParameterReader.cpp] for interaction between residues "<<myResidueNumber1.getResidueNumber()<< " and "<<myResidueNumber2.getResidueNumber() <<endl<<"trying to match :"<<endl<<
                    ","<<(myBasePairIsTwoTransformForce)<<  
                    ", myPdbResidueName1"<<(myPdbResidueName1)<<
                    ", myPdbResidueName2"<<(myPdbResidueName2)<<
                    ", myBondingEdge1"<<(myBondingEdge1)<<
                    ", myBondingEdge2"<<(myBondingEdge2)<<
                    ", dihedraltype"<<(dihedraltype)<<endl;

                } 

        return myReturnNTC_PAR_BondRow;
        if (0) cout<<"Inside getNTC_PAR_BondRow.  about to search for :"<<myPdbResidueName1<<","<<myBondingEdge1<<","<< myPdbResidueName2   <<","<<  myBondingEdge2    <<","<<    dihedraltype <<"myBasePairIsTwoTransformForce"<<myBasePairIsTwoTransformForce<<endl;
    }; 
