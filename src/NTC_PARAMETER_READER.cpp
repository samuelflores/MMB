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
#include <MMBLogger.h>
#include <iostream>
#include <istream>
#include <fstream>
#include "SimTKsimbody.h"
#include "NtC_Class_Container.h"
#include "NTC_FORCE_CLASS.h"
#include "NTC_PARAMETER_READER.h"

using namespace SimTK;
using namespace std;

static map<NTC_PAR_BondKey, NTC_PAR_BondRow, NTC_PAR_BondKeyCmp> NTC_PAR_Map;

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

    NTC_PAR_BondKey::NTC_PAR_BondKey(const NTC_PAR_BondRow &myNTC_PAR_BondRow) {
        pdbResidueName1 = myNTC_PAR_BondRow.pdbResidueName1;
        pdbResidueName2 = myNTC_PAR_BondRow.pdbResidueName2;
        bondingEdge1 = myNTC_PAR_BondRow.bondingEdge1;
        bondingEdge2 =myNTC_PAR_BondRow.bondingEdge2;
        dihedraltype = myNTC_PAR_BondRow.dihedraltype;
        isTwoTransformForce = myNTC_PAR_BondRow.isTwoTransformForce; 
    }; 	
//};

    int NTC_PAR_Class::initialize(const String &inFileName) {
        NTC_PAR_Map.clear();
        myNTC_PAR_BondMatrix.myNTC_PAR_BondRow.clear();
        ifstream inFile(inFileName.c_str(),ifstream::in);
        MMBLOG_FILE_FUNC_LINE(DEBUG, "Now checking for existence of "<<inFileName<<endl);

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
        MMBLOG_FILE_FUNC_LINE(DEBUG, "done initializing myNTC_PAR_BondMatrix"<<endl);

        return(0);
    }

    void NTC_PAR_Class::printNTC_PAR_BondRows () {
        for (size_t q =0; q < myNTC_PAR_BondMatrix.myNTC_PAR_BondRow.size(); q++) {
            MMBLOG_FILE_FUNC_LINE(
                INFO,
                (myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).pdbResidueName1
                <<(myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).pdbResidueName2
                <<(myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).bondingEdge1
                <<(myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).bondingEdge2
                <<(myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).dihedraltype << endl
            );
        }
    }


    int NTC_PAR_Class::getNTC_PAR_BondRowIndex(
        const String &myPdbResidueName1,
        const String &myPdbResidueName2,
        const String &Classtype,
        const String &dihedraltype,
        const String &myBasePairIsTwoTransformForce,
        const NTC_Classes &NTC) const {
        MMBLOG_FILE_FUNC_LINE(
            INFO,
            myNTC_PAR_BondMatrix.myNTC_PAR_BondRow.size() << " br size " << " " << Classtype << " " << dihedraltype << " " << myPdbResidueName1 << " "<< myPdbResidueName2 << endl
        );

        for (size_t q =0; q < myNTC_PAR_BondMatrix.myNTC_PAR_BondRow.size(); q++) {
            const auto &srcNTC = myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q];

            if ((srcNTC.pdbResidueName1.compare(myPdbResidueName1) == 0) &&
                (srcNTC.pdbResidueName2.compare(myPdbResidueName2) == 0) &&
                (srcNTC.bondingEdge1.compare(Classtype) == 0) &&
                (srcNTC.bondingEdge2.compare(Classtype) == 0) &&
                (srcNTC.dihedraltype.compare(dihedraltype) == 0)) {
                return q;
            }
        }

	MMBLOG_FILE_FUNC_LINE(CRITICAL, "Found no match for the above user-specified interaction.  Either add this interaction type to the parameter file, or check your spelling, syntax, or semantics.\n");
    }

