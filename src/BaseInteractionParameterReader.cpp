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
#include "BaseInteractionParameterReader.h"
const int numLeontisWesthofBondMatrixRows=(26*16+4*4);  // yes, global constants are bad.  Couldn't think of an elegant way around this one though.  This number should be exactly equal to the number of rows in the leontisWesthofBondMatrix
//const int maxParallelTorques = 1000; //max number of parallel torques to be applied.  This can be huge, minimal cost for doing that.

using namespace SimTK;
using namespace std;

/**
 * 
 * 
 * /param 
 * myPdbResidueName1,2 must be one of "A","C","G","U".
 * bondingEdge1,2 must be one of "WatsonCrick","Hoogsteen","Sugar","Bifurcated".
 * glycosidicBondOrientation must be either "Cis" or "Trans".
 *
 */

    LeontisWesthofBondKey::LeontisWesthofBondKey(String myPdbResidueName1, String myPdbResidueName2,String myBondingEdge1, String myBondingEdge2, String myGlycosidicBondOrientation, String myIsTwoTransformForce) {
        pdbResidueName1 = myPdbResidueName1;
        pdbResidueName2 = myPdbResidueName2;
        bondingEdge1 = myBondingEdge1;
        bondingEdge2 = myBondingEdge2;
        glycosidicBondOrientation = myGlycosidicBondOrientation;
        isTwoTransformForce = myIsTwoTransformForce; 
    };	
        LeontisWesthofBondKey::LeontisWesthofBondKey(LeontisWesthofBondRow myLeontisWesthofBondRow) {
        pdbResidueName1 = myLeontisWesthofBondRow.pdbResidueName1;
        pdbResidueName2 = myLeontisWesthofBondRow.pdbResidueName2;
        bondingEdge1 = myLeontisWesthofBondRow.bondingEdge1;
        bondingEdge2 =myLeontisWesthofBondRow.bondingEdge2;
        glycosidicBondOrientation = myLeontisWesthofBondRow.glycosidicBondOrientation;
        isTwoTransformForce = myLeontisWesthofBondRow.isTwoTransformForce; 
    }; 	
//};

    int  LeontisWesthofClass::initialize          ( String inFileName) {
        leontisWesthofMap.clear();
        myLeontisWesthofBondMatrix.myLeontisWesthofBondRow.clear();
        ifstream inFile(inFileName.c_str(),ifstream::in);
        MMBLOG_FILE_FUNC_LINE(INFO, "Now checking for existence of "<<inFileName<<endl);

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
             

            if ((String(s)).compare("RECORD") == 0)  { //if this is a RECORD entry
                std::getline(inFile,s,',');
            	//inFile.getline( s, 100,',' );
                LeontisWesthofBondRow tempRow;
                myLeontisWesthofBondMatrix.myLeontisWesthofBondRow.push_back(tempRow);
	        (myLeontisWesthofBondMatrix.myLeontisWesthofBondRow[q]).pdbResidueName1 = String(s);
                std::getline(inFile,s,',');
                //inFile.getline( s, 100,',' );
	        (myLeontisWesthofBondMatrix.myLeontisWesthofBondRow[q]).pdbResidueName2 = String(s);
                std::getline(inFile,s,',');
                //inFile.getline( s, 100,',' );
	        (myLeontisWesthofBondMatrix.myLeontisWesthofBondRow[q]).bondingEdge1    = String(s);
                std::getline(inFile,s,',');
                //inFile.getline( s, 100,',' );
	        (myLeontisWesthofBondMatrix.myLeontisWesthofBondRow[q]).bondingEdge2    = String(s);
                std::getline(inFile,s,',');
                //inFile.getline( s, 100,',' );
	        (myLeontisWesthofBondMatrix.myLeontisWesthofBondRow[q]).glycosidicBondOrientation = String(s);
;
                for (int r=0; r<4; r++) {
                    std::getline(inFile,s,',');
	            //inFile.getline( s, 100,',' );
	            (myLeontisWesthofBondMatrix.myLeontisWesthofBondRow[q]).residue1Atom[r] = String(s);
                } 
                for (int r=0; r<4; r++) {
                    std::getline(inFile,s,',');
	            //inFile.getline( s, 100,',' );
	            (myLeontisWesthofBondMatrix.myLeontisWesthofBondRow[q]).residue2Atom[r] = String(s);
                } 
                for (int r=0; r<4; r++) {
                    std::getline(inFile,s,',');
	            //inFile.getline( s, 100,',' );
	            (myLeontisWesthofBondMatrix.myLeontisWesthofBondRow[q]).bondLength[r] =   (double)atof(s.c_str());
                } 
                for (int r=0; r<4; r++) {
                    std::getline(inFile,s,',');
	            //inFile.getline( s, 100,',' );
	            (myLeontisWesthofBondMatrix.myLeontisWesthofBondRow[q]).springConstant[r] = (double)atof(s.c_str());
                } 
                std::getline(inFile,s,',');
	        //inFile.getline(s,100,',');
	        (myLeontisWesthofBondMatrix.myLeontisWesthofBondRow[q]).torqueConstant =  (double)atof(s.c_str());
		assert ((myLeontisWesthofBondMatrix.myLeontisWesthofBondRow[q]).torqueConstant >= 0);
                for (int r=0; r<3; r++) {
                    std::getline(inFile,s,',');
	            //inFile.getline( s, 100,',' );
	            (myLeontisWesthofBondMatrix.myLeontisWesthofBondRow[q]).attachmentPoint[r] = atof(s.c_str());
                }
                std::getline(inFile,s,',');
                //inFile.getline( s, 100,',' );
                (myLeontisWesthofBondMatrix.myLeontisWesthofBondRow[q]).rotationAngle = (double)atof(s.c_str());
                for (int r=0; r<3; r++) {
                    std::getline(inFile,s,',');
	            //inFile.getline( s, 100,',' );
	            (myLeontisWesthofBondMatrix.myLeontisWesthofBondRow[q]).rotationAxis[r] = atof(s.c_str());
                }

                std::getline(inFile,s,',');
                //inFile.getline( s, 100,',' );
                (myLeontisWesthofBondMatrix.myLeontisWesthofBondRow[q]).isTwoTransformForce= String(s);

                std::getline(inFile,s,',');
                //inFile.getline( s, 100,',' );
                (myLeontisWesthofBondMatrix.myLeontisWesthofBondRow[q]).distanceC1pC1p = (double)atof(s.c_str());

                leontisWesthofMap[LeontisWesthofBondKey(myLeontisWesthofBondMatrix.myLeontisWesthofBondRow[q])] = myLeontisWesthofBondMatrix.myLeontisWesthofBondRow[q];
	        q++;
            }
            }
        SimTK_ERRCHK_ALWAYS(leontisWesthofMap.size() ==  myLeontisWesthofBondMatrix.myLeontisWesthofBondRow.size() ,"[BaseInteractionParameterReader.cpp]","Inconsistency in number of Leontis-Westhof bond rows. This probably means that your parameter file tried to specify parameters for the same interaction twice!"); 
         
         
        //delete[] s;	
        inFile.close();
        MMBLOG_FILE_FUNC_LINE(INFO, "done initializing myLeontisWesthofBondMatrix"<<endl);

        return(0);
        };


    void LeontisWesthofClass::printLeontisWesthofBondRows () {    
        for   (int q =0; q< (int)myLeontisWesthofBondMatrix.myLeontisWesthofBondRow.size(); q++) 
            MMBLOG_FILE_FUNC_LINE(INFO, (myLeontisWesthofBondMatrix.myLeontisWesthofBondRow[q]).pdbResidueName1
                <<(myLeontisWesthofBondMatrix.myLeontisWesthofBondRow[q]).pdbResidueName2
                <<(myLeontisWesthofBondMatrix.myLeontisWesthofBondRow[q]).bondingEdge1
                <<(myLeontisWesthofBondMatrix.myLeontisWesthofBondRow[q]).bondingEdge2
                <<(myLeontisWesthofBondMatrix.myLeontisWesthofBondRow[q]).glycosidicBondOrientation
                <<(myLeontisWesthofBondMatrix.myLeontisWesthofBondRow[q]).isTwoTransformForce<<endl);
    };


    int  LeontisWesthofClass::getLeontisWesthofBondRowIndex(
        //int myResidueNumber1,
        //int myResidueNumber2,
        String myPdbResidueName1,
        String myPdbResidueName2,
        String myBondingEdge1, 
        String myBondingEdge2,
        String myGlycosidicBondOrientation,
        String myBasePairIsTwoTransformForce
        ) const {
        //if (0) { //!((myBasePairIsTwoTransformForce.compare("aromatic") == 0) || (myBasePairIsTwoTransformForce.compare("baseInteraction") == 0))) {
        //    return -11111;
        //} else  {
      
        for   (int q =0; q< (int)myLeontisWesthofBondMatrix.myLeontisWesthofBondRow.size(); q++) {
            if (
                ((((myLeontisWesthofBondMatrix.myLeontisWesthofBondRow[q]).pdbResidueName1).compare(myPdbResidueName1)) ==0)  &&
                ((((myLeontisWesthofBondMatrix.myLeontisWesthofBondRow[q]).pdbResidueName2).compare(myPdbResidueName2)) ==0)&&
                ((((myLeontisWesthofBondMatrix.myLeontisWesthofBondRow[q]).bondingEdge1).compare(myBondingEdge1))==0) &&
                ((((myLeontisWesthofBondMatrix.myLeontisWesthofBondRow[q]).bondingEdge2).compare(myBondingEdge2))==0) &&
                ((((myLeontisWesthofBondMatrix.myLeontisWesthofBondRow[q]).glycosidicBondOrientation).compare(myGlycosidicBondOrientation)) == 0) //&&
                //((((myLeontisWesthofBondMatrix.myLeontisWesthofBondRow[q]).isTwoTransformForce).compare(myBasePairIsTwoTransformForce) == (0) )  )
                )
                {
		    if (0) MMBLOG_FILE_FUNC_LINE(INFO, "found the right LeontisWesthofBondRow. residue1Atom[0], residue1Atom[1], residue1Atom[2] ,residue1Atom[3] residue2Atom[0], residue2Atom[1], residue2Atom[2] ,residue2Atom[3]  ="<<
                    (myLeontisWesthofBondMatrix.myLeontisWesthofBondRow[q]).residue1Atom[0]<<","<<
                    (myLeontisWesthofBondMatrix.myLeontisWesthofBondRow[q]).residue1Atom[1]<<","<<
                    (myLeontisWesthofBondMatrix.myLeontisWesthofBondRow[q]).residue1Atom[2]<<","<<
                    (myLeontisWesthofBondMatrix.myLeontisWesthofBondRow[q]).residue1Atom[3]<<","<<
                    (myLeontisWesthofBondMatrix.myLeontisWesthofBondRow[q]).residue2Atom[0]<<","<<
                    (myLeontisWesthofBondMatrix.myLeontisWesthofBondRow[q]).residue2Atom[1]<<","<<
                    (myLeontisWesthofBondMatrix.myLeontisWesthofBondRow[q]).residue2Atom[2]<<","<<
                    (myLeontisWesthofBondMatrix.myLeontisWesthofBondRow[q]).residue2Atom[3]<<
	            endl);
                    if (myBasePairIsTwoTransformForce.compare("contact") != 0) SimTK_ERRCHK_ALWAYS( fabs(((myLeontisWesthofBondMatrix.myLeontisWesthofBondRow[q]).rotationAxis).norm()-1.0) <.001 ,"[BaseInteractionParameterReader.cpp]","The desired interaction was found but the norm of its rotationAxis is not unity within tolerance of .001.  The interaction may be blank or incorrect.  "); 
                    
                    return q ; //myLeontisWesthofBondMatrix.myLeontisWesthofBondRow[q];
                }	 	
                    if (0) MMBLOG_FILE_FUNC_LINE(INFO, "looking at LeontisWesthofBondRow. residue1Atom[0], residue1Atom[1], residue1Atom[2] ,residue1Atom[3] residue2Atom[0], residue2Atom[1], residue2Atom[2] ,residue2Atom[3], bondingEdge1, bondingEdge2, glycosidicBondOrientation  ="<<
                    
                    ((((myLeontisWesthofBondMatrix.myLeontisWesthofBondRow[q]).pdbResidueName1)))<<","<<  
                    ((((myLeontisWesthofBondMatrix.myLeontisWesthofBondRow[q]).pdbResidueName2)))<<","<<
                    (myLeontisWesthofBondMatrix.myLeontisWesthofBondRow[q]).residue1Atom[0]<<","<<
                    (myLeontisWesthofBondMatrix.myLeontisWesthofBondRow[q]).residue1Atom[1]<<","<<
                    (myLeontisWesthofBondMatrix.myLeontisWesthofBondRow[q]).residue1Atom[2]<<","<<
                    (myLeontisWesthofBondMatrix.myLeontisWesthofBondRow[q]).residue1Atom[3]<<","<<
                    (myLeontisWesthofBondMatrix.myLeontisWesthofBondRow[q]).residue2Atom[0]<<","<<
                    (myLeontisWesthofBondMatrix.myLeontisWesthofBondRow[q]).residue2Atom[1]<<","<<
                    (myLeontisWesthofBondMatrix.myLeontisWesthofBondRow[q]).residue2Atom[2]<<","<<
                    (myLeontisWesthofBondMatrix.myLeontisWesthofBondRow[q]).residue2Atom[3]<<","<<

                    (myLeontisWesthofBondMatrix.myLeontisWesthofBondRow[q]).bondingEdge1   <<","<<
                    (myLeontisWesthofBondMatrix.myLeontisWesthofBondRow[q]).bondingEdge2   <<","<<
                    (myLeontisWesthofBondMatrix.myLeontisWesthofBondRow[q]).glycosidicBondOrientation<<","<<
	            endl);
                    if (0) MMBLOG_FILE_FUNC_LINE(INFO, "pdbResidueName1, pdbResidueName2, ="<<(myLeontisWesthofBondMatrix.myLeontisWesthofBondRow[q]).pdbResidueName1<<","<<(myLeontisWesthofBondMatrix.myLeontisWesthofBondRow[q]).pdbResidueName2<<endl);
                   
        }
	    MMBLOG_FILE_FUNC_LINE(INFO,
                "failed to match :"<<endl<<
                ","<<(myBasePairIsTwoTransformForce)<<
                ","<<(myPdbResidueName1)<<
                ","<<(myPdbResidueName2)<<
                ","<<(myBondingEdge1)<<
                ","<<(myBondingEdge2)<<
                ","<<(myGlycosidicBondOrientation)<<endl);
                //","<<(myGlycosidicBondOrientation)<<endl;  
            SimTK_ERRCHK_ALWAYS(0,"[BaseInteractionParameterReader.cpp]","Found no match for the above user-specified interaction.  Either add this interaction type to the parameter file, or check your spelling, syntax, or semantics.");

    //}
    };

    LeontisWesthofBondRow LeontisWesthofClass::getLeontisWesthofBondRow(ResidueID myResidueNumber1,ResidueID myResidueNumber2, String myPdbResidueName1, String myBondingEdge1, String myPdbResidueName2,String myBondingEdge2, String myGlycosidicBondOrientation,String myBasePairIsTwoTransformForce) const {

        static map <const LeontisWesthofBondKey, LeontisWesthofBondRow, LeontisWesthofBondKeyCmp>::iterator iter = leontisWesthofMap.begin();

        iter = leontisWesthofMap.find(LeontisWesthofBondKey(myPdbResidueName1, myPdbResidueName2, myBondingEdge1,  myBondingEdge2,  myGlycosidicBondOrientation,  myBasePairIsTwoTransformForce));
        
        LeontisWesthofBondRow myReturnLeontisWesthofBondRow;
        
        if (iter != leontisWesthofMap.end() ) 
            myReturnLeontisWesthofBondRow =  iter->second ;
        else {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "Unable to find parameters for interaction : "<< myBasePairIsTwoTransformForce<<" between residue type: \""<<myPdbResidueName1 <<"\" , and residue type \""<<myPdbResidueName2<<"\", parameter "<<myBondingEdge1<<", "<<myBondingEdge2<<", orientation "<<myGlycosidicBondOrientation<<" between residue numbers "<<myResidueNumber1.getResidueNumber()<<" and "<<myResidueNumber2.getResidueNumber()<<endl);
        }


        if 
            (!((myPdbResidueName1.compare(myReturnLeontisWesthofBondRow.pdbResidueName1) == 0) &&
            ( myPdbResidueName2.compare(myReturnLeontisWesthofBondRow.pdbResidueName2) == 0) &&            (myBondingEdge1.compare(myReturnLeontisWesthofBondRow.bondingEdge1) == 0)        &&
            (myBondingEdge2.compare(myReturnLeontisWesthofBondRow.bondingEdge2) == 0)        &&            (myGlycosidicBondOrientation.compare(myReturnLeontisWesthofBondRow.glycosidicBondOrientation) == 0) &&
            (myBasePairIsTwoTransformForce.compare(myReturnLeontisWesthofBondRow.isTwoTransformForce) == 0)))
                {

                MMBLOG_FILE_FUNC_LINE(INFO, "for interaction between residues "<<myResidueNumber1.getResidueNumber()<< " and "<<myResidueNumber2.getResidueNumber() <<endl<<"trying to match :"<<endl<<
                    ","<<(myBasePairIsTwoTransformForce)<<  
                    ", myPdbResidueName1"<<(myPdbResidueName1)<<
                    ", myPdbResidueName2"<<(myPdbResidueName2)<<
                    ", myBondingEdge1"<<(myBondingEdge1)<<
                    ", myBondingEdge2"<<(myBondingEdge2)<<
                    ", myGlycosidicBondOrientation"<<(myGlycosidicBondOrientation)<<endl);

                }

        return myReturnLeontisWesthofBondRow;
        // if (0) cout<<"Inside getLeontisWesthofBondRow.  about to search for :"<<myPdbResidueName1<<","<<myBondingEdge1<<","<< myPdbResidueName2   <<","<<  myBondingEdge2    <<","<<    myGlycosidicBondOrientation <<"myBasePairIsTwoTransformForce"<<myBasePairIsTwoTransformForce<<endl;
    };

