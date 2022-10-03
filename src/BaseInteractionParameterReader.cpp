/* -------------------------------------------------------------------------- *
 *                           MMB (MacroMoleculeBuilder)                       *
 * -------------------------------------------------------------------------- *
 *                                                                            *
 * Copyright (c) 2011-12 by the Author.                                       *
 * Author: Samuel Flores                                                      *
 *                                                                            *
 * See RNABuilder.cpp for the copyright and usage agreement.                  *
 * -------------------------------------------------------------------------- */

#include <MMBLogger.h>
#include <iostream>
#include <istream>
#include <fstream>
#include <ostream>
#include <set>
#include "SimTKsimbody.h"
#include "BaseInteractionParameterReader.h"
static constexpr int numLeontisWesthofBondMatrixRows=(26*16+4*4);  // yes, global constants are bad.  Couldn't think of an elegant way around this one though.  This number should be exactly equal to the number of rows in the leontisWesthofBondMatrix
//const int maxParallelTorques = 1000; //max number of parallel torques to be applied.  This can be huge, minimal cost for doing that.

std::ostream & operator<<(std::ostream &os, const LeontisWesthofBondKey &key) {
    os << "LeontisWesthofBondKey:\n"
       << "pdbResidueNames:\t" << "\"" << key.pdbResidueName1 << "\", \"" << key.pdbResidueName2 << "\"\n"
       << "bondingEdges:\t" << "\"" << key.bondingEdge1 << "\", \"" << key.bondingEdge2 << "\"\n"
       << "glycosidicBondOrientation:\t" << "\"" << key.glycosidicBondOrientation << "\"\n"
       << "isTwoTransformForce:\t" << "\"" << key.isTwoTransformForce << "\"\n";

    return os;
}

inline
double str_to_dbl(const std::string &s) {
    try {
        return std::stod(s);
    } catch (const std::invalid_argument &) {
        return 0.0;
    }
}

namespace std {
    template<> struct less<LeontisWesthofBondKey> {
        bool operator()(const LeontisWesthofBondKey &ti1, const LeontisWesthofBondKey &ti2) const {
            if      (ti1.pdbResidueName1 < ti2.pdbResidueName1) return 1;
            else if (ti1.pdbResidueName1 > ti2.pdbResidueName1) return 0;
            else if (ti1.pdbResidueName2 < ti2.pdbResidueName2) return 1;
            else if (ti1.pdbResidueName2 > ti2.pdbResidueName2) return 0;
            else if (ti1.bondingEdge1 < ti2.bondingEdge1) return 1;
            else if (ti1.bondingEdge1 > ti2.bondingEdge1) return 0;
            else if (ti1.bondingEdge2 < ti2.bondingEdge2) return 1;
            else if (ti1.bondingEdge2 > ti2.bondingEdge2) return 0;
            else if (ti1.glycosidicBondOrientation < ti2.glycosidicBondOrientation) return 1;
            else if (ti1.glycosidicBondOrientation > ti2.glycosidicBondOrientation) return 0;
            else if (ti1.isTwoTransformForce < ti2.isTwoTransformForce) return 1;
            else if (ti1.isTwoTransformForce > ti2.isTwoTransformForce) return 0;

	    return 0;
        }
    };
}

using namespace SimTK;
using namespace std;

bool LeontisWesthofBondKey::operator==(const LeontisWesthofBondKey &other) const noexcept {
    return
        (this->pdbResidueName1 == other.pdbResidueName1) &&
        (this->pdbResidueName2 == other.pdbResidueName2) &&
        (this->bondingEdge1    == other.bondingEdge1) &&
        (this->bondingEdge2    == other.bondingEdge2) &&
        (this->glycosidicBondOrientation == other.glycosidicBondOrientation) &&
	(this->isTwoTransformForce == other.isTwoTransformForce);
}

static map <LeontisWesthofBondKey, LeontisWesthofBondRow> leontisWesthofMap;

static const set<std::string> KnownBondingEdges{
    "",
    "Bifurcated",
    "Concentricity",
    "HelicalStackingA3",
    "HelicalStackingB3",
    "Hoogsteen",
    "ChiBondAnti",
    "InterResidueForce",
    "Parallelness",
    "PointToPlane",
    "Rigid",
    "Stacking",
    "Stacking3",
    "Stacking5",
    "SugarEdge",
    "SuperimposeC3p",
    "Superimpose",
    "WatsonCrick",
    "Weld",
    "Concentricity",
    "HelicalStackingA3",
    "HelicalStackingB3",
    "Hoogsteen",
    "ChiBondAnti",
    "InterResidueForce",
    "Parallelness",
    "PointToPlane",
    "Rigid",
    "Stacking",
    "Stacking3",
    "Stacking5",
    "SugarEdge",
    "SuperimposeC3p",
    "Superimpose",
    "WatsonCrick",
    "Weld",
    "HelicalStackingB3",
    "Superimpose",
    "WatsonCrick",
    "HelicalStackingB3",
    "Superimpose",
    "WatsonCrick",
    "HelicalStackingB3",
    "SamePolarityQuadStacking3",
    "Superimpose",
    "WatsonCrick",
    "HelicalStackingB3",
    "Superimpose",
    "WatsonCrick",
    "Concentricity",
    "HelicalStackingA3",
    "HelicalStackingB3",
    "Hoogsteen",
    "ChiBondAnti",
    "InterResidueForce",
    "Parallelness",
    "PointToPlane",
    "Rigid",
    "SamePolarityQuadStacking3",
    "SamePolarityQuadStacking5",
    "Stacking",
    "Stacking3",
    "Stacking5",
    "SugarEdge",
    "SuperimposeC3p",
    "Superimpose",
    "WatsonCrick",
    "Weld",
    "HardSphere",
    "ChiBond",
    "ProteinBackboneSterics",
    "SelectedAtoms",
    "SingleAtomNucleicAcid",
    "HelicalStackingB3",
    "WatsonCrick",
    "Concentricity",
    "HelicalStackingA3",
    "HelicalStackingA5",
    "HelicalStackingB3",
    "HelicalStackingB5",
    "Hoogsteen",
    "ChiBondAnti",
    "InterResidueForce",
    "Parallelness",
    "PointToPlane",
    "Rigid",
    "Stacking",
    "Stacking3",
    "Stacking5",
    "SugarEdge",
    "SuperimposeC3p",
    "Superimpose",
    "WatsonCrick",
    "Weld"
};

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

    void LeontisWesthofClass::initialize(const String &inFileName) {
        leontisWesthofMap.clear();
        myLeontisWesthofBondMatrix.myLeontisWesthofBondRow.clear();

        ifstream inFile(inFileName, ifstream::in);
        MMBLOG_FILE_FUNC_LINE(INFO, "Now checking for existence of "<<inFileName<<endl);

        if (!inFile.good()) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "Unable to open parameter file "<<inFileName<<endl);
        }
	MMBLOG_FILE_FUNC_LINE(INFO, "Parameter file " << inFileName << " opened sucessfully\n");

        string s;
        while (inFile.good()) {
            std::getline(inFile,s,',');

            if (s == "RECORD") { //if this is a RECORD entry
                LeontisWesthofBondRow row;

                std::getline(inFile,s,',');
	        row.pdbResidueName1 = String(s);

                std::getline(inFile,s,',');
	        row.pdbResidueName2 = String(s);

                std::getline(inFile,s,',');
		if (!isKnownBondingEdge(s)) {
		    MMBLOG_PLAIN(CRITICAL, "Bonding edge specified \"" << s << "\" is not on the list of known bonding edges");
		}
	        row.bondingEdge1 = String(s);

                std::getline(inFile,s,',');
		if (!isKnownBondingEdge(s)) {
		    MMBLOG_PLAIN(CRITICAL, "Bonding edge specified \"" << s << "\" is not on the list of known bonding edges");
		}
	        row.bondingEdge2 = String(s);

                std::getline(inFile,s,',');
	        row.glycosidicBondOrientation = String(s);

                for (int r=0; r<4; r++) {
                    std::getline(inFile,s,',');
	            row.residue1Atom[r] = String(s);
                }
                for (int r=0; r<4; r++) {
                    std::getline(inFile,s,',');
	            row.residue2Atom[r] = String(s);
                }
                for (int r=0; r<4; r++) {
                    std::getline(inFile,s,',');
	            row.bondLength[r] = str_to_dbl(s);
                }
                for (int r=0; r<4; r++) {
                    std::getline(inFile,s,',');
	            row.springConstant[r] = str_to_dbl(s);
                }

                std::getline(inFile,s,',');
	        row.torqueConstant = str_to_dbl(s);

		assert(row.torqueConstant >= 0);
                for (int r=0; r<3; r++) {
                    std::getline(inFile,s,',');
	            row.attachmentPoint[r] = str_to_dbl(s);
                }

                std::getline(inFile,s,',');
                row.rotationAngle = str_to_dbl(s);
                for (int r=0; r<3; r++) {
                    std::getline(inFile,s,',');
	            row.rotationAxis[r] = str_to_dbl(s);
                }

                std::getline(inFile,s,',');
                row.isTwoTransformForce = String(s);

                std::getline(inFile,s,',');
                row.distanceC1pC1p = str_to_dbl(s);

		myLeontisWesthofBondMatrix.myLeontisWesthofBondRow.push_back(std::move(row));
		size_t positionIdx = myLeontisWesthofBondMatrix.myLeontisWesthofBondRow.size() - 1;

		// NOTE: This creates two identical copies of the LeontisWesthofBondRow object
		auto copy = myLeontisWesthofBondMatrix.myLeontisWesthofBondRow[positionIdx];
                leontisWesthofMap[LeontisWesthofBondKey(copy)] = copy;
            }
        }

        SimTK_ERRCHK_ALWAYS(
            leontisWesthofMap.size() ==  myLeontisWesthofBondMatrix.myLeontisWesthofBondRow.size(),
            "[BaseInteractionParameterReader.cpp]","Inconsistency in number of Leontis-Westhof bond rows. This probably means that your parameter file tried to specify parameters for the same interaction twice!"
        );

        inFile.close();
        MMBLOG_FILE_FUNC_LINE(INFO, "done initializing myLeontisWesthofBondMatrix"<<endl);
    }


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

    LeontisWesthofBondRow LeontisWesthofClass::getLeontisWesthofBondRow(ResidueID myResidueNumber1, ResidueID myResidueNumber2, String myPdbResidueName1, String myBondingEdge1, String myPdbResidueName2, String myBondingEdge2, String myGlycosidicBondOrientation,String myBasePairIsTwoTransformForce) const {
        const auto iter = leontisWesthofMap.find(
            LeontisWesthofBondKey(
	        myPdbResidueName1, myPdbResidueName2,
		myBondingEdge1, myBondingEdge2,
		myGlycosidicBondOrientation, myBasePairIsTwoTransformForce
	    )
	);
        if (iter == leontisWesthofMap.end()) {
            MMBLOG_FILE_FUNC_LINE(
                CRITICAL,
                "Unable to find parameters for interaction \"" << myBasePairIsTwoTransformForce << "\" between residue type: \""<<myPdbResidueName1 <<"\" , and residue type \""<<myPdbResidueName2<<"\", parameter "<<myBondingEdge1<<", "<<myBondingEdge2<<", orientation "<<myGlycosidicBondOrientation<<" between residue numbers "<<myResidueNumber1.getResidueNumber()<<" and "<<myResidueNumber2.getResidueNumber()<<endl
            );
        }

        LeontisWesthofBondRow myReturnLeontisWesthofBondRow = iter->second;

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

bool isKnownBondingEdge(const std::string &edge) {
    return KnownBondingEdges.find(edge) != KnownBondingEdges.cend();
}
