#ifndef BaseInteractionParameterReader_H_
#define BaseInteractionParameterReader_H_
/* -------------------------------------------------------------------------- *
 *                           MMB (MacroMoleculeBuilder)                       *
 * -------------------------------------------------------------------------- *
 *                                                                            *
 * Copyright (c) 2011-12 by the Author.                                       *
 * Author: Samuel Flores                                                      *
 *                                                                            *
 * See RNABuilder.cpp for the copyright and usage agreement.                  *
 * -------------------------------------------------------------------------- */


#include "SimTKmolmodel.h"
//#include "SimTKsimbody_aux.h"
#include "SimTKsimbody.h"
#include <ios>
#include <fstream>
#include <vector>
#include <iostream>
#include <string>
#include "Utils.h"

using namespace std;
using namespace SimTK;

/**
 * 
 * 
 * /param 
 * myPdbResidueName1,2 must be one of "A","C","G","U".
 * bondingEdge1,2 must be one of "WatsonCrick","Hoogsteen","Sugar","Bifurcated".
 * glycosidicBondOrientation must be either "Cis" or "Trans".
 *
 */
struct LeontisWesthofBondRow   {
	String pdbResidueName1;
        String bondingEdge1; 
	String pdbResidueName2;
	String bondingEdge2;
	String glycosidicBondOrientation;
        String residue1Atom[4];
        String residue2Atom[4];
	double  bondLength[4];
	double  springConstant[4];
	double torqueConstant;
	Vec3   attachmentPoint;
	double  rotationAngle;
	Vec3   rotationAxis;
        String isTwoTransformForce;
        double distanceC1pC1p;
};

class LeontisWesthofBondKey {
    public:
	String pdbResidueName1;
	String pdbResidueName2;
        String bondingEdge1;
	String bondingEdge2;
	String glycosidicBondOrientation;
        String isTwoTransformForce;
        LeontisWesthofBondKey(
            String myPdbResidueName1, String myPdbResidueName2,
            String myBondingEdge1, String myBondingEdge2,
            String myGlycosidicBondOrientation, String myIsTwoTransformForce
        );
        LeontisWesthofBondKey(LeontisWesthofBondRow myLeontisWesthofBondRow);

        bool operator==(const LeontisWesthofBondKey &other) const noexcept;
};

struct LeontisWesthofBondMatrix {
	vector<LeontisWesthofBondRow> myLeontisWesthofBondRow;
};
class LeontisWesthofClass  { 
public:
    LeontisWesthofBondMatrix myLeontisWesthofBondMatrix;
    void initialize(const String &inFileName);
    Transform getLeontisWesthofTransform(LeontisWesthofBondRow myLeontisWesthofBondRow) const;

    LeontisWesthofBondRow getNearestLeontisWesthofBondRow(String myPdbResidueName1,  String myPdbResidueName2, Transform residue1Transform, Transform residue2Transform)  const;
    void printLeontisWesthofBondRows ();


    int  getLeontisWesthofBondRowIndex(
        String myPdbResidueName1,
        String myPdbResidueName2,
        String myBondingEdge1, 
        String myBondingEdge2,
        String myGlycosidicBondOrientation,
        String myBasePairIsTwoTransformForce
        ) const;
    LeontisWesthofBondRow getLeontisWesthofBondRow(ResidueID myResidueNumber1,ResidueID myResidueNumber2, String myPdbResidueName1, String myBondingEdge1, String myPdbResidueName2,String myBondingEdge2, String myGlycosidicBondOrientation,String myBasePairIsTwoTransformForce) const  ;

};

bool isKnownBondingEdge(const std::string &edge);

#endif //      BaseInteractionParameterReader_H_
