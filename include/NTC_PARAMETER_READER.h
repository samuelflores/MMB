// vim: set expandtab sw=4 ts=4 sts=4:

#ifndef NTC_PARAMETER_READER_H_
#define NTC_PARAMETER_READER_H_
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
#include <vector>
#include <string>
#include "Utils.h"
#include "NtC_Class_Container.h"

using namespace std;
using namespace SimTK;

/**
 *
 *
 * /param
 * myPdbResidueName1,2 must be one of "A","C","G","U".
 * bondingEdge1,2 must be one of "WatsonCrick","Hoogsteen","Sugar","Bifurcated".
 * dihedraltype must be either "Cis" or "Trans".
 *
 */
struct NTC_PAR_BondRow {
    String pdbResidueName1;
    String pdbResidueName2;
    String bondingEdge1;
    String bondingEdge2;
    String dihedraltype;
    String residue1Atom[4];
    String residue2Atom[4];
    int    atom_shift[4];
    String isTwoTransformForce;
    double bondLength[4];
    double springConstant[4];
    double CONFALVALUE;
    double torqueConstant;
    double distanceC1pC1p;
    double rotationAngle;
    Vec3   attachmentPoint;
    Vec3   rotationAxis;
    int    br;
};


inline
std::string ntcBondKey(
    const std::string &resName1, const std::string &resName2,
    const std::string &bondEdge1, const std::string &bondEdge2,
    const std::string &dihedralType,
    const std::string &isTwoTransformForce)
{
    return resName1 + resName2 +
	   bondEdge1 + bondEdge2 +
	   dihedralType +
	   isTwoTransformForce;
}

inline
std::string ntcBondKey(const NTC_PAR_BondRow &row) {
    return ntcBondKey(
        row.pdbResidueName1, row.pdbResidueName2,
	row.bondingEdge1, row.bondingEdge2,
	row.dihedraltype,
	row.isTwoTransformForce
    );
}

/*
class NTC_PAR_BondKey {
public:
    NTC_PAR_BondKey(String myPdbResidueName1, String myPdbResidueName2, String myBondingEdge1, String myBondingEdge2, String mydihedraltype, String myIsTwoTransformForce) noexcept;
    NTC_PAR_BondKey(const NTC_PAR_BondRow &myNTC_PAR_BondRow);

    String pdbResidueName1;
    String pdbResidueName2;
    String bondingEdge1;
    String bondingEdge2;
    String dihedraltype;
    String isTwoTransformForce;
};
*/

/*
struct NTC_PAR_BondKeyCmp {
    bool operator()(const NTC_PAR_BondKey &ti1, const NTC_PAR_BondKey &ti2) const {
        if (ti1.pdbResidueName1 < ti2.pdbResidueName1) return 1;
        else if (ti1.pdbResidueName1 > ti2.pdbResidueName1) return 0;
        else if ((ti1.pdbResidueName2 < ti2.pdbResidueName2)) return  1;
        else if ((ti1.pdbResidueName2 > ti2.pdbResidueName2)) return  0;
        else if ((ti1.bondingEdge1 < ti2.bondingEdge1)) return  1;
        else if ((ti1.bondingEdge1 > ti2.bondingEdge1)) return  0;
        else if ((ti1.bondingEdge2 < ti2.bondingEdge2)) return  1;
        else if ((ti1.bondingEdge2 > ti2.bondingEdge2)) return  0;
        else if ((ti1.dihedraltype < ti2.dihedraltype)) return  1;
        else if ((ti1.dihedraltype > ti2.dihedraltype)) return  0;
        else if ((ti1.isTwoTransformForce < ti2.isTwoTransformForce)) return  1;
        else if ((ti1.isTwoTransformForce > ti2.isTwoTransformForce)) return  0;
        else return 0;
    }
};
*/

struct NTC_PAR_BondMatrix {
	vector<NTC_PAR_BondRow> myNTC_PAR_BondRow;
};

class NTC_PAR_Class {
public:
    void initialize (const String &inFileName);
    Transform getNTC_PAR_Transform(const NTC_PAR_BondRow &myNTC_PAR_BondRow) const;
    NTC_PAR_BondRow getNearestNTC_PAR_BondRow(const String &myPdbResidueName1, const String &myPdbResidueName2, const Transform &residue1Transform, const Transform &residue2Transform)  const;
    void printNTC_PAR_BondRows();
    std::size_t getNTC_PAR_BondRowIndex(const std::string &key) const;
    std::size_t getNTC_PAR_BondRowIndex(const String &myPdbResidueName1, const String &myPdbResidueName2,
                                const String &Classtype, const String &dihedraltype, const String &myBasePairIsTwoTransformForce) const;

    NTC_PAR_BondMatrix myNTC_PAR_BondMatrix;
};

#endif //      BaseInteractionParameterReader_H_
