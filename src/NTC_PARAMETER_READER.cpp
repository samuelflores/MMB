// vim: expandtab sw=4 ts=4 sts=4 :
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
#include <fstream>
#include <stdexcept>
#include "SimTKsimbody.h"
#include "NtC_Class_Container.h"
#include "NTC_FORCE_CLASS.h"
#include "NTC_PARAMETER_READER.h"

using namespace SimTK;
using namespace std;

static map<std::string, std::size_t> NTC_PAR_Map;

/**
 *
 *
 * /param
 * myPdbResidueName1,2 must be one of "A","C","G","U".
 * bondingEdge1,2 must be one of "WatsonCrick","Hoogsteen","Sugar","Bifurcated".
 * dihedraltype must be either "Cis" or "Trans".
 *
 */

class RowInitializer {
public:
    RowInitializer(std::ifstream &ifs, NTC_PAR_BondRow &row) :
        m_ifs{ifs},
        m_row{row}
    {}

    void initField(String NTC_PAR_BondRow::*field) {
        read();
        m_row.*field = m_buf;
    }

    void initField(double NTC_PAR_BondRow::*field) {
        read();
        m_row.*field = getDouble(m_buf);
    }

    void initField(Vec3 NTC_PAR_BondRow::*field, size_t nVals = 3) {
        for (size_t idx = 0; idx < nVals; idx++) {
            read();
            (m_row.*field)[idx] = getDouble(m_buf);
        }
    }

    template <size_t Size>
    void initField(String (NTC_PAR_BondRow::*field)[Size], size_t idx) {
        read();
        (m_row.*field)[idx] = m_buf;
    }

    template <size_t Size>
    void initField(int (NTC_PAR_BondRow::*field)[Size], size_t idx) {
        read();
        (m_row.*field)[idx] = getInt(m_buf);
    }

    template <size_t Size>
    void initField(double (NTC_PAR_BondRow::*field)[Size], size_t idx) {
        read();
        (m_row.*field)[idx] = getDouble(m_buf);
    }

private:
    double getDouble(const std::string &str) {
        if (str.empty())
            return 0.0;
        try {
            return std::stod(str);
        } catch (const std::invalid_argument &) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "Cannot convert string \"" << m_buf << "\" to double" << std::endl);
        } catch (const std::out_of_range &) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "Numerical value of\"" << m_buf << "\" out of range of double" << std::endl);
        }
    }

    int getInt(const std::string &str) {
        if (str.empty())
            return 0;
        try {
            return std::stoi(str);
        } catch (const std::invalid_argument &) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "Cannot convert string \"" << m_buf << "\" to int" << std::endl);
        } catch (const std::out_of_range &) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "Numerical value of\"" << m_buf << "\" out of range of int" << std::endl);
        }
    }

    void read() {
        if (!m_ifs.good())
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "Cannot read NtC data from parameters file" << std::endl);
        std::getline(m_ifs, m_buf, ',');
    }

    std::ifstream &m_ifs;
    NTC_PAR_BondRow &m_row;
    std::string m_buf;
};

void NTC_PAR_Class::initialize(const String &inFileName) {
    if (!NTC_PAR_Map.empty() || !myNTC_PAR_BondMatrix.myNTC_PAR_BondRow.empty())
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "NTC parameters definitions are not empty. This means that parameters definition file has already been read!" << std::endl);

    ifstream inFile(inFileName.c_str(), ifstream::in);
    MMBLOG_FILE_FUNC_LINE(DEBUG, "Now checking for existence of " << inFileName << endl);

    if (!inFile.good()) {
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "Unable to open parameter file " << inFileName << endl);
    }

    /* We are reserving space for 10 more NtC classes than we
     * had as of writing this */
    myNTC_PAR_BondMatrix.myNTC_PAR_BondRow.reserve(233200);

    string s;
    while (inFile.good()) {
        std::getline(inFile, s, ',');

        if (s.compare("NTCRECORD") == 0) { // if this is a RECORD entry
            myNTC_PAR_BondMatrix.myNTC_PAR_BondRow.emplace_back(NTC_PAR_BondRow{});
            auto &row = myNTC_PAR_BondMatrix.myNTC_PAR_BondRow.back();

            RowInitializer ri{inFile, row};

            ri.initField(&NTC_PAR_BondRow::pdbResidueName1); // resname A, G, T, C
            ri.initField(&NTC_PAR_BondRow::pdbResidueName2); // resname A, G, T, C column 3
            ri.initField(&NTC_PAR_BondRow::bondingEdge1); // type NTC Class C 4
            ri.initField(&NTC_PAR_BondRow::bondingEdge2); // type NTC Class C 5
            ri.initField(&NTC_PAR_BondRow::dihedraltype); // type dihedral C 6

            for (size_t r = 0; r < 4; r++) {
                ri.initField(&NTC_PAR_BondRow::residue1Atom, r); // res atom 1 C7-C10
            }
            for (size_t r = 0; r < 4; r++) {
                ri.initField(&NTC_PAR_BondRow::atom_shift, r);
                if (!(row.atom_shift[r] == 0 || row.atom_shift[r] == 1))
                    MMBLOG_FILE_FUNC_LINE(CRITICAL, "Value of NTC atom_shift must be 0 or 1, got " << row.atom_shift[r] << std::endl);
            }
            for (size_t r = 0; r < 4; r++) {
                ri.initField(&NTC_PAR_BondRow::bondLength, r);
            }
            for (size_t r = 0; r < 4; r++) {
                ri.initField(&NTC_PAR_BondRow::springConstant, r);
            }

            ri.initField(&NTC_PAR_BondRow::torqueConstant);
            if (row.torqueConstant < 0.0)
                MMBLOG_FILE_FUNC_LINE(CRITICAL, "NtC torqueConstant must be positive but the file contains invalid value " << row.torqueConstant << std::endl);

            ri.initField(&NTC_PAR_BondRow::attachmentPoint);
            ri.initField(&NTC_PAR_BondRow::rotationAngle);
            ri.initField(&NTC_PAR_BondRow::rotationAxis, 2);
            ri.initField(&NTC_PAR_BondRow::CONFALVALUE);
            ri.initField(&NTC_PAR_BondRow::isTwoTransformForce);
            ri.initField(&NTC_PAR_BondRow::distanceC1pC1p);

            NTC_PAR_Map[ntcBondKey(row)] = myNTC_PAR_BondMatrix.myNTC_PAR_BondRow.size();
        }
    }

    if (NTC_PAR_Map.size() != myNTC_PAR_BondMatrix.myNTC_PAR_BondRow.size())
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "Inconsistency in number of NTC parameter rows. This probably means that your parameter file tried to specify parameters for the same interaction twice!" << std::endl);

    MMBLOG_FILE_FUNC_LINE(DEBUG, "Done initializing myNTC_PAR_BondMatrix" << endl);
}

void NTC_PAR_Class::printNTC_PAR_BondRows() {
  for (size_t q = 0; q < myNTC_PAR_BondMatrix.myNTC_PAR_BondRow.size(); q++) {
    MMBLOG_FILE_FUNC_LINE(
        INFO, (myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).pdbResidueName1
                  << (myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).pdbResidueName2
                  << (myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).bondingEdge1
                  << (myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).bondingEdge2
                  << (myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[q]).dihedraltype
                  << endl);
  }
}

std::size_t NTC_PAR_Class::getNTC_PAR_BondRowIndex(const std::string &key) const {
    auto it = NTC_PAR_Map.find(key);

    if (it == NTC_PAR_Map.end()) {
        MMBLOG_FILE_FUNC_LINE( CRITICAL, "Found no match for the above user-specified interaction.  Either add this interaction type to the parameter file, or check your spelling, syntax, or semantics." << std::endl);
    }

    return it->second;
}

std::size_t NTC_PAR_Class::getNTC_PAR_BondRowIndex(
    const String &myPdbResidueName1, const String &myPdbResidueName2,
    const String &Classtype, const String &dihedraltype,
    const String &myBasePairIsTwoTransformForce) const {
    MMBLOG_FILE_FUNC_LINE(INFO, myNTC_PAR_BondMatrix.myNTC_PAR_BondRow.size()
                                << " br size "
                                << " " << Classtype << " " << dihedraltype
                                << " " << myPdbResidueName1 << " "
                                << myPdbResidueName2 << endl);

    return getNTC_PAR_BondRowIndex(
        ntcBondKey(
	    myPdbResidueName1, myPdbResidueName2,
            Classtype, Classtype,
	    dihedraltype,
	    myBasePairIsTwoTransformForce
	)
    );
}
