// vim :set expandtab sw=4 ts=4 sts=4:

/* -------------------------------------------------------------------------- *
 *                           MMB (MacroMoleculeBuilder)                       *
 * -------------------------------------------------------------------------- *
 *                                                                            *
 * Copyright (c) 2011-12 by the Author.                                       *
 * Author: Samuel Flores                                                      *
 *                                                                            *
 * See RNABuilder.cpp for the copyright and usage agreement.                  *
 * -------------------------------------------------------------------------- */

#include "BiopolymerClass.h"
#include "NtC_Class_Container.h"
#include "NTC_PARAMETER_READER.h"
#include "MMBLogger.h"
#include "molmodel/internal/Compound.h"
#include <array>
#include <sstream>

static const std::array<String, 22> DIHEDRAL_TYPES = {
    "delta", "epsilon", "zeta", "alpha1", "beta1", "gamma1", "delta1", "chi", "chi1", "nccn", "nn", "cc",
    "tau0a", "tau1a", "tau2a", "tau3a", "tau4a", "tau0b", "tau1b", "tau2b", "tau3b", "tau4b"
};

void NTC_Class_Container::clear() {
    myNTC_Class_Vector.clear();
}

void NTC_Class_Container::generateAorBFormNtCs(
    BiopolymerClassContainer &myBiopolymerClassContainer, const String &chainID,
    ResidueID firstResidue, ResidueID lastResidue, double myNTCWeight,
    const NTC_PAR_Class &ntc_par_class) {
    if (myBiopolymerClassContainer.getBiopolymerClass(chainID).difference(lastResidue, firstResidue) < 1) {
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "It's not possible to apply helical stacking interactions to a run of fewer than 2 residues! " << endl);
    }

    String myNtCClassString = "XXXXZZZZ";
    if (myBiopolymerClassContainer.getBiopolymerClass(chainID).getBiopolymerType() == BiopolymerType::RNA) {
        myNtCClassString = "AA00";
    } else if (myBiopolymerClassContainer.getBiopolymerClass(chainID).getBiopolymerType() == BiopolymerType::DNA) {
        myNtCClassString = "BB00";
    } else {
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "Invalid biopolymerType: " << myBiopolymerClassContainer.getBiopolymerClass(chainID).getBiopolymerType() << endl);
    }

    add_NTC_Class(myBiopolymerClassContainer, ntc_par_class, chainID,
                  firstResidue, lastResidue, myNtCClassString, myNTCWeight, 0, 0);
}

void NTC_Class_Container::add_NTC_Class(
    BiopolymerClassContainer &myBiopolymerClassContainer,
    const NTC_PAR_Class &ntc_par_class, const String &myChain,
    const ResidueID firstNtCResidueInStretch, const ResidueID lastNtCResidueInStretch,
    const String &NtCClassString,
    const double myNtCWeight, const bool myMeta, const double myNtCWeight2) {

    MMBLOG_FILE_FUNC_LINE(ALWAYS, "Syntax: NtC <chain> <start residue> <end residue> <NtC class> "
                                  "<force constant> [meta <secondary weight>] "
                                  << endl
                                  << "For example, if (DNA) chain A, residues 1 and 2 are in a "
                                     "B-form helix helix, and you want a force constant of "
                                     "1.5, you can specify :  "
                                  << endl
                                  << "     NtC A 1 2 AA00 1.5  " << endl);

    if (myBiopolymerClassContainer.getBiopolymerClass(myChain).difference(firstNtCResidueInStretch, lastNtCResidueInStretch) != -1) {
        MMBLOG_FILE_FUNC_LINE(DEBUG, "NtCs could previously only be applied between consecutive residues." << endl);
    }

    int firstNtCResidueIndexInStretch = myBiopolymerClassContainer.getBiopolymerClass(myChain).getResidueIndex(firstNtCResidueInStretch);
    int lastNtCResidueIndexInStretch = myBiopolymerClassContainer.getBiopolymerClass(myChain).getResidueIndex(lastNtCResidueInStretch);

    for (int currentFirstResidueIndex = firstNtCResidueIndexInStretch; currentFirstResidueIndex < lastNtCResidueIndexInStretch; currentFirstResidueIndex += 1) {
        NTC_Classes NTC;
        NTC.NtC_FirstBPChain = myChain;
        NTC.NtC_Class_String = NtCClassString;
        // Convert ResidueID to index:
        // These three can be outside the loop:
        NTC.rotationCorrection1 = Rotation(0.0, UnitVec3(0, 0, 1));
        NTC.rotationCorrection2 = Rotation(0.0, UnitVec3(0, 0, 1));

        if (!(firstNtCResidueIndexInStretch < lastNtCResidueIndexInStretch)) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "Syntax error! The residue numbers must be ascending! You specified "
                                            << NTC.FirstBPResidue.outString() << " followed by "
                                            << NTC.SecondBPResidue.outString()
                                            << " . These have residue indices from "
                                            << firstNtCResidueIndexInStretch << " to "
                                            << lastNtCResidueIndexInStretch << endl);
        }

        MMBLOG_FILE_FUNC_LINE(INFO, "About to start NtC loop. Overall residue stretch is from "
                                    << firstNtCResidueInStretch.outString() << " to "
                                    << lastNtCResidueInStretch.outString() << " . " << endl);
        MMBLOG_FILE_FUNC_LINE(INFO, "Residue indices are from "
                                    << firstNtCResidueIndexInStretch << " to "
                                    << lastNtCResidueIndexInStretch << " . "
                                    << endl);

        MMBLOG_FILE_FUNC_LINE(INFO, "currentFirstResidueIndex = "
                                    << currentFirstResidueIndex << endl);
        MMBLOG_FILE_FUNC_LINE(INFO, "currentFirstResidueIndex = "
                                    << currentFirstResidueIndex << endl);
        MMBLOG_FILE_FUNC_LINE(INFO, "currentFirstResidueIndex + 1 = "
                                    << currentFirstResidueIndex + 1 << endl);

        if ((currentFirstResidueIndex + 1) > lastNtCResidueIndexInStretch) { // Don't see how this could happen, but being ultra paranoid.
            MMBLOG_FILE_FUNC_LINE(INFO, "The index: " << (currentFirstResidueIndex + 1)
                                        << " of the second residue in this NtC, is too high!"
                                        << endl);

            MMBLOG_FILE_FUNC_LINE(INFO, "Compared to the index: "
                                        << lastNtCResidueIndexInStretch
                                        << " of the last residue in the range."
                                        << endl);
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "Unexplained error! The index: "
                                            << (currentFirstResidueIndex + 1)
                                            << " of the second residue in this NtC, is too high! "
                                            << endl);
        }

        NTC.FirstBPResidue = myBiopolymerClassContainer.getBiopolymerClass(NTC.NtC_FirstBPChain).getResidueID(currentFirstResidueIndex);
        NTC.SecondBPResidue = myBiopolymerClassContainer.getBiopolymerClass(NTC.NtC_FirstBPChain).getResidueID(currentFirstResidueIndex + 1);
        NTC.NtC_step_ID = NTC.FirstBPResidue.ResidueNumber;

        MMBLOG_FILE_FUNC_LINE(INFO, "Starting NtC loop. Overall residue stretch is from "
                                    << firstNtCResidueInStretch.outString() << " to "
                                    << lastNtCResidueInStretch.outString()
                                    << " . In this round, NTC.FirstBPResidue = "
                                    << NTC.FirstBPResidue.outString()
                                    << " , NTC.SecondBPResidue = "
                                    << NTC.SecondBPResidue.outString() << " . " << endl);
        MMBLOG_FILE_FUNC_LINE(INFO, "NTC.FirstBPResidue = "
                                    << NTC.FirstBPResidue.outString() << endl);
        MMBLOG_FILE_FUNC_LINE(INFO, "NTC.SecondBPResidue = "
                                    << NTC.SecondBPResidue.outString() << endl);

        NTC.weight = myNtCWeight;
        MMBLOG_FILE_FUNC_LINE(INFO, "NTC.weight = " << NTC.weight << endl);

        if (NTC.weight > 20.0) { // Empirically found that a weight greater than 20 or so leads to strange jumpy nonconvergent behavior.
            MMBLOG_FILE_FUNC_LINE(WARNING, "You have specified an NtC weight of "
                                           << NTC.weight
                                           << " . The current NtC parameters are overconstrained, "
                                              "and unless these have been updated you may need to "
                                              "decrease the weight so as to avoid instability.  "
                                           << endl);
        }

        NTC.meta = 0;
        if (myMeta) {
            NTC.meta = 1;
            NTC.weight2 = myNtCWeight2;

            // These three do need to be in this loop:
            NTC.count = myBiopolymerClassContainer.count;
            myBiopolymerClassContainer.count++;

            MMBLOG_FILE_FUNC_LINE(INFO, NTC.count << " number of NTC meta input lines " << endl);
        }

        add_NTC_Class(myBiopolymerClassContainer, ntc_par_class, NTC);

        MMBLOG_FILE_FUNC_LINE(INFO, "At end of loop. Just added NTC with NTC.FirstBPResidue  = "
                                    << NTC.FirstBPResidue.outString()
                                    << " and NTC.SecondBPResidue = "
                                    << NTC.SecondBPResidue.outString() << endl);
        NTC.print();
    }
}

void NTC_Class_Container::add_NTC_Class(BiopolymerClassContainer &myBiopolymerClassContainer, const NTC_PAR_Class &ntc_par_class, NTC_Classes &NTC) {
    for (const auto &dht : DIHEDRAL_TYPES) {
        validate_NTC_Class(myBiopolymerClassContainer, ntc_par_class, NTC, dht);
        initMolmodelAtomIndices(myBiopolymerClassContainer, ntc_par_class, NTC);
        myNTC_Class_Vector.push_back(NTC);
    }
}

void NTC_Class_Container::validate_NTC_Class(BiopolymerClassContainer &myBiopolymerClassContainer, const NTC_PAR_Class &ntc_par_class, NTC_Classes &NTC, const String &dihedraltype) {
    int ntc2 = NTC.NtC_step_ID + 1;

    String resName1 = myBiopolymerClassContainer.getPdbResidueName(NTC.NtC_FirstBPChain, NTC.NtC_step_ID);
    String resName2 = myBiopolymerClassContainer.getPdbResidueName(NTC.NtC_FirstBPChain, ntc2);

    NTC.NTC_PAR_BondRowIndex  = ntc_par_class.getNTC_PAR_BondRowIndex(resName1,resName2,NTC.NtC_Class_String,dihedraltype,"ntcstep");

    MMBLOG_FILE_FUNC_LINE(INFO, NTC.NTC_PAR_BondRowIndex << " BOND ROW\n");
}

int NTC_Class_Container::numNTC_Torsions() { return myNTC_Class_Vector.size(); }

const NTC_Classes &NTC_Class_Container::getNTC_Class(int NTC_Class_Index) {
    return myNTC_Class_Vector[NTC_Class_Index];
}


void NTC_Class_Container::initMolmodelAtomIndices(BiopolymerClassContainer &container, const NTC_PAR_Class &ntcPars, NTC_Classes &ntc) {
    const auto bp = &container.getBiopolymerClass(ntc.NtC_FirstBPChain);
    const auto &bondRow = ntcPars.myNTC_PAR_BondMatrix.myNTC_PAR_BondRow[ntc.NTC_PAR_BondRowIndex];

    std::array<const ResidueID *const, 2> residues = { &ntc.FirstBPResidue, &ntc.SecondBPResidue };
    if (bondRow.bondLength[0] == 0.0) { // Bond of zero length is actually a torsion (I know, I know...)
        for (size_t idx = 0; idx < 4; idx++) {
            const ResidueID *resId = residues[bondRow.atom_shift[idx]];
            ntc.atomIndices[idx] = bp->atomIndex(*resId, bondRow.residue1Atom[idx]);
        }
    } else { // Non-zero bond length indicates an acutal bond
        for (size_t idx = 0; idx < 2; idx++) {
            ntc.atomIndices[idx] = bp->atomIndex(*residues[idx], bondRow.residue1Atom[idx]);
        }
    }
}
