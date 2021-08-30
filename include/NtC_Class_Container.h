/* -------------------------------------------------------------------------- *
 *                           MMB (MacroMoleculeBuilder)                       *
 * -------------------------------------------------------------------------- *
 *                                                                            *
 * Copyright (c) 2011-12 by the Author.                                       *
 * Author: Samuel Flores                                                      *
 *                                                                            *
 * See RNABuilder.cpp for the copyright and usage agreement.                  *
 * -------------------------------------------------------------------------- */

#ifndef NTC_Class_Container_H_
#define NTC_Class_Container_H_

#include "Utils.h"

class BiopolymerClassContainer;
class NTC_PAR_Class;

class MMB_EXPORT NTC_Class_Container {
public:
    void clear() ;
    void generateAorBFormNtCs(BiopolymerClassContainer &myBiopolymerClassContainer, const String &chainID, ResidueID firstResidue, ResidueID lastResidue, double myNTCWeight, const NTC_PAR_Class &ntc_par_class);
    void add_NTC_Class(BiopolymerClassContainer & myBiopolymerClassContainer,
                       const NTC_PAR_Class & ntc_par_class, const String &myChain,
                       const ResidueID firstNtCResidueInStretch, const ResidueID lastNtCResidueInStretch,
                       const String &NtCClassString, const double myNtCWeight, const bool myMeta, const double myNtCWeight2);
    void add_NTC_Class(BiopolymerClassContainer & myBiopolymerClassContainer,
                       const NTC_PAR_Class & ntc_par_class, NTC_Classes & NTC);
    void validate_NTC_Class(BiopolymerClassContainer & myBiopolymerClassContainer,
                            const NTC_PAR_Class & ntc_par_class, NTC_Classes & NTC,
                            const String &dihedraltype);
    int numNTC_Torsions();
    const NTC_Classes & getNTC_Class(int NTC_Class_Index);
    vector<NTC_Classes>	myNTC_Class_Vector;

private:
    void initMolmodelAtomIndices(const BiopolymerClassContainer &container, const NTC_PAR_Class &ntcPars, NTC_Classes &ntc);
};

#endif
