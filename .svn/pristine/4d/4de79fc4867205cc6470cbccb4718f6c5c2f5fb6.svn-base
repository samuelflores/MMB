#ifndef SimTK_MOLMODEL_PROTEIN_H_
#define SimTK_MOLMODEL_PROTEIN_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Molmodel                            *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2007 Stanford University and the Authors.           *
 * Authors: Michael Sherman, Christopher Bruns                                *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */


#include "molmodel/internal/common.h"
#include "molmodel/internal/Compound.h"
#include "molmodel/internal/CompoundSystem.h"

namespace SimTK {

/// Widely used acetyl protein N-terminal end cap
class AcetylResidue : public BiopolymerResidue {
public:
    AcetylResidue() : BiopolymerResidue("acetyl", "ACE", 'X')
    {
        setBaseCompound("methyl", MethylGroup());
        nameAtom("CH3", "methyl/C");
        nameAtom("1H", "methyl/H1");
        nameAtom("2H", "methyl/H2");
        nameAtom("3H", "methyl/H3");

        convertInboardBondCenterToOutboard();
        bondAtom(TrivalentAtom("C", Element::Carbon()), "methyl/bond", 0.15520);
        bondAtom(UnivalentAtom("O", Element::Oxygen()), "C/bond2", 0.12290);
        nameBondCenter("bondC", "C/bond3");
        nameBondCenter("bondNext", "bondC");
        
        setAtomBiotype("CH3", "Acetyl", "CH3", SimTK::Ordinality::Initial);
        setAtomBiotype("1H",  "Acetyl", "H", SimTK::Ordinality::Initial);
        setAtomBiotype("2H",  "Acetyl", "H", SimTK::Ordinality::Initial);
        setAtomBiotype("3H",  "Acetyl", "H", SimTK::Ordinality::Initial);
        setAtomBiotype("C",   "Acetyl", "C", SimTK::Ordinality::Initial);
        setAtomBiotype("O",   "Acetyl", "O", SimTK::Ordinality::Initial);    
    }
};


/// Neutral C-terminal protein cap
class NMethylAmideResidue : public BiopolymerResidue {
public:
    NMethylAmideResidue() : BiopolymerResidue("N-methyl amide", "NAC", 'X')
    {
        // TODO - set nitrogen bond angles
        setBaseAtom(TrivalentAtom("N", Element::Nitrogen()));
        bondAtom(UnivalentAtom("HN", Element::Hydrogen()), "N/bond3", 0.1010);
        bondAtom(AliphaticCarbon("CH3"), "N/bond2", 0.1449);
        bondAtom(AliphaticHydrogen("1H"), "CH3/bond2");
        bondAtom(AliphaticHydrogen("2H"), "CH3/bond3");
        bondAtom(AliphaticHydrogen("3H"), "CH3/bond4");

        nameBondCenter("bondN", "N/bond1");
        nameBondCenter("bondPrevious", "bondN");

        // Define default bonding geometry to previous residue
        setDefaultInboardBondLength(0.1335);

        // Help with automated biotype resolution
        addCompoundSynonym("N-MeAmide");

        // zero and 180 degrees are poor choices for phi angle because each eclipses either N-H or N-C
        defineDihedralAngle("phi", "bondPrevious", "CH3/bond2", 30*Deg2Rad);
        
        setAtomBiotype("N",   "N-MeAmide", "N",   SimTK::Ordinality::Final);
        setAtomBiotype("HN",  "N-MeAmide", "HN",  SimTK::Ordinality::Final);
        setAtomBiotype("CH3", "N-MeAmide", "CH3", SimTK::Ordinality::Final);
        setAtomBiotype("1H",  "N-MeAmide", "H",   SimTK::Ordinality::Final);
        setAtomBiotype("2H",  "N-MeAmide", "H",   SimTK::Ordinality::Final);
        setAtomBiotype("3H",  "N-MeAmide", "H",   SimTK::Ordinality::Final);    
    }
};

static const Angle DefaultAlphaPhiAngle = -57 * Deg2Rad;
static const Angle DefaultAlphaPsiAngle = -47 * Deg2Rad;
static const Angle DefaultParallelBetaPhiAngle = -119 * Deg2Rad;
static const Angle DefaultParallelBetaPsiAngle =  113 * Deg2Rad;
static const Angle DefaultAntiparallelBetaPhiAngle = -139 * Deg2Rad;
static const Angle DefaultAntiparallelBetaPsiAngle =  135 * Deg2Rad;

/** 
 \brief amino acid residue building block for protein polypeptide chain molecules

 AminoAcidResidue has three unsatisfied BondCenters:
  1) "bondN" at the amino nitrogen, for the preceding amino acid residue
  2) "boncC" at the carbonyl carbon, for the next amino acid residue
  3) "bondCA" at the alpha carbon, for the side chain
 
 The decision of which atoms bond to bond centers bond1, bond2, bond3 etc., 
 depends upon several criteria:
   1) chirality - in quadravalent atoms (e.g. alpha carbon), bonds 2, 3, 4
     are aranged clockwise when viewed down bond 1.
   2) inboard bond - unless explicitly changed, bond1 is the inboard bond
     for a given compound
   3) dihedral angle definition - the dihedral angle for bond1 is defined
     relative to the plane containing bond2.  bonds2,3,4... are defined
     relative to bond1

 */
class SimTK_MOLMODEL_EXPORT AminoAcidResidue : public BiopolymerResidue {
public:

    /// Factory method for creating caninical amino acid residues
    static AminoAcidResidue create(
        char oneLetterCode ///< one letter code of desired amino acid residue, e.g. "W" for tryptophan
        );


    /// Factory method for creating caninical amino acid residues
    static AminoAcidResidue create(
        const PdbResidue& ///< PdbResidue object from Pdb.h
        );
        
    explicit AminoAcidResidue(
        String name, ///< name of this type of amino acid.  e.g. "glycine"
        String threeLetterCode = "Unk", ///< three letter code.  Might appear in ResidueName field of PDB files
        char oneLetterCode = '?' ///< one letter code (use 'X' for non-canonical residues)
        )
        : BiopolymerResidue(name, threeLetterCode, oneLetterCode)
    {
        static const mdunits::Length CA_Hdistance = 0.1090;
        static const mdunits::Length N_CAdistance = 0.1449;
        static const mdunits::Length CA_Cdistance = 0.1522;
        static const mdunits::Length  C_Odistance = 0.1229;
        static const mdunits::Length  C_Ndistance = 0.1335; // peptide bond

        // Amino nitrogen
        // We want the inboard bond center to be on the main chain nitrogen peptide bond
        // But we also want the CA and main chain H to be the reference for dihedrals.
        // To reconcile these two requirements, a bit of extra work is required.
        TrivalentAtom nAtom("N", Element::Nitrogen());
        // nAtom.convertInboardBondCenterToOutboard(); // free bond1 for use with HN
        // nAtom.setInboardBondCenter("bond3");  // lowest priority bond becomes inboard

        setBaseAtom(nAtom);

        // Alpha carbon
        bondAtom(QuadrivalentAtom("CA", Element::Carbon()), "N/bond2", N_CAdistance, DefaultAntiparallelBetaPhiAngle + 180*Deg2Rad); // phi
        bondAtom(UnivalentAtom("HA", Element::Hydrogen()), "CA/bond3", CA_Hdistance);
        nameBondCenter("bondCA", "CA/bond4");

        // Carbonyl
        bondAtom(TrivalentAtom("C", Element::Carbon()), "CA/bond2", CA_Cdistance, DefaultAntiparallelBetaPsiAngle + 180*Deg2Rad); // psi
        bondAtom(UnivalentAtom("O", Element::Oxygen()), "C/bond2", C_Odistance);
        nameBondCenter("bondC", "C/bond3");
        nameBondCenter("bondNext", "bondC");

        nameBondCenter("bondN", "N/bond1");
        nameBondCenter("bondPrevious", "bondN");

        // Define default bonding geometry to previous residue
        setDefaultInboardBondLength(C_Ndistance);

        defineDihedralAngle("psi", "N", "CA", "C", "O", 180*Deg2Rad);
        setDefaultPsiAngle(DefaultAntiparallelBetaPsiAngle);

        defineDihedralAngle("phi", "N/bond3", "CA/bond2", 180*Deg2Rad);
        setDefaultPhiAngle(DefaultAntiparallelBetaPhiAngle);
    }


    /**
     * \brief Set default (initial) phi dihedral angle (rotation about N-CA bond)
     *
     * Phi and psi angles are offset 180 degrees from other dihedral definitions
     */
    AminoAcidResidue& setDefaultPhiAngle(Angle phi) {
        setDefaultDihedralAngle("phi", phi);

        return *this;
    }

    AminoAcidResidue& setDefaultPsiAngle(Angle psi) {
        setDefaultDihedralAngle("psi", psi);

        return *this;
    }

    // Make intermediate classes inner classes of AminoAcidResidue
    // to avoid registration order problems in boost.python/Py++ wrapping
    class HNAminoAcidResidue;
    class BetaUnbranchedAminoAcidResidue;
    class BetaBranchedAminoAcidResidue;

    class Alanine;
    class Cysteine;
    class Aspartate;
    class Glutamate;
    class Phenylalanine;
    class Glycine;
    class Histidine;
    class Isoleucine;
    class Lysine;
    class Leucine;
    class Methionine;
    class Asparagine;
    class Proline;
    class Glutamine;
    class Arginine;
    class Serine;
    class Threonine;
    class Valine;
    class Tryptophan;
    class Tyrosine;
};


/**
    AminoAcidResidue differs from HNAminoAcidResidue in lacking "HN" proton, so Proline can
    be derived.
 */
class  AminoAcidResidue::HNAminoAcidResidue : public AminoAcidResidue {
public:
    explicit HNAminoAcidResidue(String name, String threeLetterCode = "Unk", char oneLetterCode = '?')
        : AminoAcidResidue(name, threeLetterCode, oneLetterCode)
    {
        static const mdunits::Length  H_Ndistance = 0.1010;

        bondAtom(UnivalentAtom("H", Element::Hydrogen()), "N/bond3", H_Ndistance);
        nameAtom("HN", "H"); // synonym

        // defineDihedralAngle("phi", "H", "N", "CA", "C", 180*Deg2Rad); // already defined
    }
};


// BetaunbranchedAminoAcidResidue is base class for AminoAcidResidues
// posessing side chains that do not branch at the beta carbon
class  AminoAcidResidue::BetaUnbranchedAminoAcidResidue : public AminoAcidResidue::HNAminoAcidResidue {
public:
    BetaUnbranchedAminoAcidResidue(String name, String tlc, char olc)
        : AminoAcidResidue::HNAminoAcidResidue(name, tlc, olc) 
    {
        bondCompound("beta methylene", MethyleneGroup(), "bondCA");
        nameAtom("CB", "beta methylene/C");
        nameAtom("1HB", "beta methylene/H1");
        nameAtom("2HB", "beta methylene/H2");
        nameBondCenter("bondCB", "beta methylene/bond2");

        defineDihedralAngle("chi1", "CA/bond1", "beta methylene/bond2");
        setDefaultDihedralAngle("chi1", -60*Deg2Rad); // best for most amino acids
    }
};


// BetaunbranchedAminoAcidResidue is base class for AminoAcidResidues
// posessing side chains that do not branch at the beta carbon
class  AminoAcidResidue::BetaBranchedAminoAcidResidue : public AminoAcidResidue::HNAminoAcidResidue {
public:
    BetaBranchedAminoAcidResidue(String name, String tlc, char olc)
        : AminoAcidResidue::HNAminoAcidResidue(name, tlc, olc) 
    {
        bondAtom(AliphaticCarbon("CB"), "bondCA");

        // To get correct chirality *and* priority in isoleucine and threonine,
        // place hydrogen on bond3

        bondAtom(AliphaticHydrogen("HB"), "CB/bond3");
        nameBondCenter("bondCB1", "CB/bond2");
        nameBondCenter("bondCB2", "CB/bond4");

        defineDihedralAngle("chi1", "CA/bond1", "CB/bond2");
        setDefaultDihedralAngle("chi1", -60*Deg2Rad); // best for most amino acids
    }
};


class  AminoAcidResidue::Alanine : public AminoAcidResidue::BetaUnbranchedAminoAcidResidue {
public:
    Alanine() 
        : AminoAcidResidue::BetaUnbranchedAminoAcidResidue("alanine", "Ala", 'A')
    {
        bondAtom(AliphaticHydrogen("3HB"), "bondCB");
        
        setAtomBiotype("N",   "Alanine", "N");
        setAtomBiotype("H",   "Alanine", "HN");
        setAtomBiotype("CA",  "Alanine", "CA");
        setAtomBiotype("C",   "Alanine", "C");
        setAtomBiotype("O",   "Alanine", "O");
        setAtomBiotype("HA",  "Alanine", "HA");
        setAtomBiotype("CB",  "Alanine", "CB");
        setAtomBiotype("1HB", "Alanine", "HB");    
        setAtomBiotype("2HB", "Alanine", "HB");    
        setAtomBiotype("3HB", "Alanine", "HB");    
    }
};


class  AminoAcidResidue::Cysteine : public AminoAcidResidue::BetaUnbranchedAminoAcidResidue {
public:
    Cysteine() 
        : AminoAcidResidue::BetaUnbranchedAminoAcidResidue("cysteine", "Cys", 'C')
    {
        // From amber99 parameters in Tinker param file amber99.dat
        static const mdunits::Length H_Sdistance = 0.1336;
        static const mdunits::Length S_Cdistance = 0.1810;
        static const Angle H_S_Cangle = 96.0*Deg2Rad;
        static const Angle S_C_Cangle = 108.6*Deg2Rad;

        bondAtom(BivalentAtom("SG", Element::Sulfur(), H_S_Cangle), "bondCB", S_Cdistance, 180*Deg2Rad);
        bondAtom(UnivalentAtom("HG", Element::Hydrogen()), "SG/bond2", H_Sdistance, 180*Deg2Rad);

        defineDihedralAngle("chi2", "HG", "SG", "CB", "CA");

        setDefaultBondAngle(108.60*Deg2Rad, "SG", "CB", "CA");

        // Synonym hints to help resolve biotypes from tinker files
        // eventually, these biotypes will be pre-defined and hard coded
        addCompoundSynonym("CYS (-SH)");
        addCompoundSynonym("cysteine (-SH)");
    }
};

class  AminoAcidResidue::Aspartate : public AminoAcidResidue::BetaUnbranchedAminoAcidResidue {
public:
    Aspartate() 
        : AminoAcidResidue::BetaUnbranchedAminoAcidResidue("aspartate", "Asp", 'D')
    {
        bondCompound("sidechain carboxylate", CarboxylateGroup(), "bondCB");
        nameAtom("CG", "sidechain carboxylate/C");
        nameAtom("OD1", "sidechain carboxylate/O1");
        nameAtom("OD2", "sidechain carboxylate/O2");

        defineDihedralAngle("chi2", "OD1", "CG", "CB", "CA");

        // From Roland Dunbrack's side chain rotamer statistics
        // http://dunbrack.fccc.edu/bbdep/confanalysis.php
        setDefaultDihedralAngle("chi2", -15*Deg2Rad);

        // Synonym hints to help resolve biotypes from tinker files
        // eventually, these biotypes will be pre-defined and hard coded
        addCompoundSynonym("aspartic acid");
    }
};

class  AminoAcidResidue::Glutamate : public AminoAcidResidue::BetaUnbranchedAminoAcidResidue {
public:
    Glutamate() 
        : AminoAcidResidue::BetaUnbranchedAminoAcidResidue("glutamate", "Glu", 'E')
    {
        bondAtom(AliphaticCarbon("CG"), "bondCB");
        bondAtom(AliphaticHydrogen("1HG"), "CG/bond3");
        bondAtom(AliphaticHydrogen("2HG"), "CG/bond4");

        bondCompound("sidechain carboxylate", CarboxylateGroup(), "CG/bond2");
        nameAtom("CD", "sidechain carboxylate/C");
        nameAtom("OE1", "sidechain carboxylate/O1");
        nameAtom("OE2", "sidechain carboxylate/O2");

        defineDihedralAngle("chi2", "CD", "CG", "CB", "CA");
        defineDihedralAngle("chi3", "OE1", "CD", "CG", "CB");

        setDefaultDihedralAngle("chi3", -15*Deg2Rad);

        // Synonym hints to help resolve biotypes from tinker files
        // eventually, these biotypes will be pre-defined and hard coded
        addCompoundSynonym("glutamic acid");
    }
};


class  AminoAcidResidue::Phenylalanine : public AminoAcidResidue::BetaUnbranchedAminoAcidResidue {
public:
    Phenylalanine()
        : AminoAcidResidue::BetaUnbranchedAminoAcidResidue("phenylalanine", "Phe", 'F')
    {
        static const mdunits::Length CB_CGdistance = 0.1510;
        static const mdunits::Length C_Hdistance = 0.1080; // from tinker amber99.dat
        static const mdunits::Length C_Cdistance = 0.1400; // from tinker amber99.dat

        // default angle for TrivalentAtom is 120 degrees, which is what we want
        // for a six-membered ring
        bondAtom( TrivalentAtom("CG", Element::Carbon()), "bondCB", CB_CGdistance, 180*Deg2Rad);

        bondAtom( TrivalentAtom("CD1", Element::Carbon()), "CG/bond2", C_Cdistance, 180*Deg2Rad );
        bondAtom( UnivalentAtom("HD1", Element::Hydrogen()), "CD1/bond3", C_Hdistance );

        // Start using hidedral of zero, to loop ring back to beginning
        bondAtom( TrivalentAtom("CE1", Element::Carbon()), "CD1/bond2", C_Cdistance, 0*Deg2Rad );
        bondAtom( UnivalentAtom("HE1", Element::Hydrogen()), "CE1/bond3", C_Hdistance );

        bondAtom( TrivalentAtom("CZ", Element::Carbon()), "CE1/bond2", C_Cdistance, 0*Deg2Rad );
        bondAtom( UnivalentAtom("HZ", Element::Hydrogen()), "CZ/bond3", C_Hdistance );

        bondAtom( TrivalentAtom("CE2", Element::Carbon()), "CZ/bond2", C_Cdistance, 0*Deg2Rad );
        bondAtom( UnivalentAtom("HE2", Element::Hydrogen()), "CE2/bond3", C_Hdistance );

        bondAtom( TrivalentAtom("CD2", Element::Carbon()), "CE2/bond2", C_Cdistance, 0*Deg2Rad );
        bondAtom( UnivalentAtom("HD2", Element::Hydrogen()), "CD2/bond3", C_Hdistance );

        addRingClosingBond("CD2/bond2", "CG/bond3", C_Cdistance, 0*Deg2Rad);

        defineDihedralAngle("chi2", "CD1", "CG", "CB", "CA");

        setDefaultDihedralAngle("chi2", 90*Deg2Rad); // from Dunbrack site

        setBondMobility(BondMobility::Rigid, "CG", "CD1");
        setBondMobility(BondMobility::Rigid, "CD1", "CE1");
        setBondMobility(BondMobility::Rigid, "CE1", "CZ");
        setBondMobility(BondMobility::Rigid, "CZ", "CE2");
        setBondMobility(BondMobility::Rigid, "CE2", "CD2");
        setBondMobility(BondMobility::Rigid, "CD2", "CG");
    }
};


class  AminoAcidResidue::Glycine : public AminoAcidResidue::HNAminoAcidResidue {
public:
    Glycine() 
        : AminoAcidResidue::HNAminoAcidResidue("glycine", "Gly", 'G')
    {
        bondAtom(AliphaticHydrogen("2HA"), "bondCA");
        nameAtom("1HA", "HA");  // glycine names alpha hydrogen differently than others
    }
};


class  AminoAcidResidue::Histidine : public AminoAcidResidue::BetaUnbranchedAminoAcidResidue {
public:
    Histidine()
        : AminoAcidResidue::BetaUnbranchedAminoAcidResidue("histidine", "His", 'H')
    {
        // use average values for different protonation states
        static const mdunits::Length CB_CGdistance = 0.15040; // from tinker amber99.dat
        static const mdunits::Length CG_NDdistance = 0.13895; // average of protonation states 1385/1394
        static const mdunits::Length N_CEdistance = 0.13390; // average of protonation states 1343/1335
        static const mdunits::Length NE_CDdistance = 0.13875; // avg 1394/1381
        static const mdunits::Length CD_CGdistance = 0.1373; // avg 1371/1375
        static const mdunits::Length H_Cdistance = 0.1080;
        static const mdunits::Length H_Ndistance = 0.1010;

        // set angle to pentagonal angle 108, as amber parameters set 120 degrees as
        // "optimal" and let the forces fight it out.  Such a strategy will not
        // do for simple coordinate generation and for rigid body models.
        static const Angle ringAngle = 108.000 * Deg2Rad; // 108 degrees is exact
        static const Angle outerRingAngle = 180*Deg2Rad - 0.5*ringAngle;

        bondAtom( TrivalentAtom("CG", Element::Carbon(), outerRingAngle, outerRingAngle), "bondCB", CB_CGdistance);

        bondAtom( TrivalentAtom("ND1", Element::Nitrogen(), ringAngle, outerRingAngle), "CG/bond2", CG_NDdistance, 180*Deg2Rad);
        // optional proton depends upon protonation state
        bondAtom( UnivalentAtom("HD1", Element::Hydrogen()), "ND1/bond3", H_Ndistance);

        bondAtom( TrivalentAtom("CE1", Element::Carbon(), ringAngle, outerRingAngle), "ND1/bond2", N_CEdistance, 0*Deg2Rad);
        bondAtom( UnivalentAtom("HE1", Element::Hydrogen()), "CE1/bond3", H_Cdistance);

        bondAtom( TrivalentAtom("NE2", Element::Nitrogen(), ringAngle, outerRingAngle), "CE1/bond2", CG_NDdistance, 0*Deg2Rad);
        // optional proton depends upon protonation state
        bondAtom( UnivalentAtom("HE2", Element::Hydrogen()), "NE2/bond3", H_Ndistance);

        bondAtom( TrivalentAtom("CD2", Element::Carbon(), ringAngle, outerRingAngle), "NE2/bond2", N_CEdistance, 0*Deg2Rad);
        bondAtom( UnivalentAtom("HD2", Element::Hydrogen()), "CD2/bond3", H_Cdistance);

        addRingClosingBond("CD2/bond2", "CG/bond3", CD_CGdistance, 0*Deg2Rad);

        defineDihedralAngle("chi2", "ND1", "CG", "CB", "CA");

        setDefaultDihedralAngle("chi2", 90*Deg2Rad); // from Dunbrack site

        // Biotypes for histidine are tricky
        // Since both acidic protons are added above, this residue defaults to positive charge
        // There is more work to do to handle HIS protonation properly TODO

        // Synonym hints to help resolve biotypes from tinker files
        // eventually, these biotypes will be pre-defined and hard coded
        addCompoundSynonym("histidine (+)");

        // Make ring rigid for internal coordinate simulation
        setBondMobility(BondMobility::Rigid, "CG", "ND1");
        setBondMobility(BondMobility::Rigid, "ND1", "CE1");
        setBondMobility(BondMobility::Rigid, "CE1", "NE2");
        setBondMobility(BondMobility::Rigid, "NE2", "CD2");
        setBondMobility(BondMobility::Rigid, "CD2", "CG");
    }
};


// Unlike the case with Threonine, the chirality of Isoleucine comes out
// correctly when built up from BetaBranchedAminoAcidResidue
class  AminoAcidResidue::Isoleucine : public AminoAcidResidue::BetaBranchedAminoAcidResidue {
public:
    Isoleucine()
        : AminoAcidResidue::BetaBranchedAminoAcidResidue("isoleucine", "Ile", 'I')
    {
        bondCompound("gamma methylene", MethyleneGroup(), "bondCB1");
        nameAtom("CG1", "gamma methylene/C");
        nameAtom("1HG1", "gamma methylene/H1");
        nameAtom("2HG1", "gamma methylene/H2");

        bondCompound("delta methyl", MethylGroup(), "gamma methylene/bond2");
        nameAtom("CD", "delta methyl/C");
        nameAtom("1HD", "delta methyl/H1");
        nameAtom("2HD", "delta methyl/H2");
        nameAtom("3HD", "delta methyl/H3");

        bondCompound("gamma methyl", MethylGroup(), "bondCB2");
        nameAtom("CG2", "gamma methyl/C");
        nameAtom("1HG2", "gamma methyl/H1");
        nameAtom("2HG2", "gamma methyl/H2");
        nameAtom("3HG2", "gamma methyl/H3");

        defineDihedralAngle("chi2", "CD", "CG1", "CB", "CA");

        setDefaultDihedralAngle("chi1", -70*Deg2Rad);
        setDefaultDihedralAngle("chi2", 170*Deg2Rad);
    }
};


class  AminoAcidResidue::Lysine : public AminoAcidResidue::BetaUnbranchedAminoAcidResidue {
public:
    Lysine()
        : AminoAcidResidue::BetaUnbranchedAminoAcidResidue("lysine", "Lys", 'K') 
    {
        bondAtom(AliphaticCarbon("CG"), "bondCB");
        bondAtom(AliphaticHydrogen("1HG"), "CG/bond3");
        bondAtom(AliphaticHydrogen("2HG"), "CG/bond4");

        bondAtom(AliphaticCarbon("CD"), "CG/bond2");
        bondAtom(AliphaticHydrogen("1HD"), "CD/bond3");
        bondAtom(AliphaticHydrogen("2HD"), "CD/bond4");

        bondAtom(AliphaticCarbon("CE"), "CD/bond2");
        bondAtom(AliphaticHydrogen("1HE"), "CE/bond3");
        bondAtom(AliphaticHydrogen("2HE"), "CE/bond4");

        bondCompound("zeta amine", PrimaryAmineGroup(), "CE/bond2");
        nameAtom("NZ", "zeta amine/N");
        nameAtom("1HZ", "zeta amine/H1");
        nameAtom("2HZ", "zeta amine/H2");
        nameAtom("3HZ", "zeta amine/H3");

        defineDihedralAngle("chi2", "CD", "CG", "CB", "CA");
        defineDihedralAngle("chi3", "CE", "CD", "CG", "CB");
        defineDihedralAngle("chi4", "NZ", "CE", "CD", "CG");
        defineDihedralAngle("chi5", "1HZ", "NZ", "CE", "CD");
    }
};

    
class  AminoAcidResidue::Leucine : public AminoAcidResidue::BetaUnbranchedAminoAcidResidue {
public:
    Leucine()
        : AminoAcidResidue::BetaUnbranchedAminoAcidResidue("leucine", "Leu", 'L')
    {
        bondAtom(AliphaticCarbon("CG"), "bondCB");
        bondAtom(AliphaticHydrogen("HG"), "CG/bond4");

        // It actually takes more lines of code to add a MethylGroup and rename the 
        // atoms, than it does to add the atoms individually

        // First delta methyl group
        bondAtom(AliphaticCarbon("CD1"), "CG/bond2");
        bondAtom(AliphaticHydrogen("1HD1"), "CD1/bond2");
        bondAtom(AliphaticHydrogen("2HD1"), "CD1/bond3");
        bondAtom(AliphaticHydrogen("3HD1"), "CD1/bond4");

        // Second delta methyl group
        bondAtom(AliphaticCarbon("CD2"), "CG/bond3");
        bondAtom(AliphaticHydrogen("1HD2"), "CD2/bond2");
        bondAtom(AliphaticHydrogen("2HD2"), "CD2/bond3");
        bondAtom(AliphaticHydrogen("3HD2"), "CD2/bond4");

        defineDihedralAngle("chi2", "CD1", "CG", "CB", "CA");

        setDefaultDihedralAngle("chi1", -70*Deg2Rad);
        setDefaultDihedralAngle("chi2", 170*Deg2Rad);
    }
};


class  AminoAcidResidue::Methionine : public AminoAcidResidue::BetaUnbranchedAminoAcidResidue {
public:
    Methionine()
        : AminoAcidResidue::BetaUnbranchedAminoAcidResidue("methionine", "Met", 'M') 
    {
        bondAtom(AliphaticCarbon("CG"), "bondCB");
        bondAtom(AliphaticHydrogen("1HG"), "CG/bond3");
        bondAtom(AliphaticHydrogen("2HG"), "CG/bond4");

        bondAtom( BivalentAtom("SD", Element::Sulfur(), 98.90*Deg2Rad), "CG/bond2", 0.1810, 180*Deg2Rad);

        bondAtom(AliphaticCarbon("CE"), "SD/bond2", 0.1810, 180*Deg2Rad);
        bondAtom(AliphaticHydrogen("1HE"), "CE/bond2");
        bondAtom(AliphaticHydrogen("2HE"), "CE/bond3");
        bondAtom(AliphaticHydrogen("3HE"), "CE/bond4");

        defineDihedralAngle("chi2", "SD", "CG", "CB", "CA");
        defineDihedralAngle("chi3", "CE", "SD", "CG", "CB");
        defineDihedralAngle("chi4", "1HE", "CE", "SD", "CG");
    }
};

    
class  AminoAcidResidue::Asparagine : public AminoAcidResidue::BetaUnbranchedAminoAcidResidue {
public:
    Asparagine() 
        : AminoAcidResidue::BetaUnbranchedAminoAcidResidue("asparagine", "Asn", 'N')
    {
        bondAtom(TrivalentAtom("CG", Element::Carbon(), 120.4*Deg2Rad, 116.6*Deg2Rad), "bondCB", 0.15220);

        bondAtom(UnivalentAtom("OD1", Element::Oxygen()), "CG/bond2", 0.12290);

        bondAtom(TrivalentAtom("ND2", Element::Nitrogen()), "CG/bond3", 0.13350);
        bondAtom(UnivalentAtom("1HD2", Element::Hydrogen()), "ND2/bond2", 0.1010);
        bondAtom(UnivalentAtom("2HD2", Element::Hydrogen()), "ND2/bond3", 0.1010);

        defineDihedralAngle("chi2", "OD1", "CG", "CB", "CA");

        // From Roland Dunbrack's side chain rotamer statistics
        // http://dunbrack.fccc.edu/bbdep/confanalysis.php
        setDefaultDihedralAngle("chi2", -20*Deg2Rad);

        // For internal coordinate simulation, make amide group rigid
        setBondMobility(BondMobility::Rigid, "CG", "ND2");
    }
};


class  AminoAcidResidue::Proline : public AminoAcidResidue {
public:
    Proline() 
        : AminoAcidResidue("proline", "Pro", 'P')
    {
        // TODO check chirality of all these hydrogen names (in all residues)

        // beta methylene
        bondAtom(AliphaticCarbon("CB"), "bondCA");
        bondAtom(AliphaticHydrogen("1HB"), "CB/bond3");
        bondAtom(AliphaticHydrogen("2HB"), "CB/bond4");

        // gamma methylene
        bondAtom(AliphaticCarbon("CG"), "CB/bond2");
        bondAtom(AliphaticHydrogen("1HG"), "CG/bond3");
        bondAtom(AliphaticHydrogen("2HG"), "CG/bond4");

        // delta methylene
        bondAtom(AliphaticCarbon("CD"), "CG/bond2");
        bondAtom(AliphaticHydrogen("1HD"), "CD/bond3");
        bondAtom(AliphaticHydrogen("2HD"), "CD/bond4");

        addRingClosingBond("CD/bond2", "N/bond3", 0.14490);

        // defineDihedralAngle("phi", "CD", "N", "CA", "C", 180*Deg2Rad);

        defineDihedralAngle("chi1", "CG", "CB", "CA", "N");
        defineDihedralAngle("chi2", "CD", "CG", "CB", "CA");
        defineDihedralAngle("chi3", "N", "CD", "CG", "CB");
        defineDihedralAngle("chi4", "CA", "N", "CD", "CG"); // ring closing bond...

        // Ring strain decreases bond angles in the side chain
        // 102.9 chosen to match bond distance N-CD in beta down configuration
        // cmbruns
        setDefaultBondAngle(102.9*Deg2Rad, "N", "CA", "CB");
        setDefaultBondAngle(102.9*Deg2Rad, "CA", "CB", "CG");
        setDefaultBondAngle(102.9*Deg2Rad, "CB", "CG", "CD");
        setDefaultBondAngle(102.9*Deg2Rad, "CG", "CD", "N");

        setDefaultDownPucker();
        setDefaultPhiAngle(-70*Deg2Rad);

        // Fix phi angle for internal coordinate simulation
        setBondMobility(BondMobility::Rigid, "N", "CA");
    }

    // Up and down configurations are about equally abundant,
    // except in cis-proline, where down is preferred.
    // 
    // Ho, B.K.; Coutsias, E.A.; Seok, C.; & Dill, K.A. (2005) 
    // "The flexibility in the proline ring couples to the protein
    // backbone"  Protein Science: 14:1011-1018

    Proline& setDefaultUpPucker() { // red
        // Angles are estimated from figures in Ho et al
        // TODO get more accurate angles
        setDefaultDihedralAngle("chi1", -25*Deg2Rad);
        setDefaultDihedralAngle("chi2", 42*Deg2Rad);
        setDefaultDihedralAngle("chi3", -38*Deg2Rad);
        setDefaultDihedralAngle("chi4", 20*Deg2Rad);

        return *this;
    }

    Proline& setDefaultDownPucker() { // yellow, CG-endo
        // Angles from gaussian03 6-31G** restricted Hartree-Fock
        // minimization in beta main chain configuration
        setDefaultDihedralAngle("chi1", 32.8*Deg2Rad);
        setDefaultDihedralAngle("chi2", -35.9*Deg2Rad);
        setDefaultDihedralAngle("chi3", 25.0*Deg2Rad);
        setDefaultDihedralAngle("chi4", -16.9*Deg2Rad);

        // Angles are estimated from figures in Ho et al
        //setDefaultDihedralAngle("chi1", 30*Deg2Rad);
        //setDefaultDihedralAngle("chi2", -42*Deg2Rad);
        //setDefaultDihedralAngle("chi3", 38*Deg2Rad);
        //setDefaultDihedralAngle("chi4", -20*Deg2Rad);

        return *this;
    }

};


class  AminoAcidResidue::Glutamine : public AminoAcidResidue::BetaUnbranchedAminoAcidResidue {
public:
    Glutamine() 
        : AminoAcidResidue::BetaUnbranchedAminoAcidResidue("glutamine", "Gln", 'Q')
    {
        bondAtom( AliphaticCarbon("CG"), "bondCB" );
        bondAtom( AliphaticHydrogen("1HG"), "CG/bond3" ); 
        bondAtom( AliphaticHydrogen("2HG"), "CG/bond4" ); 

        bondAtom(TrivalentAtom("CD", Element::Carbon(), 120.4*Deg2Rad, 116.6*Deg2Rad), "CG/bond2", 0.15220);

        bondAtom(UnivalentAtom("OE1", Element::Oxygen()), "CD/bond2", 0.12290);

        bondAtom(TrivalentAtom("NE2", Element::Nitrogen()), "CD/bond3", 0.13350);
        bondAtom(UnivalentAtom("1HE2", Element::Hydrogen()), "NE2/bond2", 0.1010);
        bondAtom(UnivalentAtom("2HE2", Element::Hydrogen()), "NE2/bond3", 0.1010);

        defineDihedralAngle("chi2", "CD", "CG", "CB", "CA");
        defineDihedralAngle("chi3", "OE1", "CD", "CG", "CB");

        setDefaultDihedralAngle("chi3", -20*Deg2Rad);

        // For internal coordinate simulation, make amide group rigid
        setBondMobility(BondMobility::Rigid, "CD", "NE2");

    }
};


class  AminoAcidResidue::Arginine : public AminoAcidResidue::BetaUnbranchedAminoAcidResidue {
public:
    Arginine()
        : AminoAcidResidue::BetaUnbranchedAminoAcidResidue("arginine", "Arg", 'R') 
    {
        bondAtom(AliphaticCarbon("CG"), "bondCB");
        bondAtom(AliphaticHydrogen("1HG"), "CG/bond3");
        bondAtom(AliphaticHydrogen("2HG"), "CG/bond4");

        bondAtom(AliphaticCarbon("CD"), "CG/bond2");
        bondAtom(AliphaticHydrogen("1HD"), "CD/bond3");
        bondAtom(AliphaticHydrogen("2HD"), "CD/bond4");

        // TODO - set bond angle CB/CD/NE to 111.2

        bondAtom(TrivalentAtom("NE", Element::Nitrogen(), 123.2*Deg2Rad, 118.4*Deg2Rad), "CD/bond2", 0.1436, 180*Deg2Rad);
        bondAtom(UnivalentAtom("HE", Element::Hydrogen()), "NE/bond3", 0.1010);

        bondAtom(TrivalentAtom("CZ", Element::Carbon()), "NE/bond2", 0.1340, 180*Deg2Rad);

        bondAtom(TrivalentAtom("NH1", Element::Nitrogen()), "CZ/bond2", 0.1340, 180*Deg2Rad);
        bondAtom(UnivalentAtom("1HH1", Element::Hydrogen()), "NH1/bond2", 0.1010);
        bondAtom(UnivalentAtom("2HH1", Element::Hydrogen()), "NH1/bond3", 0.1010);

        bondAtom(TrivalentAtom("NH2", Element::Nitrogen()), "CZ/bond3", 0.1340, 180*Deg2Rad);
        bondAtom(UnivalentAtom("1HH2", Element::Hydrogen()), "NH2/bond2", 0.1010);
        bondAtom(UnivalentAtom("2HH2", Element::Hydrogen()), "NH2/bond3", 0.1010);

        defineDihedralAngle("chi2", "CD", "CG", "CB", "CA");
        defineDihedralAngle("chi3", "NE", "CD", "CG", "CB");
        defineDihedralAngle("chi4", "CZ", "NE", "CD", "CG");
        defineDihedralAngle("chi5", "NH1", "CZ", "NE", "CD");

        // TODO research default dihedral for chi4

        // Make guanidinium group rigid for internal coordinate simulation (resonance)
        setBondMobility(BondMobility::Rigid, "CZ", "NH1");
        setBondMobility(BondMobility::Rigid, "CZ", "NH2");
        setBondMobility(BondMobility::Rigid, "NE", "CZ");
    }
};

    
class  AminoAcidResidue::Serine : public AminoAcidResidue::BetaUnbranchedAminoAcidResidue {
public:
    Serine() 
        : AminoAcidResidue::BetaUnbranchedAminoAcidResidue("serine", "Ser", 'S')
    {
        // side chain hydroxyl group
        bondCompound("oh", AlcoholOHGroup(), "bondCB");
        nameAtom("OG", "oh/O", Biotype::SerineOG().getIndex() );
        nameAtom("HG", "oh/H", Biotype::SerineHG().getIndex() );

        defineDihedralAngle("chi2", "HG", "OG", "CB", "CA");

        // Set atom biotypes
        setBiotypeIndex("N", Biotype::SerineN().getIndex() );
        setBiotypeIndex("HN", Biotype::SerineHN().getIndex() );
        setBiotypeIndex("CA", Biotype::SerineCA().getIndex() );
        setBiotypeIndex("HA", Biotype::SerineHA().getIndex() );
        setBiotypeIndex("C", Biotype::SerineC().getIndex() );
        setBiotypeIndex("O", Biotype::SerineO().getIndex() );
    }
};

// Threonine has a chirality about the beta carbon, so care
// must be observed when placing groups
// Thus derive from AminoAcidResidue instead of from BetaBranchedAminoAcidResidue
class  AminoAcidResidue::Threonine : public AminoAcidResidue::BetaBranchedAminoAcidResidue {
public:
    Threonine()
        : AminoAcidResidue::BetaBranchedAminoAcidResidue("threonine", "Thr", 'T')
    {
        // hydroxyl has precedence over methyl and hydrogen
        // w.r.t. dihedral definition
        bondCompound("gamma hydroxyl", AlcoholOHGroup(), "bondCB1");
        nameAtom("OG1", "gamma hydroxyl/O");
        nameAtom("HG1", "gamma hydroxyl/H");

        // place methyl and hydrogen so as to get correct chirality
        // i.e. clockwise rotation of OG/H/CG when viewed down CA/CB bond

        bondCompound("gamma methyl", MethylGroup(), "bondCB2");
        nameAtom("CG2", "gamma methyl/C");
        nameAtom("1HG2", "gamma methyl/H1");
        nameAtom("2HG2", "gamma methyl/H2");
        nameAtom("3HG2", "gamma methyl/H3");
    }
};


// To get atom precedence correct, build custom side chain,
// rather than deriving from BetaBranchedAminoAcidResidue
class  AminoAcidResidue::Valine : public AminoAcidResidue::HNAminoAcidResidue {
public:
    Valine()
        : AminoAcidResidue::HNAminoAcidResidue("valine", "Val", 'V')
    {
        bondAtom(AliphaticCarbon("CB"), "bondCA");
        bondAtom(AliphaticHydrogen("HB"), "CB/bond4");

        bondCompound("gamma methyl 1", MethylGroup(), "CB/bond2");
        nameAtom("CG1", "gamma methyl 1/C");
        nameAtom("1HG1", "gamma methyl 1/H1");
        nameAtom("2HG1", "gamma methyl 1/H2");
        nameAtom("3HG1", "gamma methyl 1/H3");

        bondCompound("gamma methyl 2", MethylGroup(), "CB/bond3");
        nameAtom("CG2", "gamma methyl 2/C");
        nameAtom("1HG2", "gamma methyl 2/H1");
        nameAtom("2HG2", "gamma methyl 2/H2");
        nameAtom("3HG2", "gamma methyl 2/H3");

        defineDihedralAngle("chi1", "CG1", "CB", "CA", "N");

        // definition of CG1 vs. CG2 makes preferred chi1 180 instead of -60 for valine
        setDefaultDihedralAngle("chi1", 180*Deg2Rad);
    }
};


class  AminoAcidResidue::Tryptophan : public AminoAcidResidue::BetaUnbranchedAminoAcidResidue {
public:
    Tryptophan()
        : AminoAcidResidue::BetaUnbranchedAminoAcidResidue("tryptophan", "Trp", 'W')
    {
        bondAtom( TrivalentAtom("CG", Element::Carbon(), 125.0*Deg2Rad, 128.6*Deg2Rad), "bondCB", 0.14950);

        // 125.65 angle computed to place hydrogen symmetrically
        bondAtom( TrivalentAtom("CD1", Element::Carbon(), 108.7*Deg2Rad, 125.65*Deg2Rad), "CG/bond2", 0.13520, 180*Deg2Rad);
        bondAtom( UnivalentAtom("HD1", Element::Hydrogen()), "CD1/bond3", 0.1080);

        // 124.25 computed to center hydrogen
        bondAtom( TrivalentAtom("NE1", Element::Nitrogen(), 111.50*Deg2Rad, 124.25*Deg2Rad), "CD1/bond2", 0.13810, 0*Deg2Rad);
        bondAtom( UnivalentAtom("HE1", Element::Hydrogen()), "NE1/bond3", 0.1010);

        bondAtom( TrivalentAtom("CE2", Element::Carbon(), 104.4*Deg2Rad, 132.80*Deg2Rad), "NE1/bond2", 0.1380, 0*Deg2Rad);

        // 180 degrees to bulge back out for the six-membered ring
        bondAtom( TrivalentAtom("CZ2", Element::Carbon()), "CE2/bond3", 0.1400, 180*Deg2Rad);
        bondAtom( UnivalentAtom("HZ2", Element::Hydrogen()), "CZ2/bond3", 0.1080);

        bondAtom( TrivalentAtom("CH2", Element::Carbon()), "CZ2/bond2", 0.1400, 0*Deg2Rad);
        bondAtom( UnivalentAtom("HH2", Element::Hydrogen()), "CH2/bond3", 0.1080);



        bondAtom( TrivalentAtom("CD2", Element::Carbon(), 108.8*Deg2Rad, 122.7*Deg2Rad), "CE2/bond2", 0.1419, 0*Deg2Rad);
        addRingClosingBond("CD2/bond2", "CG/bond3", 0.1459);

        bondAtom( TrivalentAtom("CE3", Element::Carbon()), "CD2/bond3", 0.1404, 0*Deg2Rad);
        bondAtom( UnivalentAtom("HE3", Element::Hydrogen()), "CE3/bond3", 0.1080);

        bondAtom( TrivalentAtom("CZ3", Element::Carbon()), "CE3/bond2", 0.1400, 0*Deg2Rad);
        bondAtom( UnivalentAtom("HZ3", Element::Hydrogen()), "CZ3/bond3", 0.1080);

        addRingClosingBond("CZ3/bond2", "CH2/bond2", 0.1400);

        defineDihedralAngle("chi2", "CD1", "CG", "CB", "CA");

        setDefaultDihedralAngle("chi2", 90*Deg2Rad);

        setBondMobility(BondMobility::Rigid, "CG", "CD1");
        setBondMobility(BondMobility::Rigid, "CD1", "NE1");
        setBondMobility(BondMobility::Rigid, "NE1", "CE2");
        setBondMobility(BondMobility::Rigid, "CZ2", "CE2");
        setBondMobility(BondMobility::Rigid, "CH2", "CZ2");
        setBondMobility(BondMobility::Rigid, "CD2", "CE3");
        setBondMobility(BondMobility::Rigid, "CZ3", "CE3");
        setBondMobility(BondMobility::Rigid, "CZ3", "CH2");
    }
};


class  AminoAcidResidue::Tyrosine : public AminoAcidResidue::BetaUnbranchedAminoAcidResidue {
public:
    Tyrosine()
        : AminoAcidResidue::BetaUnbranchedAminoAcidResidue("tyrosine", "Tyr", 'Y')
    {
        static const mdunits::Length CB_CGdistance = 0.1510;
        static const mdunits::Length C_Hdistance = 0.1080; // from tinker amber99.dat
        static const mdunits::Length C_Cdistance = 0.1400; // from tinker amber99.dat
        static const mdunits::Length CZ_Odistance = 0.1364; // from tinker amber99.dat

        // default angle for TrivalentAtom is 120 degrees, which is what we want
        // for a six-membered ring
        bondAtom( TrivalentAtom("CG", Element::Carbon()), "bondCB", CB_CGdistance, 180*Deg2Rad);

        bondAtom( TrivalentAtom("CD1", Element::Carbon()), "CG/bond2", C_Cdistance, 180*Deg2Rad );
        bondAtom( UnivalentAtom("HD1", Element::Hydrogen()), "CD1/bond3", C_Hdistance );

        // Start using hidedral of zero, to loop ring back to beginning
        bondAtom( TrivalentAtom("CE1", Element::Carbon()), "CD1/bond2", C_Cdistance, 0*Deg2Rad );
        bondAtom( UnivalentAtom("HE1", Element::Hydrogen()), "CE1/bond3", C_Hdistance );

        bondAtom( TrivalentAtom("CZ", Element::Carbon()), "CE1/bond2", C_Cdistance, 180*Deg2Rad );
        bondCompound( "oh", AlcoholOHGroup(), "CZ/bond2", CZ_Odistance );
        nameAtom("OH", "oh/O");
        nameAtom("HH", "oh/H");

        bondAtom( TrivalentAtom("CE2", Element::Carbon()), "CZ/bond3", C_Cdistance, 0*Deg2Rad );
        bondAtom( UnivalentAtom("HE2", Element::Hydrogen()), "CE2/bond3", C_Hdistance );

        bondAtom( TrivalentAtom("CD2", Element::Carbon()), "CE2/bond2", C_Cdistance, 0*Deg2Rad );
        bondAtom( UnivalentAtom("HD2", Element::Hydrogen()), "CD2/bond3", C_Hdistance );

        addRingClosingBond("CD2/bond2", "CG/bond3", 0.1459);

        defineDihedralAngle("chi2", "CD1", "CG", "CB", "CA");

        setDefaultDihedralAngle("chi2", 90*Deg2Rad); // from Dunbrack site

        setBondMobility(BondMobility::Rigid, "CG", "CD1");
        setBondMobility(BondMobility::Rigid, "CD1", "CE1");
        setBondMobility(BondMobility::Rigid, "CE1", "CZ");
        setBondMobility(BondMobility::Rigid, "CZ", "CE2");
        setBondMobility(BondMobility::Rigid, "CE2", "CD2");
        setBondMobility(BondMobility::Rigid, "CD2", "CG");
    }
};


class SimTK_MOLMODEL_EXPORT Protein : public Biopolymer {
public:

    void loadFromPdbStructure(const PdbStructure& pdbStructure, SimTK::Real targetRms)
    {
        const PdbChain& pdbChain = pdbStructure.getModel(Pdb::ModelIndex(0)).getChain(Pdb::ChainIndex(0));
        loadFromPdbChain(pdbChain, targetRms);
    }
    
    void loadFromPdbChain(const PdbChain& pdbChain, SimTK::Real targetRms);
    
    explicit Protein(std::istream& pdbStream, SimTK::Real targetRms = 0.02)
    {
        PdbStructure pdbStructure(pdbStream);
        loadFromPdbStructure(pdbStructure, targetRms);
    }
    
    explicit Protein(const PdbStructure& pdbStructure, SimTK::Real targetRms = 0.02)
    {
        loadFromPdbStructure(pdbStructure, targetRms);
    }
    
    explicit Protein(const PdbChain& pdbChain, SimTK::Real targetRms)
    {
        loadFromPdbChain(pdbChain, targetRms);
    }
    
    explicit Protein(const Sequence& seq, const BondMobility::Mobility mobility= BondMobility::Rigid, bool addEndCaps = true  ) 
    {
        initialize(seq, mobility, addEndCaps);
    }
    
    explicit Protein( const Sequence& seq, CompoundSystem& system, Transform transform = Transform() ) 
    {
        initialize(seq, BondMobility::Rigid);
        system.adoptCompound(*this, transform);
    }
    
protected:
    void initialize( const Sequence& seq, const BondMobility::Mobility mobility , bool addEndCaps = true) 
    {
        Biotype::initializePopularBiotypes();
        
        // TODO - provide alternatives for end caps
        String previousResidueName = "acetyl0";
        if (addEndCaps) appendResidue(previousResidueName, AcetylResidue());

        for (int resi = 0; resi < (int)seq.size(); ++resi)
        {
            AminoAcidResidue residue = AminoAcidResidue::create(seq[resi]);
            residue.setPdbResidueNumber(resi + 1);
            residue.assignBiotypes();

            // Name residue subcompound after its number in the sequence
            // name must be unique within the protein
            String residueName(resi);
            appendResidue( residueName, residue );

            // Rigidify peptide bond 
            if (getNumResidues() > 1) // First residue lacks a preceding peptide bond
            {
                String omegaAngleName = String("omega") + String(resi);

                defineDihedralAngle(
                    omegaAngleName, 
                    previousResidueName+"/C/bond1", 
                    residueName+"/N/bond2");
                setDefaultDihedralAngle(omegaAngleName, 180*Deg2Rad);

                // Make omega torsion rigid
                String cAtomName = previousResidueName + "/" + "C";
                String nAtomName = residueName + "/" + "N";
                setBondMobility(mobility, cAtomName, nAtomName);
            }

            previousResidueName = residueName;
        }

        // Create C-terminal end cap, if we are producing end caps
        if (addEndCaps) {
            String cCapName = String("nme") + String((int)seq.size());
            appendResidue(cCapName, NMethylAmideResidue());
            // Rigidify end cap's omega angle (regardless of setting above).
            String omegaAngleName = String("omega") + String(cCapName);
            defineDihedralAngle(
                omegaAngleName, 
                previousResidueName+"/C/bond1", 
                cCapName+"/N/bond2");
            setDefaultDihedralAngle(omegaAngleName, 180*Deg2Rad);

            // Make omega torsion rigid
            String cAtomName = previousResidueName + "/" + "C";
            String nAtomName = cCapName + "/" + "N";
            setBondMobility(BondMobility::Rigid, cAtomName, nAtomName);
        }
    }
};


} // namespace SimTK

#endif // SimTK_MOLMODEL_PROTEIN_H_
