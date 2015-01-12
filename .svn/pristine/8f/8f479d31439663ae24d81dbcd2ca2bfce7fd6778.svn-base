#ifndef SimTK_MOLMODEL_RNA_H_
#define SimTK_MOLMODEL_RNA_H_

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
 * Contributors: Samuel Flores                                                *
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
#include "molmodel/internal/NA.h"       

namespace SimTK {

// Nucleoside has sugar and base but not phosphate
class SimTK_MOLMODEL_EXPORT RibonucleosideResidue : public RiboseCore {
public:
    explicit RibonucleosideResidue(String name, String threeLetterCode = "Unk", char oneLetterCode = '?')
        : RiboseCore(name, threeLetterCode, oneLetterCode)
    {
        // 2' hydroxyl group present in RNA but not in DNA
        bondAtom(AliphaticHydrogen("H2*"), "C2*/bond3");
        //#bondAtom(AliphaticHydrogen("H2*2"), "C2*/bond4");

        bondAtom(BivalentAtom("O2*", Element::Oxygen(), 108.50*Deg2Rad), "C2*/bond4", 0.1413, -60*Deg2Rad);
        bondAtom(UnivalentAtom("2HO*", Element::Hydrogen()), "O2*/bond2", 0.0960);

        // alternate atom name
        nameAtom("H2'", "H2*");

        nameAtom("O2'", "O2*");
        nameAtom("HO2'", "2HO*");
    }
    
    // Create larger residue with phosphodiester prepended
    BiopolymerResidue withPhosphodiester() const 
    {
        const RibonucleosideResidue& residue = *this;
        
        NaPhosphodiesterLinkage po2(residue.getCompoundName(), residue.getPdbResidueName(), residue.getOneLetterCode());
        po2.bondCompound("nucleoside", residue, "bondNext");
        po2.inheritAtomNames("nucleoside");
        po2.inheritCompoundSynonyms(residue);
        po2.nameBondCenter("bondNext", "nucleoside/bondNext");

        po2.setPdbResidueNumber(residue.getPdbResidueNumber());
        po2.setPdbChainId(residue.getPdbChainId());

        po2.defineDihedralAngle("alpha", "bondPrevious", "O5*/bond2");
        po2.setDefaultDihedralAngle("alpha", 295.0*Deg2Rad); // A-RNA (Schneider et al 2004)

        return po2;
    }
};

// Nucleotide has sugar and base and phosphate
class SimTK_MOLMODEL_EXPORT RibonucleotideResidue : public RibonucleosideResidue {
public:
    // Factory method
    static RibonucleotideResidue create(char oneLetterCode);
    static RibonucleotideResidue create(const PdbResidue& pdbResidue, char chainId = ' ');

    explicit RibonucleotideResidue(String name, String threeLetterCode = "Unk", char oneLetterCode = '?')
        : RibonucleosideResidue(name, threeLetterCode, oneLetterCode)
    {
        defineDihedralAngle("beta", "O5*/bond1", "C5*/bond2");
        setDefaultDihedralAngle("beta", 173.0*Deg2Rad); //  A-RNA (Schneider et al 2004)

        defineDihedralAngle("gamma", "C5*/bond1", "C4*/bond3");
        setDefaultDihedralAngle("gamma", 54.0*Deg2Rad); //  A-RNA (Schneider et al 2004)

        defineDihedralAngle("delta", "C4*/bond1", "C3*/bond2");
        setDefaultDihedralAngle("delta", 80.0*Deg2Rad); //  A-RNA (Schneider et al 2004)

        defineDihedralAngle("epsilon", "C3*/bond1", "O3*/bond2");
        setDefaultDihedralAngle("epsilon", 210.0*Deg2Rad); //  A-RNA (Schneider et al 2004)

        nameBondCenter("bondBase", "C1*/bond3");
    }


    class Adenylate;
    class Cytidylate;
    class Guanylate;
    class Uridylate;
    
    // Modified residues from tRNA
    class TwoNMethylGuanylate;
};


// PurineBaseCore represents the atoms in common among purine bases such as Adenine and Guanine
class PurineBaseCore : public Compound {
public:
    PurineBaseCore()
    {
        setBaseAtom( TrivalentAtom("N9",  Element::Nitrogen(), 125.8*Deg2Rad, 128.8*Deg2Rad));

        bondAtom( TrivalentAtom("C8",  Element::Carbon(),   113.9*Deg2Rad, 123.05*Deg2Rad),  "N9/bond3", 0.13710, 180*Deg2Rad,   BondMobility::Rigid );

        bondAtom( UnivalentAtom("H8",  Element::Hydrogen()),                                 "C8/bond3",  0.10800 );
        bondAtom(  BivalentAtom("N7",  Element::Nitrogen(), 103.8*Deg2Rad),                  "C8/bond2",  0.13040, 0*Deg2Rad,    BondMobility::Rigid );
        bondAtom( TrivalentAtom("C5",  Element::Carbon(),   110.4*Deg2Rad, 132.40*Deg2Rad),  "N7/bond2",  0.13910, 0*Deg2Rad,    BondMobility::Rigid );
        
        bondAtom( TrivalentAtom("C4",  Element::Carbon(),   126.2*Deg2Rad, 106.2*Deg2Rad),   "N9/bond2",  0.13740, 0*Deg2Rad,    BondMobility::Rigid );

        addRingClosingBond("C5/bond2",                                                       "C4/bond3",  0.13700, 0*Deg2Rad,    BondMobility::Rigid );

        // 114.30 is average between amber angle for guanine(111.30) and adenine(117.30)
        bondAtom( TrivalentAtom("C6",  Element::Carbon(), 114.30*Deg2Rad, 122.85*Deg2Rad),   "C5/bond3",  0.14040, 180*Deg2Rad,  BondMobility::Rigid );

        bondAtom(  BivalentAtom("N3",  Element::Nitrogen(), 118.60*Deg2Rad),                 "C4/bond2",  0.13540, 180*Deg2Rad,  BondMobility::Rigid );
        bondAtom( TrivalentAtom("C2",  Element::Carbon()),                                   "N3/bond2",  0.13240, 0*Deg2Rad,    BondMobility::Rigid );

        // Two bond centers open, inboard (N9), and C6
        nameBondCenter("bondC6", "C6/bond3");
        nameBondCenter("bondC2", "C2/bond3");
    }
};


class AdenineBase : public PurineBaseCore
{
public:
    AdenineBase()
    {
        // Lone pair on N1
        bondAtom(  BivalentAtom("N1",  Element::Nitrogen(), 118.60*Deg2Rad), "C6/bond2",  0.13390, 0*Deg2Rad, BondMobility::Rigid );
        addRingClosingBond("N1/bond2",                                       "C2/bond2",  0.13240, 0*Deg2Rad, BondMobility::Rigid);

        // Hydrogen on C2
        bondAtom( UnivalentAtom("H2",  Element::Hydrogen()),                 "bondC2",    0.10800,            BondMobility::Rigid );

        // Amine on C6
        bondAtom( TrivalentAtom("N6",  Element::Nitrogen()),                 "bondC6",    0.13400, 0*Deg2Rad, BondMobility::Rigid );
        bondAtom( UnivalentAtom("H61", Element::Hydrogen()),                 "N6/bond2",  0.10100,            BondMobility::Rigid );
        bondAtom( UnivalentAtom("H62", Element::Hydrogen()),                 "N6/bond3",  0.10100,            BondMobility::Rigid );
    }
};

class GuanineBase : public PurineBaseCore
{
public:
    GuanineBase()
    {
        // Hydrogen on N1
        bondAtom( TrivalentAtom("N1",  Element::Nitrogen(), 125.20*Deg2Rad, 116.8*Deg2Rad), "C6/bond2",  0.13880, 0*Deg2Rad, BondMobility::Rigid );
        bondAtom(UnivalentAtom("H1", Element::Hydrogen()),   "N1/bond3", 0.10100, BondMobility::Rigid);
        addRingClosingBond("N1/bond2",                       "C2/bond2",  0.13240, 0*Deg2Rad, BondMobility::Rigid);

        // Amine on C2
        bondAtom( TrivalentAtom("N2",  Element::Nitrogen()), "bondC2",    0.13400, 0*Deg2Rad, BondMobility::Rigid );
        bondAtom( UnivalentAtom("H21", Element::Hydrogen()), "N2/bond2",  0.10100,            BondMobility::Rigid );
        bondAtom( UnivalentAtom("H22", Element::Hydrogen()), "N2/bond3",  0.10100,            BondMobility::Rigid );

        // Carbonyl on C6
        bondAtom( UnivalentAtom("O6",  Element::Oxygen()),   "bondC6",    0.12290,            BondMobility::Rigid );

    }
};

class TwoNMethylGuanineBaseGroup : public PurineBaseCore
{
public:
	TwoNMethylGuanineBaseGroup() {
        addCompoundSynonym("Guanosine"); // for resolving biotypes KLUDGE
		// TODO
	}
};


// PyrimidineBaseCore represents the atoms in common among pyrimidine bases such as Cytosine, Uracil, and Thymine
class PyrimidineBaseCore : public Compound {
public:
    PyrimidineBaseCore()
    {
    }
};


class CytosineBase : public PyrimidineBaseCore
{
public:
    CytosineBase()
    {
        setBaseAtom(   TrivalentAtom("N1",  Element::Nitrogen(), 125.8*Deg2Rad, 128.8*Deg2Rad));
	//scf added temporarily
        //bondAtom( UnivalentAtom("HN1", Element::Hydrogen()), "N1/bond4", 0.1010, 0*Deg2Rad, BondMobility::Rigid);

        bondAtom( TrivalentAtom("C2", Element::Carbon(), 118.6*Deg2Rad, 120.9*Deg2Rad),  "N1/bond2", 0.1383, 180*Deg2Rad,   BondMobility::Rigid );
        bondAtom( UnivalentAtom("O2", Element::Oxygen()), "C2/bond3", 0.1229, 0*Deg2Rad, BondMobility::Rigid);

        bondAtom( BivalentAtom("N3", Element::Nitrogen(), 120.5*Deg2Rad), "C2/bond2", 0.1358, 0*Deg2Rad, BondMobility::Rigid);

        bondAtom( TrivalentAtom("C4", Element::Carbon(), 121.5*Deg2Rad, 119.3*Deg2Rad), "N3/bond2", 0.1339, 0*Deg2Rad, BondMobility::Rigid);
        bondAtom( TrivalentAtom("N4", Element::Nitrogen(), 120.0*Deg2Rad, 120.0*Deg2Rad), "C4/bond3", 0.1340, 0*Deg2Rad, BondMobility::Rigid);
        bondAtom( UnivalentAtom("H41", Element::Hydrogen()), "N4/bond2", 0.1010, 0*Deg2Rad, BondMobility::Rigid);
        bondAtom( UnivalentAtom("H42", Element::Hydrogen()), "N4/bond3", 0.1010, 0*Deg2Rad, BondMobility::Rigid);

        bondAtom( TrivalentAtom("C5", Element::Carbon(), 117.0*Deg2Rad, 119.70*Deg2Rad), "C4/bond2", 0.1433, 0*Deg2Rad, BondMobility::Rigid);
        bondAtom( UnivalentAtom("H5", Element::Hydrogen()), "C5/bond3", 0.1080, 0*Deg2Rad, BondMobility::Rigid);

        bondAtom( TrivalentAtom("C6", Element::Carbon(), 121.20*Deg2Rad, 119.70*Deg2Rad), "C5/bond2", 0.1350, 0*Deg2Rad, BondMobility::Rigid);
        bondAtom( UnivalentAtom("H6", Element::Hydrogen()), "C6/bond3", 0.1080, 0*Deg2Rad, BondMobility::Rigid);

        addRingClosingBond("C6/bond2", "N1/bond3", 0.1365, 0*Deg2Rad, BondMobility::Rigid);
    }
};

class UracilBase : public PyrimidineBaseCore
{
public:
    UracilBase()
    {
        setBaseAtom( TrivalentAtom("N1",  Element::Nitrogen(), 125.8*Deg2Rad, 128.8*Deg2Rad));

        bondAtom( TrivalentAtom("C2", Element::Carbon(), 118.6*Deg2Rad, 120.9*Deg2Rad),  "N1/bond2", 0.1383, 180*Deg2Rad,   BondMobility::Rigid );
        bondAtom( UnivalentAtom("O2", Element::Oxygen()), "C2/bond3", 0.1229, 0*Deg2Rad, BondMobility::Rigid);

        bondAtom(TrivalentAtom("N3", Element::Nitrogen(), 126.4*Deg2Rad, 116.8*Deg2Rad), "C2/bond2", 0.1358, 0*Deg2Rad, BondMobility::Rigid);
        bondAtom(UnivalentAtom("H3", Element::Hydrogen()), "N3/bond3", 0.1010, 0*Deg2Rad, BondMobility::Rigid);

        bondAtom( TrivalentAtom("C4", Element::Carbon(), 114.0*Deg2Rad, 120.6*Deg2Rad),  "N3/bond2", 0.1388, 0*Deg2Rad,   BondMobility::Rigid );
        bondAtom( UnivalentAtom("O4", Element::Oxygen()), "C4/bond3", 0.1229, 0*Deg2Rad, BondMobility::Rigid);

        bondAtom( TrivalentAtom("C5", Element::Carbon(), 120.7*Deg2Rad, 119.70*Deg2Rad), "C4/bond2", 0.1444, 0*Deg2Rad, BondMobility::Rigid);
        bondAtom( UnivalentAtom("H5", Element::Hydrogen()), "C5/bond3", 0.1080, 0*Deg2Rad, BondMobility::Rigid);

        bondAtom( TrivalentAtom("C6", Element::Carbon(), 121.20*Deg2Rad, 119.70*Deg2Rad), "C5/bond2", 0.1350, 0*Deg2Rad, BondMobility::Rigid);
        bondAtom( UnivalentAtom("H6", Element::Hydrogen()), "C6/bond3", 0.1080, 0*Deg2Rad, BondMobility::Rigid);

        addRingClosingBond("C6/bond2", "N1/bond3", 0.1365, 0*Deg2Rad, BondMobility::Rigid);
    }
};

// TODO - these should not be classes, e.g. Adenylate should be an instance of RibonucleotideResidue
class  RibonucleotideResidue::Adenylate : public RibonucleotideResidue {
public:
    Adenylate() : RibonucleotideResidue("adenylate", "A  ", 'A') 
    {
        bondCompound("base", AdenineBase(), "bondBase", 0.14750, 200*Deg2Rad, BondMobility::Torsion);
        inheritAtomNames("base");

        defineDihedralAngle("chi", "O4*", "C1*", "N9", "C4");
        setDefaultDihedralAngle("chi", 199.0*Deg2Rad);

        addCompoundSynonym("Adenosine"); // for resolving biotypes
    }
};

class  RibonucleotideResidue::Guanylate : public RibonucleotideResidue {
public:
    Guanylate() : RibonucleotideResidue("guanylate", "G  ", 'G') 
    {
        bondCompound("base", GuanineBase(), "bondBase", 0.14710, 200*Deg2Rad, BondMobility::Torsion);
        inheritAtomNames("base");

        defineDihedralAngle("chi", "O4*", "C1*", "N9", "C4");
        setDefaultDihedralAngle("chi", 199.0*Deg2Rad);

        addCompoundSynonym("Guanosine"); // for resolving biotypes
    }
};

class  RibonucleotideResidue::Cytidylate : public RibonucleotideResidue {
public:
    Cytidylate() : RibonucleotideResidue("cytidylate", "C  ", 'C') 
    {
        bondCompound("base", CytosineBase(), "bondBase", 0.14710, 200*Deg2Rad, BondMobility::Torsion);
        inheritAtomNames("base");

        defineDihedralAngle("chi", "O4*", "C1*", "N1", "C2");
        setDefaultDihedralAngle("chi", 199.0*Deg2Rad);

        addCompoundSynonym("Cytidine"); // for resolving biotypes
    }
};

class  RibonucleotideResidue::Uridylate : public RibonucleotideResidue {
public:
    Uridylate() : RibonucleotideResidue("uridylate", "U  ", 'U') 
    {
        bondCompound("base", UracilBase(), "bondBase", 0.14710, 200*Deg2Rad, BondMobility::Torsion);
        inheritAtomNames("base");

        defineDihedralAngle("chi", "O4*", "C1*", "N1", "C2");
        setDefaultDihedralAngle("chi", 199.0*Deg2Rad);

        addCompoundSynonym("Uridine"); // for resolving biotypes
    }
};


class TwoNMethylGuanidineGroup : public Compound {
public:
    TwoNMethylGuanidineGroup()
    {
        instantiateBiotypes();

        setPdbResidueName("2MG");

        setCompoundName("TwoNMethylGuanidineGroup");

        setBaseAtom( TrivalentAtom("N9", Element::Nitrogen(), 127.767*Deg2Rad, 112.885*Deg2Rad) );

        bondAtom( TrivalentAtom("C8", Element::Carbon(), 113.564*Deg2Rad, 120.934*Deg2Rad), "N9/bond2", 0.13740, 180.00*Deg2Rad, BondMobility::Rigid );
        bondAtom( BivalentAtom("N7", Element::Nitrogen(), 104.422*Deg2Rad), "C8/bond2", 0.12800, 0.02*Deg2Rad, BondMobility::Rigid );
        bondAtom( TrivalentAtom("C5", Element::Carbon(), 130.918*Deg2Rad, 112.885*Deg2Rad), "N7/bond2", 0.13760, 179.54*Deg2Rad, BondMobility::Rigid );
        bondAtom( TrivalentAtom("C6", Element::Carbon(), 131.213*Deg2Rad, 109.672*Deg2Rad), "C5/bond2", 0.14350, -0.17*Deg2Rad, BondMobility::Rigid );
        bondAtom( UnivalentAtom("O6", Element::Oxygen()), "C6/bond2", 0.11940 );
        bondAtom( TrivalentAtom("N1", Element::Nitrogen(), 126.495*Deg2Rad, 113.649*Deg2Rad), "C6/bond3", 0.14150, 0.32*Deg2Rad, BondMobility::Rigid );
        bondAtom( TrivalentAtom("C2", Element::Carbon(), 116.091*Deg2Rad, 123.417*Deg2Rad), "N1/bond2", 0.13620, 178.30*Deg2Rad, BondMobility::Rigid );
        bondAtom( TrivalentAtom("N2", Element::Nitrogen(), 116.703*Deg2Rad, 120.926*Deg2Rad), "C2/bond2", 0.13550, 20.71*Deg2Rad, BondMobility::Rigid );
        bondAtom( UnivalentAtom("1H2", Element::Hydrogen()), "N2/bond2", 0.09940 );
        bondAtom( QuadrivalentAtom("C10", Element::Carbon()), "N2/bond3", 0.14490, 176.13*Deg2Rad, BondMobility::Torsion );
        bondAtom( UnivalentAtom("H20", Element::Hydrogen()), "C10/bond2", 0.10810 );
        bondAtom( UnivalentAtom("H21", Element::Hydrogen()), "C10/bond3", 0.10850 );
        bondAtom( UnivalentAtom("H22", Element::Hydrogen()), "C10/bond4", 0.10790 );
        bondAtom( BivalentAtom("N3", Element::Nitrogen(), 112.885*Deg2Rad), "C2/bond3", 0.12910, 0.73*Deg2Rad, BondMobility::Rigid );
        bondAtom( TrivalentAtom("C4", Element::Carbon(), 109.368*Deg2Rad, 104.422*Deg2Rad), "N3/bond2", 0.13560, -119.37*Deg2Rad, BondMobility::Rigid );
        bondAtom( UnivalentAtom("H1", Element::Hydrogen()), "N1/bond3", 0.09980 );
        bondAtom( UnivalentAtom("H8", Element::Hydrogen()), "C8/bond3", 0.10710 );

        addRingClosingBond( "C4/bond2", "N9/bond3", 0.15);
        addRingClosingBond( "C4/bond3", "C5/bond3", 0.15);

        setBiotypeIndex( "N9", Biotype::get("TwoNMethylGuanidineGroup", "N9").getIndex());
        setBiotypeIndex( "C8", Biotype::get("TwoNMethylGuanidineGroup", "C8").getIndex());
        setBiotypeIndex( "N7", Biotype::get("TwoNMethylGuanidineGroup", "N7").getIndex());
        setBiotypeIndex( "C5", Biotype::get("TwoNMethylGuanidineGroup", "C5").getIndex());
        setBiotypeIndex( "C6", Biotype::get("TwoNMethylGuanidineGroup", "C6").getIndex());
        setBiotypeIndex( "O6", Biotype::get("TwoNMethylGuanidineGroup", "O6").getIndex());
        setBiotypeIndex( "N1", Biotype::get("TwoNMethylGuanidineGroup", "N1").getIndex());
        setBiotypeIndex( "C2", Biotype::get("TwoNMethylGuanidineGroup", "C2").getIndex());
        setBiotypeIndex( "N2", Biotype::get("TwoNMethylGuanidineGroup", "N2").getIndex());
        setBiotypeIndex( "1H2", Biotype::get("TwoNMethylGuanidineGroup", "1H2").getIndex());
        setBiotypeIndex( "C10", Biotype::get("TwoNMethylGuanidineGroup", "C10").getIndex());
        setBiotypeIndex( "H20", Biotype::get("TwoNMethylGuanidineGroup", "H20").getIndex());
        setBiotypeIndex( "H21", Biotype::get("TwoNMethylGuanidineGroup", "H21").getIndex());
        setBiotypeIndex( "H22", Biotype::get("TwoNMethylGuanidineGroup", "H22").getIndex());
        setBiotypeIndex( "N3", Biotype::get("TwoNMethylGuanidineGroup", "N3").getIndex());
        setBiotypeIndex( "C4", Biotype::get("TwoNMethylGuanidineGroup", "C4").getIndex());
        setBiotypeIndex( "H1", Biotype::get("TwoNMethylGuanidineGroup", "H1").getIndex());
        setBiotypeIndex( "H8", Biotype::get("TwoNMethylGuanidineGroup", "H8").getIndex());

        setDefaultDihedralAngle(0.020*Deg2Rad, "C5", "N7", "C8", "N9");
        setDefaultDihedralAngle(179.544*Deg2Rad, "C6", "C5", "N7", "C8");
        setDefaultDihedralAngle(-0.165*Deg2Rad, "O6", "C6", "C5", "N7");
        setDefaultDihedralAngle(-179.886*Deg2Rad, "N1", "C6", "C5", "N7");
        setDefaultDihedralAngle(0.318*Deg2Rad, "C2", "N1", "C6", "C5");
        setDefaultDihedralAngle(178.303*Deg2Rad, "N2", "C2", "N1", "C6");
        setDefaultDihedralAngle(20.715*Deg2Rad, "1H2", "N2", "C2", "N1");
        setDefaultDihedralAngle(173.145*Deg2Rad, "C10", "N2", "C2", "N1");
        setDefaultDihedralAngle(176.133*Deg2Rad, "H20", "C10", "N2", "C2");
        setDefaultDihedralAngle(-63.298*Deg2Rad, "H21", "C10", "N2", "C2");
        setDefaultDihedralAngle(57.486*Deg2Rad, "H22", "C10", "N2", "C2");
        setDefaultDihedralAngle(-0.550*Deg2Rad, "N3", "C2", "N1", "C6");
        setDefaultDihedralAngle(0.727*Deg2Rad, "C4", "N3", "C2", "N1");
        setDefaultDihedralAngle(176.522*Deg2Rad, "H1", "N1", "C6", "C5");
    } // end constructor TwoNMethylGuanidineGroup 

    static void instantiateBiotypes() {
        if (! Biotype::exists("TwoNMethylGuanidineGroup", "N9") )
            Biotype::defineBiotype(Element::Nitrogen(), 3, "TwoNMethylGuanidineGroup", "N9"); 
        if (! Biotype::exists("TwoNMethylGuanidineGroup", "C8") )
            Biotype::defineBiotype(Element::Carbon(), 3, "TwoNMethylGuanidineGroup", "C8"); 
        if (! Biotype::exists("TwoNMethylGuanidineGroup", "N7") )
            Biotype::defineBiotype(Element::Nitrogen(), 2, "TwoNMethylGuanidineGroup", "N7"); 
        if (! Biotype::exists("TwoNMethylGuanidineGroup", "C5") )
            Biotype::defineBiotype(Element::Carbon(), 3, "TwoNMethylGuanidineGroup", "C5"); 
        if (! Biotype::exists("TwoNMethylGuanidineGroup", "C6") )
            Biotype::defineBiotype(Element::Carbon(), 3, "TwoNMethylGuanidineGroup", "C6"); 
        if (! Biotype::exists("TwoNMethylGuanidineGroup", "O6") )
            Biotype::defineBiotype(Element::Oxygen(), 1, "TwoNMethylGuanidineGroup", "O6"); 
        if (! Biotype::exists("TwoNMethylGuanidineGroup", "N1") )
            Biotype::defineBiotype(Element::Nitrogen(), 3, "TwoNMethylGuanidineGroup", "N1"); 
        if (! Biotype::exists("TwoNMethylGuanidineGroup", "C2") )
            Biotype::defineBiotype(Element::Carbon(), 3, "TwoNMethylGuanidineGroup", "C2"); 
        if (! Biotype::exists("TwoNMethylGuanidineGroup", "N2") )
            Biotype::defineBiotype(Element::Nitrogen(), 3, "TwoNMethylGuanidineGroup", "N2"); 
        if (! Biotype::exists("TwoNMethylGuanidineGroup", "1H2") )
            Biotype::defineBiotype(Element::Hydrogen(), 1, "TwoNMethylGuanidineGroup", "1H2"); 
        if (! Biotype::exists("TwoNMethylGuanidineGroup", "C10") )
            Biotype::defineBiotype(Element::Carbon(), 4, "TwoNMethylGuanidineGroup", "C10"); 
        if (! Biotype::exists("TwoNMethylGuanidineGroup", "H20") )
            Biotype::defineBiotype(Element::Hydrogen(), 1, "TwoNMethylGuanidineGroup", "H20"); 
        if (! Biotype::exists("TwoNMethylGuanidineGroup", "H21") )
            Biotype::defineBiotype(Element::Hydrogen(), 1, "TwoNMethylGuanidineGroup", "H21"); 
        if (! Biotype::exists("TwoNMethylGuanidineGroup", "H22") )
            Biotype::defineBiotype(Element::Hydrogen(), 1, "TwoNMethylGuanidineGroup", "H22"); 
        if (! Biotype::exists("TwoNMethylGuanidineGroup", "N3") )
            Biotype::defineBiotype(Element::Nitrogen(), 2, "TwoNMethylGuanidineGroup", "N3"); 
        if (! Biotype::exists("TwoNMethylGuanidineGroup", "C4") )
            Biotype::defineBiotype(Element::Carbon(), 3, "TwoNMethylGuanidineGroup", "C4"); 
        if (! Biotype::exists("TwoNMethylGuanidineGroup", "H1") )
            Biotype::defineBiotype(Element::Hydrogen(), 1, "TwoNMethylGuanidineGroup", "H1"); 
        if (! Biotype::exists("TwoNMethylGuanidineGroup", "H8") )
            Biotype::defineBiotype(Element::Hydrogen(), 1, "TwoNMethylGuanidineGroup", "H8"); 
    } // end instantiateBiotypes

}; // end class TwoNMethylGuanidineGroup 


class  RibonucleotideResidue::TwoNMethylGuanylate : public RibonucleotideResidue {
public:
	TwoNMethylGuanylate() : RibonucleotideResidue("2-N-methyl-guanylate", "2MG", 'X') 
    {
        bondCompound("base", TwoNMethylGuanidineGroup(), "bondBase", 0.14750, 200*Deg2Rad, BondMobility::Torsion);
        inheritAtomNames("base");

        defineDihedralAngle("chi", "O4*", "C1*", "N9", "C4");
        setDefaultDihedralAngle("chi", 210*Deg2Rad);

        addCompoundSynonym("TwoNMethylGuanidineGroup"); // for resolving biotypes
        
        setBiotypeIndex( "N9", Biotype::get("TwoNMethylGuanidineGroup", "N9").getIndex());
        setBiotypeIndex( "C8", Biotype::get("TwoNMethylGuanidineGroup", "C8").getIndex());
        setBiotypeIndex( "N7", Biotype::get("TwoNMethylGuanidineGroup", "N7").getIndex());
        setBiotypeIndex( "C5", Biotype::get("TwoNMethylGuanidineGroup", "C5").getIndex());
        setBiotypeIndex( "C6", Biotype::get("TwoNMethylGuanidineGroup", "C6").getIndex());
        setBiotypeIndex( "O6", Biotype::get("TwoNMethylGuanidineGroup", "O6").getIndex());
        setBiotypeIndex( "N1", Biotype::get("TwoNMethylGuanidineGroup", "N1").getIndex());
        setBiotypeIndex( "C2", Biotype::get("TwoNMethylGuanidineGroup", "C2").getIndex());
        setBiotypeIndex( "N2", Biotype::get("TwoNMethylGuanidineGroup", "N2").getIndex());
        setBiotypeIndex( "1H2", Biotype::get("TwoNMethylGuanidineGroup", "1H2").getIndex());
        setBiotypeIndex( "C10", Biotype::get("TwoNMethylGuanidineGroup", "C10").getIndex());
        setBiotypeIndex( "H20", Biotype::get("TwoNMethylGuanidineGroup", "H20").getIndex());
        setBiotypeIndex( "H21", Biotype::get("TwoNMethylGuanidineGroup", "H21").getIndex());
        setBiotypeIndex( "H22", Biotype::get("TwoNMethylGuanidineGroup", "H22").getIndex());
        setBiotypeIndex( "N3", Biotype::get("TwoNMethylGuanidineGroup", "N3").getIndex());
        setBiotypeIndex( "C4", Biotype::get("TwoNMethylGuanidineGroup", "C4").getIndex());
        setBiotypeIndex( "H1", Biotype::get("TwoNMethylGuanidineGroup", "H1").getIndex());
        setBiotypeIndex( "H8", Biotype::get("TwoNMethylGuanidineGroup", "H8").getIndex());
        
    }
};






class  RNA : public Biopolymer {
public:

	explicit RNA(const PdbStructure& pdbStructure) 
	{
		// Assume first chain of first model is wanted
		const PdbChain& pdbChain = pdbStructure.getModel(Pdb::ModelIndex(0)).getChain(Pdb::ChainIndex(0));
		initializeFromPdbChain(pdbChain);
	}

	explicit RNA(const PdbChain& pdbChain)
	{
		initializeFromPdbChain(pdbChain);
	}
    /**
     * /brief First parameter is just the RNA sequence in single-letter code.  AUGC are the only residue types supported at this time.  Second parameter, when set to zero, results in an RNA with no capping hydroxyls at the termini.  Instead the terminal residues have the same complement of atoms as residues in the interior of the chain.
     *
     */
    explicit RNA(const Sequence& seq, bool useCappingHydroxyls = true) 
    {
        String previousResidueName;
        Biotype::initializePopularBiotypes();

        // TODO - create 5' end cap
        for (int resi = 0; resi < (int)seq.size(); ++resi)
        {
        	
            RibonucleotideResidue residue = RibonucleotideResidue::create(seq[resi]);
			residue.assignBiotypes();

            residue.setPdbResidueNumber(resi + 1);

            // Name residue subcompound after its number in the sequence
            // name must be unique within the protein
            String residueName(resi);

            // Cap the 3' end with a hydroxyl group
            if ((resi == seq.size() - 1) && (useCappingHydroxyls)) {
                    residue.bondCompound(
                            "3PrimeHydroxyl", 
                            ThreePrimeNaHydroxylGroup(residueName, residue.getPdbResidueName(), 'X'), 
                            "bondNext", 
                                    0.0960);
            residue.inheritAtomNames("3PrimeHydroxyl");
                    // Oxygen biotype must be changed for the amber atom types to come out correctly
            if ( (residue.hasAtom("O3*")) && ( Biotype::exists("Hydroxyl, RNA", "O3*", SimTK::Ordinality::Final)) )
                    residue.setBiotypeIndex( "O3*", Biotype::get("Hydroxyl, RNA", "O3*", SimTK::Ordinality::Final).getIndex() );
            }

 if (resi == 0) {
                    // first residue needs end cap
                    residue.convertInboardBondCenterToOutboard();

                            if (useCappingHydroxyls) { // 5' hydroxyl
                            residue.bondCompound(
                                            "5PrimeHydroxyl", 
                                            FivePrimeNaHydroxylGroup(residueName, residue.getPdbResidueName(), 'X'), 
                                            "bondPrevious", 
                                                    0.0960);
                            residue.inheritAtomNames("5PrimeHydroxyl");
                            if ( (residue.hasAtom("O5*")) && ( Biotype::exists("Hydroxyl, RNA", "O5*", SimTK::Ordinality::Initial)) )
                                    residue.setBiotypeIndex( "O5*", Biotype::get("Hydroxyl, RNA", "O5*", SimTK::Ordinality::Initial).getIndex() );
                            }

                            else { // 5' phosphate
                            residue.bondCompound(
                                            "5PrimePhosphate", 
                                            FivePrimeNaPhosphateGroup(residueName, residue.getPdbResidueName(), 'X'), 
                                            "bondPrevious");
                            residue.inheritAtomNames("5PrimePhosphate");
                            }

                    appendResidue( residueName, residue );
            }
            else {
                // non-first residue needs phosphodiester linkage
                appendResidue( residueName, residue.withPhosphodiester() );

                // Define zeta angle
                String zetaName = String("zeta") + String(resi - 1);
                defineDihedralAngle(zetaName, previousResidueName + "/O3*/bond1", residueName + "/P/bond2");
                setDefaultDihedralAngle(zetaName, 287.0*Deg2Rad); // (Schneider et al 2004)
            }

            previousResidueName = residueName;
        }
    }
 /**
  *  /brief This method sets the BondMobility::Mobility for each residue in a stretch of polynucleotide spanning residues startResidue to endResidue.  It also sets the BondMobility::Mobility for bonds connecting the included residues to the same value.
  *
  * Added by Samuel Flores
  */


/**
 * /brief This method sets the BondMobility for all bonds in a certain stretch of residues of a polynucleotide chain
 *
 * added by scf 
 */
    RNA& setRNABondMobility (BondMobility::Mobility  mobility, ResidueInfo::Index startResidue, ResidueInfo::Index endResidue){
        for (ResidueInfo::Index q =startResidue; q<=endResidue; q++) {
            setResidueBondMobility(q, mobility);
            // (updResidue(ResidueInfo::Index(q))).setCompoundBondMobility(mobility);
            if (q>startResidue) {
                std::stringstream ss1;
                ss1<<q-1<<"/O3*";
                std::stringstream ss2;
                ss2<<q<<"/P";
                setBondMobility(mobility ,ss1.str() ,ss2.str()  ); 
            }
        }
	    return *this;
    }
/**
 * /brief This method returns the base normal of a given residue.
 * the direction of the normal follows the "curl" of the atoms in the cycle of which the glycosidic nitrogen is a part, following the numbering convention in use here.
 * added by scf 
 *
 * /param residueNumber starts at 0. bodyOrGroundFrame = 0 if the normal is desired in the frame of the glycosidic nitrogen, = 1 if in the frame of Ground.
 *
 */

    Vec3 calcRNABaseNormal  (const State & state,int residueNumber, int bodyOrGroundFrame, const SimbodyMatterSubsystem & matter) const {
	    const ResidueInfo& myResidue = getResidue(ResidueInfo::Index(residueNumber));
            String myResidueName = myResidue.getPdbResidueName();
	    //cout<<"[RNA.h:calcRNABaseNormal] myResidueName ="<<myResidueName<<endl;
            //String myResidueName = (updResidue(Compound::Index(residueNumber))).getPdbResidueName();
	    assert((myResidueName.compare("A  ") == 0) || (myResidueName.compare("G  ") == 0)  || (myResidueName.compare("U  ") == 0) || (myResidueName.compare("C  ") == 0));
            String atomName1;
            String atomName1a;
            String atomName1b;
	    if ((myResidueName.compare("A  ") == 0) || (myResidueName.compare("G  ") == 0)) {
	        atomName1  = "N9";
	        atomName1a = "C4";
	        atomName1b = "C8";
	    } else {
	        atomName1  = "N1";
	        atomName1a = "C2";
	        atomName1b = "C6";
	    }
            Compound::AtomIndex myAtomIndex1 = myResidue.getAtomIndex(atomName1 );
            Compound::AtomIndex myAtomIndex1a= myResidue.getAtomIndex(atomName1a);
            Compound::AtomIndex myAtomIndex1b= myResidue.getAtomIndex(atomName1b);
	    MobilizedBodyIndex myMobilizedBodyIndex1  = getAtomMobilizedBodyIndex(myAtomIndex1 );
	    MobilizedBodyIndex myMobilizedBodyIndex1a = getAtomMobilizedBodyIndex(myAtomIndex1a);
	    MobilizedBodyIndex myMobilizedBodyIndex1b = getAtomMobilizedBodyIndex(myAtomIndex1b);
            MobilizedBody body1 = matter.getMobilizedBody(myMobilizedBodyIndex1 );
            MobilizedBody body1a= matter.getMobilizedBody(myMobilizedBodyIndex1a);
            MobilizedBody body1b= matter.getMobilizedBody(myMobilizedBodyIndex1b);
	    if (bodyOrGroundFrame == 0 ) {
                Vec3 normalInBodyFrame = (
                    (~(body1.getBodyTransform(state)) * body1a.getBodyTransform(state) * getAtomLocationInMobilizedBodyFrame(myAtomIndex1a)  
                    - getAtomLocationInMobilizedBodyFrame(myAtomIndex1)) %
                    (~(body1.getBodyTransform(state)) * body1b.getBodyTransform(state) * getAtomLocationInMobilizedBodyFrame(myAtomIndex1b)
                    - getAtomLocationInMobilizedBodyFrame(myAtomIndex1)));
		normalInBodyFrame /= normalInBodyFrame.norm();
	        return normalInBodyFrame;
	    } else if (bodyOrGroundFrame == 1) {
                Vec3 normalInGroundFrame = 
                    (    
                    body1a.getBodyTransform(state) * getAtomLocationInMobilizedBodyFrame(myAtomIndex1a)
                    - body1.getBodyTransform(state) * getAtomLocationInMobilizedBodyFrame(myAtomIndex1)
                    ) % (
                    body1b.getBodyTransform(state) * getAtomLocationInMobilizedBodyFrame(myAtomIndex1b)
                    - body1.getBodyTransform(state) * getAtomLocationInMobilizedBodyFrame(myAtomIndex1)
                    );
                normalInGroundFrame /= normalInGroundFrame.norm();
                return normalInGroundFrame;
	    }

    }

private:
    void initializeFromPdbChain(
            const PdbChain& pdbChain, 
            Compound::MatchStratagem matchStratagem = 
                    Compound::Match_TopologyOnly)
    {
        setPdbChainId( pdbChain.getChainId() );
        
        String previousResidueName("");
        
        // TODO end caps
        for (Pdb::ResidueIndex r(0); r < (Pdb::ResidueIndex)pdbChain.getNumResidues(); ++r) 
        {
            const PdbResidue& pdbResidue = pdbChain.getResidue(r);
            String residueName = pdbResidue.getName() + String(pdbResidue.getPdbResidueNumber());
            if (pdbResidue.getInsertionCode() != ' ') // include unusual insertion codes in residue name
                    residueName += pdbResidue.getInsertionCode();
            
            RibonucleotideResidue residue = RibonucleotideResidue::create(pdbResidue, pdbChain.getChainId());
            residue.assignBiotypes();
            
            if (r == 0) { // first residue does not get phosphate
                appendResidue( residueName, residue );
            }
            else {
                // non-first residue needs phosphodiester linkage
                BiopolymerResidue resWithP = residue.withPhosphodiester();
                resWithP.assignBiotypes();
                appendResidue( residueName, resWithP );

                // Define zeta angle
                String zetaName = String("zeta") + String(residue.getPdbResidueNumber());
                defineDihedralAngle(zetaName, previousResidueName + "/O3*/bond1", residueName + "/P/bond2");
                setDefaultDihedralAngle(zetaName, 287.0*Deg2Rad); // (Schneider et al 2004)
            } 
            
            previousResidueName = residueName;            
        }
        
        if (matchStratagem != Compound::Match_TopologyOnly) 
        {
            Compound::AtomTargetLocations atomTargets = 
                    createAtomTargets(pdbChain);
            
            matchDefaultConfiguration(atomTargets, matchStratagem);
        }        
    }
};


} // namespace SimTK

#endif // SimTK_MOLMODEL_RNA_H_
