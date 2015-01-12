#ifndef SimTK_MOLMODEL_NA_H_
#define SimTK_MOLMODEL_NA_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Molmodel                            *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2007 Stanford University and the Authors.           *
 * Adapted from RNA.h (original by Michael Sherman, Christopher Bruns)        *
 * by  Samuel Flores                                                          *
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

namespace SimTK {

// RiboseCore leaves 4 bond centers open:
// 1) inboard center at O5* atom, for binding 5' phosphate
// 2) 
class SimTK_MOLMODEL_EXPORT RiboseCore : public BiopolymerResidue {
public:
    RiboseCore(const String& name, const String& tlc, char olc) 
        : BiopolymerResidue(name, tlc, olc)
    {
        setBaseAtom(BivalentAtom("O5*", Element::Oxygen(), 120.90*Deg2Rad));

        bondAtom(AliphaticCarbon("C5*"), "O5*/bond2", 0.1423, 180*Deg2Rad);
        bondAtom(AliphaticHydrogen("H5*1"), "C5*/bond4");
        bondAtom(AliphaticHydrogen("H5*2"), "C5*/bond3");

        bondAtom(AliphaticCarbon("C4*"), "C5*/bond2", 0.1510); // length from (Gelbin et al 1996)
        bondAtom(AliphaticHydrogen("H4*"), "C4*/bond4");

        bondAtom(BivalentAtom("O4*", Element::Oxygen(), 112.0*Deg2Rad), "C4*/bond2", 0.1453);

        bondAtom(AliphaticCarbon("C3*"), "C4*/bond3", 0.1524);
        bondAtom(AliphaticHydrogen("H3*"), "C3*/bond4");

        bondAtom(BivalentAtom("O3*", Element::Oxygen(), 119.7*Deg2Rad), "C3*/bond2", 0.1423);
 
        bondAtom(AliphaticCarbon("C2*"), "C3*/bond3", 0.1525);
        //bondAtom(AliphaticHydrogen("H2*"), "C2*/bond3");
        //bondAtom(AliphaticHydrogen("H2*2"), "C2*/bond4");

        bondAtom(AliphaticCarbon("C1*"), "C2*/bond2", 0.1528);
        bondAtom(AliphaticHydrogen("H1*"), "C1*/bond4");

        addRingClosingBond("C1*/bond2", "O4*/bond2", 0.1414);

        // Ring strain reduces ring bond angles
        setDefaultBondAngle(109.6*Deg2Rad, "C1*", "O4*", "C4*"); 
        setDefaultBondAngle(105.5*Deg2Rad, "O4*", "C4*", "C3*"); 
        setDefaultBondAngle(102.7*Deg2Rad, "C4*", "C3*", "C2*"); 
        setDefaultBondAngle(101.5*Deg2Rad, "C3*", "C2*", "C1*"); 
        setDefaultBondAngle(106.4*Deg2Rad, "C2*", "C1*", "O4*");

        nameBondCenter("bondO5", "O5*/bond1");
        nameBondCenter("bondO3", "O3*/bond2");
        nameBondCenter("bondC1", "C1*/bond3");

        // NA backbone dihedral angles
        // defineDihedralAngle("alpha", "bondO5", ""); // TODO requires phosphate
        //defineDihedralAngle("beta", "bondO5", "C5*/bond2");
        //defineDihedralAngle("gamma", "O5*", "C5*", "C4*", "C3*");
        //defineDihedralAngle("delta", "C5*", "C4*", "C3*", "C2*");
        //defineDihedralAngle("epsilon", "bondO3", "C3*/bond1");
        // TODO zeta requires next phosphate

        // ribose ring dihedral angles
        defineDihedralAngle( "nu0", "C4*", "O4*", "C1*", "C2*" );
        defineDihedralAngle( "nu1", "O4*", "C1*", "C2*", "C3*" );
        defineDihedralAngle( "nu2", "C1*", "C2*", "C3*", "C4*" );
        defineDihedralAngle( "nu3", "C2*", "C3*", "C4*", "O4*" );
        defineDihedralAngle( "nu4", "C3*", "C4*", "O4*", "C1*" );

        // Guestimated from Figure 2 of 
        // Schneider, B; Moravek, Z; and Berman, H.M. (2004)
        // "NA conformational classes"
        // Nucleic Acids Research 32(5): 1666-1677
        // setDefaultDihedralAngle("alpha", -60*Deg2Rad);
        //setDefaultDihedralAngle("beta", 180*Deg2Rad);
        //setDefaultDihedralAngle("gamma", 60*Deg2Rad);
        //setDefaultDihedralAngle("delta", 80*Deg2Rad);
        //setDefaultDihedralAngle("epsilon", 190*Deg2Rad);
        // setDefaultDihedralAngle("zeta", -70*Deg2Rad);

        // nu angles selected from resdue U52 of 1EHZ
        setDefaultDihedralAngle("nu0", 2.2*Deg2Rad);
        setDefaultDihedralAngle("nu1", -24.5*Deg2Rad);
        setDefaultDihedralAngle("nu2", 35.8*Deg2Rad);
        setDefaultDihedralAngle("nu3", -35.9*Deg2Rad); // same bond as delta, different reference atoms
        setDefaultDihedralAngle("nu4", 21.6*Deg2Rad);

        nameBondCenter("bondNext", "bondO3");
        nameBondCenter("bondPrevious", "bondO5");

        setDefaultInboardBondLength(0.16100);

        // alternate atom names
        nameAtom("O5'", "O5*");
        nameAtom("C5'", "C5*");
        nameAtom("H5'", "H5*1");
        nameAtom("H5''", "H5*2");
        nameAtom("C4'", "C4*");
        nameAtom("H4'", "H4*");
        nameAtom("O4'", "O4*");
        nameAtom("C3'", "C3*");
        nameAtom("H3'", "H3*");
        nameAtom("O3'", "O3*");
        nameAtom("C2'", "C2*");
        //nameAtom("H2'", "H2*");
        //nameAtom("H2'1", "H2*");
        //nameAtom("H2*1", "H2*");
        nameAtom("C1'", "C1*");
        nameAtom("H1'", "H1*");
    }

};


/// Phosphate group at 5' end (beginning) of NA
class SimTK_MOLMODEL_EXPORT FivePrimeNaHydroxylGroup : public BiopolymerResidue
{
public:
	explicit FivePrimeNaHydroxylGroup(String name, String threeLetterCode = "Unk", char oneLetterCode = '?') 
    : BiopolymerResidue(name, threeLetterCode, oneLetterCode)
    {
		instantiateBiotypes();
		
		setBaseAtom( UnivalentAtom("H5T", Element::Hydrogen()) );
		
        nameBondCenter("bondNext", "H5T/bond");

        setBiotypeIndex( "H5T", Biotype::get("Hydroxyl, RNA", "H5T", SimTK::Ordinality::Initial).getIndex() );
    }
	
	static void instantiateBiotypes() 
	{
		if ( ! Biotype::exists("Hydroxyl, RNA", "H5T", SimTK::Ordinality::Initial) )
			Biotype::defineBiotype(Element::Hydrogen(), 1, "Hydroxyl, RNA", "H5T", SimTK::Ordinality::Initial);
		
		// This oxygen biotype should get applied to the 5' oxygen of the nucleotide this hydroxyl attaches to
		if ( ! Biotype::exists("Hydroxyl, RNA", "O5*", SimTK::Ordinality::Initial) )
			Biotype::defineBiotype(Element::Oxygen(), 2, "Hydroxyl, RNA", "O5*", SimTK::Ordinality::Initial);
	}
	
};


/// Phosphate group at 5' end (beginning) of NA
class SimTK_MOLMODEL_EXPORT FivePrimeNaPhosphateGroup : public BiopolymerResidue
{
public:
	explicit FivePrimeNaPhosphateGroup(String name, String threeLetterCode = "Unk", char oneLetterCode = '?') 
    : BiopolymerResidue(name, threeLetterCode, oneLetterCode)
    {
		instantiateBiotypes();
		
		setBaseAtom( QuadrivalentAtom("P", Element::Phosphorus()) );
		
        bondAtom(UnivalentAtom("OP1", Element::Oxygen()), "P/bond4", 0.14800);
        bondAtom(UnivalentAtom("OP2", Element::Oxygen()), "P/bond3", 0.14800);
        bondAtom(UnivalentAtom("OP3", Element::Oxygen()), "P/bond2", 0.14800);
		
        nameBondCenter("bondNext", "P/bond1");

        setBiotypeIndex( "P", Biotype::get("Phosphate, RNA", "P", SimTK::Ordinality::Initial).getIndex() );
        setBiotypeIndex( "OP1", Biotype::get("Phosphate, RNA", "OP", SimTK::Ordinality::Initial).getIndex() );
        setBiotypeIndex( "OP2", Biotype::get("Phosphate, RNA", "OP", SimTK::Ordinality::Initial).getIndex() );
        setBiotypeIndex( "OP3", Biotype::get("Phosphate, RNA", "OP", SimTK::Ordinality::Initial).getIndex() );
    }
	
	static void instantiateBiotypes() 
	{
		if ( ! Biotype::exists("Phosphate, RNA", "P", SimTK::Ordinality::Initial) )
			Biotype::defineBiotype(Element::Phosphorus(), 4, "Phosphate, RNA", "P", SimTK::Ordinality::Initial);
		if ( ! Biotype::exists("Phosphate, RNA", "OP", SimTK::Ordinality::Initial) )
			Biotype::defineBiotype(Element::Oxygen(), 1, "Phosphate, RNA", "OP", SimTK::Ordinality::Initial);
	}
	
};


// Phosphodiester linkage between two NA nucleotides
class SimTK_MOLMODEL_EXPORT NaPhosphodiesterLinkage : public BiopolymerResidue
{
public:
    explicit NaPhosphodiesterLinkage(String name, String threeLetterCode = "Unk", char oneLetterCode = '?') 
        : BiopolymerResidue(name, threeLetterCode, oneLetterCode)
    {
		instantiateBiotypes();

        // TODO - set bond angles: 102.6 between ether oxygens, 119.9 between lone oxygen, 108.23 between mismatched pairs
		setBaseAtom( QuadrivalentAtom("P", Element::Phosphorus() 
				   ,  104.0*Deg2Rad // 03'-P-O5' from Gelbin et al
				   ,  107.9*Deg2Rad // 03'-P-OP? from Gelbin et al
				   ,  107.9*Deg2Rad // 03'-P-OP? from Gelbin et al
				   , -108.15*Deg2Rad // to make OP1-P-OP2 angle 119.6
				   ,  108.15*Deg2Rad // to make OP1-P-OP2 angle 119.6
				   ) );

        bondAtom(UnivalentAtom("OP1", Element::Oxygen()), "P/bond3", 0.14800);
        bondAtom(UnivalentAtom("OP2", Element::Oxygen()), "P/bond4", 0.14800);

        nameBondCenter("bondPrevious", "P/bond1");
        nameBondCenter("bondNext", "P/bond2");

        setDefaultInboardBondLength(0.16100);

        addCompoundSynonym("Phosphodiester, RNA");

		setBiotypeIndex("P", Biotype::get("Phosphodiester, RNA", "P").getIndex() );
		setBiotypeIndex("OP1", Biotype::get("Phosphodiester, RNA", "OP").getIndex() );
		setBiotypeIndex("OP2", Biotype::get("Phosphodiester, RNA", "OP").getIndex() );
    }
    
    
	static void instantiateBiotypes() 
	{
		if ( ! Biotype::exists("Phosphodiester, RNA", "P") )
			Biotype::defineBiotype(Element::Phosphorus(), 4, "Phosphodiester, RNA", "P");
		if ( ! Biotype::exists("Phosphodiester, RNA", "OP") )
			Biotype::defineBiotype(Element::Oxygen(), 1, "Phosphodiester, RNA", "OP");
	}
};


/// Phosphate group at 3' end (beginning) of RNA
class SimTK_MOLMODEL_EXPORT ThreePrimeNaHydroxylGroup : public BiopolymerResidue
{
public:
	explicit ThreePrimeNaHydroxylGroup(String name, String threeLetterCode = "Unk", char oneLetterCode = '?') 
    : BiopolymerResidue(name, threeLetterCode, oneLetterCode)
    {
		instantiateBiotypes();
		
		setBaseAtom( UnivalentAtom("H3T", Element::Hydrogen()) );
		
        nameBondCenter("bondNext", "H3T/bond");

        setBiotypeIndex( "H3T", Biotype::get("Hydroxyl, RNA", "H3T", SimTK::Ordinality::Final).getIndex() );
    }
	
	static void instantiateBiotypes() 
	{
		if ( ! Biotype::exists("Hydroxyl, RNA", "H3T", SimTK::Ordinality::Final) )
			Biotype::defineBiotype(Element::Hydrogen(), 1, "Hydroxyl, RNA", "H3T", SimTK::Ordinality::Final);
		
		// This oxygen biotype should get applied to the 3' oxygen of the nucleotide this hydroxyl attaches to
		if ( ! Biotype::exists("Hydroxyl, RNA", "O3*", SimTK::Ordinality::Final) )
			Biotype::defineBiotype(Element::Oxygen(), 2, "Hydroxyl, RNA", "O3*", SimTK::Ordinality::Final);
	}
	
};

/// Phosphate group at 3' end (end) of NA
class SimTK_MOLMODEL_EXPORT ThreePrimeNaPhosphateGroup : public BiopolymerResidue
{
public:
	explicit ThreePrimeNaPhosphateGroup(String name, String threeLetterCode = "Unk", char oneLetterCode = '?') 
    : BiopolymerResidue(name, threeLetterCode, oneLetterCode)
    {
		instantiateBiotypes();
		
		setBaseAtom( QuadrivalentAtom("P", Element::Phosphorus()) );
		
        bondAtom(UnivalentAtom("OP1", Element::Oxygen()), "P/bond4", 0.14800);
        bondAtom(UnivalentAtom("OP2", Element::Oxygen()), "P/bond3", 0.14800);
        bondAtom(UnivalentAtom("OP3", Element::Oxygen()), "P/bond2", 0.14800);
		
        nameBondCenter("bondNext", "P/bond1");

        setBiotypeIndex( "P", Biotype::get("Phosphate, RNA", "P", SimTK::Ordinality::Final).getIndex() );
        setBiotypeIndex( "OP1", Biotype::get("Phosphate, RNA", "OP", SimTK::Ordinality::Final).getIndex() );
        setBiotypeIndex( "OP2", Biotype::get("Phosphate, RNA", "OP", SimTK::Ordinality::Final).getIndex() );
        setBiotypeIndex( "OP3", Biotype::get("Phosphate, RNA", "OP", SimTK::Ordinality::Final).getIndex() );
    }
	
	static void instantiateBiotypes() 
	{
		if ( ! Biotype::exists("Phosphate, RNA", "P", SimTK::Ordinality::Final) )
			Biotype::defineBiotype(Element::Phosphorus(), 4, "Phosphate, RNA", "P", SimTK::Ordinality::Final);
		if ( ! Biotype::exists("Phosphate, RNA", "OP", SimTK::Ordinality::Final) )
			Biotype::defineBiotype(Element::Oxygen(), 1, "Phosphate, RNA", "OP", SimTK::Ordinality::Final);
	}
	
};





} // namespace SimTK

#endif // SimTK_MOLMODEL_NA_H_
