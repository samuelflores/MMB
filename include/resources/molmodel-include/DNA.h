#ifndef SimTK_MOLMODEL_DNA_H_
#define SimTK_MOLMODEL_DNA_H_

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
#include "molmodel/internal/NA.h"

namespace SimTK {

// RiboseCore leaves 4 bond centers open:
// 1) inboard center at O5* atom, for binding 5' phosphate
// 2) 

// Nucleoside has sugar and base but not phosphate
class SimTK_MOLMODEL_EXPORT DeoxyribonucleosideResidue : public RiboseCore {
public:
    explicit DeoxyribonucleosideResidue(String name, String threeLetterCode = "Unk", char oneLetterCode = '?')
        : RiboseCore(name, threeLetterCode, oneLetterCode)
    {
        // 2' hydroxyl group present in DNA but not in DNA
        //bondAtom(BivalentAtom("O2*", Element::Oxygen(), 108.50*Deg2Rad), "C2*/bond4", 0.1413, -60*Deg2Rad);
        //bondAtom(UnivalentAtom("2HO*", Element::Hydrogen()), "O2*/bond2", 0.0960);
        bondAtom(AliphaticHydrogen("H2*1"), "C2*/bond3");
        //bondAtom(AliphaticHydrogen("H2*2"), "C2*/bond4");

        bondAtom(UnivalentAtom("H2*2", Element::Hydrogen()), "C2*/bond4", 0.0960);

        // alternate atom name
        //nameAtom("O2'", "O2*");
        //nameAtom("HO2'", "2HO*");
        nameAtom("H2'1", "H2*1");
        nameAtom("H2'2", "H2*2");
    }
    
    // Create larger residue with phosphodiester prepended
    BiopolymerResidue withPhosphodiester() const 
    {
        const DeoxyribonucleosideResidue& residue = *this;
        
        NaPhosphodiesterLinkage po2(residue.getCompoundName(), residue.getPdbResidueName(), residue.getOneLetterCode());
        po2.bondCompound("nucleoside", residue, "bondNext");
        po2.inheritAtomNames("nucleoside");
        po2.inheritCompoundSynonyms(residue);
        po2.nameBondCenter("bondNext", "nucleoside/bondNext");

        po2.setPdbResidueNumber(residue.getPdbResidueNumber());
        po2.setPdbChainId(residue.getPdbChainId());

        po2.defineDihedralAngle("alpha", "bondPrevious", "O5*/bond2");
        po2.setDefaultDihedralAngle("alpha", 295.0*Deg2Rad); // A-DNA (Schneider et al 2004)

        return po2;
    }
};

// Nucleotide has sugar and base and phosphate
class SimTK_MOLMODEL_EXPORT DeoxyribonucleotideResidue : public DeoxyribonucleosideResidue {
public:
    // Factory method
    static DeoxyribonucleotideResidue create(char oneLetterCode);
    static DeoxyribonucleotideResidue create(const PdbResidue& pdbResidue, char chainId = ' ');

    explicit DeoxyribonucleotideResidue(String name, String threeLetterCode = "Unk", char oneLetterCode = '?')
        : DeoxyribonucleosideResidue(name, threeLetterCode, oneLetterCode)
    {
        defineDihedralAngle("beta", "O5*/bond1", "C5*/bond2");
        setDefaultDihedralAngle("beta", 173.0*Deg2Rad); //  A-DNA (Schneider et al 2004)

        defineDihedralAngle("gamma", "C5*/bond1", "C4*/bond3");
        setDefaultDihedralAngle("gamma", 54.0*Deg2Rad); //  A-DNA (Schneider et al 2004)

        defineDihedralAngle("delta", "C4*/bond1", "C3*/bond2");
        setDefaultDihedralAngle("delta", 80.0*Deg2Rad); //  A-DNA (Schneider et al 2004)

        defineDihedralAngle("epsilon", "C3*/bond1", "O3*/bond2");
        setDefaultDihedralAngle("epsilon", 210.0*Deg2Rad); //  A-DNA (Schneider et al 2004)

        nameBondCenter("bondBase", "C1*/bond3");
    }


    class Deoxyadenosine;
    class Deoxycytidine;
    class Deoxyguanosine;
    class Deoxythymidine;
    
    // Modified residues from tDNA
    class TwoNMethylDeoxyguanosine;
};

class ThymineBase  : public PyrimidineBaseCore
{
public:
    ThymineBase()
    {
        setBaseAtom( TrivalentAtom("N1",  Element::Nitrogen(), 125.8*Deg2Rad, 128.8*Deg2Rad));

        bondAtom( TrivalentAtom("C2", Element::Carbon(), 118.6*Deg2Rad, 120.9*Deg2Rad),  "N1/bond2", 0.1383, 180*Deg2Rad,   BondMobility::Rigid );
        bondAtom( UnivalentAtom("O2", Element::Oxygen()), "C2/bond3", 0.1229, 0*Deg2Rad, BondMobility::Rigid);

        bondAtom(TrivalentAtom("N3", Element::Nitrogen(), 126.4*Deg2Rad, 116.8*Deg2Rad), "C2/bond2", 0.1358, 0*Deg2Rad, BondMobility::Rigid);
        bondAtom(UnivalentAtom("H3", Element::Hydrogen()), "N3/bond3", 0.1010, 0*Deg2Rad, BondMobility::Rigid);

        bondAtom( TrivalentAtom("C4", Element::Carbon(), 114.0*Deg2Rad, 120.6*Deg2Rad),  "N3/bond2", 0.1388, 0*Deg2Rad,   BondMobility::Rigid );
        bondAtom( UnivalentAtom("O4", Element::Oxygen()), "C4/bond3", 0.1229, 0*Deg2Rad, BondMobility::Rigid);

        bondAtom( TrivalentAtom("C5", Element::Carbon(), 120.7*Deg2Rad, 119.70*Deg2Rad), "C4/bond2", 0.1444, 0*Deg2Rad, BondMobility::Rigid);
        //bondAtom( UnivalentAtom("H5", Element::Hydrogen()), "C5/bond3", 0.1080, 0*Deg2Rad, BondMobility::Rigid);
        bondCompound("thymineMethyl",MethylGroup(),"C5/bond3",0.150 ,0*Deg2Rad,BondMobility::Rigid); // the problem is here!
        nameAtom("C7" , "thymineMethyl/C");
        //nameAtom("H7" , "thymineMethyl/H1");
        nameAtom("H71" , "thymineMethyl/H1");
        nameAtom("H72" , "thymineMethyl/H2");
        nameAtom("H73" , "thymineMethyl/H3");

        // for some reason the tinker-amber99 parameter file doesn't have parameters for H71,H72,H73.  Using the H7 parameters.
        setBiotypeIndex( "H71", Biotype::get("Deoxythymidine", "H7", SimTK::Ordinality::Any).getIndex() );
        setBiotypeIndex( "H72", Biotype::get("Deoxythymidine", "H7", SimTK::Ordinality::Any).getIndex() );
        setBiotypeIndex( "H73", Biotype::get("Deoxythymidine", "H7", SimTK::Ordinality::Any).getIndex() );

        bondAtom( TrivalentAtom("C6", Element::Carbon(), 121.20*Deg2Rad, 119.70*Deg2Rad), "C5/bond2", 0.1350, 0*Deg2Rad, BondMobility::Rigid);
        bondAtom( UnivalentAtom("H6", Element::Hydrogen()), "C6/bond3", 0.1080, 0*Deg2Rad, BondMobility::Rigid);

        addRingClosingBond("C6/bond2", "N1/bond3", 0.1365, 0*Deg2Rad, BondMobility::Rigid);
    }
};



// TODO - these should not be classes, e.g. Deoxyadenosine should be an instance of DeoxyribonucleotideResidue



class  DeoxyribonucleotideResidue::Deoxyadenosine : public DeoxyribonucleotideResidue {
public:
    Deoxyadenosine() : DeoxyribonucleotideResidue("deoxyadenylate", "A  ", 'A') 
    {
        bondCompound("base", AdenineBase(), "bondBase", 0.14750, 200*Deg2Rad, BondMobility::Torsion);
        inheritAtomNames("base");

        defineDihedralAngle("chi", "O4*", "C1*", "N9", "C4");
        setDefaultDihedralAngle("chi", 199.0*Deg2Rad);

        addCompoundSynonym("Deoxyadenosine"); // for resolving biotypes
    }
};

class  DeoxyribonucleotideResidue::Deoxyguanosine : public DeoxyribonucleotideResidue {
public:
    Deoxyguanosine() : DeoxyribonucleotideResidue("deoxyguanylate", "G  ", 'G') 
    {
        bondCompound("base", GuanineBase(), "bondBase", 0.14710, 200*Deg2Rad, BondMobility::Torsion);
        inheritAtomNames("base");

        defineDihedralAngle("chi", "O4*", "C1*", "N9", "C4");
        setDefaultDihedralAngle("chi", 199.0*Deg2Rad);

        addCompoundSynonym("Deoxyguanosine"); // for resolving biotypes
    }
};

class  DeoxyribonucleotideResidue::Deoxycytidine : public DeoxyribonucleotideResidue {
public:
    Deoxycytidine() : DeoxyribonucleotideResidue("deoxycytidylate", "C  ", 'C') 
    {
        bondCompound("base", CytosineBase(), "bondBase", 0.14710, 200*Deg2Rad, BondMobility::Torsion);
        inheritAtomNames("base");

        defineDihedralAngle("chi", "O4*", "C1*", "N1", "C2");
        setDefaultDihedralAngle("chi", 199.0*Deg2Rad);

        addCompoundSynonym("Deoxycytidine"); // for resolving biotypes
    }
};

class  DeoxyribonucleotideResidue::Deoxythymidine : public DeoxyribonucleotideResidue {
public:
    Deoxythymidine() : DeoxyribonucleotideResidue("deoxythymidine", "T  ", 'T') 
    {
        bondCompound("base", ThymineBase(), "bondBase", 0.14710, 200*Deg2Rad, BondMobility::Torsion);
        inheritAtomNames("base");

        defineDihedralAngle("chi", "O4*", "C1*", "N1", "C2");
        setDefaultDihedralAngle("chi", 199.0*Deg2Rad);

        addCompoundSynonym("Deoxythymidine"); // for resolving biotypes
    }
};


class  DNA : public Biopolymer {
public:

	explicit DNA(const PdbStructure& pdbStructure) 
	{
		// Assume first chain of first model is wanted
		const PdbChain& pdbChain = pdbStructure.getModel(Pdb::ModelIndex(0)).getChain(Pdb::ChainIndex(0));
		initializeFromPdbChain(pdbChain);
	}

	explicit DNA(const PdbChain& pdbChain)
	{
		initializeFromPdbChain(pdbChain);
	}
    /**
     * /brief First parameter is just the DNA sequence in single-letter code.  AUGC are the only residue types supported at this time.  Second parameter, when set to zero, results in an DNA with no capping hydroxyls at the termini.  Instead the terminal residues have the same complement of atoms as residues in the interior of the chain.
     *
     */
    explicit DNA(const Sequence& seq, bool useCappingHydroxyls = true) 
    {
        String previousResidueName;
        Biotype::initializePopularBiotypes();

        // TODO - create 5' end cap
        for (int resi = 0; resi < (int)seq.size(); ++resi)
        {
        	
            DeoxyribonucleotideResidue residue = DeoxyribonucleotideResidue::create(seq[resi]);
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
            if ( (residue.hasAtom("O3*")) && ( Biotype::exists("Hydroxyl, DNA", "O3*", SimTK::Ordinality::Final)) )
                    residue.setBiotypeIndex( "O3*", Biotype::get("Hydroxyl, DNA", "O3*", SimTK::Ordinality::Final).getIndex() );
            }
            if ( (residue.hasAtom("H72")) && ( Biotype::exists("Deoxythymidine", "H7", SimTK::Ordinality::Any)) ){
                    std::cout<<__FILE__<<":"<<__LINE__<<": setting biotype index to : "<< SimTK::String (Biotype::get("Deoxythymidine", "H7", SimTK::Ordinality::Any).getIndex())<<std::endl;
                    residue.setBiotypeIndex( "H72", Biotype::get("Deoxythymidine", "H7", SimTK::Ordinality::Any).getIndex() );}

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
                            if ( (residue.hasAtom("O5*")) && ( Biotype::exists("Hydroxyl, DNA", "O5*", SimTK::Ordinality::Initial)) )
                                    residue.setBiotypeIndex( "O5*", Biotype::get("Hydroxyl, DNA", "O5*", SimTK::Ordinality::Initial).getIndex() );
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
    DNA& setDNABondMobility (BondMobility::Mobility  mobility, ResidueInfo::Index startResidue, ResidueInfo::Index endResidue){
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

    Vec3 calcDNABaseNormal  (const State & state,int residueNumber, int bodyOrGroundFrame, const SimbodyMatterSubsystem & matter) const {
	    const ResidueInfo& myResidue = getResidue(ResidueInfo::Index(residueNumber));
            String myResidueName = myResidue.getPdbResidueName();
	    //cout<<"[DNA.h:calcDNABaseNormal] myResidueName ="<<myResidueName<<endl;
            //String myResidueName = (updResidue(Compound::Index(residueNumber))).getPdbResidueName();
	    assert((myResidueName.compare("A  ") == 0) || (myResidueName.compare("G  ") == 0)  || (myResidueName.compare("T  ") == 0) || (myResidueName.compare("C  ") == 0));
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
            
            DeoxyribonucleotideResidue residue = DeoxyribonucleotideResidue::create(pdbResidue, pdbChain.getChainId());
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

#endif // SimTK_MOLMODEL_DNA_H_
