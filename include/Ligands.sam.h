#include "molmodel/internal/common.h"
#include "molmodel/internal/Compound.h"
#include <iostream>
#include <fstream>
//#include <ios>

namespace SimTK 
{

//
//  This ligand is from 1UUI.pdb.  Reference is:
//  1. Davis, B. et al. Rational design of inhibitors of HIV-1 TAR RNA through the stabilisation of electrostatic "hot spots". J. Mol. Biol 336, 343-56(2004).
//

class P12 : public Compound {
public:
    //takes a TinkerDuMMForceFieldSubsystem as a reference argument so as to define all Biotypes, etc within the class.	
    P12(TinkerDuMMForceFieldSubsystem &dumm)
    {
 
        setBaseCompound("methyl", MethylGroup() );
        bondAtom(BivalentAtom("OA1",Element::Oxygen(),109.5*Deg2Rad),   "methyl/bond" ,.141);// third arg of bivalentatom is a bond angle in rads ;. angle betw bond1 and bond2.  third argument of bondAtom is a distance.
        bondAtom(TrivalentAtom("C4",   Element::Carbon()),"OA1/bond2",.141);//120*Deg2Rad);
        bondCompound("C3compound",AromaticSixMemberedCHGroup(),"C4/bond2",.141);//4th param is a distance.
        bondCompound("C5compound",AromaticSixMemberedCHGroup(),"C4/bond3",.141);
        bondCompound("C6compound",AromaticSixMemberedCHGroup(),"C5compound/C/bond2",.141);
        bondAtom(TrivalentAtom("C1",   Element::Carbon()),"C6compound/C/bond2",.141);
        bondAtom(TrivalentAtom("C2",   Element::Carbon()),"C3compound/C/bond2",.141);
	addRingClosingBond("C1/bond3",  "C2/bond3", 0.1410);
        nameAtom("C5", "C5compound/C");
        nameAtom("C6", "C6compound/C");
        defineDihedralAngle("C4C3CC2C1","C4","C3compound/C","C2"  ,"C1");//OA/bond2", "C1/bond2");
        setDefaultDihedralAngle("C4C3CC2C1", 0.00*Deg2Rad);
        bondCompound("CAcompound",MethyleneGroup(),"C2/bond2" ,.141);
        nameAtom("CA", "CAcompound/C");
        nameAtom("HA1", "CAcompound/H1");
        nameAtom("HA2", "CAcompound/H2");
        defineDihedralAngle("theta", "C4","C3compound/C","C2","CA");//OA/bond2", "C1/bond2");
        setDefaultDihedralAngle("theta",140.0 *Deg2Rad);
        bondAtom(QuadrivalentAtom("NB",Element::Nitrogen()),"CAcompound/bond2",.141);// last param is distance; set to something more real // trivalent antom takes 2 bond angles, betw bond1 and bond2, and betw bond1 and bond3.  //trivalents default to planar
        bondAtom(   UnivalentAtom("HB",Element::Hydrogen()),"NB/bond4",.141);// last param is distance; set to something more real // trivalent antom takes 2 bond angles, betw bond1 and bond2, and betw bond1 and bond3.  //trivalents default to planar
        defineDihedralAngle("omega", "C3compound/C","C2","CA","NB");//OA/bond2", "C1/bond2");
        setDefaultDihedralAngle("omega", 40.0 *Deg2Rad);
        bondAtom(BivalentAtom("OA",Element::Oxygen(),109.5*Deg2Rad),   "C1/bond2"    ,.141);

        bondCompound( "CBcompound",MethyleneGroup(),"OA/bond2" ,.141);
        nameAtom( "CB", "CBcompound/C");
        nameAtom( "HB1","CBcompound/H1");
        nameAtom( "HB2","CBcompound/H2");

        defineDihedralAngle("alpha", "C6compound/C","C1","OA","CB");//OA/bond2", "C1/bond2");
        setDefaultDihedralAngle("alpha", 40.0 *Deg2Rad);

        //bondAtom(QuadrivalentAtom("CG",Element::Carbon()),  "CBcompound/bond2" ,.141);
        bondCompound("CGcompound",MethyleneGroup(),"CBcompound/bond2" ,.141);
        nameAtom("CG", "CGcompound/C");
        nameAtom("HG1","CGcompound/H1");
        nameAtom("HG2","CGcompound/H2");
        //bondAtom(QuadrivalentAtom("CD",Element::Carbon()),  "CGcompound/bond2" ,.141);
        bondCompound( "CDcompound",MethyleneGroup(),"CGcompound/bond2" ,.141);
        nameAtom("CD", "CDcompound/C");
        nameAtom("HD1","CDcompound/H2");
        nameAtom("HD2","CDcompound/H1");
        defineDihedralAngle("zeta", "C1","OA","CB","CG");//OA/bond2", "C1/bond2");
        setDefaultDihedralAngle("zeta" , 40.0 *Deg2Rad);
        //bondAtom(QuadrivalentAtom("NE",Element::Nitrogen()),"CDcompound/bond2",.141);// last param is distance; set to something more real // trivalent antom takes 2 bond angles, betw bond1 and bond2, and betw bond1 and bond3.  //trivalents default to planar
        bondCompound("NEcompound",PrimaryAmineGroup(),"CDcompound/bond2" ,.141);
        nameAtom("NE", "NEcompound/N");
        nameAtom("HE1","NEcompound/H1");
        nameAtom("HE2","NEcompound/H2");
        nameAtom("HE3","NEcompound/H3");

        defineDihedralAngle("gamma", "C3compound/C","C2","CA","NB");//OA/bond2", "C1/bond2");
        setDefaultDihedralAngle("gamma", 140.0*Deg2Rad);
        //bondAtom(QuadrivalentAtom("CG1",Element::Carbon()),   "NB/bond3"    ,.141);
        bondCompound("CG1compound",MethyleneGroup(),"NB/bond2" ,.141);
        nameAtom("CG1", "CG1compound/C");
        nameAtom("HG11","CG1compound/H2");
        nameAtom("HG12","CG1compound/H1");
        bondCompound("CG2compound",MethyleneGroup(),"NB/bond3" ,.141);
        nameAtom("CG2", "CG2compound/C");
        nameAtom("HG21","CG2compound/H1");
        nameAtom("HG22","CG2compound/H2");
        bondCompound("CD1compound",MethyleneGroup(),"CG1compound/bond2" ,.141);
        nameAtom("CD1", "CD1compound/C");
        nameAtom("HD11","CD1compound/H1");
        nameAtom("HD12","CD1compound/H2");
        bondCompound("CD2compound",MethyleneGroup(),"CG2compound/bond2" ,.141);
        nameAtom("CD2", "CD2compound/C");
        nameAtom("HD21","CD2compound/H2");
        nameAtom("HD22","CD2compound/H1");
/**/
        bondAtom(TrivalentAtom("NE1",Element::Nitrogen()),"CD1compound/bond2",.141);
        bondAtom(   TrivalentAtom("CZ" ,Element::Carbon()),    "NE1/bond2"    ,.141);
	addRingClosingBond("CD2compound/bond2",  "NE1/bond3", 0.1410);
        bondAtom(TrivalentAtom("NH1",Element::Nitrogen()), "CZ/bond2",.141);
        bondAtom(UnivalentAtom("HH11",Element::Hydrogen()),"NH1/bond2",.141);
        bondAtom(UnivalentAtom("HH12",Element::Hydrogen()),"NH1/bond3",.141);
        bondAtom(TrivalentAtom("NH2",Element::Nitrogen()), "CZ/bond3",.141);
        bondAtom(UnivalentAtom("HH21",Element::Hydrogen()),"NH2/bond2",.141);
        bondAtom(UnivalentAtom("HH22",Element::Hydrogen()),"NH2/bond3",.141);

        defineDihedralAngle("delta", "NH1" ,"CZ","NE1","CD1");//OA/bond2", "C1/bond2");
        setDefaultDihedralAngle("delta", 70.00*Deg2Rad);

        defineDihedralAngle("beta" , "C1","OA","CB","CG");//OA/bond2", "C1/bond2");
        setDefaultDihedralAngle("beta" , 40.0 *Deg2Rad);
        
   
        

        nameAtom("CB1", "methyl/C");// just an alias or shorthand
        nameAtom("HB11", "methyl/H1");
        nameAtom("HB12", "methyl/H2");
        nameAtom("HB13", "methyl/H3");
        nameAtom("C3", "C3compound/C");
        nameAtom("H3", "C3compound/H");
        nameAtom("H5", "C5compound/H");
        nameAtom("H6", "C6compound/H");
       // nameAtom("O", "OA1");


        if (! Biotype::exists("P12", "CB1") )
            Biotype::defineBiotype(Element::Carbon(), 4, "P12", "CB1");
        if (! Biotype::exists("P12", "HC") ) // these names unrelated to the above
            Biotype::defineBiotype(Element::Hydrogen(), 1, "P12", "HC");
        if (! Biotype::exists("P12", "OA1") )
            Biotype::defineBiotype(Element::Oxygen(), 2, "P12", "OA1");
        if (! Biotype::exists("P12", "N") )
            Biotype::defineBiotype(Element::Nitrogen(), 3, "P12", "N");//this can be called for an atom that is not in use.
        if (! Biotype::exists("P12", "H") )
            Biotype::defineBiotype(Element::Hydrogen(), 3, "P12", "H");//this can be called for an atom that is not in use.
        if (! Biotype::exists("P12", "NE1") )
            Biotype::defineBiotype(Element::Nitrogen(), 3, "P12", "NE1");//this can be called for an atom that is not in use.
        if (! Biotype::exists("P12", "C4") )
            Biotype::defineBiotype(Element::Carbon(), 3, "P12", "C4");
        if (! Biotype::exists("P12", "C2") )
            Biotype::defineBiotype(Element::Carbon(), 3, "P12", "C2");
        if (! Biotype::exists("P12", "C3") )
            Biotype::defineBiotype(Element::Carbon(), 3, "P12","C3");
        if (! Biotype::exists("P12", "H3") )
            Biotype::defineBiotype(Element::Hydrogen(), 1, "P12", "H3");
        if (! Biotype::exists("P12", "CD1") )
            Biotype::defineBiotype(Element::Carbon(),4 , "P12", "CD1");
        //if (! Biotype::exists("P12", "HO") )
        //    Biotype::defineBiotype(Element::Hydrogen(), 1, "P12", "HO");

        setBiotypeIndex( "CB1", Biotype::get("P12", "CB1").getIndex() );
        setBiotypeIndex( "HB11", Biotype::get("P12", "HC").getIndex() );// use original atom name or alias.
        setBiotypeIndex( "HB12", Biotype::get("P12", "HC").getIndex() );
        setBiotypeIndex( "HB13", Biotype::get("P12", "HC").getIndex() );
        setBiotypeIndex( "HA1", Biotype::get("P12", "HC").getIndex() );
        setBiotypeIndex( "HA2", Biotype::get("P12", "HC").getIndex() );
        setBiotypeIndex( "HB1", Biotype::get("P12", "HC").getIndex() );
        setBiotypeIndex( "HB2", Biotype::get("P12", "HC").getIndex() );
        setBiotypeIndex( "HD1", Biotype::get("P12", "HC").getIndex() );
        setBiotypeIndex( "HD2", Biotype::get("P12", "HC").getIndex() );
        setBiotypeIndex( "HG1", Biotype::get("P12", "HC").getIndex() );
        setBiotypeIndex( "HG2", Biotype::get("P12", "HC").getIndex() );
        setBiotypeIndex("HG11", Biotype::get("P12", "HC").getIndex() );
        setBiotypeIndex("HG12", Biotype::get("P12", "HC").getIndex() );
        setBiotypeIndex("HG21", Biotype::get("P12", "HC").getIndex() );
        setBiotypeIndex("HG22", Biotype::get("P12", "HC").getIndex() );
        setBiotypeIndex("HD11", Biotype::get("P12", "HC").getIndex() );
        setBiotypeIndex("HD12", Biotype::get("P12", "HC").getIndex() );
        setBiotypeIndex("HD21", Biotype::get("P12", "HC").getIndex() );
        setBiotypeIndex("HD22", Biotype::get("P12", "HC").getIndex() );
        setBiotypeIndex("HE1" , Biotype::get("P12", "H").getIndex() );
        setBiotypeIndex("HE2" , Biotype::get("P12", "H").getIndex() );
        setBiotypeIndex("HE3" , Biotype::get("P12", "H").getIndex() );
        setBiotypeIndex("HB"  , Biotype::get("P12", "H").getIndex() );
        setBiotypeIndex("HH11", Biotype::get("P12", "H").getIndex() );
        setBiotypeIndex("HH12", Biotype::get("P12", "H").getIndex() );
        setBiotypeIndex("HH21", Biotype::get("P12", "H").getIndex() );
        setBiotypeIndex("HH22", Biotype::get("P12", "H").getIndex() );
        setBiotypeIndex( "OA1", Biotype::get("P12", "OA1").getIndex() );
	//cout<<"check OA"<<endl;
        setBiotypeIndex( "OA", Biotype::get("P12", "OA1").getIndex() );
        setBiotypeIndex( "NE", Biotype::get("P12", "N").getIndex() );//this can only be called for an existing atom.
        setBiotypeIndex( "NE1", Biotype::get("P12", "N").getIndex() );//this can only be called for an existing atom.
        setBiotypeIndex( "NH1", Biotype::get("P12", "N").getIndex() );//this can only be called for an existing atom.
        setBiotypeIndex( "NH2", Biotype::get("P12", "N").getIndex() );//this can only be called for an existing atom.
        setBiotypeIndex( "NB", Biotype::get("P12", "N").getIndex() );//this can only be called for an existing atom.
        setBiotypeIndex( "C4", Biotype::get("P12", "C4").getIndex() );//FIRST argument is atom name within the particular compound
        setBiotypeIndex( "C1", Biotype::get("P12", "C4").getIndex() );//FIRST argument is atom name within the particular compound
        //setBiotypeIndex( "HO", Biotype::get("P12", "HO").getIndex() );
        setBiotypeIndex( "C3", Biotype::get("P12", "C3").getIndex() );//FIRST argument is atom name within the particular compound
        setBiotypeIndex( "H3", Biotype::get("P12", "H3").getIndex() );//FIRST argument is atom name within the particular compound
        setBiotypeIndex( "C5", Biotype::get("P12", "C3").getIndex() );//FIRST argument is atom name within the particular compound
        setBiotypeIndex( "H5", Biotype::get("P12", "H3").getIndex() );//FIRST argument is atom name within the particular compound
        setBiotypeIndex( "C6", Biotype::get("P12", "C3").getIndex() );//FIRST argument is atom name within the particular compound
        setBiotypeIndex( "H6", Biotype::get("P12", "H3").getIndex() );//FIRST argument is atom name within the particular compound
        setBiotypeIndex( "C2", Biotype::get("P12", "C2").getIndex() );//FIRST argument is atom name within the particular compound
        setBiotypeIndex( "CB", Biotype::get("P12", "C2").getIndex() );//FIRST argument is atom name within the particular compound
        setBiotypeIndex( "CG", Biotype::get("P12", "C2").getIndex() );//FIRST argument is atom name within the particular compound
        setBiotypeIndex( "CD", Biotype::get("P12", "C2").getIndex() );//FIRST argument is atom name within the particular compound
        setBiotypeIndex( "CA", Biotype::get("P12", "C2").getIndex() );//FIRST argument is atom name within the particular compound
        setBiotypeIndex("CG1", Biotype::get("P12", "C2").getIndex() );//FIRST argument is atom name within the particular compound
        setBiotypeIndex("CG2", Biotype::get("P12", "C2").getIndex() );//FIRST argument is atom name within the particular compound
        setBiotypeIndex("CD1", Biotype::get("P12", "C2").getIndex() );//FIRST argument is atom name within the particular compound
        setBiotypeIndex("CD2", Biotype::get("P12", "C2").getIndex() );//FIRST argument is atom name within the particular compound
        setBiotypeIndex("CZ" , Biotype::get("P12", "C2").getIndex() );//FIRST argument is atom name within the particular compound

/*
    DuMM::AtomClassIndex amber1CTAtomClassIndex(1); 
    DuMM::AtomClassIndex amber3CAAtomClassIndex(1); 
    DuMM::AtomClassIndex amber9CBAtomClassIndex(1); 
    DuMM::AtomClassIndex amber19N2AtomClassIndex(1); 
    DuMM::AtomClassIndex amber23OSAtomClassIndex(1); 
    DuMM::AtomClassIndex amber29H1AtomClassIndex(1); 
    DuMM::AtomClassIndex amber33HAAtomClassIndex(1); 
    DuMM::AtomClassIndex amber34HCAtomClassIndex(1); 
*/

    DuMM::setBiotypeChargedAtomType(DuMM:ChargedAtomTypeIndex(37  ), Biotype::get ("P12","CB1").getIndex());
    DuMM::setBiotypeChargedAtomType(DuMM:ChargedAtomTypeIndex(38  ), Biotype::get ("P12","HB11).getIndex());
    DuMM::setBiotypeChargedAtomType(DuMM:ChargedAtomTypeIndex(34  ), Biotype::get ("P12","HC" ).getIndex());
    DuMM::setBiotypeChargedAtomType(DuMM:ChargedAtomTypeIndex(371 ), Biotype::get ("P12","HE1").getIndex());

    dumm.defineChargedAtomType(DuMM::AtomClassIndex(20),"P12 HB",DuMM::ChargedAtomTypeIndex(4000),.6)	
    DuMM::setBiotypeChargedAtomType(DuMM:ChargedAtomTypeIndex(4000 ), Biotype::get ("P12","HB").getIndex());

    DuMM::setBiotypeChargedAtomType(DuMM:ChargedAtomTypeIndex(303 ), Biotype::get ("P12","HH11).getIndex());
    DuMM::setBiotypeChargedAtomType(DuMM:ChargedAtomTypeIndex(1068), Biotype::get ("P12","OA1").getIndex());
    DuMM::setBiotypeChargedAtomType(DuMM:ChargedAtomTypeIndex(368 ), Biotype::get ("P12","NE" ).getIndex());
    DuMM::setBiotypeChargedAtomType(DuMM:ChargedAtomTypeIndex(1078), Biotype::get ("P12","NE1").getIndex());
    DuMM::setBiotypeChargedAtomType(DuMM:ChargedAtomTypeIndex(302 ), Biotype::get ("P12","NH1").getIndex());
    DuMM::setBiotypeChargedAtomType(DuMM:ChargedAtomTypeIndex(130 ), Biotype::get ("P12","C4" ).getIndex());
    DuMM::setBiotypeChargedAtomType(DuMM:ChargedAtomTypeIndex(133 ), Biotype::get ("P12","C3" ).getIndex());
    DuMM::setBiotypeChargedAtomType(DuMM:ChargedAtomTypeIndex(134 ), Biotype::get ("P12","H3" ).getIndex());
    DuMM::setBiotypeChargedAtomType(DuMM:ChargedAtomTypeIndex(1049), Biotype::get ("P12","C2" ).getIndex());
    DuMM::setBiotypeChargedAtomType(DuMM:ChargedAtomTypeIndex(33  ), Biotype::get ("P12","CD" ).getIndex());
    DuMM::setBiotypeChargedAtomType(DuMM:ChargedAtomTypeIndex(301 ), Biotype::get ("P12","CZ" ).getIndex());
 
 for (int r =0 ; r<getNBonds(); r++)    //)//(Compound::Index(q)).
    {
      	setBondMobility(BondMobility::Free ,Compound::BondIndex(r));
    }

    setPdbResidueNumber(1046);
    setPdbResidueName("P12") ;
    setPdbChainId('B');

    std::ifstream inFileStream("1UUI.pdb",std::ifstream::in);
    	
    PdbStructure pdbStructure(inFileStream);
    Compound::AtomTargetLocations atomTargets = createAtomTargets(pdbStructure);
    std::cout<<"atomtargest.szie "<<atomTargets.size()<<std::endl;
    matchDefaultBondLengths(atomTargets);
    matchDefaultBondAngles(atomTargets);
    matchDefaultAtomChirality(atomTargets);
    matchDefaultDihedralAngles(atomTargets);
    matchDefaultTopLevelTransform(atomTargets);
    //int mpResidueNumber = 1046;
    }
   /* 
*/

};
}
