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




        setAtomBiotype("CB1" ,"P12","CB1" );
        setAtomBiotype("HB11","P12","HB11");
        setAtomBiotype("HB12","P12","HB11");
        setAtomBiotype("HB13","P12","HB11");
        setAtomBiotype("HA1" ,"P12","HC"  );
        setAtomBiotype("HA2" ,"P12","HC"  );
        setAtomBiotype("HB1" ,"P12","HC"  );
        setAtomBiotype("HB2" ,"P12","HC"  );
        //setBiotypeIndex( "HD1", Biotype::get("P12", "HC").getIndex() );
        setAtomBiotype("HD1" ,"P12","HC"  );
        //setBiotypeIndex( "HD2", Biotype::get("P12", "HC").getIndex() );
        setAtomBiotype("HD2" ,"P12","HC"  );
        //setBiotypeIndex( "HG1", Biotype::get("P12", "HC").getIndex() );
        setAtomBiotype("HG1" ,"P12","HC"  );
        //setBiotypeIndex( "HG2", Biotype::get("P12", "HC").getIndex() );
        setAtomBiotype("HG2" ,"P12","HC"  );
        //setBiotypeIndex("HG11", Biotype::get("P12", "HC").getIndex() );
        setAtomBiotype("HG11","P12","HC"  );
        //setBiotypeIndex("HG12", Biotype::get("P12", "HC").getIndex() );
        setAtomBiotype("HG12","P12","HC"  );
        //setBiotypeIndex("HG21", Biotype::get("P12", "HC").getIndex() );
        setAtomBiotype("HG21","P12","HC"  );
        //setBiotypeIndex("HG22", Biotype::get("P12", "HC").getIndex() );
        setAtomBiotype("HG22","P12","HC"  );
        //setBiotypeIndex("HD11", Biotype::get("P12", "HC").getIndex() );
        setAtomBiotype("HD11","P12","HC"  );
        //setBiotypeIndex("HD12", Biotype::get("P12", "HC").getIndex() );
        setAtomBiotype("HD12","P12","HC"  );
        //setBiotypeIndex("HD21", Biotype::get("P12", "HC").getIndex() );
        setAtomBiotype("HD21","P12","HC"  );
        //setBiotypeIndex("HD22", Biotype::get("P12", "HC").getIndex() );
        setAtomBiotype("HD22","P12","HC"  );
        setAtomBiotype("HE1" ,"P12", "HE1");
        //setBiotypeIndex("HE1" , Biotype::get("P12", "HE1").getIndex() );
        //setBiotypeIndex("HE2" , Biotype::get("P12", "HE1").getIndex() );
        setAtomBiotype("HE2" ,"P12","HE1" );
        //setBiotypeIndex("HE3" , Biotype::get("P12", "HE1").getIndex() );
        setAtomBiotype("HE3" ,"P12","HE1" );
        setAtomBiotype("HA1" ,"P12","HC"  );
        //setBiotypeIndex("HH11", Biotype::get("P12", "HH11").getIndex() );
        setAtomBiotype("HH11","P12","HH11");
        //setBiotypeIndex("HH12", Biotype::get("P12", "HH11").getIndex() );
        setAtomBiotype("HH12","P12","HH11");
        //setBiotypeIndex("HH21", Biotype::get("P12", "HH11").getIndex() );
        setAtomBiotype("HH21","P12","HH11");
        //setBiotypeIndex("HH22", Biotype::get("P12", "HH11").getIndex() );
        setAtomBiotype("HH22","P12","HH11");
        //setBiotypeIndex( "OA1", Biotype::get("P12", "OA1").getIndex() );
        setAtomBiotype("OA1" ,"P12","OA1" );
	//cout<<"check OA"<<endl;
        //setBiotypeIndex( "OA", Biotype::get("P12", "OA1").getIndex() );
        setAtomBiotype("OA" ,"P12","OA1" );
        //setBiotypeIndex( "NE", Biotype::get("P12", "NE").getIndex() );//this can only be called for an existing atom.
        setAtomBiotype("NE"  ,"P12", "NE" );
        //setBiotypeIndex( "NE1", Biotype::get("P12", "NE1").getIndex() );//this can only be called for an existing atom.
        setAtomBiotype("NE1"  ,"P12", "NE1" );
        //setBiotypeIndex( "NH1", Biotype::get("P12", "NH1").getIndex() );//this can only be called for an existing atom.
        setAtomBiotype("NH1" ,"P12", "NH1");
        //setBiotypeIndex( "NH2", Biotype::get("P12", "NH1").getIndex() );//this can only be called for an existing atom.
        setAtomBiotype("NH2" ,"P12","NH1" );
        //setAtomBiotype("NE"  ,"P12", "NE" );
        //setBiotypeIndex( "NB", Biotype::get("P12", "NE").getIndex() );//this can only be called for an existing atom.
        setAtomBiotype("NB"  ,"P12","NE"  );
        //setBiotypeIndex( "C4", Biotype::get("P12", "C4").getIndex() );//FIRST argument is atom name within the particular compound
        setAtomBiotype("C4"  ,"P12","C4"  );
        //setBiotypeIndex( "C1", Biotype::get("P12", "C4").getIndex() );//FIRST argument is atom name within the particular compound
        //setBiotypeIndex( "C1","P12", "C4");
        setAtomBiotype("C1"  ,"P12","C4"  );
        setAtomBiotype("C3"  ,"P12", "C3" );
        setAtomBiotype("H3"  ,"P12", "H3" );
        //setBiotypeIndex( "C5", Biotype::get("P12", "C3").getIndex() );//FIRST argument is atom name within the particular compound
        setAtomBiotype("C5"  ,"P12","C3"  );
        //setBiotypeIndex( "H5", Biotype::get("P12", "H3").getIndex() );//FIRST argument is atom name within the particular compound
        setAtomBiotype("H5"  ,"P12","H3"  );
        //setBiotypeIndex( "C6", Biotype::get("P12", "C3").getIndex() );//FIRST argument is atom name within the particular compound
        setAtomBiotype("C6"  ,"P12","C3"  );
        //setBiotypeIndex( "H6", Biotype::get("P12", "H3").getIndex() );//FIRST argument is atom name within the particular compound
        setAtomBiotype("H6"  ,"P12","H3"  );
        //setBiotypeIndex( "C2", Biotype::get("P12", "C2").getIndex() );//FIRST argument is atom name within the particular compound
        setAtomBiotype("C2"  ,"P12", "C3" );
        setAtomBiotype("CB"  ,"P12", "CD" );//setBiotypeIndex( "CB", Biotype::get("P12", "CD").getIndex() );//FIRST argument is atom name within the particular compound
        setAtomBiotype("CD"  ,"P12", "CD" );
        //setBiotypeIndex( "CG", Biotype::get("P12", "CD").getIndex() );//FIRST argument is atom name within the particular compound
        setAtomBiotype("CG"  ,"P12","CD"  );
        //setBiotypeIndex( "CD", Biotype::get("P12", "CD").getIndex() );//FIRST argument is atom name within the particular compound
        setAtomBiotype("CD"  ,"P12","CD"  );
        setAtomBiotype("CA"  ,"P12","CD"  );
        setAtomBiotype("CG1" ,"P12","CD"  );
        setAtomBiotype("CG2" ,"P12","CD"  );
        setAtomBiotype("CD1","P12","CD");//setBiotypeIndex("CD1", Biotype::get("P12", "CD").getIndex() );//FIRST argument is atom name within the particular compound
        setAtomBiotype("CD2" ,"P12","CD"  );
        setAtomBiotype("CZ"  ,"P12", "CZ" );
        setAtomBiotype("C4"  ,"P12","C4"  );

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

    dumm.setBiotypeChargedAtomType(DuMM::ChargedAtomTypeIndex(37),Biotype::get("P12","CB1").getIndex());
    dumm.setBiotypeChargedAtomType(DuMM::ChargedAtomTypeIndex(38),Biotype::get("P12","HB11").getIndex());
    dumm.setBiotypeChargedAtomType(DuMM::ChargedAtomTypeIndex(34),Biotype::get("P12","HC").getIndex());
    dumm.setBiotypeChargedAtomType(DuMM::ChargedAtomTypeIndex(371),Biotype::get("P12","HE1").getIndex());
    cout << "setting parameters for P12 HB"<<endl;
        setAtomBiotype("HB","P12","HB");
	if (! dumm.hasChargedAtomType(DuMM::ChargedAtomTypeIndex(4000))) 
	{
    	dumm.defineChargedAtomType(DuMM::ChargedAtomTypeIndex(4000),"P12 HB",DuMM::AtomClassIndex(29),.6);	
    	dumm.setBiotypeChargedAtomType(DuMM::ChargedAtomTypeIndex(4000),Biotype::get("P12","HB").getIndex());

        }		
        /*dumm.setBiotypeChargedAtomType(DuMM::ChargedAtomTypeIndex(303),Biotype::get("P12","HB").getIndex());
    if (! Biotype::exists("P12", "HB" ) ) 
    {
	cout<<"apparently biotype P12 HB doesn't exist.."<<endl;
        //setAtomBiotype("HB","P12","HB");
        Biotype::defineBiotype(Element::Hydrogen(),1,"P12","HB");
    	dumm.defineChargedAtomType(DuMM::ChargedAtomTypeIndex(4000),"P12 HB",DuMM::AtomClassIndex(29),.6);	
    	dumm.setBiotypeChargedAtomType(DuMM::ChargedAtomTypeIndex(4000),Biotype::get("P12","HB").getIndex());
    } else 
	cout << "Biotype P12 HB exists"<<endl;
*/
    dumm.setBiotypeChargedAtomType(DuMM::ChargedAtomTypeIndex(303),Biotype::get("P12","HH11").getIndex());
    dumm.setBiotypeChargedAtomType(DuMM::ChargedAtomTypeIndex(1068),Biotype::get("P12","OA1").getIndex());
    dumm.setBiotypeChargedAtomType(DuMM::ChargedAtomTypeIndex(368),Biotype::get("P12","NE").getIndex());
    dumm.setBiotypeChargedAtomType(DuMM::ChargedAtomTypeIndex(1078),Biotype::get("P12","NE1").getIndex());
    dumm.setBiotypeChargedAtomType(DuMM::ChargedAtomTypeIndex(302),Biotype::get("P12","NH1").getIndex());
    dumm.setBiotypeChargedAtomType(DuMM::ChargedAtomTypeIndex(130),Biotype::get("P12","C4").getIndex());
    dumm.setBiotypeChargedAtomType(DuMM::ChargedAtomTypeIndex(133),Biotype::get("P12","C3").getIndex());
    dumm.setBiotypeChargedAtomType(DuMM::ChargedAtomTypeIndex(134),Biotype::get("P12","H3").getIndex());
    //dumm.setBiotypeChargedAtomType(DuMM::ChargedAtomTypeIndex(1049),Biotype::get("P12","C2").getIndex());
    dumm.setBiotypeChargedAtomType(DuMM::ChargedAtomTypeIndex(33),Biotype::get("P12","CD").getIndex());
    dumm.setBiotypeChargedAtomType(DuMM::ChargedAtomTypeIndex(301),Biotype::get("P12","CZ").getIndex());
 
 for (int r =0 ; r<getNBonds(); r++)    //)//(Compound::Index(q)).
    {
      	setBondMobility(BondMobility::Free ,Compound::BondIndex(r));
    }

    setPdbResidueNumber(1046);
    setPdbResidueName("P12") ;
    setPdbChainId('B');

    std::ifstream inFileStream("1UUI.pdb",std::ifstream::in);
    //assert(inFileStream.good()); 	
    //assert(inFileStream.bad()); 	
    //assert(0);
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
