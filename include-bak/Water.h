#include "SimTKmolmodel.h"
#include <iostream>
#include <fstream>
//#include <ios>




namespace SimTK 
{




class Water : public Compound { public:
        Water(DuMMForceFieldSubsystem &dumm)
        {   
 

     if (!dumm.hasAtomClass(DuMM::AtomClassIndex(400)))
     {	

     	dumm.defineAtomClass_KA(
         	DuMM::AtomClassIndex(400),
         	"TIP3P Hydrogen"    ,
         	1 , //element number
         	1, //expected valence
         	.0001 ,//from vdw parms
         	(.0000)
         	);
     }	
     if (! dumm.hasChargedAtomType(DuMM::ChargedAtomTypeIndex(8000)))
     {	
     	dumm.defineChargedAtomType(
         	DuMM::ChargedAtomTypeIndex(8000),
         	"TIP3P Hydrogen"    ,
         	DuMM::AtomClassIndex(400),
         	.417                   //partial charge
         	);
     }
	
     if (! Biotype::exists("TIP3P"   ,"Hydrogen" ))
         Biotype::defineBiotype(Element::Hydrogen(), 1, "TIP3P",   "Hydrogen"  ); // second arg is valence
                                                                                 // a residue name, second is an atom name
     dumm.setBiotypeChargedAtomType( DuMM::ChargedAtomTypeIndex(8000), Biotype::get("TIP3P", "Hydrogen"     ).getIndex() );
     if (!dumm.hasAtomClass(DuMM::AtomClassIndex(300)))
     {	
     	dumm.defineAtomClass_KA(
         	DuMM::AtomClassIndex(300),
         	"TIP3P Oxygen"      ,
         	16, //element number
         	2, //expected valence
         	1.7683,//from vdw parms
         	(.1520)
         	);
     }
     if (!dumm.hasChargedAtomType(DuMM::ChargedAtomTypeIndex(7000)))
     {			
     	dumm.defineChargedAtomType(
         	DuMM::ChargedAtomTypeIndex(7000),
         	"TIP3P Oxygen"      ,//magnesium",
         	DuMM::AtomClassIndex(300),
         	-.834                  //partial charge
         	);
     }
     if (! Biotype::exists("TIP3P"   ,"Oxygen"   ))
         Biotype::defineBiotype(Element::Oxygen(), 2, "TIP3P"   ,  "Oxygen"  ); // second arg is valence
                                                                                 // a residue name, second is an atom name
     dumm.setBiotypeChargedAtomType( DuMM::ChargedAtomTypeIndex(7000), Biotype::get("TIP3P" , "Oxygen"     ).getIndex() );
 

   
        setBaseAtom    (BivalentAtom("OW",Element::Oxygen()),Vec3(0));
        bondAtom    ( UnivalentAtom("HW0",Element::Hydrogen()),"OW/bond1",(.09572));
        bondAtom    ( UnivalentAtom("HW1",Element::Hydrogen()),"OW/bond2",(.09572));
        setDefaultBondAngle((104.52*Deg2Rad),"HW0","OW","HW1");

        setBiotypeIndex( "OW", Biotype::get("TIP3P"   ,"Oxygen").getIndex() );
        setBiotypeIndex("HW0", Biotype::get("TIP3P"   ,"Hydrogen").getIndex() );
        setBiotypeIndex("HW1", Biotype::get("TIP3P"   ,"Hydrogen").getIndex() );
        }   
};


}
