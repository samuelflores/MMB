

#include "MoleculeContainer.h"
#include <fstream>

#define MK_ELEMENT(name) \
    if (elementName.compare(#name) == 0) return Element::name()


void CompoundObjectMapContainer::loadCompoundMap() {
    MMBLOG_FILE_FUNC_LINE(INFO, "starting loadCompoundMap() "<<endl);
    compoundMap.insert (std::pair <const String , Compound> ("PurineBaseCore",PurineBaseCore() ));
    RibonucleotideResidue::Guanylate myGuanylate; myGuanylate.assignBiotypes();
    compoundMap.insert (std::pair <const String , Compound> ("Guanylate",myGuanylate ));
    //compoundMap.insert (std::pair <const String , Compound> ("Guanylate",RibonucleotideResidue::Guanylate() ));
    // from molmodel/include/molmodel/internal/NA.h:
    compoundMap.insert (std::pair <const String , Compound> ("NaPhosphodiesterLinkage",NaPhosphodiesterLinkage("NaPhosphodiesterLinkage") ));
    compoundMap.insert (std::pair <const String , Compound> ("FivePrimeNaPhosphateGroup",FivePrimeNaPhosphateGroup("FivePrimeNaPhosphateGroup") ));
    compoundMap.insert (std::pair <const String , Compound> ("ThreePrimeNaPhosphateGroup",ThreePrimeNaPhosphateGroup("ThreePrimeNaPhosphateGroup") ));
    compoundMap.insert (std::pair <const String , Compound> ("RibonucleosideResidue",RibonucleosideResidue("RibonucleosideResidue") ));
    compoundMap.insert (std::pair <const String , Compound> ("MethyleneGroup",MethyleneGroup() ));
    compoundMap.insert (std::pair <const String , Compound> ("methyl",MethylGroup() ));
    compoundMap.insert (std::pair <const String , Compound> ("MethylGroup",MethylGroup() ));
    compoundMap.insert (std::pair <const String , Compound> ("Methane",Methane()     ));
    compoundMap.insert (std::pair <const String , Compound> ("MagnesiumIon",MagnesiumIon()     ));
    compoundMap.insert (std::pair <const String , Compound> ("AromaticSixMemberedCHGroup",AromaticSixMemberedCHGroup() ));
    compoundMap.insert (std::pair <const String , Compound> ("AliphaticHydrogen",AliphaticHydrogen() ));

// From Compound.h:
// SingleAtom
    /*compoundMap.insert (std::pair <const String , Compound> ("UnivalentAtom",
    compoundMap.insert (std::pair <const String , Compound> ("BivalentAtom",
    compoundMap.insert (std::pair <const String , Compound> ("TrivalentAtom",
    compoundMap.insert (std::pair <const String , Compound> ("QuadrivalentAtom",*/
    compoundMap.insert (std::pair <const String , Compound> ("AlcoholOHGroup",AlcoholOHGroup()       ));
    compoundMap.insert (std::pair <const String , Compound> ("PrimaryAmineGroup",PrimaryAmineGroup() ));
    compoundMap.insert (std::pair <const String , Compound> ("CarboxylateGroup",CarboxylateGroup()   ));


    MMBLOG_FILE_FUNC_LINE(INFO, "done with loadCompoundMap() "<<endl);
}

void CompoundObjectMapContainer::printCompoundMap() {
    MMBLOG_FILE_FUNC_LINE(INFO, "Available compounds are: "<<endl);
    map <const String , Compound>::iterator compoundMapIterator = compoundMap.begin();
    for (compoundMapIterator = compoundMap.begin(); compoundMapIterator != compoundMap.end(); compoundMapIterator++) {
        std::cout<<compoundMapIterator->first<<std::endl;
    }
}

/*void CompoundObjectMapContainer::loadSingleAtomMap() {
    cout <<__FILE__<<":"<<__LINE__<<" starting loadSingleAtomMap() "<<endl;
    singleAtomMap.insert (std::pair <const String , Compound::SingleAtom> ("AliphaticHydrogen",AliphaticHydrogen() ));
    singleAtomMap.insert (std::pair <const String , Compound::SingleAtom> ("AliphaticCarbon",AliphaticCarbon() ));
    singleAtomMap.insert (std::pair <const String , Compound::SingleAtom> ("UnivalentAtom",UnivalentAtom() ));
    singleAtomMap.insert (std::pair <const String , Compound::SingleAtom> ("BivalentAtom",BivalentAtom() ));
    singleAtomMap.insert (std::pair <const String , Compound::SingleAtom> ("TrivalentAtom",TrivalentAtom() ));
    singleAtomMap.insert (std::pair <const String , Compound::SingleAtom> ("QuadrivalentAtom",QuadrivalentAtom() ));

    cout <<__FILE__<<":"<<__LINE__<<" done with loadSingleAtomMap() "<<endl;
}*/
void CompoundObjectMapContainer::loadBiotypeMap() {
    MMBLOG_FILE_FUNC_LINE(INFO, "starting loadBiotypeMap() "<<endl);
    biotypeMap.insert (std::pair <const String , Biotype> ("MethaneH",Biotype::MethaneH() ));
    biotypeMap.insert (std::pair <const String , Biotype> ("MethaneC",Biotype::MethaneC() ));
    biotypeMap.insert (std::pair <const String , Biotype> ("SerineN", Biotype::SerineN()  ));
    //biotypeMap.insert (std::pair <const String , Biotype> ("Magnesium Ion", Biotype::MagnesiumIon()  ));
    MMBLOG_FILE_FUNC_LINE(INFO, "done with loadBiotypeMap() "<<endl);
}

Compound CompoundObjectMapContainer::fetchCompound(const String compoundName)  {
    //printCompoundMap();
    if (compoundMap.find(compoundName) == compoundMap.end()) {
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "No Compound object found with name "<<compoundName <<endl);
    }
    else return compoundMap[compoundName];
}


Biotype CompoundObjectMapContainer::fetchBiotype   (const String compoundName)  {
    MMBLOG_FILE_FUNC_LINE(INFO, "You have requested biotype named : >"<<compoundName<<"< "<<endl);
    if (biotypeMap.find(compoundName) == biotypeMap.end()) {
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "No Biotype  object found with name >"<<compoundName <<"< "<<endl);
    }
    else {
        Biotype myBiotype = biotypeMap[compoundName];
        BiotypeIndex myBiotypeIndex = myBiotype.getIndex();
        if (!(myBiotype.exists(myBiotypeIndex))) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "This biotype doesn't appear to exist! >"<<compoundName <<"<, >"<<myBiotypeIndex<<endl);
        } 
        MMBLOG_FILE_FUNC_LINE(INFO, "Returning biotype with index : >"<<myBiotypeIndex<<"< "<<endl);
        return myBiotype;                
    }
}

const Element CompoundObjectMapContainer::fetchElement (String elementName) {
    MMBLOG_FILE_FUNC_LINE(INFO, "About to fetch element with name "<<elementName<<endl);
    MK_ELEMENT(Hydrogen);
    MK_ELEMENT(Deuterium);
    MK_ELEMENT(Helium);
    MK_ELEMENT(Lithium);
    MK_ELEMENT(Beryllium);
    MK_ELEMENT(Boron);
    MK_ELEMENT(Carbon);
    MK_ELEMENT(Nitrogen);
    MK_ELEMENT(Oxygen);
    MK_ELEMENT(Fluorine);
    MK_ELEMENT(Neon);
    MK_ELEMENT(Sodium);
    MK_ELEMENT(Magnesium);
    MK_ELEMENT(Aluminum);
    MK_ELEMENT(Silicon);
    MK_ELEMENT(Phosphorus);
    MK_ELEMENT(Sulfur);
    MK_ELEMENT(Chlorine);
    MK_ELEMENT(Argon);
    MK_ELEMENT(Potassium);
    MK_ELEMENT(Calcium);
    MK_ELEMENT(Scandium);
    MK_ELEMENT(Titanium);
    MK_ELEMENT(Vanadium);
    MK_ELEMENT(Chromium);
    MK_ELEMENT(Manganese);
    MK_ELEMENT(Iron);
    MK_ELEMENT(Cobalt);
    MK_ELEMENT(Nickel);
    MK_ELEMENT(Copper);
    MK_ELEMENT(Zinc);
    MK_ELEMENT(Gallium);
    MK_ELEMENT(Germanium);
    MK_ELEMENT(Arsenic);
    MK_ELEMENT(Selenium);
    MK_ELEMENT(Bromine);
    MK_ELEMENT(Krypton);
    MK_ELEMENT(Rubidium);
    MK_ELEMENT(Strontium);
    MK_ELEMENT(Yttrium);
    MK_ELEMENT(Zirconium);
    MK_ELEMENT(Niobium);
    MK_ELEMENT(Molybdenum);
    MK_ELEMENT(Technetium);
    MK_ELEMENT(Ruthenium);
    MK_ELEMENT(Rhodium);
    MK_ELEMENT(Palladium);
    MK_ELEMENT(Silver);
    MK_ELEMENT(Cadmium);
    MK_ELEMENT(Indium);
    MK_ELEMENT(Tin);
    MK_ELEMENT(Antimony);
    MK_ELEMENT(Tellurium);
    MK_ELEMENT(Iodine);
    MK_ELEMENT(Xenon);
    MK_ELEMENT(Cesium);
    MK_ELEMENT(Barium);
    MK_ELEMENT(Lanthanum);
    MK_ELEMENT(Cerium);
    MK_ELEMENT(Praseodymium);
    MK_ELEMENT(Neodymium);
    MK_ELEMENT(Promethium);
    MK_ELEMENT(Samarium);
    MK_ELEMENT(Europium);
    MK_ELEMENT(Gadolinium);
    MK_ELEMENT(Terbium);
    MK_ELEMENT(Dysprosium);
    MK_ELEMENT(Holmium);
    MK_ELEMENT(Erbium);
    MK_ELEMENT(Thulium);
    MK_ELEMENT(Ytterbium);
    MK_ELEMENT(Lutetium);
    MK_ELEMENT(Hafnium);
    MK_ELEMENT(Tantalum);
    MK_ELEMENT(Tungsten);
    MK_ELEMENT(Rhenium);
    MK_ELEMENT(Osmium);
    MK_ELEMENT(Iridium);
    MK_ELEMENT(Platinum);
    MK_ELEMENT(Gold);
    MK_ELEMENT(Mercury);
    MK_ELEMENT(Thallium);
    MK_ELEMENT(Lead);
    MK_ELEMENT(Bismuth);
    MK_ELEMENT(Polonium);
    MK_ELEMENT(Astatine);
    MK_ELEMENT(Radon);
    MK_ELEMENT(Francium);
    MK_ELEMENT(Radium);
    MK_ELEMENT(Actinium);
    MK_ELEMENT(Thorium);
    MK_ELEMENT(Protactinium);
    MK_ELEMENT(Uranium);
    MK_ELEMENT(Neptunium);
    MK_ELEMENT(Plutonium);
    MK_ELEMENT(Americium);
    MK_ELEMENT(Curium);
    MK_ELEMENT(Berkelium);
    MK_ELEMENT(Californium);
    MK_ELEMENT(Einsteinium);
    MK_ELEMENT(Fermium);
    MK_ELEMENT(Mendelevium);
    MK_ELEMENT(Nobelium);
    MK_ELEMENT(Lawrencium);
    MK_ELEMENT(Rutherfordium);
    MK_ELEMENT(Dubnium);
    MK_ELEMENT(Seaborgium);
    MK_ELEMENT(Bohrium);
    MK_ELEMENT(Hassium);
    MK_ELEMENT(Meitnerium);
    MK_ELEMENT(Darmstadtium);
    MK_ELEMENT(Roentgenium);
    MK_ELEMENT(Ununbium);
    MK_ELEMENT(Ununtrium);
    MK_ELEMENT(Ununquadium);
    MK_ELEMENT(Ununpentium);
    MK_ELEMENT(Ununhexium);
    // If we didn't return after the above, something is wrong:
    MMBLOG_FILE_FUNC_LINE(CRITICAL, "Element "<<elementName<<" not found! Did you use proper capitalization (e.g. Carbon)?"<<endl);
    //return Element::Hydrogen();
}




Compound::SingleAtom  CompoundObjectMapContainer::fetchSingleAtom(const String className, Compound::AtomName& atomName , String elementName, Angle angle1 = 180*Deg2Rad )  {
    MMBLOG_FILE_FUNC_LINE(INFO, "Fetching single atom of class >"<<className<<"< "<<endl);
    const Element myElement = fetchElement(elementName);
    MMBLOG_FILE_FUNC_LINE(INFO, "using element : >"<<myElement.getName()<<"< "<<endl);
    MMBLOG_FILE_FUNC_LINE(INFO, endl);
    if (className.compare("AliphaticHydrogen") == 0) {return AliphaticHydrogen(atomName);}
    if (className.compare("AliphaticCarbon") == 0) {return AliphaticCarbon(atomName);}
    if (className.compare("UnivalentAtom") == 0) {
        MMBLOG_FILE_FUNC_LINE(INFO, "using element : >"<<myElement.getName()<<"< "<<endl);
        return UnivalentAtom(atomName,myElement);
        MMBLOG_FILE_FUNC_LINE(INFO, endl);
    }
    if (className.compare("BivalentAtom") == 0) {return BivalentAtom(atomName,fetchElement(elementName) , angle1);} // the angle should be adjusted later
    if (className.compare("TrivalentAtom") == 0) {return TrivalentAtom(atomName,fetchElement(elementName));}
    if (className.compare("QuadrivalentAtom") == 0) {return QuadrivalentAtom(atomName,fetchElement(elementName));}
    // If we didn't return after the above, something is wrong:
    MMBLOG_FILE_FUNC_LINE(CRITICAL, "The SingleAtom subclass "<<className<<" was not found! Accepted subclasses are UnivalentAtom, BivalentAtom, BivalentAtom, TrivalentAtom, QuadrivalentAtom"<<endl);
}

CustomMolecule::CustomMolecule(vector <vector <String> > moleculeBuildCommandVector, DuMMForceFieldSubsystem & dumm) {
    CompoundObjectMapContainer compoundObjectMapContainer;
    if (moleculeBuildCommandVector.size() < 1) {
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "setBaseCompound must be the first command you issue!"<<endl);
    }
    for (size_t i = 0; i < moleculeBuildCommandVector.size(); i++) {
        MMBLOG_FILE_FUNC_LINE(INFO, "Contents of moleculeBuildCommandVector[i] : ("<<moleculeBuildCommandVector[i].size()<<" elements) " <<endl);
        for (size_t j = 0; j < moleculeBuildCommandVector[i].size(); j++) {
            MMBLOG_FILE_FUNC_LINE(INFO, ">"<<moleculeBuildCommandVector[i][j]<<"< ");
        }
        MMBLOG_FILE_FUNC_LINE(INFO, endl);
	if ((moleculeBuildCommandVector[i])[0].compare("setBaseCompound") == 0) {  // element 0 of moleculeBuildCommand is always a command
            //MMBLOG_FILE_FUNC_LINE(" Contents of moleculeBuildCommandVector[i] : >"<<moleculeBuildCommandVector[i][0]<<"<, >"<<moleculeBuildCommandVector[i][1]<<"< ."<<endl; ;
            compoundObjectMapContainer.printCompoundMap();
	    if (i > 0) {MMBLOG_FILE_FUNC_LINE(CRITICAL, "setBaseCompound must be the first command you issue!"<<endl); }
	    if ((moleculeBuildCommandVector[i]).size() != 2) {MMBLOG_FILE_FUNC_LINE(CRITICAL, "Wrong number of parameters!"<<endl);}
            MMBLOG_FILE_FUNC_LINE(INFO, "About to fetchCompound ("<<moleculeBuildCommandVector[i][1]<<") "<<endl);
	    Compound myCompound = compoundObjectMapContainer.fetchCompound(moleculeBuildCommandVector[i][1]);
	    setBaseCompound(moleculeBuildCommandVector[i][1],                                    // name of the Compound (e.g. "MethyleneGroup") is element [1] of moleculeBuildCommand
		myCompound               //followed by the corresponding Compound object
		); 
	    inheritAtomNames(moleculeBuildCommandVector[i][1]);                                  // for now I am assuming we will always want to inherit atom names of this compound.
	}
	else if ((moleculeBuildCommandVector[i])[0].compare("setBaseAtom") == 0) {  // element 0 of moleculeBuildCommand is always a command
        MMBLOG_FILE_FUNC_LINE(INFO, "Syntax: setBaseAtom <object type> <atom name> <element name> [<angle1>] "<<endl);
	    if (i > 0) {MMBLOG_FILE_FUNC_LINE(CRITICAL, "If used, this must be the first command you issue!"<<endl);}
	    if ((moleculeBuildCommandVector[i]).size() > 5) {MMBLOG_FILE_FUNC_LINE(CRITICAL, "Too many parameters!"<<endl);}
	    String objectType = moleculeBuildCommandVector[i][1]; 
	    String atomName = moleculeBuildCommandVector[i][2]; 
	    String elementname = moleculeBuildCommandVector[i][3]; 
            Angle myAngle1 = 180*Deg2Rad;
            if (moleculeBuildCommandVector[i].size() > 4)
                myAngle1 = atof(moleculeBuildCommandVector[i][4].c_str());

            MMBLOG_FILE_FUNC_LINE(INFO, endl);
	    setBaseAtom    (
		compoundObjectMapContainer.fetchSingleAtom(objectType, atomName, elementname , myAngle1)               //followed by the corresponding Compound object
		);
        MMBLOG_FILE_FUNC_LINE(INFO, endl);
	    //inheritAtomNames(moleculeBuildCommandVector[i][1]);                                  // for now I am assuming we will always want to inherit atom names of this compound.
	}

	else if ((moleculeBuildCommandVector[i])[0].compare("convertInboardBondCenterToOutboard") == 0 ) {
            //MMBLOG_FILE_FUNC_LINE(" Contents of moleculeBuildCommandVector[i] : >"<<moleculeBuildCommandVector[i][0]<<"<"<<endl;
	    if (i == 0) {MMBLOG_FILE_FUNC_LINE(CRITICAL, "This command cannot be the first to issue.  You first need to setBaseCompound!"<<endl);}
	    convertInboardBondCenterToOutboard();}

    else if ((moleculeBuildCommandVector[i])[0].compare("bondCompound") == 0 ) {
        MMBLOG_FILE_FUNC_LINE(INFO, "Syntax : bondCompound <name of added compound in parent> <compound to add> <name of parent bond at which to attach e.g. methyl/bond >  "<<endl);
            /*MMBLOG_FILE_FUNC_LINE(" Contents of moleculeBuildCommandVector[i] : "<<endl;;
            for (int j = 0; j < moleculeBuildCommandVector[i].size(); j++) {
                cout<<">"<<moleculeBuildCommandVector[i][j]<<"< ";
            }
            cout<<endl;*/
        if (moleculeBuildCommandVector[i].size  () >4){
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "Too many parameters!"<<endl);
        }
        String subcompoundNameInParent = moleculeBuildCommandVector[i][1];
        String compoundToAdd = moleculeBuildCommandVector[i][2];
        String bondName = moleculeBuildCommandVector[i][3];
        Compound myCompound = compoundObjectMapContainer.fetchCompound(compoundToAdd);
	    bondCompound(
            subcompoundNameInParent,
            myCompound,
            bondName,  // name of bond at which to attach the atom    
            .14 //hack, fix later
        );
	    // For example:
            // bondCompound("methyl2", MethylGroup(), "methyl1/bond");
            
    }
    else if ((moleculeBuildCommandVector[i])[0].compare("bondAtom") == 0 ) {
        MMBLOG_FILE_FUNC_LINE(INFO, "bondAtom <molmodel class of added atom (AliphaticHydrogen, UnivalentAtom, DivalentAtom, etc> <name of atom to be added, e.g. H1> <name of element to be added e.g. Hydrogen > <name of bond at which to attach atom e.g. methyl/bond > <bond length> [<dihedral angle, degrees>] [bond mobility: Default, Free, Torsion, or Rigid]"<<endl);
        if (moleculeBuildCommandVector[i].size  () <6){
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "Too few parameters!"<<endl);
        }
        if (moleculeBuildCommandVector[i].size  () >8){
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "Too many parameters!"<<endl);
        }
        String singleAtomSubclass = moleculeBuildCommandVector[i][1];
        String addedAtomName = moleculeBuildCommandVector[i][2];
        String addedElementName = moleculeBuildCommandVector[i][3]; 
        String bondName = moleculeBuildCommandVector[i][4];
        double bondLength = atof(moleculeBuildCommandVector[i][5].c_str()) ;
	    
        Angle myAngle1 = 180*Deg2Rad; // this is a bond angle 
        //if (moleculeBuildCommandVector[i].size() > 6)
        //    myAngle1 = atof(moleculeBuildCommandVector[i][6].c_str());
        Compound::SingleAtom addedAtom = compoundObjectMapContainer.fetchSingleAtom(singleAtomSubclass,addedAtomName,addedElementName );

        Angle myDihedral = 180*Deg2Rad;
        if (moleculeBuildCommandVector[i].size() > 6)
            myDihedral = atof(moleculeBuildCommandVector[i][6].c_str())*Deg2Rad ;

        String myBondMobilityString = "Torsion";
        if (moleculeBuildCommandVector[i].size() > 7) {
            myBondMobilityString = moleculeBuildCommandVector[i][7];}
        BondMobility::Mobility myBondMobility = stringToBondMobility(myBondMobilityString);
            
        bondAtom(
            addedAtom,
            bondName,  // name of bond at which to attach the atom    
            bondLength,
            myDihedral,          // This is the default dihedral angle parameter, defaults to 180 * Deg2Rad (i.e. pi)
            myBondMobility
        ); // Default bond length
	    // For example:
	    // bondAtom(AliphaticHydrogen("H4"), "methyl/bond", 0.1112);
	}

	else if ((moleculeBuildCommandVector[i])[0].compare("setBiotypeIndex") == 0 ) {
            MMBLOG_FILE_FUNC_LINE(INFO, "setBiotypeIndex <specific atom name e.g. H4> <Biotype e.g. MethaneH>"<<endl);
            MMBLOG_FILE_FUNC_LINE(INFO, "     or: setBiotypeIndex <specific atom name e.g. OP4> <Residue e.g. Phosphate,?RNA> <generic atom name e.g. OP> <ordinality : Initial, Any, or Final> "<<endl);
            String specificAtomName = moleculeBuildCommandVector[i][1];
            if ((moleculeBuildCommandVector[i]).size() == 3) {
                String biotypeName = moleculeBuildCommandVector[i][2];
                for (int k = 0; k < biotypeName.length(); k++) {
                    if (String(biotypeName[k]).compare("?") ==0) {
                        biotypeName[k] = (String(" "))[0];
                    }
                }
		setBiotypeIndex(specificAtomName,               //name
		    compoundObjectMapContainer.fetchBiotype(biotypeName).getIndex() //moleculeBuildCommandVector[i][2]).getIndex() // BiotypeIndex
		    );
		// For example:
		// setBiotypeIndex( "H4", Biotype::MethaneH().getIndex() );
            } else if ((moleculeBuildCommandVector[i]).size() == 5) {
		String specificAtomName = moleculeBuildCommandVector[i][1];
		String residueName      = moleculeBuildCommandVector[i][2];
                for (int k = 0; k < residueName.length(); k++) {
                    if (String(residueName[k]).compare("?") ==0) {
                        residueName[k] = (String(" "))[0];
                    }
                }
		String genericAtomName  = moleculeBuildCommandVector[i][3];
                MMBLOG_FILE_FUNC_LINE(INFO, "specificAtomName is now : >"<<specificAtomName<<"< "<<endl);
                MMBLOG_FILE_FUNC_LINE(INFO, "Residue name is now : >"<<residueName<<"< "<<endl);
                MMBLOG_FILE_FUNC_LINE(INFO, "genericAtomName is now : >"<<genericAtomName<<"< "<<endl);
                String ordinalityString = moleculeBuildCommandVector[i][4];
                enum  	ordinalityEnum { Any = 1, Initial = 2, Final = 3 };
                Ordinality::Residue myOrdinality;
                if (ordinalityString.compare("Initial") == 0) myOrdinality =  SimTK::Ordinality::Initial;
                else if (ordinalityString.compare("Final") == 0) myOrdinality =  SimTK::Ordinality::Final;
                else if (ordinalityString.compare("Any") == 0) myOrdinality =  SimTK::Ordinality::Any;
                else {
		    MMBLOG_FILE_FUNC_LINE(CRITICAL, "Invalid ordinality : >"<<ordinalityString<<"< "<<endl);
                }
		BiotypeIndex myBiotypeIndex =  Biotype::get(residueName, genericAtomName, 
                        myOrdinality
                        ).getIndex();                
                MMBLOG_FILE_FUNC_LINE(INFO, "BiotypeIndex = "<<myBiotypeIndex<<endl);
                setBiotypeIndex(specificAtomName,myBiotypeIndex);
            } else {
                MMBLOG_FILE_FUNC_LINE(CRITICAL, "Wrong number of parameters ("<<(moleculeBuildCommandVector[i]).size()<<")!"<<endl);

            }
	}
	else if ((moleculeBuildCommandVector[i])[0].compare("defineBiotype") == 0 ) {
            MMBLOG_FILE_FUNC_LINE(INFO, "Syntax: defineBiotype <element symbol (e.g. O, H, C)> <valence (integer)> <residue name> <generic atom name, e.g. Oxygen>"<<endl);
            if ((moleculeBuildCommandVector[i]).size() != 5) {
                MMBLOG_FILE_FUNC_LINE(CRITICAL, "Wrong number of parameters! "<<endl);
                }
            Biotype::defineBiotype (
                Element::getBySymbol(moleculeBuildCommandVector[i][1]),
		atoi(moleculeBuildCommandVector[i][2].c_str()),
                moleculeBuildCommandVector[i][3] ,
                moleculeBuildCommandVector[i][4]              
            );

        }

	else if ((moleculeBuildCommandVector[i])[0].compare("setCompoundName") == 0 ) {
            //MMBLOG_FILE_FUNC_LINE(" Contents of moleculeBuildCommandVector[i] : >"<<moleculeBuildCommandVector[i][0]<<"<, >"<<moleculeBuildCommandVector[i][1]<<"<, >"<<moleculeBuildCommandVector[i][2]<<"< ."<<endl;
	    if (moleculeBuildCommandVector[i][1].length() == 0 ) {MMBLOG_FILE_FUNC_LINE(CRITICAL, "Length of Compound name must be > 0"<<endl); }
        setCompoundName(moleculeBuildCommandVector[i][1]);}
	else if ((moleculeBuildCommandVector[i])[0].compare("defineAndSetChargedAtomType") == 0 ) {

            MMBLOG_FILE_FUNC_LINE(INFO, "Syntax : defineAndSetChargedAtomType <biotype name> <FF atom class index> <charge> "<<endl);
            MMBLOG_FILE_FUNC_LINE(INFO, "     or : defineAndSetChargedAtomType <residue e.g. Phosphate,?RNA> <generic atom name e.g. OP> <ordinality : Initialy, Any, or Final> <FF atom class index> <charge> "<<endl);

            //MMBLOG_FILE_FUNC_LINE(" Contents of moleculeBuildCommandVector[i] : >"<<moleculeBuildCommandVector[i][0]<<"<, >"<<moleculeBuildCommandVector[i][1]<<"<, >"<<moleculeBuildCommandVector[i][2]<<"< ."<<endl;

            if ((moleculeBuildCommandVector[i]).size() == 4) {
                MMBLOG_FILE_FUNC_LINE(INFO, endl);
		DuMM::ChargedAtomTypeIndex 	myChargedAtomTypeIndex = dumm.getNextUnusedChargedAtomTypeIndex (); 
		String myBiotypeName = moleculeBuildCommandVector[i][1];
		dumm.defineChargedAtomType(myChargedAtomTypeIndex, 
		    myBiotypeName,           // biotype name. This is actually not used.
		    DuMM::AtomClassIndex(atoi(moleculeBuildCommandVector[i][2].c_str())), // force field atom class index
		    atof(moleculeBuildCommandVector[i][3].c_str()));	              // charge
		//if (
		dumm.setBiotypeChargedAtomType(myChargedAtomTypeIndex, compoundObjectMapContainer.fetchBiotype(myBiotypeName).getIndex());
		// For example:
		//                      index just needs to be unused.    doesn't matter.Class in force field.     charge
		//defineChargedAtomType(DuMM::ChargedAtomTypeIndex(5000), "Methane C",   DuMM::AtomClassIndex(1),  0.04);
		//setBiotypeChargedAtomType(DuMM::ChargedAtomTypeIndex(5000), Biotype::MethaneC().getIndex());

            } else if ((moleculeBuildCommandVector[i]).size() == 6) {
                MMBLOG_FILE_FUNC_LINE(INFO, endl);
		DuMM::ChargedAtomTypeIndex 	myChargedAtomTypeIndex = dumm.getNextUnusedChargedAtomTypeIndex (); 
		//#String myBiotypeName = moleculeBuildCommandVector[i][1];
		dumm.defineChargedAtomType(myChargedAtomTypeIndex, 
		    String("myBiotypeName"),           // biotype name. This is actually not used.
		    DuMM::AtomClassIndex(atoi(moleculeBuildCommandVector[i][4].c_str())), // force field atom class index
		    atof(moleculeBuildCommandVector[i][5].c_str()));	              // charge
		String residueName      = moleculeBuildCommandVector[i][1];
                for (int k = 0; k < residueName.length(); k++) {
                    if (String(residueName[k]).compare("?") ==0) {
                        residueName[k] = (String(" "))[0];
                    }
                }
                MMBLOG_FILE_FUNC_LINE(INFO, "Residue name is now : >"<<residueName<<"< "<<endl);
		String genericAtomName  = moleculeBuildCommandVector[i][2];
                String ordinalityString = moleculeBuildCommandVector[i][3];
                enum  	ordinalityEnum { Any = 1, Initial = 2, Final = 3 };
                Ordinality::Residue myOrdinality;
                if (ordinalityString.compare("Initial") == 0) myOrdinality =  SimTK::Ordinality::Initial;
                else if (ordinalityString.compare("Final") == 0) myOrdinality =  SimTK::Ordinality::Final;
                else if (ordinalityString.compare("Any") == 0) myOrdinality =  SimTK::Ordinality::Any;
                else {
                    MMBLOG_FILE_FUNC_LINE(CRITICAL, "Invalid ordinality : >"<<ordinalityString<<"< "<<endl);
                }
                dumm.setBiotypeChargedAtomType (myChargedAtomTypeIndex, 
		    Biotype::get(residueName, genericAtomName, 
                        myOrdinality
                        ).getIndex()                
                    );

            } else {
                MMBLOG_FILE_FUNC_LINE(CRITICAL, "Wrong number of parameters ("<<(moleculeBuildCommandVector[i]).size()<<")!"<<endl);
            }

        }
	else if ((moleculeBuildCommandVector[i])[0].compare("setDefaultBondAngle") == 0 ) {
            MMBLOG_FILE_FUNC_LINE(INFO, "Syntax : setDefaultBondAngle <Angle (in degrees)> <atom name 1> <atom name 2 (central atom)> <atom name 3>"<<endl);
            MMBLOG_FILE_FUNC_LINE(INFO, "example: setDefaultBondAngle 104.52 HW1 OW HW2 "<<endl);
            double myAngle = atof(moleculeBuildCommandVector[i][1].c_str());
            String atomName1 = moleculeBuildCommandVector[i][2]; 
            String atomName2 = moleculeBuildCommandVector[i][3]; 
            String atomName3 = moleculeBuildCommandVector[i][4]; 
            setDefaultBondAngle(myAngle*Deg2Rad,atomName1,atomName2,atomName3);
        }
	else if ((moleculeBuildCommandVector[i])[0].compare("setBiotypeChargedAtomType") == 0 ) {
                MMBLOG_FILE_FUNC_LINE(INFO, "Syntax: setBiotypeChargedAtomType <residue name> <generic atom name> <ordinality> <FF atom class index> <charge>"<<endl);
                MMBLOG_FILE_FUNC_LINE(INFO, "Or    : setBiotypeChargedAtomType <(DuMM) chargedAtomTypeIndex> <biotypeIndex>"<<endl);
                DuMM::ChargedAtomTypeIndex      myChargedAtomTypeIndex;
                MMBLOG_FILE_FUNC_LINE(INFO, "myChargedAtomTypeIndex = >"<<myChargedAtomTypeIndex<<"< "<<endl);
                BiotypeIndex myBiotypeIndex;
                if ((moleculeBuildCommandVector[i]).size() == 6) {
			String residueName      = moleculeBuildCommandVector[i][1];
			for (int k = 0; k < residueName.length(); k++) {
			    if (String(residueName[k]).compare("?") ==0) {
				residueName[k] = (String(" "))[0];
			    } // of if 
			} // of for k
			MMBLOG_FILE_FUNC_LINE(INFO, "Residue name is now : >"<<residueName<<"< "<<endl);
			String genericAtomName  = moleculeBuildCommandVector[i][2];
			String ordinalityString = moleculeBuildCommandVector[i][3];
			DuMM::AtomClassIndex myAtomClassIndex ( atoi(moleculeBuildCommandVector[i][4].c_str())); // force field atom class index
			double charge = atof(moleculeBuildCommandVector[i][5].c_str());// charge
			enum  	ordinalityEnum { Any = 1, Initial = 2, Final = 3 };
			Ordinality::Residue myOrdinality;
			if (ordinalityString.compare("Initial") == 0) myOrdinality =  SimTK::Ordinality::Initial;
			else if (ordinalityString.compare("Final") == 0) myOrdinality =  SimTK::Ordinality::Final;
			else if (ordinalityString.compare("Any") == 0) myOrdinality =  SimTK::Ordinality::Any;
			else {
                MMBLOG_FILE_FUNC_LINE(CRITICAL, "Invalid ordinality : >"<<ordinalityString<<"< "<<endl);
			}
			myBiotypeIndex = Biotype::get(residueName, genericAtomName, 
				myOrdinality
				).getIndex();                
                        MMBLOG_FILE_FUNC_LINE(INFO, "myChargedAtomTypeIndex = >"<<myChargedAtomTypeIndex<<"< "<<endl);
			myChargedAtomTypeIndex = dumm.getNextUnusedChargedAtomTypeIndex (); 
                        MMBLOG_FILE_FUNC_LINE(INFO, "myChargedAtomTypeIndex = >"<<myChargedAtomTypeIndex<<"< "<<endl);
			String myBiotypeName = "blah";
		   
			dumm.defineChargedAtomType(myChargedAtomTypeIndex, 
				myBiotypeName,         // biotype name. This is actually not used.
				myAtomClassIndex,      // force field atom class index
			charge);	       // charge
                        dumm.setBiotypeChargedAtomType(myChargedAtomTypeIndex, myBiotypeIndex);
                } else if ((moleculeBuildCommandVector[i]).size() == 3) { MMBLOG_FILE_FUNC_LINE(INFO, endl);
                        if (dumm.hasChargedAtomType ((moleculeBuildCommandVector[i])[1])) {
                            MMBLOG_FILE_FUNC_LINE(INFO, "Found chargedAtomType for >"<<(moleculeBuildCommandVector[i])[1]<<"< "<<endl);
                        }
                        //myChargedAtomTypeIndex= atoi((moleculeBuildCommandVector[i])[1].c_str());
                        //myBiotypeIndex= atoi((moleculeBuildCommandVector[i])[2].c_str());
                        dumm.setBiotypeChargedAtomType(myChargedAtomTypeIndex, myBiotypeIndex);
                } else {
                    MMBLOG_FILE_FUNC_LINE(CRITICAL, "Incorrect number of parameters!"<<endl); }
        }
        else {
            /*MMBLOG_FILE_FUNC_LINE(" Contents of moleculeBuildCommandVector[i] : "<<endl;;
            for (int j = 0; j < moleculeBuildCommandVector[i].size(); j++) {
                cout<<">"<<moleculeBuildCommandVector[i][j]<<"< ";
            }
            cout<<endl;*/
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "Unknown command!"<<endl);
        }
    }
};

/*void MoleculeClass::setChainID(String chainID) {
    molecule.setPdbChainId(chainID);
}*/

void  MoleculeClass::setPdbResidueName() {
    if (residueName.length() > 3) {
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "residueName >"<<residueName <<"< is too long!" <<endl);
    }
    MMBLOG_FILE_FUNC_LINE(INFO, "setting residueName >"<<residueName <<"<"<<endl);
    molecule.setPdbResidueName(residueName);
};

void MoleculeClass::includeAllAtoms( DuMMForceFieldSubsystem & dumm) {
    for (Compound::AtomIndex i  = Compound::AtomIndex(0) ; i <molecule.getNumAtoms(); i++) {
	dumm.includeNonbondAtom(molecule.getDuMMAtomIndex(i));
    } 
}
void MoleculeClass::addRingClosingBond(CovalentBondClass myBond){
        MMBLOG_FILE_FUNC_LINE(INFO, "about start MoleculeClass::addRingClosingBond. molecule.getNumAtoms() = "<<molecule.getNumAtoms()<<endl);
        if (molecule.getNumAtoms() == 0) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "Molecule has  "<<molecule.getNumAtoms() <<" atoms! " <<endl);
        }
        const Compound::BondCenterPathName & centerName1 = /*String ("1/")+*/String(myBond.getAtomName1())+String('/')+String(myBond.getBondCenterName1());
        const Compound::BondCenterPathName & centerName2 = String(myBond.getAtomName2())+String('/')+String(myBond.getBondCenterName2());
        if (!(molecule.hasBondCenter(centerName1))){
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "Unable to find bond center "<<centerName1<<std::endl);
    }
        if (!(molecule.hasBondCenter(centerName2))){
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "Unable to find bond center "<<centerName2<<std::endl);
    }
    double bondLength = 111.111 ; // This doesn't matter, so I set to an absurd value
    double dihedralAngle = 0.0; // This doesn't matter either.
    
    molecule.addRingClosingBond( centerName1,    centerName2 , bondLength, dihedralAngle, SimTK::BondMobility::Rigid);
}
void MoleculeClassContainer::add(String myChainID,  String myResidueName, MoleculeClass & myMoleculeClass) {
    //myMoleculeClass.setChainID(chainID);
    myMoleculeClass.setChainID(myChainID);
    myMoleculeClass.setResidueName(myResidueName);
    if (moleculeClassMap.count(myChainID) > 0 ) {
        MMBLOG_FILE_FUNC_LINE(CRITICAL, moleculeClassMap.count(myChainID) <<" MoleculeClass's found with chain ID "<<myChainID <<endl);
    }
    moleculeClassMap.insert(std::pair <const String , MoleculeClass>  (myChainID,myMoleculeClass));
    //MMBLOG_FILE_FUNC_LINE(" Just inserted MoleculeClass with chain ID >"<<myChainID<<"< "<<endl;
    //MMBLOG_FILE_FUNC_LINE(" The inserted MoleculeClass has chain ID >"<<myMoleculeClass.getChainID()<<"< "<<endl;
    //MMBLOG_FILE_FUNC_LINE(" The retrieved MoleculeClass has chain ID >"<<updMoleculeClass(myChainID).getChainID()<<"< "<<endl;
    
}

void MoleculeClassContainer::adoptCompounds(SimTK::CompoundSystem & mySystem){
    map<const String, MoleculeClass>::iterator it;
    map<const String, MoleculeClass>::iterator next;
    //int i = 0;
    next = moleculeClassMap.begin();
    while (next != moleculeClassMap.end())
    {
       it = next;
       
       MMBLOG_FILE_FUNC_LINE(INFO, "About to adopt CustomMolecule with chain ID >"<<(it->second).getChainID()<<"<"<<endl);
       //MMBLOG_FILE_FUNC_LINE(" Top level transform before adopting: "<<(it->second.molecule).getTopLevelTransform()<<endl;
       mySystem.adoptCompound(it->second.molecule);
       //MMBLOG_FILE_FUNC_LINE(" Top level transform after adopting: "<<(it->second.molecule).getTopLevelTransform()<<endl;
       next++;
    }
};
void MoleculeClassContainer::initializeCompounds(DuMMForceFieldSubsystem & dumm){
    map<const String, MoleculeClass>::iterator it;
    map<const String, MoleculeClass>::iterator next;
    //int i = 0;
    next = moleculeClassMap.begin();
    while (next != moleculeClassMap.end())
    {
       it = next;
       (it->second).molecule = CustomMolecule((it->second).moleculeBuildCommandVector, dumm);
       //MMBLOG_FILE_FUNC_LINE(" About to set chain ID to >"<<(it->second).getChainID()<<"<"<<endl;
       (it->second).setPdbChainID(); // transfer MoleculeClass.chainID to the molecule member's PDB chain ID
       (it->second).setPdbResidueName(); // ditto for residue Name.
       (it->second).molecule.setPdbResidueNumber(1); // ditto for residue Name.
       next++;
    }
};

MoleculeClass &   MoleculeClassContainer::updMoleculeClass(String myChainID) {
    validateChainID(myChainID);
    map<const String, MoleculeClass>::iterator it;
    //if (moleculeClassMap.count(myChainID) == 1) {
    it = moleculeClassMap.find(myChainID);
    return it->second;
};

void  MoleculeClassContainer::validateChainID(String myChainID){
    if (moleculeClassMap.count(myChainID) != 1) {
        MMBLOG_FILE_FUNC_LINE(CRITICAL, moleculeClassMap.count(myChainID) <<" MoleculeClass's found with chain ID "<<myChainID <<endl);
    }
}

void MoleculeClassContainer::matchDefaultConfiguration(bool readPreviousFrameFile, String pdbFileName,bool matchExact, bool matchIdealized)
{
    MMBLOG_FILE_FUNC_LINE(INFO, "readPreviousFrameFile = "<<readPreviousFrameFile<<", pdbFileName = >"<<pdbFileName<<"< "<<endl);
    if (readPreviousFrameFile)
    {
        PdbStructure pdbStructure;
        
        //============================================ Read in PDB or CIF
        if ( pdbFileName.length() > 4 )
        {
            if ( pdbFileName.substr ( pdbFileName.length() - 4, pdbFileName.length() - 1) == ".pdb" )
            {
                std::ifstream inputFile               ( pdbFileName.c_str(), ifstream::in );
                
                if ( !inputFile.good() )
                {
                    MMBLOG_FILE_FUNC_LINE(CRITICAL, "The file " << pdbFileName << " could not be opened. If this is not the file you wanted to open, please supply the requested file name after the loadSequencesFromPdb command. Note that the supported file extensions currently are \".pdb\", \".cif\" and \".cif.gz\"." << endl);
                }
                else
                {
                    pdbStructure                      = PdbStructure (inputFile);
                }
            }
            else if ( pdbFileName.substr ( pdbFileName.length() - 4, pdbFileName.length() - 1) == ".cif" )
            {
                std::ifstream testOpen                ( pdbFileName.c_str() );
                if ( testOpen.good() )
                {
                    pdbStructure                      = PdbStructure (pdbFileName);
                }
                else
                {
                    std::string pdbFileHlp            = pdbFileName;
                    pdbFileHlp.append                 ( ".gz" );
                    
                    std::ifstream testOpen2           ( pdbFileHlp.c_str() );
                    if ( testOpen2.good() )
                    {
                        pdbStructure                  = PdbStructure ( pdbFileHlp );
                    }
                    else
                    {
                        MMBLOG_FILE_FUNC_LINE(CRITICAL, "The file " << pdbFileName << " could not be opened. If this is not the file you wanted to open, please supply the requested file name after the loadSequencesFromPdb command. Note that the supported file extensions currently are \".pdb\", \".cif\" and \".cif.gz\"." << endl);
                    }
                    testOpen2.close                   ( );
                }
                testOpen.close                        ( );
            }
            else if ( pdbFileName.length() > 7 )
            {
                if ( pdbFileName.substr ( pdbFileName.length() - 7, pdbFileName.length() - 1) == ".cif.gz" )
                {
                pdbStructure                          = PdbStructure (pdbFileName);
                }
                else
                {
                    MMBLOG_FILE_FUNC_LINE(CRITICAL, "The file " << pdbFileName << " could not be opened. If this is not the file you wanted to open, please supply the requested file name after the loadSequencesFromPdb command. Note that the supported file extensions currently are \".pdb\", \".cif\" and \".cif.gz\"." << endl);
                }
            }
            else
            {
                MMBLOG_FILE_FUNC_LINE(CRITICAL, "The file " << pdbFileName << " could not be opened. If this is not the file you wanted to open, please supply the requested file name after the loadSequencesFromPdb command. Note that the supported file extensions currently are \".pdb\", \".cif\" and \".cif.gz\"." << endl);
            }
        }
        else
        {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "The file " << pdbFileName << " could not be opened. If this is not the file you wanted to open, please supply the requested file name after the loadSequencesFromPdb command. Note that the supported file extensions currently are \".pdb\", \".cif\" and \".cif.gz\"." << endl);
        }
        
	    map<const String, MoleculeClass>::iterator it;
	    map<const String, MoleculeClass>::iterator next;
	    //int i = 0;
	    next = moleculeClassMap.begin();
	    while (next != moleculeClassMap.end())
	    {
	       it = next;
	       Compound & myCompound = (it->second).molecule;
            MMBLOG_FILE_FUNC_LINE(INFO, "About to create atom targets from file "<<pdbFileName<<endl);
            MMBLOG_FILE_FUNC_LINE(INFO, myCompound.getPdbChainId()<<endl);
	       Compound::AtomTargetLocations atomTargets = myCompound.createAtomTargets(pdbStructure);
	       map<Compound::AtomIndex, Vec3>::iterator targetIt;
	       map<Compound::AtomIndex, Vec3>::iterator targetNext;
	       targetNext = atomTargets.begin();
	       while (targetNext != atomTargets.end())
	       {
	          targetIt = targetNext;
                MMBLOG_FILE_FUNC_LINE(INFO, endl);
                MMBLOG_FILE_FUNC_LINE(INFO, targetIt->second<<endl);
                  targetNext++;
               }
	       if (matchExact)
	    	    {
                        //MMBLOG_FILE_FUNC_LINE(" Top level transform before fitting:     "<<myCompound.getTopLevelTransform()<<endl;
	    	    myCompound.matchDefaultConfiguration(atomTargets,   Compound::Match_Exact );
                        //MMBLOG_FILE_FUNC_LINE(" Top level transform after fitting:     "<<myCompound.getTopLevelTransform()<<endl;
	    	    }
	       if (matchIdealized)
	    	    {myCompound.matchDefaultConfiguration(atomTargets,   Compound::Match_Idealized );} //planarity     tolerance is in Radians, if Sherm's email is to be believed
	       next++;
        }
    }
}

void MoleculeClassContainer::includeAllAtoms( DuMMForceFieldSubsystem & dumm) {
	map<const String, MoleculeClass>::iterator it;
	map<const String, MoleculeClass>::iterator next;
	//int i = 0;
	next = moleculeClassMap.begin();
	while (next != moleculeClassMap.end())
	{
	   it = next;
           it->second.includeAllAtoms(dumm);
	   next++;
        }

}

bool MoleculeClassContainer::hasChainID( String chain) {
    MMBLOG_FILE_FUNC_LINE(INFO, "Inside MoleculeClassContainer::hasChainID("<<chain<<"). Found exactly "<<moleculeClassMap.count(chain)<<" molecules with chain ID >"<<chain<<"<."<<endl);
    if (moleculeClassMap.count(chain)== 1 ) {
        return true;}
    else if (moleculeClassMap.count(chain)== 0) {return false;} 
    else {
        MMBLOG_FILE_FUNC_LINE(CRITICAL, moleculeClassMap.count(chain) <<" MoleculeClass's found with chain ID "<<chain <<endl);
    }
}

void MoleculeClassContainer::addConstraintToGround(map<const String,double> myUserVariables,  const String chain, const String atomName, ConstraintToGroundContainer & constraintToGroundContainer){
    constraintToGroundContainer.addConstraintClassToVector(
        chain,
        ResidueID(1,(' ')), 
        atomName
        );
}



