/* -------------------------------------------------------------------------- *
 *                           MMB (MacroMoleculeBuilder)                       *
 * -------------------------------------------------------------------------- *
 *                                                                            *
 * Copyright (c) 2011-12 by the Author.                                       *
 * Author: Samuel Flores                                                      *
 *                                                                            *
 * See RNABuilder.cpp for the copyright and usage agreement.                  *
 * -------------------------------------------------------------------------- */

#ifndef AddHuntCrossleySpheres_H_
#define AddHuntCrossleySpheres_H_
#include <stdio.h>
#include <stdlib.h>
#include "Utils.h"

/*
static int AddExcludedVolume(Biopolymer & myMolecule,GeneralContactSubsystem & contacts,ContactSetIndex contactSet,HuntCrossleyForce & hc,SimbodyMatterSubsystem & matter,ParameterReader myParameterReader, string leontisWesthofInteractionType, int startResidue, int endResidue) {

    exit(1); // this is obsolete , never rigorously tested, shouldnt be used.
}
*/
/*static int AddRNABackboneSterics(Biopolymer & myMolecule,GeneralContactSubsystem & contacts,ContactSetIndex contactSet,HuntCrossleyForce & hc,SimbodyMatterSubsystem & matter,ParameterReader myParameterReader,  int startResidue, int endResidue,bool endCapsOn) {

// start adding spheres
     float huntCrossleyDissipation  =.0; 
     for (int q=startResidue; q<=endResidue;q++) {
            ResidueInfo myResidueInfo = myMolecule.updResidue(ResidueInfo::Index(q+int(endCapsOn)));
            String myPdbResidueName1 =  myResidueInfo.getPdbResidueName();            
            for (ResidueInfo::AtomIndex r(0) ; r<ResidueInfo::AtomIndex(  (myResidueInfo).getNumAtoms()); r++)
            { //loop over four possible interacting pairs of atoms.
                { //if the atom name field is not blank, do this 
                        stringstream ss3; 
                        // excluded volume radius WAS 1.25.  but now using much smaller radii. 
                        // following:
 			//Whitford PC, Noel JK, Gosavi S, Schug A, Sanbonmatsu KY & Onuchic JN, "An All-atom Structure-Based Potential for Proteins: Bridging Minimal Models with All-atom Empirical Forcefields" PROTEINS (2008) DOI: 10.1002/prot.22253. 
                        ss3<<q<<"/"<<(myResidueInfo).getAtomName(ResidueInfo::AtomIndex(r));
                       if (//myMolecule.hasAtom(ss3.str()) &&
                           ((myResidueInfo.getAtomName(r)).compare("P"  ) == 0) ||
                           ((myResidueInfo.getAtomName(r)).compare("O5*") == 0) ||
                           ((myResidueInfo.getAtomName(r)).compare("C5*") == 0) ||
                           ((myResidueInfo.getAtomName(r)).compare("C4*") == 0) ||
                           ((myResidueInfo.getAtomName(r)).compare("C3*") == 0) ||
                           ((myResidueInfo.getAtomName(r)).compare("O3*") == 0) ||
                           ((myResidueInfo.getAtomName(r)).compare("O5'") == 0) ||
                           ((myResidueInfo.getAtomName(r)).compare("C5'") == 0) ||
                           ((myResidueInfo.getAtomName(r)).compare("C4'") == 0) ||
                           ((myResidueInfo.getAtomName(r)).compare("C3'") == 0) ||
                           ((myResidueInfo.getAtomName(r)).compare("O3'") == 0) 
                          )  
                       {
                           
                           SimTK_ERRCHK1_ALWAYS(
                               (myMolecule.hasAtom(ss3.str()) ),
                               __FILE__,
                               "Failed to attach steric sphere.  Could not find specified atom: %s .",ss3.str().c_str());
                           contacts.addBody(contactSet,
                                        (matter.updMobilizedBody(myMolecule.getAtomMobilizedBodyIndex(Compound::AtomIndex(myMolecule.getAtomIndex(ss3.str()))))), 
                                        ContactGeometry::Sphere(myParameterReader.excludedVolumeRadius/10), //convert from angstrom to nanometers.
                                        (myMolecule.getAtomLocationInMobilizedBodyFrame(myMolecule.getAtomIndex(ss3.str())))
                                        );                           
                           hc.setBodyParameters(SimTK::ContactSurfaceIndex(contacts.getNumBodies(contactSet)-1),myParameterReader.excludedVolumeStiffness ,huntCrossleyDissipation, 0., 0., 0.); 
		       }
                        //hc.setBodyParameters(contacts.getNumBodies(contactSet)-1,huntCrossleyStiffness ,huntCrossleyDissipation, 0., 0., 0.); 
     
                }    
            }    
     


     
     }    
     return(0);
}*/

/*
static int AddAllHeavyAtomSterics(Biopolymer & myMolecule,GeneralContactSubsystem & contacts,ContactSetIndex contactSet,HuntCrossleyForce & hc,SimbodyMatterSubsystem & matter,ParameterReader myParameterReader,  int startResidue, int endResidue, bool endCapsOn) {

// start adding spheres
     float huntCrossleyDissipation  =.0; 
     
     for (int q=startResidue; q<=endResidue;q++) {
            ResidueInfo myResidueInfo = myMolecule.updResidue(ResidueInfo::Index(q+int(endCapsOn)));
            String myPdbResidueName1 =  myResidueInfo.getPdbResidueName();            
            for (ResidueInfo::AtomIndex r(0) ; r<ResidueInfo::AtomIndex(  myResidueInfo.getNumAtoms()); r++)
            { //loop over four possible interacting pairs of atoms.
                { //if the atom name field is not blank, do this 
                        stringstream ss3; 
                        ss3<<q<<"/"<<myResidueInfo.getAtomName(r);
                        // excluded volume radius WAS 1.25.  but now using much smaller radii. 
                        // following:
 			//Whitford PC, Noel JK, Gosavi S, Schug A, Sanbonmatsu KY & Onuchic JN, "An All-atom Structure-Based Potential for Proteins: Bridging Minimal Models with All-atom Empirical Forcefields" PROTEINS (2008) DOI: 10.1002/prot.22253. 


                        SimTK_ERRCHK_ALWAYS(
                            (myMolecule.hasAtom(ss3.str()) ),
                            __FILE__,
                            "Failed to attach steric sphere.  Could not find specified atom");//: %s .",string(ss3.str()));

                       if ((myMolecule.hasAtom(ss3.str())) &&
                          (((myMolecule.getAtomElement(myResidueInfo.getAtomIndex( r  ))).getSymbol()).compare("H") != 0))
                       {
                           contacts.addBody(contactSet,
                                        (matter.updMobilizedBody(myMolecule.getAtomMobilizedBodyIndex(Compound::AtomIndex(myMolecule.getAtomIndex(ss3.str()))))), 
                                        ContactGeometry::Sphere(myParameterReader.excludedVolumeRadius/10), //convert from angstrom to nanometers.
                                        (myMolecule.getAtomLocationInMobilizedBodyFrame(myMolecule.getAtomIndex(ss3.str())))
                                        );                           
                           hc.setBodyParameters(SimTK::ContactSurfaceIndex(contacts.getNumBodies(contactSet)-1),myParameterReader.excludedVolumeStiffness ,huntCrossleyDissipation, 0., 0., 0.); 
		       }
                }    
            }    
     }    
     return(0);
}
*/

/*
static int AddAllAtomSterics(Biopolymer & myMolecule,GeneralContactSubsystem & contacts,ContactSetIndex contactSet,HuntCrossleyForce & hc,SimbodyMatterSubsystem & matter,ParameterReader myParameterReader,  int startResidue, int endResidue, bool endCapsOn) {

// start adding spheres
     float huntCrossleyDissipation  =.0; 
     for (int q=startResidue; q<=endResidue;q++) {
            ResidueInfo myResidueInfo = myMolecule.updResidue(ResidueInfo::Index(q+int(endCapsOn)));
            String myPdbResidueName1 =  myResidueInfo.getPdbResidueName();            
            for (ResidueInfo::AtomIndex r(0) ; r<ResidueInfo::AtomIndex(  myResidueInfo.getNumAtoms()); r++)
            { //loop over four possible interacting pairs of atoms.
                { //if the atom name field is not blank, do this 
                        stringstream ss3; 
                        ss3<<q<<"/"<<myResidueInfo.getAtomName(r);
                        //cout<<ss3.str()<<endl;
                        // excluded volume radius WAS 1.25.  but now using much smaller radii. 
                        // following:
 			//Whitford PC, Noel JK, Gosavi S, Schug A, Sanbonmatsu KY & Onuchic JN, "An All-atom Structure-Based Potential for Proteins: Bridging Minimal Models with All-atom Empirical Forcefields" PROTEINS (2008) DOI: 10.1002/prot.22253. 


                        SimTK_ERRCHK_ALWAYS(
                            (myMolecule.hasAtom(ss3.str()) ),
                            __FILE__,
                            "Failed to attach steric sphere.  Could not find specified atom");//: %s .",string(ss3.str()));

                       if (myMolecule.hasAtom(ss3.str())) //&&
                       {
                           contacts.addBody(contactSet,
                                        (matter.updMobilizedBody(myMolecule.getAtomMobilizedBodyIndex(Compound::AtomIndex(myMolecule.getAtomIndex(ss3.str()))))), 
                                        ContactGeometry::Sphere(myParameterReader.excludedVolumeRadius/10), //convert from angstrom to nanometers.
                                        (myMolecule.getAtomLocationInMobilizedBodyFrame(myMolecule.getAtomIndex(ss3.str())))
                                        );                           
                           hc.setBodyParameters(SimTK::ContactSurfaceIndex(contacts.getNumBodies(contactSet)-1),myParameterReader.excludedVolumeStiffness ,huntCrossleyDissipation, 0., 0., 0.); 
		       }
                }    
            }    
     }    
     return(0);
}
*/
/*
static int AddHuntCrossleySpheres(Biopolymer & myMolecule,GeneralContactSubsystem & contacts,ContactSetIndex contactSet,HuntCrossleyForce & hc,SimbodyMatterSubsystem & matter,LeontisWesthofClass myLeontisWesthofClass,string leontisWesthofInteractionType, int startResidue, int endResidue, bool endCapsOn) {
//exit(1);
// start adding spheres
     //float huntCrossleyStiffness    =1000000;
     float huntCrossleyDissipation  =.0; 
     //float huntCrossleyStiffness    =1000;
     //float huntCrossleyDissipation  =10;      //stringstream ss3; 
     for (int q=startResidue; q<=endResidue;q++) {
            //cout<<"check 2"<<endl;
            ResidueInfo myResidueInfo = myMolecule.updResidue(ResidueInfo::Index(q+int(endCapsOn)));
            String myPdbResidueName1 = (myResidueInfo).getPdbResidueName();            
            //cout<<"check 2.3"<<endl;
            LeontisWesthofBondRow myLeontisWesthofBondRow = myLeontisWesthofClass.getLeontisWesthofBondRow(
                (myResidueInfo).getPdbResidueNumber(),
                (myResidueInfo).getPdbResidueNumber() ,
                myPdbResidueName1,
                leontisWesthofInteractionType,"",leontisWesthofInteractionType,"X","contact");
            //cout<<"check 2.35"<<endl;
            for (int r =0; r<4; r++) { //loop over four possible interacting pairs of atoms.
                if ((myLeontisWesthofBondRow.residue1Atom[r]).compare("") != 0) { //if the atom name field is not blank, do this 
                        
                        stringstream ss3; 
                        ss3<<q<<"/"<<myLeontisWesthofBondRow.residue1Atom[r]; 
                       if (myMolecule.hasAtom(ss3.str())) {
                           contacts.addBody(contactSet,
                                        (matter.updMobilizedBody(myMolecule.getAtomMobilizedBodyIndex(Compound::AtomIndex(myMolecule.getAtomIndex(ss3.str()))))), 
                                        ContactGeometry::Sphere(myLeontisWesthofBondRow.bondLength[r]/10), 
                                        (myMolecule.getAtomLocationInMobilizedBodyFrame(myMolecule.getAtomIndex(ss3.str())))
                                        );                           
                           hc.setBodyParameters(SimTK::ContactSurfaceIndex(contacts.getNumBodies(contactSet)-1),myLeontisWesthofBondRow.springConstant[r] ,huntCrossleyDissipation, 0., 0., 0.); 
		       } else if (
                           (ss3.str().compare("0/P")==0) ||
                           (ss3.str().compare("0/OP1")==0) ||
                           (ss3.str().compare("0/OP2")==0) 
                           )
                       {
                           cout<<"[AddHuntCrossleySpheres.h] Warning:  You are attempting to place a contact sphere on an atom that doesn't exist, in this case an omitted 5' Phosphate; this is pretty harmless.                 "<<endl;
                           
                       } else {
                           // LOG AND CRASH HERE
                       }
                        //hc.setBodyParameters(contacts.getNumBodies(contactSet)-1,huntCrossleyStiffness ,huntCrossleyDissipation, 0., 0., 0.); 
     
                }    
            }    
     


     
     }    
     return(0);
}
*/
/*
int AddHuntCrossleySpheres(Biopolymer & myMolecule,GeneralContactSubsystem & contacts,ContactSetIndex contactSet,HuntCrossleyForce & hc,SimbodyMatterSubsystem & matter,LeontisWesthofClass myLeontisWesthofClass,string leontisWesthofInteractionType) {

// start adding spheres
     //float huntCrossleyStiffness    =1000000;
     float huntCrossleyDissipation  =.0; 
     //float huntCrossleyStiffness    =1000;
     //float huntCrossleyDissipation  =10;      //stringstream ss3; 
     for (int q=0;q<myMolecule.getNResidues();q++) {
            String myPdbResidueName1 = (myMolecule.updResidue(ResidueInfo::Index(q))).getPdbResidueName();            LeontisWesthofBondRow myLeontisWesthofBondRow = myLeontisWesthofClass.getLeontisWesthofBondRow((myMolecule.updResidue(ResidueInfo::Index(q))).getPdbResidueNumber(),(myMolecule.updResidue(ResidueInfo::Index(q))).getPdbResidueNumber() , myPdbResidueName1,leontisWesthofInteractionType,"","","",0);
            for (int r =0; r<4; r++) { //loop over four possible interacting pairs of atoms.
                if ((myLeontisWesthofBondRow.residue1Atom[r]).compare("") != 0) { //if the atom name field is not blank, do this 
                        stringstream ss3; 
                        ss3<<q<<"/"<<myLeontisWesthofBondRow.residue1Atom[r];
                        //cout<<"about to add hard sphere to atom :"<<ss3.str()<<","<<endl;
                       if (myMolecule.hasAtom((ss3.str()))) {
                       contacts.addBody(contactSet,
                                        (matter.updMobilizedBody(myMolecule.getAtomMobilizedBodyIndex(Compound::AtomIndex(myMolecule.getAtomIndex(ss3.str()))))), 
                                        ContactGeometry::Sphere(myLeontisWesthofBondRow.bondLength[r]), 
                                        (myMolecule.getAtomLocationInMobilizedBodyFrame(myMolecule.getAtomIndex(ss3.str())))
                                        );                           
                       hc.setBodyParameters(contacts.getNumBodies(contactSet)-1,myLeontisWesthofBondRow.springConstant[r] ,huntCrossleyDissipation, 0., 0., 0.); 
                       }
                        //hc.setBodyParameters(contacts.getNumBodies(contactSet)-1,huntCrossleyStiffness ,huntCrossleyDissipation, 0., 0., 0.); 
     
                }    
            }    
     


     
     }    
     return(0);
}
*/
#endif
