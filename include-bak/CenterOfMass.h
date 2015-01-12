/* -------------------------------------------------------------------------- *
 *                           MMB (MacroMoleculeBuilder)                       *
 * -------------------------------------------------------------------------- *
 *                                                                            *
 * Copyright (c) 2011-12 by the Author.                                       *
 * Author: Samuel Flores                                                      *
 *                                                                            *
 * See RNABuilder.cpp for the copyright and usage agreement.                  *
 * -------------------------------------------------------------------------- */

#include "SimTKmolmodel.h"
#include "SimTKsimbody_aux.h"
#include <vector>
#include <list>
using namespace SimTK;
using namespace std;


/*
Vec3 GeometricCenter (Compound & myCompound, SimbodyMatterSubsystem & matter, State & state) {
    Vec3 center(0);
    for (int i = 0 ; i < myCompound.getNAtoms(); i++) {
        center += (myCompound.getAtomMobilizedBody(AtomIndex(i))).findBodyTransformInAnotherBody(state,matter.Ground()) / myCompound.getNAtoms();
    }
    return center;
}

Vec3 GeometricCenter (Biopolymer & myCompound, int firstResidue, int lastResidue, Matter & matter, State & state) { 
    Vec3 center(0);
    for (int i = firstResidue ; i <= lastResidue) {
        Compound myResidue = myCompound.updResidue(AtomIndex(i));
	for (int j = 0; j <= myResidue.getNAtoms(); j++) {
	    center += matter.getMobilizedBody(myResidue.getAtomMobilizedBodyIndex(AtomIndex(j))).findBodyTransformInAnotherBody(state,matter.Ground()) / myResidue.getNAtoms();
        }	
    }   
    return center;
}
*/

Vec3 CenterOfMass (Biopolymer & myCompound , list<int> residueList, SimbodyMatterSubsystem & matter, const State& state, TinkerDuMMForceFieldSubsystem & dumm) {
    Vec3 center(0);
    list<int>::iterator i;
    float totalMass = 0;
    for(i=residueList.begin(); i != residueList.end(); ++i) {
	//cout<<"[GeometricCenter.h : GeometricCenter] i ="<<*i<<endl;
	Compound myResidue = myCompound.updResidue(Compound::Index( *i ));
        //cout<<"[GeometricCenter.h] : myResidue.getNAtoms()"<<myResidue.getNAtoms()<<endl;
        for (int j = 0; j < myResidue.getNAtoms(); j++) 
        //int j = 0;
        { 
            //cout<<"[GeometricCenter.h] :         j ="           << j  <<endl;
            //cout<<"[GeometricCenter.h] : adding atom location ="<<myResidue.calcAtomLocationInGroundFrame(state,Compound::AtomIndex(j))<<endl;
            //cout<<"[GeometricCenter.h] : atom mass ="           << dumm.getAtomMass(DuMM::AtomIndex(j)) <<endl;
            center += myResidue.calcAtomLocationInGroundFrame(state,Compound::AtomIndex(j))*dumm.getAtomMass(DuMM::AtomIndex(j));   ///1/myResidue.getNAtoms(); //
            totalMass += dumm.getAtomMass(DuMM::AtomIndex(j));
            
	    // center += (matter.getMobilizedBody(myResidue.getAtomMobilizedBodyIndex(Compound::AtomIndex(j))).findBodyTransformInAnotherBody(state,matter.Ground())).T() / myResidue.getNAtoms();
            
        }
	

    }
    center = center/totalMass;
    return center;
}


