/* -------------------------------------------------------------------------- *
 *                           MMB (MacroMoleculeBuilder)                       *
 * -------------------------------------------------------------------------- *
 *                                                                            *
 * Copyright (c) 2011-12 by the Author.                                       *
 * Author: Samuel Flores                                                      *
 *                                                                            *
 * See RNABuilder.cpp for the copyright and usage agreement.                  *
 * -------------------------------------------------------------------------- */

#ifndef AddNASTForces_H_
#define AddNASTForces_H_

class AddNASTForces : public Biopolymer {
public:
    DuMM::AtomClassIndex nastLoopResidueIx;
    DuMM::ChargedAtomTypeIndex nastChargedAtomTypeIx;
int initialize (  DuMMForceFieldSubsystem & nastForces , ParameterReader & myParameterReader ){
    nastForces.setBondTorsionGlobalScaleFactor (myParameterReader.nastGlobalBondTorsionScaleFactor); // found that a factor of 1 was not large enough to make a floppy chain into a helix at the default temperature, at least not very quickly
    nastForces.setBondStretchGlobalScaleFactor (1);//setGlobalTorsionScaleFactor (0);
    nastForces.setBondBendGlobalScaleFactor (myParameterReader.nastGlobalBondTorsionScaleFactor);//setGlobalTorsionScaleFactor (0);
    nastForces.setGbsaGlobalScaleFactor(0);
    nastLoopResidueIx = DuMM::AtomClassIndex(1);
    nastForces.defineAtomClass(nastLoopResidueIx, "", 6, 2,
        .4,//0.5,  // Rmin, nm
        4.184); // well depth
    nastChargedAtomTypeIx = DuMM::ChargedAtomTypeIndex(1);
    nastForces.defineChargedAtomType(nastChargedAtomTypeIx, "",
    nastLoopResidueIx, 0.0);
    nastForces.defineBondStretch(nastLoopResidueIx, nastLoopResidueIx,
        779,//49.0, // bond stretch stiffness
        .578);//1.37); // bond dead length, nm
    nastForces.defineBondBend(nastLoopResidueIx,
       nastLoopResidueIx, nastLoopResidueIx,
       10,//20.0, // stiffness
       2.4*Rad2Deg);//2.6 * Rad2Deg); // radians
    // Phase must be in range 0 to 180
    nastForces.defineBondTorsion(nastLoopResidueIx,
        nastLoopResidueIx, nastLoopResidueIx, nastLoopResidueIx,
        //1, 50.0,  - 2.8*Rad2Deg); // periodicity, amplitude,
        1, 5.0,  - 2.8*Rad2Deg); // periodicity, amplitude,
        //phase (in degrees!?!)
// 19.57 = -160  original
// 160.43 = -19.74
// should be + 22
   return(0);
}
int addForces (CompoundSystem & system, SimbodyMatterSubsystem & matter , DuMMForceFieldSubsystem & nastForces, ParameterReader & myParameterReader, Biopolymer & myChain )
      {


    // I got the nonbonded, stretch, bend, and torsion parameters for
    //helical residues from Magda J. --cmb
    // This is missing some of the forces found in NAST
    // Sam has to figure out how to add the rest
    // and also add a second "atom" type for non-helical residues
 
    //myChain.setCompoundBondMobility(BondMobility::Free);
    DuMM::AtomIndex previousNastAtom;
    MobilizedBodyIndex bodyIx ;
    MobilizedBodyIndex nextBodyIx ;
    MobilizedBodyIndex oldBodyIx     = myChain.getAtomMobilizedBodyIndex(Compound::AtomIndex(0));
    DuMM::ClusterIndex clusterIx     = nastForces.createCluster("0/C3*");
    DuMM::ClusterIndex oldClusterIx  = nastForces.createCluster("0/C3*");
    stringstream ss3;
    for (ResidueInfo::Index resIx(0); resIx < (myChain).getNumResidues(); ++resIx) {
        
        if (resIx < ((myChain).getNumResidues()-1)){
            const ResidueInfo& nextResidue = (myChain).getResidue(ResidueInfo::Index(int(resIx)+1));
            Compound::AtomIndex nextC3Ix =   nextResidue.getAtomIndex("C3*");
            //Compound::AtomIndex nextC3Ix =    nextResidue.getAtomIndex("C3*");
            nextBodyIx = myChain.getAtomMobilizedBodyIndex( nextC3Ix   );// nextResidue.getAtomMobilizedBodyIndex(nextC3Ix);           
        }
        const ResidueInfo & residue = (myChain).updResidue(resIx);
        if (myParameterReader.verbose) cout<<"[AddNASTForces.h] resIx ="<<resIx<<endl;
        Compound::AtomIndex c3Ix = residue.getAtomIndex("C3*");
        if (myParameterReader.verbose) cout<<"[AddNASTForces.h] c3Ix ="<<c3Ix<<endl;
        bodyIx = myChain.getAtomMobilizedBodyIndex(c3Ix);
        //if (bodyIx == oldBodyIx)  
	  //  bodyIx ;
        if (myParameterReader.verbose) cout<<"[AddNASTForces.h] bodyIx ="<<bodyIx<<endl;
        Vec3 atomStation = myChain.getAtomLocationInMobilizedBodyFrame(c3Ix);
        if (myParameterReader.verbose) cout<<"[AddNASTForces.h] atomStation ="<<atomStation<<endl;
        if (myParameterReader.verbose) cout<<"[AddNASTForces.h] nastChargedAtomTypeIx ="<<nastChargedAtomTypeIx<<endl;
        DuMM::AtomIndex nastAtom = nastForces.addAtom(nastChargedAtomTypeIx);
        if (bodyIx !=  oldBodyIx) {
            ss3.clear();
	    ss3<<resIx<<"/C3*";
       	    clusterIx = nastForces.createCluster((ss3.str()).c_str());
            nastForces.placeAtomInCluster (nastAtom,clusterIx,atomStation);
            if ((resIx == ((myChain).getNumResidues()-1)) || (nextBodyIx != bodyIx )) {
                nastForces.attachClusterToBody(clusterIx,bodyIx,Transform(Vec3(0)));
                if (myParameterReader.verbose) cout<<"[AddNASTForces.h] attaching cluster to non-rigid body"<<endl;
            }
        } else {
            if (myParameterReader.verbose) cout<<"[AddNASTForces.h] rigid cluster detected"<<endl;
            nastForces.placeAtomInCluster (nastAtom,oldClusterIx,atomStation);
            if ((resIx == ((myChain).getNumResidues()-1)) || (nextBodyIx != bodyIx )){
                nastForces.attachClusterToBody(oldClusterIx,oldBodyIx,Transform(Vec3(0)));
                if (myParameterReader.verbose) cout<<"[AddNASTForces.h] attaching cluster to rigid body"<<endl;
            }
	}
                
        //nastForces.attachAtomToBody(nastAtom, bodyIx, atomStation);
        if (previousNastAtom.isValid() && nastAtom.isValid()) {
            nastForces.addBond(previousNastAtom, nastAtom);
        oldBodyIx = bodyIx;
        oldClusterIx = clusterIx; 
        }
        previousNastAtom = nastAtom;
   }


    return 0;
}

};
#endif
