/* -------------------------------------------------------------------------- *
 *                           MMB (MacroMoleculeBuilder)                       *
 * -------------------------------------------------------------------------- *
 *                                                                            *
 * Copyright (c) 2011-12 by the Author.                                       *
 * Author: Samuel Flores                                                      *
 *                                                                            *
 * See RNABuilder.cpp for the copyright and usage agreement.                  *
 * -------------------------------------------------------------------------- */

#ifndef AddBackboneOxygenForces_H_
#define AddBackboneOxygenForces_H_
using namespace std;
using namespace SimTK;
class AddBackboneOxygenForces : public Biopolymer {
public:
    DuMM::AtomClassIndex backboneOxygenLoopResidueIx;//(1);
    DuMM::ChargedAtomTypeIndex backboneOxygenChargedAtomTypeIx;//(1);
    DuMM::ChargedAtomTypeIndex magnesiumChargedAtomTypeIx;//(1);

int initialize (  DuMMForceFieldSubsystem & backboneOxygenForces , ParameterReader & myParameterReader ){
    //if backboneOxygenForces.loadAmber99Parameters();
        filebuf fb;
        fb.open ((myParameterReader.tinkerParameterFileName).c_str(),ios::in);
        istream is(&fb);
        if (myParameterReader.loadTinkerParameterFile)
            backboneOxygenForces.populateFromTinkerParameterFile (is);
        else 
            backboneOxygenForces.loadAmber99Parameters();
        fb.close();

    backboneOxygenForces.setBondTorsionGlobalScaleFactor (0);
    backboneOxygenForces.setBondStretchGlobalScaleFactor (0);
    backboneOxygenForces.setBondBendGlobalScaleFactor (0);
    backboneOxygenForces.setGbsaGlobalScaleFactor(0);
    backboneOxygenForces.setCoulombGlobalScaleFactor(myParameterReader.backboneOxygenGlobalCoulombScaleFactor);
    backboneOxygenForces.setVdwGlobalScaleFactor(myParameterReader.backboneOxygenGlobalVdwScaleFactor);
    backboneOxygenLoopResidueIx = DuMM::AtomClassIndex(25);
    backboneOxygenChargedAtomTypeIx = DuMM::ChargedAtomTypeIndex(1231);
    return(0); 
}


int addForces (CompoundSystem & system, SimbodyMatterSubsystem & matter , DuMMForceFieldSubsystem & backboneOxygenForces, ParameterReader & myParameterReader, MagnesiumIon myMagnesiumIon )
      {

        MobilizedBodyIndex bodyIx ;
        DuMM::ClusterIndex magnesiumClusterIx  = backboneOxygenForces.createCluster("Mg+2" );
        stringstream ss3;
        MobilizedBody myMobilizedBody;
        
        Compound::AtomIndex op1Ix(0);
        bodyIx =myMagnesiumIon.getAtomMobilizedBodyIndex(op1Ix);
        DuMM::ChargedAtomTypeIndex magnesiumIonChargedAtomTypeIx(2008);
        if (myParameterReader.verbose) cout<<"[AddBackboneOxygenForces.h] bodyIx ="<<bodyIx<<endl;
        Vec3 atomStation = myMagnesiumIon.getAtomLocationInMobilizedBodyFrame(op1Ix);
        if (myParameterReader.verbose) cout<<"[AddBackboneOxygenForces.h] atomStation ="<<atomStation<<endl;
        DuMM::AtomIndex magnesiumIonAtomIndex = backboneOxygenForces.addAtom(magnesiumIonChargedAtomTypeIx);
        backboneOxygenForces.placeAtomInCluster (magnesiumIonAtomIndex,magnesiumClusterIx,atomStation);
        backboneOxygenForces.attachClusterToBody(magnesiumClusterIx,bodyIx,Transform(Vec3(0)));
        return(0);
}

int addForces (CompoundSystem & system, SimbodyMatterSubsystem & matter , DuMMForceFieldSubsystem & backboneOxygenForces, ParameterReader & myParameterReader, Biopolymer & myChain ,int firstResidueIndex , int lastResidueIndex )
      {


    // This is missing some of the forces found in BackboneOxygen
 
    //myChain.setCompoundBondMobility(BondMobility::Free);
    DuMM::AtomIndex previousBackboneOxygenAtom;
    MobilizedBodyIndex bodyIx ;
    MobilizedBodyIndex nextBodyIx ;
    MobilizedBodyIndex oldBodyIx =myChain.getAtomMobilizedBodyIndex(myChain.getAtomIndex("1/OP1") );
    //MobilizedBodyIndex oldBodyIx = MobilizedBodyIndex(0);//myChain.getAtomMobilizedBodyIndex(myChain.getAtomIndex("0/OP1") );
    DuMM::ClusterIndex clusterIx  = backboneOxygenForces.createCluster("1/OP1");
    DuMM::ClusterIndex oldClusterIx  = backboneOxygenForces.createCluster("1/OP1");
    stringstream ss3;
    //for (Compound::Index resIx(1); resIx < (myChain).getNumResidues(); ++resIx) { //start at second residue, since the first has no phosphate
    for (ResidueInfo::Index resIx(firstResidueIndex); resIx <= ResidueInfo::Index( lastResidueIndex); ++resIx) { //start at second residue, since the first has no phosphate
        
        if (resIx < ((myChain).getNumResidues()-1)){
            const ResidueInfo & nextResidue = (myChain).getResidue(ResidueInfo::Index(int(resIx)+1));
            //const BiopolymerResidue& nextResidue = (myChain).updResidue(ResidueInfo::Index(int(resIx)+1));
            Compound::AtomIndex nextOP1Ix = nextResidue.getAtomIndex("OP1");
            nextBodyIx = myChain.getAtomMobilizedBodyIndex(nextOP1Ix);           
        }
        const ResidueInfo & residue = (myChain).getResidue(resIx);
        if (myParameterReader.verbose) cout<<"[AddBackboneOxygenForces.h] resIx ="<<resIx<<endl;
        Compound::AtomIndex op1Ix = residue.getAtomIndex("OP1");
        if (myParameterReader.verbose) cout<<"[AddBackboneOxygenForces.h] op1Ix ="<<op1Ix<<endl;
        bodyIx = myChain.getAtomMobilizedBodyIndex(op1Ix);
        if (myParameterReader.verbose) cout<<"[AddBackboneOxygenForces.h] bodyIx ="<<bodyIx<<endl;
        Vec3 atomStation = myChain.getAtomLocationInMobilizedBodyFrame(op1Ix);
        if (myParameterReader.verbose) cout<<"[AddBackboneOxygenForces.h] atomStation ="<<atomStation<<endl;
        if (myParameterReader.verbose) cout<<"[AddBackboneOxygenForces.h] backboneOxygenChargedAtomTypeIx ="<<backboneOxygenChargedAtomTypeIx<<endl;
        DuMM::AtomIndex backboneOxygenAtom = backboneOxygenForces.addAtom(backboneOxygenChargedAtomTypeIx);
       	if (int(resIx) == 1) {
            clusterIx = backboneOxygenForces.createCluster("1/OP1");
            backboneOxygenForces.placeAtomInCluster (backboneOxygenAtom,clusterIx,atomStation);
            if ((resIx == ((myChain).getNumResidues()-1)) || (nextBodyIx != bodyIx )){
                backboneOxygenForces.attachClusterToBody(oldClusterIx,oldBodyIx,Transform(Vec3(0)));
                if (myParameterReader.verbose) cout<<"[AddBackboneOxygenForces.h] attaching cluster to rigid body"<<endl;
            }
        } else 
        if (bodyIx !=  oldBodyIx) {
            ss3.clear();
	    ss3<<resIx<<"/OP1";
       	    clusterIx = backboneOxygenForces.createCluster((ss3.str()).c_str());
            backboneOxygenForces.placeAtomInCluster (backboneOxygenAtom,clusterIx,atomStation);
            if ((resIx == ((myChain).getNumResidues()-1)) || (nextBodyIx != bodyIx )) {
                backboneOxygenForces.attachClusterToBody(clusterIx,bodyIx,Transform(Vec3(0)));
                if (myParameterReader.verbose) cout<<"[AddBackboneOxygenForces.h] attaching cluster to non-rigid body"<<endl;
            }
        } else {
            if (myParameterReader.verbose) cout<<"[AddBackboneOxygenForces.h] rigid cluster detected"<<endl;
            backboneOxygenForces.placeAtomInCluster (backboneOxygenAtom,oldClusterIx,atomStation);
            if ((resIx == ((myChain).getNumResidues()-1)) || (nextBodyIx != bodyIx )){
                backboneOxygenForces.attachClusterToBody(oldClusterIx,oldBodyIx,Transform(Vec3(0)));
                if (myParameterReader.verbose) cout<<"[AddBackboneOxygenForces.h] attaching cluster to rigid body"<<endl;
            }
	}
                
        //backboneOxygenForces.attachAtomToBody(backboneOxygenAtom, bodyIx, atomStation);
        if (previousBackboneOxygenAtom.isValid() && backboneOxygenAtom.isValid()) {
            //backboneOxygenForces.addBond(previousBackboneOxygenAtom, backboneOxygenAtom);
        oldBodyIx = bodyIx;
        oldClusterIx = clusterIx; 
        }
        previousBackboneOxygenAtom = backboneOxygenAtom;
   }
   return(0);

}

};
#endif
