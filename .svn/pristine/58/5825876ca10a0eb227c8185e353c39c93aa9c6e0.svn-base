#ifndef WadleyKeatingDuartePyleTorsionForce_H_
#define WadleyKeatingDuartePyleTorsionForce_H_
/*

class SinusoidalTorsion : public DuMM::CustomBondTorsion {
public:
    SinusoidalTorsion(int periodicity, Real halfAmplitudeInKJPerMole, Real phaseInRadians) 
        : periodicity(periodicity), halfAmplitude(halfAmplitudeInKJPerMole), phase(phaseInRadians)
    {}  

    Real calcEnergy(Real torsionInRadians) const {
        return halfAmplitude * ( 1.0 + std::cos(periodicity * torsionInRadians - phase) );
    }   

    Real calcTorque(Real torsionInRadians) const {
        return periodicity * halfAmplitude * ( std::sin(periodicity * torsionInRadians - phase) );
    }   

private:
    int periodicity;
    Real halfAmplitude;
    Real phase;
};

*/



class WKDPTorsion : public DuMM::CustomBondTorsion {
public:
   WKDPTorsion(
       double wkdpGlobalBondTorsionScaleFactor
       //float amplitude1, Angle phase1, Angle stdDeviation1,
       //float amplitude2=0, Angle phase2=0, Angle stdDeviation2=0,
       //float amplitude3=0, Angle phase3=0, Angle stdDeviation3=0,
       //float amplitude4=0, Angle phase4=0, Angle stdDeviation4=0
   )
     : 
     wkdpGlobalBondTorsionScaleFactor(wkdpGlobalBondTorsionScaleFactor)
     //amplitude1(amplitude1), phase1(phase1), stdDeviation1(stdDeviation1), 
     //amplitude2(amplitude2), phase2(phase2), stdDeviation2(stdDeviation2), 
     //amplitude3(amplitude3), phase3(phase3), stdDeviation3(stdDeviation3) ,
     //amplitude4(amplitude4), phase4(phase4), stdDeviation4(stdDeviation4) 
     { }

   Real calcTorque(Angle dihedralAngle) const {
       Real   myTorque = 0;
       myTorque += 0;
       return myTorque*wkdpGlobalBondTorsionScaleFactor;
   }

   Real   calcEnergy(Angle dihedralAngle) const
   {
       float myEnergy = 0;
       myEnergy += 0;
       return  myEnergy * wkdpGlobalBondTorsionScaleFactor;
   }

private:
  double wkdpGlobalBondTorsionScaleFactor;
  /*float amplitude1;
  float amplitude2;
  float amplitude3;
  float amplitude4;
  Angle phase1;
  Angle phase2;
  Angle phase3;
  Angle phase4;
  int stdDeviation1;
  int stdDeviation2;
  int stdDeviation3;
  int stdDeviation4;*/
};


class DefineWKDPTorsion {
public:
    DefineWKDPTorsion (DuMMForceFieldSubsystem & dumm, ParameterReader & myParameterReader )

    {
       // C4' is class 1, P is class 28
       dumm.defineCustomBondTorsion(
           SimTK::DuMM::AtomClassIndex(1 )  ,  SimTK::DuMM::AtomClassIndex(28) ,  SimTK::DuMM::AtomClassIndex(1) ,  SimTK::DuMM::AtomClassIndex( 28),
           new WKDPTorsion(myParameterReader.wkdpGlobalBondTorsionScaleFactor ) 
           );
       dumm.defineCustomBondTorsion(
           SimTK::DuMM::AtomClassIndex(28)  ,  SimTK::DuMM::AtomClassIndex(1) ,  SimTK::DuMM::AtomClassIndex(28) ,  SimTK::DuMM::AtomClassIndex( 1),
           new WKDPTorsion(myParameterReader.wkdpGlobalBondTorsionScaleFactor)
           );


    }

};

class AddWKDPForces : public Biopolymer {
public:
    DuMM::AtomClassIndex wkdpLoopResidueIx;
    DuMM::ChargedAtomTypeIndex wkdpChargedAtomTypeIx;
int initialize (  DuMMForceFieldSubsystem & wkdpForces , ParameterReader & myParameterReader ){
    wkdpForces.setBondTorsionGlobalScaleFactor (myParameterReader.wkdpGlobalBondTorsionScaleFactor); // found that a factor of 1 was not large enough to make a floppy chain into a helix at the default temperature, at least not very quickly
    wkdpForces.setBondStretchGlobalScaleFactor (1);//setGlobalTorsionScaleFactor (0);
    wkdpForces.setBondBendGlobalScaleFactor (myParameterReader.wkdpGlobalBondTorsionScaleFactor);//setGlobalTorsionScaleFactor (0);
    wkdpForces.setGbsaGlobalScaleFactor(0);
    wkdpLoopResidueIx = DuMM::AtomClassIndex(1);
    wkdpForces.defineAtomClass(wkdpLoopResidueIx, "", 6, 2,
        .4,//0.5,  // Rmin, nm
        4.184); // well depth
    wkdpChargedAtomTypeIx = DuMM::ChargedAtomTypeIndex(1);
    wkdpForces.defineChargedAtomType(wkdpChargedAtomTypeIx, "",
    wkdpLoopResidueIx, 0.0);
    wkdpForces.defineBondStretch(wkdpLoopResidueIx, wkdpLoopResidueIx,
        779,//49.0, // bond stretch stiffness
        .578);//1.37); // bond dead length, nm
    wkdpForces.defineBondBend(wkdpLoopResidueIx,
       wkdpLoopResidueIx, wkdpLoopResidueIx,
       10,//20.0, // stiffness
       2.4*Rad2Deg);//2.6 * Rad2Deg); // radians
    // Phase must be in range 0 to 180
    wkdpForces.defineBondTorsion(wkdpLoopResidueIx,
        wkdpLoopResidueIx, wkdpLoopResidueIx, wkdpLoopResidueIx,
        //1, 50.0,  - 2.8*Rad2Deg); // periodicity, amplitude,
        1, 5.0,  - 2.8*Rad2Deg); // periodicity, amplitude,
        //phase (in degrees!?!)
// 19.57 = -160  original
// 160.43 = -19.74
// should be + 22
   return(0);
}
int addForces (CompoundSystem & system, SimbodyMatterSubsystem & matter , DuMMForceFieldSubsystem & wkdpForces, ParameterReader & myParameterReader, Biopolymer & myChain )
      {


    // I got the nonbonded, stretch, bend, and torsion parameters for
    //helical residues from Magda J. --cmb
    // This is missing some of the forces found in WKDP
    // Sam has to figure out how to add the rest
    // and also add a second "atom" type for non-helical residues
 
    //myChain.setCompoundBondMobility(BondMobility::Free);

/*
    DuMM::AtomIndex previousWKDPAtom;
    MobilizedBodyIndex bodyIx ;
    MobilizedBodyIndex nextBodyIx ;
    MobilizedBodyIndex oldBodyIx = myChain.getAtomMobilizedBodyIndex(Compound::AtomIndex(0));
    DuMM::ClusterIndex clusterIx  = wkdpForces.createCluster("0/C4*");
    DuMM::ClusterIndex oldClusterIx  = clusterIx;
    //stringstream ss3;
    int currentLocationResidue (int locationNumber) {
        
	return int ((locationNumber + 1) div 2);
    }
    int nextLocationResidue (int locationNumber) {
	return int ((locationNumber + 2)   / 2);
    }
    int nextAtomIx (int locationNumber) {
        if (locationNumber % 2)
           ((myChain).updResidue(Compound::Index(nextLocationResidue(locationNumber)))).getAtomIndex("P"); 
        else
           ((myChain).updResidue(Compound::Index(nextLocationResidue(locationNumber)))).getAtomIndex("C4*"); 
	
    }
    string currentLocationAtomName (int locationNumber) {
        if (!(locationNumber % 2))
           return "P"; 
        else
           return "C4*";
	
    }
    int currentLocationAtomIx (int locationNumber) {
        ((myChain).updResidue(Compound::Index(currentLocationResidue(locationNumber)))).getAtomIndex(currentLocationAtomName(locationNumber)); 
    }
    Vec3 calcCurrentLocationAtomStation(int locationNumber) {
        return ((myChain).updResidue(Compound::Index(currentLocationResidue(locationNumber)))).getAtomLocationInMobilizedBodyFrame(currentLocationAtomIx(locationNumber));
    }
    int calcNextBodyIx (int locationNumber) {
        ((myChain).updResidue(Compound::Index(nextLocationResidue(locationNumber)))).getAtomMobilizedBodyIndex(nextAtomIx(locationNumber));
    }
    int calcBodyIx(locationNumber){
        ((myChain).updResidue(Compound::Index(currentLocationResidue(locationNumber)))).getAtomMobilizedBodyIndex(currentLocationAtomIx(locationNumber));
        
    }
    
    DuMM::ClusterIndex calcClusterIx(int locationNumber) {
            stringstream ss3;
            ss3.clear();
            ss3<<currentLocationResidue (locationNumber)<<  "/"<< currentLocationAtomName(locationNumber);
            return      wkdpForces.createCluster((ss3.str()).c_str());
    } 
    for ( locationNumber = 0; locationNumber < (2*((myChain).getNumResidues())-1) ; locationNumber++)    
    {

        if (locationNumber < (2*((myChain).getNumResidues())-2))
            nextBodyIx = calcNextBodyIx(locationNumber);           
        if (calcBodyIx(locationNumber) !=  oldBodyIx) {
       	    clusterIx = calcClusterIx(locationNumber);
            wkdpForces.placeAtomInCluster (currentLocationAtomIx(locationNumber),clusterIx,calcCurrentLocationAtomStation(locationNumber));
            if ((locationNumber == (2*(myChain).getNumResidues()-2)) || (calcNextBodyIx(locationNumber) != calcBodyIx(locationNumber) )) {
                wkdpForces.attachClusterToBody(clusterIx, calcBodyIx(locationNumber),Transform(Vec3(0)));
                if (myParameterReader.verbose) cout<<"[AddWKDPForces.h] attaching cluster to non-rigid body"<<endl;
            }
        } else {
            if (myParameterReader.verbose) cout<<"[AddWKDPForces.h] rigid cluster detected"<<endl;
            wkdpForces.placeAtomInCluster (currentLocationAtomIx(locationNumber) ,oldClusterIx, calcCurrentLocationAtomStation(locationNumber));
            if ((locationNumber == (2*(myChain).getNumResidues()-2)) || (calcNextBodyIx(locationNumber != calcBodyIx(locationNumber) )){
                wkdpForces.attachClusterToBody(oldClusterIx,oldBodyIx,Transform(Vec3(0)));
                if (myParameterReader.verbose) cout<<"[AddWKDPForces.h] attaching cluster to rigid body"<<endl;
            }
	}
        if (previousWKDPAtom.isValid() && (currentLocationAtomIx(locationNumber)).isValid()) {
            wkdpForces.addBond(previousWKDPAtom,currentLocationAtomIx(locationNumber));
        oldBodyIx = bodyIx;
        oldClusterIx = clusterIx; 
        }
        previousWKDPAtom = currentLocationAtomIx(locationNumber);
    }


*/
    DuMM::AtomIndex previousWKDPAtom;
    MobilizedBodyIndex bodyIx ;
    MobilizedBodyIndex nextBodyIx ;
    MobilizedBodyIndex oldBodyIx = myChain.getAtomMobilizedBodyIndex(Compound::AtomIndex(0));
    DuMM::ClusterIndex clusterIx  = wkdpForces.createCluster("0/C3*");
    DuMM::ClusterIndex oldClusterIx  = wkdpForces.createCluster("0/C3*");
    stringstream ss3;

    for (ResidueInfo::Index resIx(0); resIx < (myChain).getNumResidues(); ++resIx) {
        
        if (resIx < ((myChain).getNumResidues()-1)){
            const ResidueInfo & nextResidue = (myChain).getResidue(ResidueInfo::Index(int(resIx)+1));
            Compound::AtomIndex nextC3Ix = nextResidue.getAtomIndex("C3*");
            nextBodyIx = myChain.getAtomMobilizedBodyIndex(nextC3Ix);           
        }
        const  ResidueInfo & residue = (myChain).updResidue(resIx);
        if (myParameterReader.verbose) cout<<"[AddWKDPForces.h] resIx ="<<resIx<<endl;
        Compound::AtomIndex c3Ix = residue.getAtomIndex("C3*");
        if (myParameterReader.verbose) cout<<"[AddWKDPForces.h] c3Ix ="<<c3Ix<<endl;
        bodyIx = myChain.getAtomMobilizedBodyIndex(c3Ix);
        //if (bodyIx == oldBodyIx)  
	  //  bodyIx ;
        if (myParameterReader.verbose) cout<<"[AddWKDPForces.h] bodyIx ="<<bodyIx<<endl;
        Vec3 atomStation =  myChain.getAtomLocationInMobilizedBodyFrame(c3Ix);
        if (myParameterReader.verbose) cout<<"[AddWKDPForces.h] atomStation ="<<atomStation<<endl;
        if (myParameterReader.verbose) cout<<"[AddWKDPForces.h] wkdpChargedAtomTypeIx ="<<wkdpChargedAtomTypeIx<<endl;
        DuMM::AtomIndex wkdpAtom = wkdpForces.addAtom(wkdpChargedAtomTypeIx);
        if (bodyIx !=  oldBodyIx) {
            ss3.clear();
	    ss3<<resIx<<"/C3*";
       	    clusterIx = wkdpForces.createCluster((ss3.str()).c_str());
            wkdpForces.placeAtomInCluster (wkdpAtom,clusterIx,atomStation);
            if ((resIx == ((myChain).getNumResidues()-1)) || (nextBodyIx != bodyIx )) {
                wkdpForces.attachClusterToBody(clusterIx,bodyIx,Transform(Vec3(0)));
                if (myParameterReader.verbose) cout<<"[AddWKDPForces.h] attaching cluster to non-rigid body"<<endl;
            }
        } else {
            if (myParameterReader.verbose) cout<<"[AddWKDPForces.h] rigid cluster detected"<<endl;
            wkdpForces.placeAtomInCluster (wkdpAtom,oldClusterIx,atomStation);
            if ((resIx == ((myChain).getNumResidues()-1)) || (nextBodyIx != bodyIx )){
                wkdpForces.attachClusterToBody(oldClusterIx,oldBodyIx,Transform(Vec3(0)));
                if (myParameterReader.verbose) cout<<"[AddWKDPForces.h] attaching cluster to rigid body"<<endl;
            }
	}
                
        //wkdpForces.attachAtomToBody(wkdpAtom, bodyIx, atomStation);
        if (previousWKDPAtom.isValid() && wkdpAtom.isValid()) {
            wkdpForces.addBond(previousWKDPAtom, wkdpAtom);
        oldBodyIx = bodyIx;
        oldClusterIx = clusterIx; 
        }
        previousWKDPAtom = wkdpAtom;
   }


    return 0;
}

};
#endif
