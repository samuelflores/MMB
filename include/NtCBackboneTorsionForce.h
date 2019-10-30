#ifndef KBBackboneTorsionForce_H_
#define KBBackboneTorsionForce_H_


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





class KBBackboneTorsion : public DuMM::CustomBondTorsion {
public:
   KBBackboneTorsion(
       float kbBackboneTorsionGlobalScaleFactor,
       float amplitude1, Angle phase1, Angle stdDeviation1,
       float amplitude2=0, Angle phase2=0, Angle stdDeviation2=0,
       float amplitude3=0, Angle phase3=0, Angle stdDeviation3=0,
       float amplitude4=0, Angle phase4=0, Angle stdDeviation4=0
   )
     : 
     kbBackboneTorsionGlobalScaleFactor(kbBackboneTorsionGlobalScaleFactor),
     amplitude1(amplitude1), phase1(phase1), stdDeviation1(stdDeviation1), 
     amplitude2(amplitude2), phase2(phase2), stdDeviation2(stdDeviation2), 
     amplitude3(amplitude3), phase3(phase3), stdDeviation3(stdDeviation3) ,
     amplitude4(amplitude4), phase4(phase4), stdDeviation4(stdDeviation4) 
     { }

   Real calcTorque(Angle dihedralAngle) const {
       Real   myTorque = 0;
       myTorque += amplitude1 * (-2*(dihedralAngle-phase1)/2/stdDeviation1/stdDeviation1) *exp (-2*(dihedralAngle-phase1)*(dihedralAngle-phase1)/2/stdDeviation1/stdDeviation1)  ;
       if (amplitude2>0) myTorque += amplitude2 * (-2*(dihedralAngle-phase2)/2/stdDeviation2/stdDeviation2) * exp (-2*(dihedralAngle-phase2)*(dihedralAngle-phase2)/2/stdDeviation2/stdDeviation2)  ;
       if (amplitude3>0) myTorque += amplitude3 * (-2*(dihedralAngle-phase3)/2/stdDeviation3/stdDeviation3) * exp (-2*(dihedralAngle-phase3)*(dihedralAngle-phase3)/2/stdDeviation3/stdDeviation3)  ;
       if (amplitude4>0) myTorque += amplitude4 * (-2*(dihedralAngle-phase4)/2/stdDeviation4/stdDeviation4) * exp (-2*(dihedralAngle-phase4)*(dihedralAngle-phase4)/2/stdDeviation4/stdDeviation4)  ;
       return myTorque*kbBackboneTorsionGlobalScaleFactor;
   }

   Real   calcEnergy(Angle dihedralAngle) const
   {
       float myEnergy = 0;
       myEnergy += amplitude1 * exp (-(dihedralAngle-phase1)*(dihedralAngle-phase1)/2/stdDeviation1/stdDeviation1);
       if (amplitude2>0) myEnergy += amplitude2 * exp (-(dihedralAngle-phase2)*(dihedralAngle-phase2)/2/stdDeviation2/stdDeviation2);
       if (amplitude3>0) myEnergy += amplitude3 * exp (-(dihedralAngle-phase3)*(dihedralAngle-phase3)/2/stdDeviation3/stdDeviation3);
       if (amplitude4>0) myEnergy += amplitude4 * exp (-(dihedralAngle-phase4)*(dihedralAngle-phase4)/2/stdDeviation4/stdDeviation4);
       return  myEnergy * kbBackboneTorsionGlobalScaleFactor;
   }

private:
  float kbBackboneTorsionGlobalScaleFactor;
  float amplitude1;
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
  int stdDeviation4;
};


class DefineKBBackboneTorsion {
public:
    DefineKBBackboneTorsion (DuMMForceFieldSubsystem & dumm, ParameterReader & myParameterReader )

    {

       dumm.defineCustomBondTorsion(
           SimTK::DuMM::AtomClassIndex(1 )  ,  SimTK::DuMM::AtomClassIndex(1) ,  SimTK::DuMM::AtomClassIndex(1) ,  SimTK::DuMM::AtomClassIndex( 23),
           new KBBackboneTorsion(myParameterReader.kbBackboneTorsionGlobalScaleFactor,1 ,25.0/Rad2Deg ,1,1,90/Rad2Deg,1,1,150/Rad2Deg,1,1,75/Rad2Deg,1 ) 
           );
       dumm.defineCustomBondTorsion(
           SimTK::DuMM::AtomClassIndex(23)  ,  SimTK::DuMM::AtomClassIndex(1) ,  SimTK::DuMM::AtomClassIndex(18) ,  SimTK::DuMM::AtomClassIndex( 2),
           new KBBackboneTorsion(myParameterReader.kbBackboneTorsionGlobalScaleFactor,1 ,30.0/Rad2Deg ,1 ,1,110./Rad2Deg,1) 
           );
       dumm.defineCustomBondTorsion(
           SimTK::DuMM::AtomClassIndex(1 )  ,  SimTK::DuMM::AtomClassIndex(1) ,  SimTK::DuMM::AtomClassIndex(23) ,  SimTK::DuMM::AtomClassIndex( 28),
           new KBBackboneTorsion(myParameterReader.kbBackboneTorsionGlobalScaleFactor,1 ,90.0/Rad2Deg ,1 ) 
           );
       dumm.defineCustomBondTorsion(
           SimTK::DuMM::AtomClassIndex(1 )  ,  SimTK::DuMM::AtomClassIndex(23) ,  SimTK::DuMM::AtomClassIndex(28) ,  SimTK::DuMM::AtomClassIndex( 23),
           new KBBackboneTorsion(myParameterReader.kbBackboneTorsionGlobalScaleFactor,1 ,140./Rad2Deg ,1,1,45/Rad2Deg,1 ) 
           );


    }

};
#endif
