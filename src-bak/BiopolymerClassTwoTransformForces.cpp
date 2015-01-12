/* -------------------------------------------------------------------------- *
 *                           MMB (MacroMoleculeBuilder)                       *
 * -------------------------------------------------------------------------- *
 *                                                                            *
 * Copyright (c) 2011-12 by the Author.                                       *
 * Author: Samuel Flores                                                      *
 *                                                                            *
 * See RNABuilder.cpp for the copyright and usage agreement.                  *
 * -------------------------------------------------------------------------- */
#include "BiopolymerClassTwoTransformForces.h"  

    AllTwoTransformLinearSprings::AllTwoTransformLinearSprings (SimbodyMatterSubsystem& matter,ParameterReader& myParameterReader,  LeontisWesthofClass& myLeontisWesthofClass, BiopolymerClassContainer & myBiopolymerClassContainer, std::ostream& outputStream ) : matter(matter),myParameterReader(myParameterReader), myLeontisWesthofClass (myLeontisWesthofClass), myBiopolymerClassContainer(myBiopolymerClassContainer), outputStream(outputStream)
        { 
    };    
    void         AllTwoTransformLinearSprings::calcAxes (const State& state,LeontisWesthofBondRow myLeontisWesthofBondRow,ResidueID residueNumber1,ResidueID residueNumber2,String chain1,String chain2,Vec3 & xAxisVector1,Vec3 & yAxisVector1, Vec3 & zAxisVector1,Vec3 & xAxisVector2,Vec3 & yAxisVector2 , Vec3 & zAxisVector2,Vec3 & glycosidicNitrogenAtom1LocationInGround,Vec3 & glycosidicNitrogenAtom2LocationInGround, Vec3 & ring1CenterLocationInGround, Vec3 & ring2CenterLocationInGround) const {
 
            glycosidicNitrogenAtom1LocationInGround = myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chain1, residueNumber1,myLeontisWesthofBondRow.residue1Atom[0]);
            glycosidicNitrogenAtom2LocationInGround = myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chain2, residueNumber2,myLeontisWesthofBondRow.residue2Atom[0]);
            Vec3 firstRingAtomvector1 = myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chain1,residueNumber1,myLeontisWesthofBondRow.residue1Atom[1])  - glycosidicNitrogenAtom1LocationInGround;
            Vec3 secondRingAtomvector1 = myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chain1,residueNumber1,myLeontisWesthofBondRow.residue1Atom[2])  - glycosidicNitrogenAtom1LocationInGround;
            Vec3 firstRingAtomvector2 = myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chain2,residueNumber2,myLeontisWesthofBondRow.residue2Atom[1])  - glycosidicNitrogenAtom2LocationInGround;
            //Vec3 firstRingAtomvector2 = myBiopolymerClass[j].calcAtomLocationInGroundFrame(state,myBiopolymerClass[j].getAtomIndex(ss2first.str()))  -glycosidicNitrogenAtom2LocationInGround;// myBiopolymerClass[j].calcAtomLocationInGroundFrame(state,myBiopolymerClass[j].getAtomIndex(ss4.str()));
            Vec3 secondRingAtomvector2 = myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chain2,residueNumber2,myLeontisWesthofBondRow.residue2Atom[2])  - glycosidicNitrogenAtom2LocationInGround;
            //Vec3 secondRingAtomvector2 = myBiopolymerClass[j].calcAtomLocationInGroundFrame(state,myBiopolymerClass[j].getAtomIndex(ss2second.str()))  - glycosidicNitrogenAtom2LocationInGround;//myBiopolymerClass[j].calcAtomLocationInGroundFrame(state,myBiopolymerClass[j].getAtomIndex(ss4.str()));
            if ((myLeontisWesthofBondRow.pdbResidueName1.compare("A  ") == 0) || (myLeontisWesthofBondRow.pdbResidueName1.compare("G  ") == 0) ) { //if purine

                xAxisVector1 =  -5.88327 * firstRingAtomvector1 - 6.13617 * secondRingAtomvector1;          
                ring1CenterLocationInGround = (myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chain1,residueNumber1,String("N3"))
                                              +myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chain1,residueNumber1,String("C6")))/2;
            }
            else if ((myLeontisWesthofBondRow.pdbResidueName1.compare("C  ") == 0)) {
                xAxisVector1 = -7.83435 * firstRingAtomvector1 -6.99265          *secondRingAtomvector1;           
                ring1CenterLocationInGround = (myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chain1,residueNumber1,String("N1"))
                                              +myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chain1,residueNumber1,String("C4")))/2;
            }
            else if ((myLeontisWesthofBondRow.pdbResidueName1.compare("U  ")) == 0) {
                xAxisVector1 = -7.3491 * firstRingAtomvector1 -6.47606 *secondRingAtomvector1;     
                ring1CenterLocationInGround = (myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chain1,residueNumber1,String("N1"))
                                              +myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chain1,residueNumber1,String("C4")))/2;
            } 
            else { cout <<__FILE__<<":"<<__LINE__<<"  Unrecognized residue type"<<endl; assert(0);} // trap errors
            if ((myLeontisWesthofBondRow.pdbResidueName2.compare("A  ")  == 0) || (myLeontisWesthofBondRow.pdbResidueName2.compare("G  ") == 0)){ //if purine
                xAxisVector2 = -5.88327 * firstRingAtomvector2 -6.13617 *secondRingAtomvector2;        
                ring2CenterLocationInGround = (myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chain2,residueNumber2,String("N3"))
                                              +myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chain2,residueNumber2,String("C6")))/2;
            }   
            else if ((myLeontisWesthofBondRow.pdbResidueName2.compare("C  ") == 0)){
                ring2CenterLocationInGround = (myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chain2,residueNumber2,String("N1"))
                                              +myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chain2,residueNumber2,String("C4")))/2;
                xAxisVector2 = -7.83435 * firstRingAtomvector2 -6.99265 *secondRingAtomvector2;           
            }
            else if ((myLeontisWesthofBondRow.pdbResidueName2.compare("U  ")) == 0) {
                xAxisVector2 = -7.3491  * firstRingAtomvector2 -6.47606 *secondRingAtomvector2;           
                ring2CenterLocationInGround = (myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chain2,residueNumber2,String("N1"))
                                              +myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chain2,residueNumber2,String("C4")))/2;
                if (myParameterReader.verbose) cout <<__FILE__<<":"<<__LINE__<<"  ring2CenterLocationInGround just computed for U: "<<ring2CenterLocationInGround<<endl;
            }
            else { cout <<__FILE__<<":"<<__LINE__<<"  Unrecognized residue type"<<endl; assert(0);} // trap errors
            zAxisVector1 = (firstRingAtomvector1%secondRingAtomvector1);
            zAxisVector1 = zAxisVector1/zAxisVector1.norm(); 
            zAxisVector2 = (firstRingAtomvector2%secondRingAtomvector2);
            zAxisVector2 = zAxisVector2/zAxisVector2.norm(); 
            yAxisVector1 = zAxisVector1%xAxisVector1;
            yAxisVector1= yAxisVector1/yAxisVector1.norm();
            yAxisVector2 = zAxisVector2%xAxisVector2;
            yAxisVector2= yAxisVector2/yAxisVector2.norm();

    };
    int AllTwoTransformLinearSprings::isThisATwoTransformForce(String myBPEdge) const {
        if  (((myBPEdge).compare("WatsonCrick") == 0) ||
            ((myBPEdge).compare("Hoogsteen"  ) == 0) ||
            ((myBPEdge).compare("SugarEdge"  ) == 0) ||
            ((myBPEdge).compare("Stacking"   ) == 0) ||
            ((myBPEdge).compare("AntiParallelStacking") == 0) ||
            ((myBPEdge).compare("HelicalStacking") == 0) ||
            ((myBPEdge).compare("Custom"     ) == 0) )
            return 1; 
        else 
            return 0;

    };
    void AllTwoTransformLinearSprings::calcForce(const State& state, Vector_<SpatialVec>& bodyForces,  
            Vector_<Vec3>& particleForces, Vector& mobilityForces) const 
        {  
        MobilizedBody body1;
        MobilizedBody body2;
        Transform transform1;
        Transform transform2;
        double forceConstant;
        double torqueConstant;
        double dutyCycle; //must be between 0 and 1.  at 1, force is applied all the time.  at 0, basically never applied.    
        double scrubberPeriod;
        double cutoffRadius;
        for (int r=0;r<myParameterReader.basePairContainer.numBasePairs();r++) 
        { 
	    if (myParameterReader.verbose) cout<<__FILE__<<":"<<__LINE__<<"  doing base pair #"<<r<<endl;	
            String chainId1=(myParameterReader.basePairContainer.getBasePair(r)).FirstBPChain;            
            String chainId2=(myParameterReader.basePairContainer.getBasePair(r)).SecondBPChain;
            String bondingEdge1=(myParameterReader.basePairContainer.getBasePair(r)).FirstBPEdge;
            String bondingEdge2=(myParameterReader.basePairContainer.getBasePair(r)).SecondBPEdge;
            String glycosidicBondOrientation=(myParameterReader.basePairContainer.getBasePair(r)).OrientationBP ;
            String basePairIsTwoTransformForce="baseInteraction";//(myParameterReader.basePairContainer.getBasePair(r)).BasePairIsTwoTransformForce ;

            ResidueID residueNumber1=(myParameterReader.basePairContainer.getBasePair(r).FirstBPResidue);
            ResidueID residueNumber2=(myParameterReader.basePairContainer.getBasePair(r).SecondBPResidue );
            String myPdbResidueName1 = ((myBiopolymerClassContainer.getPdbResidueName(chainId1,residueNumber1)));  
            String myPdbResidueName2 = ((myBiopolymerClassContainer.getPdbResidueName(chainId2,residueNumber2)));  
            LeontisWesthofBondRow myLeontisWesthofBondRow = myLeontisWesthofClass.myLeontisWesthofBondMatrix.myLeontisWesthofBondRow[(myParameterReader.basePairContainer.getBasePair(r)).leontisWesthofBondRowIndex];
            
            body1 = myBiopolymerClassContainer.updAtomMobilizedBody(matter,chainId1,residueNumber1,myLeontisWesthofBondRow.residue1Atom[0]);
            body2 = myBiopolymerClassContainer.updAtomMobilizedBody(matter,chainId2,residueNumber2,myLeontisWesthofBondRow.residue2Atom[0]);
            Vec3 station1 = myLeontisWesthofBondRow.attachmentPoint;
            Vec3 station2 = Vec3(0);
            Rotation rotation1 = Rotation(myLeontisWesthofBondRow.rotationAngle,myLeontisWesthofBondRow.rotationAxis);
            Rotation rotation2;
            rotation2.setRotationToIdentityMatrix();
            forceConstant = myLeontisWesthofBondRow.springConstant[0]; 
	    torqueConstant = myLeontisWesthofBondRow.torqueConstant;
	    dutyCycle = myParameterReader.dutyCycle;
            scrubberPeriod = myParameterReader.scrubberPeriod;
            cutoffRadius = myParameterReader.cutoffRadius; 
	    transform1 = Transform(rotation1, station1);    
            transform2 = Transform(rotation2, station2);    
            if (fmod(state.getTime(),scrubberPeriod)/scrubberPeriod >= (1-dutyCycle)) 
            { 
            Vec3 xAxisVector1 ;
	    Vec3 yAxisVector1;
	    Vec3 zAxisVector1;
 	    Vec3 xAxisVector2;
 	    Vec3 yAxisVector2;
  	    Vec3 zAxisVector2;
            Vec3 glycosidicNitrogenAtom1LocationInGround;
            Vec3 glycosidicNitrogenAtom2LocationInGround;
            Vec3 ring1CenterLocationInGround;
            Vec3 ring2CenterLocationInGround;
            Transform myX_GB1(Vec3(0));
            Transform myX_GB2(Vec3(0));
            //{ 
                //new method that uses pre-computed corrections for topology
                
            myX_GB1 =  matter.getMobilizedBody(body1).getBodyTransform(state);
            myX_GB2 =  matter.getMobilizedBody(body2).getBodyTransform(state); 
      
            myX_GB1 = Transform(~(myParameterReader.basePairContainer.getBasePair(r).rotationCorrection1*~myX_GB1.R()),myX_GB1.R()*myParameterReader.basePairContainer.getBasePair(r).translationCorrection1+myX_GB1.T());
            myX_GB2 = Transform(~(myParameterReader.basePairContainer.getBasePair(r).rotationCorrection2*~myX_GB2.R()),myX_GB2.R()*myParameterReader.basePairContainer.getBasePair(r).translationCorrection2+myX_GB2.T());
            //}
            const Transform X_GB1 = myX_GB1;
            const Transform X_GB2 = myX_GB2;
            

            const Vec3 s1_G = X_GB1.R() * station1;
            const Vec3 s2_G = X_GB2.R() * station2;
            const Rotation rot1_G = X_GB1.R() * rotation1;
            const Rotation rot2_G = X_GB2.R() * rotation2;
            const Vec3 p1_G = X_GB1.T() + s1_G; // station measured from ground origin
            const Vec3 p2_G = X_GB2.T() + s2_G;

            const Vec3 r_G = p2_G - p1_G; // vector from point1 to point2
            const Real d   = r_G.norm();  // distance between the points
            if (myParameterReader.verbose) cout<<__FILE__<<":"<<__LINE__<<" : d         = "<<d        <<endl;	
            Real  stretch   = d ; //d - x0; // + -> tension, - -> compression
	    const Vec4 rotationAngleAxis = (rot1_G*(~rot2_G)).convertRotationToAngleAxis();	
 
            Vec3 torque;

            double myFrcScalar;
            double A, B, C;
            if (myParameterReader.verbose) cout<<__FILE__<<":"<<__LINE__<<"  myParameterReader.potentialType :"<<myParameterReader.potentialType<<endl;
            for (int i = 0; i<3; i++) torque[i] = rotationAngleAxis[i+1];

            if (myLeontisWesthofBondRow.isTwoTransformForce.compare("aromatic")==0) 
                {
                cout<<__FILE__<<":"<<__LINE__<< "(rot1_G)     :  "<<((rot1_G))<<endl;
                cout<<__FILE__<<":"<<__LINE__<< "(rot2_G)     :  "<<((rot2_G))<<endl;
                cout<<__FILE__<<":"<<__LINE__<< "(rot1_G*(~rot2_G))     :  "<<(rot1_G*(~rot2_G))<<endl;
                cout<<__FILE__<<":"<<__LINE__<< "(rot1_G*(~rot2_G)).convertRotationToAngleAxis()     :  "<<(rot1_G*(~rot2_G)).convertRotationToAngleAxis()<<endl;
                cout<<__FILE__<<":"<<__LINE__<< "(rot1_G).x() :  "<<((rot1_G).x())<<endl;
                cout<<__FILE__<<":"<<__LINE__<< "(rot1_G).y() :  "<<((rot1_G).y())<<endl;
                cout<<__FILE__<<":"<<__LINE__<< "(rot1_G).z() :  "<<((rot1_G).z())<<endl;
                cout<<__FILE__<<":"<<__LINE__<< "torque       :  "<<torque<<endl;
                cout<<__FILE__<<":"<<__LINE__<< "dot(((rot1_G).x()),torque) : "<<dot(((rot1_G).x()),torque) <<endl;
                cout<<__FILE__<<":"<<__LINE__<< "dot(((rot1_G).y()),torque) : "<<dot(((rot1_G).y()),torque) <<endl;
                cout<<__FILE__<<":"<<__LINE__<< "dot(((rot1_G).z()),torque) : "<<dot(((rot1_G).z()),torque) <<endl;
                torque = torque - (dot(((rot2_G).z()),torque)*(((rot2_G).z())));   
                cout<<__FILE__<<":"<<__LINE__<< "torque       :  "<<torque<<endl;
                cout<<__FILE__<<":"<<__LINE__<< "dot(((rot1_G).z()),torque) : "<<dot(((rot1_G).z()),torque) <<endl;
                } 

            double theta = rotationAngleAxis[0];
            torque *= -theta * torqueConstant ;

            if ((myParameterReader.potentialType).compare("HarmonicInverseLorentzian") == 0)
                {if (stretch < cutoffRadius) 
                     {
                     C =  forceConstant*cutoffRadius;
                     A = -forceConstant/2/cutoffRadius/cutoffRadius;//    -C/2/cutoffRadius/cutoffRadius/cutoffRadius;
                     B = +C/cutoffRadius-A*cutoffRadius*cutoffRadius;
                     myFrcScalar = (double)(2*A*stretch); // overall positive quantity
                     myFrcScalar += (double)(theta*theta*torqueConstant*stretch/pow((1+pow((stretch/cutoffRadius),2)),2));
                     } 
                 else {
                     C = forceConstant*cutoffRadius;  
                     myFrcScalar = (double)(-C/stretch/stretch); //overall positive quantity.  Must be positive so it points towards body 2.
                     myFrcScalar += (double)(theta*theta*torqueConstant*stretch/pow((1+pow((stretch/cutoffRadius),2)),2));
                 }
                 torque *= 1/(1+pow((stretch/cutoffRadius),2));
                }
            else if ((myParameterReader.potentialType).compare("Harmonic" ) == 0) 
                {if (stretch > cutoffRadius) stretch =  0;myFrcScalar = (double)(forceConstant*stretch);  }
            else if ((myParameterReader.potentialType).compare("HarmonicLinear" ) == 0) 
                {if (stretch > cutoffRadius) stretch =  cutoffRadius;myFrcScalar = (double)(forceConstant*stretch);  }
            else if ((myParameterReader.potentialType).compare("Inverse") == 0) 
                {if (stretch < cutoffRadius) {myFrcScalar = 0;}
                 else {double C = forceConstant*cutoffRadius;  myFrcScalar = -C/stretch/stretch; }  // C should be negatve 
                }
            else if ((myParameterReader.potentialType).compare("HarmonicInverse") == 0)
                {if (stretch < cutoffRadius) 
                     {
                      C =  forceConstant*cutoffRadius;
                      A = -forceConstant/2/cutoffRadius/cutoffRadius;//    -C/2/cutoffRadius/cutoffRadius/cutoffRadius;
                      B = +C/cutoffRadius-A*cutoffRadius*cutoffRadius;
                      if (myParameterReader.verbose) cout<< __FILE__<<":"<<__LINE__  <<" A,C :"<<A<<","<<C<<endl;
                      myFrcScalar = 2*A*stretch*(1 + theta*theta*torqueConstant/2/forceConstant);
                      torque *= (A*stretch*stretch + B) /forceConstant; //*(-torqueConstant)
                     } 
                 else {
                     C = forceConstant*cutoffRadius;  
                     myFrcScalar = -C/stretch/stretch*(1 + theta*theta*torqueConstant/2/forceConstant);
                     A =  -forceConstant/2/cutoffRadius/cutoffRadius;
                     B =  +C/cutoffRadius-A*cutoffRadius*cutoffRadius;
                     if (myParameterReader.verbose) cout<< __FILE__<<":"<<__LINE__  <<" C :"<<C<<endl;

                     torque *= (C/stretch)/forceConstant;
                 }
                 if (myParameterReader.verbose) cout<< __FILE__<<":"<<__LINE__  <<" stretch, myFrcScalar: "<<stretch<<","<<myFrcScalar<<endl;
                }
            else if ((myParameterReader.potentialType).compare("HarmonicInverseCube") == 0)
                {if (stretch < cutoffRadius) 
                     {
                      C = forceConstant*cutoffRadius*cutoffRadius*cutoffRadius;
                      A   =  -3*forceConstant/2/cutoffRadius/cutoffRadius;//    -C/2/cutoffRadius/cutoffRadius/cutoffRadius;
                      //double B =      +C/cutoffRadius-A*cutoffRadius*cutoffRadius;
                      if (myParameterReader.verbose) cout<< __FILE__<<":"<<__LINE__  <<" A,C :"<<A<<","<<C<<endl;
                      myFrcScalar = 2*A*stretch;
                     } 
                 else {
                     C = forceConstant*cutoffRadius*cutoffRadius*cutoffRadius;  myFrcScalar = -3*C/stretch/stretch/stretch/stretch;
                     if (myParameterReader.verbose) cout<< __FILE__<<":"<<__LINE__  <<" C :"<<C<<endl;
                 }
                 if (myParameterReader.verbose) cout<< __FILE__<<":"<<__LINE__  <<" stretch, myFrcScalar: "<<stretch<<","<<myFrcScalar<<endl;
                }  else {assert (0);}
                
                myFrcScalar *= myParameterReader.twoTransformForceMultiplier;
                torque      *= myParameterReader.twoTransformForceMultiplier;

            const Real frcScalar = myFrcScalar;
            const Vec3 f1_G = (frcScalar/d) * r_G;
            if (myParameterReader.verbose) cout<<__FILE__<<":"<<__LINE__<<" : frcScalar = "<<frcScalar<<endl;	
	    if (myParameterReader.verbose) {cout<<__FILE__<<":"<<__LINE__<<"  torque on body 1 ="<<torque + s1_G % f1_G<<endl;}
	    if (myParameterReader.verbose) {cout<<__FILE__<<":"<<__LINE__<<"  torque on body 2 ="<<-torque - s2_G % f1_G<<endl;}
            //p1_G is station measured from ground origin 
            bodyForces[body1.getMobilizedBodyIndex()] +=  SpatialVec(torque + (-(matter.getMobilizedBody(body1).getBodyTransform(state)).T()+p1_G) % f1_G, f1_G);
            bodyForces[body2.getMobilizedBodyIndex()] -=  SpatialVec(torque + (-(matter.getMobilizedBody(body2).getBodyTransform(state)).T()+p2_G) % f1_G, f1_G);
            //cout<<__FILE__<<":"<<__LINE__<<" for mobilized body index "<<body1.getMobilizedBodyIndex()<<", just added  "<< SpatialVec(torque + (-(matter.getMobilizedBody(body1).getBodyTransform(state)).T()+p1_G) % f1_G, f1_G)<<endl;
            //std::cout<<__FILE__<<":"<<__LINE__<<": getForceIsDisabled () is  "<<   isForceDisabled( state,  index)<<  std::endl;
            }
        }
        };
    Real AllTwoTransformLinearSprings::calcPotentialEnergy(const State& state) const { 
        double energy = 0.0; 
   
        MobilizedBody body1;
        MobilizedBody body2;
        Transform transform1;
        Transform transform2;
        double forceConstant;
        double torqueConstant;
        double dutyCycle; //must be between 0 and 1.  at 1, force is applied all the time.  at 0, basically never applied.    
        double scrubberPeriod;
        double cutoffRadius;
        //ParameterReader tempParameterReader =  myParameterReader;          

        double interactingBaseRMSD=0; 
             
        int   interactingBaseRMSDCount = 0;              
        if (myParameterReader.verbose) cout<< __FILE__<<":"<<__LINE__  <<" energy before calcPotentialEnergy  :: "<< std::setiosflags(std::ios::fixed) << std::setprecision(1) << energy    <<","<<endl;
        myParameterReader.satisfiedBasePairs = 0;
        myParameterReader.unSatisfiedBasePairs = 0;
        for (int r=0;r<myParameterReader.basePairContainer.numBasePairs();r++)
        //if (( (myParameterReader.basePairContainer.getBasePair(r)).BasePairIsTwoTransformForce.compare("baseInteraction") ==0) ||
        //    ( (myParameterReader.basePairContainer.getBasePair(r)).BasePairIsTwoTransformForce.compare("aromatic"         ) ==0))   
        { 


 
            String chainId1=(myParameterReader.basePairContainer.getBasePair(r)).FirstBPChain;            
            String chainId2=(myParameterReader.basePairContainer.getBasePair(r)).SecondBPChain;
            String bondingEdge1=(myParameterReader.basePairContainer.getBasePair(r)).FirstBPEdge;
            String bondingEdge2=(myParameterReader.basePairContainer.getBasePair(r)).SecondBPEdge;
            String glycosidicBondOrientation=(myParameterReader.basePairContainer.getBasePair(r)).OrientationBP ;
            String basePairIsTwoTransformForce="baseInteraction";// (myParameterReader.basePairContainer.getBasePair(r)).BasePairIsTwoTransformForce ;


            ResidueID residueNumber1=(myParameterReader.basePairContainer.getBasePair(r)).FirstBPResidue;
            ResidueID residueNumber2=(myParameterReader.basePairContainer.getBasePair(r)).SecondBPResidue;
            String myPdbResidueName1 = myBiopolymerClassContainer.getPdbResidueName(chainId1, residueNumber1);  
            String myPdbResidueName2 = myBiopolymerClassContainer.getPdbResidueName(chainId2, residueNumber2);  
            LeontisWesthofBondRow myLeontisWesthofBondRow = myLeontisWesthofClass.getLeontisWesthofBondRow((myParameterReader.basePairContainer.getBasePair(r)).FirstBPResidue ,(myParameterReader.basePairContainer.getBasePair(r)).SecondBPResidue,myPdbResidueName1,bondingEdge1,myPdbResidueName2,bondingEdge2,glycosidicBondOrientation,  basePairIsTwoTransformForce);
            stringstream ss3;
            ss3<<residueNumber1.outString()<<"/"<<myLeontisWesthofBondRow.residue1Atom[0];
            stringstream ss4;
            ss4<<residueNumber2.outString()<<"/"<<(myLeontisWesthofBondRow.residue2Atom[0]);
            body1 = myBiopolymerClassContainer.updAtomMobilizedBody(matter,chainId1,residueNumber1,myLeontisWesthofBondRow.residue1Atom[0]);
            body2 = myBiopolymerClassContainer.updAtomMobilizedBody(matter,chainId2,residueNumber2,myLeontisWesthofBondRow.residue2Atom[0]);
            Vec3 station1 = myLeontisWesthofBondRow.attachmentPoint;
            Vec3 station2 = myBiopolymerClassContainer.getAtomLocationInMobilizedBodyFrame(chainId2, residueNumber2,myLeontisWesthofBondRow.residue2Atom[0]);//myBiopolymerClass[j].getAtomLocationInMobilizedBodyFrame(myBiopolymerClass[j].getAtomIndex(ss4.str()));
            Rotation rotation1 = Rotation(myLeontisWesthofBondRow.rotationAngle,myLeontisWesthofBondRow.rotationAxis);
            Rotation rotation2;
            rotation2.setRotationToIdentityMatrix();
            forceConstant = myLeontisWesthofBondRow.springConstant[0]; 
	    torqueConstant = myLeontisWesthofBondRow.torqueConstant;
	    dutyCycle = myParameterReader.dutyCycle;
            scrubberPeriod = myParameterReader.scrubberPeriod;
            cutoffRadius = myParameterReader.cutoffRadius; 
	    transform1 = Transform(rotation1, station1);    
            transform2 = Transform(rotation2, station2);    
	   

            { 

                // still haven't added the torsional part of the energy.  
            Vec3 station1 = myLeontisWesthofBondRow.attachmentPoint;
            Vec3 station2 = Vec3(0);
            Vec3 xAxisVector1 ;
	    Vec3 yAxisVector1;
	    Vec3 zAxisVector1;
 	    Vec3 xAxisVector2;
 	    Vec3 yAxisVector2;
  	    Vec3 zAxisVector2;
            Vec3 glycosidicNitrogenAtom1LocationInGround;
            Vec3 glycosidicNitrogenAtom2LocationInGround;
            Vec3 ring1CenterLocationInGround;
            Vec3 ring2CenterLocationInGround;

            Transform myX_GB1;
            Transform myX_GB2;
            //{
                //new method that uses pre-computed corrections for topology
            myX_GB1 =  matter.getMobilizedBody(body1).getBodyTransform(state);
            myX_GB2 =  matter.getMobilizedBody(body2).getBodyTransform(state);
            myX_GB1 = Transform(myParameterReader.basePairContainer.getBasePair(r).rotationCorrection1*myX_GB1.R(),myX_GB1.R()*myParameterReader.basePairContainer.getBasePair(r).translationCorrection1+myX_GB1.T());
            myX_GB2 = Transform(myParameterReader.basePairContainer.getBasePair(r).rotationCorrection2*myX_GB2.R(),myX_GB2.R()*myParameterReader.basePairContainer.getBasePair(r).translationCorrection2+myX_GB2.T());
            //}
            const Transform X_GB1 = myX_GB1;
            const Transform X_GB2 = myX_GB2;

            const Vec3 s1_G = X_GB1.R() * station1;
            const Vec3 s2_G = X_GB2.R() * station2;
            const Vec3 p1_G = X_GB1.T() + s1_G; // station measured from ground origin
            const Vec3 p2_G = X_GB2.T() + s2_G;
            const Vec3 r_G = p2_G - p1_G; // vector from point1 to point2
            if (myParameterReader.verbose) cout<< __FILE__<<":"<<__LINE__  <<" r_G.norm(): "<<r_G.norm()<<","<<endl;

            if (myParameterReader.verbose) cout<< __FILE__<<":"<<__LINE__  <<" energy   :: "<<energy    <<","<<endl;

            if (myParameterReader.verbose) cout<<__FILE__<<":"<<__LINE__<<"  myParameterReader.potentialType :"<<myParameterReader.potentialType<<endl;
            if ((myParameterReader.potentialType).compare("HarmonicInverseLorentzian" ) == 0) {
                {
                const Rotation rot1_G = X_GB1.R() * rotation1;
                const Rotation rot2_G = X_GB2.R() * rotation2;
                const Vec4 rotationAngleAxis = (rot1_G*(~rot2_G)).convertRotationToAngleAxis();
                double theta = rotationAngleAxis[0];

                if (r_G.norm() < cutoffRadius)
                     {
                      double C = forceConstant*cutoffRadius;
                      double A =  -forceConstant/2/cutoffRadius/cutoffRadius;//    -C/2/cutoffRadius/cutoffRadius/cutoffRadius;
                      double B =  +C/cutoffRadius-A*cutoffRadius*cutoffRadius;
                      if (myParameterReader.verbose) cout<< __FILE__<<":"<<__LINE__  <<" A,C :"<<A<<","<<C<<endl;
                      energy      += (
                          B+ A*r_G.norm()*r_G.norm()
                          +theta*theta*torqueConstant/2/(1+pow((r_G.norm()/cutoffRadius),2))           
                          )*myParameterReader.twoTransformForceMultiplier;
                     }
                 else {
                     double C = forceConstant*cutoffRadius; 
                     energy += (C/r_G.norm()
                         +theta*theta*torqueConstant/2/(1+pow((r_G.norm()/cutoffRadius),2))           
                         )*myParameterReader.twoTransformForceMultiplier;
                 }
                }
      

            }
            else if ((myParameterReader.potentialType).compare("Harmonic" ) == 0)
                {if (r_G.norm() > cutoffRadius) energy +=  0;
                 else energy += forceConstant*(-cutoffRadius*cutoffRadius  + r_G.norm()*r_G.norm())/2*myParameterReader.twoTransformForceMultiplier;}
            else if ((myParameterReader.potentialType).compare("HarmonicLinear" ) == 0)
                {if (r_G.norm() > cutoffRadius) energy +=  (r_G.norm()-cutoffRadius)*forceConstant*myParameterReader.twoTransformForceMultiplier;
                 else energy += forceConstant*(-cutoffRadius*cutoffRadius  + r_G.norm()*r_G.norm())/2*myParameterReader.twoTransformForceMultiplier;}
            else if ((myParameterReader.potentialType).compare("HarmonicInverse") == 0)
                {if (r_G.norm() <= cutoffRadius)
                     {
                      double C = forceConstant*cutoffRadius;
                      double A   =  -forceConstant/2/cutoffRadius/cutoffRadius;//    -C/2/cutoffRadius/cutoffRadius/cutoffRadius;
                      double B =  +C/cutoffRadius-A*cutoffRadius*cutoffRadius;
                      if (myParameterReader.verbose) cout<< __FILE__<<":"<<__LINE__  <<" A,C :"<<A<<","<<C<<endl;
                      energy      += (B+ A*r_G.norm()*r_G.norm())*myParameterReader.twoTransformForceMultiplier;
                     if (myParameterReader.verbose) cout<< __FILE__<<":"<<__LINE__  <<" small separation, just added this much to energy:"<< B+ A*r_G.norm()*r_G.norm()  <<endl;
                      if (myParameterReader.verbose) cout<< __FILE__<<":"<<__LINE__  <<" r_G.norm(): "<<r_G.norm()<<","<<endl;
                     }
                 else {
                     double C = forceConstant*cutoffRadius; energy += C/r_G.norm()*myParameterReader.twoTransformForceMultiplier;
                     if (myParameterReader.verbose) cout<< __FILE__<<":"<<__LINE__  <<" large separation,  energy is now:"<< energy      <<endl;
                     if (myParameterReader.verbose) cout<< __FILE__<<":"<<__LINE__  <<" C :"<<C<<endl;
                     if (myParameterReader.verbose) cout<< __FILE__<<":"<<__LINE__  <<" r_G.norm(): "<<r_G.norm()<<","<<endl;
                 }
                 if (myParameterReader.verbose) cout<< __FILE__<<":"<<__LINE__  <<" r_G.norm(): "<<r_G.norm()<<","<<endl;
                }
       

            else if ((myParameterReader.potentialType).compare("HarmonicInverseCube") == 0)
                {if (r_G.norm() < cutoffRadius) 
                     {
                      double A = -3*forceConstant/2/cutoffRadius/cutoffRadius;
                      double B = 5/2*forceConstant;
                      energy      += (B+ A*r_G.norm()*r_G.norm())*myParameterReader.twoTransformForceMultiplier;
                     } 
                 else {
                     double C = forceConstant*cutoffRadius*cutoffRadius*cutoffRadius; 
                     energy += C/r_G.norm()/r_G.norm()/r_G.norm() *myParameterReader.twoTransformForceMultiplier;
                     if (myParameterReader.verbose) cout<< __FILE__<<":"<<__LINE__  <<" C :"<<C<<endl;
                 }
                } else assert(0); 

                if (myParameterReader.verbose) cout<< __FILE__<<":"<<__LINE__  <<" energy    : "<<energy    <<","<<endl;




		
// torsional part
            Rotation rotation1 = transform1.R();
       	    Rotation rotation2 = transform2.R();
            const Rotation rot1_G = X_GB1.R() * rotation1;
            const Rotation rot2_G = X_GB2.R() * rotation2;
	    const Vec4 rotationAngleAxis = (rot1_G*(~rot2_G)).convertRotationToAngleAxis();	
            if (r_G.norm() <= cutoffRadius) {energy += 0.5*((rotationAngleAxis[0]*rotationAngleAxis[0]-180*180/Rad2Deg/Rad2Deg)*torqueConstant)*myParameterReader.twoTransformForceMultiplier; 
                if (myParameterReader.verbose) cout<< __FILE__<<":"<<__LINE__  <<" rotationAngleAxis[0] = "<<rotationAngleAxis[0]<<endl;
                if (myParameterReader.verbose) cout<< __FILE__<<":"<<__LINE__  <<" calculating torsional energy. r_G.norm(),cutoffRadius: "<<r_G.norm()<<","<<cutoffRadius<<endl;
                if (myParameterReader.verbose) cout<< __FILE__<<":"<<__LINE__  <<" energy    : "<<energy    <<","<<endl;
                
            } // add nothing if r.norm()>cutoffRadius
            else {
            }
            
            //if (((((myParameterReader.basePairContainer.getBasePair(r))).BasePairIsTwoTransformForce).compare("baseInteraction") ==0) ) 
                { 
               
                outputStream <<"REMARK Base pair : "<< myPdbResidueName1<<" " <<(myParameterReader.basePairContainer.getBasePair(r)).FirstBPResidue.outString()<<" "<< bondingEdge1 <<", "  
                                    << myPdbResidueName2<<" " <<(myParameterReader.basePairContainer.getBasePair(r)).SecondBPResidue.outString()<<" "<< bondingEdge2 
                                    <<" "<<glycosidicBondOrientation<<".  ";//<<endl;
                outputStream <<"  Base pair "<<r<<" attachment-body frame distance = "<<r_G.norm() <<endl;
                interactingBaseRMSD += r_G.norm() *  r_G.norm() ;
                interactingBaseRMSDCount++;
                if (r_G.norm() <= cutoffRadius) {
                    myParameterReader.basePairContainer.setBasePairSatisfied(r,"True") ;
                    outputStream <<"REMARK "<<__FILE__<<" :  Base pair "<<r<<" IS satisfied"<<endl;
		    myParameterReader.satisfiedBasePairs ++;	


                } else {
                    myParameterReader.basePairContainer.setBasePairSatisfied(r,"False") ;
                    outputStream  <<"REMARK "<<__FILE__ <<" :  Base pair "<<r<<" is NOT satisfied"<<endl;
		    myParameterReader.unSatisfiedBasePairs ++;
                }
                }
            }    
        }
            
        cout << __FILE__<<":"<<__LINE__<<"  Satisfied baseInteraction's:   : "<<myParameterReader.satisfiedBasePairs<< " out of : "<<myParameterReader.satisfiedBasePairs+ myParameterReader.unSatisfiedBasePairs<<endl;    
        outputStream << "REMARK [TwoTransformForces.cpp] Satisfied base interactions   : "<<myParameterReader.satisfiedBasePairs<< " out of : "<<myParameterReader.satisfiedBasePairs+ myParameterReader.unSatisfiedBasePairs<<endl;    
        outputStream <<"REMARK RMSD for interacting frames = "<<  pow((interactingBaseRMSD / interactingBaseRMSDCount), double(0.5))<<endl;
        return energy; 
    };
    bool AllTwoTransformLinearSprings::dependsOnlyOnPositions() const  { 
        return true; 
    };    
