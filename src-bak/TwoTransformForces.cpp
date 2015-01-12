#include "TwoTransformForces.h"  

    AllTwoTransformLinearSprings::AllTwoTransformLinearSprings (SimbodyMatterSubsystem& matter,ParameterReader& myParameterReader,  LeontisWesthofClass& myLeontisWesthofClass, BiopolymerClassContainer & myBiopolymerClassContainer, std::ostream& outputStream ) : matter(matter),myParameterReader(myParameterReader), myLeontisWesthofClass (myLeontisWesthofClass), myBiopolymerClassContainer(myBiopolymerClassContainer), outputStream(outputStream)
        { 
    };    

    void AllTwoTransformLinearSprings::realizeTopology(State& state) const {
	parameterReaderIndex = state.allocateDiscreteVariable(matter.getMySubsystemIndex(),  
                Stage::Dynamics, new Value<ParameterReader>(myParameterReader)); 
    };   
    void AllTwoTransformLinearSprings::setParameterReader(State& state, ParameterReader& myParameterReader) const { 
        static_cast<Value<ParameterReader>&>(state.updDiscreteVariable( matter.getMySubsystemIndex(),
                SimTK::DiscreteVariableIndex(parameterReaderIndex))).upd() = myParameterReader; 
    };

    const ParameterReader & AllTwoTransformLinearSprings::getParameterReader(const State& state) const { 
        return static_cast<const Value<ParameterReader>&>(state.getDiscreteVariable(matter.getMySubsystemIndex(), 
                SimTK::DiscreteVariableIndex(parameterReaderIndex))).get(); 
    };
    void         AllTwoTransformLinearSprings::calcAxes (const State& state,LeontisWesthofBondRow myLeontisWesthofBondRow,int residueNumber1,int residueNumber2,string chain1,string chain2,Vec3 & xAxisVector1,Vec3 & yAxisVector1, Vec3 & zAxisVector1,Vec3 & xAxisVector2,Vec3 & yAxisVector2 , Vec3 & zAxisVector2,Vec3 & glycosidicNitrogenAtom1LocationInGround,Vec3 & glycosidicNitrogenAtom2LocationInGround, Vec3 & ring1CenterLocationInGround, Vec3 & ring2CenterLocationInGround) const {
 
            stringstream ss3;
            ss3.clear();
            ss3<<"anything";
            ss3.str("");
            ss3<<residueNumber1;
            ss3<<"/";
            ss3<<myLeontisWesthofBondRow.residue1Atom[0];
            stringstream ss4(stringstream::in | stringstream::out);//(""); 
            ss4.str("");
            ss4.clear();
            ss4<<residueNumber2<<"/"<<(myLeontisWesthofBondRow.residue2Atom[0]);
          
            stringstream ss1first;
            ss1first<<residueNumber1<<"/"<<myLeontisWesthofBondRow.residue1Atom[1];
            stringstream ss1second;
            ss1second<<  residueNumber1<<"/"<<(myLeontisWesthofBondRow.residue1Atom[2]);
            stringstream ss1c1p   ;
            ss1c1p<<  residueNumber1<<"/C1*";//<<(myLeontisWesthofBondRow.residue1Atom[2]);
            stringstream ss2first;
            ss2first<<   residueNumber2<<"/"<<myLeontisWesthofBondRow.residue2Atom[1];
            stringstream ss2second;
            ss2second<<  residueNumber2<<"/"<<(myLeontisWesthofBondRow.residue2Atom[2]);
            stringstream ssRing1Atom1;
            stringstream ssRing1Atom2;
            stringstream ssRing2Atom1;
            stringstream ssRing2Atom2;
            ssRing1Atom1<<residueNumber1<<"/";
            ssRing1Atom2<<residueNumber1<<"/";
            ssRing2Atom1<<residueNumber2<<"/";
            ssRing2Atom2<<residueNumber2<<"/";
            glycosidicNitrogenAtom1LocationInGround = myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chain1,residueNumber1,myLeontisWesthofBondRow.residue1Atom[0]);
            glycosidicNitrogenAtom2LocationInGround = myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chain2,residueNumber2,myLeontisWesthofBondRow.residue2Atom[0]);

            Vec3 firstRingAtomvector1 = myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chain1, residueNumber1,myLeontisWesthofBondRow.residue1Atom[1]);  - glycosidicNitrogenAtom1LocationInGround;
            Vec3 secondRingAtomvector1 = myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chain1, residueNumber1,myLeontisWesthofBondRow.residue1Atom[2]);  - glycosidicNitrogenAtom1LocationInGround;
            Vec3 firstRingAtomvector2 = myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chain2, residueNumber2,myLeontisWesthofBondRow.residue2Atom[1]);  - glycosidicNitrogenAtom2LocationInGround;
            Vec3 secondRingAtomvector2 = myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chain2, residueNumber2,myLeontisWesthofBondRow.residue2Atom[2])  - glycosidicNitrogenAtom2LocationInGround;
            if ((myLeontisWesthofBondRow.pdbResidueName1.compare("A  ") == 0) || (myLeontisWesthofBondRow.pdbResidueName1.compare("G  ") == 0) ) { //if purine

                xAxisVector1 =  -5.88327 * firstRingAtomvector1 - 6.13617 * secondRingAtomvector1;          
                ring1CenterLocationInGround = (myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chain1,residueNumber1,string("N3"))
                                              +myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chain1,residueNumber1,string("C6")))/2;
            }
            else if ((myLeontisWesthofBondRow.pdbResidueName1.compare("C  ") == 0)) {
                xAxisVector1 = -7.83435 * firstRingAtomvector1 -6.99265          *secondRingAtomvector1;           
                ring1CenterLocationInGround = (myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chain1,residueNumber1,string("N1"))
                                              +myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chain1,residueNumber1,string("C4")))/2;
            }
            else if ((myLeontisWesthofBondRow.pdbResidueName1.compare("U  ")) == 0) {
                xAxisVector1 = -7.3491 * firstRingAtomvector1 -6.47606 *secondRingAtomvector1;     
                ring1CenterLocationInGround = (myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chain1,residueNumber1,string("N1"))
                                              +myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chain1,residueNumber1,string("C4")))/2;
            } 
            else { cout <<"[TwoTransformForces.cpp] Unrecognized residue type"<<endl; assert(0);} // trap errors
            if ((myLeontisWesthofBondRow.pdbResidueName2.compare("A  ")  == 0) || (myLeontisWesthofBondRow.pdbResidueName2.compare("G  ") == 0)){ //if purine
                xAxisVector2 = -5.88327 * firstRingAtomvector2 -6.13617 *secondRingAtomvector2;        
                ring2CenterLocationInGround = (myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chain2,residueNumber2,string("N3"))
                                              +myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chain2,residueNumber2,string("C6")))/2;
                
               
            }   
            else if ((myLeontisWesthofBondRow.pdbResidueName2.compare("C  ") == 0)){
                ring2CenterLocationInGround = (myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chain2,residueNumber2,string("N1"))
                                              +myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chain2,residueNumber2,string("C4")))/2;
                xAxisVector2 = -7.83435 * firstRingAtomvector2 -6.99265 *secondRingAtomvector2;           
            }
            else if ((myLeontisWesthofBondRow.pdbResidueName2.compare("U  ")) == 0) {
                xAxisVector2 = -7.3491  * firstRingAtomvector2 -6.47606 *secondRingAtomvector2;           
		ssRing2Atom1<<"N1";
		ssRing2Atom2<<"C4";
                if (myParameterReader.verbose) cout <<"[TwoTransformForces.cpp] ssRing2Atom1, ssRing2Atom2 : "<<ssRing2Atom1.str()<<", "<<ssRing2Atom2.str()<<                   endl;
                ring2CenterLocationInGround = (myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chain2,residueNumber2,string("N1"))
                                              +myBiopolymerClassContainer.calcAtomLocationInGroundFrame(state,chain2,residueNumber2,string("C4")))/2;
                if (myParameterReader.verbose) cout <<"[TwoTransformForces.cpp] ring2CenterLocationInGround just computed for U: "<<ring2CenterLocationInGround<<endl;
            }
            else { cout <<"[TwoTransformForces.cpp] Unrecognized residue type"<<endl; assert(0);} // trap errors
            zAxisVector1 = (firstRingAtomvector1%secondRingAtomvector1);
            zAxisVector1 = zAxisVector1/zAxisVector1.norm(); 
            zAxisVector2 = (firstRingAtomvector2%secondRingAtomvector2);
            zAxisVector2 = zAxisVector2/zAxisVector2.norm(); 
            yAxisVector1 = zAxisVector1%xAxisVector1;
            yAxisVector1= yAxisVector1/yAxisVector1.norm();
            yAxisVector2 = zAxisVector2%xAxisVector2;
            yAxisVector2= yAxisVector2/yAxisVector2.norm();

    };
    int AllTwoTransformLinearSprings::isThisATwoTransformForce(string myBPEdge) const {
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
        const ParameterReader & myParameterReader = getParameterReader(state);   
        MobilizedBody body1;
        MobilizedBody body2;
        Transform transform1;
        Transform transform2;
        float forceConstant;
        float torqueConstant;
        float dutyCycle; //must be between 0 and 1.  at 1, force is applied all the time.  at 0, basically never applied.    
        double scrubberPeriod;
        float cutoffRadius;
        for (int r=0;r<(myParameterReader.baseOperationVector).size();r++) 
        if  (((myParameterReader.baseOperationVector[r]).BasePairIsTwoTransformForce.compare("baseInteraction") == 0) ||
            ( (myParameterReader.baseOperationVector[r]).BasePairIsTwoTransformForce.compare("aromatic")          == 0))
            { 
	    if (myParameterReader.verbose) cout<<"[TwoTransformForces.cpp] doing base pair #"<<r<<endl;	
            string chainId1=(myParameterReader.baseOperationVector[r]).FirstBPChain;            
            string chainId2=(myParameterReader.baseOperationVector[r]).SecondBPChain;
            string bondingEdge1=(myParameterReader.baseOperationVector[r]).FirstBPEdge;
            string bondingEdge2=(myParameterReader.baseOperationVector[r]).SecondBPEdge;
            string glycosidicBondOrientation=(myParameterReader.baseOperationVector[r]).OrientationBP ;
            string basePairIsTwoTransformForce=(myParameterReader.baseOperationVector[r]).BasePairIsTwoTransformForce ;
            int i = -11111;
            int j = -11111;

            for (int k = 0; k < myParameterReader.numChains; k++) {
                if (strcmp((myParameterReader.chainId[k]).c_str(),chainId1.c_str()) == 0) i = k;
                if (strcmp((myParameterReader.chainId[k]).c_str(),chainId2.c_str()) == 0) j = k;
            }
            assert(i>=0);
            assert(j>=0);


            int residueNumber1=(myParameterReader.baseOperationVector[r]).FirstBPResidue-myParameterReader.getFirstResidueNumbers((myParameterReader.chainId[i]).c_str());            
            int residueNumber2=(myParameterReader.baseOperationVector[r]).SecondBPResidue-myParameterReader.getFirstResidueNumbers((myParameterReader.chainId[j]).c_str());
            string myPdbResidueName1 = myBiopolymerClassContainer.getPdbResidueName(chainId1, residueNumber1);  
            string myPdbResidueName2 = myBiopolymerClassContainer.getPdbResidueName(chainId2, residueNumber2);  
            //string myPdbResidueName2 = ((myBiopolymerClassContainer[j]).updResidue(ResidueInfo::Index(residueNumber2))).getPdbResidueName();  
            LeontisWesthofBondRow myLeontisWesthofBondRow = myLeontisWesthofClass.myLeontisWesthofBondMatrix.myLeontisWesthofBondRow[(myParameterReader.baseOperationVector[r]).leontisWesthofBondRowIndex];
            stringstream ss3;
            ss3<<residueNumber1<<"/"<<myLeontisWesthofBondRow.residue1Atom[0];
            stringstream ss4;
            ss4<<residueNumber2<<"/"<<(myLeontisWesthofBondRow.residue2Atom[0]);
            
            body1 = myBiopolymerClassContainer.updAtomMobilizedBody(matter, chainId1, residueNumber1,myLeontisWesthofBondRow.residue1Atom[0]);
            body2 = myBiopolymerClassContainer.updAtomMobilizedBody(matter, chainId2, residueNumber2,myLeontisWesthofBondRow.residue2Atom[0]);
            //body2 = matter.updMobilizedBody(myBiopolymerClassContainer[j].getAtomMobilizedBodyIndex(Compound::AtomIndex(myBiopolymerClassContainer[j].getAtomIndex(ss4.str()))));     
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
            //cout<<"residueNumber1,residueNumber2 : "<<residueNumber1<<","<<residueNumber2<<endl;
            Transform myX_GB1;
            Transform myX_GB2;
            if (myLeontisWesthofBondRow.isTwoTransformForce.compare("aromatic")==0) {
                //new method that uses pre-computed corrections for topology
            	calcAxes(state, myLeontisWesthofBondRow, residueNumber1, residueNumber2, chainId1, chainId2, xAxisVector1,yAxisVector1,zAxisVector1,xAxisVector2,yAxisVector2,zAxisVector2,glycosidicNitrogenAtom1LocationInGround,glycosidicNitrogenAtom2LocationInGround,ring1CenterLocationInGround,ring2CenterLocationInGround);
            	Rotation rotation1FromRingAtoms(Mat33(xAxisVector1,yAxisVector1,zAxisVector1));
                if (myParameterReader.verbose) cout<<__FILE__<<":"<<__LINE__<<" :xAxisVector2, yAxisVector2, zAxisVector2 = "<< xAxisVector2<<" , "<< yAxisVector2<<" , "<< zAxisVector2   <<endl;	
            	Rotation rotation2FromRingAtoms(Mat33(xAxisVector2,yAxisVector2,zAxisVector2));
                myX_GB1 =Transform(rotation1FromRingAtoms,ring1CenterLocationInGround);
                myX_GB2 =Transform(rotation2FromRingAtoms,ring2CenterLocationInGround);
            }
            else { 
                //new method that uses pre-computed corrections for topology
                
                myX_GB1 =  matter.getMobilizedBody(body1).getBodyTransform(state);
                myX_GB2 =  matter.getMobilizedBody(body2).getBodyTransform(state); 
        
                myX_GB1 = Transform(~(myParameterReader.baseOperationVector[r].rotationCorrection1*~myX_GB1.R()),myX_GB1.R()*myParameterReader.baseOperationVector[r].translationCorrection1+myX_GB1.T());
                myX_GB2 = Transform(~(myParameterReader.baseOperationVector[r].rotationCorrection2*~myX_GB2.R()),myX_GB2.R()*myParameterReader.baseOperationVector[r].translationCorrection2+myX_GB2.T());
                //old way:
                if (0) {
            	Rotation rotation1FromRingAtoms(Mat33(xAxisVector1,yAxisVector1,zAxisVector1));
            	Rotation rotation2FromRingAtoms(Mat33(xAxisVector2,yAxisVector2,zAxisVector2));
               
            	if (myParameterReader.verbose) cout<<"[TwoTransformForces.cpp] computed rotation for residue 1 based on two-vector constructor = "<<Rotation(SimTK::UnitVec3(xAxisVector1),CoordinateAxis::XCoordinateAxis(),zAxisVector1, CoordinateAxis::ZCoordinateAxis())<<endl;


                } 
            }
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

            float scale = -2* 1000;
            float xoffset = -1.2;
            float xmultiplier = .4;
            float yoffset = 3.5;
            float myFrcScalar;
            float A, B, C;
            if (myParameterReader.verbose) cout<<"[TwoTransformForces.cpp] myParameterReader.potentialType :"<<myParameterReader.potentialType<<endl;
            for (int i = 0; i<3; i++) torque[i] = rotationAngleAxis[i+1];

            if (myLeontisWesthofBondRow.isTwoTransformForce.compare("aromatic")==0) 
                {torque = torque - dot(((rot1_G).z())/(rot1_G.z()).norm(),torque)*(((rot1_G).z())/(rot1_G.z()).norm());   } 
            float theta = rotationAngleAxis[0];
            torque *= -theta * torqueConstant ;

            if ((myParameterReader.potentialType).compare("HarmonicInverseLorentzian") == 0)
                {if (stretch < cutoffRadius) 
                     {
                     C =  forceConstant*cutoffRadius;
                     A = -forceConstant/2/cutoffRadius/cutoffRadius;//    -C/2/cutoffRadius/cutoffRadius/cutoffRadius;
                     B = +C/cutoffRadius-A*cutoffRadius*cutoffRadius;
                     myFrcScalar = 2*A*stretch; // overall positive quantity
                     myFrcScalar += theta*theta*torqueConstant*stretch/pow((1+pow((stretch/cutoffRadius),2)),2);
                     } 
                 else {
                     C = forceConstant*cutoffRadius;  
                     myFrcScalar = -C/stretch/stretch; //overall positive quantity.  Must be positive so it points towards body 2.
                     myFrcScalar += theta*theta*torqueConstant*stretch/pow((1+pow((stretch/cutoffRadius),2)),2);
                 }
                 torque *= 1/(1+pow((stretch/cutoffRadius),2));
                }
            else if ((myParameterReader.potentialType).compare("Harmonic" ) == 0) 
                {if (stretch > cutoffRadius) stretch =  0;myFrcScalar = forceConstant*stretch;  }
            else if ((myParameterReader.potentialType).compare("HarmonicLinear" ) == 0) 
                {if (stretch > cutoffRadius) stretch =  cutoffRadius;myFrcScalar = forceConstant*stretch;  }
            else if ((myParameterReader.potentialType).compare("Inverse") == 0) 
                {if (stretch < cutoffRadius) {myFrcScalar = 0;}
                 else {float C = forceConstant*cutoffRadius;  myFrcScalar = -C/stretch/stretch; }  // C should be negatve 
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
                      //float B =      +C/cutoffRadius-A*cutoffRadius*cutoffRadius;
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
	    if (myParameterReader.verbose) {cout<<"[TwoTransformForces.cpp] torque on body 1 ="<<torque + s1_G % f1_G<<endl;}
	    if (myParameterReader.verbose) {cout<<"[TwoTransformForces.cpp] torque on body 2 ="<<-torque - s2_G % f1_G<<endl;}
             
            bodyForces[body1.getMobilizedBodyIndex()] +=  SpatialVec(torque + (-(matter.getMobilizedBody(body1).getBodyTransform(state)).T()+p1_G) % f1_G, f1_G);
            bodyForces[body2.getMobilizedBodyIndex()] -=  SpatialVec(torque + (-(matter.getMobilizedBody(body2).getBodyTransform(state)).T()+p2_G) % f1_G, f1_G);
            }
            }
        };
    Real AllTwoTransformLinearSprings::calcPotentialEnergy(const State& state) const { 
        double energy = 0.0; 
   
        MobilizedBody body1;
        MobilizedBody body2;
        Transform transform1;
        Transform transform2;
        float forceConstant;
        float torqueConstant;
        float dutyCycle; //must be between 0 and 1.  at 1, force is applied all the time.  at 0, basically never applied.    
        double scrubberPeriod;
        float cutoffRadius;
        ParameterReader tempParameterReader =  myParameterReader;          

        float interactingBaseRMSD=0; 
             
        int   interactingBaseRMSDCount = 0;              
        if (myParameterReader.verbose) cout<< __FILE__<<":"<<__LINE__  <<" energy before calcPotentialEnergy  :: "<< std::setiosflags(std::ios::fixed) << std::setprecision(1) << energy    <<","<<endl;
        myParameterReader.satisfiedBasePairs = 0;
        myParameterReader.unSatisfiedBasePairs = 0;
        for (int r=0;r<(myParameterReader.baseOperationVector).size();r++)
        if (( (myParameterReader.baseOperationVector[r]).BasePairIsTwoTransformForce.compare("baseInteraction") ==0) ||
            ( (myParameterReader.baseOperationVector[r]).BasePairIsTwoTransformForce.compare("aromatic"         ) ==0))   
        { 


 
            string chainId1=(myParameterReader.baseOperationVector[r]).FirstBPChain;            
            string chainId2=(myParameterReader.baseOperationVector[r]).SecondBPChain;
            string bondingEdge1=(myParameterReader.baseOperationVector[r]).FirstBPEdge;
            string bondingEdge2=(myParameterReader.baseOperationVector[r]).SecondBPEdge;
            string glycosidicBondOrientation=(myParameterReader.baseOperationVector[r]).OrientationBP ;
            string basePairIsTwoTransformForce=(myParameterReader.baseOperationVector[r]).BasePairIsTwoTransformForce ;
            /*int i = -11111;
            int j = -11111;

            for (int k = 0; k < myParameterReader.numChains; k++) {
                if (strcmp((myParameterReader.chainId[k]).c_str(),chainId1.c_str()) == 0) i = k;
                if (strcmp((myParameterReader.chainId[k]).c_str(),chainId2.c_str()) == 0) j = k;
            }
            assert             (i>=0);
            assert(j>=0);
*/

            int residueNumber1=(myParameterReader.baseOperationVector[r]).FirstBPResidue-myBiopolymerClassContainer.updBiopolymerClass(chainId1).getFirstResidueNumber();
            int residueNumber2=(myParameterReader.baseOperationVector[r]).SecondBPResidue-myBiopolymerClassContainer.updBiopolymerClass(chainId1).getFirstResidueNumber() ;
            string myPdbResidueName1 = (myBiopolymerClassContainer.getPdbResidueName(chainId1,residueNumber1));
            string myPdbResidueName2 = (myBiopolymerClassContainer.getPdbResidueName(chainId2,residueNumber2));
            LeontisWesthofBondRow myLeontisWesthofBondRow = myLeontisWesthofClass.getLeontisWesthofBondRow((myParameterReader.baseOperationVector[r]).FirstBPResidue ,(myParameterReader.baseOperationVector[r]).SecondBPResidue,myPdbResidueName1,bondingEdge1,myPdbResidueName2,bondingEdge2,glycosidicBondOrientation,  basePairIsTwoTransformForce);
            stringstream ss3;
            ss3<<residueNumber1<<"/"<<myLeontisWesthofBondRow.residue1Atom[0];
            stringstream ss4;
            ss4<<residueNumber2<<"/"<<(myLeontisWesthofBondRow.residue2Atom[0]);
            body1 = myBiopolymerClassContainer.updAtomMobilizedBody(matter,chainId1,residueNumber1,myLeontisWesthofBondRow.residue1Atom[0]);
            body2 = myBiopolymerClassContainer.updAtomMobilizedBody(matter,chainId2,residueNumber2,myLeontisWesthofBondRow.residue2Atom[0]);
            Vec3 station1 = myLeontisWesthofBondRow.attachmentPoint;
            Vec3 station2 = myBiopolymerClassContainer.getAtomLocationInMobilizedBodyFrame(chainId2,residueNumber2,myLeontisWesthofBondRow.residue2Atom[0]);
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
                float scale = -2;
                float xoffset = -1.2;
                float xmultiplier = .4;
                float yoffset = 3.5;


//scf
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
            if (myLeontisWesthofBondRow.isTwoTransformForce.compare("aromatic")==0) {
                //new method that uses pre-computed corrections for topology
                calcAxes(state, myLeontisWesthofBondRow, residueNumber1, residueNumber2, chainId1, chainId2, xAxisVector1,yAxisVector1,zAxisVector1,xAxisVector2,yAxisVector2,zAxisVector2
,glycosidicNitrogenAtom1LocationInGround,glycosidicNitrogenAtom2LocationInGround,ring1CenterLocationInGround,ring2CenterLocationInGround);
                Rotation rotation1FromRingAtoms(Mat33(xAxisVector1,yAxisVector1,zAxisVector1));
                Rotation rotation2FromRingAtoms(Mat33(xAxisVector2,yAxisVector2,zAxisVector2));
                myX_GB1 =Transform(rotation1FromRingAtoms,ring1CenterLocationInGround);
                myX_GB2 =Transform(rotation2FromRingAtoms,ring2CenterLocationInGround);
            }
            else {
                //new method that uses pre-computed corrections for topology
                myX_GB1 =  matter.getMobilizedBody(body1).getBodyTransform(state);
                myX_GB2 =  matter.getMobilizedBody(body2).getBodyTransform(state);
                myX_GB1 = Transform(myParameterReader.baseOperationVector[r].rotationCorrection1*myX_GB1.R(),myX_GB1.R()*myParameterReader.baseOperationVector[r].translationCorrection1+myX_GB1.T());
                myX_GB2 = Transform(myParameterReader.baseOperationVector[r].rotationCorrection2*myX_GB2.R(),myX_GB2.R()*myParameterReader.baseOperationVector[r].translationCorrection2+myX_GB2.T());
            }
            const Transform X_GB1 = myX_GB1;
            const Transform X_GB2 = myX_GB2;

            const Vec3 s1_G = X_GB1.R() * station1;
            const Vec3 s2_G = X_GB2.R() * station2;
            const Vec3 p1_G = X_GB1.T() + s1_G; // station measured from ground origin
            const Vec3 p2_G = X_GB2.T() + s2_G;
            const Vec3 r_G = p2_G - p1_G; // vector from point1 to point2
            if (myParameterReader.verbose) cout<< __FILE__<<":"<<__LINE__  <<" r_G.norm(): "<<r_G.norm()<<","<<endl;

            if (myParameterReader.verbose) cout<< __FILE__<<":"<<__LINE__  <<" energy   :: "<<energy    <<","<<endl;

            float A, B, C;
            if (myParameterReader.verbose) cout<<"[TwoTransformForces.cpp] myParameterReader.potentialType :"<<myParameterReader.potentialType<<endl;
            if ((myParameterReader.potentialType).compare("HarmonicInverseLorentzian" ) == 0) {
                {
                const Rotation rot1_G = X_GB1.R() * rotation1;
                const Rotation rot2_G = X_GB2.R() * rotation2;
                const Vec4 rotationAngleAxis = (rot1_G*(~rot2_G)).convertRotationToAngleAxis();
                float theta = rotationAngleAxis[0];

                if (r_G.norm() < cutoffRadius)
                     {
                      C = forceConstant*cutoffRadius;
                      A   =  -forceConstant/2/cutoffRadius/cutoffRadius;//    -C/2/cutoffRadius/cutoffRadius/cutoffRadius;
                      float B =  +C/cutoffRadius-A*cutoffRadius*cutoffRadius;
                      if (myParameterReader.verbose) cout<< __FILE__<<":"<<__LINE__  <<" A,C :"<<A<<","<<C<<endl;
                      energy      += (
                          B+ A*r_G.norm()*r_G.norm()
                          +theta*theta*torqueConstant/2/(1+pow((r_G.norm()/cutoffRadius),2))           
                          )*myParameterReader.twoTransformForceMultiplier;
                     }
                 else {
                     C = forceConstant*cutoffRadius; 
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
            else if ((myParameterReader.potentialType).compare("Inverse") == 0)
                {if (r_G.norm() <= cutoffRadius) energy += C/cutoffRadius*myParameterReader.twoTransformForceMultiplier;
                 else {float C = forceConstant*cutoffRadius;  energy += C/r_G.norm()*myParameterReader.twoTransformForceMultiplier; }  // C should be negatve
                }
            else if ((myParameterReader.potentialType).compare("HarmonicInverse") == 0)
                {if (r_G.norm() <= cutoffRadius)
                     {
                      C = forceConstant*cutoffRadius;
                      A   =  -forceConstant/2/cutoffRadius/cutoffRadius;//    -C/2/cutoffRadius/cutoffRadius/cutoffRadius;
                      float B =  +C/cutoffRadius-A*cutoffRadius*cutoffRadius;
                      if (myParameterReader.verbose) cout<< __FILE__<<":"<<__LINE__  <<" A,C :"<<A<<","<<C<<endl;
                      energy      += (B+ A*r_G.norm()*r_G.norm())*myParameterReader.twoTransformForceMultiplier;
                     if (myParameterReader.verbose) cout<< __FILE__<<":"<<__LINE__  <<" small separation, just added this much to energy:"<< B+ A*r_G.norm()*r_G.norm()  <<endl;
                      if (myParameterReader.verbose) cout<< __FILE__<<":"<<__LINE__  <<" r_G.norm(): "<<r_G.norm()<<","<<endl;
                     }
                 else {
                     C = forceConstant*cutoffRadius; energy += C/r_G.norm()*myParameterReader.twoTransformForceMultiplier;
                     if (myParameterReader.verbose) cout<< __FILE__<<":"<<__LINE__  <<" large separation,  energy is now:"<< energy      <<endl;
                     if (myParameterReader.verbose) cout<< __FILE__<<":"<<__LINE__  <<" C :"<<C<<endl;
                     if (myParameterReader.verbose) cout<< __FILE__<<":"<<__LINE__  <<" r_G.norm(): "<<r_G.norm()<<","<<endl;
                 }
                 if (myParameterReader.verbose) cout<< __FILE__<<":"<<__LINE__  <<" r_G.norm(): "<<r_G.norm()<<","<<endl;
                }
       

            else if ((myParameterReader.potentialType).compare("HarmonicInverseCube") == 0)
                {if (r_G.norm() < cutoffRadius) 
                     {
                      A   =  -3*forceConstant/2/cutoffRadius/cutoffRadius;
                      float B = 5/2*forceConstant;
                      energy      += (B+ A*r_G.norm()*r_G.norm())*myParameterReader.twoTransformForceMultiplier;
                      if (myParameterReader.verbose) cout<< __FILE__<<":"<<__LINE__  <<" A,C :"<<A<<","<<C<<endl;
                     } 
                 else {
                     C = forceConstant*cutoffRadius*cutoffRadius*cutoffRadius; // myFrcScalar = -3*C/r_G.norm()/r_G.norm()/r_G.norm()/r_G.norm();
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
                
            } // add nothing if r.norm()>cutoffRadius
            else {
            }
            
            if (((((myParameterReader.baseOperationVector)[r]).BasePairIsTwoTransformForce).compare("baseInteraction") ==0) ) { 
               
                outputStream <<"REMARK Base pair : "<< myPdbResidueName1<<" " <<(myParameterReader.baseOperationVector[r]).FirstBPResidue<<" "<< bondingEdge1 <<", "  
                                    << myPdbResidueName2<<" " <<(myParameterReader.baseOperationVector[r]).SecondBPResidue<<" "<< bondingEdge2 
                                    <<" "<<glycosidicBondOrientation<<".  ";//<<endl;
                outputStream <<"  Base pair "<<r<<" attachment-body frame distance = "<<r_G.norm() <<endl;
                interactingBaseRMSD += r_G.norm() *  r_G.norm() ;
                interactingBaseRMSDCount++;
                if (r_G.norm() <= cutoffRadius) {
                    (tempParameterReader.baseOperationVector[r]).basePairSatisfied = "True";
                    outputStream <<"REMARK "<<__FILE__<<" :  Base pair "<<r<<" IS satisfied"<<endl;
		    myParameterReader.satisfiedBasePairs ++;	


                } else {
                    ((myParameterReader.baseOperationVector)[r]).basePairSatisfied = "False";
                    outputStream  <<"REMARK "<<__FILE__ <<" :  Base pair "<<r<<" is NOT satisfied"<<endl;
		    myParameterReader.unSatisfiedBasePairs ++;
                }
                }
            }    
        }
        //State  tempState = state;
        //setParameterReader(tempState,tempParameterReader);   
            
        if (myParameterReader.verbose) cout<< __FILE__<<":"<<__LINE__  <<" energy after calcPotentialEnergy  :: "<<energy    <<","<<endl;
        //cout << "[TwoTransformForces.cpp] Angular, Linear Momentum = "<<system.calcSystemRigidBodyMomentum(state) << endl;    
        cout << "[TwoTransformForces.cpp] Satisfied baseInteraction's:   : "<<myParameterReader.satisfiedBasePairs<< " out of : "<<myParameterReader.satisfiedBasePairs+ myParameterReader.unSatisfiedBasePairs<<endl;    
        outputStream << "REMARK [TwoTransformForces.cpp] Satisfied base interactions   : "<<myParameterReader.satisfiedBasePairs<< " out of : "<<myParameterReader.satisfiedBasePairs+ myParameterReader.unSatisfiedBasePairs<<endl;    
        outputStream <<"REMARK RMSD for interacting frames = "<<  pow((interactingBaseRMSD / interactingBaseRMSDCount), float(0.5))<<endl;
        return energy; 
    };
    bool AllTwoTransformLinearSprings::dependsOnlyOnPositions() const  { 
        return true; 
    };    
