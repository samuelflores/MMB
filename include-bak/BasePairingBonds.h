/* -------------------------------------------------------------------------- *
 *                           MMB (MacroMoleculeBuilder)                       *
 * -------------------------------------------------------------------------- *
 *                                                                            *
 * Copyright (c) 2011-12 by the Author.                                       *
 * Author: Samuel Flores                                                      *
 *                                                                            *
 * See RNABuilder.cpp for the copyright and usage agreement.                  *
 * -------------------------------------------------------------------------- */
//using std::cout;
//using std::endl;

//using namespace SimTK;
//using namespace std;



class ConcentricityPoint  { 
public:
    Vec3 location2;
    Vec3 normal2InBody2Frame;
    ConcentricityPoint (
        Biopolymer & myChain,        int residueNumber2,
        LeontisWesthofBondRow & myLeontisWesthofBondRow,        SimbodyMatterSubsystem & matter,
        State & state        )
        {
                                stringstream ss2;
                                stringstream ss2a;                                stringstream ss2b;
                                ss2 <<residueNumber2<<"/"<<myLeontisWesthofBondRow.residue2Atom[0];
                                ss2a<<residueNumber2<<"/"<<myLeontisWesthofBondRow.residue2Atom[1];                                ss2b<<residueNumber2<<"/"<<myLeontisWesthofBondRow.residue2Atom[2];
                                cout<<"[Repel.h:ConcentricityPoint] ss2.str() ="<<ss2.str()<<endl;
                                Compound::AtomIndex myAtomIndex2 = myChain.getAtomIndex( ss2.str());
                                Compound::AtomIndex myAtomIndex2a= myChain.getAtomIndex(ss2a.str());
                                Compound::AtomIndex myAtomIndex2b= myChain.getAtomIndex(ss2b.str());

                                MobilizedBodyIndex myMobilizedBodyIndex2  = (myChain).getAtomMobilizedBodyIndex(myAtomIndex2 );
                                MobilizedBodyIndex myMobilizedBodyIndex2a = (myChain).getAtomMobilizedBodyIndex(myAtomIndex2a);
                                MobilizedBodyIndex myMobilizedBodyIndex2b = (myChain).getAtomMobilizedBodyIndex(myAtomIndex2b);

                                MobilizedBody body2 = matter.getMobilizedBody(myMobilizedBodyIndex2 );
                                MobilizedBody body2a= matter.getMobilizedBody(myMobilizedBodyIndex2a);
                                MobilizedBody body2b= matter.getMobilizedBody(myMobilizedBodyIndex2b);
     
                                // get normal in ground frame



                                Vec3 normal2InGroundFrame = 
                                        (    
                                        body2a.getBodyTransform(state) * myChain.getAtomLocationInMobilizedBodyFrame(myAtomIndex2a)
                                        - body2.getBodyTransform(state) * myChain.getAtomLocationInMobilizedBodyFrame(myAtomIndex2)
                                        ) % (
                                        body2b.getBodyTransform(state) * myChain.getAtomLocationInMobilizedBodyFrame(myAtomIndex2b)
                                        - body2.getBodyTransform(state) * myChain.getAtomLocationInMobilizedBodyFrame(myAtomIndex2)
                                        );   
                                cout <<"[Repel.h:ConcentricityPoint] myChain.getAtomLocationInMobilizedBodyFrame(myAtomIndex2) ="<<myChain.getAtomLocationInMobilizedBodyFrame(myAtomIndex2)<<endl;
                                cout <<"[Repel.h:ConcentricityPoint] AtomLocation(2a) ="<<  (body2.getBodyTransform(state)).invert() * body2a.getBodyTransform(state) * myChain.getAtomLocationInMobilizedBodyFrame(myAtomIndex2a)<<endl;
                                cout <<"[Repel.h:ConcentricityPoint] AtomLocation(2b) ="<<  (body2.getBodyTransform(state)).invert() * body2b.getBodyTransform(state) *  myChain.getAtomLocationInMobilizedBodyFrame(myAtomIndex2b)<<endl;
     
                                normal2InBody2Frame = (
                                    ((body2.getBodyTransform(state)).invert() * body2a.getBodyTransform(state) * myChain.getAtomLocationInMobilizedBodyFrame(myAtomIndex2a)  
                                    - myChain.getAtomLocationInMobilizedBodyFrame(myAtomIndex2)) %
                                    ((body2.getBodyTransform(state)).invert() * body2b.getBodyTransform(state) * myChain.getAtomLocationInMobilizedBodyFrame(myAtomIndex2b)                                     - myChain.getAtomLocationInMobilizedBodyFrame(myAtomIndex2)));
                                         
                                //Vec3 normal2InBody2Frame = (body2.getBodyTransform(state)).invert() * normal2InGroundFrame;
                                normal2InBody2Frame = normal2InBody2Frame/normal2InBody2Frame.norm();   
                                location2 = myChain.getAtomLocationInMobilizedBodyFrame(myAtomIndex2) - (0.4 * normal2InBody2Frame);
                                                 
        }
};



class  BasePairingBonds  : public Biopolymer { 

public:
    BasePairingBonds (
       
        //vector<Biopolymer> &  myChain,
        Biopolymer  myChain[maxChains],
	LeontisWesthofClass  & myLeontisWesthofClass ,    
	GeneralForceSubsystem & forces,
        SimbodyMatterSubsystem & matter,
        HuntCrossleyContact & myHuntCrossleyContact,
        GeneralContactSubsystem & contacts,
        ContactSetIndex contactSet,
        HuntCrossleyForce & hc,
	State & state,
        float dutyCycle,
        ParameterReader & myParameterReader
	) 
        {
	for (int r=0;r<myParameterReader.baseOperationVector.size();r++)  
        if  (((myParameterReader.firstBPEdge[r]).compare("Parallelness") != 0) &&               
	    ((myParameterReader.firstBPEdge[r]).compare("ChiBond") != 0) && 
            ((myParameterReader.firstBPEdge[r]).compare("Rigid") != 0) && 
            ((myParameterReader.firstBPEdge[r]).compare("PointToPlane") != 0))  
	{
	 
            String chainId1=myParameterReader.firstBPChain[r];
            String chainId2=myParameterReader.secondBPChain[r];
	    String bondingEdge1=myParameterReader.firstBPEdge[r];
	    String bondingEdge2=myParameterReader.secondBPEdge[r];
	    String glycosidicBondOrientation=myParameterReader.orientationBP[r] ;   	    
	    
            int i = -11111;
            int j = -11111;
            
            for (int k = 0; k < myParameterReader.numChains; k++) {
 	        if (strcmp((myParameterReader.chainId[k]).c_str(),chainId1.c_str()) == 0) i = k;
 	        if (strcmp((myParameterReader.chainId[k]).c_str(),chainId2.c_str()) == 0) j = k;
	    } 
	    assert(i>=0);
	    assert(j>=0);

            int residueNumber1=myParameterReader.firstBPResidue[r]-myParameterReader.getFirstResidueNumbers((myParameterReader.chainId[i]).c_str());
	    int residueNumber2=myParameterReader.secondBPResidue[r]-myParameterReader.getFirstResidueNumbers((myParameterReader.chainId[j]).c_str());

	    String myPdbResidueName1 = (myChain[i].updResidue(Compound::Index(residueNumber1))).getPdbResidueName();
	    String myPdbResidueName2  = (myChain[j].updResidue(Compound::Index(residueNumber2))).getPdbResidueName();
            int hardSphereCondition = ((bondingEdge1.compare("HardSphere")==0)||(bondingEdge1.compare("HardSphereSmall")==0));
	    if (! hardSphereCondition)  myPdbResidueName2 = (myChain[j].updResidue(Compound::Index(residueNumber2))).getPdbResidueName();
            Force::TwoPointLinearSpring  * myTwoPointLinearSpringArray[4];               
            Force::TwoPointLinearDamper  * myTwoPointLinearDamperArray[4];    
           
            float huntCrossleyDissipation  =.0;
	    if (myParameterReader.verbose) cout<<"[BasePairingBonds.h : 148] hardSphereCondition =" <<hardSphereCondition<<".  bondingEdge1 = "<<bondingEdge1<<endl;
            LeontisWesthofBondRow myLeontisWesthofBondRow = myLeontisWesthofClass.getLeontisWesthofBondRow(myPdbResidueName1,bondingEdge1,myPdbResidueName2,bondingEdge2,glycosidicBondOrientation);
	    if (myLeontisWesthofBondRow.residue1Atom[0].compare("") == 0) if (myParameterReader.verbose) cout << "[BasePairingBonds] No first base pairing atom.  Base pairing information is not available for this combination of bond orientation, base pairing type, and residue types.  If you wish to specify such a base pairing please edit the appropriate Leontis-Westhof bond file."<<endl; 
	    int maxAtoms;

	    if (((bondingEdge1.compare("WatsonCrick") == 0) && (bondingEdge2.compare("WatsonCrick") == 0) )|| ((bondingEdge1.compare("WatsonCrick") == 0) && (bondingEdge2.compare("Hoogsteen") == 0) )  || (bondingEdge1.compare("HelicalStacking") == 0) || (bondingEdge1.compare("Concentricity") == 0)|| (bondingEdge1.compare("Stacking") == 0) || (bondingEdge1.compare("ChiBondAnti") == 0))   maxAtoms = 1; else maxAtoms=4;
            for (int q =0; q<maxAtoms; q++) { //loop over four possible interacting pairs of atoms.
		if ((((myLeontisWesthofBondRow.residue1Atom[q]).compare("") != 0)&&((myLeontisWesthofBondRow.residue2Atom[q]).compare("") != 0 )) ||(((myLeontisWesthofBondRow.residue1Atom[q]).compare("") != 0) && hardSphereCondition )) { //if the atom name field is not blank, do this
                        Vec3 location1;
                        Rotation rotation1;
                        Vec3 location2;
                        Rotation rotation2;
			stringstream ss3;
			ss3<<residueNumber1<<"/"<<myLeontisWesthofBondRow.residue1Atom[q];
			stringstream ss3a;
			ss3a<<residueNumber1<<"/C1*";//<<myLeontisWesthofBondRow.residue1Atom[q];
			stringstream ss4;
			if (! hardSphereCondition) ss4<<residueNumber2<<"/"<<(myLeontisWesthofBondRow.residue2Atom[q]);
			stringstream ss4a;
			if (! hardSphereCondition) ss4a<<residueNumber2<<"/C1*";
			if ((bondingEdge1.compare("Concentricity") == 0) && (bondingEdge2.compare("Concentricity") == 0) ) {
				//if (myParameterReader.verbose) cout <<"[Repel.h:BasePairingBonds] before ConcentricityPoint, location2 ="<<location2<<endl; 
			        location1 = myChain[i].getAtomLocationInMobilizedBodyFrame(myChain[i].getAtomIndex(ss3.str()));
    				ConcentricityPoint myConcentricityPoint(
					myChain[j],
					residueNumber2,
					myLeontisWesthofBondRow,
        				matter,
        				state
					);
				location2 = myConcentricityPoint.location2; 
				if (myParameterReader.verbose) cout <<"[Repel.h:BasePairingBonds] after ConcentricityPoint, location2 ="<<location2<<endl; 
				if (myParameterReader.verbose) cout <<"[Repel.h:BasePairingBonds] after ConcentricityPoint, q=========="<<q        <<endl; 
				if (myParameterReader.verbose) cout <<"[Repel.h:BasePairingBonds] after ConcentricityPoint, ss4.str() ="<<ss4.str()<<endl; 
			} else if (((bondingEdge1.compare("WatsonCrick") == 0) && (bondingEdge2.compare("WatsonCrick") == 0)) || ((bondingEdge1.compare("WatsonCrick") == 0) && (bondingEdge2.compare("Hoogsteen") == 0) )  || (bondingEdge1.compare("HelicalStacking") == 0) || (bondingEdge1.compare("Stacking") == 0) || (bondingEdge1.compare("ChiBondAnti") == 0)) {
                            
			    location1 = myLeontisWesthofBondRow.attachmentPoint; 
			    rotation1 = Rotation(myLeontisWesthofBondRow.rotationAngle,myLeontisWesthofBondRow.rotationAxis);
				
			    location2 = myChain[j].getAtomLocationInMobilizedBodyFrame(myChain[j].getAtomIndex(ss4.str()));
			    rotation2.setRotationToIdentityMatrix();
  		
			}
                        else if ((bondingEdge1.compare("Rigid") == 0  ) && (bondingEdge2.compare("Rigid") == 0  )) {
			   assert (i == j);
                           RNA& rna = static_cast<RNA&>(myChain[i]);
                           if (myParameterReader.verbose) cout<<"[BasePairingBonds.h] (residueNumber1), (residueNumber2):"<<(residueNumber1)<<","<< (residueNumber2)<<endl; 
			   rna.setRNABondMobility(BondMobility::Rigid ,(residueNumber1), (residueNumber2));
			   //rna.setRNABondMobility(myParameterReader.helixBondMobility,(residueNumber1), (residueNumber2));
			   myChain[i] = rna;
			}
			else if (hardSphereCondition) {
			    location1 = myChain[i].getAtomLocationInMobilizedBodyFrame(myChain[i].getAtomIndex(ss3.str()));

			}
			else {
			    location1 = myChain[i].getAtomLocationInMobilizedBodyFrame(myChain[i].getAtomIndex(ss3.str()));
			    location2 = myChain[j].getAtomLocationInMobilizedBodyFrame(myChain[j].getAtomIndex(ss4.str()));
			}
			//if (myParameterReader.verbose) cout<<"about to Hydrogen Bond atoms :"<<ss3.str()<<","<<ss4.str()<<endl;
	               	if (myParameterReader.verbose) cout<<"[Repel.h] myLeontisWesthofBondRow.springConstant["<<q<<"] ="<<myLeontisWesthofBondRow.springConstant[q]<<endl;
			if (myParameterReader.verbose) cout<<"[BasePairingBonds.h] hardSphereCondition =" <<hardSphereCondition<<".  bondingEdge1 = "<<bondingEdge1<<endl;
			if (hardSphereCondition){
			    if (myParameterReader.verbose) cout<<"[BasePairingBonds.h] adding sphere at "<<ss3.str()<< " myLeontisWesthofBondRow.bondLength[q] ="<<myLeontisWesthofBondRow.bondLength[q]<<" myLeontisWesthofBondRow.springConstant[q] =" <<myLeontisWesthofBondRow.springConstant[q]<<   endl;	

                            contacts.addBody(contactSet,
                                        (matter.updMobilizedBody(       myChain[i].getAtomMobilizedBodyIndex(Compound::AtomIndex(myChain[i].getAtomIndex(ss3.str()))))), 
                                        ContactGeometry::Sphere(myLeontisWesthofBondRow.bondLength[q]), 
                                        (    myChain[i].getAtomLocationInMobilizedBodyFrame(       myChain[i].getAtomIndex(ss3.str())))
                                        );   
                            hc.setBodyParameters(contacts.getNumBodies(contactSet)-1,myLeontisWesthofBondRow.springConstant[q] ,huntCrossleyDissipation, 0., 0., 0.); 
                            //hc.setBodyParameters(contacts.getNumBodies(contactSet)-1,huntCrossleyStiffness ,huntCrossleyDissipation, 0., 0., 0.); 
		} 
			else {
			//for (int j = 0; j < myParameterReader.numChains; j++) assert(myParameterReader.firstResidueNumbers[(myParameterReader.chainId[j]).c_str()] > 0);
			if (myParameterReader.verbose) cout<<"[BasePairingBonds.h] adding TwoTransformLinearSpring at "<<ss3.str()<<endl;	
			if (myParameterReader.verbose) cout<<"[BasePairingBonds.h] location1,location2 ="<< location1<<","<<location2 <<endl;	
			
		        TwoTransformLinearSpring * myTwoTransformLinearSpringPointer = new TwoTransformLinearSpring 
				(
				matter,
                        	matter.updMobilizedBody(myChain[i].getAtomMobilizedBodyIndex(Compound::AtomIndex(myChain[i].getAtomIndex(ss3.str())))),  
				Transform(rotation1,location1),
                        	matter.updMobilizedBody(myChain[j].getAtomMobilizedBodyIndex(Compound::AtomIndex(myChain[j].getAtomIndex(ss4.str())))),    
				Transform(rotation2,location2),
                        	myLeontisWesthofBondRow.springConstant[q],
                        	myLeontisWesthofBondRow.torqueConstant,
				dutyCycle,
				myParameterReader.scrubberPeriod,
				myParameterReader.cutoffRadius
				);
                        Force::Custom(forces, 
			    myTwoTransformLinearSpringPointer
			);
			}	
                }
            }
        }
	}

//class  BasePairingBonds  : public Biopolymer { public:
    BasePairingBonds (
        //Biopolymer  & myChain,
        //Biopolymer  myChain[maxChains],
        vector<Biopolymer> &  myChain,
	int residueNumber1,
	int residueNumber2,
	String chainId1   ,
	String chainId2   ,
        String bondingEdge1, 
	String bondingEdge2,
	String glycosidicBondOrientation,
	LeontisWesthofClass  & myLeontisWesthofClass ,    
	GeneralForceSubsystem & forces,
        SimbodyMatterSubsystem & matter,
        //GeneralContactSubsystem & myHuntCrossleyContact,
        HuntCrossleyContact & myHuntCrossleyContact,
        GeneralContactSubsystem & contacts,
        ContactSetIndex contactSet,
        HuntCrossleyForce & hc,
	State & state,
        //MagnesiumIon & myMagnesiumIon,
        float dutyCycle,
        ParameterReader & myParameterReader
	)   
	{    
            int i = -11111;
            int j = -11111;
            
            for (int k = 0; k < myParameterReader.numChains; k++) {
 	        if (strcmp((myParameterReader.chainId[k]).c_str(),chainId1.c_str()) == 0) i = k;
 	        if (strcmp((myParameterReader.chainId[k]).c_str(),chainId2.c_str()) == 0) j = k;
	    } 
	    assert(i>=0);
	    assert(j>=0);
	    String myPdbResidueName1 = (myChain[i].updResidue(Compound::Index(residueNumber1))).getPdbResidueName();
	    String myPdbResidueName2  = (myChain[j].updResidue(Compound::Index(residueNumber2))).getPdbResidueName();
            int hardSphereCondition = ((bondingEdge1.compare("HardSphere")==0)||(bondingEdge1.compare("HardSphereSmall")==0));
	    if (! hardSphereCondition)  myPdbResidueName2 = (myChain[j].updResidue(Compound::Index(residueNumber2))).getPdbResidueName();
            Force::TwoPointLinearSpring  * myTwoPointLinearSpringArray[4];               
            Force::TwoPointLinearDamper  * myTwoPointLinearDamperArray[4];    
           
            float huntCrossleyDissipation  =.0;
	    if (myParameterReader.verbose) cout<<"[BasePairingBonds.h : 291] hardSphereCondition =" <<hardSphereCondition<<".  bondingEdge1 = "<<bondingEdge1<<endl;
            LeontisWesthofBondRow myLeontisWesthofBondRow = myLeontisWesthofClass.getLeontisWesthofBondRow(myPdbResidueName1,bondingEdge1,myPdbResidueName2,bondingEdge2,glycosidicBondOrientation);
	    if (myLeontisWesthofBondRow.residue1Atom[0].compare("") == 0) if (myParameterReader.verbose) cout << "[BasePairingBonds] No first base pairing atom.  Base pairing information is not available for this combination of bond orientation, base pairing type, and residue types.  If you wish to specify such a base pairing please edit the appropriate Leontis-Westhof bond file."<<endl; 
	    int maxAtoms;

	    if (((bondingEdge1.compare("WatsonCrick") == 0) && (bondingEdge2.compare("WatsonCrick") == 0) )|| ((bondingEdge1.compare("WatsonCrick") == 0) && (bondingEdge2.compare("Hoogsteen") == 0) )  || (bondingEdge1.compare("HelicalStacking") == 0) || (bondingEdge1.compare("Concentricity") == 0)|| (bondingEdge1.compare("Stacking") == 0) || (bondingEdge1.compare("ChiBondAnti") == 0))   maxAtoms = 1; else maxAtoms=4;
            for (int q =0; q<maxAtoms; q++) { //loop over four possible interacting pairs of atoms.
		if ((((myLeontisWesthofBondRow.residue1Atom[q]).compare("") != 0)&&((myLeontisWesthofBondRow.residue2Atom[q]).compare("") != 0 )) ||(((myLeontisWesthofBondRow.residue1Atom[q]).compare("") != 0) && hardSphereCondition )) { //if the atom name field is not blank, do this
                        Vec3 location1;
                        Rotation rotation1;
                        Vec3 location2;
                        Rotation rotation2;
			stringstream ss3;
			ss3<<residueNumber1<<"/"<<myLeontisWesthofBondRow.residue1Atom[q];
			stringstream ss3a;
			ss3a<<residueNumber1<<"/C1*";//<<myLeontisWesthofBondRow.residue1Atom[q];
			stringstream ss4;
			if (! hardSphereCondition) ss4<<residueNumber2<<"/"<<(myLeontisWesthofBondRow.residue2Atom[q]);
			stringstream ss4a;
			if (! hardSphereCondition) ss4a<<residueNumber2<<"/C1*";
			if ((bondingEdge1.compare("Concentricity") == 0) && (bondingEdge2.compare("Concentricity") == 0) ) {
				//if (myParameterReader.verbose) cout <<"[Repel.h:BasePairingBonds] before ConcentricityPoint, location2 ="<<location2<<endl; 
			        location1 = myChain[i].getAtomLocationInMobilizedBodyFrame(myChain[i].getAtomIndex(ss3.str()));
    				ConcentricityPoint myConcentricityPoint(
					myChain[j],
					residueNumber2,
					myLeontisWesthofBondRow,
        				matter,
        				state
					);
				location2 = myConcentricityPoint.location2; 
				if (myParameterReader.verbose) cout <<"[Repel.h:BasePairingBonds] after ConcentricityPoint, location2 ="<<location2<<endl; 
				if (myParameterReader.verbose) cout <<"[Repel.h:BasePairingBonds] after ConcentricityPoint, q=========="<<q        <<endl; 
				if (myParameterReader.verbose) cout <<"[Repel.h:BasePairingBonds] after ConcentricityPoint, ss4.str() ="<<ss4.str()<<endl; 
			} else if (((bondingEdge1.compare("WatsonCrick") == 0) && (bondingEdge2.compare("WatsonCrick") == 0)) || ((bondingEdge1.compare("WatsonCrick") == 0) && (bondingEdge2.compare("Hoogsteen") == 0) )  || (bondingEdge1.compare("HelicalStacking") == 0) || (bondingEdge1.compare("Stacking") == 0) || (bondingEdge1.compare("ChiBondAnti") == 0)) {
                            
			    location1 = myLeontisWesthofBondRow.attachmentPoint; 
			    rotation1 = Rotation(myLeontisWesthofBondRow.rotationAngle,myLeontisWesthofBondRow.rotationAxis);
				
			    location2 = myChain[j].getAtomLocationInMobilizedBodyFrame(myChain[j].getAtomIndex(ss4.str()));
			    rotation2.setRotationToIdentityMatrix();
  		
			}
                        else if ((bondingEdge1.compare("Rigid") == 0  ) && (bondingEdge2.compare("Rigid") == 0  )) {
			   assert (i == j);
                           RNA& rna = static_cast<RNA&>(myChain[i]);
                           if (myParameterReader.verbose) cout<<"[BasePairingBonds.h] (residueNumber1), (residueNumber2):"<<(residueNumber1)<<","<< (residueNumber2)<<endl; 
			   rna.setRNABondMobility(BondMobility::Rigid ,(residueNumber1), (residueNumber2));
			   //rna.setRNABondMobility(myParameterReader.helixBondMobility,(residueNumber1), (residueNumber2));
			   myChain[i] = rna;
			}
			else if (hardSphereCondition) {
			    location1 = myChain[i].getAtomLocationInMobilizedBodyFrame(myChain[i].getAtomIndex(ss3.str()));

			}
			else {
			    location1 = myChain[i].getAtomLocationInMobilizedBodyFrame(myChain[i].getAtomIndex(ss3.str()));
			    location2 = myChain[j].getAtomLocationInMobilizedBodyFrame(myChain[j].getAtomIndex(ss4.str()));
			}
			//if (myParameterReader.verbose) cout<<"about to Hydrogen Bond atoms :"<<ss3.str()<<","<<ss4.str()<<endl;
	               	if (myParameterReader.verbose) cout<<"[Repel.h] myLeontisWesthofBondRow.springConstant["<<q<<"] ="<<myLeontisWesthofBondRow.springConstant[q]<<endl;
			if (myParameterReader.verbose) cout<<"[BasePairingBonds.h] hardSphereCondition =" <<hardSphereCondition<<".  bondingEdge1 = "<<bondingEdge1<<endl;
			if (hardSphereCondition){
			    if (myParameterReader.verbose) cout<<"[BasePairingBonds.h] adding sphere at "<<ss3.str()<< " myLeontisWesthofBondRow.bondLength[q] ="<<myLeontisWesthofBondRow.bondLength[q]<<" myLeontisWesthofBondRow.springConstant[q] =" <<myLeontisWesthofBondRow.springConstant[q]<<   endl;	

                            contacts.addBody(contactSet,
                                        (matter.updMobilizedBody(       myChain[i].getAtomMobilizedBodyIndex(Compound::AtomIndex(myChain[i].getAtomIndex(ss3.str()))))), 
                                        ContactGeometry::Sphere(myLeontisWesthofBondRow.bondLength[q]), 
                                        (    myChain[i].getAtomLocationInMobilizedBodyFrame(       myChain[i].getAtomIndex(ss3.str())))
                                        );   
                            hc.setBodyParameters(contacts.getNumBodies(contactSet)-1,myLeontisWesthofBondRow.springConstant[q] ,huntCrossleyDissipation, 0., 0., 0.); 
                            //hc.setBodyParameters(contacts.getNumBodies(contactSet)-1,huntCrossleyStiffness ,huntCrossleyDissipation, 0., 0., 0.); 
		} 
			else {
			//for (int j = 0; j < myParameterReader.numChains; j++) assert(myParameterReader.firstResidueNumbers[(myParameterReader.chainId[j]).c_str()] > 0);
			if (myParameterReader.verbose) cout<<"[BasePairingBonds.h] adding TwoTransformLinearSpring at "<<ss3.str()<<endl;	
			if (myParameterReader.verbose) cout<<"[BasePairingBonds.h] location1,location2 ="<< location1<<","<<location2 <<endl;	
			
                        Force::Custom(forces, 
		        new TwoTransformLinearSpring 
				(
				matter,
                        	matter.updMobilizedBody(myChain[i].getAtomMobilizedBodyIndex(Compound::AtomIndex(myChain[i].getAtomIndex(ss3.str())))),  
				Transform(rotation1,location1),
                        	matter.updMobilizedBody(myChain[j].getAtomMobilizedBodyIndex(Compound::AtomIndex(myChain[j].getAtomIndex(ss4.str())))),    
				Transform(rotation2,location2),
                        	myLeontisWesthofBondRow.springConstant[q],
                        	myLeontisWesthofBondRow.torqueConstant,
				dutyCycle,
				myParameterReader.scrubberPeriod,
				myParameterReader.cutoffRadius
				)
			);
			}	
                }
            }
        }





//class  BasePairingBonds  : public Biopolymer { public:
    BasePairingBonds (
        //Biopolymer  & myChain,
        Biopolymer  myChain[maxChains],
        //vector<Biopolymer> &  myChain,
	int residueNumber1,
	int residueNumber2,
	String chainId1   ,
	String chainId2   ,
        String bondingEdge1, 
	String bondingEdge2,
	String glycosidicBondOrientation,
	LeontisWesthofClass  & myLeontisWesthofClass ,    
	GeneralForceSubsystem & forces,
        SimbodyMatterSubsystem & matter,
        //GeneralContactSubsystem & myHuntCrossleyContact,
        HuntCrossleyContact & myHuntCrossleyContact,
        GeneralContactSubsystem & contacts,
        ContactSetIndex contactSet,
        HuntCrossleyForce & hc,
	State & state,
        //MagnesiumIon & myMagnesiumIon,
        float dutyCycle,
        ParameterReader & myParameterReader
	)   
	{    
            int i = -11111;
            int j = -11111;
            
            for (int k = 0; k < myParameterReader.numChains; k++) {
 	        if (strcmp((myParameterReader.chainId[k]).c_str(),chainId1.c_str()) == 0) i = k;
 	        if (strcmp((myParameterReader.chainId[k]).c_str(),chainId2.c_str()) == 0) j = k;
	    } 
	    assert(i>=0);
	    assert(j>=0);
	    String myPdbResidueName1 = (myChain[i].updResidue(Compound::Index(residueNumber1))).getPdbResidueName();
	    String myPdbResidueName2  = (myChain[j].updResidue(Compound::Index(residueNumber2))).getPdbResidueName();
            int hardSphereCondition = ((bondingEdge1.compare("HardSphere")==0)||(bondingEdge1.compare("HardSphereSmall")==0));
	    if (! hardSphereCondition)  myPdbResidueName2 = (myChain[j].updResidue(Compound::Index(residueNumber2))).getPdbResidueName();
            Force::TwoPointLinearSpring  * myTwoPointLinearSpringArray[4];               
            Force::TwoPointLinearDamper  * myTwoPointLinearDamperArray[4];    
           
            float huntCrossleyDissipation  =.0;
	    if (myParameterReader.verbose) cout<<"[BasePairingBonds.h : 436] hardSphereCondition =" <<hardSphereCondition<<".  bondingEdge1 = "<<bondingEdge1<<endl;
            LeontisWesthofBondRow myLeontisWesthofBondRow;
            if (myParameterReader.verbose) cout<<"[BasePairingBonds.h] hardSphereCondition = "<<hardSphereCondition<<endl;
            if (hardSphereCondition) {
                if (myParameterReader.verbose) cout<<"[BasePairingBonds.h] check 0.1"<<endl;
                myLeontisWesthofBondRow = myLeontisWesthofClass.getLeontisWesthofBondRow(myPdbResidueName1,bondingEdge1,"","","");
                }
            else {
                if (myParameterReader.verbose) cout<<"[BasePairingBonds.h] check 0.2"<<endl;
                myLeontisWesthofBondRow = myLeontisWesthofClass.getLeontisWesthofBondRow(myPdbResidueName1,bondingEdge1,myPdbResidueName2,bondingEdge2,glycosidicBondOrientation);}
	    if (myParameterReader.verbose) cout<<"[BasePairingBonds.h ] check 1"<<endl;
	    if (myLeontisWesthofBondRow.residue1Atom[0].compare("") == 0) cout << "[BasePairingBonds] No first base pairing atom.  Base pairing information is not available for this combination of bond orientation, base pairing type, and residue types.  If you wish to specify such a base pairing please edit the appropriate Leontis-Westhof bond file."<<endl; 
	    int maxAtoms;
	    if (myParameterReader.verbose) cout<<"[BasePairingBonds.h ] check 2"<<endl;

	    if (((bondingEdge1.compare("WatsonCrick") == 0) && (bondingEdge2.compare("WatsonCrick") == 0) )|| ((bondingEdge1.compare("WatsonCrick") == 0) && (bondingEdge2.compare("Hoogsteen") == 0) )  || (bondingEdge1.compare("HelicalStacking") == 0) || (bondingEdge1.compare("Concentricity") == 0)|| (bondingEdge1.compare("Stacking") == 0) || (bondingEdge1.compare("ChiBondAnti") == 0))   maxAtoms = 1; else maxAtoms=4;
	    if (myParameterReader.verbose) cout<<"[BasePairingBonds.h ] check 2.2"<<endl;
            for (int q =0; q<maxAtoms; q++) { //loop over four possible interacting pairs of atoms.
	        if (myParameterReader.verbose) cout<<"[BasePairingBonds.h ] check 2.3"<<endl;
	        if (myParameterReader.verbose) cout<<"[BasePairingBonds.h ] "<< myLeontisWesthofBondRow.residue1Atom[q] <<endl;
	        if (myParameterReader.verbose) cout<<"[BasePairingBonds.h ] "<< myLeontisWesthofBondRow.residue2Atom[q] <<endl;
	        if (myParameterReader.verbose) cout<<"[BasePairingBonds.h ] "<< hardSphereCondition <<endl;
		if ((((myLeontisWesthofBondRow.residue1Atom[q]).compare("") != 0)&&((myLeontisWesthofBondRow.residue2Atom[q]).compare("") != 0 )) ||(((myLeontisWesthofBondRow.residue1Atom[q]).compare("") != 0) && hardSphereCondition )) { //if the atom name field is not blank, do this
	                if (myParameterReader.verbose) cout<<"[BasePairingBonds.h ] check 2.5"<<endl;
                        Vec3 location1;
                        Rotation rotation1;
                        Vec3 location2;
                        Rotation rotation2;
			stringstream ss3;
			ss3<<residueNumber1<<"/"<<myLeontisWesthofBondRow.residue1Atom[q];
			stringstream ss3a;
			ss3a<<residueNumber1<<"/C1*";//<<myLeontisWesthofBondRow.residue1Atom[q];
			stringstream ss4;
			if (! hardSphereCondition) ss4<<residueNumber2<<"/"<<(myLeontisWesthofBondRow.residue2Atom[q]);
			stringstream ss4a;
			if (! hardSphereCondition) ss4a<<residueNumber2<<"/C1*";
           	        if (myParameterReader.verbose) cout<<"[BasePairingBonds.h ] check 3"<<endl;
			if ((bondingEdge1.compare("Concentricity") == 0) && (bondingEdge2.compare("Concentricity") == 0) ) {
				//if (myParameterReader.verbose) cout <<"[Repel.h:BasePairingBonds] before ConcentricityPoint, location2 ="<<location2<<endl; 
			        location1 = myChain[i].getAtomLocationInMobilizedBodyFrame(myChain[i].getAtomIndex(ss3.str()));
    				ConcentricityPoint myConcentricityPoint(
					myChain[j],
					residueNumber2,
					myLeontisWesthofBondRow,
        				matter,
        				state
					);
				location2 = myConcentricityPoint.location2; 
				if (myParameterReader.verbose) cout <<"[Repel.h:BasePairingBonds] after ConcentricityPoint, location2 ="<<location2<<endl; 
				if (myParameterReader.verbose) cout <<"[Repel.h:BasePairingBonds] after ConcentricityPoint, q=========="<<q        <<endl; 
				if (myParameterReader.verbose) cout <<"[Repel.h:BasePairingBonds] after ConcentricityPoint, ss4.str() ="<<ss4.str()<<endl; 
			} else if (((bondingEdge1.compare("WatsonCrick") == 0) && (bondingEdge2.compare("WatsonCrick") == 0)) || ((bondingEdge1.compare("WatsonCrick") == 0) && (bondingEdge2.compare("Hoogsteen") == 0) )  || (bondingEdge1.compare("HelicalStacking") == 0) || (bondingEdge1.compare("Stacking") == 0) || (bondingEdge1.compare("ChiBondAnti") == 0)) {
                            
			    location1 = myLeontisWesthofBondRow.attachmentPoint; 
			    rotation1 = Rotation(myLeontisWesthofBondRow.rotationAngle,myLeontisWesthofBondRow.rotationAxis);
				
			    location2 = myChain[j].getAtomLocationInMobilizedBodyFrame(myChain[j].getAtomIndex(ss4.str()));
			    rotation2.setRotationToIdentityMatrix();
  		
			}
                        /*
                        else if ((bondingEdge1.compare("Rigid") == 0  ) && (bondingEdge2.compare("Rigid") == 0  )) {
			   assert (i == j);
                           RNA& rna = static_cast<RNA&>(myChain[i]);
                           if (myParameterReader.verbose) cout<<"[BasePairingBonds.h] (residueNumber1), (residueNumber2):"<<(residueNumber1)<<","<< (residueNumber2)<<endl; 
			   rna.setRNABondMobility(BondMobility::Rigid ,(residueNumber1), (residueNumber2));
			   //rna.setRNABondMobility(myParameterReader.helixBondMobility,(residueNumber1), (residueNumber2));
			   myChain[i] = rna;
			}
 			*/
			else if (hardSphereCondition) {
	                    if (myParameterReader.verbose) cout<<"[BasePairingBonds.h ] check 5"<<endl;
			    location1 = myChain[i].getAtomLocationInMobilizedBodyFrame(myChain[i].getAtomIndex(ss3.str()));
	                    if (myParameterReader.verbose) cout<<"[BasePairingBonds.h ] check 6"<<endl;

			}
			else {
			    location1 = myChain[i].getAtomLocationInMobilizedBodyFrame(myChain[i].getAtomIndex(ss3.str()));
			    location2 = myChain[j].getAtomLocationInMobilizedBodyFrame(myChain[j].getAtomIndex(ss4.str()));
			}
			//if (myParameterReader.verbose) cout<<"about to Hydrogen Bond atoms :"<<ss3.str()<<","<<ss4.str()<<endl;
	               	if (myParameterReader.verbose) cout<<"[Repel.h] myLeontisWesthofBondRow.springConstant["<<q<<"] ="<<myLeontisWesthofBondRow.springConstant[q]<<endl;
			if (myParameterReader.verbose) cout<<"[BasePairingBonds.h] hardSphereCondition =" <<hardSphereCondition<<".  bondingEdge1 = "<<bondingEdge1<<endl;
			if (hardSphereCondition){
			    if (myParameterReader.verbose) cout<<"[BasePairingBonds.h] adding sphere at "<<ss3.str()<< " myLeontisWesthofBondRow.bondLength[q] ="<<myLeontisWesthofBondRow.bondLength[q]<<" myLeontisWesthofBondRow.springConstant[q] =" <<myLeontisWesthofBondRow.springConstant[q]<<   endl;	

                            contacts.addBody(contactSet,
                                        (matter.updMobilizedBody(       myChain[i].getAtomMobilizedBodyIndex(Compound::AtomIndex(myChain[i].getAtomIndex(ss3.str()))))), 
                                        ContactGeometry::Sphere(myLeontisWesthofBondRow.bondLength[q]), 
                                        (    myChain[i].getAtomLocationInMobilizedBodyFrame(       myChain[i].getAtomIndex(ss3.str())))
                                        );   
                            hc.setBodyParameters(contacts.getNumBodies(contactSet)-1,myLeontisWesthofBondRow.springConstant[q] ,huntCrossleyDissipation, 0., 0., 0.); 
                            //hc.setBodyParameters(contacts.getNumBodies(contactSet)-1,huntCrossleyStiffness ,huntCrossleyDissipation, 0., 0., 0.); 

/*
                            myHuntCrossleyContact.addSphere (     
                                (matter.updMobilizedBody(myChain.getAtomMobilizedBodyIndex(Compound::AtomIndex(myChain.getAtomIndex(ss3.str()))))).getMobilizedBodyIndex(),
                                myChain.getAtomLocationInMobilizedBodyFrame(myChain.getAtomIndex(ss3.str())) ,
                                myLeontisWesthofBondRow.bondLength[q],
                                myLeontisWesthofBondRow.springConstant[q],
                                huntCrossleyDissipation
                            ) ;  			    
*/			} 
			else {
			//for (int j = 0; j < myParameterReader.numChains; j++) assert(myParameterReader.firstResidueNumbers[(myParameterReader.chainId[j]).c_str()] > 0);
			if (myParameterReader.verbose) cout<<"[BasePairingBonds.h] adding TwoTransformLinearSpring at "<<ss3.str()<<endl;	
			if (myParameterReader.verbose) cout<<"[BasePairingBonds.h] location1,location2 ="<< location1<<","<<location2 <<endl;	

                        Force::Custom(forces, 
		        new TwoTransformLinearSpring 
				(
				matter,
                        	matter.updMobilizedBody(myChain[i].getAtomMobilizedBodyIndex(Compound::AtomIndex(myChain[i].getAtomIndex(ss3.str())))),  
				Transform(rotation1,location1),
                        	matter.updMobilizedBody(myChain[j].getAtomMobilizedBodyIndex(Compound::AtomIndex(myChain[j].getAtomIndex(ss4.str())))),    
				Transform(rotation2,location2),
                        	myLeontisWesthofBondRow.springConstant[q],
                        	myLeontisWesthofBondRow.torqueConstant,
				dutyCycle,
				myParameterReader.scrubberPeriod,
				myParameterReader.cutoffRadius
				)
			);
			}	
                }
	    if (myParameterReader.verbose) cout<<"[BasePairingBonds.h ] check 10"<<endl;
            }
        }



};

