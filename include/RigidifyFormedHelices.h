#ifndef RigidifyFormedHelices_H_
#define RigidifyFormedHelices_H_
/* -------------------------------------------------------------------------- *
 *                           MMB (MacroMoleculeBuilder)                       *
 * -------------------------------------------------------------------------- *
 *                                                                            *
 * Copyright (c) 2011-12 by the Author.                                       *
 * Author: Samuel Flores                                                      *
 *                                                                            *
 * See RNABuilder.cpp for the copyright and usage agreement.                  *
 * -------------------------------------------------------------------------- */
#include <stdio.h>
#include <string.h>

//using std::cout;
//using std::endl;

//using namespace SimTK;
//using namespace std; 




class  RigidifyFormedHelices   { public:
    //this polymorphism, with one less parameter, sets up welds to rigidify the chi bond without losing the body frames of the glycosidic nitrogens.
    RigidifyFormedHelices (  
        ParameterReader & myParameterReader,
        vector<Biopolymer>   myMolecule, //Chains]),
	BondMobility::Mobility helixBondMobility
	)   
	{     
            cout <<"[RigidifyFormedHelices.h] inside helix rigidification procedure "<<endl;
            // if the endresidue of the rigid segment is partnered with a lower-numbered residue
            for (int j = 0; j<(myParameterReader.sequences).size(); j++) 
            for (int q=0;q<myParameterReader.numRigidSegments[(myParameterReader.chainId[j]).c_str()]; q++) {
            if (myParameterReader.basePairPartners[ChainResidueIndex(myParameterReader.basePairPartners[ChainResidueIndex(q,j)].rigidSegmentEndNumber,j)].BPPartner < myParameterReader.basePairPartners[ChainResidueIndex(q,j)].rigidSegmentEndNumber) {
                // find the RigidSegment that this RigidSegment is partnered with 
                // if priority is lower than current, loop over all base pairs in this RigidSegment AND ITS PARTNER and set residue BondMobility to helixBondMobility
                for (int i = 0; i<myParameterReader.baseOperationVector.size(); i++) {
                    //loop over all base pairs until we find the last/ highest-numbered base pair of the curent rigid segment   
                    if (((
			(myParameterReader.baseOperationVector[i].FirstBPResidue == myParameterReader.basePairPartners[ChainResidueIndex(q,j)].rigidSegmentEndNumber) 
			&& (myParameterReader.baseOperationVector[i].SecondBPResidue    == myParameterReader.basePairPartners[ChainResidueIndex(myParameterReader.basePairPartners[ChainResidueIndex(q,j)].rigidSegmentEndNumber-myParameterReader.getFirstResidueNumbers((myParameterReader.chainId[j]).c_str()) ,j )].BPPartner)
			) || (
			(myParameterReader.baseOperationVector[i].SecondBPResidue    == myParameterReader.basePairPartners[ChainResidueIndex(q,j)].rigidSegmentEndNumber) 
			&& (myParameterReader.baseOperationVector[i].FirstBPResidue == myParameterReader.basePairPartners[ChainResidueIndex(myParameterReader.basePairPartners[ChainResidueIndex(q,j)].rigidSegmentEndNumber-myParameterReader.getFirstResidueNumbers((myParameterReader.chainId[j]).c_str()) ,j ) ].BPPartner)
			)) && ((strcmp((myParameterReader.baseOperationVector[i].FirstBPEdge    ).c_str(),string("WatsonCrick").c_str()) == 0) ) && ((strcmp((myParameterReader.baseOperationVector[i].SecondBPEdge   ).c_str(),string("WatsonCrick").c_str()) == 0))
			&& 	((myParameterReader.basePairPartners[ChainResidueIndex(myParameterReader.basePairPartners[ChainResidueIndex(q,j)].rigidSegmentEndNumber-myParameterReader.getFirstResidueNumbers((myParameterReader.chainId[j]).c_str()) ,j )].BPPartner
                        ) < 
			myParameterReader.basePairPartners[ChainResidueIndex(q,j)].rigidSegmentEndNumber)
			)  

		    {
                    if (myParameterReader.verbose) cout<<"[RigidifyFormedHelices.h] myParameterReader.chainId[j] ="<<myParameterReader.chainId[j]<<endl;
                        assert(strcmp((myParameterReader.baseOperationVector[i].FirstBPEdge    ).c_str(),string("WatsonCrick").c_str()) == 0); // call me paranoid, want to make sure this is a watson-crick base pair.
                        assert(strcmp((myParameterReader.baseOperationVector[i].SecondBPEdge   ).c_str(),string("WatsonCrick").c_str()) == 0); // call me paranoid, want to make sure this is a watson-crick base pair.
                        // are we sure BPPartner is strictly the WatsonCrick partner?  the above should check.
                        // if the rigidSegment was completed before the current priority level, proceed with setting BondMobilities to helixBondMobility 
                        if (myParameterReader.baseOperationVector[i].BasePairPriority    < myParameterReader.priority) {
                            RNA& rna = static_cast<RNA&>(myMolecule[j]);   
                            rna.setRNABondMobility(   helixBondMobility,SimTK::ResidueInfo::Index(myParameterReader.basePairPartners[ChainResidueIndex(q,j)].rigidSegmentStartNumber-myParameterReader.getFirstResidueNumbers((myParameterReader.chainId[j]).c_str())),SimTK::ResidueInfo::Index(myParameterReader.basePairPartners[ChainResidueIndex(q,j)].rigidSegmentEndNumber-myParameterReader.getFirstResidueNumbers((myParameterReader.chainId[j]).c_str()))); 
                            rna.setRNABondMobility(   helixBondMobility,
			        SimTK::ResidueInfo::Index(myParameterReader.basePairPartners[ChainResidueIndex((myParameterReader.basePairPartners[ChainResidueIndex(q,j)].rigidSegmentEndNumber-myParameterReader.getFirstResidueNumbers((myParameterReader.chainId[j]).c_str())) ,j)].BPPartner  
			        -myParameterReader.getFirstResidueNumbers((myParameterReader.chainId[j]).c_str())

                                )  
                                , 
			        SimTK::ResidueInfo::Index(myParameterReader.basePairPartners[ChainResidueIndex((myParameterReader.basePairPartners[ChainResidueIndex(q,j)].rigidSegmentStartNumber-myParameterReader.getFirstResidueNumbers((myParameterReader.chainId[j]).c_str())) ,j )].BPPartner
			        - myParameterReader.getFirstResidueNumbers((myParameterReader.chainId[j]).c_str())  )
                            ); 
			    myMolecule[j]  = rna;
    }

                    }
                }
            }
                 //Finally, we must Weld the two RigidSegments together, since the forces will no longer work.  We'll do that further down in this program though..

        }


}

    RigidifyFormedHelices ( 
        ParameterReader myParameterReader,
        RNA & myMolecule,
	SimbodyMatterSubsystem& matter,
        CompoundSystem & system,
        State & state,
        Constraint myWeld1[maxChiWelds]
        ) 
{
        
            for (int j = 0; j<(myParameterReader.sequences).size(); j++) 
	     for (int q=0;q<myParameterReader.numRigidSegments[(myParameterReader.chainId[j]).c_str()]; q++)
	      if (myParameterReader.basePairPartners[ChainResidueIndex(myParameterReader.basePairPartners[ChainResidueIndex(q,j)].rigidSegmentEndNumber ,j )].BPPartner < myParameterReader.basePairPartners[ChainResidueIndex(q,j)].rigidSegmentEndNumber) 
                for (int i = 0; i<myParameterReader.baseOperationVector.size(); i++) {
                    //loop over all base pairs until we find the last/ highest-numbered base pair of the curent rigid segment
                    if (((
                        (myParameterReader.baseOperationVector[i].FirstBPResidue    == myParameterReader.basePairPartners[ChainResidueIndex(q,j)].rigidSegmentEndNumber)
                        && (myParameterReader.baseOperationVector[i].SecondBPResidue    == myParameterReader.basePairPartners[ChainResidueIndex(myParameterReader.basePairPartners[ChainResidueIndex(q,j)].rigidSegmentEndNumber-myParameterReader.getFirstResidueNumbers((myParameterReader.chainId[j]).c_str()) ,j )].BPPartner)
                        ) || (
                        (myParameterReader.baseOperationVector[i].SecondBPResidue    == myParameterReader.basePairPartners[ChainResidueIndex(q,j)].rigidSegmentEndNumber)
                        && (myParameterReader.baseOperationVector[i].FirstBPResidue == myParameterReader.basePairPartners[ChainResidueIndex(myParameterReader.basePairPartners[ChainResidueIndex(q,j)].rigidSegmentEndNumber-myParameterReader.getFirstResidueNumbers((myParameterReader.chainId[j]).c_str()) ,j )].BPPartner)                        )) && ((strcmp((myParameterReader.baseOperationVector[i].FirstBPEdge    ).c_str(),string("WatsonCrick").c_str()) == 0) ) && ((strcmp((myParameterReader.baseOperationVector[i].SecondBPEdge   ).c_str(),string("WatsonCrick").c_str()) == 0))                        &&      ((myParameterReader.basePairPartners[ChainResidueIndex(myParameterReader.basePairPartners[ChainResidueIndex(q,j)].rigidSegmentEndNumber-myParameterReader.getFirstResidueNumbers((myParameterReader.chainId[j]).c_str()) , j)].BPPartner) < 
                        myParameterReader.basePairPartners[ChainResidueIndex(q,j)].rigidSegmentEndNumber)
                        && (myParameterReader.baseOperationVector[i].BasePairPriority    < myParameterReader.priority)
                        )   

                                {
                                    //weld a residue from r with a residue from q
                                    stringstream ss3;
                                    ss3<<myParameterReader.basePairPartners[ChainResidueIndex(q,j)].rigidSegmentEndNumber-myParameterReader.getFirstResidueNumbers((myParameterReader.chainId[j]).c_str()) <<"/"<<"C3*";
                                    stringstream ss4;
                                    ss4<<myParameterReader.basePairPartners[ChainResidueIndex((myParameterReader.basePairPartners[ChainResidueIndex(q,j)].rigidSegmentEndNumber-myParameterReader.getFirstResidueNumbers((myParameterReader.chainId[j]).c_str())) ,j )].BPPartner - 
			 	        myParameterReader.getFirstResidueNumbers((myParameterReader.chainId[j]).c_str()) <<"/"<<"C3*";
                                            
                                    if (myParameterReader.verbose) cout<<"[RigidifyFormedHelices.h]myParameterReader.getFirstResidueNumbers((myParameterReader.chainId[j]).c_str()) ="<<myParameterReader.getFirstResidueNumbers((myParameterReader.chainId[j]).c_str())<<endl;
                                    if (myParameterReader.verbose) cout<<"[RigidifyFormedHelices.h] ss3,ss4 :"<<ss3.str()<<","<<ss4.str()<<endl;
                                    	
                                    if (myParameterReader.verbose) cout<<"[RigidifyFormedHelices.h] check .3 "<< ((myMolecule.getAtomIndex(ss3.str())))<<endl;
                                    if (myParameterReader.verbose) cout<<"[RigidifyFormedHelices.h] check .5 "<<endl;
                                    if (myParameterReader.verbose) cout<<"[RigidifyFormedHelices.h] check .5 "<< (myMolecule.getAtomMobilizedBodyIndex(myMolecule.getAtomIndex(ss3.str())))<<endl;
			            Transform      mytemp =  (matter.updMobilizedBody(myMolecule.getAtomMobilizedBodyIndex(myMolecule.getAtomIndex(ss3.str())))).getBodyTransform(state);
                                    if (myParameterReader.verbose) cout<<"[RigidifyFormedHelices.h] check 1  "<<endl;
			            MobilizedBody & mytemp2 =  (matter.updMobilizedBody(myMolecule.getAtomMobilizedBodyIndex(myMolecule.getAtomIndex(ss4.str()))));
                                    if (myParameterReader.verbose) cout<<"[RigidifyFormedHelices.h] check 1.5  "<<endl;

                                    const SimTK::Transform& mytransform1 = (matter.updMobilizedBody(myMolecule.getAtomMobilizedBodyIndex(myMolecule.getAtomIndex(ss3.str())))).findBodyTransformInAnotherBody(state,(matter.updMobilizedBody(myMolecule.getAtomMobilizedBodyIndex(myMolecule.getAtomIndex(ss4.str())))));
                                    if (myParameterReader.verbose) cout<<"[RigidifyFormedHelices.h] check 2  "<<endl;
    //Finally, we must Weld the two RigidSegments together, since the forces will no longer work.
                                    myWeld1[q] = Constraint::Weld  (
                                        matter.updMobilizedBody(myMolecule.getAtomMobilizedBodyIndex(myMolecule.getAtomIndex(ss3.str()))),
                                        Transform(Vec3(0)),
                                        (matter.updMobilizedBody(myMolecule.getAtomMobilizedBodyIndex(myMolecule.getAtomIndex(ss4.str())))),
                                        mytransform1
                                    );
				    state = system.realizeTopology();
        			    system.realize(state,Stage::Position);
	

                                }
                        }

                    }
};
#endif
