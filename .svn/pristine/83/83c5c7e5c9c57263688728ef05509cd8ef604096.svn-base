/* -------------------------------------------------------------------------- *
 *                           MMB (MacroMoleculeBuilder)                       *
 * -------------------------------------------------------------------------- *
 *                                                                            *
 * Copyright (c) 2011-12 by the Author.                                       *
 * Author: Samuel Flores                                                      *
 *                                                                            *
 * See RNABuilder.cpp for the copyright and usage agreement.                  *
 * -------------------------------------------------------------------------- */

#ifndef CalcTransformCorrection_H_
#define CalcTransformCorrection_H_

#include "ParameterReader.h"

using namespace SimTK;
using namespace std; 

    static void         newCalcAxes (const State& state,ParameterReader & myParameterReader, vector<Biopolymer> & myChain, LeontisWesthofBondRow myLeontisWesthofBondRow,int residueNumber1,int residueNumber2,int i,int j,Vec3 & xAxisVector1,Vec3 & yAxisVector1, Vec3 & zAxisVector1,Vec3 & xAxisVector2,Vec3 & yAxisVector2 , Vec3 & zAxisVector2,Vec3 & glycosidicNitrogenAtom1LocationInGround,Vec3 & glycosidicNitrogenAtom2LocationInGround, Vec3 & ring1CenterLocationInGround, Vec3 & ring2CenterLocationInGround)  {
 
            stringstream ss3;
            ss3<<residueNumber1<<"/"<<myLeontisWesthofBondRow.residue1Atom[0];
            stringstream ss4;
            ss4<<residueNumber2<<"/"<<(myLeontisWesthofBondRow.residue2Atom[0]);
          
            stringstream ss1first;
            ss1first<<residueNumber1<<"/"<<myLeontisWesthofBondRow.residue1Atom[1];
            stringstream ss1second;
            ss1second<<  residueNumber1<<"/"<<(myLeontisWesthofBondRow.residue1Atom[2]);
            stringstream ss1c1p   ;
            ss1c1p<<  residueNumber1<<"/C1*";
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
            glycosidicNitrogenAtom1LocationInGround = myChain[i].calcAtomLocationInGroundFrame(state,myChain[i].getAtomIndex(ss3.str()));
            if (myParameterReader.verbose) cout<<__FILE__<<":"<<__LINE__<<" : ss3.str() = "<<ss3.str()<<endl;	
            if (myParameterReader.verbose) cout<<__FILE__<<":"<<__LINE__<<" :myChain[i].getAtomIndex(ss3.str())  = "<<myChain[i].getAtomIndex(ss3.str())<<endl;	
            if (myParameterReader.verbose) cout<<__FILE__<<":"<<__LINE__<<" : i,j = "<<i<<" , "<<j<<endl;
            if (myParameterReader.verbose) cout<<__FILE__<<":"<<__LINE__<<" :myChain[j].getAtomIndex(ss4.str())  = "<<myChain[j].getAtomIndex(ss4.str())<<endl;	
            if (myParameterReader.verbose) cout<<__FILE__<<":"<<__LINE__<<" :myChain[j].calcAtomLocationInGroundFrame(state,myChain[j].getAtomIndex(ss4.str()))  = "<<myChain[j].calcAtomLocationInGroundFrame(state,myChain[j].getAtomIndex(ss4.str()))<<endl;	
            glycosidicNitrogenAtom2LocationInGround = myChain[j].calcAtomLocationInGroundFrame(state,myChain[j].getAtomIndex(ss4.str()));
            //cout<<__FILE__<<":"<<__LINE__<< myChain[j].calcAtomLocationInGroundFrame(state,myChain[j].getAtomIndex(ss4.str()))<<endl;
            if (myParameterReader.verbose) cout<<__FILE__<<":"<<__LINE__<<" : ss4.str() = "<<ss4.str()<<endl;	
            if (myParameterReader.verbose) cout<<__FILE__<<":"<<__LINE__<<" :myChain[j].getAtomIndex(ss4.str())  = "<<myChain[j].getAtomIndex(ss4.str())<<endl;	
            if (myParameterReader.verbose) cout<<__FILE__<<":"<<__LINE__<<" :myChain[j].calcAtomLocationInGroundFrame(state,myChain[j].getAtomIndex(ss4.str()))  = "<<myChain[j].calcAtomLocationInGroundFrame(state,myChain[j].getAtomIndex(ss4.str()))<<endl;	
            Vec3 firstRingAtomvector1 = myChain[i].calcAtomLocationInGroundFrame(state,myChain[i].getAtomIndex(ss1first.str()))  - glycosidicNitrogenAtom1LocationInGround;
            //Vec3 c1pRingAtomvector1 = myChain[i].calcAtomLocationInGroundFrame(state,myChain[i].getAtomIndex(ss1c1p.str()))  - glycosidicNitrogenAtom1LocationInGround;
            Vec3 secondRingAtomvector1 = myChain[i].calcAtomLocationInGroundFrame(state,myChain[i].getAtomIndex(ss1second.str()))  - glycosidicNitrogenAtom1LocationInGround;
            /*Vec3 c1pRingAtomvector1norm = c1pRingAtomvector1/c1pRingAtomvector1.norm();
            double alpha= asin((c1pRingAtomvector1norm%firstRingAtomvector1).norm()/c1pRingAtomvector1norm.norm()/firstRingAtomvector1.norm());
            double  beta= asin((c1pRingAtomvector1norm%secondRingAtomvector1).norm()/c1pRingAtomvector1norm.norm()/secondRingAtomvector1.norm());
            double A = sin(beta)*c1pRingAtomvector1norm.norm()/sin(180/Rad2Deg-alpha-beta)/firstRingAtomvector1.norm();
            double B = sin(alpha)*c1pRingAtomvector1norm.norm()/sin(180/Rad2Deg-alpha-beta)/secondRingAtomvector1.norm();
            */
            Vec3 firstRingAtomvector2 = myChain[j].calcAtomLocationInGroundFrame(state,myChain[j].getAtomIndex(ss2first.str()))  -glycosidicNitrogenAtom2LocationInGround;
            Vec3 secondRingAtomvector2 = myChain[j].calcAtomLocationInGroundFrame(state,myChain[j].getAtomIndex(ss2second.str()))  - glycosidicNitrogenAtom2LocationInGround;
            if (myParameterReader.verbose) cout<<__FILE__<<":"<<__LINE__<<" : firstRingAtomvector2, secondRingAtomvector2 "<<  firstRingAtomvector2<<" , "<< secondRingAtomvector2    <<endl;	
            if (myParameterReader.verbose) cout<<__FILE__<<":"<<__LINE__<<" : ss2second.str() = "<<  ss2second.str()  <<endl;	
            if (myParameterReader.verbose) cout<<__FILE__<<":"<<__LINE__<<" : myChain[j].getAtomIndex(ss2second.str())) = "<<myChain[j].getAtomIndex(ss2second.str())<<endl;
            if (myParameterReader.verbose) cout<<__FILE__<<":"<<__LINE__<<" :  myChain[j].calcAtomLocationInGroundFrame(state,myChain[j].getAtomIndex(ss2second.str()))    = "<<  myChain[j].calcAtomLocationInGroundFrame(state,myChain[j].getAtomIndex(ss2second.str())) <<endl; 
            if (myParameterReader.verbose) cout<<__FILE__<<":"<<__LINE__<<" :glycosidicNitrogenAtom2LocationInGround  = "<<glycosidicNitrogenAtom2LocationInGround   <<endl;	
            if ((myLeontisWesthofBondRow.pdbResidueName1.compare("A  ") == 0) || (myLeontisWesthofBondRow.pdbResidueName1.compare("G  ") == 0) ) { //if purine

                xAxisVector1 =  -5.88327 * firstRingAtomvector1 - 6.13617 * secondRingAtomvector1;          
                ssRing1Atom1<<"N3";
                ssRing1Atom2<<"C6";
                //cout <<"[TwoTransformForces.h] ssRing1Atom1, ssRing1Atom2 : "<<ssRing1Atom1.str()<<", "<<ssRing1Atom2.str()<< endl;
                ring1CenterLocationInGround = (myChain[i].calcAtomLocationInGroundFrame(state,myChain[i].getAtomIndex(ssRing1Atom1.str())) 
                                              +myChain[i].calcAtomLocationInGroundFrame(state,myChain[i].getAtomIndex(ssRing1Atom2.str())))/2 ;
                if (myParameterReader.verbose) cout <<"[TwoTransformForces.h] ring1CenterLocationInGround just computed for A or U:"<<ring1CenterLocationInGround<<endl;
            }
            else if ((myLeontisWesthofBondRow.pdbResidueName1.compare("C  ") == 0)) {
                xAxisVector1 = -7.83435 * firstRingAtomvector1 -6.99265          *secondRingAtomvector1;           
                ssRing1Atom1<<"N1";
                ssRing1Atom2<<"C4";
                ring1CenterLocationInGround = (myChain[i].calcAtomLocationInGroundFrame(state,myChain[i].getAtomIndex(ssRing1Atom1.str())) 
                                              +myChain[i].calcAtomLocationInGroundFrame(state,myChain[i].getAtomIndex(ssRing1Atom2.str())))/2 ;
            }
            else if ((myLeontisWesthofBondRow.pdbResidueName1.compare("U  ")) == 0) {
                xAxisVector1 = -7.3491 * firstRingAtomvector1 -6.47606 *secondRingAtomvector1;     
                ssRing1Atom1<<"N1";
                ssRing1Atom2<<"C4";
                ring1CenterLocationInGround = (myChain[i].calcAtomLocationInGroundFrame(state,myChain[i].getAtomIndex(ssRing1Atom1.str())) 
                                              +myChain[i].calcAtomLocationInGroundFrame(state,myChain[i].getAtomIndex(ssRing1Atom2.str())))/2 ;
                if (myParameterReader.verbose) cout <<"[TwoTransformForces.h] ring1CenterLocationInGround just computed for U: "<<ring1CenterLocationInGround<<endl;
            } 
            else { cout <<"[TwoTransformForces.h] Unrecognized residue type: "<<myLeontisWesthofBondRow.pdbResidueName1<<endl; assert(0);} // trap errors
            if ((myLeontisWesthofBondRow.pdbResidueName2.compare("A  ")  == 0) || (myLeontisWesthofBondRow.pdbResidueName2.compare("G  ") == 0)){ //if purine
                xAxisVector2 = -5.88327 * firstRingAtomvector2 -6.13617 *secondRingAtomvector2;        
                ssRing2Atom1<<"N3";
                ssRing2Atom2<<"C6";
                ring2CenterLocationInGround = (myChain[j].calcAtomLocationInGroundFrame(state,myChain[j].getAtomIndex(ssRing2Atom1.str())) 
                                              +myChain[j].calcAtomLocationInGroundFrame(state,myChain[j].getAtomIndex(ssRing2Atom2.str())))/2 ;
            }   
            else if ((myLeontisWesthofBondRow.pdbResidueName2.compare("C  ") == 0)){
		ssRing2Atom1<<"N1";
		ssRing2Atom2<<"C4";
		ring2CenterLocationInGround = (myChain[j].calcAtomLocationInGroundFrame(state,myChain[j].getAtomIndex(ssRing2Atom1.str())) 
					      +myChain[j].calcAtomLocationInGroundFrame(state,myChain[j].getAtomIndex(ssRing2Atom2.str())))/2 ;
                xAxisVector2 = -7.83435 * firstRingAtomvector2 -6.99265 *secondRingAtomvector2;           
            }
            else if ((myLeontisWesthofBondRow.pdbResidueName2.compare("U  ")) == 0) {
                xAxisVector2 = -7.3491  * firstRingAtomvector2 -6.47606 *secondRingAtomvector2;           
		ssRing2Atom1<<"N1";
		ssRing2Atom2<<"C4";
                if (myParameterReader.verbose) cout <<"[TwoTransformForces.h] ssRing2Atom1, ssRing2Atom2 : "<<ssRing2Atom1.str()<<", "<<ssRing2Atom2.str()<<                   endl;
		ring2CenterLocationInGround = (myChain[j].calcAtomLocationInGroundFrame(state,myChain[j].getAtomIndex(ssRing2Atom1.str())) 
					      +myChain[j].calcAtomLocationInGroundFrame(state,myChain[j].getAtomIndex(ssRing2Atom2.str())))/2 ;
                if (myParameterReader.verbose) cout <<"[TwoTransformForces.h] ring2CenterLocationInGround just computed for U: "<<ring2CenterLocationInGround<<endl;
            }
            else { cout <<"[TwoTransformForces.h] Unrecognized residue type"<<endl; assert(0);} // trap errors
            zAxisVector1 = (firstRingAtomvector1%secondRingAtomvector1);
            zAxisVector1 = zAxisVector1/zAxisVector1.norm(); 
            zAxisVector2 = (firstRingAtomvector2%secondRingAtomvector2);
            zAxisVector2 = zAxisVector2/zAxisVector2.norm(); 
            yAxisVector1 = zAxisVector1%xAxisVector1;
            yAxisVector1= yAxisVector1/yAxisVector1.norm();
            yAxisVector2 = zAxisVector2%xAxisVector2;
            yAxisVector2= yAxisVector2/yAxisVector2.norm();
            //cout<<__FILE__<<__LINE__<<xAxisVector1<<yAxisVector1<<zAxisVector1<<endl;
            //cout<<__FILE__<<__LINE__<<xAxisVector2<<yAxisVector2<<zAxisVector2<<endl;

};


static void computeCorrection (LeontisWesthofClass & myLeontisWesthofClass, ParameterReader & myParameterReader , vector<Biopolymer> & myChain, State & state, SimbodyMatterSubsystem& matter   ) {
            for (int i = 0; i<myParameterReader.baseOperationVector.size(); i++) 
            if (
                ((myParameterReader.baseOperationVector[i]).BasePairIsTwoTransformForce.compare("baseInteraction")==0 )  &&
                ((
                ((myParameterReader.sequenceTypes[myParameterReader.baseOperationVector[i].FirstBPChain]).compare("rna") == 0) ||
                ((myParameterReader.sequenceTypes[myParameterReader.baseOperationVector[i].FirstBPChain]).compare("CoarseNucleicAcid") == 0) 
                ) && (
                ((myParameterReader.sequenceTypes[myParameterReader.baseOperationVector[i].SecondBPChain]).compare("rna") == 0) || 
                ((myParameterReader.sequenceTypes[myParameterReader.baseOperationVector[i].SecondBPChain]).compare("CoarseNucleicAcid") == 0)  

                ))
               )
            {
                String myResidueName1 =       
                    getResidueName(myParameterReader.baseOperationVector[i].FirstBPChain, myParameterReader.baseOperationVector[i].FirstBPResidue,myParameterReader,myChain ); 
                String myResidueName2 =       
                    getResidueName(myParameterReader.baseOperationVector[i].SecondBPChain, myParameterReader.baseOperationVector[i].SecondBPResidue,myParameterReader,myChain ) ;
                LeontisWesthofBondRow myLeontisWesthofBondRow = myLeontisWesthofClass.getLeontisWesthofBondRow(
                    myParameterReader.baseOperationVector[i].FirstBPResidue, myParameterReader.baseOperationVector[i].SecondBPResidue, 
                    myResidueName1, //getResidueName(myParameterReader.baseOperationVector[i].FirstBPChain, myParameterReader.baseOperationVector[i].FirstBPResidue,myParameterReader,myChain ), 
                    myParameterReader.baseOperationVector[i].FirstBPEdge,   
                    myResidueName2, //getResidueName(myParameterReader.baseOperationVector[i].SecondBPChain, myParameterReader.baseOperationVector[i].SecondBPResidue,myParameterReader,myChain ) ,
                    myParameterReader.baseOperationVector[i].SecondBPEdge ,
                    myParameterReader.baseOperationVector[i].OrientationBP ,   myParameterReader.baseOperationVector[i].BasePairIsTwoTransformForce
                );
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
                ResidueInfo::Index myResidue1 (
                    myParameterReader.baseOperationVector[i].FirstBPResidue- myParameterReader.getFirstResidueNumbers(myParameterReader.baseOperationVector[i].FirstBPChain ));
                int myChain1Index = myParameterReader.getChainIndex(myParameterReader.baseOperationVector[i].FirstBPChain , myChain)  ;
                MobilizedBody body1 = matter.updMobilizedBody(myChain[
                        myChain1Index 
                        ]
                    .getAtomMobilizedBodyIndex(
                    Compound::AtomIndex(myChain[
                        myChain1Index 
                        ]
                    .getResidue(
                    myResidue1 
                    )
                    .getAtomIndex(myLeontisWesthofBondRow.residue1Atom[0]))));

                ResidueInfo::Index myResidue2 (
                    myParameterReader.baseOperationVector[i].SecondBPResidue- myParameterReader.getFirstResidueNumbers(myParameterReader.baseOperationVector[i].SecondBPChain ));
                int myChain2Index = myParameterReader.getChainIndex(myParameterReader.baseOperationVector[i].SecondBPChain , myChain)  ;
                MobilizedBody body2 = matter.updMobilizedBody(myChain[
                     myChain2Index
                     //myParameterReader.getChainIndex(myParameterReader.baseOperationVector[i].SecondBPChain , myChain)
                     ]
                    .getAtomMobilizedBodyIndex(Compound::AtomIndex(myChain[
                        myChain2Index
                        //myParameterReader.getChainIndex(myParameterReader.baseOperationVector[i].SecondBPChain , myChain)  
                        ]
                    .getResidue(
                        myResidue2
                        //ResidueInfo::Index(myParameterReader.baseOperationVector[i].SecondBPResidue- myParameterReader.getFirstResidueNumbers((myParameterReader.baseOperationVector[i].SecondBPChain) ))
                    )
                    .getAtomIndex(myLeontisWesthofBondRow.residue2Atom[0]))));

     
            	newCalcAxes(state,myParameterReader,myChain,  myLeontisWesthofBondRow, myParameterReader.baseOperationVector[i].FirstBPResidue- myParameterReader.getFirstResidueNumbers(myParameterReader.baseOperationVector[i].FirstBPChain) , myParameterReader.baseOperationVector[i].SecondBPResidue-myParameterReader.getFirstResidueNumbers(myParameterReader.baseOperationVector[i].SecondBPChain)  , 
                   myParameterReader.getChainIndex(myParameterReader.baseOperationVector[i].FirstBPChain,myChain ),
                   myParameterReader.getChainIndex(myParameterReader.baseOperationVector[i].SecondBPChain , myChain) , 
                   xAxisVector1,yAxisVector1,zAxisVector1,xAxisVector2,yAxisVector2,zAxisVector2,
                   glycosidicNitrogenAtom1LocationInGround,glycosidicNitrogenAtom2LocationInGround,
                   ring1CenterLocationInGround,ring2CenterLocationInGround);


            	Rotation rotation1FromRingAtoms(Mat33(xAxisVector1,yAxisVector1,zAxisVector1));
                if (myParameterReader.verbose) cout<<__FILE__<<":"<<__LINE__<<" :xAxisVector2, yAxisVector2, zAxisVector2 = "<< xAxisVector2<<" , "<< yAxisVector2<<" , "<< zAxisVector2   <<endl;	
            	Rotation rotation2FromRingAtoms(Mat33(xAxisVector2,yAxisVector2,zAxisVector2));
                Rotation    myRotationCorrection1 = ~rotation1FromRingAtoms * ( matter.getMobilizedBody(body1).getBodyTransform(state)).R();

 
                Rotation myRotationCorrection2 = (~rotation2FromRingAtoms * ( matter.getMobilizedBody(body2).getBodyTransform(state)).R()); 
                Vec3 myTranslationCorrection1 = (~( matter.getMobilizedBody(body1).getBodyTransform(state)).R()*(glycosidicNitrogenAtom1LocationInGround - ( matter.getMobilizedBody(body1).getBodyTransform(state)).T()  ));
                // 
                Vec3 myTranslationCorrection2 = (~( matter.getMobilizedBody(body2).getBodyTransform(state)).R()*(glycosidicNitrogenAtom2LocationInGround - ( matter.getMobilizedBody(body2).getBodyTransform(state)).T()  ));
                (myParameterReader.baseOperationVector[i]).rotationCorrection1 =myRotationCorrection1;
                (myParameterReader.baseOperationVector[i]).rotationCorrection2 = myRotationCorrection2;
                (myParameterReader.baseOperationVector[i]).translationCorrection1 = myTranslationCorrection1;
                (myParameterReader.baseOperationVector[i]).translationCorrection2 = myTranslationCorrection2;
    }
};
#endif
