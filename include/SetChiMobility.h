#ifndef SetChiMobility_H_
#define SetChiMobility_H_
/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Molmodel                            *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2006-7 Stanford University and the Authors.         *
 * Authors: Samuel Flores                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

//using std::cout;
//using std::endl;

class  SetChiMobility  : public Biopolymer { public:

    SetChiMobility (   //this polymorphism, with one less parameter, sets up welds to rigidify the chi bond without losing the body frames of the glycosidic nitrogens.
        Biopolymer & myChain,
	int residueNumber,
	LeontisWesthofClass  & myLeontisWesthofClass ,   
        SimbodyMatterSubsystem & matter,
        State& state,
        Constraint myWeld1,
        Constraint myWeld2,
        Constraint myWeld3 
        //BondMobility::Mobility mobility = BondMobility::Torsion 
	)   
	{     
            BondMobility::Mobility mobility = BondMobility::Torsion;
	    String myPdbResidueName1 = (myChain.updResidue(ResidueInfo::Index(residueNumber))).getPdbResidueName();
            LeontisWesthofBondRow myLeontisWesthofBondRow2 = myLeontisWesthofClass.getLeontisWesthofBondRow((myChain.updResidue(ResidueInfo::Index(residueNumber))).getPdbResidueNumber(), (myChain.updResidue(ResidueInfo::Index(residueNumber))).getPdbResidueNumber()  , myPdbResidueName1,"ChiBond","","","",0);	
            stringstream ss3;
            ss3<<residueNumber<<"/"<<"C1*";
	    stringstream ss4;
	    ss4<<residueNumber<<"/"<<(myLeontisWesthofBondRow2.residue1Atom[1]);
            myChain.setBondMobility(mobility,ss3.str(),ss4.str());

	    stringstream ss5;
            ss5<<residueNumber<<"/"<<"O4*";
            myChain.setBondMobility(BondMobility::Torsion,ss3.str(),ss5.str());

	    stringstream ss6;
            ss6<<residueNumber<<"/"<<"C2*";
               						  //C1*    //C2*
            myChain.setBondMobility(BondMobility::Torsion,ss3.str(),ss6.str());
 

            const SimTK::Transform& mytransform1 = (matter.updMobilizedBody(myChain.getAtomMobilizedBodyIndex(myChain.getAtomIndex(ss3.str())))).findBodyTransformInAnotherBody(state,(matter.updMobilizedBody(myChain.getAtomMobilizedBodyIndex(myChain.getAtomIndex(ss4.str())))));
            myWeld1 = Constraint::Weld  (
                matter.updMobilizedBody(myChain.getAtomMobilizedBodyIndex(myChain.getAtomIndex(ss3.str()))),
 	        Transform(Vec3(0)),
		(matter.updMobilizedBody(myChain.getAtomMobilizedBodyIndex(myChain.getAtomIndex(ss4.str())))),
                mytransform1
            );
            const SimTK::Transform& mytransform2 = (matter.updMobilizedBody(myChain.getAtomMobilizedBodyIndex(myChain.getAtomIndex(ss3.str())))).findBodyTransformInAnotherBody(state,(matter.updMobilizedBody(myChain.getAtomMobilizedBodyIndex(myChain.getAtomIndex(ss5.str())))));
            myWeld2 = Constraint::Weld (
                (matter.updMobilizedBody(myChain.getAtomMobilizedBodyIndex(myChain.getAtomIndex(ss3.str())))),
 	        Transform(Vec3(0)),
		(matter.updMobilizedBody(myChain.getAtomMobilizedBodyIndex(myChain.getAtomIndex(ss5.str())))),
                mytransform2
            );
           const SimTK::Transform& mytransform3 = (matter.updMobilizedBody(myChain.getAtomMobilizedBodyIndex(myChain.getAtomIndex(ss3.str())))).findBodyTransformInAnotherBody(state,(matter.updMobilizedBody(myChain.getAtomMobilizedBodyIndex(myChain.getAtomIndex(ss6.str())))));
            myWeld3 = Constraint::Weld (
                (matter.updMobilizedBody(myChain.getAtomMobilizedBodyIndex(myChain.getAtomIndex(ss3.str())))),
 	        Transform(Vec3(0)),
		(matter.updMobilizedBody(myChain.getAtomMobilizedBodyIndex(myChain.getAtomIndex(ss6.str())))),
                mytransform3
            );
        }
    SetChiMobility (
        Biopolymer & myChain,
	int residueNumber,
	LeontisWesthofClass  & myLeontisWesthofClass ,   
        BondMobility::Mobility mobility = BondMobility::Free 
	)   
	{     
            //std::cout<<"[SetChiMobility.h] mobility ="<<mobility<<std::endl;
	    String myPdbResidueName1 = (myChain.updResidue(ResidueInfo::Index(residueNumber))).getPdbResidueName();
            LeontisWesthofBondRow myLeontisWesthofBondRow2 = myLeontisWesthofClass.getLeontisWesthofBondRow((myChain.updResidue(ResidueInfo::Index(residueNumber))).getPdbResidueNumber(), (myChain.updResidue(ResidueInfo::Index(residueNumber))).getPdbResidueNumber(),myPdbResidueName1,"ChiBond","","","",0);	
            stringstream ss3;
            ss3<<residueNumber<<"/"<<"C1*";
	    stringstream ss4;
	    ss4<<residueNumber<<"/"<<(myLeontisWesthofBondRow2.residue1Atom[1]);
            
            myChain.setBondMobility(mobility,ss3.str(),ss4.str());
	    stringstream ss5;
            ss5<<residueNumber<<"/"<<"O4*";
            myChain.setBondMobility(BondMobility::Torsion,ss3.str(),ss5.str());
	    stringstream ss6;
            ss6<<residueNumber<<"/"<<"C2*";
            myChain.setBondMobility(BondMobility::Torsion,ss3.str(),ss6.str());

        }


};
#endif
