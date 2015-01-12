/* -------------------------------------------------------------------------- *
 *                           MMB (MacroMoleculeBuilder)                       *
 * -------------------------------------------------------------------------- *
 *                                                                            *
 * Copyright (c) 2011-12 by the Author.                                       *
 * Author: Samuel Flores                                                      *
 *                                                                            *
 * See RNABuilder.cpp for the copyright and usage agreement.                  *
 * -------------------------------------------------------------------------- */

#include "ElectrostaticPotentialGridForce.h"  
#include "DensityMap.h"
#include "time.h"
#include <stdlib.h>


ElectrostaticPotentialGridForce::ElectrostaticPotentialGridForce (SimbodyMatterSubsystem& matter,ParameterReader& myParameterReader,  DensityMap & myDensityMap ,DuMMForceFieldSubsystem & dumm , BiopolymerClassContainer & myBiopolymerClassContainer, std::ostream& outputStream ) : 
    matter(matter), myParameterReader(myParameterReader),  myDensityMap(myDensityMap) , dumm(dumm), myBiopolymerClassContainer(myBiopolymerClassContainer), outputStream(outputStream)
        {
    //myBiopolymerClassContainer.validateAtomInfoVectors(); // now done automatically upon fetching atomInfoVector
    };    

void ElectrostaticPotentialGridForce::calcForce(const State& state, Vector_<SpatialVec>& bodyForces,
            Vector_<Vec3>& particleForces, Vector& mobilityForces) const
        {
        double torque = 0.;
        //if (myParameterReader.applyHeavyAtomDensityForces)
        for (int i = 0; i < myParameterReader.electroDensityContainer.numDensityStretches(); i++) {
                String myChainID = myParameterReader.electroDensityContainer.getDensityStretch(i).getChain();
                BiopolymerClass & tempBiopolymerClass = myParameterReader.myBiopolymerClassContainer.updBiopolymerClass(myChainID );
                Biopolymer & tempBiopolymer =  myParameterReader.myBiopolymerClassContainer.updBiopolymerClass(myChainID ).updBiopolymer();
        //vector<AtomInfo> tempAtomInfoVector = tempBiopolymerClass.getAtomInfoVector();   
                vector<MMBAtomInfo> tempAtomInfoVector = tempBiopolymerClass.calcAtomInfoVector(myParameterReader.electroDensityContainer.getDensityStretch(i), matter, dumm );  
                double densitySum = 0.0; 
                for (int m = 0; m < (int)tempAtomInfoVector.size(); m++) {
                    MMBAtomInfo & tempAtomInfo = tempAtomInfoVector[m];
                    Vec3 myAtomLocation = tempBiopolymer.calcAtomLocationInGroundFrame(state, tempAtomInfo.compoundAtomIndex);
                    Vec3 myAtomForce = myDensityMap.calcInterpolatedFirstQuadrantGradient(myAtomLocation) * (myParameterReader.electroDensityForceConstant * (-tempAtomInfo.partialCharge));
                    // cout << "ElectroForce "<< tempAtomInfo.atomName << " " <<tempAtomInfo.partialCharge <<" "<< myDensityMap.getDensity(myAtomLocation) << endl;
                    densitySum += myDensityMap.getDensity(myAtomLocation);

                    bodyForces[tempAtomInfo.mobilizedBodyIndex] +=  SpatialVec(torque + (-((tempAtomInfo.mobilizedBody).getBodyTransform(state)).T()+ myAtomLocation) % myAtomForce, myAtomForce);
                } // of for m
                // cout.precision(5);
                // cout << "DensityMean "<< densitySum/(int)tempAtomInfoVector.size() << endl;
        } // of for biopolymer
        };

Real ElectrostaticPotentialGridForce::calcPotentialEnergy(const State& state) const {

        Real totalPotentialEnergy = 0;
        for (int i = 0; i < myParameterReader.electroDensityContainer.numDensityStretches(); i++) {
                String myChainID = myParameterReader.electroDensityContainer.getDensityStretch(i).getChain();

                for (ResidueID j = myParameterReader.electroDensityContainer.getDensityStretch(i).getStartResidue(); j <=  myParameterReader.electroDensityContainer.getDensityStretch(i).getEndResidue(); myParameterReader.myBiopolymerClassContainer.updBiopolymerClass(myChainID).incrementResidueID( j) ) {
                        ResidueInfo myResidueInfo = myParameterReader.myBiopolymerClassContainer.updBiopolymerClass(myChainID).updBiopolymer().updResidue(ResidueInfo::Index(myParameterReader.myBiopolymerClassContainer.updBiopolymerClass(myChainID).getResidueIndex(j) ));

                        for (ResidueInfo::AtomIndex k ( 0); k < myResidueInfo.getNumAtoms() ; k++) {
                                Compound::AtomName myAtomName = myResidueInfo.getAtomName(k);
                                Compound::AtomIndex myAtomIndex = myResidueInfo.getAtomIndex( k  );
                                DuMM::AtomIndex myDuMMAtomIndex = myParameterReader.myBiopolymerClassContainer.updBiopolymerClass(myChainID).updBiopolymer().getDuMMAtomIndex(myAtomIndex);
                                if (((myParameterReader.myBiopolymerClassContainer.updBiopolymerClass(myChainID).updBiopolymer().getAtomElement(myAtomIndex)).getSymbol()).compare("H") != 0) { 
                                        Vec3 myAtomLocation = myParameterReader.myBiopolymerClassContainer.updBiopolymerClass(myChainID).updBiopolymer().calcAtomLocationInGroundFrame(state, myAtomIndex);
                                        totalPotentialEnergy -= myDensityMap.getDensity(myAtomLocation) * myParameterReader.densityForceConstant * (-dumm.getPartialCharge(myDuMMAtomIndex));

                                }
                        }
                    if (j ==  myParameterReader.electroDensityContainer.getDensityStretch(i).getEndResidue()) break;
                }

        } 

        return totalPotentialEnergy;
    };

bool ElectrostaticPotentialGridForce::dependsOnlyOnPositions() const  { 
        return true; 
    }; 