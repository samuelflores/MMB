/* -------------------------------------------------------------------------- *
 *                           MMB (MacroMoleculeBuilder)                       *
 * -------------------------------------------------------------------------- *
 *                                                                            *
 * Copyright (c) 2011-12 by the Author.                                       *
 * Author: Samuel Flores                                                      *
 *                                                                            *
 * See RNABuilder.cpp for the copyright and usage agreement.                  *
 * -------------------------------------------------------------------------- */

#include "DensityForce.h"  
#include "time.h"
#include <stdlib.h>

DensityForce::DensityForce (SimbodyMatterSubsystem& matter,ParameterReader& myParameterReader,  DensityMap & myDensityMap ,DuMMForceFieldSubsystem & dumm , BiopolymerClassContainer & myBiopolymerClassContainer, std::ostream& outputStream ) : 
	matter(matter), myParameterReader(myParameterReader),  myDensityMap(myDensityMap) , dumm(dumm), myBiopolymerClassContainer(myBiopolymerClassContainer), outputStream(outputStream)
        {
    //myBiopolymerClassContainer.validateAtomInfoVectors(); // now done automatically upon fetching atomInfoVector
    };    

void DensityForce::calcForce(const State& state, Vector_<SpatialVec>& bodyForces,
            Vector_<Vec3>& particleForces, Vector& mobilityForces) const
        {
        double torque = 0.;
        //if (myParameterReader.applyHeavyAtomDensityForces)
        for (int i = 0; i < myParameterReader.densityContainer.numDensityStretches(); i++) {
                String myChainID = myParameterReader.densityContainer.getDensityStretch(i).getChain();
                BiopolymerClass & tempBiopolymerClass = myParameterReader.myBiopolymerClassContainer.updBiopolymerClass(myChainID );
                Biopolymer & tempBiopolymer =  myParameterReader.myBiopolymerClassContainer.updBiopolymerClass(myChainID ).updBiopolymer();
                vector<MMBAtomInfo> tempAtomInfoVector = tempBiopolymerClass.calcAtomInfoVector(myParameterReader.densityContainer.getDensityStretch(i), matter, dumm, myParameterReader.densityFitPhosphates );   
                for (int m = 0; m < (int)tempAtomInfoVector.size(); m++) {
                    MMBAtomInfo & tempAtomInfo = tempAtomInfoVector[m];
                    Vec3 myAtomLocation = tempBiopolymer.calcAtomLocationInGroundFrame(state, tempAtomInfo.compoundAtomIndex);
                    //Vec3 myAtomForce = myDensityMap.calcInterpolatedFirstQuadrantGradient(myAtomLocation) * (myParameterReader.densityForceConstant * tempAtomInfo.mass);
                    // changed to atomic number on May 30 2012, earlier was atomic mass:
                    Vec3 myAtomForce = myDensityMap.calcInterpolatedFirstQuadrantGradient(myAtomLocation) * (myParameterReader.densityForceConstant * tempAtomInfo.atomicNumber);

                    bodyForces[tempAtomInfo.mobilizedBodyIndex] +=  SpatialVec(torque + (-((tempAtomInfo.mobilizedBody).getBodyTransform(state)).T()+ myAtomLocation) % myAtomForce, myAtomForce);
                } // of for m
        } // of for biopolymer
        };

Real DensityForce::calcPotentialEnergy(const State& state) const
        {
        double totalPotentialEnergy = 0;
        for (int i = 0; i < myParameterReader.densityContainer.numDensityStretches(); i++) {
                String myChainID = myParameterReader.densityContainer.getDensityStretch(i).getChain();
                BiopolymerClass & tempBiopolymerClass = myParameterReader.myBiopolymerClassContainer.updBiopolymerClass(myChainID );
                Biopolymer & tempBiopolymer =  myParameterReader.myBiopolymerClassContainer.updBiopolymerClass(myChainID ).updBiopolymer();
                vector<MMBAtomInfo> tempAtomInfoVector = tempBiopolymerClass.calcAtomInfoVector(myParameterReader.densityContainer.getDensityStretch(i), matter, dumm, myParameterReader.densityFitPhosphates );   
                for (int m = 0; m < (int)tempAtomInfoVector.size(); m++) {
                    MMBAtomInfo & tempAtomInfo = tempAtomInfoVector[m];
                    Vec3 myAtomLocation = tempBiopolymer.calcAtomLocationInGroundFrame(state, tempAtomInfo.compoundAtomIndex);
                    // changed to atomic number on FEB 24 2021, earlier was atomic mass:
                    totalPotentialEnergy -= myDensityMap.getDensity(myAtomLocation) * myParameterReader.densityForceConstant * tempAtomInfo.atomicNumber;
                } // of for m
        } // of for biopolymer
	return totalPotentialEnergy;
        };
/*
Real DensityForce::calcPotentialEnergy(const State& state) const {

        Real totalPotentialEnergy = 0;
        for (int i = 0; i < myParameterReader.densityContainer.numDensityStretches(); i++) {
                String myChainID = myParameterReader.densityContainer.getDensityStretch(i).getChain();

                for (ResidueID j = myParameterReader.densityContainer.getDensityStretch(i).getStartResidue(); j <=  myParameterReader.densityContainer.getDensityStretch(i).getEndResidue(); myParameterReader.myBiopolymerClassContainer.updBiopolymerClass(myChainID).incrementResidueID( j) ) {
                        ResidueInfo myResidueInfo = myParameterReader.myBiopolymerClassContainer.updBiopolymerClass(myChainID).updBiopolymer().updResidue(ResidueInfo::Index(myParameterReader.myBiopolymerClassContainer.updBiopolymerClass(myChainID).getResidueIndex(j) ));

                        for (ResidueInfo::AtomIndex k ( 0); k < myResidueInfo.getNumAtoms() ; k++) {
                                //Compound::AtomName myAtomName = myResidueInfo.getAtomName(k);
                                Compound::AtomIndex myAtomIndex = myResidueInfo.getAtomIndex( k  );
                                DuMM::AtomIndex myDuMMAtomIndex = myParameterReader.myBiopolymerClassContainer.updBiopolymerClass(myChainID).updBiopolymer().getDuMMAtomIndex(myAtomIndex);
                                if (((myParameterReader.myBiopolymerClassContainer.updBiopolymerClass(myChainID).updBiopolymer().getAtomElement(myAtomIndex)).getSymbol()).compare("H") != 0) { 
                                        Vec3 myAtomLocation = myParameterReader.myBiopolymerClassContainer.updBiopolymerClass(myChainID).updBiopolymer().calcAtomLocationInGroundFrame(state, myAtomIndex);
                                        // This should be off by default. When active, it reports the density observed and the atom name. The intent is to gather statistics on density maps to judiciously choose atomicNumber overrides. 
                                        if (myParameterReader.densityReportAtEachAtomPosition) { std::cout <<__FILE__<<":"<<__LINE__<< ":" << __FUNCTION__<<" Atom_name= "<< myResidueInfo.getAtomName(k)<<" local_density= "<<myDensityMap.getDensity(myAtomLocation)<<std::endl; }
                                        totalPotentialEnergy -= myDensityMap.getDensity(myAtomLocation) * myParameterReader.densityForceConstant * dumm.getAtomMass(myDuMMAtomIndex);

                                }
                        }
                    if (j ==  myParameterReader.densityContainer.getDensityStretch(i).getEndResidue()) break;
                }

        } 
        MMBLOG_FILE_FUNC_LINE(INFO, " Total potential energy due to density fitting potential = "<<totalPotentialEnergy <<std::endl);
        return totalPotentialEnergy;
    };
*/
bool DensityForce::dependsOnlyOnPositions() const  { 
        return true; 
    };    
