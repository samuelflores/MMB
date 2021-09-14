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
                const String &myChainID = myParameterReader.densityContainer.getDensityStretch(i).getChain();
		if (myParameterReader.myBiopolymerClassContainer.hasChainID(myChainID)){
                    const BiopolymerClass &bpc = myParameterReader.myBiopolymerClassContainer.getBiopolymerClass(myChainID );
                    const Biopolymer &bp = bpc.getBiopolymer();
                    const auto atomInfoRange = bpc.calcAtomInfoVector(myParameterReader.densityContainer.getDensityStretch(i), matter, dumm,myParameterReader.densityFitPhosphates); 
                    for (auto it = atomInfoRange.first; it != atomInfoRange.second; it++) {
                        Vec3 myAtomLocation = bp.calcAtomLocationInGroundFrame(state, it->compoundAtomIndex);
                        // changed to atomic number on May 30 2012, earlier was atomic mass:
                        const Vec3 &myAtomForce = myDensityMap.calcInterpolatedFirstQuadrantGradient(myAtomLocation) * (myDensityMap.getForceConstant() * it->atomicNumber);
                        bodyForces[it->mobilizedBodyIndex] +=  SpatialVec(torque + (-((it->mobilizedBody).getBodyTransform(state)).T()+ myAtomLocation) % myAtomForce, myAtomForce);
                    } // of for m
		} // of if myBiopolymerClassContainer.hasChainID 
		else if (myParameterReader.myMonoAtomsContainer.hasChainID(myChainID)){
                        // We have not yet implemented calcForce for monoAtoms, so do nothing here.
                    }			    
                else {
                        MMBLOG_FILE_FUNC_LINE(CRITICAL, " The chain ID you specified, "<< myChainID << " does not correspond to any existing BiopolymerClass or MonoAtoms !"<<endl);
		    }
        } // of for myParameterReader.densityContainer.numDensityStretches
        };

Real DensityForce::calcPotentialEnergy(const State& state) const
        {
        double totalPotentialEnergy = 0;
        for (int i = 0; i < myParameterReader.densityContainer.numDensityStretches(); i++) {
            const String &myChainID = myParameterReader.densityContainer.getDensityStretch(i).getChain();
            if (myParameterReader.myBiopolymerClassContainer.hasChainID(myChainID)){
                    const BiopolymerClass &bpc = myParameterReader.myBiopolymerClassContainer.getBiopolymerClass(myChainID );
                    const Biopolymer &bp = bpc.getBiopolymer();
                    const auto atomInfoRange = bpc.calcAtomInfoVector(myParameterReader.densityContainer.getDensityStretch(i), matter, dumm, myParameterReader.densityFitPhosphates );   
                    for (auto it = atomInfoRange.first; it != atomInfoRange.second; it++) {
                        Vec3 myAtomLocation = bp.calcAtomLocationInGroundFrame(state, it->compoundAtomIndex);
                        MMBLOG_FILE_FUNC_LINE(DEBUG,  " myDensityMap.getDensity(Vec3(-.68, -4.79, 2.26)) = "<< myDensityMap.getDensity(Vec3(-.68, -4.79, 2.26))<<endl);
                        MMBLOG_FILE_FUNC_LINE(DEBUG,  " myDensityMap.getDensity(myAtomLocation) = "<< myDensityMap.getDensity(myAtomLocation)<<" myAtomLocation = "<<myAtomLocation<<" myDensityMap.getForceConstant()        = "<<myDensityMap.getForceConstant()<<" myAtomicNumber = "<<it->atomicNumber<<endl);
                        // changed to atomic number on FEB 24 2021, earlier was atomic mass:
                        totalPotentialEnergy -= myDensityMap.getDensity(myAtomLocation) * myDensityMap.getForceConstant() * it->atomicNumber;
                    } // of for m
            } // of if myBiopolymerClassContainer.hasChainID	
            else if (myParameterReader.myMonoAtomsContainer.hasChainID(myChainID)){
		    
                for (int i = 0; i < myParameterReader.myMonoAtomsContainer.getMonoAtoms(myChainID).getNumAtoms(); i++) {
                    Vec3 myAtomLocation = myParameterReader.myMonoAtomsContainer.getMonoAtoms(myChainID).getAtomLocationInGroundFrame(i,state);
		    //ResidueID myResidueID = myParameterReader.myMonoAtomsContainer.getMonoAtoms(myChainID).getResidueID(i);
		    // Couldn't get getAtomElement to work .. just hard coding to unity for now.
		    //SimTK::Compound::AtomIndex    myAtomIndex = myParameterReader.myMonoAtomsContainer.getMonoAtoms(myChainID).getAtomIndex(myResidueID);

		    int    myAtomicNumber = 1; //dumm.getAtomElement(myAtomIndex);
                    //MMBLOG_FILE_FUNC_LINE(DEBUG,  " myDensityMap.getDensity(Vec3(myAtomLocation[0], -4.79, 2.26)) = "<< myDensityMap.getDensity(Vec3(myAtomLocation[0], -4.79, 2.26))<<endl);
                    //MMBLOG_FILE_FUNC_LINE(DEBUG,  " myDensityMap.getDensity(Vec3(-.68,myAtomLocation[1],  2.26)) = "<< myDensityMap.getDensity(Vec3(-.68,myAtomLocation[1], 2.26))<<endl);
                    //MMBLOG_FILE_FUNC_LINE(DEBUG,  " myDensityMap.getDensity(Vec3(-.68,-4.79,myAtomLocation[2])) = "<< myDensityMap.getDensity(Vec3(-.68,-4.79,myAtomLocation[2]))<<endl);
                    //MMBLOG_FILE_FUNC_LINE(DEBUG,  " myDensityMap.getDensity(Vec3(-.68, -4.79, 2.26)) = "<< myDensityMap.getDensity(Vec3(-.68, -4.79, 2.26))<<endl);
                    MMBLOG_FILE_FUNC_LINE(DEBUG,  " myDensityMap.getDensity(myAtomLocation) = "<< myDensityMap.getDensity(myAtomLocation)<<" myAtomLocation = "<<myAtomLocation<<" myDensityMap.getForceConstant()        = "<<myDensityMap.getForceConstant()<<" myAtomicNumber = "<<myAtomicNumber<<endl);
		    totalPotentialEnergy -= myDensityMap.getDensity(myAtomLocation) * myDensityMap.getForceConstant() * myAtomicNumber;
		    
                }
            } else {
                MMBLOG_FILE_FUNC_LINE(CRITICAL, " The chain ID you specified, "<< myChainID << " does not correspond to any existing BiopolymerClass or MonoAtoms !"<<endl);
            }
        } // of for myParameterReader.densityContainer.numDensityStretches
        MMBLOG_FILE_FUNC_LINE(INFO, " Total potential energy due to density fitting potential = "<<totalPotentialEnergy <<std::endl);
        return totalPotentialEnergy;
        };

bool DensityForce::dependsOnlyOnPositions() const  { 
        return true; 
    };    
