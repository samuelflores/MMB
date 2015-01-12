/* -------------------------------------------------------------------------- *
 *                           MMB (MacroMoleculeBuilder)                       *
 * -------------------------------------------------------------------------- *
 *                                                                            *
 * Copyright (c) 2011-12 by the Author.                                       *
 * Author: Samuel Flores                                                      *
 *                                                                            *
 * See RNABuilder.cpp for the copyright and usage agreement.                  *
 * -------------------------------------------------------------------------- */

#ifndef WaterDroplet_H_
#define WaterDroplet_H_
#include <vector>
#include <string>
#include "Water.h"
#include "Utils.h"
#include "AtomSpringContainer.h"
//#include "BiopolymerClass.h"
using namespace SimTK;
using namespace std  ;

class BiopolymerClassContainer;

class WaterDroplet {
private:
        double radius;
        vector <Transform> defaultTransformVector;
public:
        WaterDroplet();
	ResidueID firstResidueNumber;
	double getRadius();
	double setRadius(double myRadius);
	Vec3 center;
	string tetherType; // options: None, ToGround, ToAtom
	double tetherStrength; 
        string chainID;
	void printPDB();
	vector<Water> waterVector;
	int validate();
	void addWaterMolecules(CompoundSystem & system, DuMMForceFieldSubsystem &dumm,  SimbodyMatterSubsystem& matter, BiopolymerClassContainer & myBiopolymerClassContainer);
        MobilizedBody & updAtomMobilizedBody(SimbodyMatterSubsystem & matter, ResidueID myResidueNumber, string myAtomName);
        void multiplySmallGroupInertia(double smallGroupInertiaMultiplier, CompoundSystem & system, SimbodyMatterSubsystem& matter, State & state) ;
	void addTethers( AtomSpringContainer &   atomSpringContainer);
        MobilizedBodyIndex getOxygenMobilizedBodyIndex(ResidueID residueNumber);
        Vec3 getOxygenLocationInMobilizedBodyFrame(ResidueID residueNumber);
        MobilizedBody getOxygenMobilizedBody(SimbodyMatterSubsystem & matter, ResidueID residueNumber);
        MobilizedBody & updOxygenMobilizedBody(SimbodyMatterSubsystem & matter, ResidueID residueNumber);
        void includeAllAtoms(DuMMForceFieldSubsystem & dumm ); // add droplet to physics region
        void validateWaterVector( ); // add droplet to physics region
        void adopt( CompoundSystem & system, bool readPreviousFrameFile);

};

class WaterDropletContainer {
private:
     	vector <WaterDroplet> waterDropletVector;

public:
	void clear();
	void printPDB();
	void addWaterMolecules(CompoundSystem & system, DuMMForceFieldSubsystem &dumm,  SimbodyMatterSubsystem& matter, BiopolymerClassContainer & myBiopolymerClassContainer );
        void add(WaterDroplet &);
        void multiplySmallGroupInertia(double smallGroupInertiaMultiplier, CompoundSystem & system, SimbodyMatterSubsystem& matter, State & state) ;
	void addTethers(AtomSpringContainer & atomSpringContainer);
        WaterDroplet getWaterDroplet(string chainID);
        WaterDroplet & updWaterDroplet(string chainID);
	bool hasChainID (string chainID);
        void includeAllAtoms(DuMMForceFieldSubsystem & dumm);
        void validateWaterVectors( ); // add droplet to physics region
        void matchDefaultConfiguration(bool readPreviousFrameFile, String pdbFileName,bool matchExact, bool matchIdealized);
        void adopt( CompoundSystem & system, bool readPreviousFrameFile);
};

#endif
