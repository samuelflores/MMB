/* -------------------------------------------------------------------------- *
 *                           MMB (MacroMoleculeBuilder)                       *
 * -------------------------------------------------------------------------- *
 *                                                                            *
 * Copyright (c) 2011-12 by the Author.                                       *
 * Author: Samuel Flores                                                      *
 *                                                                            *
 * See RNABuilder.cpp for the copyright and usage agreement.                  *
 * -------------------------------------------------------------------------- */

#ifndef MonoAtoms_H_
#define MonoAtoms_H_

#include <string>
#include "SimTKmolmodel.h"
//#include "DuMMForceFieldSubsystem.h"


using namespace std;
using namespace SimTK;

class MonoAtoms {
	private:
		String 	chainID;
    		ResidueID firstResidueID    ;
		//int 	numAtoms;
		String 	atomName;
		vector<Compound> compoundVector; 
	public:
		MonoAtoms();
		MonoAtoms(String, ResidueID, const int numAtoms, String);
		void            clear();
		void            addMonoAtom();
		void            addMonoAtom (Vec3 positionVec3 ); // Add a monoAtom at a specific default top level transform
                void            initialize(CompoundSystem & system,  bool readPreviousFrameFile,  String previousFrameFileName, bool matchExact, bool matchIdealized);
		int 		validate();  //return 0 if fine, return 1 and exit if not fine.
		void		validateResidue(ResidueID) ; //die if residue number out of range, or if atom names mismatch.
		String          getChainID();
		ResidueID	getLastResidueID ();
		ResidueID	getFirstResidueID ();
                ResidueID       getResidueID   (const int       myResidueIndex);
                int             getResidueIndex(ResidueID myResidueID);
		int             getNumAtoms() const;
		const String &  getAtomName() const;
		Compound 	getSingleCompound(ResidueID); 
		MobilizedBodyIndex getMobilizedBodyIndex(ResidueID);
		MobilizedBody &	updMobilizedBody(ResidueID, SimbodyMatterSubsystem &) ;
		bool            hasAtom(ResidueID residueNumber);
		Compound::AtomIndex getAtomIndex(ResidueID residueNumber);
		Vec3	        getAtomLocationInMobilizedBodyFrame(ResidueID residueNumber);
		Vec3	        getAtomLocationInGroundFrame       (const int       residueIndex ,const State & state);
		void		adoptCompounds(CompoundSystem & mySystem, bool readPreviousFrameFile);
		void		matchDefaultConfiguration(SimTK::PdbStructure,bool matchExact, bool matchIdealized);
                void            renumberPdbResidues();
                void            setPdbChainId(String chainID);
		String		getAtomPathName(ResidueID);
		void 		includeAllAtoms(DuMMForceFieldSubsystem & dumm);
		double          computeCurvatureSquared(const int index, const State & state);
		double          computeTotalCurvatureSquared(const State & state);

};
	
class MonoAtomsContainer  {
	private:
		map <String,MonoAtoms> monoAtomsMap;
	public:
		MonoAtomsContainer();
                void initialize(CompoundSystem &  , bool , String ,bool matchExact, bool matchIdealized  );
		//const int       getChainID(String);
		void            remove(String myChainID) ; 
		bool            hasChainID(const String &) const;
		MonoAtoms 	getMonoAtoms(const String &);
		const MonoAtoms & getMonoAtoms(const String &) const;
		void 		addMonoAtoms(MonoAtoms) ; 
		void		adoptCompounds(CompoundSystem & mySystem, bool readPreviousFrameFile);
		void 		matchDefaultConfiguration(SimTK::PdbStructure ,bool matchExact, bool matchIdealized );
		String		getAtomPathName(String,ResidueID);
		void 		clear();
		void 		includeAllAtoms(DuMMForceFieldSubsystem & dumm);
		double          computeTotalCurvatureSquared(const State & state);
};
#endif

